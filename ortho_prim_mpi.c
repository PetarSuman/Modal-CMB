#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include "global.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <hdf5.h>


#define FILE_NAME "/nfs/st01/hpc-gr-epss/ps792/PlanckV3/C/gamma.h5"
#define DATASET_NAME "gamma_matrix"

/*
Primordial basis Q_n(k) isn't (necessarily) orthonormal. The interior product
between different basis elements is given by the \gamma matrix.

ortho_prim_mpi.c calculates          < Q_n(k) Q_m(k) >  :=  \gamma_{nm} (k) 

Its job submission script is in the wmpijob files.
*/

int main( int argc, char *argv[] ){
	
// Initialisation
	double time1, time2, time3, time4, time5, duration, duration1, duration2;

	if (argc < 2 || argc > 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	//	load relevant parameters from the parameter.ini file or alike
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);

	set_terms_prim();		//set primordial ordering n <---> ijk relating Q_n to q_iq_jq_k
	set_terms_late();		//set late time ordering n <---> ijk 
	
	int i,j,k,l,m,n;

	// MPI variables
	int rank, nproc, onext, enext;

	// MPI initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	double result=0;
	
	//	Load bessel and transfer functions from files specified by parameter.ini
	load_bessel();
	load_transfer();
	
	double kmax = get_kmax();	//calculate kmax of the targeted survey

	/*
	Alpha points (set in parameter.ini files) define the 
	number of points used in integration. The xvec array partitions the k-space
	into alpha_points number of intervals:
	*/
	double xvec[alpha_points+1];
	for(i=0;i<alpha_points+1;i++){
		xvec[i] = kmax*(double)i/alpha_points;
	}

	int ortho_size = get_terms_prim();

	/*	Print out all the mappings n <----> {ijk},
	there will be ortho_size number of them.
	*/
	if(rank==0){
		printf("kmax: %e\n", kmax);
		printf("pmax: %d\n", get_pmax_prim());
		printf("Primordial ordering:\n");
		printf("n\ti\tj\tk\n");
		for(n=0;n<ortho_size;n++){
			find_perm_prim(n,&i,&j,&k);
			printf("%d\t%d\t%d\t%d\n",n,i,j,k);
		}
	}

	/*	Create the primordial basis, defined by resolution, 
 	number of terms in the expansion and the order_prim
	*/	
	create_basis_prim(alpha_points+1,kmax,xvec);
	gsl_matrix *cholesky = gsl_matrix_alloc(ortho_size,ortho_size);
	gsl_matrix *matrixin = gsl_matrix_alloc(ortho_size,ortho_size);
	gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(ortho_size);
	gsl_vector *eigen = gsl_vector_alloc(ortho_size);
	gsl_matrix *eigenv = gsl_matrix_alloc(ortho_size,ortho_size);
	double x;
		
	create_ortho();
	int ortho_total = (ortho_size)*(ortho_size+1)/2;
    int nthreads = omp_get_max_threads();
	tetrapyd_limits* block;
	block = (tetrapyd_limits*) malloc(nthreads * sizeof(tetrapyd_limits));
	MPI_Barrier(MPI_COMM_WORLD);
	decompose_tetrapyd_prim(block, rank, nproc, nthreads, 0, alpha_points);
	
	for(i=0;i<nthreads;i++){
		printf("[%d:%d]\t%d\t%d\t%d\t%d\n",rank,i,block[i].i_bgn,block[i].j_bgn,block[i].k_bgn,block[i].loops);
	}
	double **ortho = (double **)create_array(ortho_size,ortho_size);
	double *ortho_send = (double *)malloc( ortho_total * sizeof(double) );
	double* ortho_recv = (double *)malloc( ortho_total * sizeof(double) );

	for (i=0;i<ortho_size;i++) {
		for (j=0;j<ortho_size;j++) {
			ortho[i][j] = 0e0;
		}
	}
	for(i=0;i<ortho_total;i++){
		ortho_send[i] = 0e0;
		ortho_recv[i] = 0e0;
	}
	
	// Calculate ortho computes <q_n q_m> = \gamma_nm
	calculate_ortho(0,alpha_points,block,ortho_send);
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(ortho_send,ortho_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		n=0;
		for (i=0;i<ortho_size;i++) {
			for (j=i;j<ortho_size;j++) {
				ortho[i][j] = ortho_recv[n];
				if(i!=j)ortho[j][i] = ortho[i][j];	//symmetrize
				n++;
				printf("%d\t%d\t%e\n",i,j,ortho[i][j]);
			}
		}
		/*
		//--------------- HDF5 OUTPUT CODE ----------------------------------------
		
		hid_t file_id, dataset_id, dataspace_id;
	  	int ROWS = ortho_size;
	  	int COLS = ortho_size;
  	  	hsize_t dims[2] = {ROWS, COLS};
  		int i, j;

	  	// Allocate memory for the 2D array
  		// Allocate memory for the flattened 1D array
  		double* data1d = (double*)malloc(ortho_size * ortho_size * sizeof(double));

  		// Flatten the 2D array into a 1D array
  		for (i = 0; i < ortho_size; i++) {
    		for (j = 0; j < ortho_size; j++) {
      			data1d[i * ortho_size + j] = ortho[i][j];
    		}
  		}

  		// Create a new HDF5 file
 		file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  		// Create the data space for the 2D array
  		dataspace_id = H5Screate_simple(2, dims, NULL);

  		// Create the dataset
  		dataset_id = H5Dcreate(file_id, DATASET_NAME, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  		// Write the flattened 1D array to the dataset
  		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);

  		// Close the HDF5 objects
  		H5Dclose(dataset_id);
  		H5Sclose(dataspace_id);
  		H5Fclose(file_id);

  		free(data1d);

		//---------------- END OF HDF5 CODE ---------------------------
		*/
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	double *results_o =  malloc( sizeof(double)*3);
	results_o[0] = 0.0;
	results_o[1] = 0.0;
	results_o[2] = 0.0;
		
	create_orthol();
	if (rank == 0) printf("Ortho l created.\n");
	create_lambdal();
	if (rank == 0) printf("Lambda l created.\n");

	if(rank==0){

		gsl_matrix_set_zero(cholesky);
		gsl_matrix_set_zero(matrixin);
		gsl_vector_set_zero(eigen);
		double norm[ortho_size];
			
		for (i=0;i<ortho_size;i++) {
			norm[i] = sqrt(ortho[i][i]);
			//printf("%d\t%e\n",i,norm[i]);
		}
	
		double **dummy = (double **)create_array(ortho_size,ortho_size);
		
		for (i=0;i<ortho_size;i++) {
			for (j=0;j<ortho_size;j++) {
				x = ortho[i][j]/(norm[i]*norm[j]);
				dummy[i][j] = x;
				gsl_matrix_set(cholesky,i,j,x);
				gsl_matrix_set(matrixin,i,j,x);
			}
		}
		
		printf("GSL Matrix set.\n");
		
		gsl_eigen_symmv(matrixin,eigen,eigenv,workspace);
		gsl_eigen_symmv_sort(eigen,eigenv,GSL_EIGEN_SORT_VAL_DESC);

		printf("GSL Eigen Symmv.\n");

		/*for(i=0;i<ortho_size;i++){
			printf("%d\t%e\n",i,gsl_vector_get(eigen,i));
		}*/

		gsl_linalg_cholesky_decomp(cholesky);

		printf("GSL Cholesky decomposition.\n");


		for(i=0;i<ortho_size;i++){
			for(j=0;j<ortho_size;j++){
				results_o[0] = (double)i;
				results_o[1] = (double)j;
	 			results_o[2] = norm[i]*gsl_matrix_get(cholesky, i, j);

				update_ortho(results_o);
			}
		}
		printf("Want to output ortho...\n");
		output_ortho();
		printf("Ortho updated and output.\n");
		
		
		gsl_matrix *upper = gsl_matrix_alloc(ortho_size,ortho_size);
		gsl_matrix *inverse = gsl_matrix_alloc(ortho_size,ortho_size);
		gsl_matrix_set_zero(upper);

		printf("GSL Upper and Inverse allocated.\n");

		
		for(i=0;i<ortho_size;i++){
			for(j=i;j<ortho_size;j++){
				x = gsl_matrix_get(cholesky, i, j);
				gsl_matrix_set(upper,i,j,x);
			}
		}
		printf("GSL Matrix set upper.\n");
		
		int s;
		gsl_permutation *perm = gsl_permutation_calloc(ortho_size);

		gsl_linalg_LU_invert(upper,perm,inverse);
		
		create_lambda();
		printf("Lamabda created.\n");
		
		for(i=0;i<ortho_size;i++){
			for(j=0;j<ortho_size;j++){
				results_o[0] = (double)i;
				results_o[1] = (double)j;
				results_o[2] = gsl_matrix_get(inverse, j, i)/norm[j];
				update_lambda(results_o);
				//printf("%d\t%d\t%e\n",i, j, results_o[2]);
			}
		}
		output_lambda();
	}

	printf("Code finished\n");
	
	MPI_Finalize();
}