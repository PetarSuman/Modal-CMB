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

int main( int argc, char *argv[] ){
	
// **1**

	double time1, time2, time3, time4, time5, duration, duration1, duration2;


	if (argc < 2 || argc > 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	
	set_terms_tri_prim();
	set_terms_tri_late();
	
	int i,j,k,l,m,n;

// mpi vars
	int rank, nproc, onext, enext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	double result=0;
	
	load_bessel();
	load_transfer();
	
	double kmax = get_kmax_cut();

	double xvec[alpha_points+1];
	for(i=0;i<alpha_points+1;i++){
		xvec[i] = kmax*(double)i/alpha_points;
	}

	int ortho_size = get_terms_prim();
	if(rank==0){
		printf("kmax: %e\n", kmax);
		printf("pmax: %d\n", get_pmax_prim());
		for(n=0;n<ortho_size;n++){
			find_perm_tri_prim(n,&i,&j,&k,&l);
			printf("%d\t%d\t%d\t%d\t%d\n",n,i,j,k,l);
		}
	}

	create_basis_tri_prim(alpha_points+1,kmax,xvec);
	
// 	if(rank==0){
// 		for(i=0;i<alpha_points+1;i++){
// 			printf("%e",xvec[i]);
// 			for(j=0;j<get_pmax_prim()+1;j++){
// 				printf("\t%e",get_basis_prim(i,j));
// 			}
// 			printf("\n");
// 		}
// 	}
	
	int loops;
	int auxloop ;
	int start_loop;
	int end_loop;
	
	gsl_matrix *cholesky = gsl_matrix_alloc(ortho_size,ortho_size);
	gsl_matrix *matrixin = gsl_matrix_alloc(ortho_size,ortho_size);
	gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(ortho_size);
	gsl_vector *eigen = gsl_vector_alloc(ortho_size);
	gsl_matrix *eigenv = gsl_matrix_alloc(ortho_size,ortho_size);
	double x;
		
	create_ortho();

	int ortho_total = (ortho_size)*(ortho_size+1)/2;
	int ortho_pairs[ortho_total][2];
	double **ortho = (double **)create_array(ortho_size,ortho_size);
		
	n=0;
	for (i=0;i<ortho_size;i++) {
		for (j=i;j<ortho_size;j++) {
			ortho_pairs[n][0] = i;
			ortho_pairs[n][1] = j;
			n++;
		}
	}
	
	loops=ortho_total/nproc;
	auxloop = fmod(ortho_total,nproc);
	start_loop = rank*loops;
	end_loop = (rank+1)*loops-1;
	
	if (auxloop != 0){
		if (rank < auxloop){
			start_loop = start_loop + rank;
			end_loop = end_loop + rank + 1;
		}else{
			start_loop = start_loop + auxloop;
			end_loop = end_loop + auxloop;
		}
	}
		
	printf("mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
	
	for (i=0;i<ortho_size;i++) {
		for (j=0;j<ortho_size;j++) {
			ortho[i][j] = 0.0;
		}
	}
	
	time1 = csecond();
	for(n=start_loop;n<end_loop+1;n++){
		
		i = ortho_pairs[n][0];
		j = ortho_pairs[n][1];
		ortho[i][j] = calculate_ortho_tri(i,j);
			
		duration1 = csecond() - time1;
		printf("[%d]\t%d\t%d\t%e\t%e\t%e\n", rank, i, j, ortho[i][j], x, duration1);
		time1 = csecond();
		
	} // end of MPI loop
	printf("[%d] finished\n", rank);
		
	double* ortho_send = (double *)malloc( ortho_total * sizeof(double) );
	double* ortho_recv = (double *)malloc( ortho_total * sizeof(double) );
		
	MPI_Barrier(MPI_COMM_WORLD);
	
	n=0;
	for (i=0;i<ortho_size;i++) {
		for (j=i;j<ortho_size;j++) {
			ortho_send[n] = ortho[i][j];
			n++;
		}
	}
	
	MPI_Reduce(ortho_send,ortho_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		n=0;
		for (i=0;i<ortho_size;i++) {
			for (j=i;j<ortho_size;j++) {
				ortho[i][j] = ortho_recv[n];
				if(i!=j)ortho[j][i] = ortho[i][j];
				n++;
			}
		}
		
// 		char filename[100] = "/home/cosmos/tmp/jf334/DX9/master_ortho_prim_5_6000.unf";
// 		int big_size = ortho_size*ortho_size;
// 		array_write(&big_size, filename, &ortho[0][0]);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	double *results_o =  malloc( sizeof(double)*3);
	results_o[0] = 0.0;
	results_o[1] = 0.0;
	results_o[2] = 0.0;
		
	create_ortho_tri();
	create_lambda_tri();

	if(rank==0){

		gsl_matrix_set_zero(cholesky);
		gsl_matrix_set_zero(matrixin);
		gsl_vector_set_zero(eigen);
		double norm[ortho_size];
			
		for (i=0;i<ortho_size;i++) {
			norm[i] = sqrt(ortho[i][i]);
			printf("%d\t%e\n",i,norm[i]);
		}
	
		for (i=0;i<ortho_size;i++) {
			for (j=0;j<ortho_size;j++) {
// 				x = ortho[i][j]/(norm[i]*norm[j]);
				x = ortho[i][j];
				gsl_matrix_set(cholesky,i,j,x);
				gsl_matrix_set(matrixin,i,j,x);
			}
		}
		
		gsl_eigen_symmv(matrixin,eigen,eigenv,workspace);
		gsl_eigen_symmv_sort(eigen,eigenv,GSL_EIGEN_SORT_VAL_DESC);

		for(i=0;i<ortho_size;i++){
			printf("%d\t%e\n",i,gsl_vector_get(eigen,i));
		}

		gsl_linalg_cholesky_decomp(cholesky);

		for(i=0;i<ortho_size;i++){
			for(j=0;j<ortho_size;j++){
				results_o[0] = (double)i;
				results_o[1] = (double)j;
// 	 			results_o[2] = norm[i]*gsl_matrix_get(cholesky, i, j);
	 			results_o[2] = gsl_matrix_get(cholesky, i, j);
				update_ortho_tri(results_o);
			}
		}
		output_ortho_tri();
		
		
		gsl_matrix *upper = gsl_matrix_alloc(ortho_size,ortho_size);
		gsl_matrix *inverse = gsl_matrix_alloc(ortho_size,ortho_size);
		gsl_matrix_set_zero(upper);

		
		for(i=0;i<ortho_size;i++){
			for(j=i;j<ortho_size;j++){
				x = gsl_matrix_get(cholesky, i, j);
				gsl_matrix_set(upper,i,j,x);
			}
		}
		
		int s;
		gsl_permutation *perm = gsl_permutation_calloc(ortho_size);

		gsl_linalg_LU_invert(upper,perm,inverse);
		
		create_lambda_tri();
		
		for(i=0;i<ortho_size;i++){
			for(j=0;j<ortho_size;j++){
				results_o[0] = (double)i;
				results_o[1] = (double)j;
// 				results_o[2] = gsl_matrix_get(inverse, j, i)/norm[j];
				results_o[2] = gsl_matrix_get(inverse, j, i);
				update_lambda_tri(results_o);
			}
		}
		output_lambda_tri();
		
	}
	
	MPI_Finalize();
}

