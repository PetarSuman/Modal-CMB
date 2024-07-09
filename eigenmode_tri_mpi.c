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


int main( int argc, char *argv[] ){
	
// **1**

	printf("start\n");
	double time1, time2, time3, time4, time5, duration, duration1, duration2;


	if (argc < 3 || argc > 5) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile model (kstar) (phase)\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);

	set_terms_tri_prim();
	set_terms_tri_late();
	
	model = atoi(argv[2]);
	printf("model: %d\n",model);
	
	if(argc>3){
		kstar = atof(argv[3]);
	}else{
		kstar = 1e0;
	}
	printf("kstar: %e\n",kstar);
	
	
	if(argc>4){
		phase = atof(argv[4]);
	}else{
		phase = 0e0;
	}
	printf("phase: %e\n",phase);
	
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
		xvec[i] = kmax*(double)i/(double)alpha_points;
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
		
 	MPI_Barrier(MPI_COMM_WORLD);
 	
	read_ortho_tri();
	read_lambda_tri();
	
	double *results_e =  malloc( sizeof(double)*2);
	results_e[0] = 0.0;
	results_e[1] = 0.0;
	double c1,c2,c3,c4,x;
	double result2;
	
	gsl_matrix *cholesky = gsl_matrix_alloc(ortho_size,ortho_size);

	int eigen_size = ortho_size;
	double *eigen = (double *)create_vector(eigen_size);
	
	for (i=0;i<eigen_size;i++) {
		eigen[i] = 0.0;
	}
	
	int loops;
	int auxloop ;
	int start_loop;
	int end_loop;
	
	loops=eigen_size/nproc;
	auxloop = fmod(eigen_size,nproc);
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

// 	MPI_Barrier(MPI_COMM_WORLD);
	for(n=start_loop;n<end_loop+1;n++){
		time1 = csecond();
		if(eflag_order_prim!=5){
			eigen[n] = calculate_eigen_tri(enext);
		}else{
// 			eigen[n] = calculate_eigen_mc(enext);
			eigen[n] = calculate_eigen_tri(enext);
		}
		results_e[0] = (double)enext;
		results_e[1] = result;

		duration1 = csecond() - time1;
		printf("[%d]\t%d\t%e\t%e\t%e\n", rank, enext, result, result2, duration1);

	}

	// end of MPI loop
	printf("[%d] finished\n", rank);
		
	double* eigen_send = (double *)malloc( eigen_size * sizeof(double) );
	double* eigen_recv = (double *)malloc( eigen_size * sizeof(double) );
		
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (i=0;i<eigen_size;i++) {
		eigen_send[i] = eigen[i];
	}
	
	MPI_Reduce(eigen_send,eigen_recv,eigen_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		
		gsl_vector *eigenin = gsl_vector_alloc(eigen_size);
		gsl_vector *eigenout = gsl_vector_alloc(eigen_size);

		gsl_vector_set_zero(eigenin);
		gsl_vector_set_zero(eigenout);

		for(i=0;i<eigen_size;i++){
			x = eigen_recv[i];
			gsl_vector_set(eigenin, i, x);
		}

		for(i=0;i<ortho_size;i++){
			for(j=0;j<ortho_size;j++){
				x = get_ortho(i,j);
	 			gsl_matrix_set(cholesky, i, j, x);
			}
		}

		gsl_linalg_cholesky_solve(cholesky, eigenin, eigenout);

		printf("Alpha Q\n");
		for(i=0;i<eigen_size;i++){
			results_e[0] = (double)i;
			results_e[1] = gsl_vector_get(eigenout, i);
			update_eigen_tri(results_e);
			printf("[%d]\t%d\t%e\n", rank, i, results_e[1]);
		}

		output_eigen_tri();

		printf("Alpha R\n");
		for(i=0;i<eigen_size;i++){
			x=0.0;
			for(j=i;j<eigen_size;j++){
				x = x + get_ortho_tri(i,j)*get_eigen_tri(j);
			}
			results_e[0] = (double)i;
			results_e[1] = x;
			printf("[%d]\t%d\t%e\n", rank, i, results_e[1]);
		}

	}
		
	if(rank==0){
		create_modes_tri();
		read_gamma_tri();
		read_eigen_tri();

// 		printf("Alpha Q/R\n");
// 		for(i=0;i<ortho_size;i++){
// 			x=0.0;
// 			for(j=i;j<ortho_size;j++){
// 				x = x + get_ortho_tri(i,j)*get_eigen_tri(j);
// 			}
// 			printf("[%d]\t%d\t%e\t%e\n", rank, i, get_eigen_tri(i), x);
// 		}
		
		printf("Modes Q\n");
		for(i=0;i<ortho_size;i++){
			results_e[0] = (double)i;
			results_e[1] = 0.0;
			for(j=0;j<ortho_size;j++){
				x = get_gamma_tri(i,j);
				results_e[1] += x*get_eigen_tri(j);
			}
			update_modes_tri(results_e);
			printf("%d\t%e\t%e\n", i, get_modes_tri(i),get_eigen_tri(i));
		}

		output_modes_tri();

		read_orthol_tri();

		printf("Modes R\n");
		for(i=0;i<ortho_size;i++){
			x=0.0;
			for(j=0;j<ortho_size;j++){
				x += get_orthol_tri(j,i)*get_modes_tri(j);
			}
			printf("%d\t%e\n", i, x);
		}
	}
		
	MPI_Finalize();
}

