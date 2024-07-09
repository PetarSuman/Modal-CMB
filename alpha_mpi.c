#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <nagmk21.h>
#include "global.h"

int main( int argc, char *argv[] ){

	double pi = 3.141592653589793;
	
// **1**

	double time1, time2, time3, time4, time5, duration;


	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s data\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");
	
	load_transfer(directory);
	
// ini file
	char inifile[MAXLEN];
	int i,j,k,n,m;
	double result;

// mpi vars
	int rank, nproc, lnext, knext, anext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	printf("[%d] Process[%d]  \n", rank, getpid());

// Calculate alpha

	int alpha_total = (alpha_max+1)*(alpha_max+2)*(alpha_max+3)/6;
	int alpha_triples[alpha_total][3];
		
	n=0;	
	for (i=0;i<alpha_max+1;i++) {
		for (j=i;j<alpha_max+1;j++) {
			for (k=j;k<alpha_max+1;k++) {
				alpha_triples[n][0] = i;
				alpha_triples[n][1] = j;
				alpha_triples[n][2] = k;
				n++;
			}
		}
	}
	
	double *results =  malloc( sizeof(double)*4);
	results[0] = 0.0;
	results[1] = 0.0;
	results[2] = 0.0;
	results[3] = 0.0;
		
	double kmax = get_kmax();
	
	int ksize = alpha_points + 1;
	double kvec[ksize];
	
	for (i=0;i<ksize;i++){
		kvec[i] = (kmax)*(double)i/alpha_points;
		//printf("%e\n", kvec[i]);
	}
	
	create_basis(ksize, kvec);
	printf("[%d] Created basis\n",rank);
	
	int spl_size = ksize + 4;
	int wrk_size = 6 * ksize + 16;
	
	double x[ksize];
	double y[ksize];
	double spl_k[spl_size];
	double spl_c[spl_size];
	double wrk[wrk_size];
	int ifail;
	
	for (i=0;i<ksize;i++) {
		x[i] = kvec[i]/kmax;
	}
	
	y[0] = 1;
	
	if ( rank == 0 ){
	double sum;
		for (i=0;i<alpha_max+1;i++) {
			sum = 0;
			for (j=0;j<ksize;j++) {
				y[j] = get_basis(j,i)*get_basis(j,i)*(2*i+1);
			}
			e01baf_(&ksize,x,y,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
			e02bdf_(&spl_size,spl_k,spl_c,&sum,&ifail);
			printf("%d\t%e\n", i, 1.0-sum);
		}
	}
	
	
	if ( rank == 0 ) record_tasks(alpha_total);
	
	create_alpha();
	
	MPI_Barrier(MPI_COMM_WORLD);
	time1 = csecond();
	
	while ( (anext = get_next_task(1,4,rank,results)) < alpha_total ) {
		
		i = alpha_triples[anext][0];
		j = alpha_triples[anext][1];
		k = alpha_triples[anext][2];
		
		result = calculate_alpha(rank, kvec, i, j, k);
		results[0] = (double)i;
		results[1] = (double)j;
		results[2] = (double)k;
 		results[3] = result;
		
		duration = csecond() - time1;
		//printf("[%d] atask: %d result: %d %d %d %e time: %e\n", rank, anext, i, j, k, result, duration);
		printf("%d %d %d %e %e\n", i, j, k, result, duration);
		time1 = csecond();
		
		
	} // end of MPI loop
	
	if (rank == 0){
		signal_end_tasks(1,4);
		output_alpha(directory);
	}
	
	destroy_basis(ksize);
	ksize = get_ksize();
	double kvalues[ksize];
	get_kvec(kvalues);
	create_basis(ksize, kvalues);
	printf("[%d] Created basis\n",rank);
	
	double prim;
	double kval, orig;
	
	MPI_Barrier(MPI_COMM_WORLD);	
	
	if (rank == 0){
		printf("\nCentre\n");
		for (n=0;n<ksize;n++){
			prim = 0;
			for (i=0;i<alpha_max+1;i++){
				for (j=0;j<alpha_max+1;j++){
					for (k=0;k<alpha_max+1;k++){
					prim += get_alpha(i, j, k) * get_basis(n,i) * get_basis(n,j) * get_basis(n,k);				
					}
				}
			}
			kval = kvalues[n];
			orig = shape3(kval, kval, kval);
			result = prim / orig;
			printf("%e %e\n", kval, result);	
		}
		
		for (n=0;kvalues[n]<kmax/2.0;n++) m=n;
		
// 		printf("\nDiagonal\n");
// 		
// 		for (n=0;kvalues[n]<ksize;n++){
// 			prim = 0;
// 			for (i=0;i<alpha_max+1;i++){
// 				for (j=0;j<alpha_max+1;j++){
// 					for (k=0;k<alpha_max+1;k++){
// 					prim += get_alpha(i, j, k) * get_basis(n,i) * get_basis(m,j) * get_basis(m,k);				
// 					}
// 				}
// 			}
// 			kval = kvalues[n];
// 			orig = shape3(kval, kvalues[m], kvalues[m]);
// 			result = prim / orig;
// 			if(kval>2*kvalues[m])break;
// 			printf("%e %e\n", kval, result);	
// 		}
	}
	
	printf("[%d] finished\n", rank);
	
	MPI_Finalize();
	
	// End of code
	return 0;
}