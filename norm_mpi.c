#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
// #include <nag.h>
// #include <nag_stdlib.h>
// #include <nagf02.h>
// #include <nagf07.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include "global.h"

int main( int argc, char *argv[] ){

// set up timing
	double time1, time2, time3, time4, time5, duration;
	
	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s dir\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	time1 = csecond();

// mpi vars
	int rank, nproc, lnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);


// read in data
	
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");
	
	set_terms_prim();
	set_terms_late();

	read_beta(directory);
	int terms = get_terms_prim();
	int lsize = get_b_lsize();
	int *lvec = create_ivector(lsize);
	get_b_lvec(lvec);
	create_cl(lsize);
	create_beam(lsize);
	create_noise(lsize);
	load_cl(directory, lsize, lvec);
	load_BN(directory, lsize, lvec);
	int lmax = lvec[lsize-1];
	init_lmax(lmax);
		
	int i,j,k,r,s,t,l,m,n;
	double sum1,sum2,sum3,sum4;
	double x1,x2,x3,x4;
	
	if(rank==0){
		printf("lmax: %d\n", (int)lmax);
		printf("prim pmax: %d\n", get_pmax_prim());
		printf("late pmax: %d\n", get_pmax_late());
		printf("prim terms: %d\n", get_terms_prim());
		printf("late terms: %d\n", get_terms_late());
		printf("xsize: %d\n",get_b_xsize());
		for(n=0;n<terms;n++){
			find_perm_prim(n,&i,&j,&k);
			find_perm_late(n,&r,&s,&t);
			printf("%d\t\t%d\t%d\t%d\t\t%d\t%d\t%d\n",n,i,j,k,r,s,t);
		}
	}
	
	int xsize = get_b_xsize();
	double *xvec = create_vector(xsize);
	get_b_xvec(xvec);
	double *yvec = create_vector(xsize);
	init_norm_glint();

	int loops=xsize/nproc;
	int auxloop = fmod(xsize,nproc);
	int start_loop = rank*loops;
	int end_loop = (rank+1)*loops-1;

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

	for (i=0;i<xsize;i++) {
		yvec[i] = 0.0;
	}
	
	time1 = csecond();
	for(n=start_loop;n<end_loop+1;n++){
		
		yvec[n] = calculate_norm(0,n);
		
		duration = csecond() - time1;
		printf("[%d]\t%d\t%e\t%e\t%e\n", rank, n, xvec[n], yvec[n], duration);
		time1 = csecond();
		
	} // end of MPI loop
	printf("[%d] finished\n", rank);
	
	double* norm_send = (double *)malloc( xsize * sizeof(double) );
	double* norm_recv = (double *)malloc( xsize * sizeof(double) );
	
	MPI_Barrier(MPI_COMM_WORLD);

	for (i=0;i<xsize;i++) {
		norm_send[i] = yvec[i];
	}
	
        MPI_Reduce(norm_send,norm_recv,xsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		for (i=0;i<xsize;i++) {
			yvec[i] = norm_recv[i];
			printf("[%d]\t%d\t%e\t%e\n", rank, i, xvec[i], yvec[i]);
		}

		double xmin = xvec[0];
		double xmax = xvec[xsize-1];
		double norm;
	
		gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, xsize);
		gsl_interp_accel* acc = gsl_interp_accel_alloc();
		gsl_spline_init(sp,xvec,yvec,xsize);
		norm = gsl_spline_eval_integ(sp,xmin,xmax,acc);
		gsl_spline_free(sp);
		gsl_interp_accel_free(acc);

		printf("Result: %e\n",norm);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

