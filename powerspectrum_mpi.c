#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <gsl/gsl_spline.h>
#include "global.h"

int main( int argc, char *argv[] ){

	double pi = 3.141592653589793;
	
// **1**

	double time1, time2, time3, time4, time5, duration;


	
	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s dir\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	
	load_transfer();

// ini file
	int  lcol = 0;
	int restart = 0;
	int lsize, ltotal;
	int i,j,k,n,l;
	double area, result;
	double x1,x2,x3,x4,x5,x6;

// mpi vars
	int rank, nproc, lnext, bnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	printf("[%d] Process[%d]  \n", rank, getpid());


// load l's from ini file or cli 
	
	int ksize = get_ksize();
	double kvalues[ksize];
	double transferT[ksize];
	double transferE[ksize];
	get_kvec(kvalues);
	double kmin,kmax;
	kmin = kvalues[0];
	kmax = kvalues[ksize-1];

	printf("kmax %e\n",kmax);
	
	double integrandTT[ksize];
	double integrandTE[ksize];
	double integrandEE[ksize];
	double pt;
	l=0;
	double tau0 = 14.6e3;
	
	gsl_spline* splTT =  gsl_spline_alloc (gsl_interp_cspline, ksize);
	gsl_spline* splTE =  gsl_spline_alloc (gsl_interp_cspline, ksize);
	gsl_spline* splEE =  gsl_spline_alloc (gsl_interp_cspline, ksize);
	gsl_interp_accel* accTT = gsl_interp_accel_alloc();
	gsl_interp_accel* accTE = gsl_interp_accel_alloc();
	gsl_interp_accel* accEE = gsl_interp_accel_alloc();
	
	integrandTT[0] = 0;
	integrandTE[0] = 0;
	integrandEE[0] = 0;
	// for(i=0;i<lsize;i++){
		// l = lvalues[i];
	for(l=2;l<5001;l++){
		get_tvec_T(l, transferT);
		get_tvec_E(l, transferE);
		// x4 = sqrt((l+2e0)*(l+1e0)*l*(l-1e0));
		x4 = 1e0;
		for(j=1;j<ksize;j++){
// 			get_tvec(l, transfer);
			pt = kpivot* pow(kvalues[j]/kpivot,2.0-nscalar);
// 			pt = pow(kvalues[j],2.0-nscalar);
// 			if(kvalues[j] > 0.*l/tau0){
			
			integrandTT[j] = transferT[j] * transferT[j] / (pt);
			integrandTE[j] = x4*transferT[j] * transferE[j] / (pt);
			integrandEE[j] = x4*x4*transferE[j] * transferE[j] / (pt);
// 			}else{
// 				integrand[j]=0.0;
// 			}
// 			integrand[j] = pow(spherical_bessel(l,tau0*kvalues[j]),2) /(9.0*pt);
			// if(l==20)printf("20\t%e\t%e\t%e\n",kvalues[j],transferT[j],x4*transferE[j]);
			// if(l==200)printf("200\t%e\t%e\t%e\n",kvalues[j],transferT[j],x4*transferE[j]);
			// if(l==2000)printf("2000\t%e\t%e\t%e\n",kvalues[j],transferT[j],x4*transferE[j]);
		}
		result=0.0;
		gsl_spline_init(splTT,kvalues,integrandTT,ksize);
		gsl_spline_init(splTE,kvalues,integrandTE,ksize);
		gsl_spline_init(splEE,kvalues,integrandEE,ksize);
// 		result = M_2_PI * deltaphi * pow(kpivot,1.0-nscalar)* gsl_spline_eval_integ(spl,kmin,kmax,acc);
// 		result = M_2_PI * 1.55e-8* 1.11 * gsl_spline_eval_integ(spl,kmin,kmax,acc);
		x1 = M_2_PI * gsl_spline_eval_integ(splTT,kmin,kmax,accTT);
		x2 = M_2_PI * gsl_spline_eval_integ(splTE,kmin,kmax,accTE);
		x3 = M_2_PI * gsl_spline_eval_integ(splEE,kmin,kmax,accEE);
// 		result = M_2_PI * gsl_spline_eval_integ(spl,kmin,kmax,acc);
		printf("%d\t%e\t%e\t%e\n",l,x1,x2,x3);
	}
	
	MPI_Finalize();
	
	// End of code
	return 0;
}