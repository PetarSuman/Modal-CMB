#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include "global.h"

int main( int argc, char *argv[] ){

	double pi = 3.141592653589793;
	
// **1**

	double time1, time2, time3, time4, time5, duration;


	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	
	load_bessel();
	load_transfer();

// ini file
	int* lvalues = malloc( sizeof(int)*MAXLINES);
	int* lceil = malloc( sizeof(int));
	int  lcol = 0;
	int restart = 0;
	int lsize, ltotal;
	int i,j,k,l,n;
	double area, result;

// mpi vars
	int rank, nproc, lnext, bnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	printf("[%d] Process[%d]  \n", rank, getpid());


// load l's from ini file or cli 


	if(use_l_file==1){
		if ( load_one(inifile, lvalues, lceil ) ) {			
			printf("[%d] Ini load complete\n", rank);
		} else {
			printf("[%d] Ini file error\n", rank);
			MPI_Finalize();
			exit(1);
		}
		lsize = *lceil;
	}else{
		lsize = eflag_lmax+1;
		for(i=0;i<lsize;i++){
			lvalues[i] = i;
		}
	}
	

	set_terms_tri_prim();
	int lmax = lvalues[lsize-1];
	double kmax = get_kmax_cut();
	int ksize = get_ksize();
	double kvalues[ksize];
	get_kvec(kvalues);
	create_basis_tri_prim(ksize, kmax, kvalues);
	int pmax = get_pmax_prim();
	double prim;
	double kval, orig;
	
// 	int bbxsize = get_xsize();
// 	int bbxsize = ksize;
// 	double bbx[bbxsize];
// 	double bbb[bbxsize];
// 	get_xvec(bbx);
	
	if(rank==0){
// 		for (i=1990;i<lsize;i++){
// 			get_tvec(i, bbb);
// 			printf("%d\n",i);
// 			for(j=0;j<bbxsize;j++){
// 				printf("%e\t%e\n",kvalues[j],bbb[j]);
// 			}
// 			printf("\n");
// 		}
	
		printf("kmax_cut: %e\n", kmax);
		for(n=0;n<get_terms_prim();n++){
			find_perm_tri_prim(n,&i,&j,&k,&l);
			printf("%d\t%d\t%d\t%d\t%d\n",n,i,j,k,l);
		}
		
		for(i=0;i<ksize;i++){
			printf("%e",kvalues[i]);
			for(n=0;n<pmax;n++){
				printf("\t%e",get_basis_tri_prim(i,n));
			}
			printf("\n");
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
// calculate k integrals
	
	double xstep;

	double step3 = 200.0;
	double step2 = 50.0;
	double step1 = 10.0;
	double step0 = 5.0;

// 	double step3 = 60.0;
// 	double step2 = 30.0;
// 	double step1 = 15.0;
// 	double step0 = 3.0;
	
	double cut1 = 13000.0;
	double cut2 = 13600.0;
	double cut3 = 14200.0;
	double cut4 = 15000.0;

	double x=0.0;
	double xmax = 2*get_tau0();
	
	for(n=0;x<=xmax;n++){
		if(x<=cut1){
			xstep = step3;
		}else if (x<=cut2){
			xstep = step2;
		}else if (x<=cut3){
			xstep = step1;
		}else if (x<=cut4){
			xstep = step2;
		}else{
			xstep = step3;
		}
		x+=xstep;
	}
	int xsize=n;
	printf("[%d] xsize %d\n",rank,xsize);
	double xvec[xsize];
	x=0.0;
	xvec[0]=0.0;
	for(n=0;x<=xmax;n++){
		if(x<=cut1){
			xstep = step3;
		}else if (x<=cut2){
			xstep = step2;
		}else if (x<=cut3){
			xstep = step1;
		}else if (x<=cut4){
			xstep = step2;
		}else{
			xstep = step3;
		}
		x+=xstep;
		xvec[n]=x;
	}
	
	double *results =  malloc( sizeof(double)*(2+xsize));
	for(i=0;i<2+xsize;i++) results[i] = 0.0;
	
	int ktotal = (pmax+1) * lsize;
	double ktriples[ktotal][3];
	n=0;
	for (j=0;j<lsize;j++){
		for (k=0;k<pmax+1;k++){
			ktriples[n][0] = (double)lvalues[j];
			ktriples[n][1] = (double)j;
			ktriples[n][2] = (double)k;
			n++;
		}
	}
	
	int alpha;

	// master record total number of tasks	
	if ( rank == 0 ) record_tasks(ktotal);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	create_beta_tri(lsize, xsize);

	// mpi loop over the kvalues
	
	time1 = csecond();
	while ( (bnext = get_next_task(12,2+xsize,rank,results)) < ktotal ) {
		
		time1 = csecond();
		
		l = (int)ktriples[bnext][0];
		alpha = (int)ktriples[bnext][2];
		results[0] = ktriples[bnext][1];
		results[1] = ktriples[bnext][2];
		
		for(i=0;i<xsize;i++){
			x = xvec[i];
			if(l>=2){
				results[2+i] = M_2_PI * calculate_kint_tri(l, alpha, x);
			}else{
				results[2+i] = 0.0;
			}
		}
		
		duration = csecond() - time1;
		
		printf("[%d] Done:\t%d\t%d\tin: %es\n", rank, l, alpha, duration);
		
	} // end of MPI loop
	
	if (rank == 0){
		signal_end_tasks(12,3+xsize);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank == 0){
		output_beta_tri(xsize,xvec,lsize,lvalues);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	destroy_beta_tri();
	MPI_Barrier(MPI_COMM_WORLD);
	
	time1 = csecond();
	read_beta_tri();
	time2 = csecond();
	printf("[%d] Read beta in %e\n", rank, time2-time1);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
	
		int size = get_terms_prim();
		for (j=0;j<lsize;j++){
			printf("%d\t",j);
			x = deltaphi*deltaphi*calculate_xint_tri(j, j, j, j, 0, xsize, xvec)/(4.0*M_PI);
			printf("%e\n",x);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	printf("[%d] Finished code\n", rank);

	MPI_Finalize();
	
	// End of code
	return 0;
}