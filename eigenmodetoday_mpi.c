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
// #include <nagmk21.h>
// #include <nag.h>
// #include <nagf03.h>
#include "tetrapyd_tools.h"

int main( int argc, char *argv[] ){

	double time1, time2, time3, time4, time5, duration;

 	if (argc != 3) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s data mode\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	
	// MPI Initialization
	int rank, nproc, onext, mnext;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);


	//	Obtain the mode: 0 = ISW, 1 = point sources
	int mode = atoi(argv[2]);
	if (rank == 0) printf("mode: %d\n",mode);

	int i,j,k,l,m,n,r,s;
	int l_size = eflag_T_lmax+1;
	int *l_values = create_ivector(l_size);
	for(i=0;i<l_size;i++){
		l_values[i] = i;
	}
	
	double lmax = (double)l_values[l_size-1];
	init_lmax((int)lmax);
	
	long int l_size_long = l_size;
	set_terms_prim();
	// find_perm_prim(4,&i,&j,&k);
	// printf("prim: %d\t%d\t%d\n",i,j,k);
	set_terms_late();
	load_cl();
	load_BN();
	load_lens();
	if(eflag_order_late==5)load_TL(l_size, l_values);
	// int xsize = (int)lmax+1;

	double Tmax = (double)eflag_T_lmax;
	double Tmin = (double)eflag_T_lmin;
	int Tsize = Tmax-Tmin+1;  
	double Tvec[Tsize];
	for(i=0;i<Tsize;i++)Tvec[i] = (double)i+Tmin;
	
	double Emax = (double)eflag_E_lmax;
	double Emin = (double)eflag_E_lmin;
	int Esize = Emax-Emin+1;
	double Evec[Esize];
	for(i=0;i<Esize;i++)Evec[i] = (double)i+Emin;
	
	create_basis_late(Tsize, Esize, Tmin, Tmax, Emin, Emax, Tvec, Evec);
	
	double result=0;
	double x,x1,x2,x3,x4;
	
	int ortho_size = get_terms_late();
	
	if(rank==0){
		printf("lmax: %d\n", (int)lmax);
		printf("pmax: %d\n", get_pmax_late_T());
		printf("terms: %d\n", get_terms_late());
		for(n=0;n<ortho_size;n++){
			find_perm_late_TTT(n,&i,&j,&k);
			printf("%d\t%d\t%d\t%d\n",n,i,j,k);
		}
		// for(n=0;n<eflag_T_lmax-eflag_T_lmin+1;n++){
		// 	x1 = get_basis_late_T(n,0);
		// 	x2 = get_basis_late_T(n,1);
		// 	x3 = get_basis_late_T(n,2);
		// 	x4 = get_basis_late_T(n,3);
		// 	printf("%d\t%e\t%e\t%e\t%e\n",n,x1,x2,x3,x4);
		// }
		// for(n=0;n<eflag_T_lmax+1;n++){
		// 	x1 = plijk_TTT(0,n,n,n);
		// 	x2 = plijk_TTT(1,n,n,n);
		// 	x3 = plijk_TTT(2,n,n,n);
		// 	x4 = plijk_TTT(3,n,n,n);
		// 	printf("%d\t%e\t%e\t%e\t%e\n",n,x1,x2,x3,x4);
		// }
	}
	int loops;
	int auxloop;
	int start_loop;
	int end_loop;

 	MPI_Barrier(MPI_COMM_WORLD);

	read_lambdal();
	read_qtilde();
	read_orthol();
	// int xsize = get_qtilde_xsize();
	//     double *xvec = create_vector(xsize);
	// get_qtilde_xvec(xvec);
	// for(i=0;i<xsize;i++)printf("xvec\t%d\t%e\n",i,xvec[i]);
			
	double **orthoinv = (double **)create_array(ortho_size,ortho_size);
	double **lambda = (double **)create_array(ortho_size,ortho_size);
	double **lambdainv = (double **)create_array(ortho_size,ortho_size);

	/*
	for(i=0;i<ortho_size;i++){
		for(j=0;j<ortho_size;j++){
			if(i<j){
				lambdainv[i][j] = 0.0;
			}
			else{
				lambdainv[i][j] = get_orthol_TTT(i,j);
			}
		}
	}
	for(i=0;i<ortho_size;i++){
		for(j=0;j<ortho_size;j++){
			orthoinv[i][j]=0.0;
			for(k=0;k<ortho_size;k++){
				orthoinv[i][j] += get_lambdal_TTT(k,i)*get_lambdal_TTT(k,j);
			}
		}
	}
	*/

	//	Petar's attempt to speed up by iterating the loops only once:
	for(i=0;i<ortho_size;i++){
		for(j=0;j<ortho_size;j++){
			if(i<j){
				lambdainv[i][j] = 0.0;
			}
			else{
				lambdainv[i][j] = get_orthol_TTT(i,j);
			}

			orthoinv[i][j]=0.0;
			for(k=0;k<ortho_size;k++){
				orthoinv[i][j] += get_lambdal_TTT(k,i)*get_lambdal_TTT(k,j);
			}
		}
	}	

	int modes_size = ortho_size;
	create_modes();
	
	double *results_m =  malloc( sizeof(double)*2);
	results_m[0] = 0.0;
	results_m[1] = 0.0;
	
	init_ISW_glint();
	init_PS_glint();
	init_3D(mode);

	// update_modes(results_m);
	//
	// if ( rank == 0 ) record_tasks(modes_size);
	//
	// time1 = csecond();
	//
	// double c1,c2,c3,c4;
	//
	//
	// printf("[%d] Starting\n", rank);
	// MPI_Barrier(MPI_COMM_WORLD);
	//
	// MPI_Barrier(MPI_COMM_WORLD);
	//
	// loops=modes_size/nproc;
	// auxloop = fmod(modes_size,nproc);
	// start_loop = rank*loops;
	// end_loop = (rank+1)*loops-1;
	//
	// if (auxloop != 0){
	// 	if (rank < auxloop){
	// 		start_loop = start_loop + rank;
	// 		end_loop = end_loop + rank + 1;
	// 	}else{
	// 		start_loop = start_loop + auxloop;
	// 		end_loop = end_loop + auxloop;
	// 	}
	// }
	//
	// double* modes_send = (double *)malloc( modes_size * sizeof(double) );
	// double* modes_recv = (double *)malloc( modes_size * sizeof(double) );
	//
	// for(i=0;i<modes_size;i++){
	// 	modes_send[i] = 0.0;
	// 	modes_recv[i] = 0.0;
	// }
	//
	// time1 = csecond();
	// for(n=start_loop;n<end_loop+1;n++){
	//
	// 	if(mode==0){
	// 		c2 = calculate_ISW_3D(n);
	// 	}elseif(mode==1){
	// 		c2 = calculate_PS_3D(n);
	// 	}else{
	// 		c2 = calculate_mode_3D(n);
	// 	}
	// 	modes_send[n]=c2;
	//
	// 	duration = csecond() - time1;
	// 	printf("[%d]\t%d\t%e\t%e\t%e\n", rank, n, c1, c2, duration);
	// 	time1 = csecond();
	//
	// }
	//
	// printf("[%d] Finished\n", rank);
	// MPI_Reduce(modes_send,modes_recv,modes_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);


	
    int nthreads = omp_get_max_threads();
	tetrapyd_limits* block;
	block = (tetrapyd_limits*) malloc(nthreads * sizeof(tetrapyd_limits));
	decompose_tetrapyd_XXX(block, rank, nproc, nthreads, 2, lmax);
	
	for(i=0;i<nthreads;i++){
		printf("[%d:%d]\t%d\t%d\t%d\t%d\n",rank,i,block[i].i_bgn,block[i].j_bgn,block[i].k_bgn,block[i].loops);
	}
	
	double *modes = (double *)create_vector(modes_size);
	
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Check 7\n");	
	if(mode==0){ 
		if (rank == 0) printf("Calculate ISW\n");
		calculate_ISW_3D(2,lmax,block,modes);
	}
	else if (mode==1){
		if (rank == 0) printf("Calculate Point Source\n");
		calculate_PS_3D(2,lmax,block,modes);
		// for (i=0;i<modes_size;i++) {
			// printf("[%d] modes 1:\t%d\t%e\n",rank,i,modes[i]);
		// }
	}
	else{
		calculate_mode_3D(2,lmax,block,modes);
	}
	
	double* modes_send = (double *)malloc( modes_size * sizeof(double) );
	double* modes_recv = (double *)malloc( modes_size * sizeof(double) );

	MPI_Barrier(MPI_COMM_WORLD);
	printf("Check 8\n");
	for (i=0;i<modes_size;i++) {
		modes_send[i] = modes[i];
		// printf("[%d] modes send:\t%d\t%e\t%e\n",rank,i,modes_send[i],modes[i]);
	}

	MPI_Reduce(modes_send,modes_recv,modes_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		if(mode==0){
			model = 999;
		}else if(mode==1){
			model = 998;
		}else{
			model = mode;
			// model = 6;
		}
		
		double *modes_tmp = (double *)create_vector(ortho_size);
		for(i=0;i<modes_size;i++){
			modes_tmp[i] = 0.0;
			for(j=0;j<modes_size;j++){
				modes_tmp[i] += orthoinv[i][j]*modes_recv[j];
			}
			printf("[%d] modes revc:\t%d\t%e\t%e\n",rank,i,modes_tmp[i],modes_recv[i]);
		}
		for(i=0;i<modes_size;i++){
			results_m[0] = (double)i;
			results_m[1] = modes_tmp[i];
			update_modes(results_m);
			printf("[%d]\t%d\t%e\n", rank, i, results_m[1]);
		}
		
		output_modes(0,0,0);
		
		printf("Modes R\n");

		create_modesR();
		for(i=0;i<modes_size;i++){
			x=0.0;
			for(j=0;j<modes_size;j++){
				x += lambdainv[j][i]*modes_tmp[j];
				// if(i>1995)printf("ortho: %d\t%d\t%e\n",j,i,get_orthol_TTT(j,i));
			}
			results_m[0] = (double)i;
			results_m[1] = x;
			update_modesR(results_m);
			printf("[%d]\t%d\t%e\n", rank, i, results_m[1]);
		}
	
	}
	printf("Check FINAL\n");
	
	MPI_Finalize();
}

