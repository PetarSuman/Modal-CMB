#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include "global.h"

/*	
	The idea is to project primordial (k-space) basis to late time (l-space).

	This function (split_mpi.c) projects the primordial basis polynomials q(k) onto 
	their late-time value qtilde (i,l,x). 
	The projection is given by theoretical formula:

		qtilde_i(l,x) = \int \Delta(l,k) * j_l(kx) * q_i(k)

	where:	Delta(l,k) = transfer function
			j_l(kx) = Bessel function
			q_i(k) = primordial basis polynomial

	This is the crucial step needed to later obtain the actual late time basis via the
	late-time-projection of the primordial basis
*/

int main( int argc, char *argv[] ){

	double pi = 3.141592653589793;

	double time1, time2, time3, time4, time5, duration;

	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
	}

	//	Initialise parameters and read from file
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);

	load_bessel();
	load_transfer();

	int* lvalues_T = malloc( sizeof(int)*MAXLINES);
	int* lvalues_E = malloc( sizeof(int)*MAXLINES);
	int* lceil = malloc( sizeof(int));
	int lcol = 0;
	int restart = 0;
	int lsize_T;
	int lsize_E;
	int ltotal;
	int p,l,i,j,k,n;
	double area, result;

	//	MPI Variables and Initialization
	int rank, nproc, lnext, bnext;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);


	//	Load l's from ini file or cli 
	if(use_l_file==1){
		if ( load_one(inifile, lvalues_T, lceil ) ) {			
			printf("[%d] Ini load complete\n", rank);
		}
		else{
			printf("[%d] Ini file error\n", rank);
			MPI_Finalize();
			exit(1);
		}
		lsize_T = *lceil;
		if(do_polarisation==1){
			lsize_E = lsize_T;
			for (i=0;i<lsize_T;i++){
				lvalues_E[i] = lvalues_T[i];
			}
		}
	}
	else{
		lsize_T = eflag_T_lmax+1;
		for(i=0;i<lsize_T;i++){
			lvalues_T[i] = i;
		}
		if(do_polarisation==1){
			lsize_E = eflag_E_lmax+1;
			for(i=0;i<lsize_E;i++){
				lvalues_E[i] = i;
			}
		}
	}

	//	Set up primordial basis and late-time ordering
	set_terms_prim();
	set_terms_late();
	double kmax = get_kmax();
	int ksize = get_ksize();
	double kvalues[ksize];
	get_kvec(kvalues);
	create_basis_prim(ksize, kmax, kvalues);
	int pmax = get_pmax_prim();
	double prim;
	double kval, orig;

	//	Informative print statements	
	if(rank==0){						
		printf("kmax_cut: %e\n", kmax);
		printf("Primordial basis ordering:\nn\ti\tj\tk\n");
		for(n=0;n<get_terms_prim();n++){
			find_perm_prim(n,&i,&j,&k);
			printf("%d\t%d\t%d\t%d\n",n,i,j,k);
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	//	Calculate k integrals
	double xstep;
	double step3 = 200.0;
	double step2 = 50.0;
	double step1 = 10.0;
	double step0 = 5.0;
	
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

	int xsize=n;			//	i.e. xsize = xmax
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
	
	int elements_T = lsize_T * xsize;
	int element_list_T[elements_T][2];
	double **array_T = (double **)create_array(lsize_T,xsize);
	n=0;
	for (i=0;i<lsize_T;i++){
		for (j=0;j<xsize;j++){
			element_list_T[n][0] = i;
			element_list_T[n][1] = j;
			n++;
			array_T[i][j] = 0e0;
		}
	}
		
	int loops;
	int auxloop ;
	int start_loop_T;
	int end_loop_T;
	
	loops=elements_T/nproc;
	auxloop = fmod(elements_T,nproc);
	start_loop_T = rank*loops;
	end_loop_T = (rank+1)*loops-1;
	
	if (auxloop != 0){
		if (rank < auxloop){
			start_loop_T = start_loop_T + rank;
			end_loop_T = end_loop_T + rank + 1;
		}else{
			start_loop_T = start_loop_T + auxloop;
			end_loop_T = end_loop_T + auxloop;
		}
	}

	double* array_T_send = (double *)malloc( elements_T * sizeof(double) );
	double* array_T_recv = (double *)malloc( elements_T * sizeof(double) );
	
	create_qtilde(lsize_T,lsize_E,xsize);
	
	time1 = csecond();
	for(p=0;p<pmax+1;p++){
		//for (i=0;i<lsize_T;i++){
		//	for (j=0;j<xsize;j++){
		//		array_T[i][j]=0e0;
		//	}
			//		This is already computed before...
		//}
		for(n=start_loop_T;n<end_loop_T+1;n++){
			i = element_list_T[n][0];
			j = element_list_T[n][1];
			l = lvalues_T[i];
			x = xvec[j];
			if(l>=2){
				array_T[i][j] = M_2_PI * calculate_kint_T(p, l, x);
			}else{
				array_T[i][j] = 0e0;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		n=0;	
		for (i=0;i<lsize_T;i++){
			for (j=0;j<xsize;j++){
				array_T_send[n] = array_T[i][j];
				n++;
			}
		}
		MPI_Reduce(array_T_send,array_T_recv,elements_T,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(rank==0){
			n=0;	
			for (i=0;i<lsize_T;i++){
				for (j=0;j<xsize;j++){
					x = array_T_recv[n];
					update_qtilde_T(p,i,j,x);
					printf("qtilde[p = %d][l = %d][x = %d] = %e\n", p,i,j,x);
					n++;
				}
			}
			duration = csecond() - time1;
			printf("[%d] Done temperature projectionin rank %d\t time = %es\n",rank,p,duration);
			time1 = csecond();
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	// this part of the code is attempt to split code into T-only and E-only part.
	// otherwise we risk running out of memory

	//create_qtilde(lsize_T,lsize_E,xsize);

	// end of the modified code

	
	if(do_polarisation==1){						//	Repeat same procedure for E-mode
		printf("[%d] Polarization included\n", rank);
		int elements_E = lsize_E * xsize;
		int element_list_E[elements_E][2];
		double **array_E = (double **)create_array(lsize_E,xsize);
		n=0;
		if (rank == 0) printf("Check 1\n");
		for (i=0;i<lsize_E;i++){
			for (j=0;j<xsize;j++){
				element_list_E[n][0] = i;
				element_list_E[n][1] = j;
				n++;
				array_E[i][j] = 0e0;
			}
		}
		if (rank == 0) printf("Check 2\n");
		int loops;
		int auxloop ;
		int start_loop_E;
		int end_loop_E;
		loops=elements_E/nproc;
		auxloop = fmod(elements_E,nproc);
		start_loop_E = rank*loops;
		end_loop_E = (rank+1)*loops-1;
		if (rank == 0) printf("Check 3\n");
		if (auxloop != 0){
			if (rank < auxloop){
				start_loop_E = start_loop_E + rank;
				end_loop_E = end_loop_E + rank + 1;
			}else{
				start_loop_E = start_loop_E + auxloop;
				end_loop_E = end_loop_E + auxloop;
			}
		}
		if (rank == 0) printf("Check 4\n");
		double* array_E_send = (double *)malloc( elements_E * sizeof(double) );
		double* array_E_recv = (double *)malloc( elements_E * sizeof(double) );
		
		time1 = csecond();
		for(p=0;p<pmax+1;p++){
			for (i=0;i<lsize_E;i++){
				for (j=0;j<xsize;j++){
					array_E[i][j]=0e0;
				}
			}
			if (rank == 0) printf("Check 4a\n");
			for(n=start_loop_E;n<end_loop_E+1;n++){
				i = element_list_E[n][0];
				j = element_list_E[n][1];
				l = lvalues_E[i];
				x = xvec[j];
				if(l>=2){
					if (rank == 0) printf("Check 4b\n");
					array_E[i][j] = M_2_PI * calculate_kint_E(p, l, x);
				}else{
					array_E[i][j] = 0e0;
				}
			}
			if (rank == 0) printf("Check 4c\n");
			MPI_Barrier(MPI_COMM_WORLD);
			
			n=0;	
			for (i=0;i<lsize_E;i++){
				for (j=0;j<xsize;j++){
					array_E_send[n] = array_E[i][j];
					n++;
				}
			}
			if (rank == 0) printf("Check 5\n");
			MPI_Reduce(array_E_send,array_E_recv,elements_E,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			if(rank==0){
				n=0;
				if (rank == 0) printf("Check 6\n");	
				for (i=0;i<lsize_E;i++){
					for (j=0;j<xsize;j++){
						x = array_E_recv[n];
						update_qtilde_E(p,i,j,x);
						n++;
					}
				}
				duration = csecond() - time1;
				printf("[%d] Done polarisation projection in rank %d\t time = %es\n",rank,p,duration);
				time1 = csecond();
			}	
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank == 0){
		output_qtilde(xsize,xvec,lsize_T,lvalues_T,lsize_E,lvalues_E);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	/*
	time1 = csecond();
	printf("Read qtilde:------------------------------------------\n");
	read_qtilde();
	time2 = csecond();
	printf("[%d] Read qtilde in %e s\n", rank, time2-time1);
	*/

	//...... Original part of code (incompatible with Alex's part)

	/*double lmax_T;
	double *Tvec;
	printf("[%d] Want to create vector(lsize_T)\n", rank);
	Tvec = create_vector(lsize_T);
	printf("[%d] Vector T created.\n", rank);
	for(i=0;i<lsize_T;i++)Tvec[i] = (double)lvalues_T[i];
	lmax_T = Tvec[lsize_T-1];

	double lmax_E;
	double *Evec;
	//	lsize_E = 1;	//just so that I don't have malloc(0)
	printf("[%d] Want to create vector(lsize_E)\n", rank);
	Evec = create_vector(lsize_E);
	printf("[%d] Vector E created.\n", rank);
	for(i=0;i<lsize_E;i++)Evec[i] = (double)lvalues_E[i];
	lmax_E = Evec[lsize_E-1];
	printf("[%d] lsize_E = %d\n", rank, lsize_E);

	printf("[%d] Want to create basis late...\n", rank);
	create_basis_late(lsize_T, lsize_E, 2, lmax_T, 2, lmax_E, Tvec, Evec);
	printf("[%d] Successfully created basis late...\n", rank);*/

/*
	//...........this is Alex's version of the code 
	printf("[%d] Start of Alex's part......\n", rank);
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

	printf("[%d] Want to create basis late...... and lsize_E = %d\n", rank, Esize);
    create_basis_late(Tsize, Esize, Tmin, Tmax, Emin, Emax, Tvec, Evec);

	// ...................    end of ALex's version

	int start_loop;
	int end_loop;
	double x1,x2,x3,x4;
	model = 62;
	//read_eigen(1,1,1);
	int size = get_terms_prim();

	n=0;
	for (i=0;i<eflag_T_lmax+1;i+=10){
		for (j=i;j<eflag_T_lmax+1;j+=10){
			for (k=j;k<eflag_T_lmax+1;k+=10){
				if(i+j>=k)n++;
			}
		}
	}
	ltotal = n+1;
	printf("[%d] ltotal: %d\n",rank, ltotal);

	int **ltriplesXXX = create_iarray(ltotal,3);

	printf("[%d] tasks TTT EEE %d\n",rank,ltotal);

	n=0;

	for (i=0;i<eflag_T_lmax+1;i+=10){
		for (j=i;j<eflag_T_lmax+1;j+=10){
			for (k=j;k<eflag_T_lmax+1;k+=10){
				if(i+j>=k){
					ltriplesXXX[n][0] = i;
					ltriplesXXX[n][1] = j;
					ltriplesXXX[n][2] = k;
					n++;
				}
			}
		}
	}
	
	loops=(ltotal-1)/nproc;
	auxloop = fmod(ltotal,nproc);
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
	MPI_Barrier(MPI_COMM_WORLD);

	printf("l1\tl2\tl3\tQ_tilde\tS(l1,l2,l3)\n");

	for(n=start_loop;n<end_loop+1;n++){
		time1 = csecond();

		l1 = ltriplesXXX[n][0];
		l2 = ltriplesXXX[n][1];
		l3 = ltriplesXXX[n][2];

		if(l1==0)l1=2;
		if(l2==0)l2=2;		//change this to T_lmin
		if(l3==0)l3=2;
		
		x3 = 0e0;
		x4 = 0e0;

		for (i=0;i<size;i++){
			x1 = calculate_xint_TTT(l1, l2, l3, i, xsize, xvec);
			x3 += get_eigen(i) * x1;
		}

		duration = csecond() - time1;
		//if (rank == 0) printf("%d\t%d\t%d\t%e\t%e\n",l1,l2,l3,x1,x3);

	} 
	MPI_Barrier(MPI_COMM_WORLD);*/
	
	printf("[%d] Finished code\n", rank);

	MPI_Finalize();
	
	return 0;
}