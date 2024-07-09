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

	int* lvalues_T = malloc( sizeof(int)*MAXLINES);
	int* lvalues_E = malloc( sizeof(int)*MAXLINES);
	int* lceil = malloc( sizeof(int));
	int  lcol = 0;
	int restart = 0;
	int lsize_T;
	int lsize_E = 0;
	int ltotal;
	int p,l,i,j,k,n;
	double area, result;

    int rank, nproc, lnext, bnext;

    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	if(use_l_file==1){
		if ( load_one(inifile, lvalues_T, lceil ) ) {			
			printf("[%d] Ini load complete\n", rank);
		} else {
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

 	double lmax_T;
	double *Tvec;
	printf("[%d] Want to create vector(lsize_T)\n", rank);
	Tvec = create_vector(lsize_T);
	printf("[%d] Vector T created.\n", rank);
	for(i=0;i<lsize_T;i++)Tvec[i] = (double)lvalues_T[i];
	lmax_T = Tvec[lsize_T-1];


	double lmax_E;
	double *Evec;
	lsize_E = 1;	//just so that I don't have malloc(0)
	printf("[%d] Want to create vector(lsize_E)\n", rank);
	Evec = create_vector(lsize_E);
	printf("[%d] Vector E created.\n", rank);
	for(i=0;i<lsize_E;i++)Evec[i] = (double)lvalues_E[i];
	lmax_E = Evec[lsize_E-1];

	//create_basis_late(lsize_T, lsize_E, 2, lmax_T, 2, lmax_E, Tvec, Evec);
    printf("%d = \n", get_cl_TT(2));

    
	MPI_Finalize();    
    return 0;


}