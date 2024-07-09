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
		printf( "Usage is >:%s file\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char file[100] = "";
	strcat(file, argv[1]);
	
	int i,j,k,l,m,n;

// mpi vars
	int rank, nproc, lnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

// 	printf("[%d] Process[%d]  \n", rank, getpid());
	
	double* data_raw = malloc( sizeof(double)*MAXLINES*6);
	int* size = malloc( sizeof(int));
	
	load_six(file, data_raw, size);
	
	int length = *size;
	
	double** data = create_array(6,length);
	
	j = 0;
	for(i=0;i<length;i++){
		data[0][i] = data_raw[j++];
		data[1][i] = data_raw[j++];
		data[2][i] = data_raw[j++];
		data[3][i] = data_raw[j++];
		data[4][i] = data_raw[j++];
		data[5][i] = data_raw[j++];
		if(data[3][i] != 0.0 && data[0][i]<=2000.0 && data[1][i]<=2000.0 && data[2][i]<=2000.0 && data[0][i]+data[1][i]>=data[2][i] && data[1][i]+data[2][i]>=data[0][i] && data[2][i]+data[0][i]>=data[1][i])printf("%d\t%d\t%d\t%e\n", (int)data[0][i], (int)data[1][i], (int)data[2][i], data[3][i]);
	}
	
	int lsize =1;
	for(i=1;data[2][i]>data[2][i-1] && data[2][i]<=2000.0;i++) lsize++;
	int lvalues[lsize];
	for(i=0;i<lsize;i++) lvalues[i] = (int)data[2][i];
	
	
	create_bispectrum(lsize);
	
	double results[4];
	for(n=0;n<length;n++){
		if(data[3][n] != 0.0 && data[0][n]<=2000.0 && data[1][n]<=2000.0 && data[2][n]<=2000.0 && data[0][n]+data[1][n]>=data[2][n] && data[1][n]+data[2][n]>=data[0][n] && data[2][n]+data[0][n]>=data[1][n]) {
			results[0] = round((data[0][n] - 9.0) / 19.335);
			results[1] = round((data[1][n] - 9.0) / 19.335);
			results[2] = round((data[2][n] - 9.0) / 19.335);
			results[3] = data[3][n];
		}
	}
	
// 	for(i=0;i<length;i++){
// 		if (data[0][i] == data[1][i] && data[1][i] == data[2][i] && data[3][i] != 0.0){
// 			printf("%d\t%d\t%d\t%e\t%e\n", (int)data[0][i], (int)data[1][i], (int)data[2][i], data[3][i], data[4][i]);
// 		}
// 	}
// 	int test;
// 	for(j=1;j<lsize;j++){
// 		printf("\n");
// 		for(i=0;i<length;i++){
// 			test = (data[0][i] == lvalues[j] && data[1][i] == data[2][i]);
// 			test = test || (data[1][i] == lvalues[j] && data[2][i] == data[0][i]);
// 			test = test || (data[2][i] == lvalues[j] && data[0][i] == data[1][i]);
// 			test = test && data[0][i]+data[1][i]>=data[2][i];
// 			if (test){
// 				printf("%d\t%d\t%d\t%e\t%e\n", (int)data[0][i], (int)data[1][i], (int)data[2][i], data[3][i], data[4][i]);
// 			}
// 		}
// 	}
	
	MPI_Finalize();
	
	// End of code
	return 0;
}