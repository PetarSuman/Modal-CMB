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


	if (argc != 3) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s directory file\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	int i,j,k,n,m;
	
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");
	
	char filename[100] = "";
	strcat(filename, directory);
	strcat(filename, argv[2]);
	
// ini file
	char inifile[MAXLEN];
	double result;

// mpi vars
	int rank, nproc;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	printf("[%d] Process[%d]  \n", rank, getpid());
	
// load alpha

	double* values = malloc( sizeof(int)*4*MAXLINES);
	int* size = malloc( sizeof(int));

	load_four(filename, values, size);
	
	printf("%s\t%d\n",filename,*size);

	create_alpha();
	
	double results[4];
	
	j=0;
	for(i=0;i<*size;i++){
		results[0] = values[j++];
		results[1] = values[j++];
		results[2] = values[j++];
		results[3] = values[j++];
		if(results[3]!=0)printf("%d\t%d\t%d\t%e\n",(int)results[0],(int)results[1],(int)results[2],results[3]);
		update_alpha(results);
	}

	output_alpha(directory);
	
	printf("[%d] finished\n", rank);

	MPI_Finalize();
	
	// End of code
	return 0;
}