#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include "global.h"

static int glbcounter = 0; // next l value in the list
static int maxcounter;   // total number of tasks 
static int npdone = 0; // processes done 
/*
1 = alpha
2 = beta
3 = bispectrum
4 = decompose
5 = eigen
6 = ortho
7 = modes
8 = orthol
9 = modesR
*/

int update(int mode, int size, double *results){

  switch(mode){
  	case 1:
  		update_alpha(results);
  	break;
  	case 2:
//   		update_beta(size, results);
  	break;
  	case 3:
  		update_bispectrum(results);
  	break;
  	case 4:
  		update_decompose(results);
  	break;
  	case 5:
  		update_eigen(results);
  	break;
  	case 6:
  		update_ortho(results);
  	break;
  	case 7:
  		update_modesR(results);
  	break;
  	case 8:
  		update_orthol(results);
  	break;
  	case 9:
  		update_modes(results);
  	break;
  	case 10:
  		update_ortho_tri(results);
  	break;
  	case 11:
  		update_eigen_tri(results);
  	break;
  	case 12:
  		update_beta_tri(size, results);
  	break;
  	case 13:
  		update_trispectrum(results);
  	break;
  	case 14:
  		update_orthol_tri(results);
  	break;
  	case 15:
  		update_modesR_tri(results);
  	break;
  	case 16:
  		update_modes_tri(results);
  	break;
  	case 17:
  		update_transfer(size, results);
  	break;
  	case 18:
  		update_fisher_tri(results);
  	break;
  	case 19:
  		update_gamma_tri(results);
  	break;
  	case 20:
  		update_gamma(results);
  	break;
  	case 21:
//   		update_beta2(size,results);
  	break;
  	default:
  		printf("invalid mode specified for mpi comm\n");
  	break;
  	
	free(results);
  }
  
  return 0;
}

// master poll
int sync_tasks(int mode, int size){
	int flag;
	
	MPI_Request s;
	MPI_Status status;
	
	MPI_Iprobe( MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
	
	if ( flag ) {
// 		printf("ST: found message\n");
		double *results =  malloc( sizeof(double)*size);
		MPI_Irecv(results, size, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &s);
		MPI_Wait(&s, &status);
// 		MPI_Request_free(&s);
// 		printf("ST: recieved message\n");
		update(mode, size, results);
// 		printf("ST: updated message\n");
		
		MPI_Isend(&glbcounter, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &s);
// 		printf("ST: sent reply\n");
		MPI_Wait(&s, &status);
// 		printf("ST: reply recieved\n");
// 		MPI_Request_free(&s);
		if ( glbcounter >= maxcounter ) npdone++;
		glbcounter++;
		
		return 1;
	}
	
	return 0;
}


// worker polls & sends results
int get_next_task(int mode, int size, int rank, double *results){
	int lc, ln;
	int one, i, np;
	MPI_Request s;
	MPI_Status status;
	double t0, t1;
	int flag;
// MPI sync
	if ( rank == 0 ) {
		if ( glbcounter == 0 ) {
// 			printf("GNT: Starting\n");
			sleep(1);
			MPI_Comm_size(MPI_COMM_WORLD, &np);
			for ( i=0; i<np ; i++ ) {
// 			for ( i=0; i<np/2 ; i++ ) {
				sync_tasks(mode, size);
// 				printf("GNT: sent start value\n");
			}
		}
		update(mode, size, results);
		ln = glbcounter;
		glbcounter++;
		while (sync_tasks(mode, size)) {}
		
	} else {
		MPI_Iprobe( 0, 0, MPI_COMM_WORLD, &flag, &status);
		
		if ( !flag ){
			MPI_Isend(results, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &s);
			MPI_Wait(&s, &status);
// 			MPI_Request_free(&s);
		}
		
		MPI_Recv(&ln, 1, MPI_INT, 0,0, MPI_COMM_WORLD, &status);
		MPI_Wait(&s, &status);
// 		MPI_Request_free(&s);
	}
	
	return ln;
}

// master to workers
int signal_end_tasks(int mode, int size){
  int np,i;
  MPI_Request s;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  
  npdone++; // master also done
  while ( npdone < np ) sync_tasks(mode, size);
//   while ( npdone < np/2 ) sync_tasks(mode, size);
  printf("[0] All done. Exiting. \n");
  return 1;
}

// master record tasks
int record_tasks(int ntasks){
	glbcounter = 0;
	npdone = 0;
  maxcounter = ntasks;
  printf("[0] tasks = %d \n", ntasks);
  return 0;
}