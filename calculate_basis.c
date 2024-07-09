#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "global.h"

double calculate_baseq_tri(int i,int j, double max){
 
 	double integral=0;
 	
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	double k1,k2,k3,k2,di;

	int r,q1,q2,q3,q4;

	for(r=0;r<alpha_points*2;r++){
		for(q1=0;q1<alpha_points;q1++){
			for(q2=0;q2<alpha_points;q2++){
				for(q3=0;q3<alpha_points;q3++){
					for(q4=0;q4<alpha_points;q4++){
						if(tri(r,q1,q2) && tri(r,q3,q4)){
							
						}
					}
				}	
			}
		}
		// MPI sync
		if ( myrank == 0 ) sync_tasks(5,2);
	}
	return integral;
}
double calculate_baser_tri(int i,int j, double max){
 
 	double integral=0;
 	double cube_size = pow((double)alpha_points,-3);
 	
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	int p1,p2,p3;
	
// 	printf("model_eigen: %d\n",model);
	
	double k1,k2,k3;

	int r,q1,q2,q3,q4;

	for(r=0;r<alpha_points*2;r++){
		for(q1=0;q1<alpha_points;q1++){
			for(q2=0;q2<alpha_points;q2++){
				for(q3=0;q3<alpha_points;q3++){
					for(q4=0;q4<alpha_points;q4++){
						if(tri(r,q1,q2) && tri(r,q3,q4)){
							
						}
					}
				}	
			}
		}
		// MPI sync
		if ( myrank == 0 ) sync_tasks(5,2);
	}
	return integral;
}