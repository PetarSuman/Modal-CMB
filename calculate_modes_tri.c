#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "global.h"

double calculate_modes_tri(int r, int size, int* vec, double**** trispectrum){
	
	double perm;

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	int i,j,k,l;
	double result;
	double c[size];
	
	for (i=0; i<size; i++){
		c[i] = 0;
		if (get_cl(i)!=0 && vec[i]!=0) c[i] = get_beam(i) / (pow(2.0*vec[i]+1.0, 1.0/6.0) * sqrt(get_beam(i)*get_beam(i)*get_cl(i)+get_noise(i)));
	}
	
	result = 0.0;
	for(i=0;i<size-1;i++){
		for(j=i;j<size-1;j++){
			for(k=j;k<size-1;k++){
				for(l=k;l<size-1;l++){
					
					if(l<=i+j+k+2){
						if(i==j){
							if(j==k){
								if(k==l){
									perm=1.0*c[i]*c[j]*c[k]*c[l];
								}else{
									perm=4.0*c[i]*c[j]*c[k]*c[l];
								}
							}else{
								if(k==l){
									perm=6.0*c[i]*c[j]*c[k]*c[l];
								}else{
									perm=12.0*c[i]*c[j]*c[k]*c[l];
								}
							}
						}else{
							if(j==k){
								if(k==l){
									perm=4.0*c[i]*c[j]*c[k]*c[l];
								}else{
									perm=12.0*c[i]*c[j]*c[k]*c[l];
								}
							}else{
								if(k==l){
									perm=12.0*c[i]*c[j]*c[k]*c[l];
								}else{
									perm=24.0*c[i]*c[j]*c[k]*c[l];
								}
							}
						}
						
						if(vec[i]+vec[j]+vec[k]>=vec[l]&&vec[i]>1){
							result += perm * calculate_geometric_tri(vec[i], vec[j], vec[k], vec[l]) *trilR(r,i,j,k,l)*trispectrum[i][j][k][l];
						}
					}
				}	
			}
		}
// 		printf("result %d\t%e\n",r,result);
		// MPI sync
		if ( myrank == 0 ) sync_tasks(15,3);
	}
	
	return result;
}
