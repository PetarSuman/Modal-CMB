#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "global.h"

void modify_bispectrum(int size, int* vec, double*** bispectrum){

	int i,j,k;
	double c[size];
	double time, duration;
	int test;
	double x1;
	
	for (i=0; i<size; i++){
		c[i] = 0;
// 		x1 = calculate_ISW(size,vec,2,i,i);
// 		x1 = get_beam(i)*get_beam(i)*get_cl(i)+get_noise(i);
// 		printf("%d\t%e\n",i,x1);
		if (get_cl_TT(i)!=0 && vec[i]!=0) c[i] = get_beam_TT(i) / (pow(2.0*vec[i]+1.0, 1.0/6.0) * sqrt(get_beam_TT(i)*get_beam_TT(i)*get_cl_TT(i)+get_noise_TT(i)));
	}
	
	for(i=0;i<size;i++){
		l1 = vec[i];
// 		time = csecond();
		for(j=0;j<size;j++){
			l2 = vec[j];
			for(k=0;k<size;k++){
				l3 = vec[k];
				test = (l1+l2+l3) - 2*((l1+l2+l3)/2);
				if(!test && l1+l2>=l3 && l2+l3>=l1 && l3+l1>=l2){
// 					x1 = calculate_ISW(size,vec,i,j,k);
					x1 = bispectrum[i][j][k]*1.55e-8*1.55e-8;
					bispectrum[i][j][k] = c[i]*c[j]*c[k] * calculate_geometric(l1, l2, l3) * x1;
// 					if((i==2||fmod(i,10)==0)&&(j==2||fmod(j,10)==0)&&(k==2||fmod(k,10)==0)&&i<=j&&j<=k){
// 						printf("%d\t%d\t%d\t%e\n",l1,l2,l3,x1);
// 					}
				} else {
					bispectrum[i][j][k] = 0.0;
				}
			}
		}
// 		duration = csecond() - time;
// 		printf("done: %d\t%e\n",i,duration);
	}
}

double calculate_modes(int r, int size, double*** bispectrum){

	double factor;
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int i,j,k;
	double result;
	double time, duration;
	
	result = 0.0;
	for(i=0;i<size;i++){
// 		time = csecond();
		for(j=i;j<size;j++){
			for(k=j;k<size;k++){
				if(bispectrum[i][j][k]!=0.0){
					if(i==j&&j==k){
						factor = 1.0;
					} else if (i==j || j==k){
						factor = 3.0;
					} else {
						factor = 6.0;
					}
					result += factor * plijk_TTT(r,i,j,k) * bispectrum[i][j][k];
				}
			}
		}
// 		duration = csecond() - time;
// 		printf("done: %d\t%e\n",i,duration);
		// MPI sync
		if ( myrank == 0 ) sync_tasks(9,2);
	}
	return result;
}

// double calculate_modes(int r, int size, int* vec, double*** bispectrum){
// 
// 	double factor;
// 	int myrank;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
// 	int i,j,k;
// 	double result;
// 	double c[size];
// 	double time, duration;
// 	
// 	for (i=0; i<size; i++){
// 		c[i] = 0;
// 		if (get_cl(i)!=0 && vec[i]!=0) c[i] = get_beam(i) / (pow(2.0*vec[i]+1.0, 1.0/6.0) * sqrt(get_beam(i)*get_beam(i)*get_cl(i)+get_noise(i)));
// 	}
// 	result = 0.0;
// 	for(i=0;i<size;i++){
// 		time = csecond();
// 		for(j=i;j<size;j++){
// 			for(k=j;k<size;k++){
// 				if(i==j&&j==k){
// 					factor = 1.0*c[i]*c[j]*c[k];
// 				} else if (i==j || j==k){
// 					factor = 2.0*c[i]*c[j]*c[k];
// 				} else {
// 					factor = 6.0*c[i]*c[j]*c[k];
// 				}
// 				if(vec[i]+vec[j]>=vec[k]){
// 					result += factor * calculate_geometric(vec[i], vec[j], vec[k]) * pijk(r,i,j,k) * bispectrum[i][j][k];
// 				}
// 			}
// 		}
// 		duration = csecond() - time;
// 		printf("done: %d\t%e\n",i,duration);
// 		// MPI sync
// 		if ( myrank == 0 ) sync_tasks(7,2);
// 	}
// 	return result;
// }

double*** reconstruct_b(int size, int* vec){

	int i,j,k;
	double value;
	double factor;
	
	long int l_size_long = size;
	double ***bispectrum = create_3Darray_long(l_size_long,l_size_long,l_size_long);

	double c[size];
	
	for (i=0; i<size; i++){
		c[i] = 0;
		if (get_cl_TT(i)!=0 && vec[i]!=0) c[i] = sqrt(get_beam_TT(i)*get_beam_TT(i)*get_cl_TT(i)+get_noise_TT(i)) / (get_beam_TT(i) * pow(2.0*vec[i]+1.0, 1.0/6.0));
// 		if (get_cl(i)!=0 && vec[i]!=0) c[i] = sqrt(get_cl(i)) / pow(2.0*vec[i]+1.0, 1.0/6.0);
	}
	
	for(i=0;i<size;i++){
		for(j=i;j<size;j++){
			for(k=j;k<size;k++){
				if(k<i+j+1){
					
					factor = c[i]*c[j]*c[k];
					
					value = factor*bijk_TTT(i,j,k);
					bispectrum[i][j][k] = value;
					bispectrum[j][k][i] = value;
					bispectrum[k][i][j] = value;
					bispectrum[k][j][i] = value;
					bispectrum[j][i][k] = value;
					bispectrum[i][k][j] = value;
				} else {
					bispectrum[i][j][k] = 0;
					bispectrum[j][k][i] = 0;
					bispectrum[k][i][j] = 0;
					bispectrum[k][j][i] = 0;
					bispectrum[j][i][k] = 0;
					bispectrum[i][k][j] = 0;
				}
			}
		}
	}
	
	return bispectrum;
}

double*** reconstruct_bn(int max, int size, int* vec){

	double value;
	double factor;
	
	long int l_size_long = size;
	double ***bispectrum = create_3Darray_long(l_size_long,l_size_long,l_size_long);

	double c[size];
	
	
	int terms=get_terms_late();
	double* modes = (double *)create_vector(max+1);
	double* modesR = (double *)create_vector(max+1);
	double** ortho = (double **)create_array(max+1,max+1);
	
	double k1,k2,k3,ksum;

	int i,j,k,l,m,n,p;
	
	for(n=0;n<max+1;n++){
		modesR[n] = get_modesR_TTT(n);
	}
	for(i=0;i<max+1;i++){
		for(j=0;j<max+1;j++){
			ortho[i][j] = get_orthol_TTT(i,j);
		}
	}
	
	for(i=0;i<max+1;i++){
		modes[i] = 0.0;
		for( j=i;j<max+1;j++){
			modes[i] += ortho[j][i]*modesR[j];
		}
	}

// 	for(n=0;n<max+1;n++){
// 		modes[n] = get_modes(n);
// 	}
	
	for (i=0; i<size; i++){
		c[i] = 0;
		if (get_cl_TT(i)!=0 && vec[i]!=0) c[i] = sqrt(get_beam_TT(i)*get_beam_TT(i)*get_cl_TT(i)+get_noise_TT(i)) / (get_beam_TT(i) * pow(2.0*vec[i]+1.0, 1.0/6.0));
	}
	
	for(i=0;i<size;i++){
		for(j=i;j<size;j++){
			for(k=j;k<size;k++){
				if(k<i+j+1){
// 					printf("here %d %d %d\n",i,j,k);
					factor = c[i]*c[j]*c[k];
					
					for(n=0;n<max+1;n++){
						value += modes[n]*plijk_TTT(n,i,j,k);
					}
					
					value *= c[i]*c[j]*c[k];
					
					bispectrum[i][j][k] = value;
					bispectrum[j][k][i] = value;
					bispectrum[k][i][j] = value;
					bispectrum[k][j][i] = value;
					bispectrum[j][i][k] = value;
					bispectrum[i][k][j] = value;
				} else {
					bispectrum[i][j][k] = 0;
					bispectrum[j][k][i] = 0;
					bispectrum[k][i][j] = 0;
					bispectrum[k][j][i] = 0;
					bispectrum[j][i][k] = 0;
					bispectrum[i][k][j] = 0;
				}
			}
		}
	}
	
	printf("here %d\n", max);
	return bispectrum;
}
