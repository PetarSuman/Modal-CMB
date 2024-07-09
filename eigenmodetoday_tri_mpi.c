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
	
// **1**

	double time1, time2, time3, time4, time5, duration;


	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s data\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}

// mpi vars
	int rank, nproc, onext, mnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");

	char l_number_file[100] = "";
	strcat(l_number_file, directory);
	strcat(l_number_file, l_size_file);

	char l_file[100] = "";
	strcat(l_file, directory);
	strcat(l_file, l_data_file);
	
	char tri_file[100] = "";
	strcat(tri_file, directory);
	strcat(tri_file, trispectrum_file);
	
	int number = 1;
	
	int *l_data_size = create_ivector(number);
	ivector_read(&number, l_number_file, &l_data_size[0]);

	int l_size = l_data_size[0];
	
	int *l_values = create_ivector(l_size);
	ivector_read(&l_size, l_file, &l_values[0]);
	
	double lmax = (double)l_values[l_size-1];
	
	set_terms_tri();
	create_cl(l_size);
	create_beam(l_size);
	create_noise(l_size);
	load_cl(directory, l_size, l_values);
	load_BN(directory, l_size, l_values);
	int i,j,k,l,m,n;
		
	double result=0;
	
	int modes_size = get_terms_tri();
	
	if(rank==0){
		printf("lmax: %d\n", (int)lmax);
		printf("pmax: %d\n", get_qmax());
		printf("terms: %d\n", get_terms_tri());
		for(n=0;n<modes_size;n++){
			find_perm_tri(n,&i,&j,&k,&l);
			printf("%d\t%d\t%d\t%d\t%d\n",n,i,j,k,l);
		}
	}
	
	double xvec[l_size];
	for(i=0;i<l_size;i++){
		xvec[i] = (double)l_values[i] / lmax;
	}

	create_basis_tri(l_size,xvec);
	read_orthol_tri(directory);
	read_gamma_tri(directory);
	read_eigen_tri(directory);

	create_modes_tri();

	double *results_m =  malloc( sizeof(double)*2);
	results_m[0] = 0.0;
	results_m[1] = 0.0;
	double *a1 = create_vector(modes_size);
	double *a2 = create_vector(modes_size);
	double *a3 = create_vector(modes_size);
	
	for(i=0;i<modes_size;i++){
		a1[i] = 0.0;
		for(j=0;j<modes_size;j++){
			a1[i] += get_gamma_tri(j,i)*get_eigen_tri(j);
			printf("%d\t%d\t%e\n", i, j, get_gamma_tri(i,j));
		}
		printf("[1]\t%d\t%e\n", i, a1[i]);
	}
	
	for(i=0;i<modes_size;i++){
		a2[i] = 0.0;
		for(j=0;j<modes_size;j++){
			a2[i] += get_orthol_tri(i,j)*a1[j];
		}
		printf("[2]\t%d\t%e\n", i, a2[i]);
	}
	
	for(i=0;i<modes_size;i++){
		a3[i] = 0.0;
		for(j=0;j<modes_size;j++){
			a3[i] += get_orthol_tri(j,i)*a2[j];
		}
		printf("[3]\t%d\t%e\n", i, a3[i]);
		results_m[0] = (double)i;
		results_m[1] = a3[i];
		update_modes_tri(results_m);
	}
	
// 	for(i=0;i<modes_size;i++){
// 		sum1 = 0.0;
// 		for(j=0;j<modes_size;j++){
// 			sum2 = 0.0;
// 			for(k=0;k<modes_size;k++){
// 				sum3 = 0.0;
// 				for(l=0;l<modes_size;l++){
// 					sum3 += get_gamma_tri(l,k)*get_eigen_tri(l);
// 				}
// 				sum2 += get_orthol_tri(j,k)*sum3;
// 				printf("[3]\t%d\t%e\n", k, sum3);
// 			}
// 			sum1 += get_orthol_tri(j,i)*sum2;
// 				printf("[2]\t%d\t%e\n", j, sum2);
// 		}
// 		results_m[0] = (double)i;
// 		results_m[1] = sum1;
// 		update_modes_tri(results_m);
// 		printf("[%d]\t%d\t%e\t%e\n", rank, i, get_eigen_tri(i), sum1);
// 	}
	
	output_modes_tri(directory);
	double sum1,sum2;
	for(i=0;i<l_size;i++){
		sum1 = 2.0*l_values[i]+1.0;
		sum2 = get_beam(i)*get_beam(i)*get_cl(i)+get_noise(i);
		for(n=0;n<modes_size;n++){
			result += a3[n]*triQ(n,i,i,i,i);
		}
		printf("%d\t%e\n", l_values[i], sum2*sum2*result/sum1);
	}
	
	MPI_Finalize();
}

