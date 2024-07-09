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

// set up timing
	double time1, time2, time3, time4, time5, duration;
	
	time1 = csecond();

// mpi vars
	int rank, nproc, lnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);


// read in data
	
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");
	
// 	load_flat(directory);
	
// 	double tau0 = get_ftau();
	
	int i,j,k;
	
	int number = 1;
	int lmax = 800;
	int lsize = lmax+1;
	int long lsize_long = lsize;
	
	double ***bispectrum = create_3Darray_long(lsize_long,lsize_long,lsize_long);
	
	double k1,k2,k3,s;
	
	for (i=0;i<lsize;i++) {
		k1 = (double)i;
		for (j=i;j<lsize;j++) {
			k2 = (double)j;
			for (k=j;k<lsize;k++) {
				k3 = (double)k;
				s = 0.0;
				if (k1+k2>=k3) s = smooth3(k1,k2,k3);
				bispectrum[i][j][k] = s;
				bispectrum[i][k][j] = s;
				bispectrum[j][i][k] = s;
				bispectrum[j][k][i] = s;
				bispectrum[k][i][j] = s;
				bispectrum[k][j][i] = s;
				
				if (k1==k2 && k2==k3) printf("%d\t%d\t%d\t%e\n",i,j,k,bispectrum[i][j][k]);
				
			}
		}
	}
	


	char l_number_file_int[100] = "";
	strcat(l_number_file_int, directory);
	strcat(l_number_file_int, interpolated_l_size_file);

	char l_file_int[100] = "";
	strcat(l_file_int, directory);
	strcat(l_file_int, interpolated_l_data_file);
	
	char reduced_bi_file_int[100] = "";
	strcat(reduced_bi_file_int, directory);
	strcat(reduced_bi_file_int, interpolated_bispectrum_file);
	
	int *l_values = create_ivector(lsize);
	
	for (i=0;i<lsize;i++) l_values[i] = i;
	
	long int bi_size = lsize*lsize*lsize;
	ivector_write(&number, l_number_file_int, &lsize);
	ivector_write(&lsize, l_file_int, &l_values[0]);
	array_write_long(&bi_size, reduced_bi_file_int, &bispectrum[0][0][0]);
	
	MPI_Finalize();
	
	return 0;	
}

