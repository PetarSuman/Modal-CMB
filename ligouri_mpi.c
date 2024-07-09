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

	int flag = 2;
// ****	flag = 0: source **** flag = 1: transfer **** flag = 2: source CMBFast ****
// set up timing
	double time1, time2, time3, time4, time5, duration;
	
	time1 = csecond();

// read in data
	
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");
	
// ini file
	char inifile[MAXLEN];
	
// mpi vars
	int rank, nproc, lnext, bnext;
	
// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	int i,j,k,r,s,t,m,n,l;

	if (flag == 0){
		
		char k_file[100] = "";
		strcat(k_file, directory);
		strcat(k_file, "KVAL.out");
		
		char tau_file[100] = "";
		strcat(tau_file, directory);
		strcat(tau_file, "TAUVAL.out");
		
		char src_file[100] = "";
		strcat(src_file, directory);
		strcat(src_file, "SOURCE.out");	
		
		double* kval = malloc( sizeof(double)*MAXLINES);
		double* tauval = malloc( sizeof(double)*MAXLINES);
		int* size = malloc( sizeof(int));
		
		int ksize, tausize, ksize_long, tausize_long;
		
		load_one_double(k_file, kval, size);
		ksize_long = *size;
		load_one_double(tau_file, tauval, size);
		tausize_long = *size;
		
		printf("%d %d\n", ksize_long, tausize_long);
		
		double* srcraw = malloc( sizeof(double)*MAXLINES);
		load_one_double(src_file, srcraw, size);
		
		ksize = ksize_long+1;
		for(i=1;i<ksize_long;i++){
			printf("%d %e\n", i, kval[i]);
			if(kval[i]==0){
				ksize = i+1;
				break;
			}
		}
		
		tausize = tausize_long;
		for(i=1;i<tausize_long;i++){
			printf("%d %e\n", i, tauval[i]);
			if(tauval[i]==0){
				tausize = i;
				break;
			}
		}
 		
		printf("%d %d\n", ksize, tausize);
		
		double **s_data = create_array(tausize, ksize);
		s_data[0][0] = srcraw[0];
		for(i=1;i<ksize;i++) s_data[0][i] = kval[i-1];
		for(i=1;i<tausize;i++) s_data[i][0] = tauval[i];
		
		for(j=1;j<ksize;j++){
			for(i=1;i<tausize;i++){
			 k = (j-1)*tausize_long + i + 1;
			 s_data[i][j] = srcraw[k];
			}
		}
		
		for(j=0;j<ksize;j++){
			for(i=0;i<tausize;i++){
			 printf("%e\t",s_data[i][j]);
			}
			printf("\n");
		}
		
		char s_file[100] = "";
		strcat(s_file, directory);
		strcat(s_file, source_data_file);
 		
		char s_size_file[100] = "";
		strcat(s_size_file, directory);
		strcat(s_size_file, source_size_file);
		
		int number = 2;
		
		int *source_sizes = create_ivector(number);
		source_sizes[0] = tausize;	
		source_sizes[1] = ksize;
    
    int ssize = ksize * tausize; 
    
    ivector_write(&number, s_size_file, &source_sizes[0]);
    destroy_ivector(source_sizes);
    
    array_write(&ssize, s_file, &s_data[0][0]);
    destroy_array(s_data, tausize);
		
	} else if (flag == 1) {
		
		char l_file[100] = "";
		strcat(l_file, directory);
		strcat(l_file, "lval.dat");	
		
		int* lvalues = malloc( sizeof(int)*MAXLINES);
		int* size = malloc( sizeof(int));
		
		load_one(l_file, lvalues, size);
		
		int l_size = *size;
		
		l = lvalues[0];
		char lval[4] = "";
		sprintf(lval, "%d", l);
		char l4[4] = "";
		if (l<10){
			strcat(l4, "000");
			strcat(l4, lval);
		} else if (l<100) {
			strcat(l4, "00");
			strcat(l4, lval);
		} else if (l<1000) {
			strcat(l4, "0");
			strcat(l4, lval);
		} else {
			strcat(l4, lval);
		}
		
		double* transfer = malloc( sizeof(double)*MAXLINES*2);
		
		char data_file[100] = "";
		strcat(data_file, directory);
		strcat(data_file, "tf_");
		strcat(data_file, l4);
		strcat(data_file, ".dat");
		
		load_two(data_file, transfer, size);
		
		int k_size = *size;
		double kvec[k_size];
		
		double **data_array = create_array(l_size+1,k_size+1);
		data_array[0][0] = 0;
		data_array[1][0] = l;
		
		j=0;
		for(k=0; k<k_size; k++){
			data_array[0][k+1] = transfer[j++];
			data_array[1][k+1] = transfer[j++];
		}
		
		for (i=1; i<l_size; i++){
			l = lvalues[i];
			char lval[4] = "";
			sprintf(lval, "%d", l);
			char l4[4] = "";
			if (l<10){
				strcat(l4, "000");
				strcat(l4, lval);
			} else if (l<100) {
				strcat(l4, "00");
				strcat(l4, lval);
			} else if (l<1000) {
				strcat(l4, "0");
				strcat(l4, lval);
			} else {
				strcat(l4, lval);
			}
			
			
			char data_file[100] = "";
			strcat(data_file, directory);
			strcat(data_file, "tf_");
			strcat(data_file, l4);
			strcat(data_file, ".dat");
			
			load_two(data_file, transfer, size);
			data_array[i+1][0] = l;
			
			j=0;
			for(k=0; k<k_size; k++){
				if (data_array[0][k+1] == transfer[j++])	data_array[i+1][k+1] = transfer[j++];
			}
		}
		
		for(l=0; l<l_size; l++){
			printf("%d\n", (int)data_array[l+1][0]);
		}
		
		
		
 		
		char t_file[100] = "";
		strcat(t_file, argv[1]);
		strcat(t_file, "/");
		strcat(t_file, transfer_data_file);
	 
		char t_size_file[100] = "";
		strcat(t_size_file, argv[1]);
		strcat(t_size_file, "/");
		strcat(t_size_file, transfer_size_file);
		
		int t_data_size = (l_size+1) * (k_size+1);
		array_write(&t_data_size, t_file, &data_array[0][0]);
		destroy_array(data_array, l_size);
		double **data_array2 = create_array(l_size+1,k_size+1);
		array_read(&t_data_size, t_file, &data_array2[0][0]);
		
		for(l=0; l<l_size; l++){
			printf("%d\n", (int)data_array[l+1][0]);
		}
		
		int number = 2;
		
		int *transfer_sizes = create_ivector(number);
		transfer_sizes[0] = l_size+1;
		transfer_sizes[1] = k_size+1;	
		ivector_write(&number, t_size_file, &transfer_sizes[0]);
		
		destroy_array(data_array, l_size);
		destroy_ivector(transfer_sizes);
	} else if (flag == 2) {
		
		char data_file[100] = "/home/cosmos/ccc-cam/jf334/output.dat";
		
		double* data = malloc( sizeof(double)*MAXLINES);
		long int* size = malloc( sizeof(long int));
		
		load_one_double_long(data_file, data, size);
		
		int k_size;
		int t_size;
		
		k_size = (int)data[0]+1;
		t_size = (int)data[1];
		
		
		double **source = create_array(t_size, k_size);
		
		source[0][0] = data[2];
		
		printf("k_size %d\n",k_size);
		printf("t_size %d\n",t_size);
		printf("tau_R  %e\n",source[0][0]);
		printf("\nK\n");
		
		long int m,n,tl,kl;
		
		tl = (long int)t_size;
		kl = (long int)k_size;
		
		for (i=1;i<k_size;i++) {
			source[0][i] = data[i+2];
// 			printf("%e\n",source[0][i]);
		}
		
		printf("\nTAU\n");
		
		for (i=1;i<t_size;i++) {
			source[i][0] = data[i+2+k_size];
// 			printf("%e\n",source[i][0]);
		}
		
		
		for (m=1;m<kl;m++) {
			for (n=1;n<tl;n++) {
				
				source[n][m] = data[n+m*tl+2+kl+tl]*(5.0/3.0);
// 				printf("%e\t",source[n][m]);
				
			}
			printf("\n");
		}
		
	 
		char s_file[100] = "";
		strcat(s_file, directory);
		strcat(s_file, source_data_file);
	 
		char s_size_file[100] = "";
		strcat(s_size_file, directory);
		strcat(s_size_file, source_size_file);
		
		int two = 2;
		int *source_sizes = create_ivector(two);
		source_sizes[0] = t_size;
		source_sizes[1] = k_size;
		ivector_write(&two, s_size_file, &source_sizes[0]);
		
		long int s_data_size = tl * kl;	
		array_write_long(&s_data_size, s_file, &source[0][0]);
		
	} else {
		printf("invalid flag\n");
	}
	
	MPI_Finalize();
	
	return 0;	
}

