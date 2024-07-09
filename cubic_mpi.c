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
// 	printf("%s\n", argv[2]);
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");
	
	int i,j,k,r,s,t,m,n;

	char l_number_file[100] = "";
	strcat(l_number_file, directory);
	strcat(l_number_file, l_size_file);

	char l_file[100] = "";
	strcat(l_file, directory);
	strcat(l_file, l_data_file);
	
	char reduced_bi_file[100] = "";
	strcat(reduced_bi_file, directory);
	strcat(reduced_bi_file, bispectrum_file);
	strcat(reduced_bi_file, "_f");
	strcat(reduced_bi_file, argv[2]);
	strcat(reduced_bi_file, "_");
	strcat(reduced_bi_file, argv[3]);
	
	printf("%s\n",reduced_bi_file);
	int number = 1;
	
	int *l_data_size = create_ivector(number);
	ivector_read(&number, l_number_file, &l_data_size[0]);

	int l_size_raw = l_data_size[0];
	
	int *l_values_raw = create_ivector(l_size_raw);
	ivector_read(&l_size_raw, l_file, &l_values_raw[0]);
	
	long int b_size_raw = l_size_raw;
	double ***bispectrum_raw = create_3Darray_long(b_size_raw,b_size_raw,b_size_raw);
	long int bi_size_raw = b_size_raw*b_size_raw*b_size_raw;
	array_read_long(&bi_size_raw, reduced_bi_file, &bispectrum_raw[0][0][0]);
	
	// Cut down to maximum.
	
	int l_max = 500;
	int limit = 0;
	for (i=0; i<l_size_raw; i++){
		if(l_values_raw[i]>l_max) break;
		limit = i+1;
	}
	l_size_raw = limit;
// 	for (i=0;i<l_size_raw;i++) printf("%d %d\n", i, l_values_raw[i]);
	
	long int b_size_raw_spun = l_size_raw;
	double ***bispectrum_raw_spun = create_3Darray_long(b_size_raw_spun,b_size_raw_spun,b_size_raw_spun);
	double test;
	
	//printf("\n%d\n\n", b_size_raw_spun);
	
	time2 = csecond();
	
	for (i = 0; i < l_size_raw; i++){
		if(bispectrum_raw[i][i][i]!=0.0)printf("%d\t%e\n", l_values_raw[i], bispectrum_raw[i][i][i]);
		for (j = 0; j < l_size_raw; j++){
			for (k = 0; k < l_size_raw; k++){
// 				if(bispectrum_raw[i][j][k]!=0.0)printf("%d\t%d\t%d\t%e\n", l_values_raw[i], l_values_raw[j], l_values_raw[k], bispectrum_raw[i][j][k]);
				test = (double)(i+j+k)/2;
				test = test-(int)test;
				//printf("%d %d %d %e\n", (j+k-i)/2, (k+i-j)/2, (i+j-k)/2, test);
				if (test == 0 && j+k-i>=0 && i+k-j>=0 && i+j-k>=0) {
					bispectrum_raw_spun[(j+k-i)/2][(k+i-j)/2][(i+j-k)/2] = bispectrum_raw[i][j][k];
// 					printf("%d\t%d\t%d\t%e\n", i, j, k);
				}
			}
		}
	}
	
	duration = csecond() - time2;
	printf("Finished rotation:\t%e\n", duration);
	time2 = csecond();
	
	
// interpolate data
	step *= 2;
	int start;
	double factor=1,gap,l1,l2,l3;
	int l_size = (l_size_raw - 1)*step + 1;
	double *l_values = create_vector(l_size);
	long int b_size = (l_size_raw - 1)*step/2 + 1;
	long int b_size_spun = l_size;
	double ***bispectrum_spun = create_3Darray_long(b_size_spun,b_size_spun,b_size_spun);
	double ***in = create_3Darray(4,4,4);
	double ***out = create_3Darray(step+1,step+1,step+1);
	int *cell = create_ivector(3);
	
	l_values[0] = 0;
	
	for (i = 1; i < l_size_raw; i++){
		l_values[i*step] = (double)l_values_raw[i];
		gap = l_values[i*step] - l_values[(i-1)*step];
		start = (i-1)*step;
		for (j=1; j<step+1; j++){
			l_values[start+j] = l_values[start] + j* (gap / step);
		}
	}
	
	int pt1, pt2, pt3, lpt1, lpt2, lpt3;

	for (i = 0; i < b_size_raw_spun-1; i++){
		time3 = csecond();
		for (j = i; j < b_size_raw_spun-1; j++){
			for (k = j; k < b_size_raw_spun-1; k++){
				
				cell[0] = 0;
				cell[1] = 0;
				cell[2] = 0;
				
				if (i==0) {cell[0] = -1;}
				if (j==0) {cell[1] = -1;}
				if (k==0) {cell[2] = -1;}
				if (i==b_size_raw_spun-2) {cell[0] = 1;}
				if (j==b_size_raw_spun-2) {cell[1] = 1;}
				if (k==b_size_raw_spun-2) {cell[2] = 1;}
				
				for(r=0;r<4;r++){
					pt1 = i-1+r-cell[0];
					for(s=0;s<4;s++){
						pt2 = j-1+s-cell[1];
						for(t=0;t<4;t++){
							pt3 = k-1+t-cell[2];
							
							lpt1 = pt2+pt3;
							lpt2 = pt3+pt1;
							lpt3 = pt1+pt2;
							
							factor = 0;
							if (lpt1 < l_size_raw && lpt2 < l_size_raw && lpt3 < l_size_raw ){
								l1 = l_values_raw[lpt1]*(l_values_raw[lpt1]+1);
								l2 = l_values_raw[lpt2]*(l_values_raw[lpt2]+1);
								l3 = l_values_raw[lpt3]*(l_values_raw[lpt3]+1);
								
								if(l1+l2+l3!=0)factor = (3*3*3*3*M_PI*M_PI*(l1*l2*l3))/(l1+l2+l3);
							}
							
							in[r][s][t] = bispectrum_raw_spun[pt1][pt2][pt3]*factor;
							
						}
					}
				}
				
				cubic_interpolation(in, out, cell);
				
				for(r=0;r<step+1;r++){
					pt1 = i*step+r;
					for(s=0;s<step+1;s++){
						pt2 = j*step+s;
						for(t=0;t<step+1;t++){
							pt3 = k*step+t;
							
							lpt1 = pt2+pt3;
							lpt2 = pt3+pt1;
							lpt3 = pt1+pt2;
							
							if (lpt1 < l_size && lpt2 < l_size && lpt3 < l_size ){
								l1 = l_values[lpt1]*(l_values[lpt1]+1);
								l2 = l_values[lpt2]*(l_values[lpt2]+1);
								l3 = l_values[lpt3]*(l_values[lpt3]+1);
							} else {
								l1=l2=l3=0;
							}
							
							factor = 0;
							if (l1*l2*l3!=0)factor = (l1+l2+l3)/(3*3*3*3*M_PI*M_PI*(l1*l2*l3));
							out[r][s][t] *= factor;
							
							bispectrum_spun[pt1][pt2][pt3] = out[r][s][t];
							bispectrum_spun[pt2][pt3][pt1] = out[r][s][t];
							bispectrum_spun[pt3][pt1][pt2] = out[r][s][t];
							bispectrum_spun[pt3][pt2][pt1] = out[r][s][t];
							bispectrum_spun[pt2][pt1][pt3] = out[r][s][t];
							bispectrum_spun[pt1][pt3][pt2] = out[r][s][t];
							
						}
					}
				}
				
			}
		}
		duration = csecond() - time3;
		printf("Finished %d\t%e\n", i, duration);
	}

	duration = csecond() - time2;
	printf("Finished interpolation:\t%e\n", duration);
	time2 = csecond();

	double ***bispectrum = create_3Darray_long(b_size,b_size,b_size);
	for (i = 0; i < b_size; i++){
		for (j = 0; j < b_size; j++){
			for (k = 0; k < b_size; k++){
				if (j+k-i>=0 && i+k-j>=0 && i+j-k>=0) {
					bispectrum[i][j][k] = bispectrum_spun[j+k-i][i+k-j][i+j-k];
				} else {
					bispectrum[i][j][k] = 0;
				}
			}
		}
	}
	time2 = csecond();

	duration = csecond() - time2;
	printf("Finished rotation back:\t%e\n", duration);
	time2 = csecond();
// output dense bispectrum file


	step = step/2;
	for(i = 0; i < l_size_raw; i++) printf("%d\t%e\t%e\n", l_values_raw[i], bispectrum_raw[i][i][i], bispectrum[i*step][i*step][i*step]);


	char l_number_file_int[100] = "";
	strcat(l_number_file_int, directory);
	strcat(l_number_file_int, interpolated_l_size_file);

	char l_file_int[100] = "";
	strcat(l_file_int, directory);
	strcat(l_file_int, interpolated_l_data_file);
	
	char reduced_bi_file_int[100] = "";
	strcat(reduced_bi_file_int, directory);
	strcat(reduced_bi_file_int, interpolated_bispectrum_file);
// 	strcat(reduced_bi_file_int, "_f");
// 	strcat(reduced_bi_file_int, argv[2]);
// 	strcat(reduced_bi_file_int, "_");
// 	strcat(reduced_bi_file_int, argv[3]);
	
	int l_size2 = b_size;
	int *l_values2 = create_ivector(l_size2);
	
	for (i=0;i<l_size2;i++) l_values2[i] = (int)l_values[2*i];
	
	long int bi_size = b_size*b_size*b_size;
	ivector_write(&number, l_number_file_int, &l_size2);
	ivector_write(&l_size2, l_file_int, &l_values2[0]);
	array_write_long(&bi_size, reduced_bi_file_int, &bispectrum[0][0][0]);


// 	int number9 = 1;
	
// 	int *l_data_size9 = create_ivector(number9);
// 	ivector_read(&number9, l_number_file_int, &l_data_size9[0]);
	
// 	int l_size9 = l_data_size9[0];
	
// 	int *l_values9 = create_ivector(l_size9);
// 	ivector_read(&l_size9, l_file_int, &l_values9[0]);
	
// 	long int l_size_long9 = l_size9;
// 	double ***bispectrum9 = create_3Darray(l_size_long9,l_size_long9,l_size_long9);
// 	long int bi_size9 = l_size_long9*l_size_long9*l_size_long9;
// 	printf("%d\n%ld\n%ld\n", l_size9, l_size_long9, bi_size9);
// 	array_read_long(&bi_size9, reduced_bi_file_int, &bispectrum9[0][0][0]);

// 	for(i = 0; i < l_size_raw; i++) printf("%d\t%e\t%e\n", l_values_raw[i], bispectrum_raw[i][i][i], bispectrum[step*i][step*i][step*i]);
	
// 	for (i=0;i<l_size_raw;i++){
// 		for (j=i;j<l_size_raw;j++){
// 			for (k=j;k<l_size_raw;k++){
// 				test = (double)(i+j+k)/2;
// 				test = test-(int)test;
// 				if(i+j>=k && test ==0.0)printf("%d\t%d\t%d\t%e\n", l_values_raw[i], l_values_raw[j], l_values_raw[k], bispectrum_raw[i][j][k]);
// 			}
// 		}
// 	}

// 	for (i=2;i<b_size;i++){
// 		for (j=i;j<b_size;j++){
// 			for (k=j;k<b_size;k++){
// 				test = (double)(i+j+k)/2;
// 				test = test-(int)test;
// 				if(i+j>=k && test ==0.0)printf("%d\t%d\t%d\t%e\n", l_values2[i], l_values2[j], l_values2[k], bispectrum[i][j][k]);
// 			}
// 		}
// 	}

	MPI_Finalize();
	
	duration = csecond() - time1;
	printf("Finished:\t%e\n", duration);
	
	return 0;	
}

