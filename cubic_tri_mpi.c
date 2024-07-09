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
	
	int i,j,k,l,r,s,t,u,m,n;

	char l_number_file[100] = "";
	strcat(l_number_file, directory);
	strcat(l_number_file, l_size_file);

	char l_file[100] = "";
	strcat(l_file, directory);
	strcat(l_file, l_data_file);
	
	char reduced_tri_file[100] = "";
	strcat(reduced_tri_file, directory);
	strcat(reduced_tri_file, trispectrum_file);
	
	printf("%s\n",reduced_tri_file);
	int number = 1;
	
	int *l_data_size = create_ivector(number);
	ivector_read(&number, l_number_file, &l_data_size[0]);

	int l_size_raw = l_data_size[0];
	
	int *l_values_raw = create_ivector(l_size_raw);
	ivector_read(&l_size_raw, l_file, &l_values_raw[0]);
	
	long int t_size_raw = l_size_raw;
	double ****trispectrum_raw = create_4Darray_long(t_size_raw,t_size_raw,t_size_raw,t_size_raw);
	long int tri_size_raw = t_size_raw*t_size_raw*t_size_raw*t_size_raw;
	array_read_long(&tri_size_raw, reduced_tri_file, &trispectrum_raw[0][0][0][0]);
	
	// Cut down to maximum.
	
	int l_max = 500;
	int limit = 0;
	for (i=0; i<l_size_raw; i++){
		if(l_values_raw[i]>l_max) break;
		limit = i+1;
	}
	l_size_raw = limit;
	for (i=0;i<l_size_raw;i++) printf("%d %d\n", i, l_values_raw[i]);
	
	long int t_size_raw_spun = l_size_raw;
	double ****trispectrum_raw_spun = create_4Darray_long(t_size_raw_spun,t_size_raw_spun,t_size_raw_spun,t_size_raw_spun);
	double test;
	
	//printf("\n%d\n\n", t_size_raw_spun);
	
	time2 = csecond();
	
	for (i = 0; i < l_size_raw; i++){
		for (j = 0; j < l_size_raw; j++){
			for (k = 0; k < l_size_raw; k++){
				for (l = 0; l < l_size_raw; l++){
				test = (double)(i+j+k+l)/2;
				test = test-(int)test;
				if (test == 0 && j+k+l-i>=0 && i+k+l-j>=0 && i+j+l-k>=0 && i+j+k-l>=0) {
					trispectrum_raw_spun[(j+k+l-i)/2][(k+i+l-j)/2][(i+j+l-k)/2][(i+j+k-l)/2] = trispectrum_raw[i][j][k][l];
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
	double factor=1,gap,l1,l2,l3,l4;
	int l_size = (l_size_raw - 1)*step + 1;
	double *l_values = create_vector(l_size);
	long int t_size = (l_size_raw - 1)*step/2 + 1;
	long int t_size_spun = l_size;
	double ****trispectrum_spun = create_4Darray_long(t_size_spun,t_size_spun,t_size_spun,t_size_spun);
	double ****in = create_4Darray(4,4,4,4);
	double ****out = create_4Darray(step+1,step+1,step+1,step+1);
	int *cell = create_ivector(4);
	
	l_values[0] = 0;
	
	for (i = 1; i < l_size_raw; i++){
		l_values[i*step] = (double)l_values_raw[i];
		gap = l_values[i*step] - l_values[(i-1)*step];
		start = (i-1)*step;
		for (j=1; j<step+1; j++){
			l_values[start+j] = l_values[start] + j* (gap / step);
		}
	}
	
	int pt1, pt2, pt3, pt4, lpt1, lpt2, lpt3, lpt4;

	for (i = 0; i < t_size_raw_spun-1; i++){
		time3 = csecond();
		for (j = i; j < t_size_raw_spun-1; j++){
			for (k = j; k < t_size_raw_spun-1; k++){
				for (l = k; l < t_size_raw_spun-1; l++){
					
					cell[0] = 0;
					cell[1] = 0;
					cell[2] = 0;
					cell[3] = 0;
					
					if (i==0) {cell[0] = -1;}
					if (j==0) {cell[1] = -1;}
					if (k==0) {cell[2] = -1;}
					if (l==0) {cell[3] = -1;}
					if (i==t_size_raw_spun-2) {cell[0] = 1;}
					if (j==t_size_raw_spun-2) {cell[1] = 1;}
					if (k==t_size_raw_spun-2) {cell[2] = 1;}
					if (l==t_size_raw_spun-2) {cell[3] = 1;}
					
					for(r=0;r<4;r++){
						pt1 = i-1+r-cell[0];
						for(s=0;s<4;s++){
							pt2 = j-1+s-cell[1];
							for(t=0;t<4;t++){
								pt3 = k-1+t-cell[2];
								for(u=0;u<4;u++){
									pt4 = l-1+u-cell[3];
									
									lpt1 = pt2+pt3+pt4;
									lpt2 = pt3+pt1+pt4;
									lpt3 = pt1+pt2+pt4;
									lpt4 = pt1+pt2+pt3;
									
									factor = 0;
									if (lpt1 < l_size_raw && lpt2 < l_size_raw && lpt3 < l_size_raw && lpt4 < l_size_raw ){
										l1 = l_values_raw[lpt1]*(l_values_raw[lpt1]+1);
										l2 = l_values_raw[lpt2]*(l_values_raw[lpt2]+1);
										l3 = l_values_raw[lpt3]*(l_values_raw[lpt3]+1);
										l4 = l_values_raw[lpt4]*(l_values_raw[lpt4]+1);
										
										if(l1+l2+l3+l4!=0)factor = (3*3*3*3*M_PI*M_PI*(l1*l2*l3))/(l1+l2+l3);
									}
									
									in[r][s][t][u] = trispectrum_raw_spun[pt1][pt2][pt3][pt4]*factor;
								}	
							}
						}
					}
					
					cubic_interpolation_tri(in, out, cell);
					
					for(r=0;r<step+1;r++){
						pt1 = i*step+r;
						for(s=0;s<step+1;s++){
							pt2 = j*step+s;
							for(t=0;t<step+1;t++){
								pt3 = k*step+t;
								for(u=0;u<step+1;u++){
									pt4 = l*step+u;
									
									lpt1 = pt2+pt3+pt4;
									lpt2 = pt3+pt1+pt4;
									lpt3 = pt1+pt2+pt4;
									lpt4 = pt1+pt2+pt3;
									
									if (lpt1 < l_size && lpt2 < l_size && lpt3 < l_size && lpt4 < l_size){
										l1 = l_values[lpt1]*(l_values[lpt1]+1);
										l2 = l_values[lpt2]*(l_values[lpt2]+1);
										l3 = l_values[lpt3]*(l_values[lpt3]+1);
										l4 = l_values[lpt4]*(l_values[lpt4]+1);
									} else {
										l1=l2=l3=0;
									}
									
									factor = 0;
									if (l1*l2*l3*l4!=0)factor = (l1+l2+l3)/(3*3*3*3*M_PI*M_PI*(l1*l2*l3));
									out[r][s][t][u] *= factor;
									
									trispectrum_spun[pt1][pt2][pt3][pt4] = out[r][s][t][u];
									trispectrum_spun[pt1][pt2][pt4][pt3] = out[r][s][t][u];
									trispectrum_spun[pt1][pt3][pt2][pt4] = out[r][s][t][u];
									trispectrum_spun[pt1][pt3][pt4][pt2] = out[r][s][t][u];
									trispectrum_spun[pt1][pt4][pt2][pt3] = out[r][s][t][u];
									trispectrum_spun[pt1][pt4][pt3][pt2] = out[r][s][t][u];
									trispectrum_spun[pt2][pt1][pt3][pt4] = out[r][s][t][u];
									trispectrum_spun[pt2][pt1][pt4][pt3] = out[r][s][t][u];
									trispectrum_spun[pt2][pt3][pt2][pt4] = out[r][s][t][u];
									trispectrum_spun[pt2][pt3][pt4][pt2] = out[r][s][t][u];
									trispectrum_spun[pt2][pt4][pt1][pt3] = out[r][s][t][u];
									trispectrum_spun[pt2][pt4][pt3][pt1] = out[r][s][t][u];
									trispectrum_spun[pt3][pt1][pt2][pt4] = out[r][s][t][u];
									trispectrum_spun[pt3][pt1][pt4][pt2] = out[r][s][t][u];
									trispectrum_spun[pt3][pt2][pt1][pt4] = out[r][s][t][u];
									trispectrum_spun[pt3][pt2][pt4][pt1] = out[r][s][t][u];
									trispectrum_spun[pt3][pt4][pt1][pt2] = out[r][s][t][u];
									trispectrum_spun[pt3][pt4][pt2][pt1] = out[r][s][t][u];
									trispectrum_spun[pt4][pt1][pt2][pt3] = out[r][s][t][u];
									trispectrum_spun[pt4][pt1][pt3][pt2] = out[r][s][t][u];
									trispectrum_spun[pt4][pt2][pt1][pt3] = out[r][s][t][u];
									trispectrum_spun[pt4][pt2][pt3][pt1] = out[r][s][t][u];
									trispectrum_spun[pt4][pt3][pt3][pt2] = out[r][s][t][u];
									trispectrum_spun[pt4][pt3][pt2][pt1] = out[r][s][t][u];
								}
							}
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

	double ****trispectrum = create_4Darray_long(t_size,t_size,t_size,t_size);
	for (i = 0; i < t_size; i++){
		for (j = 0; j < t_size; j++){
			for (k = 0; k < t_size; k++){
				for (l = 0; l < t_size; l++){
					if (j+k+l-i>=0 && i+k+l-j>=0 && i+j+l-k>=0 && i+j+k-l>=0){
						trispectrum[i][j][k][l]= trispectrum_spun[j+k+l-i][i+k+l-j][i+j+l-k][i+j+k-l];
					} else {
						trispectrum[i][j][k][l] = 0;
					}
				}
			}
		}
	}
	time2 = csecond();

	duration = csecond() - time2;
	printf("Finished rotation back:\t%e\n", duration);
	time2 = csecond();
	
// output dense trispectrum file

	step = step/2;
	for(i = 0; i < l_size_raw; i++) printf("%d\t%e\t%e\n", l_values_raw[i], trispectrum_raw[i][i][i][i], trispectrum[i*step][i*step][i*step][i*step]);


	char l_number_file_int[100] = "";
	strcat(l_number_file_int, directory);
	strcat(l_number_file_int, interpolated_l_size_file);

	char l_file_int[100] = "";
	strcat(l_file_int, directory);
	strcat(l_file_int, interpolated_l_data_file);
	
	char reduced_tri_file_int[100] = "";
	strcat(reduced_tri_file_int, directory);
	strcat(reduced_tri_file_int, interpolated_trispectrum_file);
	
	int l_size2 = t_size;
	int *l_values2 = create_ivector(l_size2);
	
	for (i=0;i<l_size;i++) l_values2[i] = (int)l_values[2*i];
	
	long int tri_size = t_size*t_size*t_size*t_size;
	ivector_write(&number, l_number_file_int, &l_size2);
	ivector_write(&l_size2, l_file_int, &l_values2[0]);
	array_write_long(&tri_size, reduced_tri_file_int, &trispectrum[0][0][0][0]);

	MPI_Finalize();
	
	duration = csecond() - time1;
	printf("Finished:\t%e\n", duration);
	
	return 0;	
}

