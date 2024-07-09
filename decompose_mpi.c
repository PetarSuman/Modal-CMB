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
	
	int i,j,k,r,s,t,m,n;

	char l_number_file[100] = "";
	strcat(l_number_file, directory);
	strcat(l_number_file, interpolated_l_size_file);

	char l_file[100] = "";
	strcat(l_file, directory);
	strcat(l_file, interpolated_l_data_file);
	
	char bi_file[100] = "";
	strcat(bi_file, directory);
	strcat(bi_file, interpolated_bispectrum_file);
	
	int number = 1;
	
	int *l_data_size = create_ivector(number);
	ivector_read(&number, l_number_file, &l_data_size[0]);

	int l_size = l_data_size[0];
	
	int *l_values = create_ivector(l_size);
	ivector_read(&l_size, l_file, &l_values[0]);
	
	long int l_size_long = l_size;
	double ***bispectrum = create_3Darray_long(l_size_long,l_size_long,l_size_long);
	long int bi_size = l_size_long*l_size_long*l_size_long;
	array_read_long(&bi_size, bi_file, &bispectrum[0][0][0]);
	
// 	for(i = 0; i < l_size9; i++) printf("%d\t%e\n", l_values[i], bispectrum[i][i][i]);
	
// generate basis

	double lvec[l_size];
	for(i=0;i<l_size;i++) lvec[i] = (double)l_values[i];
	create_basis(l_size, lvec);

	if (rank == 0){
		create_decompose();
	}
	
// calculate delta coefficents	
	double bi,p1,p2,p3;
	
	double l_max = l_values[l_size-1];
	
	
	int spl_size = l_size + 4;
	int wrk_size = 6 * l_size + 16;
	
	double x[l_size];
	double y1[l_size];
	double y2[l_size];
	double y3[l_size];
	double spl_k[spl_size];
	double spl_c[spl_size];
	double wrk[wrk_size];
	int ifail;
	double result;
	double factor;
	
	for (r=0;r<l_size;r++) {
		x[r] = l_values[r]/l_max;
	}
	
	y1[0] = 1;
	y2[0] = 1;
	y3[0] = 1;

	if ( rank == 0 ){
	double sum;
		for (i=0;i<alpha_max+1;i++) {
			sum = 0;
			for (r=0;r<l_size;r++) {
				y1[r] = get_basis(r,i)*get_basis(r,i)*(2*i+1);
			}
			e01baf_(&l_size,x,y1,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
			e02bdf_(&spl_size,spl_k,spl_c,&sum,&ifail);
			printf("%d\t%e\n", i, sum);
		}
	}
	
	int delta_total = (alpha_max+1)*(alpha_max+2)*(alpha_max+3)/6;
	int delta_triples[delta_total][3];
	
	n=0;	
	for (i=0;i<alpha_max+1;i++) {
		for (j=i;j<alpha_max+1;j++) {
			for (k=j;k<alpha_max+1;k++) {
				delta_triples[n][0] = i;
				delta_triples[n][1] = j;
				delta_triples[n][2] = k;
				n++;
			}
		}
	}
	
	double *results =  malloc( sizeof(double)*4);
	results[0] = 0.0;
	results[1] = 0.0;
	results[2] = 0.0;
	results[3] = 0.0;

	// master record total number of tasks	
	if ( rank == 0 ) record_tasks(delta_total);

	create_decompose();

	int s1,s2,s3;
	double scale;

	// mpi loop over the lvalues 
	while ( (lnext = get_next_task(4,4,rank,results)) < delta_total ) {
		
		time1 = csecond();
		 
		i = delta_triples[lnext][0];
		j = delta_triples[lnext][1];
		k = delta_triples[lnext][2];
		
		factor = (2*i+1)*(2*j+1)*(2*k+1);
		
		for (r=0;r<l_size;r++) {
			
			l1 = l_values[r];
// 			s1 = l1*(l1+1);
			s1 = l1*l1;
			p1 = get_basis(r,i);
			
			for (s=0;s<l_size;s++) {
				
				l2 = l_values[s];
// 				s2 = l2*(l2+1);
				s2 = l2*l2;
				p2 = get_basis(s,j);
				
				for (t=0;t<l_size;t++) {
					
					l3 = l_values[t];
// 					s3 = l3*(l3+1);
					s3 = l3*l3;
					p3 = get_basis(t,k);
					
					scale = 0.0;
// 					if(s1+s2+s3!=0) scale = (double)(s1*s2*s3) / (double)(s1+s2+s3);
					scale = (double)s1*(double)s2*(double)s3;
					
					bi = bispectrum[r][s][t] * scale;
// 					if (l1==l2 && l2==l3) printf("%d\t%d\t%d\t%e\n",l1,l2,l3,bi);
					y3[t] = bi*p3;
				}
				
				e01baf_(&l_size,x,y3,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
				e02bdf_(&spl_size,spl_k,spl_c,&result,&ifail);
				
				y2[s] = result*p2;
			}
			
			e01baf_(&l_size,x,y2,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
			e02bdf_(&spl_size,spl_k,spl_c,&result,&ifail);
			
			y1[r] = result*p1;
		}
		
		e01baf_(&l_size,x,y1,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
		e02bdf_(&spl_size,spl_k,spl_c,&result,&ifail);
		
		result *= factor;
		
		results[0] = (double)i;
		results[1] = (double)j;
		results[2] = (double)k;
		results[3] = result;
		
		time2 = csecond();
		
		printf("[%d] task %d complete for %d %d %d %e in %e\n", rank, lnext, i, j, k, result, time2-time1);
		
	}
	
// end of MPI loop
	
	if (rank == 0){
		signal_end_tasks(4,4);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
// 	destroy_3Darray_long(bispectrum, l_size_long, l_size_long);
	
	time3 = csecond();
	
// reconstruct bispectrum

	if (rank == 0){
		output_decompose(directory);
		
		if (dflag_out == 1){
			printf("Starting\n");
			double sum1, sum2, sum3, delta;
			
// 			read_decompose(directory);
			
			double ***bispectrum_rec = create_3Darray(l_size,l_size,l_size);
			for (r=0;r<l_size;r++) {
				time1 = csecond();
				l1 = l_values[r];
// 				s1 = l1*(l1+1);
				s1 = l1*l1;
				for (s=r;s<l_size;s++) {
					l2 = l_values[s];
// 					s2 = l2*(l2+1);
					s2 = l2*l2;
					for (t=s;t<l_size;t++) {
						l3 = l_values[t];
// 						s3 = l3*(l3+1);
						s3 = l3*l3;
						
						scale = 0.0;
// 						if(s1*s2*s3!=0) scale = (double)(s1+s2+s3) / (s1*s2*s3);
						if((double)s1*(double)s2*(double)s3!=0) scale = 1.0 / ((double)s1*(double)s2*(double)s3);
						
						sum1 = 0;
						for (i=0;i<dflag_dmax+1;i++) {
							sum2 = 0;
							p1 = get_basis(r,i);
							for (j=0;j<dflag_dmax+1;j++) {
								sum3 = 0;
								p2 = get_basis(s,j);
								for (k=0;k<dflag_dmax+1;k++) {
									p3 = get_basis(t,k);
									delta = get_decompose(i,j,k);
									sum3 += delta * p3;	
								}
								sum2 += sum3 * p2;
							}
							sum1 += sum2 * p1;
						}
						sum1 *= scale;
						
						bispectrum_rec[r][s][t] = sum1;
						bispectrum_rec[s][t][r] = sum1;
						bispectrum_rec[t][r][s] = sum1;
						bispectrum_rec[t][s][r] = sum1;
						bispectrum_rec[s][r][t] = sum1;
						bispectrum_rec[r][t][s] = sum1;
						
						if (l1==l2 && l2==l3) printf("%d\t%d\t%d\t%e\t%e\n",r,s,t,bispectrum[r][s][t],bispectrum_rec[r][s][t]);
						
					}
				}
				time2 = csecond();
// 				printf("Done %d in %e\n",l1,time2-time1);
			}
			
			char dmaxs[1];
			sprintf(dmaxs,"%d",dflag_dmax);
			char bi_fileout[100] = "";
			strcat(bi_fileout, directory);
			strcat(bi_fileout, "bispectrum_int");
			strcat(bi_fileout, dmaxs);
			array_write_long(&bi_size, bi_fileout, &bispectrum_rec[0][0][0]);
		}
	}

	printf("[%d] Done.\n", rank);
	
	MPI_Finalize();
	
	return 0;	
}

