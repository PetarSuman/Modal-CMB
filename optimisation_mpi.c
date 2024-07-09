#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include "global.h"

int main( int argc, char *argv[] ){

// set up timing
	double time1, time2, time3, time4, time5, duration;
	
	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s dir\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
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
	
	set_terms_tri();
	read_eigen_tri(directory);
	read_beta_tri(directory);
	int lsize = get_bt_lsize();
	int *lvec = create_ivector(lsize);
	get_bt_lvec(lvec);
	create_cl(lsize);
	load_cl(directory, lsize, lvec);
	
	int xsize = get_bt_xsize();
	double *xvec = create_vector(xsize);
	get_bt_xvec(xvec);
	double **fisher = create_array(xsize,xsize);
	double **fishtmp = create_array(xsize,xsize);
	double *weight = create_vector(xsize);
	double distance;
	double sigma;

	int i,j,k,l,m,n;
	double sum1,sum2,sum3,sum4;
	double x1,x2,x3,x4;
	
// 	for(i=0;i<xsize;i++){
// 		for(j=0;j<lsize;j++){
// 			printf("%d\t%d\n",j,lvec[j]);
// 		}
// 			printf("\n");
// 	}
	
	time1 = csecond();
	
	double *results =  malloc( sizeof(double)*3);
	results[0] = 0.0;
	results[1] = 0.0;
	results[2] = 0.0;

	int entries = xsize*(xsize+1)/2;
	int jobs[entries][2];

	n=0;
	for(i=0;i<xsize;i++){
		for(j=i;j<xsize;j++){
			jobs[n][0]=i;
			jobs[n][1]=j;
			n++;
		}
	}
	create_fisher_tri(xsize);
	if ( rank == 0 ) record_tasks(entries);

	double c1,c2,c3,c4;
	int next = 0;
	double result;

	MPI_Barrier(MPI_COMM_WORLD);
	time1 = csecond();

	if (rank==0){
		for (l=0;l<lsize;l++){
			c1 = 0;
			for (n=0;n<get_terms_tri();n++){
				c1 += get_eigen_tri(n)*calculate_xint_tri(l, l, l, l, n, xsize);
			}
			printf("%d\t%e\n",l,c1);
			
		}
	}
	printf("[%d] starting\n", rank);
	while ( (next = get_next_task(18,3,rank,results)) < entries ) {

		i = jobs[next][0];
		j = jobs[next][1];

		result = calculate_fisher_tri(i,j);

		results[0] = (double)i;
		results[1] = (double)j;
		results[2] = result;

		duration = csecond() - time1;
// 		if(i==j-2)printf("[%d]\t%d\t%d\t%e\t%e\n", rank, i, j, result, duration);
		printf("%d\t%d\t%e\n", i, j, result);
		time1 = csecond();

	}
	
// 	end of MPI loop
	printf("[%d] finished\n", rank);

	if(rank==0){
		signal_end_tasks(18, 3);
	}
	MPI_Barrier(MPI_COMM_WORLD);

// 	xsize=187;
// 	int fisher_size=34969;
// 	double **fisher = create_array(xsize,xsize);
// 	double **fishtmp = create_array(xsize,xsize);
// 	double *weight = create_vector(xsize);
// 	array_read(&fisher_size,"/home/cosmos-tmp/jf334/Optimisation/Fisher.out",&fisher[0][0]);
	
	if (rank==0){
		sigma = 0.0;
		for(i=0;i<xsize;i++){
			for(j=0;j<xsize;j++){
				fisher[i][j] = get_fisher_tri(i,j);
// 				fisher[i][j] = i+j+1;
// 				printf("%d\t%d\t%e\n",i,j,fisher[i][j]);
				sigma += fisher[i][j];
// 				printf("Fisher: %d\t%d\t%e\n",i,j,fisher[i][j]);
			}
		}
		printf("Sigma: %e\n",sigma);
		
// 		for(i=0;i<xsize;i++){
// 			for(j=0;j<xsize;j++){
// 				fisher[i][j] = xsize*xsize*fisher[i][j]/sigma;
// 			}
// 		}
		
// 		sigma = 0.0;
// 		for(i=0;i<xsize;i++){
// 			for(j=0;j<xsize;j++){
// 				sigma += fisher[i][j];
// 			}
// 		}
		
		distance = 1.0;
		
		
		int index[xsize];
		int indtmp[xsize];
		for(i=0;i<xsize;i++){
			index[i] = -1;
			indtmp[i] = i;
		}
		
		
		for(n=0;n<xsize-1;n++){
		
			x1=0.0;
			for(i=n;i<xsize;i++){
			
				sum1=0.0;
				for(j=n;j<xsize;j++){
					sum1+=fisher[i][j];
				}
				x2 = sum1*sum1 / fisher[i][i];
				
				if(x2>x1){
					m = i;
					x1=x2;
				}
// 				printf("%d\t%e\t%e\n",i,x2,x1);
			}
			
			index[n]=indtmp[m];
			i = indtmp[n];
			indtmp[n]=indtmp[m];
			indtmp[m]=i;
			
			distance = distance - x1/sigma;
// 			distance = distance - x1;
			printf("Distance: %d\t%d\t%e\n",n,index[n],distance);
	
			for(i=0;i<xsize;i++){
				for(j=0;j<xsize;j++){
					fishtmp[i][j] = fisher[i][j];
				}
			}
			
//			Permute
			for(i=0;i<xsize;i++){
				x1 = fishtmp[n][i];
				fishtmp[n][i] = fishtmp[index[n]][i];
				fishtmp[index[n]][i] = x1;
			}
			for(i=0;i<xsize;i++){
				x1 = fishtmp[i][n];
				fishtmp[i][n] = fishtmp[i][index[n]];
				fishtmp[i][index[n]] = x1;
			}
			
			for(i=0;i<xsize;i++){
				for(j=1;j<xsize;j++){
					fisher[i][j] = fishtmp[i][j] - fishtmp[i][n]*fishtmp[n][j]/fishtmp[n][n];
				}
			}
			for(i=0;i<xsize;i++){
				fisher[i][n] = -fishtmp[i][n]/fishtmp[n][n];
			}
			for(i=0;i<xsize;i++){
				fisher[n][i] = -fishtmp[n][i]/fishtmp[n][n];
			}
			fisher[n][n] = -1.0/fishtmp[n][n];
			
			sum1 = 0.0;
			for(i=n+1;i<xsize;i++){
				for(j=n+1;j<xsize;j++){
					sum1 += fisher[i][j];
				}
			}
			
			if(distance<1e-10){
				printf("Weights are:\n");
				for(i=0;i<n+1;i++){
					weight[i] = 1.0;
					for(j=n+1;j<xsize;j++){
						weight[i] += -fisher[i][j];
					}
					printf("%d\t%e\n",index[i]+1,weight[i]);
				}
				sum1=0.0;
				sum2=0.0;
				for(i=0;i<n+1;i++){
					for(j=0;j<n+1;j++){
						sum1 += weight[i]*weight[j]*get_fisher_tri(index[i],index[j]);
// 						sum1 += weight[i]*weight[j]*(index[i]+index[j]+1);
					}
					for(j=0;j<xsize;j++){
						sum2 += weight[i]*get_fisher_tri(index[i],j);
// 						sum2 += weight[i]*(index[i]+j+1);
					}
				}
				
				printf("Correlator: %e\t%e\t%e\t%e\n",sigma,sum1,sum2,sum2/sqrt(sigma*sum1));
				
				double **data = create_array(lsize, n+1);
				int two = 2;

				char size_file[100] = "";
				strcat(size_file, directory);
				strcat(size_file, tri_proj_size_file);
				strcat(size_file, "*");

				int *sizes = create_ivector(two);
				sizes[0] = lsize;
				sizes[1] = n+1;
				ivector_write(&two, size_file, &sizes[0]);

				int data_size = lsize * (n+1);
				printf("data_size = %d\n",data_size);
// 				for (i=0; i<n+1; i++) data[0][i] = xvec[index[i]];
// 				for (l=0; l<lsize; l++) data[l][0] = (double)lvec[l];

				int qmax = get_qmax();
				printf("qmax = %d\n",qmax);
				for(m=0;m<qmax+1;m++){
				printf("m = %d\n",m);
					char suffix[2] = "";
					sprintf(suffix, "%d*", m);

					char filename[100] = "";
					strcat(filename, directory);
					strcat(filename, tri_proj_data_file);
					strcat(filename, suffix);

					data[0][0] = (double)m;
					for(i=0; i<n+1; i++){
						x1 = sqrt(xvec[index[i]]);
						x2 = copysign(pow(fabs(weight[i]),0.25),weight[i]);
						if(index[i]==0){
							x3 = xvec[index[i]+1] - xvec[index[i]];
						}else if(index[i]==xsize-1){
							x3 = xvec[index[i]] - xvec[index[i]-1];
						}else{
							x3 = (xvec[index[i]+1]-xvec[index[i]-1])/2.0;
						}
						x3 = pow(x3,0.25);
						
						for(l=0; l<lsize; l++){
							data[l][i] = x1*x2*x3*get_beta_tri(m,l,index[i]);
							if(m==0)printf("%d\t%d\t%e\n",i,lvec[l],data[l][i]);
						}
					}

					array_write(&data_size, filename, &data[0][0]);

				}
				
				break;
			}
		}
	}
	
	
	MPI_Finalize();
	
	return 0;	
}

