#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include "global.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
// #include <nagmk21.h>
// #include <nag.h>
// #include <nagf03.h>

int main( int argc, char *argv[] ){
	
// **1**

	double time1, time2, time3, time4, time5, duration;


//  	if (argc != 2) {
	if (argc != 4) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s data filename models\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");
	
	char filename[100] = "";
	strcat(filename, argv[2]);
	
	int models = atoi(argv[3]);
	
	int number = 1;

	int i,j,k,l,m,n,r,s;
	
	int l_size = 2001;
	int *l_values = create_ivector(l_size);
	for(i=0;i<l_size;i++){
		l_values[i] = i;
	}
	
	double lmax = (double)l_values[l_size-1];
	init_lmax((int)lmax);
	
	long int l_size_long = l_size;
	set_terms_late(directory);
	create_cl(l_size);
	create_beam(l_size);
	create_noise(l_size);
	create_t_wgt(l_size);
	create_lens(l_size);
	load_cl(directory, l_size, l_values);
	load_BN(directory, l_size, l_values);
	load_TL(directory, l_size, l_values);
	load_lens(directory, l_size, l_values);
	int xsize = (int)lmax+1;
	
// mpi vars
	int rank, nproc, onext, mnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	if(rank==0)printf("data file: %s\n", filename);
	
	double result=0;
	double x;

	double xvec[xsize];
	for(i=0;i<xsize;i++){
		xvec[i] = (double)i;
	}
	
	double xmax = (double)lmax;
	
	create_basis_late(xsize, xmax, xvec);
	
	int ortho_size = get_terms_late();
	
	if(rank==0){
		printf("lmax: %d\n", (int)lmax);
		printf("pmax: %d\n", get_pmax_late());
		printf("terms: %d\n", get_terms_late());
		for(n=0;n<ortho_size;n++){
			find_perm_late(n,&i,&j,&k);
			printf("%d\t%d\t%d\t%d\n",n,i,j,k);
		}
	}
	
	int loops;
	int auxloop ;
	int start_loop;
	int end_loop;
	
	if(oflag_load==0){
		
		int ortho_total = (ortho_size)*(ortho_size+1)/2;
		int ortho_pairs[ortho_total][2];
		double **ortho = (double **)create_array(ortho_size,ortho_size);
		
		n=0;	
		for (i=0;i<ortho_size;i++) {
			for (j=i;j<ortho_size;j++) {
				ortho_pairs[n][0] = i;
				ortho_pairs[n][1] = j;
				n++;
			}
		}
		
		if(eflag_order_late!=6){
			loops=ortho_total/nproc;
			auxloop = fmod(ortho_total,nproc);
			start_loop = rank*loops;
			end_loop = (rank+1)*loops-1;
		
			if (auxloop != 0){
				if (rank < auxloop){
					start_loop = start_loop + rank;
					end_loop = end_loop + rank + 1;
				}else{
					start_loop = start_loop + auxloop;
					end_loop = end_loop + auxloop;
				}
			}
 		
			printf("mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
			
			for (i=0;i<ortho_size;i++) {
				for (j=0;j<ortho_size;j++) {
					ortho[i][j] = 0.0;
				}
			}
				
			init_orthol_glint();
		
			time1 = csecond();
			for(n=start_loop;n<end_loop+1;n++){
				
				i = ortho_pairs[n][0];
				j = ortho_pairs[n][1];
				
				ortho[i][j] = calculate_orthol(i,j);
				
				duration = csecond() - time1;
				printf("QQ\t%d\t%d\t%e\t%e\n", i, j, ortho[i][j], duration);
				time1 = csecond();
				
			} // end of MPI loop
			printf("[%d] finished\n", rank);
			
			double* ortho_send = (double *)malloc( ortho_total * sizeof(double) );
			double* ortho_recv = (double *)malloc( ortho_total * sizeof(double) );
			
			MPI_Barrier(MPI_COMM_WORLD);
		
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
					ortho_send[n] = ortho[i][j];
					n++;
				}
			}
			
	        	MPI_Reduce(ortho_send,ortho_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
			
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==0){
				n=0;
				for (i=0;i<ortho_size;i++) {
					for (j=i;j<ortho_size;j++) {
						ortho[i][j] = ortho_recv[n];
						if(i!=j)ortho[j][i] = ortho[i][j];
						n++;
					}
				}
			}
		}else{
		
			m=0;
			n=0;
			for (i=0;i<ortho_total;i++) {
				r=ortho_pairs[i][0];
				s=ortho_pairs[i][1];
				if((r>1&&r<5)||(s>1&&s<5)){
					m++;
				}else{
					n++;
				}
			}
			
			int ortho_3D = m;
			int ortho_pairs_3D[ortho_3D][2];
				
			int ortho_GL = n;
			int ortho_pairs_GL[ortho_GL][2];

			m=0;
			n=0;
			for (i=0;i<ortho_total;i++) {
				r=ortho_pairs[i][0];
				s=ortho_pairs[i][1];
				if((r>1&&r<5)||(s>1&&s<5)){
					ortho_pairs_3D[m][0] = r;
					ortho_pairs_3D[m][1] = s;
					m++;
				}else{
					ortho_pairs_GL[n][0] = r;
					ortho_pairs_GL[n][1] = s;
					n++;
				}
			}
			
			for (i=0;i<ortho_size;i++) {
				for (j=0;j<ortho_size;j++) {
					ortho[i][j] = 0.0;
				}
			}
			
			loops=ortho_3D/nproc;
			auxloop = fmod(ortho_3D,nproc);
			start_loop = rank*loops;
			end_loop = (rank+1)*loops-1;
		
			if (auxloop != 0){
				if (rank < auxloop){
					start_loop = start_loop + rank;
					end_loop = end_loop + rank + 1;
				}else{
					start_loop = start_loop + auxloop;
					end_loop = end_loop + auxloop;
				}
			}
 		
			printf("3D mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
			MPI_Barrier(MPI_COMM_WORLD);
		
			time1 = csecond();
			for(n=start_loop;n<end_loop+1;n++){
				
				i = ortho_pairs_3D[n][0];
				j = ortho_pairs_3D[n][1];
				
				ortho[i][j] = calculate_orthol_3D(i,j);
				
				duration = csecond() - time1;
				printf("QQ\t%d\t%d\t%e\t%e\n", i, j, ortho[i][j], duration);
				time1 = csecond();
				
			} // end of MPI loop
			printf("[%d] finished\n", rank);
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			loops=ortho_GL/nproc;
			auxloop = fmod(ortho_GL,nproc);
			start_loop = rank*loops;
			end_loop = (rank+1)*loops-1;
		
			if (auxloop != 0){
				if (rank < auxloop){
					start_loop = start_loop + rank;
					end_loop = end_loop + rank + 1;
				}else{
					start_loop = start_loop + auxloop;
					end_loop = end_loop + auxloop;
				}
			}
 		
			printf("GL mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
			MPI_Barrier(MPI_COMM_WORLD);
				
			init_orthol_glint();
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			time1 = csecond();
			for(n=start_loop;n<end_loop+1;n++){
				
				i = ortho_pairs_GL[n][0];
				j = ortho_pairs_GL[n][1];
				
				ortho[i][j] = calculate_orthol(i,j);
				
				duration = csecond() - time1;
				printf("QQ\t%d\t%d\t%e\t%e\n", i, j, ortho[i][j], duration);
				time1 = csecond();
				
			} // end of MPI loop
			printf("[%d] finished\n", rank);
			
			double* ortho_send = (double *)malloc( ortho_total * sizeof(double) );
			double* ortho_recv = (double *)malloc( ortho_total * sizeof(double) );
			
			MPI_Barrier(MPI_COMM_WORLD);
		
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
					ortho_send[n] = ortho[i][j];
					n++;
				}
			}
			
	        	MPI_Reduce(ortho_send,ortho_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
			
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==0){
				n=0;
				for (i=0;i<ortho_size;i++) {
					for (j=i;j<ortho_size;j++) {
						ortho[i][j] = ortho_recv[n];
						if(i!=j)ortho[j][i] = ortho[i][j];
						n++;
					}
				}
			}
		}
		
		double *results_o =  malloc( sizeof(double)*3);
		results_o[0] = 0.0;
		results_o[1] = 0.0;
		results_o[2] = 0.0;
		
		create_orthol();
		create_lambdal();
			
		if(rank==0){
			gsl_matrix *cholesky = gsl_matrix_alloc(ortho_size,ortho_size);
			gsl_matrix *matrixin = gsl_matrix_alloc(ortho_size,ortho_size);
			gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(ortho_size);
			gsl_vector *eigen = gsl_vector_alloc(ortho_size);
			gsl_matrix *eigenv = gsl_matrix_alloc(ortho_size,ortho_size);
			double x1;

			gsl_matrix_set_zero(cholesky);
			gsl_matrix_set_zero(matrixin);
			gsl_vector_set_zero(eigen);
			gsl_matrix_set_zero(eigenv);
			double norm[ortho_size];
				
			for (i=0;i<ortho_size;i++) {
				norm[i] = sqrt(ortho[i][i]);
				printf("%d\t%e\n",i,norm[i]);
			}
	
			for (i=0;i<ortho_size;i++) {
				for (j=0;j<ortho_size;j++) {
					x1 = ortho[i][j]/(norm[i]*norm[j]);
					gsl_matrix_set(cholesky,i,j,x1);
					gsl_matrix_set(matrixin,i,j,x1);
				}
			}
				
			gsl_eigen_symmv(matrixin,eigen,eigenv,workspace);
			gsl_eigen_symmv_sort(eigen,eigenv,GSL_EIGEN_SORT_VAL_DESC);

			for(i=0;i<ortho_size;i++){
				printf("%d\t%e\n",i,gsl_vector_get(eigen,i));
				
			}

			if(gflag_pca==0||gsl_vector_get(eigen,ortho_size-1)<0){
				
				int r_size;
				for(i=0;i<ortho_size && gsl_vector_get(eigen,i)>0.0;i++){
					r_size = i+1;
				}
				
				double **lambda = (double **)create_array(ortho_size,ortho_size);
				double **lambdainv = (double **)create_array(ortho_size,ortho_size);
				
				for(i=0;i<ortho_size;i++){
					for(j=0;j<ortho_size;j++){
						if(i<r_size){
							lambda[i][j] = gsl_matrix_get(eigenv,j,i)/(norm[j]*sqrt(gsl_vector_get(eigen,i)));
						}else{
							lambda[i][j] = 0.0;
						}
						if(j<r_size){
							lambdainv[i][j] = norm[i]*gsl_matrix_get(eigenv,i,j)*sqrt(gsl_vector_get(eigen,j));
						}else{
							lambdainv[i][j] = 0.0;
						}
					}
				}
				
				for(i=0;i<ortho_size;i++){
					for(j=0;j<ortho_size;j++){
						results_o[0] = (double)i;
						results_o[1] = (double)j;
			 			results_o[2] = lambdainv[i][j];
						update_orthol(results_o);
					}
				}
				
				for(i=0;i<ortho_size;i++){
					for(j=0;j<ortho_size;j++){
						results_o[0] = (double)i;
						results_o[1] = (double)j;
			 			results_o[2] = lambda[i][j];
		 				update_lambdal(results_o);
					}
				}
				
			}else{
				
				gsl_linalg_cholesky_decomp(cholesky);
			
				for(i=0;i<ortho_size;i++){
					for(j=0;j<ortho_size;j++){
						results_o[0] = (double)i;
						results_o[1] = (double)j;
		 				if(i>=j){
		 					results_o[2] = norm[i]*gsl_matrix_get(cholesky, i, j);
		 				} else {
		 					results_o[2] = 0;
		 				}
						update_orthol(results_o);
					}
				}
				
				gsl_matrix *upper = gsl_matrix_alloc(ortho_size,ortho_size);
				gsl_matrix *inverse = gsl_matrix_alloc(ortho_size,ortho_size);
				gsl_matrix_set_zero(upper);
				
				
				for(i=0;i<ortho_size;i++){
					for(j=i;j<ortho_size;j++){
			 			x1 = gsl_matrix_get(cholesky, i, j);
						gsl_matrix_set(upper,i,j,x1);
		 			}
		 		}
				
				int s;
				gsl_permutation *perm = gsl_permutation_calloc(ortho_size);
				gsl_linalg_LU_invert(upper,perm,inverse);
				
				for(i=0;i<ortho_size;i++){
					for(j=0;j<ortho_size;j++){
						results_o[0] = (double)i;
						results_o[1] = (double)j;
						results_o[2] = gsl_matrix_get(inverse, j, i)/norm[j];
						update_lambdal(results_o);
					}
				}
			
			}
			output_lambdal(directory);
			output_orthol(directory);
		}
		
 	}
 	
 	MPI_Barrier(MPI_COMM_WORLD);

	read_orthol(directory);
	read_lambdal(directory);
		
	double **orthoinv = (double **)create_array(ortho_size,ortho_size);
	double **lambda = (double **)create_array(ortho_size,ortho_size);
	double **lambdainv = (double **)create_array(ortho_size,ortho_size);

	for(i=0;i<ortho_size;i++){
		for(j=0;j<ortho_size;j++){
			if(i<j){
				lambdainv[i][j] = 0.0;
			}else{
				lambdainv[i][j] = get_orthol(i,j);
			}
		}
	}

	for(i=0;i<ortho_size;i++){
		for(j=0;j<ortho_size;j++){
			orthoinv[i][j]=0.0;
			for(k=0;k<ortho_size;k++){
				orthoinv[i][j] += get_lambdal(k,i)*get_lambdal(k,j);
			}
		}
	}
	
	int modes_size = ortho_size;
	create_modes();
	
	double *results_m =  malloc( sizeof(double)*2);
	results_m[0] = 0.0;
	results_m[1] = 0.0;
	
	time1 = csecond();
	
	double c1,c2,c3,c4;
	
	printf("[%d] Starting\n", rank);
	
	double **modearray = (double **)create_array(models,modes_size);
	
	calculate_modes_file(filename, models, modearray);

	double *modes_tmp = (double *)create_vector(ortho_size);
	
	for (n=0;n<models;n++){
// 		model = n;
		if(n==0)model=6;
		if(n==1)model=2;
		if(n==2)model=16;
		if(n==3)model=18;
		if(n==4)model=3;
		if(n==5)model=14;
		if(n==6)model=1;
		if(n==7)model=5;
		if(n==8)model=7;
		if(n==9)model=8;
		if(n==10)model=9;
		
		for(i=0;i<modes_size;i++){
			modes_tmp[i] = 0.0;
			for(j=0;j<modes_size;j++){
				modes_tmp[i] += orthoinv[i][j]*modearray[n][j]*deltaphi*deltaphi;
			}
		}
		
		printf("Modes Q\n");
		
		for(i=0;i<modes_size;i++){
			results_m[0] = (double)i;
			results_m[1] = modes_tmp[i];
			update_modes(results_m);
			printf("[%d]\t%d\t%e\n", n, i, results_m[1]);
		}
		
		output_modes(directory);
		
		printf("Modes R\n");

		for(i=0;i<modes_size;i++){
			x=0.0;
			for(j=0;j<modes_size;j++){
				x += get_orthol(j,i)*modes_tmp[j];
			}
			printf("[%d]\t%d\t%e\n", n, i, x);
		}
	
	}
	
	MPI_Finalize();
}

