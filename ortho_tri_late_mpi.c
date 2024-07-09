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


 	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);

// mpi vars
	int rank, nproc, onext, mnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	set_terms_tri_prim();
	set_terms_tri_late();

	int i,j,k,l,m,n,r,s;
	int l_size = eflag_lmax+1;
	int *l_values = create_ivector(l_size);
	for(i=0;i<l_size;i++){
		l_values[i] = i;
	}
	
	double lmax = (double)l_values[l_size-1];
	init_lmax((int)lmax);
	
	long int l_size_long = l_size;
	
	create_cl(l_size);
	create_beam(l_size);
	create_noise(l_size);
	create_t_wgt(l_size);
	create_lens(l_size);
	load_cl(l_size, l_values);
	load_BN(l_size, l_values);
	load_TL(l_size, l_values);
	load_lens(l_size, l_values);
	int xsize = (int)lmax+1;
// 	printf("%e\n",get_cl(1));
// 	for(i=0;i<l_size;i++){
// 		printf("%d\t%e\n",l_values[i],get_cl(i)*get_beam(i)*get_beam(i)+get_noise(i));
// 	}
	
	double result=0;
	double x;

	double xvec[xsize];
	for(i=0;i<xsize;i++){
// 		xvec[i] = (double)i / lmax;
		xvec[i] = (double)i;
	}
// 	
	double xmax = (double)lmax;
	
	create_basis_tri_late(xsize, xmax, xvec);
	
	int ortho_size = get_terms_late();
	
	if(rank==0){
		printf("lmax: %d\n", (int)lmax);
		printf("pmax: %d\n", get_pmax_late());
		printf("terms: %d\n", get_terms_late());
		for(n=0;n<ortho_size;n++){
			find_perm_tri_late(n,&i,&j,&k,&l);
			printf("%d\t%d\t%d\t%d\t%d\n",n,i,j,k,l);
		}
	}
		
// 	if(rank==0){
// 		for(i=0;i<xsize;i++){
// 			printf("%e",xvec[i]);
// 			for(j=0;j<get_pmax_late()+1;j++){
// 				printf("\t%e",get_basis_late(i,j));
// 			}
// 			printf("\n");
// 		}
// 	}
	
	int loops;
	int auxloop ;
	int start_loop;
	int end_loop;
/*
 	gsl_matrix *cholesky = gsl_matrix_alloc(ortho_size,ortho_size);
	gsl_matrix *matrixin = gsl_matrix_alloc(ortho_size,ortho_size);
	gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(ortho_size);
	gsl_vector *eigen = gsl_vector_alloc(ortho_size);
	gsl_matrix *eigenv = gsl_matrix_alloc(ortho_size,ortho_size);
	double x1
*/	
		
	long int ortho_total = (ortho_size)*(ortho_size+1)/2;
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

	MPI_Barrier(MPI_COMM_WORLD);
	
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
			
		init_orthol_tri_glint();
		
		time1 = csecond();
		
		for(n=start_loop;n<end_loop+1;n++){
			
			i = ortho_pairs[n][0];
			j = ortho_pairs[n][1];
			ortho[i][j] = calculate_orthol_tri(i,j);
			
			duration = csecond() - time1;
			printf("[%d] QQ\t%d\t%d\t%d\t%e\t%e\n", rank, n, i, j, ortho[i][j], duration);
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
// 					printf("QQ\t%d\t%d\t%e\n",i,j,ortho[i][j]);
				}
// 				printf("\n");
			}
			
			
			double **dummy = (double **)create_array(ortho_size,ortho_size);
			for (i=0;i<ortho_size;i++) {
				for (j=0;j<ortho_size;j++) {
					dummy[i][j] = ortho[i][j] / sqrt(ortho[i][i]*ortho[j][j]);
				}
			}
		
			char suffix1[5],suffix2[5],suffix3[5];
			suffix1[0] = '\0';
			sprintf(suffix1, "%d", (int)lmax);
			suffix2[0] = '\0';
			sprintf(suffix2, "%d", eflag_order_late);
			suffix3[0] = '\0';
			sprintf(suffix3, "%d", ortho_size);
			char filename[100] = "/home/cosmos/tmp/jf334/DX9/master_orthol_tri_";
			strcat(filename, suffix1);
			strcat(filename, "_");
			strcat(filename, suffix2);
			strcat(filename, "_");
			strcat(filename, suffix3);
			strcat(filename, ".unf");

			int big_size = ortho_size*ortho_size;
			array_write(&big_size, filename, &dummy[0][0]);
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
// 				if(rank==0)printf("3D\t%d\t%d\n",r,s);
			}else{
				ortho_pairs_GL[n][0] = r;
				ortho_pairs_GL[n][1] = s;
				n++;
// 				if(rank==0)printf("GL\t%d\t%d\n",r,s);
			}
		}
		
		for (i=0;i<ortho_size;i++) {
			for (j=0;j<ortho_size;j++) {
				ortho[i][j] = 0.0;
			}
		}
		
// 		MPI_Barrier(MPI_COMM_WORLD);
		
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
			
// 			ortho[i][j] = calculate_orthol_tri_3D(i,j);
			
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
			
		init_orthol_tri_glint();
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		time1 = csecond();
		for(n=start_loop;n<end_loop+1;n++){
			
			i = ortho_pairs_GL[n][0];
			j = ortho_pairs_GL[n][1];
			
			ortho[i][j] = calculate_orthol_tri(i,j);
			
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
			
// 			char filename[100] = "/home/cosmos/tmp/jf334/DX9/master_gamma_3_1000.unf";
// 			int big_size = ortho_size*ortho_size;
// 			array_write(&big_size, filename, &ortho[0][0]);
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
		
	double *results_o =  malloc( sizeof(double)*3);
	results_o[0] = 0.0;
	results_o[1] = 0.0;
	results_o[2] = 0.0;
		
	create_orthol_tri();
	create_lambdal_tri();
			
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
// 				x1 = ortho[i][j];
				gsl_matrix_set(cholesky,i,j,x1);
				gsl_matrix_set(matrixin,i,j,x1);
// 				printf("%d\t%d\t%e\t%e\n",i,j,x1,ortho[i][j]);
			}
		}
		
		gsl_eigen_symmv(matrixin,eigen,eigenv,workspace);
		gsl_eigen_symmv_sort(eigen,eigenv,GSL_EIGEN_SORT_VAL_DESC);

		for(i=0;i<ortho_size;i++){
			printf("%d\t%e\n",i,gsl_vector_get(eigen,i));
		}
			
// 		for(i=0;i<ortho_size;i++){
// 			for(j=0;j<ortho_size;j++){
// 				printf("%d\t%d\t%e\n",i,j,gsl_matrix_get(eigenv,i,j));
// 			}
// 		}
			

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
					update_orthol_tri(results_o);
				}
			}
			
			for(i=0;i<ortho_size;i++){
				for(j=0;j<ortho_size;j++){
					results_o[0] = (double)i;
					results_o[1] = (double)j;
		 			results_o[2] = lambda[i][j];
	 				update_lambdal_tri(results_o);
// 					if(results_o[2]!=0)printf("Lambda:\t%d\t%d\t%e\n", i, j, results_o[2]);
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
// 	 					printf("Choleski:%d\t%d\t%e\n",i,j,results_o[2]);
	 				} else {
	 					results_o[2] = 0;
	 				}
					update_orthol_tri(results_o);
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
					update_lambdal_tri(results_o);
// 					if(results_o[2]!=0)printf("Lambda:\t%d\t%d\t%e\n", i, j, results_o[2]);
				}
			}
		
		}
		output_lambdal_tri();
		output_orthol_tri();
	}
 	
 	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Finalize();
}
