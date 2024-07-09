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
//	if (argc != 4) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);

	set_terms_tri_prim();
	set_terms_tri_late();

	int i,j,k,l,m,n,r,s,t,u;
	
	int rank, nproc;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	double x;
	
	int ortho_size = get_terms_late();
	int big_size = ortho_size*ortho_size;
	int lmax = eflag_lmax;
	
	if(rank==0){
		printf("pmax: %d\n", get_pmax_late());
		printf("terms: %d\n", get_terms_late());
// 		for(n=0;n<ortho_size;n++){
// 			find_perm_tri_late(n,&i,&j,&k,&l);
// 			printf("%d\t%d\t%d\t%d\t%d\n",n,i,j,k,l);
// 		}
	}

	
	double **ortho = (double **)create_array(ortho_size,ortho_size);
	
	
	char suffix1[5],suffix2[5],suffix3[5];
	suffix1[0] = '\0';
	sprintf(suffix1, "%d", lmax);
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
	
 	array_read(&big_size, filename, &ortho[0][0]);

	double norm[ortho_size];
		
	for (i=0;i<ortho_size;i++) {
		norm[i] = sqrt(ortho[i][i]);
	}
	
	int index, trueindex;
	double test;
		
	int* list = (int *)malloc( ortho_size * sizeof(int) );
	double* listout = (double *)malloc( ortho_size * sizeof(double) );
	double* listin = (double *)malloc( ortho_size * sizeof(double) );

	char filename2[100] = "/home/cosmos/ccc-cam/jf334/Cosmos2/Test1/faillist_tri_";
	strcat(filename2, suffix1);
	strcat(filename2, "_");
	strcat(filename2, suffix2);
	strcat(filename2, "_");
	strcat(filename2, suffix3);
	strcat(filename2, ".txt");
	int *faillist = malloc(sizeof(int)*ortho_size);
	int *failsize = malloc(sizeof(int));
	load_txt_int(filename2, 1, faillist, failsize);
	
	int list_size = ortho_size - *failsize;
	
	if(rank==0)printf("modes complete\n");
	n=0;
	for (i=0;i<ortho_size;i++) {
		m=0;
		for(j=0;j<*failsize;j++){
			if(faillist[j]==i){
				m=1;
				if(rank==0){
					find_perm_tri_late(i,&r,&s,&t,&u);
					printf("%d\t%d\t%d\t%d\t%d\n",i,r,s,t,u);
				}
				break;
			}
		}
		if(m==0){
			list[n] = i;
			n++;
		}
	}
	
	if(rank==0){
		printf("modes active\n");
		for(n=0;n<list_size;n++){
			s = list[n];
			find_perm_tri_late(s,&i,&j,&k,&l);
			printf("%d\t%d\t%d\t%d\t%d\n",s,i,j,k,l);
		}
		printf("starting reduction\n");
	}
	
	double **array = (double **)create_array(list_size,list_size);
	
	
	for (i=0;i<list_size;i++) {
		for (j=0;j<list_size;j++) {
			array[i][j] = ortho[list[i]][list[j]]/(norm[list[i]]*norm[list[j]]);
// 			if(rank==0)printf("%d\t%d\t%e\n",i,j,array[i][j]);
		}
	}
		
	if(rank==0)time1 = csecond();
	for (n=list_size-1;n>0;n--) {
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		gsl_matrix *matrix = gsl_matrix_alloc(n,n);
		gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(n);
		gsl_vector *eigen = gsl_vector_alloc(n);
		gsl_matrix *eigenv = gsl_matrix_alloc(n,n);
		
		
		for (i=0;i<n+1;i++){
			listout[i] = 0.0;
			listin[i] = 0.0;
		}
		
		int loops=(n+1)/nproc;
		int auxloop = fmod((n+1),nproc);
		int start_loop = rank*loops;
		int end_loop = (rank+1)*loops-1;
		
		if (auxloop != 0){
			if (rank < auxloop){
				start_loop = start_loop + rank;
				end_loop = end_loop + rank + 1;
			}else{
				start_loop = start_loop + auxloop;
				end_loop = end_loop + auxloop;
			}
		}
		
// 		printf("mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
// 		MPI_Barrier(MPI_COMM_WORLD);
		
		for(m=start_loop;m<end_loop+1;m++){
			for (i=0;i<m;i++) {
				for (j=i;j<m;j++) {
					x = array[i][j];
					gsl_matrix_set(matrix,i,j,x);
					gsl_matrix_set(matrix,j,i,x);
				}
				for (j=m+1;j<n+1;j++) {
					x = array[i][j];
					gsl_matrix_set(matrix,i,j-1,x);
					gsl_matrix_set(matrix,j-1,i,x);
				}
			}
			for (i=m+1;i<n+1;i++) {
				for (j=i;j<n+1;j++) {
					x = array[i][j];
					gsl_matrix_set(matrix,i-1,j-1,x);
					gsl_matrix_set(matrix,j-1,i-1,x);
				}
			}
				
			gsl_eigen_symmv(matrix,eigen,eigenv,workspace);
			gsl_eigen_symmv_sort(eigen,eigenv,GSL_EIGEN_SORT_VAL_DESC);
				
			listout[m] = gsl_vector_get(eigen,n-1);
// 			listout[m] = m;
		}
		
		MPI_Allreduce(listout,listin,n+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
		
		test = listin[0];
		index = 0;
		for (i=1;i<n+1;i++){
			x = listin[i];
			if(x>test){
				test = x;
				index = i;
			}
		}
		
		trueindex = list[index];
		
		for (i=index+1;i<n+1;i++) {
			list[i-1] = list[i];
		}
		
		for (i=0;i<n+1;i++) {
			for (j=index+1;j<n+1;j++) {
				array[i][j-1] = array[i][j];
			}
		}
		for (i=index+1;i<n+1;i++) {
			for (j=0;j<n+1;j++) {
				array[i-1][j] = array[i][j];
			}
		}
		
// 		s = -1;
// 		for(i=0;i<list_size;i++){
// 			if(list[i]!=-1){
// 				s++;
// 				if(s==index){
// 					trueindex = i;
// 					list[i] = -1;
// 					break;
// 				}
// 			}
// 		}
		
		find_perm_tri_late(trueindex,&i,&j,&k,&l);
		if(rank==0){
			duration = csecond() - time1;
			printf("%d\t%d\t%d\t%d\t%d\t%e\t%e\n",trueindex,i,j,k,l,test, duration);
			time1 = csecond();
		}
		gsl_matrix_free(matrix);
		gsl_eigen_symmv_free(workspace);
		gsl_vector_free(eigen);
		gsl_matrix_free(eigenv);
	}
	if(rank==0){
		n = list[0];
		find_perm_tri_late(n,&i,&j,&k,&l);
		printf("%d\t%d\t%d\t%d\t%d\t%e\t%e\n",n,i,j,k,l,1e0, 0e0);
	}

	MPI_Finalize();
}

