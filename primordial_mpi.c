#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "global.h"


/*
1. set up timing variables
2. read in command line arguments for directories and l values
3. variable definitions
4. zero triangle and cell list to stop odd behaviour
5. read in Bessel and transfer data from file
6. create sparse grid of points over the triangle
7. calculate rough estimate of integral for input to area calculations
8. start adaptive algorithm for triangles
9. print result
*/

int main( int argc, char *argv[] ){

	double pi = 3.141592653589793;
	
// **1**

	double time1, time2, time3, time4, time5, duration;


	if (argc < 3 || argc > 6) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile model (p1) (p2) (p3)\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	
	model = atoi(argv[2]);
	shape_params params;
	
	if(argc>3){
		params.a1 = atoi(argv[3]);
	}else{
		params.a1 = 0;
	}
	
	if(argc>4){
		params.a2 = atof(argv[4]);
	}else{
		params.a2 = 0e0;
	}
	
	if(argc>5){
		params.a3 = atof(argv[5]);
	}else{
		params.a3 = 0e0;
	}
	
	int i,j,k,l,n;
	
	if(model==44)loadnicola();

// mpi vars
	int rank, nproc, lnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

// 	char lambdafile[100] = "/home/cosmos/tmp/jf334/Test2/lambda_143_DX8";
// 	char gammafile[100] = "/home/cosmos/tmp/jf334/Test2/gamma_143_DX8";
// 	char txtfile[100] = "/home/cosmos/users/jf334/Test2/michele.txt";

	set_terms_prim();
	set_terms_late();
	
	load_bessel();
	load_transfer();
	
	int size = get_terms_prim();
	
	double kmax = get_kmax();
	double x1,x2,x3;
	double s1,s2;
	
	double xvec[alpha_points+1];

	for(i=0;i<alpha_points+1;i++){
		xvec[i] = kmax*(double)i/(double)alpha_points;
	}

	
	create_basis_prim(alpha_points+1,kmax,xvec);
	int p1,p2,p3;
	p1 = params.a1;
	p2 = params.a2;
	p3 = params.a3;
	read_eigen(p1,p2,p3);
	
	for(i=0;i<alpha_points+1;i++){
		x1 = kmax*(double)i/alpha_points;
		s1 = shape3(x1,x1,x1,params);
		s2 = 0e0;
		for(n=0;n<size;n++){
			s2 += get_eigen(n)*pijk(n,i,i,i);
		}
		printf("%d\t%e\t%e\n",i,s1,s2);
	}
	
	for(i=0;i<alpha_points+1;i+=10){
		x1 = kmax*(double)i/alpha_points;
		for(j=i;j<alpha_points+1;j+=10){
			x2 = kmax*(double)j/alpha_points;
			for(k=j;k<min(i+j+1,alpha_points+1);k+=10){
				x3 = kmax*(double)k/alpha_points;
				s1 = shape3(x1,x2,x3,params);
				s2 = 0e0;
				for(n=0;n<size;n++){
					s2 += get_eigen(n)*pijk(n,i,j,k);
				}
				printf("%d\t%d\t%d\t%e\t%e\n",i,j,k,s1,s2);
			}
		}
	}
	
	// create_basis_prim(alpha_points+1,kmax,xvec);
	//
	// read_ortho();
	// read_lambda();
/*	
	double **lambda = create_array(size,size);
	double **gamma = create_array(size,size);
	
	double *results =  malloc( sizeof(double)*3);
	
	array_read(&size2, lambdafile, &lambda[0][0]);
	array_read(&size2, gammafile, &gamma[0][0]);

// 	double* list = malloc(sizeof(int)*4*MAXLINES);
// 	int* size_tmp = malloc(sizeof(int));
	
// 	load_txt_dbl(txtfile,4,list,size_tmp);
// 	printf("%d\n",*size_tmp);
	
// 	for (n=0;n<*size_tmp;n+=4) {
// 		i = (int)list[n];
// 		j = (int)list[n+1];
// 		lambda[i][j] = list[n+2];
// 		gamma[i][j] = list[n+3];
// 		printf("%d\t%d\t%e\t%e\n",i,j,lambda[i][j],gamma[i][j]);
// 	}

	gsl_matrix *cholesky = gsl_matrix_alloc(size,size);
	gsl_matrix *matrixin = gsl_matrix_alloc(size,size);
	gsl_eigen_symm_workspace *workspace = gsl_eigen_symm_alloc(size);
	gsl_vector *eigen = gsl_vector_alloc(size);
	double x;
	
	create_orthol();

	gsl_matrix_set_zero(cholesky);
	gsl_matrix_set_zero(matrixin);
	gsl_vector_set_zero(eigen);

	for (i=0;i<size;i++) {
		for (j=0;j<size;j++) {
			x = gamma[i][j];
// 			x = lambda[i][j];
			gsl_matrix_set(cholesky,i,j,x);
			gsl_matrix_set(matrixin,i,j,x);
		}
	}
	
	gsl_eigen_symm(matrixin,eigen,workspace);

	for(i=0;i<size;i++){
		printf("%d\t%e\n",i,gsl_vector_get(eigen,i));
	}
	
	gsl_linalg_cholesky_decomp(cholesky);
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			results[0] = (double)i;
			results[1] = (double)j;
			results[2] = gsl_matrix_get(cholesky, i, j);
			update_orthol(results);
		}
	}
	
	output_orthol(directory);
	
	gsl_matrix *upper = gsl_matrix_alloc(size,size);
	gsl_matrix *inverse = gsl_matrix_alloc(size,size);
	gsl_matrix_set_zero(upper);

		
	for(i=0;i<size;i++){
		for(j=i;j<size;j++){
			x = gsl_matrix_get(cholesky, i, j);
			gsl_matrix_set(upper,i,j,x);
		}
	}
		
	int s;
	gsl_permutation *perm = gsl_permutation_calloc(size);

	gsl_linalg_LU_invert(upper,perm,inverse);
		
	create_lambdal();
		
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			results[0] = (double)i;
			results[1] = (double)j;
			results[2] = gsl_matrix_get(inverse, j, i);
			update_lambdal(results);
			printf("Lambda:\t%d\t%d\t%e\n", i, j, results[2]);
		}
	}
	output_lambdal(directory);
*/
	//printf("[%d] Process[%d]  \n", rank, getpid());
	//
	// int tri_size = initial_grid_size + 1;
	// double k1, k2, k3, k4, x, x1, x2, y;
	// double height;

// 	double kmax = 2000.0/get_tau0();
// 	double kmin = 2.0/get_tau0();
/*	
// 	double kmax = 1.0;
// 	double kmin = 0.0;
	int size = 301;
	int fnlsize = 201;
	double fnl_input = 100.0;
	double nu_input = 1.0;
	int nu_index;
	int loadflag=1;
	int mcloops=1;
	
	for ( i = 0; i < size; i++) {
		x1 = 1.5*(double)i/(size-1);
		if(x1==nu_input)nu_index = i;
	}
// 	printf("%d\n",nu_index);

	double **covariance = create_array(size,size);
	double **correlation = create_array(size,size);
	double **correlation2 = create_array(size,size);

	double *norm = create_vector(size);
	y = calculate_localnorm(kmin, kmax, 1.0, 1.0);
// 	printf("%e\n",norm);

	char filename[100] = "/home/cosmos/tmp/jf334/Test2/XingangCov";
	int total_size = size*size;
	if(loadflag==0){
		for ( i = 0; i < size; i++) {
			x1 = 1.5*(double)i/(size-1);
			for ( j = i; j < size; j++) {
				x2 = 1.5*(double)j/(size-1);
				covariance[i][j] = calculate_correlation(kmin, kmax, x1, x2)/y;
				covariance[j][i] = covariance[i][j];
			}
		}
		array_write(&total_size, filename, &covariance[0][0]);
	}else{
		array_read(&total_size, filename, &covariance[0][0]);
	}
	
	for ( i = 0; i < size; i++) {
		norm[i]=1.0/sqrt(covariance[i][i]);
// 		printf("%d\t%e\n",i,norm[i]*5.0);
	}
	
// 	fnl_input /= norm[nu_index];
// 	printf("%e\n",fnl_input);
	
	gsl_matrix *A = gsl_matrix_alloc(size,size);
	gsl_matrix *V = gsl_matrix_alloc(size,size);
	gsl_vector *S = gsl_vector_alloc(size);
	gsl_vector *W = gsl_vector_alloc(size);

	for ( i = 0; i < size; i++) {
		for ( j = 0; j < size; j++) {
			y = covariance[i][j]/sqrt(covariance[i][i]*covariance[j][j]);
			correlation[i][j] = y;
			gsl_matrix_set(A,i,j,y);
// 			printf("%e\t",correlation[i][j]);
		}
// 		printf("\n");
	}
	
// 	printf("\n");
 
// 	for (i=0;i<size;i++) {
// 		for (j=0;j<size;j++) {
// 			y=0.0;
// 			x1=0.0;
// 			x2=0.0;
// 			for (k=0;k<size;k++) {
// 				y += correlation[i][k]*correlation[j][k];
// 				x1 += correlation[i][k]*correlation[i][k];
// 				x2 += correlation[j][k]*correlation[j][k];
// 			}
// 			correlation2[i][j] = x2*x1-y*y;
// 			printf("%e\t",correlation2[i][j]);
// 		}
// 		printf("\n");
// 	}
	
// 	for ( i = 0; i < size; i++) {
// 		for ( j = 0; j < size; j++) {
// 			y = 0.0;
// 			x1 = 0.0;
// 			x2 = 0.0;
// 			for ( k = 0; k < size; k++) {
// 				y += correlation[i][k]*correlation[j][k];
// 				x1 += correlation[i][k]*correlation[i][k];
// 				x2 += correlation[j][k]*correlation[j][k];
// 			}
// 			correlation2[i][j] = y/sqrt(x1*x2);
// 			printf("%e\t",correlation2[i][j]);
// 		}
// 		printf("\n");
// 	}
// 	printf("\n");
	
	gsl_linalg_SV_decomp(A, V, S, W);
	
	double *eigenvalue = create_vector(size);
	double **eigenvector = create_array(size,size);
// 	printf("\n");
	
	for (i=0;i<size;i++) {
		eigenvalue[i] = gsl_vector_get(S,i);
// 		printf("%d\t%e\n",i,eigenvalue[i]);
	}
// 	printf("\n");
	
	for (i=0;i<size;i++) {
		for (j=0;j<size;j++) {
			eigenvector[i][j] = gsl_matrix_get(V,i,j);
// 			printf("%e\t",eigenvector[i][j]);
		}
// 		printf("\n");
	}
	
// 	printf("done eigenvector\n");
	double **CV = create_array(size,size);
	
	for (i=0;i<size;i++) {
		for (j=0;j<size;j++) {
			CV[i][j] = eigenvector[i][j]*sqrt(eigenvalue[j]);
// 			CV[i][j] = eigenvector[i][j];
			if(j<10)printf("%e\t",CV[i][j]);
		}
// 		printf("%e\n",CV[i][1]/CV[i][0]);
		printf("\n");
	}
// 	printf("\n");

// 	for (i=0;i<size;i++) {
// 		x1 = 2.0*(double)i/(size-1)-1.0;
// 		for (j=0;j<size;j++) {
// 			x2 = 2.0*(double)j/(size-1)-1.0;
// 			x=pow((CV[0][1]/CV[0][0]-x2/x1),2);
// 			nu_index=0;
// 			for (k=1;k<size;k++) {
// 				y=pow((CV[k][1]/CV[k][0]-x2/x1),2);
// 				if(y<x){
// 					x=y;
// 					nu_index=k;
// 					
// 				}
// 			}
// 			printf("%e\t",1.5*(double)nu_index/(size-1));
// 		}
// 		printf("\n");
// 	}

// 	for (i=0;i<size;i++) {
// 		for (j=0;j<size;j++) {
// 			x1 = 0;
// 			for (k=0;k<size;k++) {
// 				x1 += CV[i][k]*eigenvector[j][k];
// 			}
// 			printf("%e\t",x1);
// 		}
// 		printf("\n");
// 	}

// 	printf("done CV\n");
	double *data = create_vector(size);
	double *sim = create_vector(size);
	
	double fnl, error;
	int nu;

	
	int **result = create_iarray(fnlsize,size);
		
	for (i=0;i<fnlsize;i++) {
		for (j=0;j<size;j++) {
			result[i][j] = 0;
		}
	}
	
	init_rng(10*(rank+1));
	
// 	printf("starting\n");

	int loops=mcloops/nproc;
	int auxloop = fmod(mcloops,nproc);
	int start_loop = rank*loops;
	int end_loop = (rank+1)*loops-1;
	int result_total = fnlsize*size;

	if (auxloop != 0){
		if (rank < auxloop){
			start_loop = start_loop + rank;
			end_loop = end_loop + rank + 1;
		}else{
			start_loop = start_loop + auxloop;
			end_loop = end_loop + auxloop;
		}
	}
	
	for(n=start_loop;n<end_loop+1;n++){
		
// 		y=0.0;
// 		for (i=0;i<size;i++) {
// 			y+=CV[nu_index][i]*CV[nu_index][i];
// 		}
// 		y=sqrt(y);
// 		printf("%e\n",y);
		
		for (i=0;i<size;i++) {
			data[i] = 0.0;
			sim[i] = fnl_input*CV[nu_index][i] + 5.0*get_rng_gauss();
// 			printf("%d\t%es\t%e\n",i,y,CV[nu_index][i]);
		}
		
 		x1=0.0;
 		x2=0.0;
 		for (j=0;j<size;j++) {
 			x1 += sim[j]*CV[0][j];
 			x2 += CV[0][j]*CV[0][j];
 		}

 		y = x1/x2;

		l = lround(y);
		if(l>=0&&l<fnlsize)result[l][0] += 1;

 		x2 = 0.0;
 		for (j=0;j<size;j++) {
 			x1 = sim[j]-y*CV[0][j];
 			x2 += x1*x1;
 		}

 		error = x2;
 		fnl = y;
 		nu = 0;


 		for (i=1;i<size;i++) {
 			x1=0.0;
 			x2=0.0;
 			for (j=0;j<size;j++) {
 				x1 += sim[j]*CV[i][j];
 				x2 += CV[i][j]*CV[i][j];
 			}
 			y = x1/x2;

			l = lround(y);
			if(l>=0&&l<fnlsize)result[l][i] += 1;

 			x2 = 0.0;
 			for (j=0;j<size;j++) {
 				x1 = sim[j]-y*CV[i][j];
 				x2 += x1*x1;
 			}
 			if(x2<error){
 				error = x2;
 				fnl = y;
 				nu = i;
 			}
 		}

		
// 		for (i=size-1;i>=0;i--) {
// 			for (j=size-1;j>=0;j--) {
// 				data[i] += sim[j]*eigenvector[i][j];
// 			}
// 			l = lround(data[i]);
// 			if(l>=0&&l<fnlsize)result[l][i] += 1;
// 		}
// 		printf("\n");
//
// 		for (i=0;i<size;i++) {
// 			printf("%d\t%e\n",i,data[i]);
// 		}
// 		printf("\n");

// 		x1=0.0;
// 		x2=0.0;
// 		for (j=0;j<size;j++) {
// 			x1 += data[j]*correlation[0][j];
// 			x2 += correlation[0][j]*correlation[0][j];
// 		}
// 		y = x1/x2;
// 
// 		x2 = 0.0;
// 		for (j=0;j<size;j++) {
// 			x1 = data[j]-y*correlation[0][j];
// 			x2 += x1*x1;
// 		}
// 
// 		error = x2;
// 		fnl = y;
// 		nu = 0;
// 
// 		for (i=1;i<size;i++) {
// 			x1=0.0;
// 			x2=0.0;
// 			for (j=0;j<size;j++) {
// 				x1 += data[j]*correlation[i][j];
// 				x2 += correlation[i][j]*correlation[i][j];
// 			}
// 			y = x1/x2;
// 
// 			x2 = 0.0;
// 			for (j=0;j<size;j++) {
// 				x1 = data[j]-y*correlation[i][j];
// 				x2 += x1*x1;
// 			}
// 			if(x2<error){
// 				error = x2;
// 				fnl = y;
// 				nu = i;
// 			}
// 		}
		
// 		l = lround(fnl*norm[nu]);
// 		if(l>=0&&l<fnlsize)result[l][nu] += 1;
		
	}
	
	int* result_send = (int *)malloc( result_total * sizeof(int) );
	int* result_recv = (int *)malloc( result_total * sizeof(int) );
	
	MPI_Barrier(MPI_COMM_WORLD);

	n=0;
	for (i=0;i<fnlsize;i++) {
		for (j=0;j<size;j++) {
			result_send[n] = result[i][j];
			n++;
		}
	}
	
        MPI_Reduce(result_send,result_recv,result_total,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		n=0;
		for (i=0;i<fnlsize;i++) {
			for (j=0;j<size;j++) {
				result[i][j] = result_recv[n];
				n++;
			}
		}
	}
	
	if(rank==0){
		for (i=0;i<fnlsize;i++) {
			for (j=0;j<size;j++) {
				printf("%d\t",result[i][j]);
			}
			printf("\n");
		}
	}




	
	height = slice(200, 300, 250);


*/
/*	
// 	double kmax = 1.0;
// 	double kmin = 0.0;
	int terms = size;
	size = 151;
	int fnlsize = 201;
// 	double fnl_input = 100.0;
// 	double nu_input = 1.0;
	int nu_index;
	int loadflag=1;
	int mcloops=2000000000;
	double **modelist = create_array(size,terms);
	
	for ( i = 0; i < size; i++) {
		x1 = 1.5*(double)i/(size-1);
		kstar = x1;
		read_modes(directory);
		for ( j = 0; j < size; j++) {
			modelist[i][j] = get_modes(j);
		}
// 		if(x1==nu_input)nu_index = i;
	}
// 	printf("%d\n",nu_index);

	double **covariance = create_array(size,size);
	double **correlation = create_array(size,size);
	double **correlation2 = create_array(size,size);
	double enorm = 0e0;
	double *norm = create_vector(size);
// 	y = calculate_localnorm(kmin, kmax, 1.0, 1.0);
// 	printf("%e\n",norm);

// 	char filename[100] = "/home/cosmos/tmp/jf334/Test2/XingangCov";
// 	int total_size = size*size;
// 	if(loadflag==0){
		for ( i = 0; i < size; i++) {
// 			x1 = 1.5*(double)i/(size-1);
			for ( j = i; j < size; j++) {
				x1 = 0e0;
				for ( n = 0; n < terms; n++) {
					x1 += modelist[i][n]*modelist[j][n];
				}
// 				x2 = 1.5*(double)j/(size-1);
				covariance[i][j] = x1;
				covariance[j][i] = covariance[i][j];
			}
		}
// 		array_write(&total_size, filename, &covariance[0][0]);
// 	}else{
// 		array_read(&total_size, filename, &covariance[0][0]);
// 	}
	
	for ( i = 0; i < size; i++) {
		norm[i]=1.0/sqrt(covariance[i][i]);
// 		printf("%d\t%e\n",i,norm[i]);
	}
	
// 	fnl_input /= norm[nu_index];
// 	printf("%e\n",fnl_input);
	
	gsl_matrix *A = gsl_matrix_alloc(size,size);
	gsl_matrix *V = gsl_matrix_alloc(size,size);
	gsl_vector *S = gsl_vector_alloc(size);
	gsl_eigen_symmv_workspace *W = gsl_eigen_symmv_alloc(size);

	for ( i = 0; i < size; i++) {
		for ( j = 0; j < size; j++) {
			y = covariance[i][j]/sqrt(covariance[i][i]*covariance[j][j]);
			correlation[i][j] = y;
			gsl_matrix_set(A,i,j,y);
// 			printf("%e\t",correlation[i][j]);
		}
// 		printf("\n");
	}
	
// 	printf("\n");
 
// 	for (i=0;i<size;i++) {
// 		for (j=0;j<size;j++) {
// 			y=0.0;
// 			x1=0.0;
// 			x2=0.0;
// 			for (k=0;k<size;k++) {
// 				y += correlation[i][k]*correlation[j][k];
// 				x1 += correlation[i][k]*correlation[i][k];
// 				x2 += correlation[j][k]*correlation[j][k];
// 			}
// 			correlation2[i][j] = x2*x1-y*y;
// 			printf("%e\t",correlation2[i][j]);
// 		}
// 		printf("\n");
// 	}
	
// 	for ( i = 0; i < size; i++) {
// 		for ( j = 0; j < size; j++) {
// 			y = 0.0;
// 			x1 = 0.0;
// 			x2 = 0.0;
// 			for ( k = 0; k < size; k++) {
// 				y += correlation[i][k]*correlation[j][k];
// 				x1 += correlation[i][k]*correlation[i][k];
// 				x2 += correlation[j][k]*correlation[j][k];
// 			}
// 			correlation2[i][j] = y/sqrt(x1*x2);
// 			printf("%e\t",correlation2[i][j]);
// 		}
// 		printf("\n");
// 	}
// 	printf("\n");
	
// 	gsl_linalg_SV_decomp(A, V, S, W);
	gsl_eigen_symmv(A,S,V,W);
	gsl_eigen_symmv_sort(S,V,GSL_EIGEN_SORT_VAL_DESC);
	
	double *eigenvalue = create_vector(size);
	double **eigenvector = create_array(size,size);
// 	printf("\n");
	
	for (i=0;i<size;i++) {
		eigenvalue[i] = gsl_vector_get(S,i);
		enorm += eigenvalue[i]*eigenvalue[i];
// 		printf("%d\t%e\n",i,eigenvalue[i]);
	}
// 	printf("\n");
	
	for (i=0;i<size;i++) {
		for (j=0;j<size;j++) {
			eigenvector[i][j] = gsl_matrix_get(V,i,j);
// 			printf("%e\t",eigenvector[i][j]);
		}
// 		printf("\n");
	}
	int size2 = 5;
// 	printf("done eigenvector\n");
	double **CV = create_array(size,size2);
	
	for (i=0;i<size;i++) {
		for (j=0;j<size2;j++) {
			CV[i][j] = eigenvector[i][j]*sqrt(eigenvalue[j]);
// 			CV[i][j] = eigenvector[i][j];
// 			printf("%e\t",CV[i][j]);
		}
// 		printf("%e\n",CV[i][1]/CV[i][0]);
// 		printf("\n");
	}
*/
// 	printf("\n");

// 	for (i=0;i<size;i++) {
// 		x1 = 2.0*(double)i/(size-1)-1.0;
// 		for (j=0;j<size;j++) {
// 			x2 = 2.0*(double)j/(size-1)-1.0;
// 			x=pow((CV[0][1]/CV[0][0]-x2/x1),2);
// 			nu_index=0;
// 			for (k=1;k<size;k++) {
// 				y=pow((CV[k][1]/CV[k][0]-x2/x1),2);
// 				if(y<x){
// 					x=y;
// 					nu_index=k;
// 					
// 				}
// 			}
// 			printf("%e\t",1.5*(double)nu_index/(size-1));
// 		}
// 		printf("\n");
// 	}

// 	for (i=0;i<size;i++) {
// 		for (j=0;j<size;j++) {
// 			x1 = 0;
// 			for (k=0;k<size;k++) {
// 				x1 += CV[i][k]*eigenvector[j][k];
// 			}
// 			printf("%e\t",x1);
// 		}
// 		printf("\n");
// 	}

// 	printf("done CV\n");
// 	double *data = create_vector(size2);
// 	double *sim = create_vector(size2);
// 	
// 	double fnl, error;
// 	int nu;
// 
// 	
// 	int **result = create_iarray(fnlsize,size);
// 		
// 	for (i=0;i<fnlsize;i++) {
// 		for (j=0;j<size;j++) {
// 			result[i][j] = 0;
// 		}
// 	}
// 	
// 	init_rng(10*(rank+1));
// 	
// // 	printf("starting\n");
// 
// 	int loops=mcloops/nproc;
// 	int auxloop = fmod(mcloops,nproc);
// 	int start_loop = rank*loops;
// 	int end_loop = (rank+1)*loops-1;
// 	int result_total = fnlsize*size;
// 
// 	if (auxloop != 0){
// 		if (rank < auxloop){
// 			start_loop = start_loop + rank;
// 			end_loop = end_loop + rank + 1;
// 		}else{
// 			start_loop = start_loop + auxloop;
// 			end_loop = end_loop + auxloop;
// 		}
// 	}
// 	
// 	for(n=start_loop;n<end_loop+1;n++){
// 		
// // 		y=0.0;
// // 		for (i=0;i<size;i++) {
// // 			y+=CV[nu_index][i]*CV[nu_index][i];
// // 		}
// // 		y=sqrt(y);
// // 		printf("%e\n",y);
// 		
// 		for (i=0;i<size2;i++) {
// 			data[i] = 0.0;
// // 			sim[i] = (4.79e0/norm[150])*CV[150][i] + sqrt(6e0*0.72)*get_rng_gauss();
// 			sim[i] = (4.79e0/norm[150])*CV[150][i] + sqrt(6e0*0.72)*get_rng_gauss()/enorm;
// // 			printf("%d\t%es\t%e\n",i,y,CV[nu_index][i]);
// 		}
// 		
//  		x1=0.0;
//  		x2=0.0;
//  		for (j=0;j<size2;j++) {
//  			x1 += sim[j]*CV[0][j];
//  			x2 += CV[0][j]*CV[0][j];
//  		}
// 
//  		y = x1/x2;
// 
// // 		l = lround(y)+100;
// // 		if(l>=0&&l<fnlsize)result[l][0] += 1;
// 
//  		x2 = 0.0;
//  		for (j=0;j<size2;j++) {
//  			x1 = sim[j]-y*CV[0][j];
//  			x2 += x1*x1;
//  		}
// 
//  		error = x2;
//  		fnl = y;
//  		nu = 0;
// //  		printf("%e\t%e\t%d\n",error,fnl,nu);
// 
// 
//  		for (i=1;i<size;i++) {
//  			x1=0.0;
//  			x2=0.0;
//  			for (j=0;j<size2;j++) {
//  				x1 += sim[j]*CV[i][j];
//  				x2 += CV[i][j]*CV[i][j];
//  			}
//  			y = x1/x2;
// // 			if(i==150)printf("%e\n",y);
// // 			l = lround(y)+100;
// // 			if(l>=0&&l<fnlsize)result[l][i] += 1;
// 
//  			x2 = 0.0;
//  			for (j=0;j<size2;j++) {
//  				x1 = sim[j]-y*CV[i][j];
//  				x2 += x1*x1;
//  			}
//  			if(x2<error){
//  				error = x2;
//  				fnl = y;
//  				nu = i;
// //  				printf("%e\t%e\t%d\n",error,fnl,nu);
//  			}
//  		}
// 
// // 		data = 0e0
// // 		for (i=size-1;i>=0;i--) {
// // 			for (j=size-1;j>=0;j--) {
// // 				data[i] += sim[j]*eigenvector[i][j];
// // 			}
// // 			l = lround(data[i])+100;
// // 			if(l>=0&&l<fnlsize)result[l][i] += 1;
// // 		}
// // 		printf("\n");
// // 
// // 		for (i=0;i<size;i++) {
// // 			printf("%d\t%e\n",i,data[i]);
// // 		}
// // 		printf("\n");
// // 
// // 		x1=0.0;
// // 		x2=0.0;
// // 		for (j=0;j<size;j++) {
// // 			x1 += data[j]*correlation[0][j];
// // 			x2 += correlation[0][j]*correlation[0][j];
// // 		}
// // 		y = x1/x2;
// // 
// // 		x2 = 0.0;
// // 		for (j=0;j<size;j++) {
// // 			x1 = data[j]-y*correlation[0][j];
// // 			x2 += x1*x1;
// // 		}
// // 
// // 		error = x2;
// // 		fnl = y;
// // 		nu = 0;
// // 
// // 		for (i=1;i<size;i++) {
// // 			x1=0.0;
// // 			x2=0.0;
// // 			for (j=0;j<size;j++) {
// // 				x1 += data[j]*correlation[i][j];
// // 				x2 += correlation[i][j]*correlation[i][j];
// // 			}
// // 			y = x1/x2;
// // 
// // 			x2 = 0.0;
// // 			for (j=0;j<size;j++) {
// // 				x1 = data[j]-y*correlation[i][j];
// // 				x2 += x1*x1;
// // 			}
// // 			if(x2<error){
// // 				error = x2;
// // 				fnl = y;
// // 				nu = i;
// // 			}
// // 		}
// // 		printf("%d\t%e\n",nu,fnl*norm[nu]);
// 		l = lround(fnl*norm[nu])+100;
// 		if(l>=0&&l<fnlsize)result[l][nu] += 1;
// 		
// 	}
// 	
// 	int* result_send = (int *)malloc( result_total * sizeof(int) );
// 	int* result_recv = (int *)malloc( result_total * sizeof(int) );
// 	
// 	MPI_Barrier(MPI_COMM_WORLD);
// 
// 	n=0;
// 	for (i=0;i<fnlsize;i++) {
// 		for (j=0;j<size;j++) {
// 			result_send[n] = result[i][j];
// 			n++;
// 		}
// 	}
// 	
//         MPI_Reduce(result_send,result_recv,result_total,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
// 	
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	
// 	if(rank==0){
// 		n=0;
// 		for (i=0;i<fnlsize;i++) {
// 			for (j=0;j<size;j++) {
// 				result[i][j] = result_recv[n];
// 				n++;
// 			}
// 		}
// 	}
// 	
// // 	printf("\n");
// 	if(rank==0){
// 		for (i=0;i<fnlsize;i++) {
// // 			printf("%d",i-100);
// 			for (j=0;j<size;j++) {
// 				printf("\t%d",result[i][j]);
// 			}
// 			printf("\n");
// 		}
// 	}

// 	printf("\n");

// 	double test = shape(1.0, 1.0, 1.0);
//  	printf("Normalisation: %e\n",test);
 	
// 	test = correlation_prim();
//  	printf("Correlation: %e\n",test);

// 	for ( i = 0; i < tri_size; i++) {
// 		x = (double)i / (tri_size - 1);
// 		for ( j = 0; j < tri_size - i; j++) {
// 			y = (double)(i+2*j - tri_size + 1) / (tri_size - 1);
// 			k1 = (1-x);
// 			k2 = (1+x+y)/2;
// 			k3 = (1+x-y)/2;
// // 			height = pow(k1*k2*k3,2)*xingang(0.0,k1, k2, k3);
// 			height = slice(k1, k2, k3);
//  			printf("%e\t%e\t%e\n",x,y,height);
// 		}
// 	}
// 	printf("[%d] Done.\n", rank);

// 	FILE *outfile;
// 	double k1,k2,k3,height,x,y;
//
// 	if (rank == 0){
//
// 		for ( i = 1; i < 200; i++) {
// 			k1 = get_kmax()*i/150e0;
// 			height = scale(k1);
// 			printf("%e\t%e\n",k1,height);
// 		}
//
// 		char modelnum[2] = "";
// 		sprintf(modelnum,"_%d",model);
//
// 		char filename[100] = "";
// 		strcat(filename, "/home/cosmos/users/jf334/Matlab/Data/Tri_Shape");
// 		strcat(filename, modelnum);
//
// 		printf("%s\n", filename);
// 		outfile = fopen(filename, "w");
//
// 		double scale = get_kmax()/2.0;
//
// 		for ( i = 0; i < 101; i++) {
// 			x = (double)i / 1e2;
// 			for ( j = 0; j < 101 - i; j++) {
// 				y = (double)(i+2e0*j - 1e2) / 1e2;
// 				k1 = scale*(1e0-x);
// 				k2 = scale*(1e0+x+y)/2e0;
// 				k3 = scale*(1e0+x-y)/2e0;
// // 				if(k1>k2+k3||k2>k1+k3||k3>k1+k2){
// // 					height=0.0;
// // 				}else{
// 					height = slice(k1, k2, k3);
// // 				}
// 	 			fprintf(outfile,"%e\t%e\t%e\n",x,y,height);
// 			}
// 		}
//
// 		fclose(outfile);

// 		char modelnum[2] = "";
// 		sprintf(modelnum,"_%d",model);
// 
// 		char filename[200] = "/home/cosmos/ccc-cam/jf334/Cosmos2/MapsDX9/beta_DDX9_1500_7_601_smica_1e4.txt";
// 		double* values = malloc( sizeof(double)*2*MAXLINES);
// 		double modesR[size];
// 		double modes[size];
// 		int* loadsize = malloc( sizeof(int));
// 		load_txt_dbl(filename, 2, values, loadsize);
// 		j=0;
// 		for(i=0;i<*loadsize;i++){
// 			j++;
// 			modesR[i] = values[j++];
// 		}
// 		
// 		for(i=0;i<size;i++){
// 			modes[i] = 0e0;
// 			for(j=i;j<size;j++){
// 				modes[i] += get_lambda(i,j)*modesR[j];
// 			}
// 		}
// 
// 		double result;
// 		double k1,k2,k3;
// 		
// 		char filename2[200] = "/home/cosmos/ccc-cam/jf334/Matlab/Data/beta_prim_1500_7_601_smica_1e4.txt";
// 		
// 		outfile = fopen(filename2, "w");
// 		
// 		
// 		for (i=0;i<alpha_points+1;i++){
// 			k1 = xvec[i];
// 			for (j=i;j<alpha_points+1;j++){
// 				k2 = xvec[j];
// 				for (k=j;k<alpha_points+1;k++){
// 					k3 = xvec[k];
// 					if(k<=i+j){
// 						result = 0e0;
// 						for(n=0;n<get_terms_prim();n++){
// 							result += modes[n]*pijk(n,i,j,k);
// 						}
// 						fprintf(outfile,"%d\t%d\t%d\t%e\n",i,j,k,result);
// 					}
// 				}
// 			}
// 		}
// 		
// 		fclose(outfile);
// 
	// }
// 	double a,b,c;
// 	double step=2.0/(double)initial_grid_size;
// 	for ( i = 0; i < tri_size; i++) {
// 		for ( j = 0; j < tri_size; j++) {
// 			for ( k = 0; k < tri_size; k++) {
// 				a=-0.5+i*step;
// 				b=-1.0+j*step;
// 				c=-1.0+k*step;
// 				k1=2.0*(1.0+a+sqrt(2)*b+sqrt(6)*c)/3.0;
// 				k2=2.0*(1.0+a+sqrt(2)*b-sqrt(6)*c)/3.0;
// 				k3=2.0*(1.0+a-2.0*sqrt(2)*b)/3.0;
// 				k4=2.0*(1.0-a);
// 				if(0.0<k1 && k1<=k2+k3+k4 && 0.0<k2 && k2<=k1+k3+k4 && 0.0<k3 && k3<=k2+k1+k4 && 0.0<k4 && k4<=k2+k3+k1){
// 					height = shape_tri(k1, k2, k3, k4);
//  					printf("%d\t%d\t%d\t%e\n",i,j,k,height);
//  				} else {
// 
//  				}
//  			}
// 		}
	// }

	MPI_Finalize();
	
	// End of code
	return 0;
}
