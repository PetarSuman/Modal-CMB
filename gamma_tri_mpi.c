#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
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
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	
	set_terms_tri_prim();
	set_terms_tri_late();

	int i,j,k,r,s,t,l,m,n,u;
	
	read_beta_tri();
	int terms = get_terms_prim();
	int lsize = get_bt_lsize();
	int *lvec = create_ivector(lsize);
	get_bt_lvec(lvec);
	
	int lmax = lvec[lsize-1];
	int xsize = lmax+1;
	init_lmax(lmax);
	create_cl(lsize);
	create_beam(lsize);
	create_noise(lsize);
	create_t_wgt(lsize);
	create_lens(lsize);
	load_cl(lsize, lvec);
	load_BN(lsize, lvec);
	load_TL(lsize, lvec);
	load_lens(lsize, lvec);
		
	double sum1,sum2,sum3,sum4;
	double x1,x2,x3,x4;
	
	int ortho_size = terms;
	
	if(rank==0){
		printf("lmax: %d\n", (int)lmax);
		printf("prim pmax: %d\n", get_pmax_prim());
		printf("late pmax: %d\n", get_pmax_late());
		printf("prim terms: %d\n", get_terms_prim());
		printf("late terms: %d\n", get_terms_late());
		printf("xsize: %d\n",get_bt_xsize());
		for(n=0;n<ortho_size;n++){
			find_perm_tri_prim(n,&i,&j,&k,&l);
			find_perm_tri_late(n,&r,&s,&t,&u);
			printf("%d\t(%d,%d,%d,%d)\t(%d,%d,%d,%d)\n",n,i,j,k,l,r,s,t,u);
		}
// 		printf("lvec\n");
// 		for(l=0;l<lsize;l++){
// 			printf("%d\t%e\t%e\t%e\t%e\n",l,get_cl(l),get_beam(l),get_noise(l),get_beam(l)*get_beam(l)*get_cl(l)+get_noise(l));
// 		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	double xvec[lmax+1];
	for(i=0;i<xsize;i++){
// 		xvec[i] = (double)i / lmax;
		xvec[i] = (double)i;
	}
// 	
	double xmax = (double)lmax;
	create_basis_tri_late(xsize, xmax, xvec);
	
	int loops;
	int auxloop ;
	int start_loop;
	int end_loop;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	read_orthol_tri();
	read_lambdal_tri();
	double **orthoinv = (double **)create_array(ortho_size,ortho_size);

	for(i=0;i<ortho_size;i++){
		for(j=0;j<ortho_size;j++){
			orthoinv[i][j]=0.0;
			for(k=0;k<ortho_size;k++){
				orthoinv[i][j] += get_lambdal_tri(k,i)*get_lambdal_tri(k,j);
			}
// 			printf("%d\t%d\t%e\n",i,j,orthoinv[i][j]);
		}
	}
	
	create_gamma_tri();
	int gamma_size = terms;
	int gamma_total = terms*terms;
	double **gamma = (double **)create_array(gamma_size,gamma_size);

	if(eflag_order_late==6 || rflag_do3D == 0){
		
		int gamma_4D = (lmax-1)*terms;
		int gamma_pairs_4D[gamma_4D][2];
		m=0;
		
		for(n=0;n<gamma_size;n++){
			for(i=2;i<lmax+1;i++){
				gamma_pairs_4D[m][0] = n;
				gamma_pairs_4D[m][1] = i;
				m++;
			}
		}
	
		loops=gamma_4D/nproc;
		auxloop = fmod(gamma_4D,nproc);
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
 		
		printf("4D mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
		MPI_Barrier(MPI_COMM_WORLD);
		
		double *mvec = (double *)create_vector(gamma_size);
		
		time1 = csecond();
		for(n=start_loop;n<end_loop+1;n++){
			
			i = gamma_pairs_4D[n][0];
			j = gamma_pairs_4D[n][1];
// 			calculate_gamma_tri_4D(i,j,mvec);
			
			for (r=0;r<gamma_size;r++) {
				gamma[r][i] += mvec[r];
			}
			
			duration = csecond() - time1;
			printf("[%d]\t4D QQ\t%d\t%d\t%d\t%e\n",rank,n, i, j, duration);
			time1 = csecond();
			
		}
		printf("[%d] finished\n", rank);
		
	}else{
	
		int gamma_total = terms*terms;
		int gamma_pairs[gamma_total][2];
		n=0;
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				gamma_pairs[n][0] = i;
				gamma_pairs[n][1] = j;
				n++;
			}
		}
		
		init_gamma_tri_glint();
		
		loops=gamma_total/nproc;
		auxloop = fmod(gamma_total,nproc);
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
	 
		printf("Gl mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
		
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				gamma[i][j] = 0.0;
			}
		}
		
		time1 = csecond();
		for(n=start_loop;n<end_loop+1;n++){
			
			i = gamma_pairs[n][0];
			j = gamma_pairs[n][1];
			
			gamma[i][j] = calculate_gamma_tri(i,j);
			
			duration = csecond() - time1;
			printf("[%d]\t%d\t%d\t%e\t%e\n", rank, i, j, gamma[i][j], duration);
			time1 = csecond();
			
		} // end of MPI loop
		printf("[%d] finished\n", rank);
	}

	double* gamma_send = (double *)malloc( gamma_total * sizeof(double) );
	double* gamma_recv = (double *)malloc( gamma_total * sizeof(double) );
	
	MPI_Barrier(MPI_COMM_WORLD);

	n=0;
	for (i=0;i<gamma_size;i++) {
		for (j=0;j<gamma_size;j++) {
			gamma_send[n] = gamma[i][j];
			n++;
		}
	}
	
        MPI_Reduce(gamma_send,gamma_recv,gamma_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		n=0;
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				gamma[i][j] = gamma_recv[n];
				n++;
				printf("%e\t",gamma[i][j]);
			}
			printf("\n");
		}
	}
	
	if(rank==0){
		
		double **gamma_temp = (double **)create_array(gamma_size,gamma_size);
		
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
		 		gamma_temp[i][j] = 0.0;
				for (n=0;n<gamma_size;n++) {
		 			gamma_temp[i][j] +=orthoinv[i][n]*gamma[n][j];
		 		}
			}
		}
		
		printf("\n");
			
		double *results_g =  malloc( sizeof(double)*3);
		results_g[0] = 0.0;
		results_g[1] = 0.0;
		results_g[2] = 0.0;

		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				results_g[0] = (double)i;
				results_g[1] = (double)j;
				results_g[2] = gamma_temp[i][j];
				update_gamma_tri(results_g);
// 				printf("Gamma:\t%d\t%d\t%e\n", i, j, results_g[2]);
				printf("%e\t",results_g[2]);
			}
			printf("\n");
		}
		output_gamma_tri();
		printf("[%d] Written Gamma\n", rank);
		
		printf("\n");

		for (n=0;n<gamma_size;n++) {
			for (i=0;i<gamma_size;i++) {
				gamma[i][n] = 0.0;
				for (j=i;j<gamma_size;j++) {
					gamma[i][n] += gamma_temp[j][n]*get_orthol_tri(j,i);
				}
				printf("%e\t",gamma[i][n]);
			}
			printf("\n");
		}
		
		printf("\n");
		
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				x1=x2=x3=0.0;
				for (n=0;n<gamma_size;n++) {
					x1 += gamma[n][i]*gamma[n][i];
					x2 += gamma[n][j]*gamma[n][j];
					x3 += gamma[n][i]*gamma[n][j];
				}
				x4 = x3 / sqrt(x1*x2);
				printf("%e\t",x4);
			}
			printf("\n");
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}