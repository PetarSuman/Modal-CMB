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
#include <hdf5.h>


#define FILE_NAME "/nfs/st01/hpc-gr-epss/ps792/PlanckV3/C/gamma.h5"
#define DATASET_NAME "gamma_matrix"

/*
Late time basis Q_n(l) isn't (necessarily) orthonormal. The interior product
between different basis elements is given by the \gamma matrix.

ortho_late_mpi.c calculates          < Q_n(l) Q_m(l) >  :=  \gamma_{nm} (l) 

Its job submission script is in the ympijob files.
*/

int main( int argc, char *argv[] ){

//	Initialisation

	double time1, time2, time3, time4, time5, duration;

 	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}

//	MPI variables
	int rank, nproc, onext, mnext;

//	MPI Initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

//	load relevant parameters from the parameter.ini file or alike
	if (rank == 0) printf("Initialise\n");
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);

	int i,j,k,l,m,n,r,s;
	set_terms_late();
	if (rank == 0) printf("Terms set late...\n");
	
	double x;

	double Tmax = (double)eflag_T_lmax;
	double Emax = (double)eflag_E_lmax;
	double Tmin = (double)eflag_T_lmin;
	double Emin = (double)eflag_E_lmin;
	int Tsize = Tmax-Tmin+1;
	int Esize = Emax-Emin+1;
	double Tvec[Tsize];
	double Evec[Esize];
	for(i=0;i<Tsize;i++)Tvec[i] = (double)i+Tmin;
	for(i=0;i<Esize;i++)Evec[i] = (double)i+Emin;
	
	int pmax_T = get_pmax_late_T();
	int pmax_E = get_pmax_late_E();
	int psize_T, psize_E;
	psize_T = get_pmax_late_T()+1;
	if(do_polarisation==1) psize_E = get_pmax_late_E()+1;
	
	int ortho_size = get_terms_late();
	
	if(rank==0){
		printf("lmax T: %d\n", eflag_T_lmax);
		printf("lmin T: %d\n", eflag_T_lmin);
		printf("pmax T: %d\n", get_pmax_late_T());
		if(do_polarisation==1){
			printf("lmax E: %d\n", eflag_E_lmax);
			printf("lmin E: %d\n", eflag_E_lmin);
			printf("pmax E: %d\n", get_pmax_late_E());
		}
		printf("terms: %d\n", get_terms_late());
		if(do_polarisation==1){
			for(n=0;n<ortho_size;n++){
				find_perm_late_TTT(n,&i,&j,&k);
				printf("TTT\t%d\t%d\t%d\t%d\n",n,i,j,k);
			}
			for(n=0;n<ortho_size;n++){
				find_perm_late_TTE(n,&i,&j,&k);
				printf("TTE\t%d\t%d\t%d\t%d\n",n,i,j,k);
			}
			for(n=0;n<ortho_size;n++){
				find_perm_late_TEE(n,&i,&j,&k);
				printf("TEE\t%d\t%d\t%d\t%d\n",n,i,j,k);
			}
			for(n=0;n<ortho_size;n++){
				find_perm_late_EEE(n,&i,&j,&k);
				printf("EEE\t%d\t%d\t%d\t%d\n",n,i,j,k);
			}
		}else{
			for(n=0;n<ortho_size;n++){
				find_perm_late_TTT(n,&i,&j,&k);
				printf("TTT\t%d\t%d\t%d\t%d\n",n,i,j,k);
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) printf("Want to create basis late...\n");
	create_basis_late(Tsize, Esize, Tmin, Tmax, Emin, Emax, Tvec, Evec);

	if (rank == 0) printf("Basis successfully created.\n");
	
	double* result;

	if(do_polarisation==1){
		result = (double *)malloc( 4 * sizeof(double) );
	}else{
		result = (double *)malloc( 1 * sizeof(double) );
	}

	//printf("Check 1\n");
	
	int loops;
	int auxloop ;
	int start_loop;
	int end_loop;
		
	int ortho_total = (ortho_size)*(ortho_size+1)/2;
	//int ortho_pairs[ortho_total][2];
	int **ortho_pairs = (int **)create_array(ortho_total, 2);
	double **ortho_TTT = (double **)create_array(ortho_size,ortho_size);
	double **ortho_TTE;
	double **ortho_TEE;
	double **ortho_EEE;

	//printf("Check 2\n");

	//if (rank == 0) printf("Ortho_XXX created.\n");
	 
	if(do_polarisation==1){
		ortho_TTE = (double **)create_array(ortho_size,ortho_size);
		ortho_TEE = (double **)create_array(ortho_size,ortho_size);
		ortho_EEE = (double **)create_array(ortho_size,ortho_size);
	}
	//printf("Check 3\n");
	n=0;
	for (i=0;i<ortho_size;i++) {
		for (j=i;j<ortho_size;j++) {
			ortho_pairs[n][0] = i;
			ortho_pairs[n][1] = j;
			n++;
	// 		printf("[%d]\t%d\t%d\t%d\n",rank,n,i,j);
		}
	}
	//printf("Check 4\n");
	if(eflag_order_late!=6){
			
		init_orthol_glint();
		int Mij_total;
		if(do_polarisation==1){
			Mij_total = psize_T*psize_T+psize_E*psize_E;
		}else{
			Mij_total = psize_T*psize_T;
		}

		int Mij_pairs[Mij_total][3];
		n=0;
		for (i=0;i<psize_T;i++) {
			for (j=0;j<psize_T;j++) {
				Mij_pairs[n][0] = 0;
				Mij_pairs[n][1] = i;
				Mij_pairs[n][2] = j;
				n++;
			}
		}	
		if(do_polarisation==1){
			for (i=0;i<psize_E;i++) {
				for (j=0;j<psize_E;j++) {
					Mij_pairs[n][0] = 1;
					Mij_pairs[n][1] = i;
					Mij_pairs[n][2] = j;
					n++;
				}
			}
		}
		printf("Check 5\n");
		loops=Mij_total/nproc;
		auxloop = fmod(Mij_total,nproc);
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
		
		time1 = csecond();
		for(n=start_loop;n<end_loop+1;n++){
			
			i = Mij_pairs[n][1];
			j = Mij_pairs[n][2];
			if(Mij_pairs[n][0]==0){
				calculate_orthoMij_T(i,j);
			}else{
				calculate_orthoMij_E(i,j);
			}
			if((n-start_loop+1)%10==0){
				duration = csecond() - time1;
				//printf("[%d] done 10 loops M_ij in %e\n",rank,duration);
				time1 = csecond();
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		sync_orthoMij();
		MPI_Barrier(MPI_COMM_WORLD);
		
		for (i=0;i<ortho_size;i++) {
			for (j=0;j<ortho_size;j++) {
				ortho_TTT[i][j] = 0.0;
				if(do_polarisation==1){
					ortho_TTE[i][j] = 0.0;
					ortho_TEE[i][j] = 0.0;
					ortho_EEE[i][j] = 0.0;
				}
			}
		}
		//printf("Check 6\n");
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
 	
		//printf("mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
	
		time1 = csecond();
		int count = 0;
		for(n=start_loop;n<end_loop+1;n++){
			
			i = ortho_pairs[n][0];
			j = ortho_pairs[n][1];
			
			calculate_orthol(i,j,result);
			ortho_TTT[i][j] = result[0];
			//printf("ortho_TTT[%d][%d] = %f\n", i, j, result[0]);
			if(do_polarisation==1){
				ortho_TTE[i][j] = result[1];
				ortho_TEE[i][j] = result[2];
				ortho_EEE[i][j] = result[3];
			}


			
			if((n-start_loop+1)%100==0){
				duration = csecond() - time1;
				if(do_polarisation==1){
					//printf("QQ\t%d\t%d\t%e\t%e\t%e\t%e\t%e\n", i, j, ortho_TTT[i][j], ortho_TTE[i][j], ortho_TEE[i][j], ortho_EEE[i][j], duration);
				}else{
					//printf("QQ\t%d\t%d\t%e\t%e\n", i, j, ortho_TTT[i][j], duration);
				}
				time1 = csecond();
			}
			
		} // end of MPI loop
		//printf("[%d] finished\n", rank);
		
		double* ortho_TTT_send = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TTT_recv = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TTE_send = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TTE_recv = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TEE_send = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TEE_recv = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_EEE_send = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_EEE_recv = (double *)malloc( ortho_total * sizeof(double) );
		
		MPI_Barrier(MPI_COMM_WORLD);
	
		n=0;
		for (i=0;i<ortho_size;i++) {
			for (j=i;j<ortho_size;j++) {
				ortho_TTT_send[n] = ortho_TTT[i][j];
				n++;
				printf("ortho_TTT[%d][%d] = %f\n", i, j, ortho_TTT[i][j]);
			}
		}
		if(do_polarisation==1){
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
					ortho_TTE_send[n] = ortho_TTE[i][j];
					n++;
					printf("ortho_TTE[%d][%d] = %f\n", i, j, ortho_TTE[i][j]);
				}
			}
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
					ortho_TEE_send[n] = ortho_TEE[i][j];
					n++;
					printf("ortho_TEE[%d][%d] = %f\n", i, j, ortho_TEE[i][j]);
				}
			}
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
					ortho_EEE_send[n] = ortho_EEE[i][j];
					n++;
					printf("ortho_EEE[%d][%d] = %f\n", i, j, ortho_EEE[i][j]);
				}
			}
		}
		
		MPI_Reduce(ortho_TTT_send,ortho_TTT_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
		if(do_polarisation==1){
			MPI_Reduce(ortho_TTE_send,ortho_TTE_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(ortho_TEE_send,ortho_TEE_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(ortho_EEE_send,ortho_EEE_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0){
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
				
					ortho_TTT[i][j] = ortho_TTT_recv[n];
					if(i!=j)ortho_TTT[j][i] = ortho_TTT[i][j];
					
					if(do_polarisation==1){
						ortho_TTE[i][j] = ortho_TTE_recv[n];
						if(i!=j)ortho_TTE[j][i] = ortho_TTE[i][j];
						ortho_TEE[i][j] = ortho_TEE_recv[n];
						if(i!=j)ortho_TEE[j][i] = ortho_TEE[i][j];
						ortho_EEE[i][j] = ortho_EEE_recv[n];
						if(i!=j)ortho_EEE[j][i] = ortho_EEE[i][j];
					}
					
					n++;
				}
			}
			
			double **dummy_TTT = (double **)create_array(ortho_size,ortho_size);
			double **dummy_TTE = (double **)create_array(ortho_size,ortho_size);
			double **dummy_TEE = (double **)create_array(ortho_size,ortho_size);
			double **dummy_EEE = (double **)create_array(ortho_size,ortho_size);
			for (i=0;i<ortho_size;i++) {
				for (j=0;j<ortho_size;j++) {
					dummy_TTT[i][j] = ortho_TTT[i][j] / sqrt(ortho_TTT[i][i]*ortho_TTT[j][j]);
				}
			}
			
			if(do_polarisation==1){
				for (i=0;i<ortho_size;i++) {
					for (j=0;j<ortho_size;j++) {
						dummy_TTE[i][j] = ortho_TTE[i][j] / sqrt(ortho_TTE[i][i]*ortho_TTE[j][j]);
						dummy_TEE[i][j] = ortho_TEE[i][j] / sqrt(ortho_TEE[i][i]*ortho_TEE[j][j]);
						dummy_EEE[i][j] = ortho_EEE[i][j] / sqrt(ortho_EEE[i][i]*ortho_EEE[j][j]);
					}
				}
			}
			/*
			// It is for saving a large Gamma for reduction to remove degenerate modes when creating a new basis
			// If you are using existing bases, then you don’t need it

			char suffix1[5],suffix2[5],suffix3[5],suffix4[5];
			suffix1[0] = '\0';
			sprintf(suffix1, "%d", (int)eflag_T_lmax);
			suffix2[0] = '\0';
			sprintf(suffix2, "%d", (int)eflag_E_lmax);
			suffix3[0] = '\0';
			sprintf(suffix3, "%d", eflag_order_late);
			suffix4[0] = '\0';
			sprintf(suffix4, "%d", ortho_size);

			char filename1[200] = "/fast/space/projects/planck/ps792.private/Master/master_orthol_TTT_";
			strcat(filename1, suffix1);
			strcat(filename1, "_");
			strcat(filename1, suffix2);
			strcat(filename1, "_");
			strcat(filename1, suffix3);
			strcat(filename1, "_");
			strcat(filename1, suffix4);
			strcat(filename1, ".unf");
			char filename2[200] = "/fast/space/projects/planck/ps792.private/Master/master_orthol_TTE_";
			strcat(filename2, suffix1);
			strcat(filename2, "_");
			strcat(filename2, suffix2);
			strcat(filename2, "_");
			strcat(filename2, suffix3);
			strcat(filename2, "_");
			strcat(filename2, suffix4);
			strcat(filename2, ".unf");
			char filename3[200] = "/fast/space/projects/planck/ps792.private/Master/master_orthol_TEE_";
			strcat(filename3, suffix1);
			strcat(filename3, "_");
			strcat(filename3, suffix2);
			strcat(filename3, "_");
			strcat(filename3, suffix3);
			strcat(filename3, "_");
			strcat(filename3, suffix4);
			strcat(filename3, ".unf");
			char filename4[200] = "/fast/space/projects/planck/ps792.private/Master/master_orthol_EEE_";
			strcat(filename4, suffix1);
			strcat(filename4, "_");
			strcat(filename4, suffix2);
			strcat(filename4, "_");
			strcat(filename4, suffix3);
			strcat(filename4, "_");
			strcat(filename4, suffix4);
			strcat(filename4, ".unf");

			int big_size = ortho_size*ortho_size;
			array_write(&big_size, filename1, &dummy_TTT[0][0]);
			array_write(&big_size, filename2, &dummy_TTE[0][0]);
			array_write(&big_size, filename3, &dummy_TEE[0][0]);
			array_write(&big_size, filename4, &dummy_EEE[0][0]);
			*/
		}
	}
	else{
		
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
		// 		if(rank==0)printf("3D\t%d\t%d\n",r,s);
			}else{
				ortho_pairs_GL[n][0] = r;
				ortho_pairs_GL[n][1] = s;
				n++;
		// 		if(rank==0)printf("GL\t%d\t%d\n",r,s);
			}
		}
		
		for (i=0;i<ortho_size;i++) {
			for (j=0;j<ortho_size;j++) {
				ortho_TTT[i][j] = 0.0;
				if(do_polarisation==1){
					ortho_TTE[i][j] = 0.0;
					ortho_TEE[i][j] = 0.0;
					ortho_EEE[i][j] = 0.0;
				}
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
		
		//printf("3D mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
		MPI_Barrier(MPI_COMM_WORLD);
		
		time1 = csecond();
		for(n=start_loop;n<end_loop+1;n++){
			
			i = ortho_pairs_3D[n][0];
			j = ortho_pairs_3D[n][1];
			
			calculate_orthol_3D(i,j,result);
			ortho_TTT[i][j] = result[0];
			if(do_polarisation==1){
				ortho_TTE[i][j] = result[1];
				ortho_TEE[i][j] = result[2];
				ortho_EEE[i][j] = result[3];
			}

			
			if((n-start_loop+1)%100==0){
				duration = csecond() - time1;
				if(do_polarisation==1){
					//printf("QQ\t%d\t%d\t%e\t%e\t%e\t%e\t%e\n", i, j, ortho_TTT[i][j], ortho_TTE[i][j], ortho_TEE[i][j], ortho_EEE[i][j], duration);
				}else{
					//printf("QQ\t%d\t%d\t%e\t%e\n", i, j, ortho_TTT[i][j], duration);
				}
				time1 = csecond();
			}
			
		} // end of MPI loop
		
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
 		
		//printf("GL mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
		MPI_Barrier(MPI_COMM_WORLD);
			
		init_orthol_glint();
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		time1 = csecond();
		for(n=start_loop;n<end_loop+1;n++){
			
			i = ortho_pairs_GL[n][0];
			j = ortho_pairs_GL[n][1];
			
			calculate_orthol(i,j,result);
			ortho_TTT[i][j] = result[0];
			if(do_polarisation==1){
				ortho_TTE[i][j] = result[1];
				ortho_TEE[i][j] = result[2];
				ortho_EEE[i][j] = result[3];
			}
			
			if((n-start_loop+1)%100==0){
				duration = csecond() - time1;
				if(do_polarisation==1){
					//printf("QQ\t%d\t%d\t%e\t%e\t%e\t%e\t%e\n", i, j, ortho_TTT[i][j], ortho_TTE[i][j], ortho_TEE[i][j], ortho_EEE[i][j], duration);
				}else{
					//printf("QQ\t%d\t%d\t%e\t%e\n", i, j, ortho_TTT[i][j], duration);
				}
				time1 = csecond();
			}
			
		} // end of MPI loop
		
		double* ortho_TTT_send = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TTT_recv = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TTE_send = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TTE_recv = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TEE_send = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_TEE_recv = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_EEE_send = (double *)malloc( ortho_total * sizeof(double) );
		double* ortho_EEE_recv = (double *)malloc( ortho_total * sizeof(double) );
		
		MPI_Barrier(MPI_COMM_WORLD);
	
		n=0;
		for (i=0;i<ortho_size;i++) {
			for (j=i;j<ortho_size;j++) {
				ortho_TTT_send[n] = ortho_TTT[i][j];
				n++;
			}
		}
		if(do_polarisation==1){
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
					ortho_TTE_send[n] = ortho_TTE[i][j];
					n++;
				}
			}
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
					ortho_TEE_send[n] = ortho_TEE[i][j];
					n++;
				}
			}
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
					ortho_EEE_send[n] = ortho_EEE[i][j];
					n++;
				}
			}
		}
		
		MPI_Reduce(ortho_TTT_send,ortho_TTT_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
		if(do_polarisation==1){
			MPI_Reduce(ortho_TTE_send,ortho_TTE_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(ortho_TEE_send,ortho_TEE_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(ortho_EEE_send,ortho_EEE_recv,ortho_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(rank==0){
			n=0;
			for (i=0;i<ortho_size;i++) {
				for (j=i;j<ortho_size;j++) {
				
					ortho_TTT[i][j] = ortho_TTT_recv[n];
					if(i!=j)ortho_TTT[j][i] = ortho_TTT[i][j];
					
					if(do_polarisation==1){
						ortho_TTE[i][j] = ortho_TTE_recv[n];
						if(i!=j)ortho_TTE[j][i] = ortho_TTE[i][j];
						ortho_TEE[i][j] = ortho_TEE_recv[n];
						if(i!=j)ortho_TEE[j][i] = ortho_TEE[i][j];
						ortho_EEE[i][j] = ortho_EEE_recv[n];
						if(i!=j)ortho_EEE[j][i] = ortho_EEE[i][j];
					}
					
					n++;
				}
			}
		}
	}
	
	/*if (rank == 0){
		//--------------- HDF5 OUTPUT CODE ----------------------------------------

		hid_t file_id, dataset_id, dataspace_id;
	  	int ROWS = ortho_size;
	  	int COLS = ortho_size;
  	  	hsize_t dims[2] = {ROWS, COLS};
  		int i, j;

	  	// Allocate memory for the 2D array
  		// Allocate memory for the flattened 1D array
  		double* data1d = (double*)malloc(ortho_size * ortho_size * sizeof(double));

  		// Flatten the 2D array into a 1D array
  		for (i = 0; i < ortho_size; i++) {
    		for (j = 0; j < ortho_size; j++) {
      			data1d[i * ortho_size + j] = ortho_TTT[i][j];
    		}
  		}

  		// Create a new HDF5 file
 		file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  		// Create the data space for the 2D array
  		dataspace_id = H5Screate_simple(2, dims, NULL);

  		// Create the dataset
  		dataset_id = H5Dcreate(file_id, DATASET_NAME, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  		// Write the flattened 1D array to the dataset
  		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);

  		// Close the HDF5 objects
  		H5Dclose(dataset_id);
  		H5Sclose(dataspace_id);
  		H5Fclose(file_id);

  		free(data1d);

		//---------------- END OF HDF5 CODE ---------------------------
	}*/


	double *results;
	if(do_polarisation==1){
		results =  malloc( sizeof(double)*6);
	}else{
		results =  malloc( sizeof(double)*3);
	}
	results[0] = 0.0;
	results[1] = 0.0;
	results[2] = 0.0;
	if(do_polarisation==1){
		results[3] = 0.0;
		results[4] = 0.0;
		results[5] = 0.0;
	}
		
	create_orthol();
	create_lambdal();
			
	if(rank==0){
		int npol=1;
		if(do_polarisation==1)npol=4;
		double ***ortho = (double ***)create_3Darray(npol,ortho_size,ortho_size);
		double ***lambda = (double ***)create_3Darray(npol,ortho_size,ortho_size);
		double ***lambdainv = (double ***)create_3Darray(npol,ortho_size,ortho_size);
		double norm[npol][ortho_size];
		
		for (i=0;i<ortho_size;i++) {
			for (j=0;j<ortho_size;j++) {
				ortho[0][i][j] = ortho_TTT[i][j];
				if(do_polarisation==1){
					ortho[1][i][j] = ortho_TTE[i][j];
					ortho[2][i][j] = ortho_TEE[i][j];
					ortho[3][i][j] = ortho_EEE[i][j];
				}

				//printf("ortho[0][%d][%d] = %f\n", i, j, ortho[0][i][j]);
			}
		}
	
		
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
			
		//printf("Normalisation\n");
		for (i=0;i<ortho_size;i++) {
			for (n=0;n<npol;n++) {
				norm[n][i] = sqrt(ortho[n][i][i]);
				//printf("Normalized ortho(%d)\t%e\n",i, norm[n][i]);
			}
		}
	
		for (n=0;n<npol;n++) {
			for (i=0;i<ortho_size;i++) {
				for (j=0;j<ortho_size;j++) {
					x1 = ortho[n][i][j]/(norm[n][i]*norm[n][j]);
					gsl_matrix_set(cholesky,i,j,x1);
					gsl_matrix_set(matrixin,i,j,x1);
				}
			}
			
			gsl_eigen_symmv(matrixin,eigen,eigenv,workspace);
			gsl_eigen_symmv_sort(eigen,eigenv,GSL_EIGEN_SORT_VAL_DESC);
	
			//printf("[%d]\tEigenvalues\n",n);
			/*for(i=0;i<ortho_size;i++){
				printf("%d\t%e\n",i,gsl_vector_get(eigen,i));
			}*/
				
			
			if(gflag_pca==0||gsl_vector_get(eigen,ortho_size-1)<0){
				
				int r_size;
				for(i=0;i<ortho_size && gsl_vector_get(eigen,i)>0.0;i++){
					r_size = i+1;
				}
				
				for(i=0;i<ortho_size;i++){
					for(j=0;j<ortho_size;j++){
						if(i<r_size){
							lambda[n][i][j] = gsl_matrix_get(eigenv,j,i)/(norm[n][j]*sqrt(gsl_vector_get(eigen,i)));
						}else{
							lambda[n][i][j] = 0.0;
						}
						if(j<r_size){
							lambdainv[n][i][j] = norm[n][i]*gsl_matrix_get(eigenv,i,j)*sqrt(gsl_vector_get(eigen,j));
						}else{
							lambdainv[n][i][j] = 0.0;
						}
					}
				}
				
			}else{
				
				gsl_linalg_cholesky_decomp(cholesky);
			
				for(i=0;i<ortho_size;i++){
					for(j=0;j<ortho_size;j++){
		 				if(i>=j){
		 					lambdainv[n][i][j] = norm[n][i]*gsl_matrix_get(cholesky, i, j);
		 				} else {
		 					lambdainv[n][i][j] = 0e0;
		 				}
						update_orthol(results);
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
						lambda[n][i][j] = gsl_matrix_get(inverse, j, i)/norm[n][j];
					}
				}
			
			}
		}
	
			
		for(i=0;i<ortho_size;i++){
			for(j=0;j<ortho_size;j++){
				results[0] = (double)i;
				results[1] = (double)j;
				results[2] = lambdainv[0][i][j];
				if(do_polarisation==1){
					results[3] = lambdainv[1][i][j];
					results[4] = lambdainv[2][i][j];
					results[5] = lambdainv[3][i][j];
				}
				update_orthol(results);
			}
		}
		
		for(i=0;i<ortho_size;i++){
			for(j=0;j<ortho_size;j++){
				results[0] = (double)i;
				results[1] = (double)j;
				results[2] = lambda[0][i][j];
				if(do_polarisation==1){
					results[3] = lambda[1][i][j];
					results[4] = lambda[2][i][j];
					results[5] = lambda[3][i][j];
				}
				update_lambdal(results);
			}
		}
		output_lambdal();
		output_orthol();
	}
	// if (rank == 0){

	// 	int lmax = get_lmax();
	// 	double xvec[alpha_points+1];

	// 	for(i=0;i<alpha_points+1;i++){
	// 		xvec[i] = lmax*(double)i/(double)alpha_points;
    // 	}

	// 	create_basis_late(alpha_points+1,kmax,xvec);


	// }
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Finalize();
}
