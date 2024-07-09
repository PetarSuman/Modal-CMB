#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>
#include "global.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_bessel.h>
#include <hdf5.h>


#define FILE_NAME "/nfs/st01/hpc-gr-epss/ps792/PlanckV3/C/corr_alphas_legendre.h5"
#define DATASET_NAME "corr_alpha"

//#define FILE_NAME "/nfs/st01/hpc-gr-epss/ps792/PlanckV3/C/bispectrum.h5"
//#define DATASET_NAME "primordial_bispectrum"

int main( int argc, char *argv[] ){
	
	double time1, time2, time3, time4, time5, duration, duration1, duration2;

	bool slice_late = false;
	bool slice_prim = false;
	bool do_corr_prim = false;
	bool parallel_corr = false;
	bool plot_bispectrum = false;
	bool print_shapes = false;
	bool brute_corr = false;

	if (argc < 3 || argc > 3) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile model\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);

	// mpi vars
	int rank, nproc, onext, enext;

	// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	set_terms_prim();
	//set_terms_late();

	int i,j,k,l,m,n;
	double result=0;
	load_bessel();
	load_transfer();
	if(model==44)loadnicola();
    int nthreads = omp_get_max_threads();
	double kmax = get_kmax();
	int ortho_size = get_terms_prim();

	double xvec[alpha_points+1];
	for(i=0;i<alpha_points+1;i++){
		xvec[i] = kmax*(double)i/(double)alpha_points;
	}
	create_basis_prim(alpha_points+1,kmax,xvec);
    read_ortho();
	read_lambda();

	
	double c1, c2,x;
	double corr_local, corr_equi, corr_ortho;
	double sum, sum1, sum2, sum3, sum4, sum5, sum6, k1, k2, k3;
	int p1,p2,p3,r,s,t;

	gsl_matrix *cholesky = gsl_matrix_alloc(ortho_size,ortho_size);

	int eigen_size = ortho_size;
	double *eigen = (double *)create_vector(eigen_size);
	//double *eigen_R_model = (double *)create_vector(eigen_size);
	//double *eigen_R_qsf = (double *)create_vector(eigen_size);
	//double *eigen_R_local = (double *)create_vector(eigen_size);
	//double *eigen_R_equi = (double *)create_vector(eigen_size);
	//double *eigen_R_ortho = (double *)create_vector(eigen_size);
	//double *eigen_R_const = (double *)create_vector(eigen_size);
	//model = 2;
	shape_params params;
	params.a1 = 30;
	params.a2 = 0e0;
	params.a3 = 0e0;
	
	model = atoi(argv[2]);
	create_eigen();
	//printf("test1\n");
	read_eigen(30, 0, 0);
	//printf("test2\n");
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < ortho_size; i++){
		x = get_eigen(i);
		printf("%d %e\n", i, x);
	}

	if (rank == 0){
		double c3;
		for(i=0;i<alpha_points+1;i+=1){
			k1 = xvec[i];
			c1 = 0e0;
			c3 = 0e0;
			for(n=0;n<ortho_size;n++){
				//find_perm_prim(n,&q1,&q2,&q3);
				//sum1 = get_basis_prim(i,q1)*get_basis_prim(j,q2)*get_basis_prim(k,q3);
				//sum2 = get_basis_prim(i,q2)*get_basis_prim(j,q3)*get_basis_prim(k,q1);
				//sum3 = get_basis_prim(i,q3)*get_basis_prim(j,q1)*get_basis_prim(k,q2);
				//sum4 = get_basis_prim(i,q3)*get_basis_prim(j,q2)*get_basis_prim(k,q1);
				//sum5 = get_basis_prim(i,q2)*get_basis_prim(j,q1)*get_basis_prim(k,q3);
				//sum6 = get_basis_prim(i,q1)*get_basis_prim(j,q3)*get_basis_prim(k,q2);
				//c1 += get_eigen(n) * (sum1+sum2+sum3+sum4+sum5+sum6)/6e0;
				c3 += get_eigen(n) * pijk(n, i, i, i);
			}
			c2 = shape3(k1,k1,k1,params);
			printf("%e\t%e\t%e\n", k1, c3, c2);
		}
	}


	if (parallel_corr){
		// Calculate the subrange for the current rank
		int total_values = (slim.p1_end - slim.p1_bgn) / slim.p1_stp + 1;
		int subrange_size = total_values / nproc;
		int remainder = total_values % nproc;
		int start_value = slim.p1_bgn + rank * subrange_size * slim.p1_stp;
		int end_value = start_value + subrange_size * slim.p1_stp - slim.p1_stp;

		// Ensure the last rank handles any remaining values
		if (rank == nproc - 1) {
			end_value += remainder * slim.p1_stp;
		}
		for (p1 = start_value; p1 <= end_value; p1 += slim.p1_stp){
			//printf("%d\n", p1);
			/*for(p2=slim.p2_bgn;p2<slim.p2_end+1;p2+=slim.p2_stp){
				for(p3=slim.p3_bgn;p3<slim.p3_end+1;p3+=slim.p3_stp){
					params.a1 = (double)p1;
					params.a2 = (double)p2;
					params.a3 = (double)p3;

					corr1 = 0.0;
					corr2 = 0.0;
					corr3 = 0.0;
					corr = 0.0;
					read_eigen(p1,p2,p3);
					for(i=0;i<alpha_points+1;i+=1){
						k1 = xvec[i];
						for(j=i;j<alpha_points+1;j+=1){
							k2 = xvec[j];
							for(k=j;k<alpha_points+1;k+=1){
								k3 = xvec[k];
								if((k1+k2>=k3) && (k1+k3>=k2) && (k2+k3>=k1)){
									c1 = 0e0;
									for(n=0;n<ortho_size;n++){
										c1 += get_eigen(n) * pijk(n,i,j,k);
									}
									x = volumes[i][j][k] * cube_size;
									c2 = shape3(k1,k2,k3,params);
									ksum = k1 + k2 + k3;
									if (ksum>1e-10) x = x/ksum;
									else x = 0e0;
									corr1 += c1 * c2 * x;
									corr2 += c1 * c1 * x;
									corr3 += c2 * c2 * x;
									//printf("%e %e %e %e %e\n", k1, k2, k3, c1, c2);
								}
							}
						}
					}
					printf("%e %e\n", corr2, corr3);
					corr2 = 1e0/sqrt(corr2);
					corr3 = 1e0/sqrt(corr3);
					corr = corr1 * corr2 * corr3;
					printf("%d\t%d\t%d\t%e\n", p1, p2, p3, corr);
				}
			}*/
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0){

			/*
			model = 2;	//equilateral
			read_modes(0,0,0);
			sum1 = 0.0;
			for (i = 0; i < eigen_size; i++){
				eigen_R_equi[i] = 0.0;
			}
			for (i = 0; i < eigen_size; i++){
				for (j = 0; j < eigen_size; j++){
					eigen_R_equi[i] += get_orthol_TTT(j,i)*get_modes_TTT(j);
				}
				sum1 += eigen_R_equi[i] * eigen_R_equi[i];
			}
			//printf("%e\n", sum1);
			//------------------------------------------------------------------------------
			model = 6;	//local
			read_modes(0,0,0);
			sum2 = 0.0;
			for (i = 0; i < eigen_size; i++){
				eigen_R_local[i] = 0.0;
			}
			for (i = 0; i < eigen_size; i++){
				for (j = 0; j < eigen_size; j++){
					eigen_R_local[i] += get_orthol_TTT(j,i)*get_modes_TTT(j);
				}
				sum2 += eigen_R_local[i] * eigen_R_local[i];
			}
			//printf("%e\n", sum2);
			//------------------------------------------------------------------------------
			model = 16;	//orthogonal
			read_modes(0,0,0);
			sum3 = 0.0;
			for (i = 0; i < eigen_size; i++){
				eigen_R_ortho[i] = 0.0;
			}
			for (i = 0; i < eigen_size; i++){
				for (j = 0; j < eigen_size; j++){
					eigen_R_ortho[i] += get_orthol_TTT(j,i)*get_modes_TTT(j);
				}
				sum3 += eigen_R_ortho[i] * eigen_R_ortho[i];
			}
			//printf("%e\n", sum3);
			//------------------------------------------------------------------------------
			model = 14;	//constant
			read_modes(0,0,0);
			sum4 = 0.0;
			for (i = 0; i < eigen_size; i++){
				eigen_R_const[i] = 0.0;
			}
			for (i = 0; i < eigen_size; i++){
				for (j = 0; j < eigen_size; j++){
					eigen_R_const[i] += get_orthol_TTT(j,i)*get_modes_TTT(j);
				}
				sum4 += eigen_R_const[i] * eigen_R_const[i];
			}
			//printf("%e\n", sum4);
			//------------------------------------------------------------------------------
			equi_local = 0.0;
			equi_ortho = 0.0;
			equi_const = 0.0;
			local_ortho = 0.0;
			local_const = 0.0;
			ortho_const = 0.0;
			for (i = 0; i < eigen_size; i++){
				equi_local += eigen_R_equi[i] * eigen_R_local[i];
				equi_ortho += eigen_R_equi[i] * eigen_R_ortho[i];
				equi_const += eigen_R_equi[i] * eigen_R_const[i];
				local_ortho += eigen_R_ortho[i] * eigen_R_local[i];
				local_const += eigen_R_local[i] * eigen_R_const[i];
				ortho_const += eigen_R_ortho[i] * eigen_R_const[i];

			}
			equi_local = equi_local / (sqrt(sum1) * sqrt(sum2));
			equi_ortho = equi_ortho / (sqrt(sum1) * sqrt(sum3));
			equi_const = equi_const / (sqrt(sum1) * sqrt(sum4));
			local_ortho = local_ortho / (sqrt(sum2) * sqrt(sum3));
			local_const = local_const / (sqrt(sum2) * sqrt(sum4));
			ortho_const = ortho_const / (sqrt(sum3) * sqrt(sum4));
			printf("********** LATE TIME **************\n");
			printf("%e %e %e\n", equi_local, equi_ortho, equi_const);
			printf("%e %e\n", local_ortho, local_const);
			printf("%e\n", ortho_const);
			*/

		printf("________ SLICE	PRIM ________\n");
		if (slice_prim){
			model = 61;
			double K = 0.5*kmax;
			p1 = 30;
			p2 = 0;
			p3 = 0;
			read_eigen(p1,p2,p3);
			printf("k1\tk2\tk3\tReconstructed\tTheory\n");
			for(i=0;i<alpha_points+1;i+=1){
				k1 = xvec[i];
				// printf("%d\t%e\n",i, get_basis_prim(i,1));
				for(j=0;j<alpha_points+1;j+=1){
					k2 = xvec[j];
					for(k=0;k<alpha_points+1;k+=1){
						k3 = xvec[k];
						if((k1+k2>=k3) && (k1+k3>=k2) && (k2+k3>=k1) && (k1 + k2 + k3 < 0.51*kmax) && (k1 + k2 + k3 > 0.49*kmax)){
							c1 = 0e0;
							for(n=0;n<ortho_size;n++){
								//find_perm_prim(n,&q1,&q2,&q3);
								//sum1 = get_basis_prim(i,q1)*get_basis_prim(j,q2)*get_basis_prim(k,q3);
								//sum2 = get_basis_prim(i,q2)*get_basis_prim(j,q3)*get_basis_prim(k,q1);
								//sum3 = get_basis_prim(i,q3)*get_basis_prim(j,q1)*get_basis_prim(k,q2);
								//sum4 = get_basis_prim(i,q3)*get_basis_prim(j,q2)*get_basis_prim(k,q1);
								//sum5 = get_basis_prim(i,q2)*get_basis_prim(j,q1)*get_basis_prim(k,q3);
								//sum6 = get_basis_prim(i,q1)*get_basis_prim(j,q3)*get_basis_prim(k,q2);
								//c1 += get_eigen(n) * (sum1+sum2+sum3+sum4+sum5+sum6)/6e0;
								c1 += get_eigen(n) * pijk(n, i, j, k);
							}
							c2 = shape3(k1,k2,k3,params);
							printf("%e\t%e\t%e\t%e\t%e\n", k1, k2, k3, c1, c2);
						}
					}
				}
			}
		}
		printf("________ SLICE	LATE ________\n");
		if (slice_late){
			model = 61;
			read_modes(80,0,0);
			double Tmax = (double)eflag_T_lmax;
			double Tmin = (double)eflag_T_lmin;
			int Tsize = Tmax-Tmin+1;  
			double Tvec[Tsize];
			for(i=0;i<Tsize;i++)Tvec[i] = (double)i+Tmin;

			double Emax = (double)eflag_E_lmax;
			double Emin = (double)eflag_E_lmin;
			int Esize = Emax-Emin+1;
			double Evec[Esize];
			for(i=0;i<Esize;i++)Evec[i] = (double)i+Emin;

			create_basis_late(Tsize, Esize, Tmin, Tmax, Emin, Emax, Tvec, Evec);
			load_cl();
			load_BN();
			load_lens();

			double *weight = (double *)create_vector(Tmax+1);
			int pmaxT = get_pmax_late_T();
			double **QTT = (double **)create_array(pmaxT+1,Tmax+1);
			for(r=0;r<Tmax+1;r++){
				weight[r] = sqrt(get_cl_TT(r)+get_noise_TT(r)/(get_beam_TT(r)*get_beam_TT(r)))* pow(2e0*r+1e0,-1e0/6e0);
				//printf("%d\t%e\t%e\t%e\t%e\n",r,get_beam_TT(r),get_cl_TT(r),get_noise_TT(r),weight[r]);
			}
			for(i=0;i<pmaxT+1;i++){
				for(r=Tmin;r<Tmax+1;r++){
					// QTT[i][r] = get_basis_late_T(r-Tmin,i);
					QTT[i][r] = get_basis_late_T(r-Tmin,i) * weight[r];
						// if(i>pmaxT-3)printf("%d\t%e\t%e\t%e\n",r,QTT[i][r],get_basis_late_T(r-Tmin,i),weight[r]);
				}
			}
			int l1,l2,l3,L;
			double x,x1,x2,x3,x4,x5,x6;
			printf("L\tl1\tl2\tl3\tB(l1,l2,l3)\n");
			for(r=0;r<eflag_T_lmax+1;r+=1){
				for(s=0;s<eflag_T_lmax+1;s+=1){
					for(t=0;t<eflag_T_lmax+1;t+=1){
						l1 = r;
						l2 = s;
						l3 = t;
						if(l1==0)l1=eflag_T_lmin;
						if(l2==0)l2=eflag_T_lmin;
						if(l3==0)l3=eflag_T_lmin;
						x=0e0;
						L = l1 + l2 + l3;
						if(l1+l2>=l3 && l2+l3>=l1 && l1+l3>=l2 && L%2==0 && (L == 500 || L == 1000 || L == 1500 || L==2000)){
							for(n=0;n<ortho_size;n++){
								find_perm_late_TTT(n,&i,&j,&k);
								//printf("i j k = %d %d %d\n", i, j, k);
								x1 = QTT[i][l1]*QTT[j][l2]*QTT[k][l3];
								x2 = QTT[j][l1]*QTT[k][l2]*QTT[i][l3];
								x3 = QTT[k][l1]*QTT[i][l2]*QTT[j][l3];
								x4 = QTT[k][l1]*QTT[j][l2]*QTT[i][l3];
								x5 = QTT[j][l1]*QTT[i][l2]*QTT[k][l3];
								x6 = QTT[i][l1]*QTT[k][l2]*QTT[j][l3];
								x += get_modes_TTT(n)*(x1+x2+x3+x4+x5+x6)/6e0;
							}
							printf("%d\t%d\t%d\t%d\t%e\n", L, l1, l2, l3, x);
						}
					}
				}
			}
		}
		printf("________ FINISHED CODE ________\n");
    }
	MPI_Finalize();
}