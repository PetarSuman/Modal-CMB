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
	set_terms_late();
	
	model = atoi(argv[2]);
	int i,j,k,l,m,n;
	double result=0;
	load_bessel();
	load_transfer();
	
	double kmax = get_kmax();
	if(model==44)loadnicola();
	
	double xvec[alpha_points+1];
	for(i=0;i<alpha_points+1;i++){
		xvec[i] = kmax*(double)i/(double)alpha_points;
	}
	double ***volumes = create_3Darray(alpha_points + 1, alpha_points + 1, alpha_points + 1);
    int nthreads = omp_get_max_threads();
	tetrapyd_limits* block;
	block = (tetrapyd_limits*) malloc(nthreads * sizeof(tetrapyd_limits));
	
	int ortho_size = get_terms_prim();
	if(rank==0){
		printf("model: %d\n",model);
		printf("kmax: %e\n", kmax);
		printf("pmax: %d\n", get_pmax_prim());
		/*for(n=0;n<ortho_size;n++){
			find_perm_prim(n,&i,&j,&k);
			printf("%d\t%d\t%d\t%d\n",n,i,j,k);
		}*/
	}
	if(slim.p1_end!=0){
		if(rank==0)printf("p1:\t%e\t%e\t%e\n",(double) slim.p1_bgn,(double) slim.p1_stp,(double)slim.p1_end);
	}else{
		slim.p1_bgn = 0;
		slim.p1_stp = 1;
		slim.p1_end = 0;
	}
	if(slim.p2_end!=0){
		if(rank==0)printf("p2:\t%e\t%e\t%e\n",(double)slim.p2_bgn,(double)slim.p2_stp,(double)slim.p2_end);
	}else{
		slim.p2_bgn = 0;
		slim.p2_stp = 1;
		slim.p2_end = 0;
	}
	if(slim.p3_end!=0){
		if(rank==0)printf("p3:\t%e\t%e\t%e\n",(double)slim.p3_bgn,(double)slim.p3_stp,(double)slim.p3_end);
	}else{
		slim.p3_bgn = 0;
		slim.p3_stp = 1;
		slim.p3_end = 0;
	}
	//printf("Check 1\n");
	create_basis_prim(alpha_points+1,kmax,xvec);
	MPI_Barrier(MPI_COMM_WORLD);
	//printf("Check 2\n");
	decompose_tetrapyd_prim(block, rank, nproc, nthreads, 0, alpha_points);
	//printf("Check 3\n");
	for(i=0;i<nthreads;i++){
		printf("[%d:%d]\t%d\t%d\t%d\t%d\n",rank,i,block[i].i_bgn,block[i].j_bgn,block[i].k_bgn,block[i].loops);
	}
	//printf("Check 4\n");
	read_ortho();
	read_lambda();
	//printf("Check 5\n");
	double *results_e =  malloc( sizeof(double)*2);
	results_e[0] = 0.0;
	results_e[1] = 0.0;
	double c1,c2,c3,c4,x;
	double k1,k2,k3, sum, sum1,sum2,sum3,sum4,sum5,sum6;
	double ksum, grid, factor, cube_size;
	int q1,q2,q3;
	double result2;

	ksum;
	grid = 1e0/((double)alpha_points);
	factor = kmax*grid;	// = \delta k
	cube_size = factor*factor*factor;	// = \delta k^3	

	if (rank == 0){
		for (i = 0; i < alpha_points + 1; i++){
			for (j = 0; j < alpha_points + 1; j++){
				for (k = 0; k < alpha_points + 1; k++){
					k1 = xvec[i];
					k2 = xvec[j];
					k3 = xvec[k];
					volumes[i][j][k] = calculate_weight(0, alpha_points, k1, k2, k3);
				}
			}
		}
	}

	gsl_matrix *cholesky = gsl_matrix_alloc(ortho_size,ortho_size);

	double y;
	double y_send = 0.0;
	double y_recv = 0.0;
	int eigen_size = ortho_size;
	double *eigen_R_model = (double *)create_vector(eigen_size);
	double *eigen = (double *)create_vector(eigen_size);
	double *pre_eigen = (double *)create_vector(eigen_size);

	int p1,p2,p3;
	double corr1, corr2, corr3, corr;
	//printf("Check 5\n");
	double* eigen_send = (double *)malloc( eigen_size * sizeof(double) );
	double* eigen_recv = (double *)malloc( eigen_size * sizeof(double) );

	//double* pre_eigen_send = (double *)malloc( eigen_size * sizeof(double) );
	//double* pre_eigen_recv = (double *)malloc( eigen_size * sizeof(double) );
	//printf("Check 6\n");
	shape_params params;
	params.a1 = 0e0;
	params.a2 = 0e0;
	params.a3 = 0e0;
	
	create_eigen();

	bool do_correlation = true;
	
	time1 = MPI_Wtime();
	for(p1=slim.p1_bgn;p1<slim.p1_end+1;p1+=slim.p1_stp){
		for(p2=slim.p2_bgn;p2<slim.p2_end+1;p2+=slim.p2_stp){
			for(p3=slim.p3_bgn;p3<slim.p3_end+1;p3+=slim.p3_stp){
				params.a1 = (double)p1;
				params.a2 = (double)p2;
				params.a3 = (double)p3;
				for (i=0;i<eigen_size;i++) {
					eigen[i] = 0.0;
				}
				// Compute   < Q(k1,k2,3), S(k1,k2,k3) >
				calculate_eigen(0,alpha_points,params,block,eigen);
				// end of MPI loop
				// printf("[%d] finished\n", rank);
				// duration = MPI_Wtime() - time1;
				// printf("[%d] done eigen: %e\n",rank,duration);
				time1 = MPI_Wtime();
				MPI_Barrier(MPI_COMM_WORLD);
				for (i=0;i<eigen_size;i++){
					eigen_send[i] = eigen[i];
				}
				MPI_Reduce(eigen_send,eigen_recv,eigen_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				duration = MPI_Wtime() - time1;
				// printf("[%d] done reduction: %e\n",rank,duration);
				time1 = MPI_Wtime();
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank==0){
					for(i=0;i<eigen_size;i++){
						x = eigen_recv[i];
					}
					gsl_vector *eigenin = gsl_vector_alloc(eigen_size);
					gsl_vector *eigenout = gsl_vector_alloc(eigen_size);
					gsl_vector_set_zero(eigenin);
					gsl_vector_set_zero(eigenout);
					for(i=0;i<eigen_size;i++){
						x = eigen_recv[i];
						gsl_vector_set(eigenin, i, x);
					}
					for(i=0;i<ortho_size;i++){
						for(j=0;j<ortho_size;j++){
							x = get_ortho(i,j);
				 			gsl_matrix_set(cholesky, i, j, x);
						}
					}
					gsl_linalg_cholesky_solve(cholesky, eigenin, eigenout);
					//printf("Alpha Q\n");
					for(i=0;i<eigen_size;i++){
						results_e[0] = (double)i;
						results_e[1] = gsl_vector_get(eigenout, i);
						update_eigen(results_e);
						//printf("[%d]\t%d\t%e\t%e\n", rank, i, results_e[1], get_eigen(i));
					}
					output_eigen(p1,p2,p3);
					duration = MPI_Wtime() - time1;
					//printf("done eigen:\t%e\t%e\t%e\tin\t%e\n",(double)p1,(double)p2,(double)p3,duration);
					time1 = MPI_Wtime();
					if (do_correlation){
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
						corr2 = 1e0/sqrt(corr2);
						corr3 = 1e0/sqrt(corr3);
						corr = corr1 * corr2 * corr3;
						printf("%d\t%d\t%d\t%e\n", p1, p2, p3, corr);
					}
				}
			}
		}
	}
	//---------------------------------------------------------------------------------------------------------------
	bool plot_bispectrum = false;
	bool print_shapes = false;
	bool brute_corr = false;



	if (rank == 0 ){
		double ksum;
		double grid = 1e0/((double)alpha_points);
		double factor = kmax*grid;	// = \delta k
		double cube_size = factor*factor*factor;	// = \delta k^3	
		//params.a1 = 1.0;
		//params.a2 = 0.0;
		read_eigen(p1,p2,p3);
		if (brute_corr){
			corr1 = 0e0;
			corr2 = 0e0;
			corr3 = 0e0;
			corr = 0e0;
			// maybe set beginning of i from 1 to exclude (0,k2,k3) points
			for(i=0;i<alpha_points+1;i+=1){
				k1 = xvec[i];
				// printf("%d\t%e\n",i, get_basis_prim(i,1));
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
			printf("----------------------------------------------------------------------\n");
			printf("Rec_Temp = %e\nRec_Rec = %e\nTemp_Temp = %e\n", corr1, corr2, corr3);
			corr2 = 1e0/sqrt(corr2);
			corr3 = 1e0/sqrt(corr3);
			corr = corr1 * corr2 * corr3;
			printf("Shape cosine = %e\n", corr);
		}
		if (print_shapes){
			double K;
			corr1 = 0e0;
			corr2 = 0e0;
			corr3 = 0e0;
			printf("K\tx\tk1\tk2\tk3\tReconstructed\tTheory\n");
			for(i=0;i<alpha_points+1;i+=1){
				k1 = xvec[i];
				// printf("%d\t%e\n",i, get_basis_prim(i,1));
				for(j=0;j<alpha_points+1;j+=1){
					k2 = xvec[j];
					for(k=0;k<alpha_points+1;k+=1){
						k3 = xvec[k];
						K = k1 + k2 + k3;
						// slices: (((K > 0.49*kmax) && (K < 0.51 * kmax)) || ((K > 0.29*kmax) && (K < 0.31 * kmax)) || ((K > 0.74*kmax) && (K < 0.76 * kmax)))
						// bulk: ((K > 1.5 * kmax) && (K < 2.5 * kmax))
						// edges: ((K < 0.3 * kmax) || (K > 2.8 * kmax))
						// interior: ((K > 0.3 * kmax) || (K < 2.8 * kmax))
						// bottom: (K < 0.3 * kmax)
						// top: (K > 2.8 * kmax)
						if((k1+k2>=k3) && (k1+k3>=k2) && (k2+k3>=k1)){
							c1 = 0e0;
							x = volumes[i][j][k] * cube_size;
							if (K > 1e-10) x = x / K;
							else x = 0e0;
							for(n=0;n<ortho_size;n++){
								c1 += get_eigen(n) * pijk(n,i,j,k);
							} 
							c2 = shape3(k1,k2,k3,params);
							corr1 += c1 * c2 * x;
							corr2 += c1 * c1 * x;
							corr3 += c2 * c2 * x;

							//printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n", K, x, k1, k2, k3, c1, c2);
						}
					}
				}
			}
			printf("Rec_Temp = %e\nRec_Rec = %e\nTemp_Temp = %e\n", corr1, corr2, corr3);
			corr2 = 1e0/sqrt(corr2);
			corr3 = 1e0/sqrt(corr3);
			corr = corr1 * corr2 * corr3;
			printf("Shape cosine = %e\n", corr);
		}
		if (plot_bispectrum){
			printf("Primordial bispectrum reconstruction\n--------------------------------------------------------------\n");
            double weight;
            double *** bisp_recon_array = (double ***)create_3Darray(alpha_points+1, alpha_points+1, alpha_points+1);
			double *** bisp_true_array = (double ***)create_3Darray(alpha_points+1,alpha_points+1,alpha_points+1);
			for(i=0;i<alpha_points+1;i++){
				k1 = xvec[i];
				// printf("%d\t%e\n",i, get_basis_prim(i,1));
				for(j=0;j<alpha_points+1;j++){
					k2 = xvec[j];
					for(k=0;k<alpha_points+1;k++){
						k3= xvec[k];
						if((k1>=k2) && (k1+k2>=k3) && (k1+k3>=k2) && (k2+k3>=k1)){
							c1 = 0e0;
							for(n=0;n<ortho_size;n++){
								find_perm_prim(n,&q1,&q2,&q3);
								sum1 = get_basis_prim(i,q1)*get_basis_prim(j,q2)*get_basis_prim(k,q3);
								sum2 = get_basis_prim(i,q2)*get_basis_prim(j,q3)*get_basis_prim(k,q1);
								sum3 = get_basis_prim(i,q3)*get_basis_prim(j,q1)*get_basis_prim(k,q2);
								sum4 = get_basis_prim(i,q3)*get_basis_prim(j,q2)*get_basis_prim(k,q1);
								sum5 = get_basis_prim(i,q2)*get_basis_prim(j,q1)*get_basis_prim(k,q3);
								sum6 = get_basis_prim(i,q1)*get_basis_prim(j,q3)*get_basis_prim(k,q2);
								c1 += get_eigen(n) * (sum1+sum2+sum3+sum4+sum5+sum6)/6e0;
							}
							c2 = shape3(k1,k2,k3,params);
							bisp_recon_array[i][j][k] = c1;
							bisp_true_array[i][j][k] = c2;
                        }
						else{
							bisp_recon_array[i][j][k] = NAN;
                            bisp_true_array[i][j][k] = NAN;
						}
                    }
                }
            }
			//bisp_recon_array[alpha_points+1][alpha_points+1][alpha_points+1] = bisp_recon_array[alpha_points+1][alpha_points+1][alpha_points+1];    
			//bisp_recon_array[alpha_points-10][alpha_points-10][alpha_points-9] = bisp_recon_array[alpha_points-10][alpha_points-10][alpha_points-8];

			hid_t h5file, h5file2, dataspace, dataspace2, dataset, dataset2;
			herr_t status;
			hsize_t dims[3];
			char bispectrum_h5file[100],bispectrum_h5file2[100], modelstr[50], termsstr[50], modelstr2[50], termsstr2[50]; 
			sprintf(modelstr,"%i.h5", model);
			sprintf(termsstr,"%i_",alpha_max);
			strcpy(bispectrum_h5file,"/nfs/st01/hpc-gr-epss/ps792/PlanckV3/C/recon_prim_bispectrum");
			strcat(bispectrum_h5file,termsstr);
			strcat(bispectrum_h5file,modelstr);

			sprintf(modelstr2,"%i.h5", model);
			sprintf(termsstr2,"%i_",alpha_max);
			strcpy(bispectrum_h5file2,"/nfs/st01/hpc-gr-epss/ps792/PlanckV3/C/th_prim_bispectrum");
			strcat(bispectrum_h5file2,termsstr);
			strcat(bispectrum_h5file2,modelstr);
			printf(bispectrum_h5file2);
			printf("\n");

			h5file = H5Fcreate(bispectrum_h5file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			if(h5file<0){
				printf("Failed to create file! \n");
				exit(1);
			}
			dims[0] = alpha_points+1;
			dims[1] = alpha_points+1;
			dims[2] = alpha_points+1;
			
			dataspace = H5Screate_simple(3, dims, NULL);
			if(dataspace<0){
				printf("Failed to create dataset! \n");
                exit(1);
            }

            dataset = H5Dcreate(h5file, "bispectrum", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (dataset < 0){
                printf("Failed to create dataset!\n");
                exit(1);
            }
            H5Sclose(dataspace);
            // Write the array data to the dataset
            status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bisp_recon_array[0][0][0]);
            if (status < 0){
                printf("Failed to write data to dataset!\n");
                exit(1);
            }
            // Close the dataset, dataspace, and file
			status = H5Dclose(dataset);
			if (status < 0) {
					printf("Failed to close dataset!\n");
					exit(1);
			}
			status = H5Fclose(h5file);
			if (status < 0) {
					printf("Failed to close file!\n");
					exit(1);
			}
			
			printf("Success to write recon_bispectrum to dataset.\n");

			h5file2 = H5Fcreate(bispectrum_h5file2, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			dims[0] = alpha_points+1;
			dims[1] = alpha_points+1;
			dims[2] = alpha_points+1;
			dataspace2 = H5Screate_simple(3, dims, NULL);
			dataset2 = H5Dcreate(h5file2, "bispectrum", H5T_NATIVE_DOUBLE, dataspace2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			H5Sclose(dataspace2);
			status = H5Dwrite(dataset2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bisp_true_array[0][0][0]);
			H5Dclose(dataset2);
			H5Fclose(h5file2);

			printf("Success to write th_bispectrum to dataset.\n");
		}
	}

	MPI_Finalize();
}