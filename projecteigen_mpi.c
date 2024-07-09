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

int main( int argc, char *argv[] ){
	
// **1**

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
	
	int i,j,k,l,m,n,r,s,t;
	double result=0;

	int ortho_size = get_terms_prim();
	if(rank==0){
		printf("model: %d\n",model);
		printf("pmax: %d\n", get_pmax_prim());
		printf("lmax = %d\n", eflag_T_lmax);
		/*printf("TTT\n");
		for(n=0;n<ortho_size;n++){
			find_perm_late_TTT(n,&i,&j,&k);
			printf("%d\t%d\t%d\t%d\n",n,i,j,k);
		}
		if (do_polarisation){
			printf("TTE\n");
			for(n=0;n<ortho_size;n++){
				find_perm_late_TTE(n,&i,&j,&k);
				printf("%d\t%d\t%d\t%d\n",n,i,j,k);
			}
			printf("TEE\n");
			for(n=0;n<ortho_size;n++){
				find_perm_late_TEE(n,&i,&j,&k);
				printf("%d\t%d\t%d\t%d\n",n,i,j,k);
			}
			printf("EEE\n");
			for(n=0;n<ortho_size;n++){
				find_perm_late_EEE(n,&i,&j,&k);
				printf("%d\t%d\t%d\t%d\n",n,i,j,k);
			}
		}
		*/
	}
		
	MPI_Barrier(MPI_COMM_WORLD);
	
	read_ortho();
	read_lambda();
	read_orthol();
	read_lambdal();
	
	double *results_e =  malloc( sizeof(double)*2);
	results_e[0] = 0.0;
	results_e[1] = 0.0;
	double c1,c2,c3,c4,x,x1,x2,x3,x4;
	double result2;

	int eigen_size = ortho_size;
	double *eigen = (double *)create_vector(eigen_size);
	
	for (i=0;i<eigen_size;i++) {
		eigen[i] = 0.0;
	}
	
	double *results_m;
	double normT;
	double normTE;

	read_qtilde();
	int ell = get_qtilde_E_lsize();
	int exx = get_qtilde_xsize();
	int p1 = 0;
	int p2 = 1;
	int p3 = 2;
	int p4 = 3;
	/*for(int i = 0; i < ell; i++){
		for (int j = 0; j < exx; j++){
			printf("%e\t%e\t%e\t%e\n", get_qtilde_E(p1,i,j),get_qtilde_E(p2,i,j),get_qtilde_E(p3,i,j),get_qtilde_E(p4,i,j));
		}
	}*/

	
	if(rank==0){
		int p1,p2,p3;
		if(do_polarisation==1){
			results_m =  (double *)malloc( sizeof(double)*5);
		}else{
			results_m =  (double *)malloc( sizeof(double)*2);
		}
		create_modes();
		read_gamma();
		for(p1=slim.p1_bgn;p1<slim.p1_end+1;p1+=slim.p1_stp){
			for(p2=slim.p2_bgn;p2<slim.p2_end+1;p2+=slim.p2_stp){
				for(p3=slim.p3_bgn;p3<slim.p3_end+1;p3+=slim.p3_stp){
					normT = 0.0;
					normTE = 0.0;
					read_eigen(p1,p2,p3);
					//printf("Modes Q\n");
					for(i=0;i<ortho_size;i++){
						results_m[0] = (double)i;
						results_m[1] = 0.0;
						if(do_polarisation==1){
							results_m[2] = 0.0;
							results_m[3] = 0.0;
							results_m[4] = 0.0;
						}
						for(j=0;j<ortho_size;j++){
							x = get_eigen(j);
							c1 = get_gamma_TTT(i,j);
							//printf("gamma_TTT = %e\t", c1);
							//printf("eigen = %e\n",x);
							if(do_polarisation==1){
								c2 = get_gamma_TTE(i,j);
								c3 = get_gamma_TEE(i,j);
								c4 = get_gamma_EEE(i,j);
								//printf("gamma_TTE\tgamma_TEE\tgamma_EEE = %e\t%e\t%e\n", c2,c3,c4);
							}
							results_m[1] += c1*x;
							if(do_polarisation==1){
								results_m[2] += c2*x;
								results_m[3] += c3*x;
								results_m[4] += c4*x;
							}
						}
						update_modes(results_m);
						if(do_polarisation==1){
							//printf("%d\t%e\t%e\t%e\t%e\t%e\n", i, get_modes_TTT(i), get_modes_TTE(i), get_modes_TEE(i), get_modes_EEE(i), get_eigen(i));
						}else{
							//printf("%d\t%e\t%e\n", i, get_modes_TTT(i), get_eigen(i));
						}
					}
					output_modes(p1,p2,p3);
					x1 = 0e0;
					x2 = 0e0;
					x3 = 0e0;
					x4 = 0e0;
					//printf("Modes R\n");
					for(i=0;i<ortho_size;i++){
						//printf("%d iteration\n", i);
						c1 = 0e0;
						c2 = 0e0;
						c3 = 0e0;
						c4 = 0e0;
						for(j=0;j<ortho_size;j++){
							c1 += get_modes_TTT(j)*get_orthol_TTT(j,i);
							//printf("c1 = %f\n", c1);
							if (do_polarisation==1){
								c2 += get_modes_TTE(j)*get_orthol_TTE(j,i);
								//printf("c2 = %f\n", c2);
								c3 += get_modes_TEE(j)*get_orthol_TEE(j,i);
								//printf("c3 = %f\n", c3);
								c4 += get_modes_EEE(j)*get_orthol_EEE(j,i);
								//printf("c4 = %f\n", c4);
							}
						}
						x1 += c1*c1;
						x2 += c2*c2;
						x3 += c3*c3;
						x4 += c4*c4;
						//printf("CHECK 1\n");
						//printf("%d\t%e\t%e\t%e\t%e\n", i, c1, c2, c3, c4);
					}
					//printf("WANT X4 AND X1.....\n");
					normT = sqrt(6.0 / x1);
					normTE = sqrt(6.0 / (x1 + 3e0*x2 + 3e0*x3 + x4));
					//x4 = 7.76e-1*(x4+3e0*x3+3e0*x2+x1)/6e0;
					//x1 = 7.76e-1*x1/6e0;
					//printf("X4 = %e\tX3 = %e\n", x4, x1);
					//printf("%d\t%d\t%d\t%e\t%e\n",p1,p2,p3,1e0/sqrt(x1),1e0/sqrt(x4));
					printf("%d\t%d\t%d\t%e\t%e\n",p1,p2,p3,normT,normTE);
				}
			}
		}
		printf("Want to destroy modes...\n");
		destroy_modes();
		
		printf("Want to read modes.\n");
		read_modes(0,0,0);
		printf("Modes successfully read!\n");
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
 				QTT[i][r] = get_basis_late_T(r-Tmin,i);
 				//QTT[i][r] = get_basis_late_T(r-Tmin,i) * weight[r];
			        // if(i>pmaxT-3)printf("%d\t%e\t%e\t%e\n",r,QTT[i][r],get_basis_late_T(r-Tmin,i),weight[r]);
			}
		}
		
		int l1,l2,l3,L;
		double x,x1,x2,x3,x4,x5,x6;
		double y,y1,y2,y3,y4,y5,y6;
		
		printf("l1\tl2\tl3\tS(l1,l2,l3)\tb(l1,l2,l3)\n");
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
					y=0e0;
					L = l1 + l2 + l3;
					if(l1+l2>=l3 && l2+l3>=l1 && l1+l3>=l2 && L%2==0 && (L == 500 || L == 1000 || L == 1500 || L==2000 || L==2500)){
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
							y += x * weight[l1] * weight[l2] * weight[l3];
						}
						printf("%d\t%d\t%d\t%d\t%e\t%e\n", L, l1, l2, l3, x, y);
					}
				}
			}
		}
		
		
		// x = 0e0;
		// for(n=0;n<ortho_size;n++){
		// 	find_perm_late_TTT(n,&i,&j,&k);
		// 	x1 = QTT[i][2]*QTT[j][2]*QTT[k][2];
		// 	x2 = QTT[j][2]*QTT[k][2]*QTT[i][2];
		// 	x3 = QTT[k][2]*QTT[i][2]*QTT[j][2];
		// 	x4 = QTT[k][2]*QTT[j][2]*QTT[i][2];
		// 	x5 = QTT[j][2]*QTT[i][2]*QTT[k][2];
		// 	x6 = QTT[i][2]*QTT[k][2]*QTT[j][2];
		// 	x += get_modes_TTT(n)*(x1+x2+x3+x4+x5+x6)/6e0;
		// 	printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",n,x1,x2,x3,x4,x5,x6,get_modes_TTT(n),x);
		// }
		/*
		bool plot_bispectrum = false;
	
		// Reconstruct late time bispectrum
		if (plot_bispectrum){
			int l_step = 10;
			int l_start = 10;
			int l_range = ((eflag_T_lmax - l_start)/l_step);
			printf("Late time bispectrum reconstruction\n------------------------------------------------\n");
			//double *** bisp_recon_array = (double ***)create_3Darray(eflag_T_lmax+1, eflag_T_lmax+1, eflag_T_lmax+1);
			double *** bisp_recon_array = (double ***)create_3Darray(l_range+1, l_range+1, l_range+1);
			for(r=l_start;r<eflag_T_lmax+1;r+=l_step){
				for(s=l_start;s<eflag_T_lmax+1;s+=l_step){
					for(t=l_start;t<eflag_T_lmax+1;t+=l_step){
						l1 = r;
						l2 = s;
						l3 = t;
						//if(l1==0)bisp_recon_array[l1][l2][l3] = NAN;
						//if(l2==0)bisp_recon_array[l1][l2][l3] = NAN;
						//if(l3==0)bisp_recon_array[l1][l2][l3] = NAN;

						if(l1==0)l1=eflag_T_lmin;
						if(l2==0)l2=eflag_T_lmin;
						if(l3==0)l3=eflag_T_lmin;
						x=0e0;
						L = l1 + l2 + l3;
						if(l1+l2>=l3 && l2+l3>=l1 && l1+l3>=l2 && L%2==0){
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
							bisp_recon_array[l1/l_step-1][l2/l_step-1][l3/l_step-1] = x;
							//printf("%d\t%d\t%d\t%e\n",l1,l2,l3,x);
						}
						else{
							bisp_recon_array[l1/l_step-1][l2/l_step-1][l3/l_step-1] = NAN;
						}
					}
				}
			}
			printf("Triple for loop passed.\n");
			hid_t h5file, dataspace, dataset;
			herr_t status;
			hsize_t dims[3];
			char bispectrum_h5file[50], modelstr[50], termsstr[50];
			sprintf(modelstr,"%i.h5", model);
			sprintf(termsstr,"%e_",Tmax);
			strcpy(bispectrum_h5file,"/nfs/st01/hpc-gr-epss/ps792/PlanckV3/C/recon_late_bispectrum_9_");
			strcat(bispectrum_h5file,termsstr);
			strcat(bispectrum_h5file,modelstr);
			printf("\n");

			h5file = H5Fcreate(bispectrum_h5file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			if(h5file<0){
				printf("Failed to create file! \n");
				exit(1);
			}
			dims[0] = l_range+1;
			dims[1] = l_range+1;
			dims[2] = l_range+1;
			
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
			//printf("File name = %c\n", bispectrum_h5file);
			printf("Success to write recon_bispectrum to dataset.\n");
		}
		*/
		// read_orthol();
		//
		// printf("Modes R\n");
		// for(i=0;i<ortho_size;i++){
		// 	c1=0.0;
		// 	if(do_polarisation==1){
		// 		c2=0e0;
		// 		c3=0e0;
		// 		c4=0e0;
		// 	}
		// 	for(j=0;j<ortho_size;j++){
		// 		c1 += get_orthol_TTT(j,i)*get_modes_TTT(j);
		// 		if(do_polarisation==1){
		// 			c2 += get_orthol_TTE(j,i)*get_modes_TTE(j);
		// 			c3 += get_orthol_TEE(j,i)*get_modes_TEE(j);
		// 			c4 += get_orthol_EEE(j,i)*get_modes_EEE(j);
		// 		}
		// 	}
		// 	if(do_polarisation==1){
		// 		printf("%d\t%e\t%e\t%e\t%e\n", i, c1, c2, c3, c4);
		// 	}else{
		// 		printf("%d\t%e\n", i, c1);
		// 	}
		// }
			
	printf("Code finished.");
	}
	
	MPI_Finalize();
	return 0;
}