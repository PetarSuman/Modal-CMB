#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "global.h"

struct my_int_params{int p1; int p2; int p3; int p4;};
double integrand_tri_QS(double *k, size_t dim, void *params);

double calculate_eigen_tri(int r){
 
 	double integral=0;
	double ****points = create_4Darray(2,2,2,2);
 	
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	double perm;
	int i,j,k,l,m,n,o,p;
	
	double k1,k2,k3,k4,ksum;
	double grid = 1e0/((double)alpha_points);
	double max = get_kmax_cut();
	double factor = max*grid;
 	double cube_size = factor*factor*factor;
	
	for(i=0;i<alpha_points;i++){
		for(j=i;j<alpha_points;j++){
			for(k=j;k<alpha_points;k++){
				for(l=k;l<alpha_points;l++){
					
					if(l<=i+j+k+2){
						if(i==j){
							if(j==k){
								if(k==l){
									perm=1.0;
								}else{
									perm=4.0;
								}
							}else{
								if(k==l){
									perm=6.0;
								}else{
									perm=12.0;
								}
							}
						}else{
							if(j==k){
								if(k==l){
									perm=4.0;
								}else{
									perm=12.0;
								}
							}else{
								if(k==l){
									perm=12.0;
								}else{
									perm=24.0;
								}
							}
						}
						for(m=0;m<2;m++){
 							k1 = factor*(i+m);
							for(n=0;n<2;n++){
 								k2 = factor*(j+n);
								for(o=0;o<2;o++){
 									k3 = factor*(k+o);
									for(p=0;p<2;p++){
 										k4 = factor*(l+p);
 										
										points[m][n][o][p] = triQ(r,i+m,j+n,k+o,l+p)*shape_tri3(k1,k2,k3,k4);
 										
//  										ksum = k1+k2+k3+k4;
//  										if (ksum>1e-10) ksum = 1e0/ksum;
// 										points[m][n][o][p] = triQ(r,i+m,j+n,k+o,l+p)*shape_tri3(k1,k2,k3,k4)*ksum;
									}
								}
							}
						}
						integral += perm * calculate_volume_tri(i,j,k,l,points) * cube_size;
					}
				}
				
			}
		}
	}
	return integral;
}

 double calculate_eigen_tri_mc(int r){

	double res1, err1;
	double res2, err2;
	double res3, err3;
	double kmin = 2e0/get_tau0();
	double kmax = get_kmax_cut();
	double xl[4] = {kmin,kmin,kmin,kmin};
	double xu[4] = {kmax,kmax,kmax,kmax};
	
	const gsl_rng_type *T;
	gsl_rng *rng;
	
	int p1,p2,p3,p4;
	find_perm_tri_prim(r,&p1,&p2,&p3,&p4);

	struct my_int_params params = {p1,p2,p3,p4};
	gsl_monte_function QS = { &integrand_tri_QS, 4, &params};
	
	size_t calls = 20000000;
	
	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);
	
	gsl_monte_plain_state *state = gsl_monte_plain_alloc(4);
	gsl_monte_plain_integrate (&QS, xl, xu, 4, calls, rng, state, &res1, &err1);
	gsl_monte_plain_free(state);
	gsl_rng_free (rng);
	
	return res1;
}

double integrand_tri_QS(double *k, size_t dim, void *params){
	struct my_int_params *p = (struct my_int_params *)params;

	double k1 = k[0];
	double k2 = k[1];
	double k3 = k[2];
	double k4 = k[3];
	int p1,p2,p3,p4;
	p1 = p->p1;
	p2 = p->p2;
	p3 = p->p3;
	p4 = p->p4;
	
	double kinv = 1e0/(k1+k2+k3+k4);
	
	double result;
	double b1,b2,b3,b4,b5,b6;
	double part1,part2;
	double sixinv = 1e0/6e0;
	
	if(k1>k2+k3+k4||k2>k3+k4+k1||k3>k4+k1+k2||k4>k1+k2+k3){
		result = 0e0;
	}else{
 		b1 = sinlog_pt(k1,p1)*sinlog_pt(k2,p2)*sinlog_pt(k3,p3)*sinlog_pt(k4,p4);
 		b2 = sinlog_pt(k1,p1)*sinlog_pt(k2,p3)*sinlog_pt(k3,p4)*sinlog_pt(k4,p2);
 		b3 = sinlog_pt(k1,p1)*sinlog_pt(k2,p4)*sinlog_pt(k3,p2)*sinlog_pt(k4,p3);
 		b4 = sinlog_pt(k1,p1)*sinlog_pt(k2,p4)*sinlog_pt(k3,p3)*sinlog_pt(k4,p2);
 		b5 = sinlog_pt(k1,p1)*sinlog_pt(k2,p3)*sinlog_pt(k3,p2)*sinlog_pt(k4,p4);
 		b6 = sinlog_pt(k1,p1)*sinlog_pt(k2,p2)*sinlog_pt(k3,p4)*sinlog_pt(k4,p3);
 		
 		b1 += sinlog_pt(k1,p2)*sinlog_pt(k2,p1)*sinlog_pt(k3,p3)*sinlog_pt(k4,p4);
 		b2 += sinlog_pt(k1,p2)*sinlog_pt(k2,p3)*sinlog_pt(k3,p4)*sinlog_pt(k4,p1);
 		b3 += sinlog_pt(k1,p2)*sinlog_pt(k2,p4)*sinlog_pt(k3,p1)*sinlog_pt(k4,p3);
 		b4 += sinlog_pt(k1,p2)*sinlog_pt(k2,p4)*sinlog_pt(k3,p3)*sinlog_pt(k4,p1);
 		b5 += sinlog_pt(k1,p2)*sinlog_pt(k2,p3)*sinlog_pt(k3,p1)*sinlog_pt(k4,p4);
 		b6 += sinlog_pt(k1,p2)*sinlog_pt(k2,p1)*sinlog_pt(k3,p4)*sinlog_pt(k4,p3);
 		
 		b1 += sinlog_pt(k1,p3)*sinlog_pt(k2,p2)*sinlog_pt(k3,p1)*sinlog_pt(k4,p4);
 		b2 += sinlog_pt(k1,p3)*sinlog_pt(k2,p1)*sinlog_pt(k3,p4)*sinlog_pt(k4,p2);
 		b3 += sinlog_pt(k1,p3)*sinlog_pt(k2,p4)*sinlog_pt(k3,p2)*sinlog_pt(k4,p1);
 		b4 += sinlog_pt(k1,p3)*sinlog_pt(k2,p4)*sinlog_pt(k3,p1)*sinlog_pt(k4,p2);
 		b5 += sinlog_pt(k1,p3)*sinlog_pt(k2,p1)*sinlog_pt(k3,p2)*sinlog_pt(k4,p4);
 		b6 += sinlog_pt(k1,p3)*sinlog_pt(k2,p2)*sinlog_pt(k3,p4)*sinlog_pt(k4,p1);
 		
 		b1 += sinlog_pt(k1,p4)*sinlog_pt(k2,p2)*sinlog_pt(k3,p3)*sinlog_pt(k4,p1);
 		b2 += sinlog_pt(k1,p4)*sinlog_pt(k2,p3)*sinlog_pt(k3,p1)*sinlog_pt(k4,p2);
 		b3 += sinlog_pt(k1,p4)*sinlog_pt(k2,p1)*sinlog_pt(k3,p2)*sinlog_pt(k4,p3);
 		b4 += sinlog_pt(k1,p4)*sinlog_pt(k2,p1)*sinlog_pt(k3,p3)*sinlog_pt(k4,p2);
 		b5 += sinlog_pt(k1,p4)*sinlog_pt(k2,p3)*sinlog_pt(k3,p2)*sinlog_pt(k4,p1);
 		b6 += sinlog_pt(k1,p4)*sinlog_pt(k2,p2)*sinlog_pt(k3,p1)*sinlog_pt(k4,p3);
 		
 		part1 = (b1+b2+b3+b4+b5+b6)*sixinv;

 		part2 = shape_tri3(k1,k2,k3,k4);
 		
		result = kinv*part1*part2;
	}
	return result;
}

double check_tri(double kmax){
 	double integral1=0;
 	double integral2=0;
 	double integral3=0;
	double ****points1 = create_4Darray(2,2,2,2);
	double ****points2 = create_4Darray(2,2,2,2);
	double ****points3 = create_4Darray(2,2,2,2);
	
	int alpha_cut = (int)( (double)(eflag_lmax*alpha_points) / (kmax*get_tau0()) );
	if(alpha_cut>alpha_points) alpha_cut = alpha_points;
	printf("alpha_cut:\t%d\n",alpha_cut);
 	double cube_size = pow((double)alpha_cut,-4);
 	double result;
	
	double k1,k2,k3,k4,ksum;

	int i,j,k,l,m,n,o,p;
	double perm;

	for(i=0;i<alpha_cut;i++){
		for(j=i;j<alpha_cut;j++){
			for(k=j;k<alpha_cut;k++){
				for(l=k;l<alpha_cut;l++){
					
					if(l<=i+j+k+2){
						if(i==j){
							if(j==k){
								if(k==l){
									perm=1.0;
								}else{
									perm=4.0;
								}
							}else{
								if(k==l){
									perm=6.0;
								}else{
									perm=12.0;
								}
							}
						}else{
							if(j==k){
								if(k==l){
									perm=4.0;
								}else{
									perm=12.0;
								}
							}else{
								if(k==l){
									perm=12.0;
								}else{
									perm=24.0;
								}
							}
						}
						for(m=0;m<2;m++){
							for(n=0;n<2;n++){
								for(o=0;o<2;o++){
									for(p=0;p<2;p++){
 										k1 = kmax*(double)(i+m)/alpha_points;
 										k2 = kmax*(double)(j+n)/alpha_points;
 										k3 = kmax*(double)(k+o)/alpha_points;
 										k4 = kmax*(double)(l+p)/alpha_points;
										points1[m][n][o][p] = triS(i+m,j+n,k+o,l+p)*shape_tri3(k1,k2,k3,k4);
										points2[m][n][o][p] = triS(i+m,j+n,k+o,l+p)*triS(i+m,j+n,k+o,l+p);
										points3[m][n][o][p] = shape_tri3(k1,k2,k3,k4)*shape_tri3(k1,k2,k3,k4);
									}
								}
							}
						}
						integral1 += perm * calculate_volume_tri(i,j,k,l,points1) * cube_size;
						integral2 += perm * calculate_volume_tri(i,j,k,l,points2) * cube_size;
						integral3 += perm * calculate_volume_tri(i,j,k,l,points3) * cube_size;
// 						printf("%e\t%e\t%e\n",integral1,integral2,integral3);
					}
				}
			}
		}
	}
	printf("%e\t%e\t%e\n",integral1,integral2,integral3);
	result = integral1 / sqrt(integral2*integral3);
// 	printf("%e\n",result);
	return result;
}

double check_eigen_tri(int num, double max){

 	double integral1=0;
 	double integral2=0;
 	double integral3=0;
	double ****points1 = create_4Darray(2,2,2,2);
	double ****points2 = create_4Darray(2,2,2,2);
	double ****points3 = create_4Darray(2,2,2,2);
 	double cube_size = pow((double)alpha_points,-4);
 	double recon,factor;
	int terms=get_terms_prim();
	double* eigen = (double *)create_vector(terms);
	double* eigenR = (double *)create_vector(terms);
	double** ortho = (double **)create_array(terms,terms);
	double** ortho_inv = (double **)create_array(terms,terms);
	double** dummy = (double **)create_array(terms,terms);
	double perm;
	
	double k1,k2,k3,k4;

	int i,j,k,l,m,n,o,p,q;
	
	for(n=0;n<terms;n++){
		eigen[n] = get_eigen_tri(n);
	}
	
	for(i=0;i<terms;i++){
		for(j=0;j<terms;j++){
			ortho[i][j] = get_ortho_tri(i,j);
			dummy[i][j] = get_ortho_tri(i,j);
			ortho_inv[i][j] = 0.0;
		}
	}
	
// 	Calculate alpha inverse	
	
	for(i=0;i<terms;i++){
		ortho_inv[i][i] = 1.0 / dummy[i][i];
		for(j=0;j<terms;j++){
			if(j>i){
				dummy[i][j] = 0.0;
				ortho[i][j] = 0.0;
			} else {
				dummy[i][j] = dummy[i][j] / dummy[i][i];
			}
		}
	}
	
	for(k=1;k<terms;k++){
		for(m=terms-1;m>=k;m--){
// 			m = terms-i;
			factor = -(dummy[m][m-k]/dummy[m-k][m-k]);
			for(j=0;j<=m-k;j++){
				ortho_inv[m][j] = ortho_inv[m][j] + factor*ortho_inv[m-k][j];
				dummy[m][j] = dummy[m][j] + factor*dummy[m-k][j];
			}
		}
	}
	
	for(i=0;i<terms;i++){
		for( j=0;j<terms;j++){
			eigenR[i] = eigenR[i] + ortho_inv[j][i]*eigen[j];
		}
	}

	for(i=0;i<alpha_points;i++){
		for(j=i;j<alpha_points;j++){
			for(k=j;k<alpha_points;k++){
				for(l=k;l<alpha_points;l++){
					
					if(l<=i+j+k+2){
						if(i==j){
							if(j==k){
								if(k==l){
									perm=1.0;
								}else{
									perm=4.0;
								}
							}else{
								if(k==l){
									perm=6.0;
								}else{
									perm=12.0;
								}
							}
						}else{
							if(j==k){
								if(k==l){
									perm=4.0;
								}else{
									perm=12.0;
								}
							}else{
								if(k==l){
									perm=12.0;
								}else{
									perm=24.0;
								}
							}
						}
						
						for(m=0;m<2;m++){
							for(n=0;n<2;n++){
								for(o=0;o<2;o++){
									for(p=0;p<2;p++){
										
										recon = 0.0;
										for(q=0;q<num+1;q++){
											recon += eigenR[q] * triR(q,i+m,j+n,k+o,l+p);
										}
										
 										k1 = max*(double)(i+m)/alpha_points;
 										k2 = max*(double)(j+n)/alpha_points;
 										k3 = max*(double)(k+o)/alpha_points;
 										k4 = max*(double)(l+p)/alpha_points;
 										
										points1[m][n][o][p] = recon*shape_tri3(k1,k2,k3,k4);
										points2[m][n][o][p] = recon*recon;
										points3[m][n][o][p] = shape_tri3(k1,k2,k3,k4)*shape_tri3(k1,k2,k3,k4);
										
									}
								}
							}
						}
						
						integral1 += perm * calculate_volume_tri(i,j,k,l,points1) * cube_size;
						integral2 += perm * calculate_volume_tri(i,j,k,l,points2) * cube_size;
						integral3 += perm * calculate_volume_tri(i,j,k,l,points3) * cube_size;
					}
				}
			}
		}
	}
// 	printf("%e\t%e\t%e\n",integral1,integral2,integral3);
	return integral1 / sqrt(integral2*integral3);
}

double check_norm_tri(double kmax){
 	double integral1=0;
 	double integral2=0;
	double ****points1 = create_4Darray(2,2,2,2);
	double ****points2 = create_4Darray(2,2,2,2);
	
	int alpha_cut = (int)( (double)(eflag_lmax*alpha_points) / (kmax*get_tau0()) );
	if(alpha_cut>alpha_points) alpha_cut = alpha_points;
	printf("alpha_cut:\t%d\n",alpha_cut);
 	double cube_size = pow((double)alpha_cut,-4);
 	double result;
	
	double k1,k2,k3,k4,ksum;

	int i,j,k,l,m,n,o,p;
	double perm;

	for(i=0;i<alpha_cut;i++){
		for(j=i;j<alpha_cut;j++){
			for(k=j;k<alpha_cut;k++){
				for(l=k;l<alpha_cut;l++){
					
					if(l<=i+j+k+2){
						if(i==j){
							if(j==k){
								if(k==l){
									perm=1.0;
								}else{
									perm=4.0;
								}
							}else{
								if(k==l){
									perm=6.0;
								}else{
									perm=12.0;
								}
							}
						}else{
							if(j==k){
								if(k==l){
									perm=4.0;
								}else{
									perm=12.0;
								}
							}else{
								if(k==l){
									perm=12.0;
								}else{
									perm=24.0;
								}
							}
						}
						for(m=0;m<2;m++){
							for(n=0;n<2;n++){
								for(o=0;o<2;o++){
									for(p=0;p<2;p++){
 										k1 = kmax*(double)(i+m)/alpha_points;
 										k2 = kmax*(double)(j+n)/alpha_points;
 										k3 = kmax*(double)(k+o)/alpha_points;
 										k4 = kmax*(double)(l+p)/alpha_points;
										points1[m][n][o][p] = shape_tri3(k1,k2,k3,k4);
										points2[m][n][o][p] = triS(i+m,j+n,k+o,l+p);
									}
								}
							}
						}
						integral1 += perm * calculate_volume_tri(i,j,k,l,points1) * cube_size;
						integral2 += perm * calculate_volume_tri(i,j,k,l,points2) * cube_size;
					}
				}
			}
		}
	}
	printf("%e\t%e\n",integral1,integral2);
	result = integral2 / integral1;
// 	printf("%e\n",result);
	return result;
}