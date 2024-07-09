#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "global.h"

struct my_int_params{int p1; int p2; int p3;};
double integrand_QS(double *k, size_t dim, void *params);

void integrate_template2(int min, int max, shape_params params, tetrapyd_limits* block, double* full_res){
	// This function computes the self-inner product of the theoretical template:
	//   it computes <S_th, S_th> in primordial space

	int n,r,t,i,j,k;
	double x,y;
	//y = 0.0;
	//full_res = 0.0;
	
	double k1,k2,k3,ksum;
	double grid = 1e0/((double)max);
	double kmax = get_kmax();
	double factor = kmax*grid;	// = \delta k
 	double cube_size = factor*factor*factor;	// = \delta k^3	

	#pragma omp parallel default(none) private(n,r,t,i,j,k,x,y,k1,k2,k3,ksum) shared(full_res, params, block, factor, min, max, cube_size)
	{
		t = omp_get_thread_num();
		i = block[t].i_bgn;
		j = block[t].j_bgn;
		k = block[t].k_bgn;

		y = 0.0;
		
		/*for(r=0;r<terms;r++){
			alpha_private[r] = 0e0;
		}*/
		
		for(n=0;n<block[t].loops;n++){
			// dV_k = calculate_weight computes what fraction of a cube lies inside the tetrapyd
			x = calculate_weight(min, max, i, j, k)*cube_size;
			//printf("CUBE SIZE = %e\t", cube_size);
			//printf("CALCULATE WEIGHT = %e\n", x);
			k1 = factor*(i);
			k2 = factor*(j);
			k3 = factor*(k);

			ksum = k1+k2+k3;
			if (ksum>1e-10) ksum = 1e0/ksum;

			//printf("ksum = %e\n", ksum);

			x *= pow(shape3(k1,k2,k3,params),2) * ksum;

			//printf("x = %e\n", x);

			y += x;
			//printf("Integral = %e\n", y);
			//printf("XSHAPE = %e\n",x);

			/*for(r=0;r<terms;r++){
				y = pijk(r,i,j,k);
				alpha_private[r] += x*y;
			}*/
			//printf("Y = PIJK = %e\t", y);
			get_ijk_next_prim(max,&i,&j,&k);
		}
		
		#pragma omp critical
		{
			/*for(r=0;r<terms;r++){
				alpha[r] += alpha_private[r];
				//printf("ALPHA = %e\t", alpha[r]);
				//printf("ALPHA PRIVATE = %e\n", alpha_private[r]);
			}*/
			*full_res += y;
			//printf("Full res = %e\n", *full_res);
		}
	}
	return;
}


void calculate_eigen(int min, int max, shape_params params, tetrapyd_limits* block, double* alpha){
	int n,r,t,i,j,k;
	double x,y;
	
	double k1,k2,k3,ksum;
	double grid = 1e0/((double)max);
	double kmax = get_kmax();
	double factor = kmax*grid;	// = \delta k
 	double cube_size = factor*factor*factor;	// = \delta k^3
	int terms = get_terms_prim();
	double alpha_private[terms];

	#pragma omp parallel default(none) private(n,r,t,i,j,k,x,y,k1,k2,k3,ksum,alpha_private) shared(alpha, params, block, factor, min, max, cube_size, terms)
	{
		t = omp_get_thread_num();
		i = block[t].i_bgn;
		j = block[t].j_bgn;
		k = block[t].k_bgn;
		
		for(r=0;r<terms;r++){
			alpha_private[r] = 0e0;
		}
		
		for(n=0;n<block[t].loops;n++){
			// dV_k = calculate_weight computes what fraction of a cube lies inside the tetrapyd
			x = calculate_weight(min, max, i, j, k)*cube_size;
			//printf("CUBE SIZE = %e\t", cube_size);
			//printf("CALCULATE WEIGHT = %e\n", x);
			k1 = factor*(i);
			k2 = factor*(j);
			k3 = factor*(k);

			ksum = k1+k2+k3;
			if (ksum>1e-10) ksum = 1e0/ksum;

			//printf("ksum = %e\n", ksum);

			x *= shape3(k1,k2,k3,params)*ksum;

			//printf("XSHAPE = %e\n",x);

			for(r=0;r<terms;r++){
				y = pijk(r,i,j,k);
				alpha_private[r] += x*y;
			}
			//printf("Y = PIJK = %e\t", y);
			get_ijk_next_prim(max,&i,&j,&k);
		}
		
		#pragma omp critical
		{
			for(r=0;r<terms;r++){
				alpha[r] += alpha_private[r];
				//printf("ALPHA = %e\t", alpha[r]);
				//printf("ALPHA PRIVATE = %e\n", alpha_private[r]);
			}
		}
	}
	return;
}

// double calculate_eigen_old(int r, shape_params params){
//
//  	double integral=0;
// 	double ***points = create_3Darray(2,2,2);
//
// 	int myrank;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//
// 	double perm;
// 	int i,j,k,l,m,n;
//
// 	double k1,k2,k3,ksum;
// 	double grid = 1e0/((double)alpha_points);
// 	double max = get_kmax();
// 	double factor = max*grid;
//  	double cube_size = factor*factor*factor;
//
// 	for(i=0;i<alpha_points;i++){
// 		for(j=i;j<alpha_points;j++){
// 			for(k=j;k<alpha_points;k++){
//
// 				if(k<=i+j+2){
// 					if(i==j){
// 						if(j==k){
// 							perm=1.0;
// 						}else{
// 							perm=3.0;
// 						}
// 					}else{
// 						if(j==k){
// 							perm=3.0;
// 						}else{
// 							perm=6.0;
// 						}
// 					}
//
// 					for(l=0;l<2;l++){
// 						for(m=0;m<2;m++){
// 							for(n=0;n<2;n++){
//
//  								k1 = factor*(i+l);
//  								k2 = factor*(j+m);
//  								k3 = factor*(k+n);
//
//  								ksum = k1+k2+k3;
//  								if (ksum>1e-10) ksum = 1e0/ksum;
// 								points[l][m][n] = pijk(r,i+l,j+m,k+n)*shape3(k1,k2,k3,params)*ksum;
// // 								points[l][m][n] = pijk(r,i+l,j+m,k+n)*pijk(1,i+l,j+m,k+n)*ksum;
// // 								points[l][m][n] = pijk(r,i+l,j+m,k+n)*shape3(k1,k2,k3);
//
// 							}
// 						}
// 					}
//
// 					integral += perm * calculate_volume(i,j,k,points) * cube_size;
//
//
// 				}
//
// 			}
// 		}
// 		// MPI sync
// // 		if ( myrank == 0 ) sync_tasks(5,2);
// 	}
// 	return integral;
// }
//
// double calculate_eigen_mc(int r){
//
// 	double res1, err1;
// 	double res2, err2;
// 	double res3, err3;
// 	double kmin = 2e0/get_tau0();
// 	double kmax = get_kmax();
// 	double xl[3] = {kmin,kmin,kmin};
// 	double xu[3] = {kmax,kmax,kmax};
//
// 	const gsl_rng_type *T;
// 	gsl_rng *rng;
//
// 	int p1,p2,p3;
// 	find_perm_prim(r,&p1,&p2,&p3);
//
// 	struct my_int_params params = {p1,p2,p3};
// 	gsl_monte_function QS = { &integrand_QS, 3, &params};
//
// 	size_t calls = 20000000;
//
// 	gsl_rng_env_setup();
//
// 	T = gsl_rng_default;
// 	rng = gsl_rng_alloc(T);
//
// 	gsl_monte_plain_state *state = gsl_monte_plain_alloc(3);
// 	gsl_monte_plain_integrate (&QS, xl, xu, 3, calls, rng, state, &res1, &err1);
// 	gsl_monte_plain_free(state);
// 	gsl_rng_free (rng);
//
// 	return res1;
// }
//
// double integrand_QS(double *k, size_t dim, void *params){
// 	struct my_int_params *p = (struct my_int_params *)params;
//
// 	double k1 = k[0];
// 	double k2 = k[1];
// 	double k3 = k[2];
// 	int p1,p2,p3;
// 	p1 = p->p1;
// 	p2 = p->p2;
// 	p3 = p->p3;
//
// 	double kinv = 1e0/(k1+k2+k3);
//
// 	double result;
// 	double b1,b2,b3,b4,b5,b6;
// 	double part1,part2;
// 	double sixinv = 1e0/6e0;
//
// 	if(k1>k2+k3||k2>k3+k1||k3>k1+k2){
// 		result = 0e0;
// 	}else{
//  		b1 = sinlog_pt(k1,p1)*sinlog_pt(k2,p2)*sinlog_pt(k3,p3);
//  		b2 = sinlog_pt(k1,p2)*sinlog_pt(k2,p3)*sinlog_pt(k3,p1);
//  		b3 = sinlog_pt(k1,p3)*sinlog_pt(k2,p1)*sinlog_pt(k3,p2);
//  		b4 = sinlog_pt(k1,p3)*sinlog_pt(k2,p2)*sinlog_pt(k3,p1);
//  		b5 = sinlog_pt(k1,p2)*sinlog_pt(k2,p1)*sinlog_pt(k3,p3);
//  		b6 = sinlog_pt(k1,p1)*sinlog_pt(k2,p3)*sinlog_pt(k3,p2);
//  		part1 = (b1+b2+b3+b4+b5+b6)*sixinv;
//
//  		// part2 = shape3(k1,k2,k3);
//
// 		result = kinv*part1*part2;
// 	}
// 	return result;
// }
//
// double check(shape_params params){
//  	double integral1=0;
//  	double integral2=0;
//  	double integral3=0;
// 	double ***points1 = create_3Darray(2,2,2);
// 	double ***points2 = create_3Darray(2,2,2);
// 	double ***points3 = create_3Darray(2,2,2);
//  	double cube_size = pow((double)alpha_points,-3);
//
// 	int p1,p2,p3;
//
// 	double k1,k2,k3,ksum;
// 	double max = get_kmax();
//
// 	int i,j,k,l,m,n;
//
// 	for(i=0;i<alpha_points;i++){
// 		for(j=0;j<alpha_points;j++){
// 			for(k=0;k<alpha_points;k++){
// 				if(i<=j+k+1 && j<=k+i+1 && k<=i+j+1){
//
// 					for(l=0;l<2;l++){
// 						for(m=0;m<2;m++){
// 							for(n=0;n<2;n++){
//  								k1 = max*(double)(i+l)/alpha_points;
//  								k2 = max*(double)(j+m)/alpha_points;
//  								k3 = max*(double)(k+n)/alpha_points;
//
//  								ksum = k1+k2+k3;
//  								if (ksum!=0) ksum = 1.0/ksum;
// // 								ksum = 1.0;
//
// 								points1[l][m][n] = sijk(i+l,j+m,k+n)*shape3(k1,k2,k3,params)*ksum;
// 								points2[l][m][n] = sijk(i+l,j+m,k+n)*sijk(i+l,j+m,k+n)*ksum;
// 								points3[l][m][n] = shape3(k1,k2,k3,params)*shape3(k1,k2,k3,params)*ksum;
// 							}
// 						}
// 					}
//
// 					integral1 += calculate_volume(i,j,k,points1) * cube_size;
// 					integral2 += calculate_volume(i,j,k,points2) * cube_size;
// 					integral3 += calculate_volume(i,j,k,points3) * cube_size;
//
// 				}
//
// 			}
// 		}
// 	}
// // 	printf("%e\t%e\t%e\n",integral1,integral2,integral3);
// 	return integral1 / sqrt(integral2*integral3);
//
// }
//
double checkmode(int num, shape_params params){
 	double integral1=0;
  	double integral2=0;
  	double integral3=0;
 	double ***points1 = create_3Darray(2,2,2);
 	double ***points2 = create_3Darray(2,2,2);
 	double ***points3 = create_3Darray(2,2,2);
  	double cube_size = pow((double)alpha_points,-3);
  	double recon,factor;
 	int terms=get_terms_prim();
 	double* eigen = (double *)create_vector(num+1);
 	int p1,p2,p3;
 	double k1,k2,k3,ksum,x;
 	double max = get_kmax();

 	int i,j,k,l,m,n,p;

 	for(i=0;i<num+1;i++){
 		eigen[i] = 0.0;
 		//for( j=i;j<num+1;j++){
 		//	eigen[i] += get_lambda(j,i)*get_eigenR(j);
 		//}
		eigen[i] = get_eigen(i);
 	}
	printf("Checkmode\n");
 	for(i=0;i<alpha_points;i++){
 		for(j=0;j<alpha_points;j++){
 			for(k=0;k<alpha_points;k++){
 				if(i<=j+k+1 && j<=k+i+1 && k<=i+j+1){

 					for(l=0;l<2;l++){
 						for(m=0;m<2;m++){
 							for(n=0;n<2;n++){
 								recon = 0.0;
 								for(p=0;p<num+1;p++){
 									recon += eigen[p]*pijk(p,i+l,j+m,k+n);
 								}

  								k1 = max*(double)(i+l)/alpha_points;
  								k2 = max*(double)(j+m)/alpha_points;
  								k3 = max*(double)(k+n)/alpha_points;
  								ksum = k1+k2+k3;
  								if (ksum > 1e-15) ksum = 1.0/ksum;
  								//ksum =1.0;
 								points1[l][m][n] = recon*shape3(k1,k2,k3,params)*ksum;
 								points2[l][m][n] = recon*recon*ksum;
 								points3[l][m][n] = shape3(k1,k2,k3,params)*shape3(k1,k2,k3,params)*ksum;
 							}
 						}
 					}

 					integral1 += calculate_volume(i,j,k,points1) * cube_size;
 					integral2 += calculate_volume(i,j,k,points2) * cube_size;
 					integral3 += calculate_volume(i,j,k,points3) * cube_size;

 				}

 			}
 		}
 	}
 	return integral1 / sqrt(integral2*integral3);
}

//
// double correlation_prim(shape_params params1, shape_params params2){
//  	double integral1=0;
//  	double integral2=0;
//  	double integral3=0;
// 	double ***points1 = create_3Darray(2,2,2);
// 	double ***points2 = create_3Darray(2,2,2);
// 	double ***points3 = create_3Darray(2,2,2);
//  	double cube_size = pow((double)alpha_points,-3);
//
// 	int p1,p2,p3;
//
// 	double k1,k2,k3,ksum;
// 	double max = get_kmax();
//
// 	int i,j,k,l,m,n;
//
// 	for(i=0;i<alpha_points;i++){
// 		for(j=0;j<alpha_points;j++){
// 			for(k=0;k<alpha_points;k++){
// 				if(i<=j+k+1 && j<=k+i+1 && k<=i+j+1){
//
// 					for(l=0;l<2;l++){
// 						for(m=0;m<2;m++){
// 							for(n=0;n<2;n++){
//  								k1 = max*(double)(i+l)/alpha_points;
//  								k2 = max*(double)(j+m)/alpha_points;
//  								k3 = max*(double)(k+n)/alpha_points;
//
//  								ksum = k1+k2+k3;
//  								if (ksum!=0) ksum = 1.0/ksum;
// // 								ksum = 1.0;
//  								ksum *= k1*k1*k2*k2*k3*k3*k1*k1*k2*k2*k3*k3;
// 								points1[l][m][n] = (single(k1,k2,k3,params1)+nicola529(k1,k2,k3,params1))*equilateral(k1,k2,k3,params2)*ksum;
// 								points2[l][m][n] = (single(k1,k2,k3,params1)+nicola529(k1,k2,k3,params1))*(single(k1,k2,k3,params1)+nicola529(k1,k2,k3,params1))*ksum;
// 								points3[l][m][n] = equilateral(k1,k2,k3,params2)*equilateral(k1,k2,k3,params2)*ksum;
// 							}
// 						}
// 					}
//
// 					integral1 += calculate_volume(i,j,k,points1) * cube_size;
// 					integral2 += calculate_volume(i,j,k,points2) * cube_size;
// 					integral3 += calculate_volume(i,j,k,points3) * cube_size;
//
// 				}
//
// 			}
// 		}
// 	}
// // 	printf("%e\t%e\t%e\n",integral1,integral2,integral3);
// 	return integral1 / sqrt(integral2*integral3);
//
// }
