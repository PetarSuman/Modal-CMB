#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "global.h"

struct my_int_params{int p1; int p2; int p3; int q1; int q2; int q3;};
double integrand_QQ(double *k, size_t dim, void *params);


void calculate_ortho(int min, int max, tetrapyd_limits* block, double* ortho){
	int m,n,r,s,t,i,j,k;
	double x,y;
	double k1,k2,k3,ksum;
	double grid = 1e0/((double)max);
	double kmax = get_kmax();
	double factor = kmax*grid;
 	double cube_size = factor*factor*factor;

	int terms = get_terms_prim();
	int vsize = (terms)*(terms+1)/2;
	
	/*
	double ortho_private[vsize];
	double basis[terms];
	*/

	double * ortho_private = (double *)malloc(vsize*sizeof(double));
	double * basis = (double *)malloc(terms*sizeof(double));

	//#pragma omp parallel default(none) private(m,n,r,s,t,i,j,k,x,y,k1,k2,k3,ksum,basis,ortho_private) shared(ortho, block, factor, min, max, cube_size, terms, vsize)
	{
		//t = omp_get_thread_num();
		t = 0;
		i = block[t].i_bgn;
		j = block[t].j_bgn;
		k = block[t].k_bgn;
		//printf("t i j k = %d %d %d %d\n", t, i, j, k);
        m = 0;
		for(r=0;r<terms;r++){
			basis[r] = 0e0;
            for(s=r;s<terms;s++){
	            ortho_private[m] = 0.0;
	            m++;
            }
		}
		for(n=0;n<block[t].loops;n++){
			x = calculate_weight(min, max, i, j, k)*cube_size;
			k1 = factor*(i);
			k2 = factor*(j);
			k3 = factor*(k);

			ksum = k1+k2+k3;
			if (ksum>1e-10) ksum = 1e0/ksum;

			x *= ksum;

			for(r=0;r<terms;r++){
				basis[r] = pijk(r,i,j,k);
			}
			//printf("check\n");
			m = 0;
			for(r=0;r<terms;r++){
				for(s=r;s<terms;s++){
					ortho_private[m] += x*basis[r]*basis[s];
					m++;
				}
			}
			get_ijk_next_prim(max,&i,&j,&k);
		}
		//printf("check here\n");
		
		//#pragma omp critical
		{
			for(r=0;r<vsize;r++){
				ortho[r] += ortho_private[r];
			}
		}
		
	}
	
	return;
}

// double calculate_ortho(int r,int s){
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
//
// 								points[l][m][n] = pijk(r,i+l,j+m,k+n)*pijk(s,i+l,j+m,k+n)*ksum;
// // 								points[l][m][n] = pijk(r,i+l,j+m,k+n)*pijk(s,i+l,j+m,k+n);
//
// 							}
// 						}
// 					}
//
// 					integral += perm*calculate_volume(i,j,k,points) * cube_size;
//
// 				}
//
// 			}
// 		}
// 		// MPI sync
// // 		if ( myrank == 0 ) sync_tasks(6,3);
// 	}
// 	return integral;
// }

double calculate_ortho_mc(int r,int s){

	double res1, err1;
	double res2, err2;
	double res3, err3;
	double kmin = 2e0/get_tau0();
	double kmax = get_kmax();
	double xl[3] = {kmin,kmin,kmin};
	double xu[3] = {kmax,kmax,kmax};
	
	const gsl_rng_type *T;
	gsl_rng *rng;
	
	int p1,p2,p3,q1,q2,q3;
	find_perm_prim(r,&p1,&p2,&p3);
	find_perm_prim(s,&q1,&q2,&q3);

	struct my_int_params params = {p1,p2,p3,q1,q2,q3};
	gsl_monte_function QQ = { &integrand_QQ, 3, &params};
	
	size_t calls = 20000000;
	
	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);
	
	gsl_monte_plain_state *state = gsl_monte_plain_alloc(3);
	gsl_monte_plain_integrate (&QQ, xl, xu, 3, calls, rng, state, &res1, &err1);
	gsl_monte_plain_free(state);
	gsl_rng_free(rng);
	
	return res1;
}

double integrand_QQ(double *k, size_t dim, void *params){
	struct my_int_params *p = (struct my_int_params *)params;

	double k1 = k[0];
	double k2 = k[1];
	double k3 = k[2];
	int p1,p2,p3,q1,q2,q3;
	p1 = p->p1;
	p2 = p->p2;
	p3 = p->p3;
	q1 = p->q1;
	q2 = p->q2;
	q3 = p->q3;
	
	double kinv = 1e0/(k1+k2+k3);
	
	double result;
	double b1,b2,b3,b4,b5,b6;
	double part1,part2;
	double sixinv = 1e0/6e0;
	
	if(k1>k2+k3||k2>k3+k1||k3>k1+k2){
		result = 0e0;
	}else{
 		b1 = sinlog_pt(k1,p1)*sinlog_pt(k2,p2)*sinlog_pt(k3,p3);
 		b2 = sinlog_pt(k1,p2)*sinlog_pt(k2,p3)*sinlog_pt(k3,p1);
 		b3 = sinlog_pt(k1,p3)*sinlog_pt(k2,p1)*sinlog_pt(k3,p2);
 		b4 = sinlog_pt(k1,p3)*sinlog_pt(k2,p2)*sinlog_pt(k3,p1);
 		b5 = sinlog_pt(k1,p2)*sinlog_pt(k2,p1)*sinlog_pt(k3,p3);
 		b6 = sinlog_pt(k1,p1)*sinlog_pt(k2,p3)*sinlog_pt(k3,p2);
 		part1 = (b1+b2+b3+b4+b5+b6)*sixinv;
 		
 		b1 = sinlog_pt(k1,q1)*sinlog_pt(k2,q2)*sinlog_pt(k3,q3);
 		b2 = sinlog_pt(k1,q2)*sinlog_pt(k2,q3)*sinlog_pt(k3,q1);
 		b3 = sinlog_pt(k1,q3)*sinlog_pt(k2,q1)*sinlog_pt(k3,q2);
 		b4 = sinlog_pt(k1,q3)*sinlog_pt(k2,q2)*sinlog_pt(k3,q1);
 		b5 = sinlog_pt(k1,q2)*sinlog_pt(k2,q1)*sinlog_pt(k3,q3);
 		b6 = sinlog_pt(k1,q1)*sinlog_pt(k2,q3)*sinlog_pt(k3,q2);
 		part2 = (b1+b2+b3+b4+b5+b6)*sixinv;
 		
		result = kinv*part1*part2;
	}
	return result;
}

double calculate_delta(int r,int s, int size){
 
 	double integral=0;
	double ***points = create_3Darray(2,2,2);
 	double cube_size = pow((double)(size-1),-3);
	
	int p1,p2,p3;

	int i,j,k,l,m,n;
	double k1,k2,k3,ksum;
	double max = get_kmax();

	for(i=0;i<size-1;i++){
		for(j=0;j<size-1;j++){
			for(k=0;k<size-1;k++){
				if(i<=j+k+1 && j<=k+i+1 && k<=i+j+1){
					
					for(l=0;l<2;l++){
						for(m=0;m<2;m++){
							for(n=0;n<2;n++){
							
 								k1 = max*(double)(i+l)/alpha_points;
 								k2 = max*(double)(j+m)/alpha_points;
 								k3 = max*(double)(k+n)/alpha_points;
 								
 								ksum = k1+k2+k3;
 								if (ksum!=0) ksum = 1.0/ksum;
 								
								points[l][m][n] = rijk(r,i+l,j+m,k+n)*rijk(s,i+l,j+m,k+n)*ksum;
							}
						}
					}
					
					integral += calculate_volume(i,j,k,points) * cube_size;
					
				}
				
			}
		}
	}
	
	return integral;
	
}

double pijk(int r,int i,int j,int k){
	
	int p1,p2,p3;
	find_perm_prim(r,&p1,&p2,&p3);
	
	double b1,b2,b3,b4,b5,b6;
 	b1 = get_basis_prim(i,p1)*get_basis_prim(j,p2)*get_basis_prim(k,p3);
 	b2 = get_basis_prim(i,p2)*get_basis_prim(j,p3)*get_basis_prim(k,p1);
 	b3 = get_basis_prim(i,p3)*get_basis_prim(j,p1)*get_basis_prim(k,p2);
 	b4 = get_basis_prim(i,p3)*get_basis_prim(j,p2)*get_basis_prim(k,p1);
 	b5 = get_basis_prim(i,p2)*get_basis_prim(j,p1)*get_basis_prim(k,p3);
 	b6 = get_basis_prim(i,p1)*get_basis_prim(j,p3)*get_basis_prim(k,p2);
	
 	double result = (b1+b2+b3+b4+b5+b6)/(6.0);
 	
	return result;
}

double rijk(int r,int i,int j,int k){

	int n;
	double result = 0;
	
	for(n=0;n<r+1;n++){
		result += get_lambda(r,n)*pijk(n,i,j,k);
	}
 	
	return result;
}

double sijk(int i,int j,int k){

	int n;
	double result = 0;
	int s=get_terms_prim();
	
	for(n=0;n<s;n++){
		result += get_eigen(n)*pijk(n,i,j,k);
	}
 	
	return result;
}

