#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "global.h"

// Declare external functions
static int npt;
static double *gl_pts;
static double *gl_wgt;
static double **gl_Pl;
static double *lweight;
static double **basis;
static double **power;

double big_b(int i,int j,int k);
double little_b(int i);

double init_PS_glint(){

	int lmax = eflag_T_lmax;
	npt = 3*lmax/2 + 1;
		
	int i,j,l;
	
	gsl_integration_glfixed_table *table =  gsl_integration_glfixed_table_alloc(npt);
	
	gl_pts = (double *)create_vector(npt);
	gl_wgt = (double *)create_vector(npt);
	
	for(i=0;i<npt;i++){
		gsl_integration_glfixed_point (-1.0, 1.0, i, &gl_pts[i], &gl_wgt[i], table);
	}
	
	gl_Pl = (double **)create_array(npt,lmax+1);
	
	for(l=0;l<lmax+1;l++){
		for(i=0;i<npt;i++){
			gl_Pl[i][l] = gsl_sf_legendre_Pl(l,gl_pts[i]);
		}
	}
		
	lweight = (double *)create_vector(lmax+1);
	double s1,s2,s3;
		
	
	for(l=0;l<lmax+1;l++){
		if(l > 1){
			s1 = get_beam_TT(l)*pow(2.0*l+1.0,5.0/6.0);
			s2 = get_beam_TT(l)*get_beam_TT(l)*get_cl_TT(l)+get_noise_TT(l);
			s3 = 1.0 / sqrt(s2);
			lweight[l] =  s1*s3;
		} else {
			lweight[l] = 0.0;
		}
	}

	int pmax = get_pmax_late_T();
	basis = (double **)create_array(pmax+1,lmax+1);
	for(l=0;l<eflag_T_lmin;l++){
		for(i=0;i<pmax+1;i++){
			basis[i][l] = 0e0;
		}
	}
	for(l=eflag_T_lmin;l<lmax+1;l++){
		for(i=0;i<pmax+1;i++){
			basis[i][l] = get_basis_late_T(l-eflag_T_lmin,i);
		}
	}

	return 0;
}

double calculate_PS(int n) {

	int l,r,s,j;
	int pmax = get_pmax_prim();
	int lsize = eflag_T_lmax+1;
	double *factor = create_vector(lsize);
	double **Nmap = create_array(3,1);
	int p1,p2,p3,q1,q2,q3;
	double s1,s2,s3;
	double sum1;
	double sum2;
	
	double result = 0.0;
	double *scale = create_vector(lsize);

	find_perm_late_TTT(n,&p1,&p2,&p3);

	int pvec[3];

	pvec[0] = p1;
	pvec[1] = p2;
	pvec[2] = p3;
		
	for(l=0;l<lsize;l++){
		if(l > 1){
			s1 = pow(2.0*l+1.0,5.0/6.0);
			s2 = get_cl_TT(l)+get_noise_TT(l)/(get_beam_TT(l)*get_beam_TT(l));
			s3 = 1.0/sqrt(s2);
			scale[l] =  s1*s3;
		} else {
			scale[l] = 0.0;
		}
	}
	
	for(j=0;j<npt;j++){
		
		for(l=0;l<lsize;l++){
			factor[l] =  scale[l]*gl_Pl[j][l];
			
		}
		
		for(r=0;r<3;r++){
			sum1 = 0.0;
			for(l=2;l<lsize;l++){
				sum1 += factor[l]*basis[pvec[r]][l]*little_b(l);
			}
			Nmap[r][0] = sum1;
		}
		
		sum2  = Nmap[0][0]*Nmap[1][0]*Nmap[2][0];

		result += gl_wgt[j]*sum2;
	}
	result /= (8.0*M_PI);
	
	destroy_array(Nmap);
	destroy_vector(factor);
	return result;
}

double calculate_PS_norm(){

	int l,r,s,j,t1,t2,t3;
	int pmax = get_pmax_prim();
	int lsize = eflag_T_lmax+1;
	double *factor = create_vector(lsize);
	
	int p1,p2,p3,q1,q2,q3;
	double sum1;
	double sum2;
	double s1,s2,s3;
	
	double result = 0.0;
	
	double *scale = create_vector(lsize);
	
		
	for(l=0;l<lsize;l++){
		if(l > 1){
			s1 = (2.0*l+1.0);
			s2 = get_cl_TT(l)+get_noise_TT(l)/(get_beam_TT(l)*get_beam_TT(l));
			s3 = 1.0/s2;
			scale[l] =  s1*s3;
		} else {
			scale[l] = 0.0;
		}
	}
	
	for(j=0;j<npt;j++){
		
		for(l=0;l<lsize;l++){
			factor[l] =  scale[l]*gl_Pl[j][l];
			
		}
				
		sum1 = 0.0;
		for(l=2;l<lsize;l++){
			sum1 += factor[l]*little_b(l)*little_b(l);
		}
		
		result += gl_wgt[j]*sum1*sum1*sum1;
	}
	
	result *= 1.0/(8.0*M_PI);
	destroy_vector(factor);
	return result;
}

double calculate_PS_norm_3D(){

	int i,j,k,t1,t2,t3;
	double x;
	int lsize = eflag_T_lmax+1;
	int *lvec = create_ivector(lsize);
	for(i=0;i<lsize;i++){
		lvec[i] = i;
	}
	double result = 0.0;
	double s1,s2,s3;
	double y;
// 	int test;
		
	for(i=2;i<lsize;i++){
		s1 = get_cl_TT(i)+get_noise_TT(i)/(get_beam_TT(i)*get_beam_TT(i));
		for(j=i;j<lsize;j++){
			s2 = get_cl_TT(j)+get_noise_TT(j)/(get_beam_TT(j)*get_beam_TT(j));
			t1 = i+j;
			t2 = j%2;
			t3 = t1%2;
			if(t1>eflag_T_lmax)t1=eflag_T_lmax;
			if(t2==0&&t3==0){
				for(k=j;k<t1+1;k+=2){
					y = big_b(i,j,k);
					s3 = get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k));
					result += permsix(i,j,k)*calculate_geometric(i,j,k)*y*y/(s1 * s2 * s3);
				}
			}else if(t2==0&&t3==1){
				for(k=j+1;k<t1+1;k+=2){
					y = big_b(i,j,k);
					s3 = get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k));
					result += permsix(i,j,k)*calculate_geometric(i,j,k)*y*y/(s1 * s2 * s3);
				}
			}else if(t2==1&&t3==0){
				for(k=j+1;k<t1+1;k+=2){
					y = big_b(i,j,k);
					s3 = get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k));
					result += permsix(i,j,k)*calculate_geometric(i,j,k)*y*y/(s1 * s2 * s3);
				}
			}else if(t2==1&&t3==1){
				for(k=j;k<t1+1;k+=2){
					y = big_b(i,j,k);
					s3 = get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k));
					result += permsix(i,j,k)*calculate_geometric(i,j,k)*y*y/(s1 * s2 * s3);
				}
			}
		}
	}
	return result;
}

// double calculate_PS_3D(int n){
//
// 	int i,j,k,t1,t2,t3;
// 	double x,y;
// 	int lsize = eflag_T_lmax+1;
// 	int *lvec = create_ivector(lsize);
// 	for(i=0;i<lsize;i++){
// 		lvec[i] = i;
// 	}
// 	double result = 0.0;
// 	double s1,s2,s3;
//
// 	for(i=2;i<lsize;i++){
// 		s1 = pow(2.0*i+1.0,1.0/3.0)*(get_cl_TT(i)+get_noise_TT(i)/(get_beam_TT(i)*get_beam_TT(i)));
// 		for(j=i;j<lsize;j++){
// 			s2 = pow(2.0*j+1.0,1.0/3.0)*(get_cl_TT(j)+get_noise_TT(j)/(get_beam_TT(j)*get_beam_TT(j)));
// 			t1 = i+j;
// 			t2 = j%2;
// 			t3 = t1%2;
// 			if(t1>eflag_T_lmax)t1=eflag_T_lmax;
// 			if(t2==0&&t3==0){
// 				for(k=j;k<t1+1;k+=2){
// 					x = big_b(i,j,k);
// 					y = plijk_TTT(n,i,j,k);
// 					s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k)));
// 					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*y/sqrt(s1 * s2 * s3);
// 				}
// 			}else if(t2==0&&t3==1){
// 				for(k=j+1;k<t1+1;k+=2){
// 					x = big_b(i,j,k);
// 					y = plijk_TTT(n,i,j,k);
// 					s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k)));
// 					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*y/sqrt(s1 * s2 * s3);
// 				}
// 			}else if(t2==1&&t3==0){
// 				for(k=j+1;k<t1+1;k+=2){
// 					x = big_b(i,j,k);
// 					y = plijk_TTT(n,i,j,k);
// 					s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k)));
// 					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*y/sqrt(s1 * s2 * s3);
// 				}
// 			}else if(t2==1&&t3==1){
// 				for(k=j;k<t1+1;k+=2){
// 					x = big_b(i,j,k);
// 					y = plijk_TTT(n,i,j,k);
// 					s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k)));
// 					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*y/sqrt(s1 * s2 * s3);
// 				}
// 			}
// 		}
// 	}
// 	return result;
// }

double big_b(int i,int j,int k){
	return little_b(i)*little_b(j)*little_b(k);
}

double little_b(int l){
// 	double s1,s2,s3;
// 	s1 = get_beam(l)*pow(2.0*l+1.0,1.0/6.0);
// 	s2 = get_beam(l)*get_beam(l)*get_cl(l)+get_noise(l);
// 	s3 = sqrt(s2)/s1;
	
// 	return s3;
	return 1e0;
}