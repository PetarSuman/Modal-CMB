#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "global.h"

static int npt;
static double *lweight;
static double *gl_pts;
static double *gl_wgt;
static double **gl_Pl;
static double ***beta;

double init_norm_glint(){

	int lsize = get_b_lsize();
	int *lvec = create_ivector(lsize);
	get_b_lvec(lvec);
	int xsize = get_b_xsize();
	int lmax = lvec[lsize-1];
	npt = 3*lmax/2 + 1;
	
	lweight = (double *)create_vector(lsize);
	double s1,s2,s3;
		
	int i,j,l;
	
	for(l=0;l<lsize;l++){
		if(get_cl(l)>0.0 && lvec[l] > 1){
			s1 = get_beam(l)*get_beam(l)*(2.0*lvec[l]+1.0);
			s2 = get_beam(l)*get_beam(l)*get_cl(l)+get_noise(l);
			lweight[l] =  s1/s2;
		} else {
			lweight[l] = 0.0;
		}
	}
	gsl_integration_glfixed_table *table =  gsl_integration_glfixed_table_alloc(npt);
	
	gl_pts = (double *)create_vector(npt);
	gl_wgt = (double *)create_vector(npt);
	
	for(i=0;i<npt;i++){
		gsl_integration_glfixed_point (-1.0, 1.0, i, &gl_pts[i], &gl_wgt[i], table);
	}
	
	gl_Pl = (double **)create_array(npt,lsize);
	for(l=0;l<lsize;l++){
		for(i=0;i<npt;i++){
			gl_Pl[i][l] = gsl_sf_legendre_Pl(lvec[l],gl_pts[i]);
		}
	}
	
	int pmax = get_pmax_prim();
	beta = (double ***)create_3Darray(xsize,pmax+1,lsize);
	for(l=0;l<lsize;l++){
		for(i=0;i<xsize;i++){
			for(j=0;j<pmax+1;j++){
				beta[i][j][l] = get_beta(l,j,i);
			}
		}
	}
	
	destroy_beta();
	
	return 0;
}

double norm_pt(int n, int i, int j){

	int l,r,s,k;
	int pmax = get_pmax_prim();
	int lsize = get_b_lsize();
	double *factor = create_vector(lsize);
	double **Nmap = create_array(3,3);
	int p1,p2,p3,q1,q2,q3;
	double sum1;
	double s1,s2,s3,s4,s5,s6;
	
	double result = 0.0;

	find_perm_late(n,&p1,&p2,&p3);

	int pvec[3];

	pvec[0] = p1;
	pvec[1] = p2;
	pvec[2] = p3;

	for(k=0;k<npt;k++){
		
		for(l=0;l<lsize;l++){
			factor[l] =  lweight[l]*gl_Pl[k][l];
		}
		
		for(r=0;r<3;r++){
			for(s=0;s<3;s++){
				sum1 = 0.0;
				for(l=0;l<lsize;l++){
					sum1 += factor[l]*beta[i][pvec[r]][l]*beta[j][pvec[s]][l];
				}
				Nmap[r][s] = sum1;
			}
		}


		s1  = Nmap[0][0]*Nmap[1][1]*Nmap[2][2];
		s2  = Nmap[0][0]*Nmap[1][2]*Nmap[2][1];
		s3  = Nmap[0][1]*Nmap[1][0]*Nmap[2][2];
		s4  = Nmap[0][1]*Nmap[1][2]*Nmap[2][0];
		s5  = Nmap[0][2]*Nmap[1][0]*Nmap[2][1];
		s6  = Nmap[0][2]*Nmap[1][1]*Nmap[2][0];
		
		sum1 = s1+s2+s3+s4+s5+s6;
		result += gl_wgt[k]*sum1;
	}
	
	destroy_array(Nmap);
	destroy_vector(factor);
	return result;
}

double calculate_norm(int n, int i) {

	double result = 0.0;
	int j;
	
	int xsize = get_b_xsize();
	double *xvec = create_vector(xsize);
	get_b_xvec(xvec);
	double *yvec = create_vector(xsize);

	double xmin = xvec[0];
	double xmax = xvec[xsize-1];
	
	for(j=0;j<xsize;j++){
		yvec[j] = xvec[j]*xvec[j]*norm_pt(n,i,j);
	}
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, xsize);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	gsl_spline_init(sp,xvec,yvec,xsize);
	result = gsl_spline_eval_integ(sp,xmin,xmax,acc);
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
	result *= 3.0*xvec[i]*xvec[i]*pow(1.55e-8,4)/(4.0*M_PI);
	
 	return result;

}