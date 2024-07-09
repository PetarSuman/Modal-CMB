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
static double **basis;
static double ***beta;

double init_gamma_tri_glint(){

	int lsize = get_bt_lsize();
	int *lvec = create_ivector(lsize);
	get_bt_lvec(lvec);
	int xsize = get_bt_xsize();
	int lmax = lvec[lsize-1];
	npt = 3*lmax/2 + 1;
	
	lweight = (double *)create_vector(lsize);
	double s1,s2,s3;	
	int i,j,l;
	
	for(i=0;i<lsize;i++){
		l = lvec[i];
		if(l > 1){
			s1 = pow(2.0*l+1.0,3.0/4.0);
			s2 = get_cl(l)+get_noise(l)/(get_beam(l)*get_beam(l));
			s3 = 1.0 / sqrt(s2);
			lweight[l] =  s1*s3;
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

	asy(gl_pts,gl_wgt,npt);

	gl_Pl = (double **)create_array(npt,lsize);
	
	for(l=0;l<lsize;l++){
		for(i=0;i<npt;i++){
			gl_Pl[i][l] = gsl_sf_legendre_Pl(lvec[l],gl_pts[i]);
		}
	}

	int pmax = get_pmax_late();

	basis = (double **)create_array(pmax+1,lsize);
	for(l=0;l<lsize;l++){
		for(i=0;i<pmax+1;i++){
			basis[i][l] = get_basis_tri_late(lvec[l],i);
		}
	}
	
	
	pmax = get_pmax_prim();
	beta = (double ***)create_3Darray(xsize,pmax+1,lsize);
	for(l=0;l<lsize;l++){
		for(i=0;i<xsize;i++){
			for(j=0;j<pmax+1;j++){
				beta[i][j][l] = get_beta_tri(l,j,i);
			}
		}
	}	
	
	destroy_beta();
	
	return 0;
}

double gamma_tri_pt(int m, int n, int i){

	int l,r,s,j;
	int pmax = get_pmax_prim();
	int lsize = get_bt_lsize();
	double *factor = create_vector(lsize);
	double **Nmap = create_array(4,4);
	int p1,p2,p3,p4,q1,q2,q3,q4;
	double sum1;
	double s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24;
	
	double result = 0.0;

	find_perm_tri_late(m,&p1,&p2,&p3,&p4);
	find_perm_tri_prim(n,&q1,&q2,&q3,&q4);

	int pvec[4];
	int qvec[4];

	pvec[0] = p1;
	pvec[1] = p2;
	pvec[2] = p3;
	pvec[3] = p4;

	qvec[0] = q1;
	qvec[1] = q2;
	qvec[2] = q3;
	qvec[3] = q4;
	
	for(j=0;j<npt;j++){
		
		for(l=0;l<lsize;l++){
			factor[l] =  lweight[l]*gl_Pl[j][l];
		}
		
		for(r=0;r<4;r++){
			for(s=0;s<4;s++){
				sum1 = 0.0;
				for(l=0;l<lsize;l++){
					sum1 += factor[l]*basis[pvec[r]][l]*beta[i][qvec[s]][l];
				}
				Nmap[r][s] = sum1;
			}
		}
		
		s1  = Nmap[0][0]*Nmap[1][1]*Nmap[2][2]*Nmap[3][3];
		s2  = Nmap[0][0]*Nmap[1][1]*Nmap[2][3]*Nmap[3][2];
		s3  = Nmap[0][0]*Nmap[1][2]*Nmap[2][1]*Nmap[3][3];
		s4  = Nmap[0][0]*Nmap[1][2]*Nmap[2][3]*Nmap[3][1];
		s5  = Nmap[0][0]*Nmap[1][3]*Nmap[2][1]*Nmap[3][2];
		s6  = Nmap[0][0]*Nmap[1][3]*Nmap[2][2]*Nmap[3][1];
		
		s7  = Nmap[0][1]*Nmap[1][0]*Nmap[2][2]*Nmap[3][3];
		s8  = Nmap[0][1]*Nmap[1][0]*Nmap[2][3]*Nmap[3][2];
		s9  = Nmap[0][1]*Nmap[1][2]*Nmap[2][0]*Nmap[3][3];
		s10 = Nmap[0][1]*Nmap[1][2]*Nmap[2][3]*Nmap[3][0];
		s11 = Nmap[0][1]*Nmap[1][3]*Nmap[2][0]*Nmap[3][2];
		s12 = Nmap[0][1]*Nmap[1][3]*Nmap[2][2]*Nmap[3][0];
		
		s13 = Nmap[0][2]*Nmap[1][0]*Nmap[2][1]*Nmap[3][3];
		s14 = Nmap[0][2]*Nmap[1][0]*Nmap[2][3]*Nmap[3][1];
		s15 = Nmap[0][2]*Nmap[1][1]*Nmap[2][0]*Nmap[3][3];
		s16 = Nmap[0][2]*Nmap[1][1]*Nmap[2][3]*Nmap[3][0];
		s17 = Nmap[0][2]*Nmap[1][3]*Nmap[2][0]*Nmap[3][1];
		s18 = Nmap[0][2]*Nmap[1][3]*Nmap[2][1]*Nmap[3][0];
		
		s19 = Nmap[0][3]*Nmap[1][0]*Nmap[2][1]*Nmap[3][2];
		s20 = Nmap[0][3]*Nmap[1][0]*Nmap[2][2]*Nmap[3][1];
		s21 = Nmap[0][3]*Nmap[1][1]*Nmap[2][0]*Nmap[3][2];
		s22 = Nmap[0][3]*Nmap[1][1]*Nmap[2][2]*Nmap[3][0];
		s23 = Nmap[0][3]*Nmap[1][2]*Nmap[2][0]*Nmap[3][1];
		s24 = Nmap[0][3]*Nmap[1][2]*Nmap[2][1]*Nmap[3][0];
		
		sum1 =s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19+s20+s21+s22+s23+s24;
		result += gl_wgt[j]*sum1;
	}
	
	destroy_array(Nmap);
	destroy_vector(factor);
	return result;
}

double calculate_gamma_tri(int m, int n) {

	double result = 0.0;
	int i;
	
	int xsize = get_bt_xsize();
	double *xvec = create_vector(xsize);
	get_bt_xvec(xvec);
	double *yvec = create_vector(xsize);

	double xmin = xvec[0];
	double xmax = xvec[xsize-1];
	
	for(i=0;i<xsize;i++){
		yvec[i] = xvec[i]*xvec[i]*gamma_tri_pt(m,n,i);
		if(i<10&&fabs(yvec[i])>1e-18) yvec[i]=0;
		
	}
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, xsize);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	gsl_spline_init(sp,xvec,yvec,xsize);
	result = gsl_spline_eval_integ(sp,xmin,xmax,acc);
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
	result *= deltaphi*deltaphi*deltaphi/(32.0*M_PI*M_PI);
	
 	return result;

}