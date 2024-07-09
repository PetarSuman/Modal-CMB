#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "global.h"

static int npt;
static double *gl_pts;
static double *gl_wgt;
static double **gl_Pl;
static double **basis;

double init_orthol_tri_glint(){

	int lmax = get_lmax();
	npt = 3*lmax/2 + 1;
		
	int i,l;
	
	gsl_integration_glfixed_table *table =  gsl_integration_glfixed_table_alloc(npt);
	
	gl_pts = (double *)create_vector(npt);
	gl_wgt = (double *)create_vector(npt);
	
//	for(i=0;i<npt;i++){
//		gsl_integration_glfixed_point (-1.0, 1.0, i, &gl_pts[i], &gl_wgt[i], table);
//	}
	asy(gl_pts,gl_wgt,npt);	
	
	gl_Pl = (double **)create_array(npt,lmax+1);
	
	for(l=0;l<lmax+1;l++){
		for(i=0;i<npt;i++){
			gl_Pl[i][l] = gsl_sf_legendre_Pl(l,gl_pts[i]);
		}
	}

	int pmax = get_pmax_late();
	printf("pmax = %d\n",pmax);
	basis = (double **)create_array(pmax+1,lmax+1);
	for(l=0;l<lmax+1;l++){
		for(i=0;i<pmax+1;i++){
			basis[i][l] = get_basis_tri_late(l,i);
		}
	}

	return 0;
}


double calculate_orthol_tri(int m, int n) {

	int l,i,r,s,a,b;
	int pmax = get_pmax_late();
	int lsize = get_lmax()+1;
	double *factor = create_vector(lsize);
	double **Nmap = create_array(4,4);
	int p1,p2,p3,p4,q1,q2,q3,q4;
	double sum1;
	double s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24;
	
	double result = 0.0;

	find_perm_tri_late(m,&p1,&p2,&p3,&p4);
	find_perm_tri_late(n,&q1,&q2,&q3,&q4);

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
	
	for(i=0;i<npt;i++){
		for(l=0;l<lsize;l++){
			s1 = sqrt(2e0*l+1e0);
			factor[l] =  s1*gl_Pl[i][l];
		}

		for(r=0;r<4;r++){
			a = pvec[r];
// 			printf("%d\t%d\t%d\t%d\n",m,n,i,a);
			for(s=r;s<4;s++){
				b = qvec[s];
				sum1 = 0e0;
				for(l=2;l<lsize;l++){
					sum1 += factor[l]*basis[a][l]*basis[b][l];
				}
				Nmap[r][s] = sum1;
				Nmap[s][r] = sum1;
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
		result += gl_wgt[i]*sum1/(768e0*M_PI*M_PI);
	}
	destroy_array(Nmap);
	destroy_vector(factor);

	return result;
}

double trilR(int r,int i,int j,int k,int l){

	int n;
	double result = 0;
	
	for(n=0;n<r+1;n++){
		result += get_orthol_tri(r,n)*triQ(n,i,j,k,l);
	}
 	
	return result;
}