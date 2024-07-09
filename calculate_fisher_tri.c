#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "global.h"

struct fish_params {int i; int j; int m; int n;};

double fisher_tri_pt(double x, void *p){

	struct fish_params *params = (struct fish_params *)p;
	int i = (params->i);
	int j = (params->j);
	int m = (params->m);
	int n = (params->n);

	int l,r,s;
	int qmax = get_qmax();
	int lsize = get_bt_lsize();
	int *lvec = create_ivector(lsize);
	double *factor = create_vector(lsize);
	get_bt_lvec(lvec);
	double **Nmap = create_array(qmax+1,qmax+1);
	int p1,p2,p3,p4,q1,q2,q3,q4;
	double sum1;
	double s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24;
	
	double result = 0.0;
	
	for(l=0;l<lsize;l++){
	
		printf("%d\t%d\n",l,lvec[l]);
		
		if(get_cl(l)>0.0){
			factor[l] = (2.0*lvec[l]+1.0)*gsl_sf_legendre_Pl(lvec[l],x)/(sqrt(4.0*M_PI)*get_cl(l));
		} else {
			factor[l] = 0.0;
		}
	}

	for(r=0;r<qmax+1;r++){
		for(s=r;s<qmax+1;s++){
			sum1 = 0.0;
			for(l=2;l<lsize;l++){
				sum1 += factor[l]*get_beta_tri(r,l,i)*get_beta_tri(s,l,j);
			}
			Nmap[r][s] = sum1;
			Nmap[s][r] = sum1;
		}
	}

	find_perm_tri(m,&p1,&p2,&p3,&p4);
	find_perm_tri(n,&q1,&q2,&q3,&q4);
	
	s1  = Nmap[p1][q1]*Nmap[p2][q2]*Nmap[p3][q3]*Nmap[p4][q4];
	s2  = Nmap[p1][q1]*Nmap[p2][q2]*Nmap[p3][q4]*Nmap[p4][q3];
	s3  = Nmap[p1][q1]*Nmap[p2][q3]*Nmap[p3][q2]*Nmap[p4][q4];
	s4  = Nmap[p1][q1]*Nmap[p2][q3]*Nmap[p3][q4]*Nmap[p4][q2];
	s5  = Nmap[p1][q1]*Nmap[p2][q4]*Nmap[p3][q2]*Nmap[p4][q3];
	s6  = Nmap[p1][q1]*Nmap[p2][q4]*Nmap[p3][q3]*Nmap[p4][q2];
	
	s7  = Nmap[p1][q2]*Nmap[p2][q1]*Nmap[p3][q3]*Nmap[p4][q4];
	s8  = Nmap[p1][q2]*Nmap[p2][q1]*Nmap[p3][q4]*Nmap[p4][q3];
	s9  = Nmap[p1][q2]*Nmap[p2][q3]*Nmap[p3][q1]*Nmap[p4][q4];
	s10 = Nmap[p1][q2]*Nmap[p2][q3]*Nmap[p3][q4]*Nmap[p4][q1];
	s11 = Nmap[p1][q2]*Nmap[p2][q4]*Nmap[p3][q1]*Nmap[p4][q3];
	s12 = Nmap[p1][q2]*Nmap[p2][q4]*Nmap[p3][q3]*Nmap[p4][q1];
	
	s13 = Nmap[p1][q3]*Nmap[p2][q1]*Nmap[p3][q2]*Nmap[p4][q4];
	s14 = Nmap[p1][q3]*Nmap[p2][q1]*Nmap[p3][q4]*Nmap[p4][q2];
	s15 = Nmap[p1][q3]*Nmap[p2][q2]*Nmap[p3][q1]*Nmap[p4][q4];
	s16 = Nmap[p1][q3]*Nmap[p2][q2]*Nmap[p3][q4]*Nmap[p4][q1];
	s17 = Nmap[p1][q3]*Nmap[p2][q4]*Nmap[p3][q1]*Nmap[p4][q2];
	s18 = Nmap[p1][q3]*Nmap[p2][q4]*Nmap[p3][q2]*Nmap[p4][q1];
	
	s19 = Nmap[p1][q4]*Nmap[p2][q1]*Nmap[p3][q2]*Nmap[p4][q3];
	s20 = Nmap[p1][q4]*Nmap[p2][q1]*Nmap[p3][q3]*Nmap[p4][q2];
	s21 = Nmap[p1][q4]*Nmap[p2][q2]*Nmap[p3][q1]*Nmap[p4][q3];
	s22 = Nmap[p1][q4]*Nmap[p2][q2]*Nmap[p3][q3]*Nmap[p4][q1];
	s23 = Nmap[p1][q4]*Nmap[p2][q3]*Nmap[p3][q1]*Nmap[p4][q2];
	s24 = Nmap[p1][q4]*Nmap[p2][q3]*Nmap[p3][q2]*Nmap[p4][q1];
	
	sum1 =s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19+s20+s21+s22+s23+s24;
	result = sum1/1152.0;
	
	destroy_array(Nmap,qmax);
	destroy_vector(factor);
	destroy_ivector(lvec);
	return result;
}

double calculate_fisher_tri(int i, int j) {

	double result = 0.0;
	int k,m,n;
	int lsize = get_bt_lsize();
	int *lvec = create_ivector(lsize);
	get_bt_lvec(lvec);
	int lmax = lvec[lsize-1];
	int npt = 2*lmax + 1;
	int terms = get_terms_tri();
	
	gsl_function F;
	struct fish_params params;
	F.function = &fisher_tri_pt;
	F.params = &params;
	
	params.i = i;
	params.j = j;
	
	gsl_integration_glfixed_table *table =  gsl_integration_glfixed_table_alloc(npt);
	for(m=0;m<terms;m++){
		for(n=m;m<terms;m++){
			params.m = m;
			params.n = n;
			
			if(m==n){
				result += get_eigen_tri(m)*get_eigen_tri(n)*gsl_integration_glfixed(&F, -1.0, 1.0, table);
			}else{
				result += 2.0*get_eigen_tri(m)*get_eigen_tri(n)*gsl_integration_glfixed(&F, -1.0, 1.0, table);
			}
		}
	}
	
	result *= 9.0*pow(1.53E-8,6);
 	return result;

}