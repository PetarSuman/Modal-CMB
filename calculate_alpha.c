#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <nagmk21.h>
#include "global.h"

double calculate_alpha(int rank, double *x, int i, int j, int k){

	int r,s,t,n,m;

	double kmax = get_kmax();
	int ksize = alpha_points + 1;
	
	int spl_size = ksize + 4;
	int wrk_size = 6 * ksize + 16;
	
	double k1,k2,k3,p1,p2,p3,p4;
	
	
	double y1[ksize];
	double y2[ksize];
	double y3[ksize];
	double spl_k[spl_size];
	double spl_c[spl_size];
	double wrk[wrk_size];
	int ifail;
	double factor, result=0;
	
	factor = (2*i+1)*(2*j+1)*(2*k+1)/(kmax * kmax * kmax);
	
	for (r=0;r<ksize;r++) {
		
		k1 = x[r];
		p1 = get_basis(r,i);
		
		for (s=0;s<ksize;s++) {
			
			k2 = x[s];
			p2 = get_basis(s,j);
			
			for (t=0;t<ksize;t++) {
				
				k3 = x[t];
				p3 = get_basis(t,k);
				
				p4 = shape3(k1,k2,k3);
				
				y3[t] = p4*p3;
			}
			
			e01baf_(&ksize,x,y3,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
			e02bdf_(&spl_size,spl_k,spl_c,&result,&ifail);
			
			y2[s] = result*p2;
		}
		
		e01baf_(&ksize,x,y2,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
		e02bdf_(&spl_size,spl_k,spl_c,&result,&ifail);
		
		y1[r] = result*p1;
	}
	
	e01baf_(&ksize,x,y1,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
	e02bdf_(&spl_size,spl_k,spl_c,&result,&ifail);
	
	result *= factor;
	
	return result;
}
