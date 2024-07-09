#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "global.h"

// Declare external functions
static int xsize;
static double *xvec;
static int terms;
static double *modes;
static int lsize;
static int *lvec;
static double *weight;
static double ***bispectrum;

struct record
	{
		int l1,l2,l3;
		float bisp;
	};

double get_weight(int i, int j, int k);

void init_3D(int mode){
	
	int i;
	xsize = get_qtilde_xsize();
    xvec = create_vector(xsize);
	get_qtilde_xvec(xvec);
	terms = get_terms_prim();

	if(mode>1000){
	    modes = create_vector(terms);
		model = mode-1000;

		read_eigen(0,0,0);
		for(i=0;i<terms;i++){
			modes[i] = get_eigen(i);
		}
	}
	
	lsize = eflag_T_lmax+1;
	lvec = create_ivector(lsize);
	for(i=0;i<lsize;i++){
		lvec[i] = i;
	}
	
	weight = create_vector(lsize);
	weight[0] = 0e0;
	weight[1] = 0e0;
	
	double x;
	for(i=2;i<lsize;i++){
		 x = (get_cl_TT(i)+get_noise_TT(i)/(get_beam_TT(i)*get_beam_TT(i)));
		 weight[i] = pow(2.0*i+1.0,1.0/6.0)/sqrt(x);
		 // printf("%d\t%e\n",i,weight[i]);
	}
	


	
	int L,L1,L2,L3;
	double twopi = sqrt(2e0/3.1415927);
	double thrd = 1e0/3e0;
	double sxth = 1e0/6e0;
	double p1,p2,p3,p4,p5,p6,p7,p8;
	// double bisp,factor;
	// bispectrum = create_3Darray(eflag_T_lmax+1,eflag_T_lmax+1,eflag_T_lmax+1);
	//
	// FILE *ptr_myfile;
	// long int size = 334831501;
	// struct record my_record[size];
	//
	// ptr_myfile=fopen("/fast/space/projects/planck/jf334.private/TSS_Bispectrum/Blll_TSS.unf","rb");
	// if (!ptr_myfile){
	// 	printf("Unable to open file!");
	// }
	// fread(&my_record, size*sizeof(struct record), 1, ptr_myfile);
	// fclose(ptr_myfile);
	//
	// // #pragma omp parallel for private(l1,l2,l3,bisp,L,L1,L2,L3,p1,p2,p3,p4,p5,p6,factor)
	// for(i=0;i<size;i++){
	// 	l1 = my_record[i].l1;
	// 	l2 = my_record[i].l2;
	// 	l3 = my_record[i].l3;
	// 	bisp = my_record[i].bisp;
	//
	// 	if(fmin(fmin(l1,l2),l3)<2 || fmax(fmax(l1,l2),l3)>eflag_T_lmax) continue;
	//
	// 	L = l1+l2+l3;
	// 	L1 = L-2*l1;
	// 	L2 = L-2*l2;
	// 	L3 = L-2*l3;
	// 	p1 = sqrt((2e0*l1+1e0)*(2e0*l2+1e0)*(2e0*l2+1e0)/4e0*3.1415927e0);
	// 	p2 = pow(-1,(L/2)%2)*twopi*pow(L+1e0,-5e-1);
	// 	p3 = pow(L+thrd,5e-1)*pow(L+sxth,-2.5e-1);
	// 	p4 = pow(L1+sxth,2.5e-1)*pow(L1+thrd,-5e-1);
	// 	p5 = pow(L2+sxth,2.5e-1)*pow(L2+thrd,-5e-1);
	// 	p6 = pow(L3+sxth,2.5e-1)*pow(L3+thrd,-5e-1);
	// 	factor = p1*p2*p3*p4*p5*p6;
	//
	// 	// if(l1+l2+l3==1002)printf("%d\t%d\t%d\t%e\n",l1,l2,l3,bisp/factor);
	// 	bispectrum[l1][l2][l3] = bisp/factor;
	//
	// }
	
	// model = 6;
	// read_eigen(0,0,0);
	
}

void calculate_ISW_3D(int min, int max, tetrapyd_limits* block, double* alpha){

	int n,r,t,i,j,k;
	double x,y;
	int size = get_terms_late();
	double alpha_private[size];
	//printf("Alpha allocated\n");
	for(r=0;r<size;r++){
		alpha[r] = 0e0;
	}
	//printf("Alpha initialized\n");
	#pragma omp parallel default(none) private(n,r,t,i,j,k,x,y,alpha_private) shared(alpha, block, min, max, size, weight)
	{
		t = omp_get_thread_num();
		i = block[t].i_bgn;
		j = block[t].j_bgn;
		k = block[t].k_bgn;

		for(r=0;r<size;r++){
			alpha_private[r] = 0e0;
		}
		//printf("Alpha private initialized\n");
		for(n=0;n<block[t].loops;n++){
			//printf("Try to compute x");
			//printf("get wieght = %e\n", get_weight(i,j,k));
			//printf("weight = %e\n", weight[i]);
			//printf("ISW Bispectrum = %e\n", calculate_ISW_bispectrum(i,j,k));
			//x = get_weight(i,j,k)*plijk_TTT(5,i,j,k);
			x = get_weight(i,j,k)*calculate_ISW_bispectrum(i,j,k)*weight[i]*weight[j]*weight[k];
			//printf("x computed\n");
			for(r=0;r<size;r++){
				y = plijk_TTT(r,i,j,k);
				alpha_private[r] += x*y;
				//printf("y = %e and alpha_private computed\n", y);
			}
			get_ijk_next_XXX(max,&i,&j,&k);
			//printf("get ijk next XXX\n");
		}

		#pragma omp critical
		{
			for(r=0;r<size;r++){
				alpha[r] += alpha_private[r];
				//printf("Alpha update\n");
			}
		}


	}
	
	return;
}

void calculate_PS_3D(int min, int max, tetrapyd_limits* block, double* alpha){

	int n,r,t,i,j,k;
	double x,y;
	int size = get_terms_late();
	int test = 0;
	double alpha_private[size];
	
	for(r=0;r<size;r++){
		alpha[r] = 0e0;
	}

	#pragma omp parallel default(none) private(n,r,t,i,j,k,x,y,test,alpha_private) shared(alpha, block, min, max, size, weight)
	{
		t = omp_get_thread_num();
		i = block[t].i_bgn;
		j = block[t].j_bgn;
		k = block[t].k_bgn;

		for(r=0;r<size;r++){
			alpha_private[r] = 0e0;
		}

		for(n=0;n<block[t].loops;n++){
			x = get_weight(i,j,k)*(weight[i]*weight[j]*weight[k]/sqrt(calculate_geometric(i,j,k)));
			
			for(r=0;r<size;r++){
				y = plijk_TTT(r,i,j,k);
				alpha_private[r] += x*y;
			}
			get_ijk_next_XXX(max,&i,&j,&k);
			test = 1;
		}

		#pragma omp critical
		{
			for(r=0;r<size;r++){
				alpha[r] += alpha_private[r];
			}
		}


	}
	return;
}
	

void calculate_mode_3D(int min, int max, tetrapyd_limits* block, double* alpha){

	int n,r,t,i,j,k;
	double x,y,z;
	int size = get_terms_late();
	double alpha_private[size];
	
	for(r=0;r<size;r++){
		alpha[r] = 0e0;
	}

	#pragma omp parallel default(none) private(n,r,t,i,j,k,x,y,alpha_private) shared(alpha, block, min, max, size, weight, bispectrum)
	{
		t = omp_get_thread_num();
		i = block[t].i_bgn;
		j = block[t].j_bgn;
		k = block[t].k_bgn;

		for(r=0;r<size;r++){
			alpha_private[r] = 0e0;
		}

		for(n=0;n<block[t].loops;n++){
			x = get_weight(i,j,k)*calculate_3D_bispectrum(i,j,k) * weight[i] * weight[j] * weight[k];
			// x = bispectrum[i][j][k] * weight[i] * weight[j] * weight[k];
			// x = (i+j+k) * weight[i] * weight[j] * weight[k];
			// if((i==2 || i%10==0) && (j==2 || j%10==0) && (k==2 || k%10==0)) printf("%d\t%d\t%d\t%e\n",i,j,k,x);
			// x = get_weight(i,j,k);
			for(r=0;r<size;r++){
				y = plijk_TTT(r,i,j,k);
				alpha_private[r] += x*y;
			}
			
			get_ijk_next_XXX(max,&i,&j,&k);
		}

		#pragma omp critical
		{
			for(r=0;r<size;r++){
				alpha[r] += alpha_private[r];
			}
		}


	}
	
	return;
}

double calculate_3D_bispectrum(int l1, int l2, int l3){
	
	int i;
	double area = 0e0;
	double result = 0e0;
	double x = 0e0;
	int terms = get_terms_late();
	for(i=0;i<terms;i++){
		x = calculate_xint_TTT(l1, l2, l3, i, xsize, xvec);
		result += modes[i] * x;
	}
	
	if((l1==2 || l1%10==0) && (l2==2 || l2%10==0) && (l3==2 || l3%10==0)) printf("%d\t%d\t%d\t%e\n",l1,l2,l3,result);
	
	return result;
}

double get_weight(int i, int j, int k){
	
	double result = 0e0;
	
	result = permsix(i,j,k)*calculate_geometric(i,j,k)*pow(2.0*i+1.0,-1.0/3.0)*pow(2.0*j+1.0,-1.0/3.0)*pow(2.0*k+1.0,-1.0/3.0);
		
	return result;
}
			
			
			
			
			
			
			
			
			
			
			
			
			