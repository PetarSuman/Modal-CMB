#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_spline.h>
#include "global.h"

double calculate_kint_T(int alpha, int l, double p){
	double result;
	int i,j;
	
	double xmax = get_xmax();
	double kmax = get_kmax();
	
	int xsize = get_xsize();
	//double x[xsize];
	//double b[xsize];
	double *x = (double *)create_vector(xsize);
	double *b = (double *)create_vector(xsize);
	get_xvec(x);
	get_bvec(l, b);
	
	
	int ksize = get_ksize();
	//double k[ksize];
	//double t[ksize];
	//double y[ksize];
	double *k = (double *)create_vector(ksize);
	double *t = (double *)create_vector(ksize);
	double *y = (double *)create_vector(ksize);
	get_kvec(k);
	get_tvec_T(l, t);
	
	double kmin = k[0];
	double pt;
	
	gsl_spline* spx =  gsl_spline_alloc (gsl_interp_cspline, xsize);
	gsl_interp_accel* accx = gsl_interp_accel_alloc();
	
	gsl_spline_init(spx,x,b,xsize);
	y[0] = 0;
	
	for(i=1;i<ksize;i++){
		pt = p*k[i];
		if(pt>xmax){
			y[i] = 0.0;
		}else{
			y[i] = t[i] * get_basis_prim(i,alpha) * gsl_spline_eval(spx,pt,accx);
		}
 		//if(alpha==0&&l==20&&p==200.0)printf("%d\t%e\t%d\t%e\t%e\n",alpha,p,i,k[i],y[i]);
	}
	
	if(alpha==0&&l==20&&p==200.0)printf("%e\t%e\n",kmin,kmax);
	
	gsl_spline_free(spx);
	gsl_interp_accel_free(accx);
	
	gsl_spline* spk =  gsl_spline_alloc (gsl_interp_cspline, ksize);
	gsl_interp_accel* acck = gsl_interp_accel_alloc();
	gsl_spline_init(spk,k,y,ksize);
	result = gsl_spline_eval_integ(spk,kmin,kmax,acck);
	gsl_spline_free(spk);
	gsl_interp_accel_free(acck);
		
	if(alpha==0&&l==20&&p==200.0)printf("result\t%e\n",result);
		
	if(p<1000&&l>500&&(result>1e-5||result<-1e-5))result = 0.0;
		
	destroy_vector(x);
	destroy_vector(b);
	destroy_vector(k);
	destroy_vector(t);
	destroy_vector(y);

	return result;
	
}

double calculate_kint_E(int alpha, int l, double p){
	double result;
	int i,j;
	
	double xmax = get_xmax();
	double kmax = get_kmax();
	//printf("Kmax\n");
	int xsize = get_xsize();
	//double x[xsize];
	//double b[xsize];

	double *x = (double *)create_vector(xsize);
	double *b = (double *)create_vector(xsize);	

	get_xvec(x);
	get_bvec(l, b);
	//printf("bvec\n");
	
	int ksize = get_ksize();
	//printf("ksize\n");
	double *k = (double *)create_vector(ksize);
	//double k[ksize];
	//printf("double k\n");
	double *t = (double *)create_vector(ksize);
	//double t[ksize];
	//printf("double t\n");
	//double y[ksize];
	double *y = (double *)create_vector(ksize);
	//printf("double y\n");
	get_kvec(k);
	//printf("Try to get tvec E\n");
	get_tvec_E(l, t);
	
	double kmin = k[0];
	double pt;
	//printf("Evec passed\n");
	gsl_spline* spx =  gsl_spline_alloc (gsl_interp_cspline, xsize);
	gsl_interp_accel* accx = gsl_interp_accel_alloc();
	
	gsl_spline_init(spx,x,b,xsize);
	y[0] = 0;
	//printf("For loop:\n");
	for(i=1;i<ksize;i++){
		pt = p*k[i];
		if(pt>xmax){
			y[i] = 0.0;
		}else{
			y[i] = t[i] * get_basis_prim(i,alpha) * gsl_spline_eval(spx,pt,accx);
		}
 		//if(l==1990&&p==200.0)printf("%d\t%e\t%d\t%e\t%e\n",alpha,p,i,x[i],y[i]);
	}
	//printf("Free...\n");
	gsl_spline_free(spx);
	gsl_interp_accel_free(accx);
	
	gsl_spline* spk =  gsl_spline_alloc (gsl_interp_cspline, ksize);
	gsl_interp_accel* acck = gsl_interp_accel_alloc();
	gsl_spline_init(spk,k,y,ksize);
	result = gsl_spline_eval_integ(spk,kmin,kmax,acck);
	gsl_spline_free(spk);
	gsl_interp_accel_free(acck);
		
		
	if(p<1000&&l>500&&(result>1e-5||result<-1e-5))result = 0.0;
		
	destroy_vector(x);
	destroy_vector(b);
	destroy_vector(k);
	destroy_vector(t);
	destroy_vector(y);	
	return result;
	
}

double calculate_xint_TTT(int l1, int l2, int l3, int n, int xsize, double *xvec) {

	double result;
	int i,j,k;
	int a1,a2,a3;
	find_perm_prim(n,&a1,&a2,&a3);
	
	int size = xsize;
	double x[size];
	double y[size];
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	for(i=0;i<size;i++){
		x[i] = xvec[i];
	}

    // double* xdiff = (double*) _mm_malloc(xsize * sizeof(double), 64);
    // double* ixdiff = (double*) _mm_malloc(xsize * sizeof(double), 64);
    // for (i = 0; i < xsize-1; i++)
    // {
    //     xdiff[i] = xvec[i+1] - xvec[i];
    //     ixdiff[i] = 1.0 / xdiff[i];
    // }
	
	double min = x[0];
	double max = x[size-1];
	
	double p1,p2,p3,p4,p5,p6,p7,p8,p9;
	double t1,t2,t3,t4,t5,t6;
	
	if(l1<eflag_T_lmin || l1>eflag_T_lmax || l2<eflag_T_lmin || l2>eflag_T_lmax ||l3<eflag_T_lmin || l3>eflag_T_lmax){
		printf("[%d] error l's out of range: %d\t%d\t%d\t%d\n",rank,n,l1,l2,l3);
		result = 0e0;
	}else{
	
		for (i=0;i<size;i++){
	
			p1 = get_qtilde_T(a1,l1,i);
			p2 = get_qtilde_T(a2,l1,i);
			p3 = get_qtilde_T(a3,l1,i);
			p4 = get_qtilde_T(a1,l2,i);
			p5 = get_qtilde_T(a2,l2,i);
			p6 = get_qtilde_T(a3,l2,i);
			p7 = get_qtilde_T(a1,l3,i);
			p8 = get_qtilde_T(a2,l3,i);
			p9 = get_qtilde_T(a3,l3,i);
		
			t1 = p1*p5*p9;
			t2 = p2*p6*p7;
			t3 = p3*p4*p8;
			t4 = p3*p5*p7;
			t5 = p2*p4*p9;
			t6 = p1*p6*p8;
		
			y[i] = x[i] * x[i] * (t1+t2+t3+t4+t5+t6)/6.0;
			if(isnan(y[i])||isinf(y[i]))printf("%d\t%d\t%d\t%d\t%d\t%e\t%e\n",i,l1,l2,l3,n,x[i],y[i]);
			if( i<10 && fabs(y[i])>1e-18 )y[i]=0;
		}
	
		gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, size);
		gsl_interp_accel* acc = gsl_interp_accel_alloc();

		gsl_spline_init(sp,x,y,size);
		result = gsl_spline_eval_integ(sp,min,max,acc);

		gsl_spline_free(sp);
		gsl_interp_accel_free(acc);
	
	}
	
// 	return result*pow(kmax,-3);
    // result = integrate(size, xdiff, ixdiff, y);
    // _mm_free(xdiff);
    // _mm_free(ixdiff);
	
	return result;
}

double calculate_xint_TTE(int l1, int l2, int l3, int n, int xsize, double *xvec) {

	double result;
	int i,j,k;
	int a1,a2,a3;
	find_perm_prim(n,&a1,&a2,&a3);
	
	int size = xsize;
	double x[size];
	double y[size];
	
	for(i=0;i<size;i++){
		x[i] = xvec[i];
	}
	
	double min = x[0];
	double max = x[size-1];
	
	double p1,p2,p3,p4,p5,p6,p7,p8,p9;
	double t1,t2,t3,t4,t5,t6;
	
	for (i=0;i<size;i++){
	
		p1 = get_qtilde_T(a1,l1,i);
		p2 = get_qtilde_T(a2,l1,i);
		p3 = get_qtilde_T(a3,l1,i);
		p4 = get_qtilde_T(a1,l2,i);
		p5 = get_qtilde_T(a2,l2,i);
		p6 = get_qtilde_T(a3,l2,i);
		p7 = get_qtilde_E(a1,l3,i);
		p8 = get_qtilde_E(a2,l3,i);
		p9 = get_qtilde_E(a3,l3,i);
		
		t1 = p1*p5*p9;
		t2 = p2*p6*p7;
		t3 = p3*p4*p8;
		t4 = p3*p5*p7;
		t5 = p2*p4*p9;
		t6 = p1*p6*p8;
		
		y[i] = x[i] * x[i] * (t1+t2+t3+t4+t5+t6)/6.0;
		if( i<10 && fabs(y[i])>1e-18 )y[i]=0;
	}
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, size);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	
	gsl_spline_init(sp,x,y,size);
	result = gsl_spline_eval_integ(sp,min,max,acc);
	
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
// 	return result*pow(kmax,-3);
	return result;
}

double calculate_xint_TEE(int l1, int l2, int l3, int n, int xsize, double *xvec) {

	double result;
	int i,j,k;
	int a1,a2,a3;
	find_perm_prim(n,&a1,&a2,&a3);
	
	int size = xsize;
	double x[size];
	double y[size];
	
	for(i=0;i<size;i++){
		x[i] = xvec[i];
	}
	
	double min = x[0];
	double max = x[size-1];
	
	double p1,p2,p3,p4,p5,p6,p7,p8,p9;
	double t1,t2,t3,t4,t5,t6;
	
	for (i=0;i<size;i++){
	
		p1 = get_qtilde_E(a1,l1,i);
		p2 = get_qtilde_E(a2,l1,i);
		p3 = get_qtilde_E(a3,l1,i);
		p4 = get_qtilde_E(a1,l2,i);
		p5 = get_qtilde_E(a2,l2,i);
		p6 = get_qtilde_E(a3,l2,i);
		p7 = get_qtilde_T(a1,l3,i);
		p8 = get_qtilde_T(a2,l3,i);
		p9 = get_qtilde_T(a3,l3,i);
		
		t1 = p1*p5*p9;
		t2 = p2*p6*p7;
		t3 = p3*p4*p8;
		t4 = p3*p5*p7;
		t5 = p2*p4*p9;
		t6 = p1*p6*p8;
		
		y[i] = x[i] * x[i] * (t1+t2+t3+t4+t5+t6)/6.0;
		if( i<10 && fabs(y[i])>1e-18 )y[i]=0;
	}
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, size);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	
	gsl_spline_init(sp,x,y,size);
	result = gsl_spline_eval_integ(sp,min,max,acc);
	
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
// 	return result*pow(kmax,-3);
	return result;
}

double calculate_xint_EEE(int l1, int l2, int l3, int n, int xsize, double *xvec) {

	double result;
	int i,j,k;
	int a1,a2,a3;
	find_perm_prim(n,&a1,&a2,&a3);
	
	int size = xsize;
	double x[size];
	double y[size];
	
	for(i=0;i<size;i++){
		x[i] = xvec[i];
	}
	
	double min = x[0];
	double max = x[size-1];
	
	double p1,p2,p3,p4,p5,p6,p7,p8,p9;
	double t1,t2,t3,t4,t5,t6;
	
	for (i=0;i<size;i++){
	
		p1 = get_qtilde_E(a1,l1,i);
		p2 = get_qtilde_E(a2,l1,i);
		p3 = get_qtilde_E(a3,l1,i);
		p4 = get_qtilde_E(a1,l2,i);
		p5 = get_qtilde_E(a2,l2,i);
		p6 = get_qtilde_E(a3,l2,i);
		p7 = get_qtilde_E(a1,l3,i);
		p8 = get_qtilde_E(a2,l3,i);
		p9 = get_qtilde_E(a3,l3,i);
		
		t1 = p1*p5*p9;
		t2 = p2*p6*p7;
		t3 = p3*p4*p8;
		t4 = p3*p5*p7;
		t5 = p2*p4*p9;
		t6 = p1*p6*p8;
		
		y[i] = x[i] * x[i] * (t1+t2+t3+t4+t5+t6)/6.0;
		if( i<10 && fabs(y[i])>1e-18 )y[i]=0;
	}
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, size);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	
	gsl_spline_init(sp,x,y,size);
	result = gsl_spline_eval_integ(sp,min,max,acc);
	
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
// 	return result*pow(kmax,-3);
	return result;
}

double calculate_kint_tri(int l, int alpha, double p){

	double result;
	int i,j;
	
	double xmax = get_xmax();
	double kmax = get_kmax();
	
	int xsize = get_xsize();
	double x[xsize];
	double b[xsize];
	get_xvec(x);
	get_bvec(l, b);
	
	
	int ksize = get_ksize();
	double k[ksize];
	double t[ksize];
	double y[ksize];
	get_kvec(k);
	get_tvec_T(l, t);
	
	double kmin = k[0];
	double pt,pt2;	
	
	gsl_spline* spx =  gsl_spline_alloc (gsl_interp_cspline, xsize);
	gsl_interp_accel* accx = gsl_interp_accel_alloc();
	
	gsl_spline_init(spx,x,b,xsize);
	y[0] = 0;
	
	for(i=1;i<ksize;i++){
		pt = p*k[i];
		pt2 = pow(k[i],-0.25e0);
		if(pt>xmax){
			y[i] = 0.0;
		}else{
			y[i] = t[i] * get_basis_tri_prim(i,alpha) * gsl_spline_eval(spx,pt,accx)*pt2;
		}
	}
	
	gsl_spline_free(spx);
	gsl_interp_accel_free(accx);
	
	gsl_spline* spk =  gsl_spline_alloc (gsl_interp_cspline, ksize);
	gsl_interp_accel* acck = gsl_interp_accel_alloc();
	gsl_spline_init(spk,k,y,ksize);
	result = gsl_spline_eval_integ(spk,kmin,kmax,acck);
	gsl_spline_free(spk);
	gsl_interp_accel_free(acck);
		
		
	if(p<1000&&l>500&&(result>1e-5||result<-1e-5))result = 0.0;
	
	return result;
	
}

double calculate_xint_tri(int l1, int l2, int l3, int l4, int n, int xsize, double *xvec) {

	double result;
	int i,j,k;
	int a1,a2,a3,a4;
	find_perm_tri_prim(n,&a1,&a2,&a3,&a4);
	
	int size = xsize;
	double x[size];
	double y[size];
	
	for(i=0;i<size;i++){
		x[i] = xvec[i];
	}
	
	double min = x[0];
	double max = x[size-1];
	
	double p[4][4];
	double t[24];
	
	for (i=0;i<size;i++){
		p[0][0] = get_beta_tri(a1,l1,i);
		p[0][1] = get_beta_tri(a1,l2,i);
		p[0][2] = get_beta_tri(a1,l3,i);
		p[0][3] = get_beta_tri(a1,l4,i);
		p[1][0] = get_beta_tri(a2,l1,i);
		p[1][1] = get_beta_tri(a2,l2,i);
		p[1][2] = get_beta_tri(a2,l3,i);
		p[1][3] = get_beta_tri(a2,l4,i);
		p[2][0] = get_beta_tri(a3,l1,i);
		p[2][1] = get_beta_tri(a3,l2,i);
		p[2][2] = get_beta_tri(a3,l3,i);
		p[2][3] = get_beta_tri(a3,l4,i);
		p[3][0] = get_beta_tri(a4,l1,i);
		p[3][1] = get_beta_tri(a4,l2,i);
		p[3][2] = get_beta_tri(a4,l3,i);
		p[3][3] = get_beta_tri(a4,l4,i);

		t[0] = p[0][0]*p[1][1]*p[2][2]*p[3][3];
		t[1] = p[0][0]*p[1][1]*p[2][3]*p[3][2];
		t[2] = p[0][0]*p[1][2]*p[2][1]*p[3][3];
		t[3] = p[0][0]*p[1][2]*p[2][3]*p[3][1];
		t[4] = p[0][0]*p[1][3]*p[2][1]*p[3][2];
		t[5] = p[0][0]*p[1][3]*p[2][2]*p[3][1];
		t[6] = p[0][1]*p[1][0]*p[2][2]*p[3][3];
		t[7] = p[0][1]*p[1][0]*p[2][3]*p[3][2];
		t[8] = p[0][1]*p[1][2]*p[2][0]*p[3][3];
		t[9] = p[0][1]*p[1][2]*p[2][3]*p[3][0];
		t[10] = p[0][1]*p[1][3]*p[2][0]*p[3][2];
		t[11] = p[0][1]*p[1][3]*p[2][2]*p[3][0];
		t[12] = p[0][2]*p[1][1]*p[2][0]*p[3][3];
		t[13] = p[0][2]*p[1][1]*p[2][3]*p[3][0];
		t[14] = p[0][2]*p[1][0]*p[2][1]*p[3][3];
		t[15] = p[0][2]*p[1][0]*p[2][3]*p[3][1];
		t[16] = p[0][2]*p[1][3]*p[2][1]*p[3][0];
		t[17] = p[0][2]*p[1][3]*p[2][0]*p[3][1];
		t[18] = p[0][3]*p[1][1]*p[2][2]*p[3][0];
		t[19] = p[0][3]*p[1][1]*p[2][0]*p[3][2];
		t[20] = p[0][3]*p[1][2]*p[2][1]*p[3][0];
		t[21] = p[0][3]*p[1][2]*p[2][0]*p[3][1];
		t[22] = p[0][3]*p[1][0]*p[2][1]*p[3][2];
		t[23] = p[0][3]*p[1][0]*p[2][2]*p[3][1];
		
		y[i]=0.0;
		for(j=0;j<24;j++){
			y[i] += t[j];
		}
		
		if( i<10 && fabs(y[i])>1e-18 )y[i]=0;
		
		y[i] /= 24.0;
	}

	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, size);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	
	gsl_spline_init(sp,x,y,size);
	result = gsl_spline_eval_integ(sp,min,max,acc);
	
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
	return result;
}
