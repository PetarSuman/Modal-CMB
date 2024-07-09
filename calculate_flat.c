#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include "global.h"
#include <gsl/gsl_spline.h>

double calculate_flat(int l1, int l2, int l3){
	
	int f_size = get_fsize();
	
	double total_area = 0;
	
	int i,j,k,n;
	
// 	double tau = 1.42e4;
	double tau = get_ftau();
// 	printf("tauR:\t%e\n",tau);
	
	double x[f_size];
	double f1[f_size];
	double f2[f_size];
	double f3[f_size];
	
	get_fkvec(x);
	get_fvec(l1,f1);
	get_fvec(l2,f2);
	get_fvec(l3,f3);
// 	printf("%e\n",f1[10]);
	double factor = tau*tau*M_1_PI*M_1_PI/2.0;
	double area1 = 0;
	double area2 = 0;
	
	double k1, k2, k3, k4, bi;
	
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	gsl_interp_accel* acc1 = gsl_interp_accel_alloc();
	gsl_interp_accel* acc2 = gsl_interp_accel_alloc();
	gsl_interp_accel* acc3 = gsl_interp_accel_alloc();
	gsl_interp_accel* acc4 = gsl_interp_accel_alloc();
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, f_size);
	gsl_spline* spint1 =  gsl_spline_alloc (gsl_interp_cspline, f_size);
	gsl_spline* spint2 =  gsl_spline_alloc (gsl_interp_cspline, f_size);
	gsl_spline* spint3 =  gsl_spline_alloc (gsl_interp_cspline, f_size);
	gsl_spline* spint4 =  gsl_spline_alloc (gsl_interp_cspline, f_size);
	
	double x1 = 0;
	double x2 = 0;
	double x3 = 0;
	double x4 = 0;
	double y1 = 0;
	double y2 = 0;
	double y3 = 0;
	double y4 = 0;
	double step = 0;
	double min = x[0];
	double max = get_fmax();
	int ifail;
	
	double integrand1[f_size];
	double integrand2[f_size];
	double integrand3[f_size];
	double integrand4[f_size];
	integrand1[0] = 0;
	integrand2[0] = 0;
	integrand3[0] = 0;
	integrand4[0] = 0;
	
	shape_params params;
	params.a1 = 0e0;
	params.a2 = 0e0;
	params.a3 = 0e0;
	
//	****** Spline ******
	
	gsl_spline_init (sp,x,f3,f_size);
	
	for ( i = 1; i <  f_size; i++ ) {
		integrand2[i] = 0;
		integrand4[i] = 0;
		
		x1 = x[i];
		k1 = sqrt(x1*x1 + (l1*l1 / (tau*tau)));
		y1 = f1[i];
		
		if (y1 != 0) {
			for ( j = 1; j <  f_size; j++ ) {
				integrand1[j] = 0;
				integrand3[j] = 0;
				
				x2 = x[j];
				k2 = sqrt(x2*x2 + (l2*l2 / (tau*tau)));
				y2 = f2[j];
				
				x3 = x1 + x2;
				k3 = sqrt(x3*x3 + (l3*l3 / (tau*tau)));
				
				x4 = fabs(x1 - x2);
				k4 = sqrt(x4*x4 + (l3*l3 / (tau*tau)));
				
				if (y2 != 0 && x3 < max) {
					y3 = gsl_spline_eval (sp,x3,acc);
					bi =  shape(k1, k2, k3, params);
					integrand1[j] = y1 * y2 * y3 * bi;
				}
				if (y2 != 0 && x4 < max) {
					y4 = gsl_spline_eval (sp,x4,acc);
					bi =  shape(k1, k2, k4, params);
					integrand3[j] = y1 * y2 * y4 * bi;
				}
			}
			
			gsl_spline_init(spint1,x,integrand1,f_size);
			integrand2[i] = gsl_spline_eval_integ(spint1,min,max,acc1);
			
			gsl_spline_init(spint3,x,integrand3,f_size);
			integrand4[i] = gsl_spline_eval_integ(spint3,min,max,acc3);
			
		}
	}

	gsl_spline_init(spint2,x,integrand2,f_size);
	area1 = gsl_spline_eval_integ(spint2,min,max,acc2);

	gsl_spline_init(spint4,x,integrand4,f_size);
	area2 = gsl_spline_eval_integ(spint4,min,max,acc4);
	
	total_area = factor*(area1 + area2);
	
	return total_area;
}
