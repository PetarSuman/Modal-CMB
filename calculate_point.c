#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_spline.h>
#include "global.h"

double calculate_integral(double k1, double k2, double k3, double power, int size, double *x, gsl_spline* sp1, gsl_interp_accel* acc1, gsl_spline* sp2, gsl_interp_accel* acc2, gsl_spline* sp3, gsl_interp_accel* acc3);

double calculate_point(double **array, double x, double y, int b_size, double *b_x, gsl_spline* b_sp1, gsl_interp_accel* b_acc1, gsl_spline* b_sp2, gsl_interp_accel* b_acc2, gsl_spline* b_sp3, gsl_interp_accel* b_acc3, int t_size, double *t_x, gsl_spline* t_sp1, gsl_interp_accel* t_acc1, gsl_spline* t_sp2, gsl_interp_accel* t_acc2, gsl_spline* t_sp3, gsl_interp_accel* t_acc3){

// 	MPI sync
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if ( myrank == 0 ) sync_tasks(3,4);
	
	double result;
	double tolerance = 1/ (initial_grid_size * pow(2,depth + 1));
	double b_power = 2;
	double t_power = -1;
	
	double part1;
	double part2;
	double part3;
	
	double k1 = (1-x);
	double k2 = (1+x+y)/2;
	double k3 = (1+x-y)/2;
// **1**
	if (k1 > 1) {k1 = 1;}
	if (k2 > 1) {k2 = 1;}
	if (k3 > 1) {k3 = 1;}
// **2**
	
	shape_params params;
	params.a1 = 0e0;
	params.a2 = 0e0;
	params.a3 = 0e0;

	if (k1 < tolerance || k2 < tolerance ||  k3 < tolerance ) {
		result = 0;
	} else {
		part1 = slice(k1, k2, k3, params);
		if (part1 == 0) {
			part2 = 0;
			part3 = 0;
		} else {
			part2 = calculate_integral(k1, k2, k3, t_power, t_size, t_x, t_sp1, t_acc1, t_sp2, t_acc2, t_sp3, t_acc3);
			if (part2 == 0) {
				part3 = 0;
			} else {
				part3 = calculate_integral(k1, k2, k3, b_power, b_size, b_x, b_sp1, b_acc1, b_sp2, b_acc2, b_sp3, b_acc3);
			}
		}
		result = part1 * part2 * part3;
	}

// **4**
	if ( 1 - k1 < tolerance || 1 - k2 < tolerance ||  1 - k3 < tolerance ) {
		result = 2 * result;
	}
	
	return  result;
}

double calculate_integral(double k1, double k2, double k3, double power, int size, double *x, gsl_spline* sp1, gsl_interp_accel* acc1, gsl_spline* sp2, gsl_interp_accel* acc2, gsl_spline* sp3, gsl_interp_accel* acc3){

	int n;
	double z;
	double z1;
	double z2;
	double z3;
	
	double factor;
	double result = 0;
	double part1;
	double part2;
	double part3;
	
	double min = x[0];
	double max = x[size-1];
	
	double y[size];
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, size);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	
	y[0] = 0;
	
	for (n = 1; n < size; n++) {
		
		z = x[n];
		
		z1 = z*k1;
		z2 = z*k2;
		z3 = z*k3;
		
		if (power == 2) {
			factor = z * z;
		} else if (power == -1) {
			factor = scale(2.0*z/3.0) / z;
		} else {
			factor = pow(z,power);
		}
		
		part1 = gsl_spline_eval(sp1,z1,acc1);
		part2 = gsl_spline_eval(sp2,z2,acc2);
		part3 = gsl_spline_eval(sp3,z3,acc3);
		
		y[n] = factor * part1 * part2 * part3;
		
	}
	
	gsl_spline_init(sp,x,y,size);
	result = gsl_spline_eval_integ(sp,min,max,acc);
	
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
	return result;
}
