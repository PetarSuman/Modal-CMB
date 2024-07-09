#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "global.h"

struct my_int_params{double nu1; double nu2;};
double integrand_ll(double *k, size_t dim, void *params);
double integrand_lx(double *k, size_t dim, void *params);
double integrand_xx(double *k, size_t dim, void *params);

gsl_rng *rng;

double calculate_correlation(double kmin, double kmax, double nu1, double nu2){

	double res1, err1;
	double res2, err2;
	double res3, err3;
	double xl[3] = {kmin,kmin,kmin};
	double xu[3] = {kmax,kmax,kmax};
	
	const gsl_rng_type *T;
	gsl_rng *r;

	struct my_int_params params = {nu1,nu2};
// 	gsl_monte_function LL = { &integrand_ll, 3, &params};
// 	gsl_monte_function LX = { &integrand_lx, 3, &params};
	gsl_monte_function XX = { &integrand_xx, 3, &params};
	
	size_t calls = 100000;
	
	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	gsl_monte_plain_state *s = gsl_monte_plain_alloc(3);
// 	gsl_monte_plain_integrate (&LL, xl, xu, 3, calls, r, s, &res1, &err1);
// 	gsl_monte_plain_init(s);
// 	gsl_monte_plain_integrate (&LX, xl, xu, 3, calls, r, s, &res2, &err2);
// 	gsl_monte_plain_init(s);
	gsl_monte_plain_integrate (&XX, xl, xu, 3, calls, r, s, &res1, &err1);
	gsl_monte_plain_free(s);
	gsl_rng_free (r);
	
// 	printf("Nu:\t%e\t%e\t%e\t%e\n",nu, res1, res2, res3);
	
// 	return res2/sqrt(res1*res3);
	return res1;
}

double calculate_localnorm(double kmin, double kmax, double nu1, double nu2){

	double res1, err1;
	double res2, err2;
	double res3, err3;
	double xl[3] = {kmin,kmin,kmin};
	double xu[3] = {kmax,kmax,kmax};
	
	const gsl_rng_type *T;
	gsl_rng *r;

	struct my_int_params params = {nu1,nu2};
	gsl_monte_function LL = { &integrand_ll, 3, &params};
// 	gsl_monte_function LX = { &integrand_lx, 3, &params};
// 	gsl_monte_function XX = { &integrand_xx, 3, &params};
	
	size_t calls = 100000;
	
	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	gsl_monte_plain_state *s = gsl_monte_plain_alloc(3);
	gsl_monte_plain_integrate (&LL, xl, xu, 3, calls, r, s, &res1, &err1);
	gsl_monte_plain_init(s);
// 	gsl_monte_plain_integrate (&LX, xl, xu, 3, calls, r, s, &res2, &err2);
// 	gsl_monte_plain_init(s);
// 	gsl_monte_plain_integrate (&XX, xl, xu, 3, calls, r, s, &res1, &err1);
// 	gsl_monte_plain_free(s);
	gsl_rng_free(r);
	
// 	printf("Nu:\t%e\t%e\t%e\t%e\n",nu, res1, res2, res3);
	
// 	return res2/sqrt(res1*res3);
	return res1;
}

double integrand_ll(double *k, size_t dim, void *params){
	struct my_int_params *p = (struct my_int_params *)params;

	double k1 = k[0];
	double k2 = k[1];
	double k3 = k[2];
	double nu1 = p->nu1;
	double nu2 = p->nu2;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p444 = pow(p111,4);
	double result;
	shape_params params2;
	params2.a1 = 0e0;
	params2.a2 = 0e0;
	params2.a3 = 0e0;
	
	if(k1>k2+k3||k2>k3+k1||k3>k1+k2||p111==0.0){
		result = 0.0;
	}else{
		result = p444*local(k1,k2,k3,params2)*local(k1,k2,k3,params2)/p1;
	}
	return result;
}

double integrand_lx(double *k, size_t dim, void *params){
	struct my_int_params *p = (struct my_int_params *)params;

	double k1 = k[0];
	double k2 = k[1];
	double k3 = k[2];
	double nu1 = p->nu1;
	double nu2 = p->nu2;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p444 = pow(p111,4);
	double result;
	
	if(k1>k2+k3||k2>k3+k1||k3>k1+k2||p111==0.0){
		result = 0.0;
	}else{
// 		result = p444*local(k1,k2,k3)*xingang(nu1,k1,k2,k3)/p1;
	}
	return result;
}

double integrand_xx(double *k, size_t dim, void *params){
	struct my_int_params *p = (struct my_int_params *)params;

	double k1 = k[0];
	double k2 = k[1];
	double k3 = k[2];
	double nu1 = p->nu1;
	double nu2 = p->nu2;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p444 = pow(p111,4);
	double result;
	
	if(k1>k2+k3||k2>k3+k1||k3>k1+k2||p111==0.0){
		result = 0.0;
	}else{
// 		result = p444*xingang(nu1,k1,k2,k3)*xingang(nu2,k1,k2,k3)/p1;
	}
	return result;
}

void init_rng(int seed){
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,seed);
	return;
}

double get_rng_gauss(){
	double x1,x2,z1;
	x1 = 1e0-gsl_rng_uniform(rng);
	x2 = 1e0-gsl_rng_uniform(rng);
	z1 = sqrt(-2e0*log(x1))*cos(2e0*M_PI*x2);
	return z1;
}