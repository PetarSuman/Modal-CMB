#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include "global.h"

#define CMD "python ./Python/SRFT_call.py"
#define MAX_LINE_LEN 500

double function1(double x);
double function2(double x);
double function3(double ksum, double kmin2);
double adshead_power(double x, shape_params params);
double adshead_i0(double x, shape_params params);
double adshead_i1(double x, shape_params params);
double adshead_i2(double x, shape_params params);
double adshead_D(double x, shape_params params);
static double **NicolaShape44;

int cs_feature_prior_check(shape_params params);
void cs_feature_init(shape_params params, int *n_total, double **k, double **dpp);

double shape(double k1, double k2, double k3, shape_params params) {

	double result;

	if(nscalar!=1.0 && (model==2||model==6||model==14||model==16)){
		k1 = kpivot*pow(k1/kpivot,(4.0-nscalar)/3.0);
		k2 = kpivot*pow(k2/kpivot,(4.0-nscalar)/3.0);
		k3 = kpivot*pow(k3/kpivot,(4.0-nscalar)/3.0);
	}
	switch (model) {
		case 1: result = dbi(k1, k2, k3, params); break;
		case 2: result = equilateral(k1, k2, k3, params); break;
		case 3: result = smooth2(k1, k2, k3, params); break;
		case 4: result = feature(k1, k2, k3, params); break;
		case 5: result = ghost(k1, k2, k3, params); break;
		case 6: result = local(k1, k2, k3, params); break;
		case 7: result = single(k1, k2, k3, params); break;
		case 8: result = warm(k1, k2, k3, params); break;
		case 9: result = smooth1(k1, k2, k3, params); break;
		case 10: result = nonlocal(k1, k2, k3, params); break;
		case 11: result = maldacena(k1, k2, k3, params); break;
		case 12: result = string(k1, k2, k3, params); break;
		case 13: result = smooth3(k1, k2, k3, params); break;
		case 14: result = constant(k1, k2, k3, params); break;
		case 15: result = chen(k1, k2, k3, params); break;
		case 16: result = orthogonal(k1, k2, k3, params); break;
		case 17: result = mode1(k1, k2, k3, params); break;
		case 18: result = flat(k1, k2, k3, params); break;
		case 19: result = neil1(k1, k2, k3, params); break;
		case 20: result = neil2(k1, k2, k3, params); break;
		case 21: result = NBD_ivan_cos1(k1, k2, k3, params); break;
		case 22: result = NBD_ivan_cos2(k1, k2, k3, params); break;
		case 23: result = NBD_xingang(k1, k2, k3, params); break;
		case 24: result = sinlog(k1, k2, k3, params); break;
		case 25: result = NBD_sinlog(k1, k2, k3, params); break;
		case 26: result = EFT1(k1, k2, k3, params); break;
		case 27: result = EFT2(k1, k2, k3, params); break;
		case 28: result = featureENV(k1, k2, k3, params); break;
		case 29: result = featureENV2(k1, k2, k3, params); break;
		case 30: result = featureENV3(k1, k2, k3, params); break;
		case 31: result = featureENV4(k1, k2, k3, params); break;
		case 32: result = featureENV5(k1, k2, k3, params); break;
		case 33: result = featureENV6(k1, k2, k3, params); break;
		case 34: result = featureENV7(k1, k2, k3, params); break;
		case 35: result = QSF(k1, k2, k3, params); break;
		case 36: result = adshead(k1, k2, k3, params); break;
		case 37: result = NBD_sin(k1, k2, k3, params); break;
		case 38: result = NBD_ivan_sin1(k1, k2, k3, params); break;
		case 39: result = NBD_ivan_sin2(k1, k2, k3, params); break;
		case 40: result = sinfit(k1, k2, k3, params); break;
		case 41: result = ribero(k1, k2, k3, params); break;
		case 42: result = sinflat(k1, k2, k3, params); break;
		case 43: result = nicola(k1, k2, k3, params); break;
		case 44: result = nicolanum(k1, k2, k3, params); break;
		case 45: result = ksin(k1, k2, k3, params); break;
		case 46: result = sinlogequi(k1, k2, k3, params); break;
		case 47: result = sinlogflat(k1, k2, k3, params); break;
		case 48: result = sinlogtest(k1, k2, k3, params); break;
		case 49: result = sinloglocal(k1, k2, k3, params); break;
		case 50: result = expsin(k1, k2, k3, params); break;
		case 51: result = expcos(k1, k2, k3, params); break;
		case 52: result = equiexpsin(k1, k2, k3, params); break;
		case 53: result = equiexpcos(k1, k2, k3, params); break;
		case 54: result = flatexpsin(k1, k2, k3, params); break;
		case 55: result = flatexpcos(k1, k2, k3, params); break;
		case 56: result = NBD_DBI(k1, k2, k3, params); break;
		case 57: result = cs_feature(k1, k2, k3, params); break;
		case 58: result = cosmo_collider(k1, k2, k3, params); break;
		case 60: result = freqdrift(k1, k2, k3, params); break;
		case 61: result = equicollider(k1, k2, k3, params); break;
		case 62: result = multispeed(k1, k2, k3, params); break;
		case 63: result = quasilocalfield(k1, k2, k3, params); break;
	}
	return result;
}

double shape2(double k1, double k2, double k3, shape_params params) {

	double result;
	
	if(nscalar!=1.0 && (model==2||model==6||model==14||model==16)){
		k1 = pow(k1,(4.0-nscalar)/3.0);
		k2 = pow(k2,(4.0-nscalar)/3.0);
		k3 = pow(k3,(4.0-nscalar)/3.0);
	}
	
	switch (model2) {
		case 1: result = dbi(k1, k2, k3, params); break;
		case 2: result = equilateral(k1, k2, k3, params); break;
		case 3: result = smooth2(k1, k2, k3, params); break;
		case 4: result = feature(k1, k2, k3, params); break;
		case 5: result = ghost(k1, k2, k3, params); break;
		case 6: result = local(k1, k2, k3, params); break;
		case 7: result = single(k1, k2, k3, params); break;
		case 8: result = warm(k1, k2, k3, params); break;
		case 9: result = smooth1(k1, k2, k3, params); break;
		case 10: result = nonlocal(k1, k2, k3, params); break;
		case 11: result = maldacena(k1, k2, k3, params); break;
		case 12: result = string(k1, k2, k3, params); break;
		case 13: result = smooth3(k1, k2, k3, params); break;
		case 14: result = constant(k1, k2, k3, params); break;
		case 15: result = chen(k1, k2, k3, params); break;
		case 16: result = orthogonal(k1, k2, k3, params); break;
		case 17: result = mode1(k1, k2, k3, params); break;
		case 18: result = flat(k1, k2, k3, params); break;
		case 19: result = neil1(k1, k2, k3, params); break;
		case 20: result = neil2(k1, k2, k3, params); break;
		case 21: result = NBD_ivan_cos1(k1, k2, k3, params); break;
		case 22: result = NBD_ivan_cos2(k1, k2, k3, params); break;
		case 23: result = NBD_xingang(k1, k2, k3, params); break;
		case 24: result = sinlog(k1, k2, k3, params); break;
		case 25: result = NBD_sinlog(k1, k2, k3, params); break;
		case 26: result = EFT1(k1, k2, k3, params); break;
		case 27: result = EFT2(k1, k2, k3, params); break;
		case 28: result = featureENV(k1, k2, k3, params); break;
		case 29: result = featureENV2(k1, k2, k3, params); break;
		case 30: result = featureENV3(k1, k2, k3, params); break;
		case 31: result = featureENV4(k1, k2, k3, params); break;
		case 32: result = featureENV5(k1, k2, k3, params); break;
		case 33: result = featureENV6(k1, k2, k3, params); break;
		case 34: result = featureENV7(k1, k2, k3, params); break;
		case 35: result = QSF(k1, k2, k3, params); break;
		case 36: result = adshead(k1, k2, k3, params); break;
		case 37: result = NBD_sin(k1, k2, k3, params); break;
		case 38: result = NBD_ivan_sin1(k1, k2, k3, params); break;
		case 39: result = NBD_ivan_sin2(k1, k2, k3, params); break;
		case 40: result = sinfit(k1, k2, k3, params); break;
		case 41: result = ribero(k1, k2, k3, params); break;
		case 42: result = sinflat(k1, k2, k3, params); break;
		case 43: result = nicola(k1, k2, k3, params); break;
		case 44: result = nicolanum(k1, k2, k3, params); break;
		case 45: result = ksin(k1, k2, k3, params); break;
		case 46: result = sinlogequi(k1, k2, k3, params); break;
		case 47: result = sinlogflat(k1, k2, k3, params); break;
		case 48: result = sinlogtest(k1, k2, k3, params); break;
		case 49: result = sinloglocal(k1, k2, k3, params); break;
		case 50: result = expsin(k1, k2, k3, params); break;
		case 51: result = expcos(k1, k2, k3, params); break;
		case 52: result = equiexpsin(k1, k2, k3, params); break;
		case 53: result = equiexpcos(k1, k2, k3, params); break;
		case 54: result = flatexpsin(k1, k2, k3, params); break;
		case 55: result = flatexpcos(k1, k2, k3, params); break;
		case 56: result = NBD_DBI(k1, k2, k3, params); break;
		case 57: result = cs_feature(k1, k2, k3, params); break;
		case 60: result = freqdrift(k1, k2, k3, params); break;
	}

	return result;
}

double scale(double k, shape_params params) {

	double k6 = pow(k,6);
	double result = 1.0;
	if(k>1e-5)result = k6*shape(k, k, k, params);
	return result;
}

double slice(double k1, double k2, double k3, shape_params params) {

	double k = (k1+k2+k3)/3e0;
	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p222 = k12*k22*k32;

	double result = 1e0;
	if(scale(k, params)>1e-15)result =  p222*shape(k1, k2, k3, params)/scale(k, params);
	return result;
}

double shape3(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p222 = k12*k22*k32;
	
 	//double k13 = k12*k1;
 	//double k23 = k22*k2;
 	//double k33 = k32*k3;
 	//double p3 = k13+k23+k33;
 	//double p333 = k13*k23*k33;
	//double p111 = k1*k2*k3;
	
	double result;

	result = p222 * shape(k1,k2,k3,params);

	return result;

 	//if(p111==0){
 	//	result = 0.0;
 	//}else{
 	//	result = p222 * shape(k1, k2, k3, params);
 	//}
	/*if (model == 2 || model == 6 || model == 16){
		result = shape(k1, k2, k3, params);
	}
	else result = p222 * shape(k1, k2, k3, params);
	return result;*/
}

double quasilocalfield(double k1, double k2, double k3, shape_params params){
	// Sqsf = 3\sqrt(3kappa) N_v(8kappa)/N_v(8/27) with kappa = k1k2k3/K
	double K = k1 + k2 + k3;
	double k = k1 * k2 * k3;
	double K3inv = 0.0;
	double k2inv = 0.0;
	double result = 0.0;
	double nu = 1.0 * params.a1; //change this as needed to scan different parameter space
	double kappa;
	double arg827 = 8.0 / 27.0;
	if (k > 1e-15){
		k2inv = 1e0 / (k * k);
	}
	if (K > 1e-10){
		K3inv = 1e0 / (K * K * K);
	}
	kappa = k * K3inv;
	// next: try cutoff kappa > 1e-4 and maybe 2 or 5 e-3 and compare...
	if (kappa > 1e-15){
		result = 3 * sqrt(3 * kappa) * gsl_sf_bessel_Yn(nu, 8 * kappa) / gsl_sf_bessel_Yn(nu, arg827);
	}
	result = result * k2inv;
	return result;
	}

double equicollider(double k1, double k2, double k3, shape_params params){
	double phase_step = 10.0;
	double freq_step = 10.0;
	double phase = 3.1415926 / phase_step;
	double u = params.a1 / freq_step;		// this is the frequency $\mu$
	double d = params.a2;		// this will determine the overall phase as a single parameter phi=u*ln(c_sigma/(2*c_s)+delta)
	double phi = phase * d;
	double result;
	double p111 = k1 * k2 * k3;
	double pinv = 0e0;
	double p12, p13, p23, p12inv, p13inv, p23inv;
	double result1, result2, result3;

	p12 = k1 + k2;
	p13 = k1 + k3;
	p23 = k2 + k3;

	if (p111 > 1e-15){
		pinv = 1e0 / (p111 * p111);

		p12inv = 1e0 / p12;
		p13inv = 1e0 / p13;
		p23inv = 1e0 / p23;

		result1 = k1*k2*p12inv*p12inv*sqrt(k3)*sqrt(p12inv) * cos(u * log(k3 * p12inv) + phi);
		result2 = k1*k3*p13inv*p13inv*sqrt(k2)*sqrt(p13inv) * cos(u * log(k2 * p13inv) + phi);
		result3 = k2*k3*p23inv*p23inv*sqrt(k1)*sqrt(p23inv) * cos(u * log(k1 * p23inv) + phi);

		result = result1 + result2 + result3;
	}
	result = result * pinv;
	return result;
}


/*double equicollider_old(double k1, double k2, double k3, shape_params params){

	double phase = 2 * 3.1415926 / 20;
	double u = params.a1;	// this is the factor, \mu
	double d = params.a2;	// this is the phase, delta
	double c = 0.05;	// c_{\sigma} / c_s will be a single parameter == c
	
	double result;
	double p111 = k1 * k2 * k3;
	double pinv = 0e0;
	if (p111 > 1e-15){
		pinv = 1e0 / (p111 * p111);
	}

	double p12 = k1 + k2;
	double p13 = k1 + k3;
	double p23 = k2 + k3;

	double p12inv = 0e0;
	double p13inv = 0e0;
	double p23inv = 0e0;

	double result1 = 0.0;
	double result2 = 0.0;
	double result3 = 0.0;

	if (k1 > 1e-15 && k2 > 1e-15 && k3 > 1e-15){
		if (p12 < 1e-15){
			result = 0.0;
			return result;
		}
		else {
			p12inv = 1e0 / p12;
			//result1 = k1*k2*p12inv*p12inv*sqrt(k3)*sqrt(p12inv) * cos(u * log(0.5 * c * k3 * p12inv) + d*phase);
		}
		if (p13 < 1e-15){
			result = 0.0;
			return result;
		}
		else{
			p13inv = 1e0 / p13;
			//result2 = k1*k3*p13inv*p13inv*sqrt(k2)*sqrt(p13inv) * cos(u * log(0.5 * c * k2 * p13inv) + d*phase);
		} 
		if (p23 < 1e-15){
			result = 0.0;
			return result;
		}
		else{
			p23inv = 1e0 / p23;
			//result3 = k2*k3*p23inv*p23inv*sqrt(k1)*sqrt(p23inv) * cos(u * log(0.5 * c * k1 * p23inv) + d*phase);
		} 

		result1 = k1*k2*p12inv*p12inv*sqrt(k3)*sqrt(p12inv) * cos(u * log(0.5 * c * k3 * p12inv) + d*phase);
		result2 = k1*k3*p13inv*p13inv*sqrt(k2)*sqrt(p13inv) * cos(u * log(0.5 * c * k2 * p13inv) + d*phase);
		result3 = k2*k3*p23inv*p23inv*sqrt(k1)*sqrt(p23inv) * cos(u * log(0.5 * c * k1 * p23inv) + d*phase);

		result = result1 + result2 + result3;
	}
	
	result = result * pinv;
	return result;
}
*/

double multispeed(double k1, double k2, double k3, shape_params params){

	double c1 = params.a1;
	double c2 = params.a2;
	double c3 = params.a3;

	
	if(params.a3 == 0){
		// Effectively we are interested in c1/c3 and c2/c3, so we can just set c3=1
		c3 = 1e0;
	}

	double cubic1inv = 0e0;
	double cubic2inv = 0e0;
	double cubic3inv = 0e0;
	double cubic4inv = 0e0;
	double cubic5inv = 0e0;
	double cubic6inv = 0e0;

	double cubic1 = c1*k1 + c2*k2 + c3*k3;
	cubic1 = cubic1 * cubic1 * cubic1;
	if (cubic1 > 1e-10) cubic1inv = 1e0 / cubic1;

	double cubic2 = c1*k2 + c2*k1 + c3*k3;
	cubic2 = cubic2 * cubic2 * cubic2;
	if (cubic2 > 1e-10) cubic2inv = 1e0 / cubic2;

	double cubic3 = c1*k3 + c2*k2 + c3*k1;
	cubic3 = cubic3 * cubic3 * cubic3;
	if (cubic3 > 1e-10) cubic3inv = 1e0 / cubic3;

	double cubic4 = c1*k1 + c2*k3 + c3*k2;
	cubic4 = cubic4 * cubic4 * cubic4;
	if (cubic4 > 1e-10) cubic4inv = 1e0 / cubic4;

	double cubic5 = c2*k1 + c3*k2 + c1*k3;
	cubic5 = cubic5 * cubic5 * cubic5;
	if (cubic5 > 1e-10) cubic5inv = 1e0 / cubic5;

	double cubic6 = c3*k1 + c1*k2 + c2*k3;
	cubic6 = cubic6 * cubic6 * cubic6;
	if (cubic6 > 1e-10) cubic6inv = 1e0 / cubic6;

	double cubicinv = cubic1inv + cubic2inv + cubic3inv + cubic4inv + cubic5inv + cubic6inv;

	double p111 = k1 * k2 * k3;
	double pinv = 0e0;
	if (p111 > 1e-10) pinv = 1e0 / p111;

	double result = pinv * cubicinv;
	return result;
}


double cosmo_collider(double k1, double k2, double k3, shape_params params){

	double K = k1 + k2 + k3;
	double omega = params.a1;
	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p222 = k12*k22*k32;

	double pinv = 0e0;
	if (p222>1e-10) pinv = 1/p222;

	double result = pinv * sin(K*omega);

	return result;
}

double smooth1_in(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	
	double p222 = k12*k22*k32;
	
	double result = p222*warm_cut(k1, k2, k3, params);
	return result;
}

double smooth2_in(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	
	double p222 = k12*k22*k32;
	
	double result = p222*NBD_cut(k1, k2, k3, params);

	return result;
}

double smooth3_in(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double p1 = (k1+k2+k3)/2.0;
	double p222 = k12*k22*k32;
	double p333 = k13*k23*k33;
	double result;

	if (k1t<=50.1 || k2t<=50.1 || k3t<=50.1) {
		result = 0.0;
	} else {
		result = string(k1, k2, k3, params);
	}

// 	result = 1.0;

	return result;
}

double smooth1(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	
	double p222 = k12*k22*k32;
	double pinv = 0e0;
	if (p222>1e-10) pinv = 1/p222;

	double kk = (k1+k2+k3)/2.0;
	double step = 0.01*kk;
	
	double cn = smooth1_in(k1,k2,k3,params);

	double p1 = 0e0;
	double p2 = 0e0;
	double p3 = 0e0;
	double p4 = 0e0;
	double p5 = 0e0;
	double p6 = 0e0;
	if (k1>step && k1+k3>=k2+2*step) p1 = smooth1_in(k1-step,k2+step,k3,params);
	if (k1>step && k1+k2>=k3+2*step) p2 = smooth1_in(k1-step,k2,k3+step,params);
	if (k2>step && k1+k2>=k3+2*step) p3 = smooth1_in(k1,k2-step,k3+step,params);
	if (k2>step && k2+k3>=k1+2*step) p4 = smooth1_in(k1+step,k2-step,k3,params);
	if (k3>step && k2+k3>=k1+2*step) p5 = smooth1_in(k1+step,k2,k3-step,params);
	if (k3>step && k1+k3>=k2+2*step) p6 = smooth1_in(k1,k2+step,k3-step,params);
	
	double s1 = 0e0;
	double s2 = 0e0;
	double s3 = 0e0;
	double s4 = 0e0;
	double s5 = 0e0;
	double s6 = 0e0;
	if (k1>2*step && (k1+k2>=k3+2*step) && (k1+k3>=k2+2*step) ) s1 = smooth1_in(k1-2*step,k2+step,k3+step,params);
	if (k1>step && k2>step && k1+k2>=k3+4*step) s2 = smooth1_in(k1-step,k2-step,k3+2*step,params);
	if (k2>2*step && (k2+k1>=k3+2*step) && (k2+k3>=k1+2*step) ) s3 = smooth1_in(k1+step,k2-2*step,k3+step,params);
	if (k2>step && k3>step && k2+k3>=k1+4*step) s4 = smooth1_in(k1+2*step,k2-step,k3-step,params);
	if (k3>2*step && (k3+k2>=k1+2*step) && (k3+k1>=k2+2*step) ) s5 = smooth1_in(k1+step,k2+step,k3-2*step,params);
	if (k1>step && k3>step && k1+k3>=k2+4*step) s6 = smooth1_in(k1-step,k2+2*step,k3-step,params);
	
	double r1 = 0e0;
	double r2 = 0e0;
	double r3 = 0e0;
	double r4 = 0e0;
	double r5 = 0e0;
	double r6 = 0e0;
	if (k1>2*step && k1+k3>=k2+4*step) r1 = smooth1_in(k1-2*step,k2+2*step,k3,params);
	if (k1>2*step && k1+k2>=k3+4*step) r2 = smooth1_in(k1-2*step,k2,k3+2*step,params);
	if (k2>2*step && k1+k2>=k3+4*step) r3 = smooth1_in(k1,k2-2*step,k3+2*step,params);
	if (k2>2*step && k2+k3>=k1+4*step) r4 = smooth1_in(k1+2*step,k2-2*step,k3,params);
	if (k3>2*step && k2+k3>=k1+4*step) r5 = smooth1_in(k1+2*step,k2,k3-2*step,params);
	if (k3>2*step && k1+k3>=k2+4*step) r6 = smooth1_in(k1,k2+2*step,k3-2*step,params);
	
	double result = (cn + 0.607*(p1+p2+p3+p4+p5+p6) + 0.223*(s1+s2+s3+s4+s5+s6) + 0.135*(r1+r2+r3+r4+r5+r6))/6.79;
	
 	result *= pinv;

	//printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
	//cn, p1,p2,p3,p4,p5,p6,s1,s2,s3,s4,s5,s6,r1,r2,r3,r4,r5,r6);
	printf("RESULT = %e\n", result);
	return result;
	
}

double smooth2(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	
	double p222 = k12*k22*k32;
	double pinv = 0e0;
	if (p222>1e-10) pinv = 1/p222;

	double kk = (k1+k2+k3)/2.0;
	double step = 0.005*kk;
	
	double cn = smooth2_in(k1,k2,k3,params);

	double p1 = 0e0;
	double p2 = 0e0;
	double p3 = 0e0;
	double p4 = 0e0;
	double p5 = 0e0;
	double p6 = 0e0;
	if (k1>step && k1+k3>=k2+2*step) p1 = smooth2_in(k1-step,k2+step,k3,params);
	if (k1>step && k1+k2>=k3+2*step) p2 = smooth2_in(k1-step,k2,k3+step,params);
	if (k2>step && k1+k2>=k3+2*step) p3 = smooth2_in(k1,k2-step,k3+step,params);
	if (k2>step && k2+k3>=k1+2*step) p4 = smooth2_in(k1+step,k2-step,k3,params);
	if (k3>step && k2+k3>=k1+2*step) p5 = smooth2_in(k1+step,k2,k3-step,params);
	if (k3>step && k1+k3>=k2+2*step) p6 = smooth2_in(k1,k2+step,k3-step,params);
	
	double s1 = 0e0;
	double s2 = 0e0;
	double s3 = 0e0;
	double s4 = 0e0;
	double s5 = 0e0;
	double s6 = 0e0;
	if (k1>2*step && (k1+k2>=k3+2*step) && (k1+k3>=k2+2*step) ) s1 = smooth2_in(k1-2*step,k2+step,k3+step,params);
	if (k1>step && k2>step && k1+k2>=k3+4*step) s2 = smooth2_in(k1-step,k2-step,k3+2*step,params);
	if (k2>2*step && (k2+k1>=k3+2*step) && (k2+k3>=k1+2*step) ) s3 = smooth2_in(k1+step,k2-2*step,k3+step,params);
	if (k2>step && k3>step && k2+k3>=k1+4*step) s4 = smooth2_in(k1+2*step,k2-step,k3-step,params);
	if (k3>2*step && (k3+k2>=k1+2*step) && (k3+k1>=k2+2*step) ) s5 = smooth2_in(k1+step,k2+step,k3-2*step,params);
	if (k1>step && k3>step && k1+k3>=k2+4*step) s6 = smooth2_in(k1-step,k2+2*step,k3-step,params);
	
	double r1 = 0e0;
	double r2 = 0e0;
	double r3 = 0e0;
	double r4 = 0e0;
	double r5 = 0e0;
	double r6 = 0e0;
	if (k1>2*step && k1+k3>=k2+4*step) r1 = smooth2_in(k1-2*step,k2+2*step,k3,params);
	if (k1>2*step && k1+k2>=k3+4*step) r2 = smooth2_in(k1-2*step,k2,k3+2*step,params);
	if (k2>2*step && k1+k2>=k3+4*step) r3 = smooth2_in(k1,k2-2*step,k3+2*step,params);
	if (k2>2*step && k2+k3>=k1+4*step) r4 = smooth2_in(k1+2*step,k2-2*step,k3,params);
	if (k3>2*step && k2+k3>=k1+4*step) r5 = smooth2_in(k1+2*step,k2,k3-2*step,params);
	if (k3>2*step && k1+k3>=k2+4*step) r6 = smooth2_in(k1,k2+2*step,k3-2*step,params);
	
	double result = (cn + 0.607*(p1+p2+p3+p4+p5+p6) + 0.223*(s1+s2+s3+s4+s5+s6) + 0.135*(r1+r2+r3+r4+r5+r6))/6.79;
	
 	result *= pinv;
 
	return result;
}

double smooth3(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	
	double p222 = k12*k22*k32;
	double pinv = 0e0;
	if (p222>1e-10) pinv = 1/p222;
	
	double step = 5;
	
	int max = 10;
	int i,j,k;
	double s,w;
	double sum = 0.0;
	double weight =0.0;
	
	double p1,p2,p3;
	
	for (i = -max; i<max+1; i++){
		for (j = -max; j<max+1; j++){
			for (k = -max; k<max+1; k++){
				
				w = exp(-(i*i+j*j+k*k)/40.0);
				s = 0.0;
				
				p1 = k1 + i*step;
				p2 = k2 + j*step;
				p3 = k3 + k*step;
				if(p1>0e0 && p2>0e0 && p3>0e0 && p1+p2>=p3 && p2+p3>=p1 && p3+p1>=p2) s=w*smooth3_in(p1,p2,p3,params);
				
// 				printf("%d\t%d\t%d\t%e\t%e\n",i,j,k,w,s);
				
				weight += w;
				sum += s;
			}
		}
	}
	
 	double result = sum * pinv / weight;
 
	return result;
}

double dbi(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double k14 = k13*k1;
	double k24 = k23*k2;
	double k34 = k33*k3;
	double k15 = k14*k1;
	double k25 = k24*k2;
	double k35 = k34*k3;
	double p1 = k1+k2+k3;
	double p5 = k15+k25+k35;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double p113 = k13*k2*k3 + k23*k3*k1 + k33*k1*k2;
	double p122 = k12*k22*k3 + k22*k32*k1 + k32*k12*k2;
	double p1inv = 0e0;
	if(p1>1e-10) p1inv = 1e0/p1;
	
	double d1 = k14*(k2+k3) + k24*(k3+k1) + k34*(k1+k2);
	double d2 = k13*(k22+k32) + k23*(k32+k12) + k33*(k12+k22);
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p333;
	
	double result = -3.0*factor*p1inv*p1inv*(p5 + 2*d1 - 3*d2 + 2*p113 - 8*p122) / 7.0;
	return result;
}

double equilateral(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p333;
	double result = factor*k1t*k2t*k3t;

	/*if (p111 > 1e-15){
		result = k1/k2 + k2/k1 + k1/k3 + k3/k1 + k2/k3 + k3/k2 - k12/(k2*k3) - k22/(k1*k3) - k32/(k1*k2) - 2;
	}*/
	return result;
}

double NBD(double k1, double k2, double k3, shape_params params) {

	double result;
	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p1 = k1+k2+k3;
	double p2 = k12+k22+k32;
	double p11t1 = k3*k1+k1*k2-k2*k3;
	double p11t2 = k1*k2+k2*k3-k3*k1;
	double p11t3 = k2*k3+k3*k1-k1*k2;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double k1tinv = 0e0;
	if(k1t>1e-10) k1tinv = 1/k1t;
	double k2tinv = 0e0;
	if(k2t>1e-10) k2tinv = 1/k2t;
	double k3tinv = 0e0;
	if(k3t>1e-10) k3tinv = 1/k3t;

	double i21 = p1*(p1*p1-2*p2);
	double i22 = p111*(3 + 4*p11t1*k1tinv*k1tinv + 4*p11t2*k2tinv*k2tinv + 4*p11t3*k3tinv*k3tinv);
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1/p333;
	
	result = factor*(i21+2*i22)/39.0;
	
	return result;
}

double NBD_cut(double k1, double k2, double k3, shape_params params) {

	double result;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double p1 = (k1+k2+k3)/2.0;

	if (k1t<0.021*p1 || k2t<0.021*p1 || k3t<0.021*p1){
		result = 0e0;
	} else {
		result = flat(k1, k2, k3, params);
	}
	
	return result;
}

double ghost(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double k12t = k22+k32-k12;
	double k22t = k32+k12-k22;
	double k32t = k12+k22-k32;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p333 = k13*k23*k33;
	
	double f1 = 0e0;
	if(p1>1e-10)  f1 = 2*k1/p1;
	double f2 = 0e0;
	if(p1>1e-10)  f2 = 2*k2/p1;
	double f3 = 0e0;
	if(p1>1e-10)  f3 = 2*k3/p1;
	
	double g1 = 0e0;
	if(k1>1e-10) g1 = fabs(k2-k3)/k1;
	double g2 = 0e0;
	if(k2>1e-10) g2 = fabs(k1-k3)/k2;
	double g3 = 0e0;
	if(k3>1e-10) g3 = fabs(k2-k1)/k3;
	
	double int1 = function1(f1)*(1-function1(f1)*function2(g1));
	double int2 = function1(f2)*(1-function1(f2)*function2(g2));
	double int3 = function1(f3)*(1-function1(f3)*function2(g3));
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p333;
	
	double result = factor*(int1*k1*k12t + int2*k2*k22t + int3*k3*k32t)/1.584319e-01;
	
	return result;
}

double local(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p3 = k13+k23+k33;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double factor = 0e0;
	//double result = 0e0;
	if(p111>1e-15) {
		factor = 1e0/p333;
		//result = (k12/(k2*k3) + k22/(k1*k3) + k32/(k1*k2))/3e0;
	}
	double result = factor*p3/3e0;
	return result;
}

double mode1(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p12 = k12*(k2+k3)+k22*(k1+k3)+k32*(k2+k1);
// 	double p3 = k13+k23+k33;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p333;
	
	double result = factor*p12/6.0;
	return result;
}

double single(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p1inv = 0e0;
	if(p1>1e-10) p1inv = 1/p1;
	double p111inv = 0e0;
	if(p111>1e-15) p111inv = 1e0/p111;

	double result = 27e0*p111inv*p1inv*p1inv*p1inv;
	return result;
}

double nicola(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p1inv = 0e0;
	if(p1>1e-10) p1inv = 1/p1;
	double p111inv = 0e0;
	if(p111>1e-15) p111inv = 1/p111;

	
	double result = 3e0*p111inv*p111inv*p1inv*pow(p111,1e0/3e0);
	return result;
}

double warm(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double k12t = k22+k32-k12;
	double k22t = k32+k12-k22;
	double k32t = k12+k22-k32;
	
	double w1 = k13*(k22+k32)*k12t;
	double w2 = k23*(k12+k32)*k22t;
	double w3 = k33*(k12+k22)*k32t;
	double w4 = 0e0;
	if(k22*k32>1e-10) w4 = 1e0/(k22*k32);
	double w5 = 0e0;
	if(k12*k32>1e-10) w5 = 1e0/(k12*k32);
	double w6 = 0e0;
	if(k12*k22>1e-10) w6 = 1e0/(k12*k22);
	
	double factor = 0e0;
	if(p111>1e-10) factor = 1e0/p333;
	
	double result = factor*(w1*w4 + w2*w5 + w3*w6)/6e0;
	
	return result;
}

double warm_cut(double k1, double k2, double k3, shape_params params) {

	double result;
	double p1 = (k1+k2+k3)/2.0;

	if (k1<0.021*p1 || k2<0.021*p1 || k3<0.021*p1){
		result = 0e0;
	} else {
		result = warm(k1, k2, k3, params);
	}

	return result;
}

double nonlocal(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p1 = k1+k2+k3;
	double p3 = k13+k23+k33;
	double p11 = k1*k2 + k2*k3 + k3*k1;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p333;
	
	double result = factor*(60e0*p3 + 3e0*p1*p11-4e0*p1*p1*p1/3e0)/171e0;
	return result;
}

double feature(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;

	double factor = 0e0;
	if(p111>1e-15) factor = 1e0 / (p111*p111);
	
	double s1 = params.a1*p1;
	double s2 = 2e0*3.1415927*params.a2*1e-2;
	double result = factor*sin(s1+s2);
	return result;

}

double feature_slice(double k1, double k2, double k3, shape_params params) {

	return 1.0;

}

double feature_scale(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	
	double s1 = params.a1*p1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double result = sin(s1+s2);
	return result;

}

double sinlog(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = (k1+k2+k3)/3e0;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;

	double s1,s2;
	double result = 0.0;
	if(p1>1e-10){
		s1 = params.a1;
		s2 = 2e0*3.1415927e0*params.a2*1e-2;
		result = factor*sin(s1*log(p1) + s2);
	}
	return result;

}

double sinlog_slice(double k1, double k2, double k3, shape_params params) {

	return 1.0;
	
}

double sinlog_scale(double k1, double k2, double k3, shape_params params) {

	double p1 = (k1+k2+k3)/3e0;
	double s1,s2;
	double result = 0.0;
	if(p1>1e-10){
		s1 = params.a1;
		s2 = 2e0*3.1415927e0*params.a2*1e-2;
		result = sin(s1*log(p1) + s2);
	}
	return result;

}

double sinlogequi(double k1, double k2, double k3, shape_params params) {
	
	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;
	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = (k1+k2+k3)/3e0;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double p333 = p222*p111;
	
	
	double factor = 0e0;
	if(p111>1e-15) factor = (k1t*k2t*k3t)/p333;

	double s1,s2;
	double result = 0.0;
	if(p1>1e-10){
		s1 = params.a1;
		s2 = 2e0*3.1415927e0*params.a2*1e-2;
		result = factor*sin(s1*log(p1) + s2);
	}
	return result;

}

double sinlogflat(double k1, double k2, double k3, shape_params params) {

	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;
	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = (k1+k2+k3)/3e0;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double p333 = p222*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = (1e0 - k1t*k2t*k3t/p111)/p222;

	double s1,s2;
	double result = 0.0;
	if(p1>1e-10){
		s1 = params.a1;
		s2 = 2e0*3.1415927e0*params.a2*1e-2;
		result = factor*sin(s1*log(p1) + s2);
	}
	return result;

}

double sinlogtest(double k1, double k2, double k3, shape_params params) {
	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;
	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = (k1+k2+k3)/3e0;
	double p12 = (k1*k2*(k1+k2) + k1*k3*(k1+k3) + k3*k2*(k3+k2))/6e0;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double p333 = p222*p111;
	
	double factorE = 0e0;
	double factorF = 0e0;
	double factorL = 0e0;
	double omegaE = 35e0;
	double omegaF = 12e0;
	double omegaL = 15e0;
	double ampliE = 525e0;
	double ampliF = -220e0;
	double ampliL = 81.4e0;
	double phaseE = 5e-2*2e0*3.1415927;
	double phaseF = 2.5e-1*2e0*3.1415927;
	double phaseL = 2e-1*2e0*3.1415927;
	if(p111>1e-15){
		factorE = ampliE*(k1t*k2t*k3t)/p333;
		factorF = ampliF*(1e0 - k1t*k2t*k3t/p111)/p222;
		factorL = ampliL*(p12/p111)/p222;
	} 

	double s1,s2,s3;
	double result = 0.0;
	s1 = (p1-0.05)*(p1-0.05)/(1.8e-3);
	if(p1>1e-10){
		// result = exp(-s1)*(factorE*(2e0*sin(36*log(p1)+phase) + sin(18*log(p1))+phase) + factorF*cos(12*log(p1))+phase);
		// result = exp(-s1)*(factorE*sin(omegaE*log(p1)+phaseE) + factorF*sin(omegaF*log(p1)+phaseF) + factorL*sin(omegaL*log(p1))+phaseL);
		result = exp(-s1)*(factorE*sin(omegaE*log(p1)+phaseE) + factorF*sin(omegaF*log(p1)+phaseF));  
		
	}
	return result;
}

double sinloglocal(double k1, double k2, double k3, shape_params params) {

	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;
	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = (k1+k2+k3)/3e0;
	double p12 = (k1*k2*(k1+k2) + k1*k3*(k1+k3) + k3*k2*(k3+k2))/6e0;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double p333 = p222*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = (p12/p111)/p222;

	double s1,s2;
	double result = 0.0;
	if(p1>1e-10){
		s1 = params.a1;
		s2 = 2e0*3.1415927e0*params.a2*1e-2;
		result = factor*sin(s1*log(p1) + s2);
	}
	return result;

}

double maldacena(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p1 = k1+k2+k3;
	double p3 = k13+k23+k33;
	double p12 = k1*k22+k2*k12+k2*k32+k3*k22+k1*k32+k3*k12;
	double p22 = k12*k22+k22*k32+k12*k32;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	
	double p1inv = 0e0;
	if(p1>1e-10) p1inv = 1/p1;
	
	double a = 1;
	double b = 1;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1/p333;
	
	double result = factor*(a*p3 + b*p12 + 8e0*b*p22*p1inv)/14e0;
	return result;

}

double string(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double k14 = k13*k1;
	double k24 = k23*k2;
	double k34 = k33*k3;
	
	double p1 = k1+k2+k3;
	double p4 = k14+k24+k34;
	double p22 = k12*k22+k22*k32+k32*k12;
	double k12t = k22+k32-k12;
	double k22t = k32+k12-k22;
	double k32t = k12+k22-k32;
	double p222 = k12*k22*k32;
	double adelta = 0.0;
	double adelta1 = 0.0;
	double adelta2 = 0.0;
	double adelta3 = 0.0;
	
	
	if(p22/2.0 - p4/4.0>0.0) adelta = sqrt(p22/2.0 - p4/4.0);
// 	if(k1>1e-10) adelta1 = adelta / k1;
// 	if(k2>1e-10) adelta2 = adelta / k2;
// 	if(k3>1e-10) adelta3 = adelta / k3;
	
	double kmin;
	double kmin2;
	
	if(k1>k2 && k1>k3){
		kmin = k1;
	} else if(k2>k1 && k2>k3) {
		kmin = k2;
	} else {
		kmin = k3;
	}
	kmin2 = 500.0;
	
	if(kmin < 500.0){
		kmin2 = kmin;
	}
	
	double zeta;
	
	if(kmin<500.0){
		zeta = 1.0/500.0;
	} else {
		zeta = 1.0/kmin;
	}
	
	if(k1>1e-10) adelta1 = adelta * zeta / k1;
	if(k2>1e-10) adelta2 = adelta * zeta / k2;
	if(k3>1e-10) adelta3 = adelta * zeta / k3;
	
	double s1,s2,s3;
	
	s1 = k12t * function3(adelta1, kmin2) * erf(0.324*zeta*k1);
	s2 = k22t * function3(adelta2, kmin2) * erf(0.324*zeta*k2);
	s3 = k32t * function3(adelta3, kmin2) * erf(0.324*zeta*k3);
	
	double result;
	
	double factor = 0e0;
	if(p222>1e-10) factor = 1/p222;
	
	if(adelta!=0.0) result = -factor*(s1+s2+s3) / (zeta*zeta*adelta);
	
	return result;
	
}

double constant(double k1, double k2, double k3, shape_params params) {

	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;
	
	double result = factor;
	
	return result;
	
}

double chen(double k1, double k2, double k3, shape_params params) {

	double result = 0.0;
	
	return result;
	
}

double QSF(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double x1,x2,x3,x4,x5;
	double result = 0.0;
	
	if(p111>1e-15){
		x1 = 8.0*p111/pow(p1,3);
		x2 = pow(p111*p1/3.0,-1.5);
		x3 = gsl_sf_bessel_Ynu(1e0/params.a1, x1);
		result = x3*x2;
	}
	
	result /= gsl_sf_bessel_Ynu(1e0/params.a1, 8e0/27e0);
	
	return result;
}

double orthogonal(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p3 = (k13+k23+k33);
	double p12 = k1*k22+k2*k12+k2*k32+k3*k22+k1*k32+k3*k12;
	double p111 = k1*k2*k3;
	double p222 = k12*k22*k32;
	double d111 = 0.0;
	double d222 = 0.0;

	if (p111>1e-15){
		d111 = 1e0/p111;
		d222 = 1e0/p222;
	}

	//double result = d222*(3e0*(p12-p3)*d111 - 8e0);
// 	double result = 3.0*equilateral(k1,k2,k3,params)-2.0*constant(k1,k2,k3,params);
	
	double result = (3e0*(p12-p3)*d111 - 8e0);
	return result;
	
}

double flat(double k1, double k2, double k3, shape_params params) {

	double result = constant(k1,k2,k3,params)-equilateral(k1,k2,k3,params);
	
	return result;
	
}

double neil1(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k14 = k12*k12;
	double k24 = k22*k22;
	double k34 = k32*k32;
	double k12t = k22+k32-k12;
	double k22t = k32+k12-k22;
	double k32t = k12+k22-k32;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double p444 = p222*p222;
	double d444 = 0.0;

	if (p111>1e-15){
		d444 = 1e0/p444;
	}

	double result = d444*(k14*k12t+k24*k22t+k34*k32t)/3e0;
	
	return result;
	
}

double neil2(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k15 = k12*k12*k1;
	double k25 = k22*k22*k2;
	double k35 = k32*k32*k3;
	double k12t = k22+k32-k12;
	double k22t = k32+k12-k22;
	double k32t = k12+k22-k32;
	double k12t2 = k12t*k12t;
	double k22t2 = k22t*k22t;
	double k32t2 = k32t*k32t;
	double p111 = k1*k2*k3;;
	double p555 = p111*p111*p111*p111*p111;
	double d555 = 0.0;

	if (p111>1e-15){
		d555 = 1e0/p555;
	}

	double result = 9e0*d555*(k15*k12t2+k25*k22t2+k35*k32t2)/4.0 - local(k1,k2,k3,params);
	result *= 4e0/23e0;
	
	return result;
	
}

double NBD_ivan_cos1(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double d333 = 0.0;

	if (p111>1e-15){
		d333 = 1e0/p333;
	}
	
	double s1,s2,s3;
	
	if(k1t<1e-3){
		s1 = 0.0;
	}else{
		s1 = k22*k32*(1.0-cos(k1t*params.a1))/k1t;
	}
	if(k2t<1e-3){
		s2 = 0.0;
	}else{
		s2 = k12*k32*(1.0-cos(k2t*params.a1))/k2t;
	}
	if(k3t<1e-3){
		s3 = 0.0;
	}else{
		s3 = k22*k12*(1e0-cos(k3t*params.a1))/k3t;
	}
	
	double result = d333*(s1+s2+s3)/3e0;
	
	return result;
	
}

double NBD_ivan_cos2(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double d333 = 0.0;

	if (p111>1e-15){
		d333 = 1.0/p333;
	}
	
	double s1,s2,s3;
	
	if(k1t<1e-3){
		s1 = 0.0;
	}else{
		s1 = k12*(k22+k32)*(1.0-cos(k1t*params.a1))/k1t;
	}
	if(k2t<1e-3){
		s2 = 0.0;
	}else{
		s2 = k22*(k32+k12)*(1.0-cos(k2t*params.a1))/k2t;
	}
	if(k3t<1e-3){
		s3 = 0.0;
	}else{
		s3 = k32*(k12+k22)*(1.0-cos(k3t*params.a1))/k3t;
	}
	
	double result = d333*(s1+s2+s3)/6.0;
	
	return result;
	
}

double NBD_ivan_sin1(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double d333 = 0.0;
	double c1 = 2e0*3.1415927*params.a1;
	double c2 = 2e0*3.1415927*params.a2*1e-2;

	if (p111>1e-15){
		d333 = 1e0/p333;
	}
	
	double s1,s2,s3;
	
	if(k1t<1e-3){
		s1 = 0e0;
	}else{
		s1 = k22*k32*(sin(k1t*c1+c2))/k1t;
	}
	if(k2t<1e-3){
		s2 = 0e0;
	}else{
		s2 = k12*k32*(sin(k2t*c1+c2))/k2t;
	}
	if(k3t<1e-3){
		s3 = 0e0;
	}else{
		s3 = k22*k12*(sin(k3t*c1+c2))/k3t;
	}
	
	double result = d333*(s1+s2+s3)/3e0;
	
	return result;
	
}

double NBD_ivan_sin2(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double d333 = 0.0;
	double c1 = 2.0*3.1415927*params.a1;
	double c2 = 2.0*3.1415927*params.a2*1e-2;

	if (p111>1e-15){
		d333 = 1.0/p333;
	}
	
	double s1,s2,s3;
	
	if(k1t<1e-3){
		s1 = 0.0;
	}else{
		s1 = k12*(k22+k32)*(sin(k1t*c1+c2))/k1t;
	}
	if(k2t<1e-3){
		s2 = 0.0;
	}else{
		s2 = k22*(k32+k12)*(sin(k2t*c1+c2))/k2t;
	}
	if(k3t<1e-3){
		s3 = 0.0;
	}else{
		s3 = k32*(k12+k22)*(sin(k3t*c1+c2))/k3t;
	}
	
	double result = d333*(s1+s2+s3)/3.0;
	
	return result;
	
}

double NBD_xingang(double k1, double k2, double k3, shape_params params) {

	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double p111 = k1*k2*k3;
	double d111 = 0.0;

	if (p111>1e-15){
		d111 = 1e0/p111;
	}
	
	double s1 = k1t / pow(k1t+0.001,4);
	double s2 = k2t / pow(k2t+0.001,4);
	double s3 = k3t / pow(k3t+0.001,4);
	
	double result = d111*(s1+s2+s3)/3.0;
	
	return result;
	
}

double NBD_sinlog(double k1, double k2, double k3, shape_params params) {

	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double d222 = 0.0;
	double p1t = 0.0;
	double p2t = 0.0;
	double p3t = 0.0;
	double lk1 = 0.0;
	double lk2 = 0.0;
	double lk3 = 0.0;
	double c1 = -pow(params.a1,-0.6)/2e0;
	double c2 = 5e-1/params.a1;
	double c3 = 2e0*3.1415927*params.a2*1e-2;

	if(k1>1e-3){
		p1t = k1t/k1;
		lk1 = 2e0*log(k1);
	}

	if(k2>1e-3){
		p2t = k2t/k2;
		lk2 = 2e0*log(k2);
	}

	if(k3>1e-3){
		p3t = k3t/k3;
		lk3 = 2e0*log(k3);
	}

	if (p111>1e-15){
		d222 = 1e0/p222;
	}
	
	double s1 = exp(c1*p1t)*sin(c2*(p1t+lk1)+c3);
	double s2 = exp(c1*p2t)*sin(c2*(p2t+lk2)+c3);
	double s3 = exp(c1*p3t)*sin(c2*(p3t+lk3)+c3);
	
	double result = d222*(s1+s2+s3)/3e0;
	
	return result;
	
}

double NBD_sin(double k1, double k2, double k3, shape_params params) {

	double k1t = k1-k2-k3;
	double k2t = k2-k1-k3;
	double k3t = k3-k1-k2;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double s1 = params.a1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double d222;

	if (p111>1e-15){
		d222 = 1e0/p222;
	}
	
	double e1 = exp(s1*k1t);
	double e2 = exp(s1*k2t);
	double e3 = exp(s1*k3t);
	
	double result = d222*(e1+e2+e3)*sin(s1*p1+s2)/3.0;
	
	return result;
	
}

double sinfit(double k1, double k2, double k3, shape_params params) {

	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	
	double equi = 0e0;
	if(p111>1e-15){
		equi = (k1t*k2t*k3t)/(k1*k2*k3);
	}

	double factor = 0e0;
	if(p111>1e-15) factor = 1e0 / (p111*p111);
	
	double s1 = params.a1*p1;
	double s2 = 2e0*3.1415927*params.a2*1e-2;
	
	double result = factor*equi*sin(s1+s2);
	
	return result;
}

double sinflat(double k1, double k2, double k3, shape_params params) {


	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double s1 = params.a1;
	double s2 = 2e0*3.1415927e0*params.a2*1e-2;
	
	double d111 = 0e0;
	double d222 = 0e0;
	if(p111>1e-15){
		d111 = 1e0/p111;
		d222 = 1e0/p222;
	}
	
  	double result = d222*(1e0 - d111*k1t*k2t*k3t)*sin(s1*p1+s2);
	
	return result;
}

double EFT1(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p11 = k1*k2+k2*k3+k3*k1;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	double p333 = p222*p111;
	double d3 = 0.0;
	double d333 = 0.0;

	if (p111>1e-15){
		d3 = 1e0/(p1*p1*p1);
		d333 = 1e0/p333;
	}
	
	double result = d333*d3*(24e0*p222+p1*(-8e0*p11*p111+p1*(-8e0*p11*p11+p1*(22e0*p111+p1*(-6e0*p11+2e0*p1*p1)))));
	result *= -9e0/34e0;
	return result;
	
}

double EFT2(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double d1 = 0.0;
	double d111 = 0.0;

	if (p111>1e-15){
		d1 = 1.0/p1;
		d111 = 1.0/p111;
	}
	
	double result = 27.0*d111*d1*d1*d1;
	
	return result;
	
}

double featureENV(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;
	
	double s1 = 2.0*3.1415927*p1*params.a1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double result = factor*exp(-(p1-4.5e-2)*(p1-4.5e-2)/(4e-3))*sin(s1+s2);
	return result;

}

double featureENV2(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;
	
	double s1 = 2.0*3.1415927*p1*params.a1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double result = factor*exp(-(p1-4.5e-2)*(p1-4.5e-2)/(7.11e-3))*sin(s1+s2);
	return result;

}

double featureENV3(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;
	
	double s1 = 2.0*3.1415927*p1*params.a1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double result = factor*exp(-(p1-4.5e-2)*(p1-4.5e-2)/(1.11e-2))*sin(s1+s2);
	return result;

}

double featureENV4(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;
	
	double s1 = 2.0*3.1415927*p1*params.a1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double result = factor*exp(-(p1-4.5e-2)*(p1-4.5e-2)/(1.6e-2))*sin(s1+s2);
	return result;

}

double featureENV5(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;
	
	double s1 = 2.0*3.1415927*p1*params.a1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double result = factor*exp(-(p1-4.5e-2)*(p1-4.5e-2)/(2.18e-2))*sin(s1+s2);
	return result;

}

double featureENV6(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;
	
	double s1 = 2.0*3.1415927*p1*params.a1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double result = factor*exp(-(p1-4.5e-2)*(p1-4.5e-2)/(2.84e-2))*sin(s1+s2);
	return result;

}

double featureENV7(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double factor = 0e0;
	if(p111>1e-15) factor = 1e0/p222;
	
	double s1 = 2.0*3.1415927*p1*params.a1;
	double s2 = 2.0*3.1415927*params.a2*1e-2;
	double result = factor*exp(-(p1-4.5e-2)*(p1-4.5e-2)/(3.6e-2))*sin(s1+s2);
	return result;

}

double nicola529(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p1 = k1+k2+k3;
	double p3 = k13+k23+k33;
	double p22 = k12*k22+k22*k32+k32*k12;
	double p23 = k12*(k23+k33)+k22*(k33+k13)+k32*(k13+k23);
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double d1 = 0.0;
	
	double factor = 0e0;
	if(p111>1e-15){
		factor = 1e0/p333;
		d1 = 1e0/p1;
	}
	
	double result = -24e0*factor*(p3/8e0 + d1*d1*p23/2e0 - d1*p22)/31e0;
	return result;

}

double ribero(double k1, double k2, double k3, shape_params params) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double d1 = 0.0;
	
	double factor = 0e0;
	if(p111>1e-15){
		factor = 1e0/p333;
// 		d1 = 3e0/p1;
		d1 = 2e1;
	}
	
	double s1 = (k23+k33)*pow(k1*d1,params.a1);
	double s2 = (k13+k33)*pow(k2*d1,params.a1);
	double s3 = (k33+k13)*pow(k3*d1,params.a1);
	
	double result = factor*(s1+s2+s3)/6e0;
	return result;

}

double adshead(double k1, double k2, double k3, shape_params params) {
//  equation 48 of 1110.3050
	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double p1 = k1+k2+k3;
	double p2 = k12+k22+k32;
	double p12 = k1*(k22+k32)+k2*(k32+k12)+k3*(k12+k22);
	double p111 = k1*k2*k3;
	double p333 = p111*p111*p111;
	double d333 = 0e0;

	double s1=0e0;
	double s2=0e0;
	double s3=0e0;
	double s4=0e0;
	
	double t1=params.a1;
	
	if(p111>1e-15){
		d333 = 1e0/p333;
		s1 = adshead_power(k1*t1,params)*adshead_power(k2*t1,params)*adshead_power(k3*t1,params);
		s2 = adshead_i0(p1*t1,params)*p111;
		s3 = adshead_i1(p1*t1,params)*p12;
		s4 = adshead_i2(p1*t1,params)*p1*p2;
	}

	return d333*s1*(-s2-s3+s4);
// 	return d333*s4;
}

double adshead_power(double x,shape_params params){
	double s1,s2,s3;
	double A = 0.0509/6e0;
	s1 = ((-3e0*x*x+9e0)/(x*x))*cos(2e0*x);
	s2 = ((15e0*x*x-9e0)/(x*x*x))*sin(2e0*x);
	s3 = adshead_D(2e0*x,params);
	return exp(A*(s1+s2)*s3);
}

double adshead_i0(double x,shape_params params){
	double s1,s2,s3;
	s1 = -((x*x*x*x-9e0*x*x+54e0)/(x*x))*cos(x);
	s2 = ((2e0*x*x*x*x-27e0*x*x+54e0)/(x*x*x))*sin(x);
	s3 = adshead_D(x,params);
	return (s1+s2)*s3;
}

double adshead_i1(double x,shape_params params){
	double s1,s2,s3;
	s1 = (3e0*(x*x-6e0)/(x*x))*cos(x);
	s2 = ((x*x-6e0)*(x*x-3e0)/(x*x*x))*sin(x);
	s3 = adshead_D(x,params);
	return (s1+s2)*s3;
}

double adshead_i2(double x,shape_params params){
	double s1,s2,s3;
	s1 = -((x*x-9e0)/(x*x))*cos(x);
	s2 = ((4e0*x*x-9e0)/(x*x*x))*sin(x);
	s3 = adshead_D(x,params);
	return (s1+s2)*s3;
}

double adshead_D(double x,shape_params params){
	double y;
	double result = 1e0;
	if(params.a2>0e0){
		y = (x*params.a2)/(params.a1);
		result = y/sinh(y);
	}
	return result;
}

double function1(double x){

	double x2 = x*x;
	double x3 = x2*x;
	double x4 = x3*x;
	double x5 = x4*x;
	
	double c1 = 7841.0/2000000.0;
	double c2 = 195971.0/2400000.0;
	double c3 = 1049.0/4800.0;
	double c4 = 6457.0/19200.0;
	double c5 = 71.0/600.0;
	
	return c1*x + c2*x2 + c3*x3 - c4*x4 + c5*x5;
}

double function2(double x){

	double x2 = x*x;
	double x3 = x2*x;
	double x4 = x3*x;
	double x5 = x4*x;
	
	double c1 = 4278907.0/12000000.0;
	double c2 = 34691.0/4800000.0;
	double c3 = 4509229.0/384000.0;
	double c4 = 3424133.0/192000.0;
	double c5 = 216413.0/25600.0;
	
	return c1*x + c2*x2 + c3*x3 - c4*x4 + c5*x5;
}

double function3(double ksum, double kmin2){

	double result;
	
	if(ksum<1.0e-5){
		result = 0.0;
// 	} else if(ksum>1.0) {
// 		result = sqrt(M_PI)*erf(ksum*0.324);
	} else {
		result = (0.4881*ksum + 0.01946 / (ksum)) * sqrt(kmin2/500.0);
	}
	
	return result;
}

void loadnicola(){
	
	double* data = malloc( sizeof(double)*MAXLINES*3);
	int* len = malloc( sizeof(int));
	char filename[100];
	filename[0] = '\0';
	strcat(filename, "/fast/space/projects/planck/jf334.private/Polarisation/NicolaShape44.txt");
	load_txt_dbl(filename, 3, data, len);
	int size = *len;
	int i,j,m,n;
	double x1,x2,x3;
	
	NicolaShape44 = (double **)create_array(41,41);

	for (i=0; i<41; i++){
		for (j=0; j<41; j++){
			NicolaShape44[i][j] = 0e0;
		}
	}
	
	j=0;
	for (i=0; i<size; i++){
		x1 = data[j++];
		x2 = data[j++];
		x3 = data[j++];
		m = (int)round(40e0*x1);
		n = (int)round(40e0*x2);
		NicolaShape44[m][n] = x3;
		NicolaShape44[n][m] = x3;
	}
	
	for (i=0; i<40; i++){
		NicolaShape44[i][39-i] = NicolaShape44[i][40-i]+NicolaShape44[i+1][39-i]-NicolaShape44[i+1][40-i];
	}
	
	for (i=0; i<39; i++){
		NicolaShape44[i][38-i] = NicolaShape44[i][39-i]+NicolaShape44[i+1][38-i]-NicolaShape44[i+1][39-i];
	}
	
}

double nicolanum(double k1, double k2, double k3, shape_params params) {

	double p111 = k1*k2*k3;
	double p222 = p111*p111;
	
	double x1,x2;
	int i,j;
	double c1,c2,c3,c4;
	double result;
	
	double factor = 0e0;
	if(p111>1e-15){
		factor = 1e0/p222;
	
		if(k1>=k2 && k1>=k3){
			if(k2>=k3){
				x1 = k2/k1;
				x2 = k3/k1;
			}else{
				x1 = k3/k1;
				x2 = k2/k1;
			}
		}else if (k2>=k3){
			if(k1>=k3){
				x1 = k1/k2;
				x2 = k3/k2;			
			}else{
				x1 = k3/k2;
				x2 = k1/k2;			
			}
		}else{
			if(k1>=k2){
				x1 = k1/k3;
				x2 = k2/k3;			
			}else{
				x1 = k2/k3;
				x2 = k1/k3;			
			}
		}
	
		i = (int)floor(40e0*x1);
		j = (int)floor(40e0*x2);
	
		if(i==40)i=39;
		if(j==40)j=39;
		
		c1 = NicolaShape44[i][j];
		c2 = NicolaShape44[i+1][j];
		c3 = NicolaShape44[i][j+1];
		c4 = NicolaShape44[i+1][j+1];
	
		c2 = c2-c1;
		c3 = c3-c1;
		c4 = c4-c1-c2-c3;
		
		x1 = 40e0*x1 - 1e0*i;
		x2 = 40e0*x2 - 1e0*j;
	
		result = c1+c2*x1+c3*x2+c4*x1*x2;
		result *= factor;
	}else{
		result = 0e0;
	}
	
	return result;

}

double ksin(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;

	double factor = 0e0;
	if(p111>1e-15) factor = 1e0 / (p111*p111);
	
	double window = 1e0;
	double y;
	if(params.a2>0e0&&p1>0e0){ 
		y = params.a2*p1;
		window = y/sinh(y);
	}
	
	double s1 = params.a1*p1;
	double result = factor*p1*sin(s1)*window;
	return result;

}


double expsin(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;

	double factor = 0e0;
	if(p111>1e-15) factor = 1e0 / (p111*p111);
	
	double s1 = params.a1*(p1-3e-1*params.a2);
	double window = exp(-s1*s1/32e0);
	
	double result = factor*sin(s1)*window;
	return result;

}

double expcos(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;

	double factor = 0e0;
	if(p111>1e-15) factor = 1e0 / (p111*p111);
	
	double s1 = params.a1*(p1-3e-1*params.a2);
	double window = exp(-s1*s1/32e0);
	
	double result = factor*cos(s1)*window;
	return result;

}

double equiexpsin(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;

	double factor = 0e0;
	double cross = 0e0;
	if(p111>1e-15){
		factor = 1e0 / (p111*p111);
		cross = (k1t*k2t*k3t)/(p111);
	}
	
	double s1 = params.a1*(p1-3e-1*params.a2);
	double window = exp(-s1*s1/32e0);
	
	double result = factor*cross*sin(s1)*window;
	return result;

}

double equiexpcos(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;

	double factor = 0e0;
	double cross = 0e0;
	if(p111>1e-15){
		factor = 1e0 / (p111*p111);
		cross = (k1t*k2t*k3t)/(p111);
	}
	
	double s1 = params.a1*(p1-3e-1*params.a2);
	double window = exp(-s1*s1/32e0);
	
	double result = factor*cross*cos(s1)*window;
	return result;

}

double flatexpsin(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;

	double factor = 0e0;
	double cross = 0e0;
	if(p111>1e-15){
		factor = 1e0 / (p111*p111);
		cross = 1e0 - (k1t*k2t*k3t)/(p111);
	}
	
	double s1 = params.a1*(p1-3e-1*params.a2);
	double window = exp(-s1*s1/32e0);
	
	double result = factor*cross*sin(s1)*window;
	return result;

}

double flatexpcos(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double k1t = -k1+k2+k3;
	double k2t = -k2+k1+k3;
	double k3t = -k3+k1+k2;

	double factor = 0e0;
	double cross = 0e0;
	if(p111>1e-15){
		factor = 1e0 / (p111*p111);
		cross = 1e0 - (k1t*k2t*k3t)/(p111);
	}
	
	double s1 = params.a1*(p1-3e-1*params.a2);
	double window = exp(-s1*s1/32e0);
	
	double result = factor*cross*cos(s1)*window;
	return result;

}

double freqdrift(double k1, double k2, double k3, shape_params params) {

	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double phase = 2.0*3.1415927e0*params.a2/1e2;

	double factor = 0e0;
	if(p111>1e-15){
		factor = 1e0 / (p111*p111);
	}
	
	double s1 = params.a1*(p1 + 2e-1*(params.a3-30)*p1*p1) + phase;
	
	double result = factor*sin(s1);
	return result;

}

double NBD_DBI(double k1, double k2, double k3, shape_params params) {
	// arXiv:0910.4986
	double p1 = k1+k2+k3;
	double p111 = k1*k2*k3;
	double kt1 = -k1+k2+k3;
	double kt2 =  k1-k2+k3;
	double kt3 =  k1+k2-k3;
	double p = 2.0*3.1415927e0*params.a2/1e2;
	double a;
	
	double factor = 0e0;
	if(p111>1e-15){
		factor = 1e0 / (p111);
	}
	
	double part1 = 0e0;
	a = 3e0*params.a1*kt1;
	if(a>1e-5){
		part1 = cos(a + p)/(2e0*a) - sin(a+p)/(a*a) + (cos(p)-cos(a+p))/(a*a*a);
	}else{
		part1 = -sin(p)/2e0;
	}
	
	double part2 = 0e0;
	a = 3e0*params.a1*kt2;
	if(a>1e-5){
		part2 = cos(a + p)/(2e0*a) - sin(a+p)/(a*a) + (cos(p)-cos(a+p))/(a*a*a);
	}else{
		part2 = -sin(p)/2e0;
	}
	
	double part3 = 0e0;
	a = 3e0*params.a1*kt3;
	if(a>1e-5){
		part3 = cos(a + p)/(2e0*a) - sin(a+p)/(a*a) + (cos(p)-cos(a+p))/(a*a*a);
	}else{
		part3 = -sin(p)/2e0;
	}
	

	double result = factor*(part1+part2+part3);
	return result;

}

// ******************************************************** Jesus ********************************************************

int cs_feature_prior_check(shape_params params) {

  // Get the parameters
  double log10u = params.a1;
  double log10s = params.a2;
  double tau0 = params.a3;
  double aux, logbeta, minsr;
  aux = pow(10., log10u)/(exp(0.5)+pow(10., log10u));
  logbeta = log(2.*pow((pow(10., log10s)/aux), 2));
  minsr = log10(0.03);

  if (
      // hard limits
      log10u < -4 || log10u > 0 ||
      log10s < -2 || log10s > 0 ||
      tau0 > 0 ||
      // slow-roll limit
      (log10u < minsr && log10s < minsr) ||
      // beta limit
      logbeta < 0 || logbeta > 14) {
    return 0;
  };
  return 1;
}


// Initialisation: get the vector of samples of deltaP/P
void cs_feature_init(shape_params params, int *n_total, double **k, double **dpp) {

  // Get the parameters
  double log10u = params.a1;
  double log10s = params.a2;
  double tau0 = params.a3;
  
  // Prepare the command
  char cmd_args[MAX_LINE_LEN];
  sprintf(cmd_args, "%s %g %g %g", CMD, log10u, log10s, tau0);

  // Prepare the vars and arrays to track and store the output
  int i = 0, status;
  char line[MAX_LINE_LEN];
  int n_guess = 100; // Initial guess for the output size.
                     // If too small, arrays are doubled on the fly.
  double *tmp_k, *tmp_dpp, k_i, dpp_i;
  *k   = (double *)malloc(n_guess*sizeof(double));
  *dpp = (double *)malloc(n_guess*sizeof(double));
  // Launch the process
  FILE *process;
  printf("# Calling command: %s\n", cmd_args);
  process = popen(cmd_args, "r");
  if (process == NULL) {
  printf("ERROR: failed to set up the environment when running '%s'\n", cmd_args);
  }

  // Read the output from the buffer
  while (fgets(line, sizeof(line)-1, process) != NULL) {
      sscanf(line, "%lf %lf", &k_i, &dpp_i);
      if ((i+1) > n_guess) {
        n_guess *= 2;
        tmp_k   = (double *)realloc(*k,   n_guess*sizeof(double));
        tmp_dpp = (double *)realloc(*dpp, n_guess*sizeof(double));
        if (tmp_k == NULL || tmp_dpp == NULL) {
            printf("ERROR: allocation failed!");
        };
        *k   = tmp_k;
        *dpp = tmp_dpp;
      };
      (*k  )[i] = k_i;
      (*dpp)[i] = dpp_i;
      i++;
  };
  *n_total = i;

  // Close the process (or check if it failed)
  status = pclose(process);
  if (status != 0) {
    printf("ERROR: Something failed along the process.\n");
    printf("Run it by hand to test it: '%s'",cmd_args);
  };

  // Uncomment for tests.
  // Output must be the same as the python script when called directly
  /*
  int j;
  for (j=0; j<i;j++) {printf("%d : %15.8e %15.8e\n", j, (*k)[j], (*dpp)[j]);};
  */
  
};


// Get B(k1,k2,k3). First call initialises the shape 
double cs_feature(double k1, double k2, double k3, shape_params params) {

  // Static set of parameters and spline helpers
  static shape_params params_current = {.a1 = 0, .a2 = 0, .a3 = 0};
  static gsl_interp_accel *acc = NULL;
  static gsl_spline *spline = NULL;
  static double k_min_spline, k_max_spline;
  double result;
  
  double log10u = params.a1;
  double log10s = params.a2;
  double tau0 = params.a3;
  params.a1 = log10u/100e0 - 4e0;
  params.a2 = log10s/100e0 - 2e0;
  params.a3 = -tau0;
  
  // If the parameters changed or this is the 1st call, initialise
  if (params.a1 != params_current.a1 ||
      params.a2 != params_current.a2 ||
      params.a3 != params_current.a3) {
    printf("# Params changed or 1st call -- Initialising...\n");
    // Check the prior constraints, and fail if they fall outside the prior
    if (cs_feature_prior_check(params) == 0) {
      printf("ERROR: Outside prior bounds!\n");
    };
    // Release the memory used by the interpolation, if 1st call
    if (params_current.a1 != 0 &&
        params_current.a2 != 0 &&
        params_current.a3 != 0) {
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
    };
    // Define the new interpolating function
    int n_samples;
    double *k, *dpp;
    cs_feature_init(params, &n_samples, &k, &dpp);
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc(gsl_interp_cspline, n_samples);
    gsl_spline_init(spline, k, dpp, n_samples);
    // Finally, save the parameters so we know the spline is initialised
    params_current = params;
    k_min_spline = k[0];
    k_max_spline = k[n_samples-1];
    free(k);
    free(dpp);
    printf("# Interpolating spline created!\n");
  };
  
  // Evaluate the shape at the given k's
  if(k1<1e-3||k2<1e-3||k3<1e-3){
	  result = 0e0;
  }else{
	  double kt, kto2, ki2, ki12, ki123;
	  kt   = k1+k2+k3;
	  kto2 = kt/2.;  
	  if (kto2 <= k_min_spline || kto2 >= k_max_spline) {
	    return 0;
	  };
	  ki2   = pow(k1,2)+pow(k2,2)+pow(k3,2);
	  ki12  = k1*k2 + k1*k3 + k2*k3;
	  ki123 = k1*k2*k3;
	  double c0, c1, c2;
	  c0 = pow(kt,-2) * (-1.)   * ki12 +
	       pow(kt,-1) * ( 1./4.) * (pow(k1,4)+pow(k2,4)+pow(k3,4))/ki123 +
	       pow(kt,-1) * (-3./2.) * (k1*k2/k3+k1*k3/k2+k2*k3/k1) +
	       1.         * (-5./4.) +
	       pow(kt, 1) * ( 1./4.) * ki12/ki123;
	  c1 = pow(kt,-1) * ( 1./2.) * ki12 +
	       1.         * (-1./8.) * ((ki2-pow(k1,2))/k1 + (ki2-pow(k2,2))/k2 + (ki2-pow(k3,2))/k3);
	  c2 = ki2        * (1./16.);

	  result = c0 * gsl_spline_eval(spline, kto2, acc) + c1 * gsl_spline_eval_deriv(spline, kto2, acc) + c2 * gsl_spline_eval_deriv2(spline, kto2, acc);
  }
  
  return result;    
};





