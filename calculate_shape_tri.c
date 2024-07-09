
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "global.h"

double shape_tri(double k1, double k2, double k3, double k4) {

	double result;

	k1 = kpivot*pow(k1/kpivot,(4.0-nscalar)/3.0);
	k2 = kpivot*pow(k2/kpivot,(4.0-nscalar)/3.0);
	k3 = kpivot*pow(k3/kpivot,(4.0-nscalar)/3.0);
	k4 = kpivot*pow(k4/kpivot,(4.0-nscalar)/3.0);
	
	switch (model) {
		case 1: result = const_tri(k1, k2, k3, k4); break;
		case 2: result = equi_tri(k1, k2, k3, k4); break;
		case 3: result = localg_tri(k1, k2, k3, k4); break;
		case 4: result = localt_tri(k1, k2, k3, k4); break;
	}

// 	printf("shape\t%d\t%e\t%e\n",model,k1,result);
	return result;
}

double shape_tri3(double k1, double k2, double k3, double k4) {
	
	double result = pow(k1*k2*k3*k4,2.25) * shape_tri(k1, k2, k3, k4);
// 	printf("shape\t%d\t%e\t%e\n",model,k1,result);
	return result;
	
}

double const_tri(double k1, double k2, double k3, double k4){

	double result = 0.0;
	if(k1*k2*k3*k4!=0) result = 1.0/pow(k1*k2*k3*k4,2.25);
	return result;
}

double equi_tri(double k1, double k2, double k3, double k4) {

	double p1 = (k1+k2+k3+k4)/4.0;
	
	double result = 0;
	if(k1*k2*k3*k4!=0) result = 1.0/(k1*k2*k3*k4*pow(p1,5));
	
	return result;
}

double localg_tri(double k1, double k2, double k3, double k4) {

	double k12 = k1*k1;
	double k22 = k2*k2;
	double k32 = k3*k3;
	double k42 = k4*k4;
	double k13 = k12*k1;
	double k23 = k22*k2;
	double k33 = k32*k3;
	double k43 = k42*k4;
	
	double factor = 0;
	if(k1*k2*k3*k4!=0) factor = 1.0/(k13*k23*k33*k43);
	
	double result = factor*(k13+k33+k23+k43)/4.0;
	
	return result;
}

double localt_tri(double k1, double k2, double k3, double k4){

	double result=0.0;
	double k13=k1*k1*k1;
	double k23=k2*k2*k2;
	double k33=k3*k3*k3;
	double k43=k4*k4*k4;
	
	if(k1!=k2&&k1!=k3&&k1!=k4&&k2!=k3&&k2!=k4&&k3!=k4){
		double k13m3=1/(k13*k33);
		double k12m3=1/(k13*k23);
		double k14m3=1/(k13*k43);
		double k23m3=1/(k23*k33);
		double k24m3=1/(k23*k43);
		double k34m3=1/(k33*k43);
		double combi1=k13m3+k23m3+k14m3+k24m3;
		double combi2=k12m3+k23m3+k14m3+k34m3;
		double combi3=k12m3+k24m3+k13m3+k34m3;
		
		double pref1=(1/(k1*k2))*(1/sqrt(k1*k1+k2*k2-2*k1*k2)-1/(k1+k2))+(1/(k3*k4))*(1/sqrt(k3*k3+k4*k4-2*k3*k4)-1/(k3+k4));
		double pref2=(1/(k1*k3))*(1/sqrt(k1*k1+k3*k3-2*k1*k3)-1/(k1+k3))+(1/(k2*k4))*(1/sqrt(k2*k2+k4*k4-2*k2*k4)-1/(k2+k4));
		double pref3=(1/(k1*k4))*(1/sqrt(k1*k1+k4*k4-2*k1*k4)-1/(k1+k4))+(1/(k2*k3))*(1/sqrt(k2*k2+k3*k3-2*k2*k3)-1/(k2+k3));
		
		result=0.25*(pref1*combi1+pref2*combi2+pref3*combi3);//we've 24 terms in this!!!
	}
	return (result/24.0);
}