#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "global.h"

struct geo_params{int geo_l1; int geo_l2; int geo_l3;};

// Declare external functions

double geo_pt(double x, void *p){

	struct geo_params *params = (struct geo_params *)p;
	int geo_l1 = (params->geo_l1);
	int geo_l2 = (params->geo_l2);
	int geo_l3 = (params->geo_l3);
	
	double result = 0.0;
	double p1,p2,p3;
	
	p1 = gsl_sf_legendre_Pl(geo_l1,x);
	p2 = gsl_sf_legendre_Pl(geo_l2,x);
	p3 = gsl_sf_legendre_Pl(geo_l3,x);
	
	result = p1*p2*p3;
	return result;
}

// double calculate_geometric(int l1, int l2, int l3){
// 
// 	double result = 0.0;
// 	int npt = (l1+l2+l3)/2 + 1;
// 	gsl_function F;
// 	struct geo_params params;
// 	F.function = &geo_pt;
// 	F.params = &params;
// 	params.geo_l1 = l1;
// 	params.geo_l2 = l2;
// 	params.geo_l3 = l3;
// 
// 	int test;
// 
// 	test = (l1+l2+l3) - 2*((l1+l2+l3)/2);
// 
// 	if(!test && l1+l2>=l3 && l2+l3>=l1 && l3+l1>=l2){
// 		gsl_integration_glfixed_table *table =  gsl_integration_glfixed_table_alloc(npt);
// 		result = gsl_integration_glfixed(&F, -1.0, 1.0, table);
// 		gsl_integration_glfixed_table_free(table);
// 	} else {
// 		result = 0.0;
// 	}
// 
// 	result = (2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)*result / 8.0*M_PI;
// 
// 	return result;
// 
// }


double calculate_geometric(int l1, int l2, int l3){

	double result;
	double s1,s3,s6,s13,s23,s33,s16,s26,s36;
	int sign;
	double fact,p1,p2,p3;
	int test;

	test = (l1+l2+l3)%2;

	if(!test && l1+l2>=l3 && l2+l3>=l1 && l3+l1>=l2){
		s1 = l1+l2+l3+1.0;
		s3 = l1+l2+l3+1.0/3.0;
		s6 = l1+l2+l3+1.0/6.0;
		s13 = l2+l3-l1+1.0/3.0;
		s23 = l3+l1-l2+1.0/3.0;
		s33 = l1+l2-l3+1.0/3.0;
		s16 = l2+l3-l1+1.0/6.0;
		s26 = l3+l1-l2+1.0/6.0;
		s36 = l1+l2-l3+1.0/6.0;

		fact = M_1_PI * M_1_PI / 2.0;
		p1 = fact*((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)) / s1;
		p2 = s3 / (s13*s23*s33);
		p3 = sqrt(s16*s26*s36/s6);

		result = p1*p2*p3;
	} else {
		result = 0.0;
	}

	return result;

}

double calculate_geometric_tri(int l1, int l2, int l3, int l4){

	double result;
	int i;

	result = 0.0;
	for (i=0;i<l1+l2;i++){
		result+= (calculate_geometric(l1,l2,i)*calculate_geometric(l3,l4,i)) / (double)(2*i+1);
	}
	for (i=0;i<l1+l3;i++){
		result+= (calculate_geometric(l1,l3,i)*calculate_geometric(l2,l4,i)) / (double)(2*i+1);
	}
	if(l1+l4>l2+l3){
		for (i=0;i<l2+l3;i++){
			result+= (calculate_geometric(l2,l3,i)*calculate_geometric(l1,l4,i)) / (double)(2*i+1);
		}
	}else{
		for (i=0;i<l1+l4;i++){
			result+= (calculate_geometric(l1,l4,i)*calculate_geometric(l2,l3,i)) / (double)(2*i+1);
		}
	}
	
	
	return result/3.0;
	
}