#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nagmk21.h>
#include "integrate.h"

int legendre(double **array, int l_size, int *l_values, int alpha_max ){

	int i,j,k,ifail=0;
	
	int mode = 1;
	int m = 0;
	double l_max = l_values[l_size - 1];
	double p[alpha_max+1];
	double x = 0;
	
	for (i=0;i<l_size;i++){
		x = 2*(l_values[i]/l_max) - 1;
		s22aaf_(&mode,&x,&m,&alpha_max,p,&ifail);
		for (j=0;j<alpha_max+1;j++){
			array[i][j] = p[j];
		}
	}
	return 0;
}
