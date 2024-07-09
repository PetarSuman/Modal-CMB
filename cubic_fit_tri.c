#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"

double function(double ****coefficients, double w, double x, double y, double z);

int cubic_interpolation_tri(double ****in, double ****out, int *cell){

	double ****coefficients = create_4Darray(4,4,4,4);
	int i,j,k;
	double x,y,z;
	
	coefficients[0][0][0][0] = in[1][1][1][1];
	
	for (i=0;i<step+1;i++){
		x = cell[0]+(double)i/(double)step;
		for (j=0;j<step+1;j++){
			y = cell[1]+(double)j/(double)step;
			for (k=0;k<step+1;k++){
				z = cell[2]+(double)k/(double)step;
				out[i][j][k] = function(coefficients, w, x, y, z);
			}
		}
	}
	
	
	return 0;	
}

double function(double ***coefficients, double w, double x, double y, double z){

	double pt = 0.0;
	double wvec[4];
	double xvec[4];
	double yvec[4];
	double zvec[4];
	int i,j,k,l;
	
	wvec[0] = 1;
	wvec[1] = w;
	wvec[2] = w*w;
	wvec[3] = w*w*w;
	
	xvec[0] = 1;
	xvec[1] = x;
	xvec[2] = x*x;
	xvec[3] = x*x*x;
	
	yvec[0] = 1;
	yvec[1] = y;
	yvec[2] = y*y;
	yvec[3] = y*y*y;
	
	zvec[0] = 1;
	zvec[1] = z;
	zvec[2] = z*z;
	zvec[3] = z*z*z;
	
	for (i=0;i<4;i++){
		for (j=0;j<4;j++){
			for (k=0;k<4;k++){
				for (l=0;l<4;l++){
					pt += coefficients[i][j][k][l]*wvec[i]*xvec[j]*yvec[k]*zvec[l];
				}
			}
		}
	}
	
	return pt;
}
