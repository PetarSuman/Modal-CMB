#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
// #include <nagmk21.h>
#include "global.h"

// double calculate_limber3(int l1, int l2, int l3){
// 	
// 	int i,j,k,m,n,l;
// 	
// 	double k1,k2,k3;
// 	double total_area;
// 	
// 	int size = get_tausize();
// 	
// 	double x[size];
// 	double y[size];
// 	double y1[size];
// 	double y2[size];
// 	double y3[size];
// 	
// 	//get_xvec(x);
// 	
// 	get_tauvec(x);
// 	get_lmvec(l1,y1);
// 	get_lmvec(l2,y2);
// 	get_lmvec(l3,y3);
// 	
// 	x[0] = 0;
// 	y[0] = 0;
// 	
// 	double tau0 = x[size-1];
// 	double xpt;
// 	
// 	
// 	for(i=0;i<size;i++){
// 		xpt = x[i];
// 		if (xpt == 0) {
// 			y[i] = 0;
// 		} else {
// 			k1 = (double)l1 / xpt;
// 			k2 = (double)l2 / xpt;
// 			k3 = (double)l3 / xpt;
// 			y[i] = shape(k1,k2,k3) * y1[i] * y2[i] * y3[i] / (xpt*xpt*xpt*xpt);
// 		}
// 	}
// 
// 	int spl_size = size + 4;
// 	int wrk_size = 6 * size + 16;
// 	double spl_k[spl_size];
// 	double spl_c[spl_size];
// 	double wrk[wrk_size];
// 	int ifail;
// 	
// // 	e01baf_(&size,x,y,spl_k,spl_c,&spl_size,wrk,&wrk_size,&ifail);
// // 	e02bdf_(&spl_size,spl_k,spl_c,&total_area,&ifail);
// 	
// 	return total_area;
// 
// }
// 
// double calculate_limber2(int l1, int l2, int l3){
// 
// 	double total_area = 0;
// 	return total_area;
// 
// }
// 
// double calculate_limber1(int l1, int l2, int l3){
// 
// 	double total_area = 0;
// 	return total_area;
// 	
// }
// 
