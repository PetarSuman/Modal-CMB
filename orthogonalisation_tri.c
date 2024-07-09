#include <math.h>
#include <stdio.h>
#include <mpi.h>
// #include <nagmk21.h>
#include "global.h"

double calculate_ortho_tri(int r,int s){
 
 	double integral=0;
	double ****points = create_4Darray(2,2,2,2);
 	
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	double perm;
	int i,j,k,l,m,n,o,p;
	
	double k1,k2,k3,k4,ksum;
	double grid = 1e0/((double)alpha_points);
	double max = get_kmax_cut();
	double factor = max*grid;
 	double cube_size = factor*factor*factor;
	
	for(i=0;i<alpha_points;i++){
		for(j=i;j<alpha_points;j++){
			for(k=j;k<alpha_points;k++){
				for(l=k;l<alpha_points;l++){
					
					if(l<=i+j+k+2){
						if(i==j){
							if(j==k){
								if(k==l){
									perm=1.0;
								}else{
									perm=4.0;
								}
							}else{
								if(k==l){
									perm=6.0;
								}else{
									perm=12.0;
								}
							}
						}else{
							if(j==k){
								if(k==l){
									perm=4.0;
								}else{
									perm=12.0;
								}
							}else{
								if(k==l){
									perm=12.0;
								}else{
									perm=24.0;
								}
							}
						}
						for(m=0;m<2;m++){
							k1 = factor*(i+m);
							for(n=0;n<2;n++){
								k2 = factor*(j+n);
								for(o=0;o<2;o++){
									k3 = factor*(k+o);
									for(p=0;p<2;p++){
										k4 = factor*(l+p);
										
										points[m][n][o][p] = triQ(r,i+m,j+n,k+o,l+p)*triQ(s,i+m,j+n,k+o,l+p);
										
// 										ksum = k1+k2+k3+k4;
// 										if (ksum>1e-10) ksum = 1e0/ksum;
// 										points[m][n][o][p] = triQ(r,i+m,j+n,k+o,l+p)*triQ(s,i+m,j+n,k+o,l+p)*ksum;
									}
								}
							}
						}
						integral += perm * calculate_volume_tri(i,j,k,l,points) * cube_size;
					}
				}
				
			}
		}
	}
	destroy_4Darray(points);
	return integral;
}

double calculate_check_tri(int r,int s){
 
 	double integral=0;
	double ****points = create_4Darray(2,2,2,2);
 	
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	double perm;
	int i,j,k,l,m,n,o,p;
	
	double k1,k2,k3,k4,ksum;
	double grid = 1e0/((double)alpha_points);
	double max = get_kmax_cut();
	double factor = max*grid;
 	double cube_size = factor*factor*factor;
	
	for(i=0;i<alpha_points;i++){
		for(j=i;j<alpha_points;j++){
			for(k=j;k<alpha_points;k++){
				for(l=k;l<alpha_points;l++){
					
					if(l<=i+j+k+2){
						if(i==j){
							if(j==k){
								if(k==l){
									perm=1.0;
								}else{
									perm=4.0;
								}
							}else{
								if(k==l){
									perm=6.0;
								}else{
									perm=12.0;
								}
							}
						}else{
							if(j==k){
								if(k==l){
									perm=4.0;
								}else{
									perm=12.0;
								}
							}else{
								if(k==l){
									perm=12.0;
								}else{
									perm=24.0;
								}
							}
						}
						for(m=0;m<2;m++){
 							k1 = factor*(i+m);
							for(n=0;n<2;n++){
 								k2 = factor*(j+n);
								for(o=0;o<2;o++){
 									k3 = factor*(k+o);
									for(p=0;p<2;p++){
 										k4 = factor*(l+p);
 										
 										ksum = k1+k2+k3+k4;
 										if (ksum>1e-10) ksum = 1e0/ksum;
 										
										points[m][n][o][p] = triR(r,i+m,j+n,k+o,l+p)*triR(s,i+m,j+n,k+o,l+p)*ksum;
									}
								}
							}
						}
						integral += perm * calculate_volume_tri(i,j,k,l,points) * cube_size;
					}
				}
				
			}
		}
	}
	destroy_4Darray(points);
	return integral;
}

double calculate_orthol_tri_old(int r,int s, int size, int* vec){
 
 	double integral=0;
 	double perm;
	double c[size];

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	int i,j,k,l,m,n,o,p;
 	
	for (i=0; i<size; i++){
		c[i] = 0;
// 		printf("here %e %d\n",get_cl(i),vec[i]);
		if (vec[i]!=0){
			c[i] = pow(2.0*vec[i]+1.0, -1.0/4.0);
		}
	}

	for(i=0;i<size-1;i++){
		for(j=i;j<size-1;j++){
			for(k=j;k<size-1;k++){
				for(l=k;l<size-1;l++){
					
					if(l<=i+j+k+2){
						if(i==j){
							if(j==k){
								if(k==l){
									perm=1.0*c[i]*c[j]*c[k]*c[l];
								}else{
									perm=4.0*c[i]*c[j]*c[k]*c[l];
								}
							}else{
								if(k==l){
									perm=6.0*c[i]*c[j]*c[k]*c[l];
								}else{
									perm=12.0*c[i]*c[j]*c[k]*c[l];
								}
							}
						}else{
							if(j==k){
								if(k==l){
									perm=4.0*c[i]*c[j]*c[k]*c[l];
								}else{
									perm=12.0*c[i]*c[j]*c[k]*c[l];
								}
							}else{
								if(k==l){
									perm=12.0*c[i]*c[j]*c[k]*c[l];
								}else{
									perm=24.0*c[i]*c[j]*c[k]*c[l];
								}
							}
						}
						
						if(vec[i]+vec[j]+vec[k]>=vec[l]&&vec[i]>1){
							integral += perm * calculate_geometric_tri(vec[i], vec[j], vec[k], vec[l]) *triQ(r,i,j,k,l)*triQ(s,i,j,k,l);
						}
					}
				}	
			}
		}
	}
	return integral;
}

double triQ(int n,int q1,int q2,int q3,int q4){
	
	int p1,p2,p3,p4;
	find_perm_tri_prim(n,&p1,&p2,&p3,&p4);
	double result;
 	result = get_basis_tri_prim(q1,p1)*get_basis_tri_prim(q2,p2)*get_basis_tri_prim(q3,p3)*get_basis_tri_prim(q4,p4);
 	result += get_basis_tri_prim(q1,p1)*get_basis_tri_prim(q2,p2)*get_basis_tri_prim(q3,p4)*get_basis_tri_prim(q4,p3);
 	result += get_basis_tri_prim(q1,p1)*get_basis_tri_prim(q2,p3)*get_basis_tri_prim(q3,p2)*get_basis_tri_prim(q4,p4);
 	result += get_basis_tri_prim(q1,p1)*get_basis_tri_prim(q2,p3)*get_basis_tri_prim(q3,p4)*get_basis_tri_prim(q4,p2);
 	result += get_basis_tri_prim(q1,p1)*get_basis_tri_prim(q2,p4)*get_basis_tri_prim(q3,p2)*get_basis_tri_prim(q4,p3);
 	result += get_basis_tri_prim(q1,p1)*get_basis_tri_prim(q2,p4)*get_basis_tri_prim(q3,p3)*get_basis_tri_prim(q4,p2);
 	
 	result += get_basis_tri_prim(q1,p2)*get_basis_tri_prim(q2,p1)*get_basis_tri_prim(q3,p3)*get_basis_tri_prim(q4,p4);
 	result += get_basis_tri_prim(q1,p2)*get_basis_tri_prim(q2,p1)*get_basis_tri_prim(q3,p4)*get_basis_tri_prim(q4,p3);
 	result += get_basis_tri_prim(q1,p2)*get_basis_tri_prim(q2,p3)*get_basis_tri_prim(q3,p1)*get_basis_tri_prim(q4,p4);
 	result += get_basis_tri_prim(q1,p2)*get_basis_tri_prim(q2,p3)*get_basis_tri_prim(q3,p4)*get_basis_tri_prim(q4,p1);
 	result += get_basis_tri_prim(q1,p2)*get_basis_tri_prim(q2,p4)*get_basis_tri_prim(q3,p1)*get_basis_tri_prim(q4,p3);
 	result += get_basis_tri_prim(q1,p2)*get_basis_tri_prim(q2,p4)*get_basis_tri_prim(q3,p3)*get_basis_tri_prim(q4,p1);
 	
 	result += get_basis_tri_prim(q1,p3)*get_basis_tri_prim(q2,p2)*get_basis_tri_prim(q3,p1)*get_basis_tri_prim(q4,p4);
 	result += get_basis_tri_prim(q1,p3)*get_basis_tri_prim(q2,p2)*get_basis_tri_prim(q3,p4)*get_basis_tri_prim(q4,p1);
 	result += get_basis_tri_prim(q1,p3)*get_basis_tri_prim(q2,p1)*get_basis_tri_prim(q3,p2)*get_basis_tri_prim(q4,p4);
 	result += get_basis_tri_prim(q1,p3)*get_basis_tri_prim(q2,p1)*get_basis_tri_prim(q3,p4)*get_basis_tri_prim(q4,p2);
 	result += get_basis_tri_prim(q1,p3)*get_basis_tri_prim(q2,p4)*get_basis_tri_prim(q3,p2)*get_basis_tri_prim(q4,p1);
 	result += get_basis_tri_prim(q1,p3)*get_basis_tri_prim(q2,p4)*get_basis_tri_prim(q3,p1)*get_basis_tri_prim(q4,p2);
 	
 	result += get_basis_tri_prim(q1,p4)*get_basis_tri_prim(q2,p2)*get_basis_tri_prim(q3,p3)*get_basis_tri_prim(q4,p1);
 	result += get_basis_tri_prim(q1,p4)*get_basis_tri_prim(q2,p2)*get_basis_tri_prim(q3,p1)*get_basis_tri_prim(q4,p3);
 	result += get_basis_tri_prim(q1,p4)*get_basis_tri_prim(q2,p3)*get_basis_tri_prim(q3,p2)*get_basis_tri_prim(q4,p1);
 	result += get_basis_tri_prim(q1,p4)*get_basis_tri_prim(q2,p3)*get_basis_tri_prim(q3,p1)*get_basis_tri_prim(q4,p2);
 	result += get_basis_tri_prim(q1,p4)*get_basis_tri_prim(q2,p1)*get_basis_tri_prim(q3,p2)*get_basis_tri_prim(q4,p3);
 	result += get_basis_tri_prim(q1,p4)*get_basis_tri_prim(q2,p1)*get_basis_tri_prim(q3,p3)*get_basis_tri_prim(q4,p2);

	result = result / 24.0;

	return result;
}

double triR(int r,int i,int j,int k,int l){

	int n;
	double result = 0;
	
	for(n=0;n<r+1;n++){
		result += get_ortho_tri(r,n)*triQ(n,i,j,k,l);
	}
 	
	return result;
}

double triS(int i,int j,int k,int l){

	int n;
	double result = 0;
	int s=get_terms_prim();
	
	for(n=0;n<s;n++){
		result += get_eigen_tri(n)*triQ(n,i,j,k,l);
	}
 	
	return result;
}

double triT(int i,int j,int k,int l){

	int n;
	double result = 0;
	int s=get_terms_late();
	
	for(n=0;n<s;n++){
		result += get_modes_tri(n)*triQ(n,i,j,k,l);
	}
 	
	return result;
}

