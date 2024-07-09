#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "global.h"

// Declare external functions
static int npt;
static double *gl_pts;
static double *gl_wgt;
static double *sp_x;
static double **gl_Pl;
static double *lweight;
static double **basis;
static double **power;

double init_ISW_glint(){

	int lmax = eflag_T_lmax;
	npt = 3*lmax/2 + 1;
//	npt = 10001;
		
	int i,j,l;
	double sp_x[npt];
	
//	gsl_integration_glfixed_table *table =  gsl_integration_glfixed_table_alloc(npt);
	
	gl_pts = (double *)create_vector(npt);
	gl_wgt = (double *)create_vector(npt);

	asy(gl_pts,gl_wgt,npt);	
//	for(i=0;i<npt;i++){
//		gsl_integration_glfixed_point(-1.0, 1.0, i, &gl_pts[i], &gl_wgt[i], table);
//	}
	
	gl_Pl = (double **)create_array(npt,lmax+1);
	
//	for(i=0;i<npt;i++){
//		gl_pts[i] = -1.0 + 2.0*i/(npt-1.0);
//	}
	
	for(l=0;l<lmax+1;l++){
		for(i=0;i<npt;i++){
			gl_Pl[i][l] = gsl_sf_legendre_Pl(l,gl_pts[i]);
		}
	}
	
	lweight = (double *)create_vector(lmax+1);
	double s1,s2,s3;
		
	
	for(l=0;l<lmax+1;l++){
		if(l > 1){
			s1 = get_beam_TT(l)*pow(2.0*l+1.0,5.0/6.0);
			s2 = get_beam_TT(l)*get_beam_TT(l)*get_cl_TT(l)+get_noise_TT(l);
			s3 = 1.0 / sqrt(s2);
			lweight[l] =  s1*s3;
		} else {
			lweight[l] = 0.0;
		}
	}

	int pmax = get_pmax_late_T();
	basis = (double **)create_array(pmax+1,lmax+1);
	for(l=0;l<eflag_T_lmin;l++){
		for(i=0;i<pmax+1;i++){
			basis[i][l] = 0e0;
		}
	}
	for(l=eflag_T_lmin;l<lmax+1;l++){
		for(i=0;i<pmax+1;i++){
			basis[i][l] = get_basis_late_T(l-eflag_T_lmin,i);
		}
	}

	power = (double **)create_array(6,lmax+1);
	for(l=0;l<lmax+1;l++){
		power[0][l] = 1.0;
		power[1][l] = l*(l+1.0);
		power[2][l] = get_cl_TT(l);
		power[3][l] = get_cl_TP(l);
		power[4][l] = l*(l+1.0)*get_cl_TT(l);
		power[5][l] = l*(l+1.0)*get_cl_TP(l);
// 		printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\n%",l,power[0][l],power[1][l],power[2][l],power[3][l],power[4][l],power[5][l]);
	}

	return 0;
}

double calculate_ISW(int n) {

	int i,l,r,s,j;
	int pmax = get_pmax_prim();
	int lsize = eflag_T_lmax+1;
	double *factor = create_vector(lsize);
	double **Nmap = create_array(3,6);
	int p1,p2,p3,q1,q2,q3;
	double s1,s2,s3;
	double sum1;
	double sum2;
	double x;
	
	double result = 0.0;
	double *scale = create_vector(lsize);
	double sp_x[npt];
	double sp_y[npt];
	
	for(i=0;i<npt;i++){
		sp_x[i] = gl_pts[i];
	}
	
//	gsl_spline* spline =  gsl_spline_alloc (gsl_interp_cspline, npt);
//	gsl_interp_accel* acc = gsl_interp_accel_alloc();

	find_perm_late_TTT(n,&p1,&p2,&p3);

	int pvec[3];

	pvec[0] = p1;
	pvec[1] = p2;
	pvec[2] = p3;
		
// 	for(l=0;l<lsize;l++){
// 		if(l > 1){
// 			s1 = pow(2.0*l+1.0,5.0/6.0);
// 			s2 = get_cl(l)+get_noise(l)/(get_beam(l)*get_beam(l));
// 			s3 = 1.0/sqrt(s2);
// 			scale[l] =  s1*s3;
// 		} else {
// 			scale[l] = 0.0;
// 		}
// 	}
	
	for(j=0;j<npt;j++){
		
		for(l=0;l<lsize;l++){
			factor[l] =  lweight[l]*gl_Pl[j][l];
			
		}
		
		for(r=0;r<3;r++){
			for(s=0;s<6;s++){
				sum1 = 0.0;
				for(l=0;l<lsize;l++){
					x = factor[l]*basis[pvec[r]][l]*power[s][l];
					sum1 += x;
// 					if(n==1&&(j==0||j==1500)&&r!=1)printf("%d\t%d\t%d\t%e\t%e\n",r,s,l,x,sum1);
				}
				Nmap[r][s] = sum1;
			}
		}
		
		sum1 = 0.0;
		
		sum1 += Nmap[0][0]*Nmap[1][2]*Nmap[2][5];
		sum1 += Nmap[1][0]*Nmap[2][2]*Nmap[0][5];
		sum1 += Nmap[2][0]*Nmap[0][2]*Nmap[1][5];
		sum1 += Nmap[2][0]*Nmap[1][2]*Nmap[0][5];
		sum1 += Nmap[1][0]*Nmap[0][2]*Nmap[2][5];
		sum1 += Nmap[0][0]*Nmap[2][2]*Nmap[1][5];
		
		sum1 -= Nmap[0][1]*Nmap[1][2]*Nmap[2][3];
		sum1 -= Nmap[1][1]*Nmap[2][2]*Nmap[0][3];
		sum1 -= Nmap[2][1]*Nmap[0][2]*Nmap[1][3];
		sum1 -= Nmap[2][1]*Nmap[1][2]*Nmap[0][3];
		sum1 -= Nmap[1][1]*Nmap[0][2]*Nmap[2][3];
		sum1 -= Nmap[0][1]*Nmap[2][2]*Nmap[1][3];
		
		sum1 += Nmap[0][0]*Nmap[1][3]*Nmap[2][4];
		sum1 += Nmap[1][0]*Nmap[2][3]*Nmap[0][4];
		sum1 += Nmap[2][0]*Nmap[0][3]*Nmap[1][4];
		sum1 += Nmap[2][0]*Nmap[1][3]*Nmap[0][4];
		sum1 += Nmap[1][0]*Nmap[0][3]*Nmap[2][4];
		sum1 += Nmap[0][0]*Nmap[2][3]*Nmap[1][4];

		result += gl_wgt[j]*sum1;
//		sp_y[j] = sum1;
//		if(n==0)printf("%e\t%e\n",gl_pts[j],sum1);
	}
	
//	for(j=1;j<npt;j++){
//		result += (sp_x[j]-sp_x[j-1])*(sp_y[j]+sp_y[j-1])/2.0;
//	}
	
//	gsl_spline_init(spline,sp_x,sp_y,npt);
//	result = gsl_spline_eval_integ(spline,sp_x[0],sp_x[npt-1],acc);
	
//	gsl_spline_free(spline);
//	gsl_interp_accel_free(acc);
		
	result /= (16.0*M_PI);
	destroy_array(Nmap);
	destroy_vector(factor);
	return result;
}

double calculate_ISW_norm(){

	int l,r,s,j;
	int pmax = get_pmax_prim();
	int lsize = eflag_T_lmax+1;
	double *factor = create_vector(lsize);
	double **Nmap = create_array(6,6);
	int p1,p2,p3,q1,q2,q3;
	double sum1;
	double sum2;
	double s1,s2,s3;
	
	double result = 0.0;
	
	double *scale = create_vector(lsize);
	
		
	for(l=0;l<lsize;l++){
		if(l > 1){
			s1 = (2.0*l+1.0);
			s2 = get_cl_TT(l)+get_noise_TT(l)/(get_beam_TT(l)*get_beam_TT(l));
			s3 = 1.0/s2;
			scale[l] =  s1*s3;
		} else {
			scale[l] = 0.0;
		}
	}
	
	for(j=0;j<npt;j++){
		
		for(l=0;l<lsize;l++){
			factor[l] =  scale[l]*gl_Pl[j][l];
			
		}
		
		for(r=0;r<6;r++){
			for(s=0;s<6;s++){
				sum1 = 0.0;
				for(l=0;l<lsize;l++){
					sum1 += factor[l]*power[r][l]*power[s][l];
				}
				Nmap[r][s] = sum1;
			}
		}

		sum1 = 0.0;
		
		sum1 += Nmap[0][0]*Nmap[3][3]*Nmap[4][4];
		sum1 += Nmap[3][0]*Nmap[4][3]*Nmap[0][4];
		sum1 += Nmap[4][0]*Nmap[0][3]*Nmap[3][4];
		sum1 += Nmap[4][0]*Nmap[3][3]*Nmap[0][4];
		sum1 += Nmap[3][0]*Nmap[0][3]*Nmap[4][4];
		sum1 += Nmap[0][0]*Nmap[4][3]*Nmap[3][4];
		
		sum1 += Nmap[0][0]*Nmap[2][2]*Nmap[5][5];
		sum1 += Nmap[2][0]*Nmap[5][2]*Nmap[0][5];
		sum1 += Nmap[5][0]*Nmap[0][2]*Nmap[2][5];
		sum1 += Nmap[5][0]*Nmap[2][2]*Nmap[0][5];
		sum1 += Nmap[2][0]*Nmap[0][2]*Nmap[5][5];
		sum1 += Nmap[0][0]*Nmap[5][2]*Nmap[2][5];
		
		sum1 += Nmap[1][1]*Nmap[2][2]*Nmap[3][3];
		sum1 += Nmap[2][1]*Nmap[3][2]*Nmap[1][3];
		sum1 += Nmap[3][1]*Nmap[1][2]*Nmap[2][3];
		sum1 += Nmap[3][1]*Nmap[2][2]*Nmap[1][3];
		sum1 += Nmap[2][1]*Nmap[1][2]*Nmap[3][3];
		sum1 += Nmap[1][1]*Nmap[3][2]*Nmap[2][3];
		
		sum1 -= 2*Nmap[0][1]*Nmap[3][2]*Nmap[4][3];
		sum1 -= 2*Nmap[3][1]*Nmap[4][2]*Nmap[0][3];
		sum1 -= 2*Nmap[4][1]*Nmap[0][2]*Nmap[3][3];
		sum1 -= 2*Nmap[4][1]*Nmap[3][2]*Nmap[0][3];
		sum1 -= 2*Nmap[3][1]*Nmap[0][2]*Nmap[4][3];
		sum1 -= 2*Nmap[0][1]*Nmap[4][2]*Nmap[3][3];
		
		sum1 -= 2*Nmap[0][1]*Nmap[2][2]*Nmap[5][3];
		sum1 -= 2*Nmap[2][1]*Nmap[5][2]*Nmap[0][3];
		sum1 -= 2*Nmap[5][1]*Nmap[0][2]*Nmap[2][3];
		sum1 -= 2*Nmap[5][1]*Nmap[2][2]*Nmap[0][3];
		sum1 -= 2*Nmap[2][1]*Nmap[0][2]*Nmap[5][3];
		sum1 -= 2*Nmap[0][1]*Nmap[5][2]*Nmap[2][3];
		
		sum1 += 2*Nmap[0][0]*Nmap[2][3]*Nmap[5][4];
		sum1 += 2*Nmap[2][0]*Nmap[5][3]*Nmap[0][4];
		sum1 += 2*Nmap[5][0]*Nmap[0][3]*Nmap[2][4];
		sum1 += 2*Nmap[5][0]*Nmap[2][3]*Nmap[0][4];
		sum1 += 2*Nmap[2][0]*Nmap[0][3]*Nmap[5][4];
		sum1 += 2*Nmap[0][0]*Nmap[5][3]*Nmap[2][4];
		
		result += gl_wgt[j]*sum1;
	}
	
	result *= 3.0/(8.0*M_PI);
	destroy_array(Nmap);
	destroy_vector(factor);
	return result;
}

double calculate_ISW_bispectrum(int i, int j, int k){

        double result = 0.0;

        double p1 = (double)i*(i+1e0);
        double p2 = (double)j*(j+1e0);
        double p3 = (double)k*(k+1e0);
        double tt1 = get_cl_TT(i);
        double tt2 = get_cl_TT(j);
        double tt3 = get_cl_TT(k);
        double tp1 = get_cl_TP(i);
        double tp2 = get_cl_TP(j);
        double tp3 = get_cl_TP(k);
        double f1 = p2+p3-p1;
        double f2 = p1+p3-p2;
        double f3 = p1+p2-p3;

        result = f1*(tp2*tt3+tt2*tp3) + f2*(tp1*tt3+tt1*tp3) + f3*(tp1*tt2+tt1*tp2);

        return result/2.0;
}

/*
	In computations of ISW effect, such as calculate_ISW_3D, the function called is calculate_ISW_bispectrum(i,j,k),
	so I am going to update the function to only depend on these three arguments. The original function with 5 arguments
	is only called by calculate_ISW_norm_3D(), and this function is never called in the pipeline, as far as I can tell,
	so we are safe to do so.

double calculate_ISW_bispectrum(int size, int* vec, int i, int j, int k){

	double result = 0.0;
	
	double p1 = vec[i]*(vec[i]+1);
	double p2 = vec[j]*(vec[j]+1);
	double p3 = vec[k]*(vec[k]+1);
	double tt1 = get_cl_TT(i);
	double tt2 = get_cl_TT(j);
	double tt3 = get_cl_TT(k);
	double tp1 = get_cl_TP(i);
	double tp2 = get_cl_TP(j);
	double tp3 = get_cl_TP(k);
	double f1 = p2+p3-p1;
	double f2 = p1+p3-p2;
	double f3 = p1+p2-p3;

	result = f1*(tp2*tt3+tt2*tp3) + f2*(tp1*tt3+tt1*tp3) + f3*(tp1*tt2+tt1*tp2);

	return result/2.0;
}*/

//===================================================================================================
//===================================================================================================



// This function seem to be obsolete as it is never called in the pipeline...
/*double calculate_ISW_norm_3D(){

	int i,j,k,t1,t2,t3;
	double x;
	int lsize = eflag_T_lmax+1;
	int *lvec = create_ivector(lsize);
	for(i=0;i<lsize;i++){
		lvec[i] = i;
	}
	double result = 0.0;
	double s1,s2,s3;
		
	for(i=2;i<lsize;i++){
		s1 = get_cl_TT(i)+get_noise_TT(i)/(get_beam_TT(i)*get_beam_TT(i));
		for(j=i;j<lsize;j++){
			s2 = get_cl_TT(j)+get_noise_TT(j)/(get_beam_TT(j)*get_beam_TT(j));
			t1 = i+j;
			t2 = j%2;
			t3 = t1%2;
			if(t1>eflag_T_lmax)t1=eflag_T_lmax;
			if(t2==0&&t3==0){
				for(k=j;k<t1+1;k+=2){
					x = calculate_ISW_bispectrum(lsize,lvec,i,j,k);
					s3 = get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k));
					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*x/(s1 * s2 * s3);
				}
			}else if(t2==0&&t3==1){
				for(k=j+1;k<t1+1;k+=2){
					x = calculate_ISW_bispectrum(lsize,lvec,i,j,k);
					s3 = get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k));
					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*x/(s1 * s2 * s3);
				}
			}else if(t2==1&&t3==0){
				for(k=j+1;k<t1+1;k+=2){
					x = calculate_ISW_bispectrum(lsize,lvec,i,j,k);
					s3 = get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k));
					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*x/(s1 * s2 * s3);
				}
			}else if(t2==1&&t3==1){
				for(k=j;k<t1+1;k+=2){
					x = calculate_ISW_bispectrum(lsize,lvec,i,j,k);
					s3 = get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k));
					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*x/(s1 * s2 * s3);
				}
			}
		}
	}
	return result;
}*/

// double calculate_ISW_3D(int n){
//
// 	int i,j,k,t1,t2,t3;
// 	double x,y;
// 	int lsize = eflag_T_lmax+1;
// 	int *lvec = create_ivector(lsize);
// 	for(i=0;i<lsize;i++){
// 		lvec[i] = i;
// 	}
// 	double result = 0.0;
// 	double s1,s2,s3;
//
// 	for(i=2;i<lsize;i++){
// 		s1 = pow(2.0*i+1.0,1.0/3.0)*(get_cl_TT(i)+get_noise_TT(i)/(get_beam_TT(i)*get_beam_TT(i)));
// 		for(j=i;j<lsize;j++){
// 			s2 = pow(2.0*j+1.0,1.0/3.0)*(get_cl_TT(j)+get_noise_TT(j)/(get_beam_TT(j)*get_beam_TT(j)));
// 			t1 = i+j;
// 			t2 = j%2;
// 			t3 = t1%2;
// 			if(t1>eflag_T_lmax)t1=eflag_T_lmax;
// 			if(t2==0&&t3==0){
// 				for(k=j;k<t1+1;k+=2){
// 					x = calculate_ISW_bispectrum(lsize,lvec,i,j,k);
// 					y = plijk_TTT(n,i,j,k);
// 					s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k)));
// 					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*y/sqrt(s1 * s2 * s3);
// 				}
// 			}else if(t2==0&&t3==1){
// 				for(k=j+1;k<t1+1;k+=2){
// 					x = calculate_ISW_bispectrum(lsize,lvec,i,j,k);
// 					y = plijk_TTT(n,i,j,k);
// 					s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k)));
// 					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*y/sqrt(s1 * s2 * s3);
// 				}
// 			}else if(t2==1&&t3==0){
// 				for(k=j+1;k<t1+1;k+=2){
// 					x = calculate_ISW_bispectrum(lsize,lvec,i,j,k);
// 					y = plijk_TTT(n,i,j,k);
// 					s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k)));
// 					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*y/sqrt(s1 * s2 * s3);
// 				}
// 			}else if(t2==1&&t3==1){
// 				for(k=j;k<t1+1;k+=2){
// 					x = calculate_ISW_bispectrum(lsize,lvec,i,j,k);
// 					y = plijk_TTT(n,i,j,k);
// 					s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl_TT(k)+get_noise_TT(k)/(get_beam_TT(k)*get_beam_TT(k)));
// 					result += permsix(i,j,k)*calculate_geometric(i,j,k)*x*y/sqrt(s1 * s2 * s3);
// 				}
// 			}
// 		}
// 	}
// 	return result;
// }