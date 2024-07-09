#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "global.h"

static int npt;
static double *gl_pts;
static double *gl_wgt;
static double **gl_Pl;
static double **basis_T;
static double **basis_E;
static double ***Mij_T;
static double ***Mij_E;
double plijk_TTT(int r, int i, int j, int k);
double plijk_TTE(int r, int i, int j, int k);
double plijk_TEE(int r, int i, int j, int k);
double plijk_EEE(int r, int i, int j, int k);

double init_orthol_glint(){

    //int npt; (this was here written by James) already defined as static int

	int lmax_T = eflag_T_lmax;
	int lmin_T = eflag_T_lmin;
	int lmax_E = eflag_E_lmax;
	int lmin_E = eflag_E_lmin;
	int lsize_T = lmax_T - lmin_T + 1;
	int lsize_E = lmax_E - lmin_E + 1;
	int lmax = lmax_T;
	if(do_polarisation==1){
		if(lmax_E>lmax_T)lmax = lmax_E;
	}
	npt = 3*lmax/2 + 1;
	int i,j,l;
	
	gsl_integration_glfixed_table *table =  gsl_integration_glfixed_table_alloc(npt);
	
	gl_pts = (double *)create_vector(npt);
	gl_wgt = (double *)create_vector(npt);
	
	asy(gl_pts,gl_wgt,npt);
	
	gl_Pl = (double **)create_array(npt,lmax+1);
	double x1;
	for(l=0;l<lmax+1;l++){
		x1 = pow(2e0*l+1e0,2e0/3e0);
		for(i=0;i<npt;i++){
			gl_Pl[i][l] = x1*gsl_sf_legendre_Pl(l,gl_pts[i]);
		}
	}

	int pmax_T = get_pmax_late_T();
	int pmax_E = get_pmax_late_E();
	
	basis_T = (double **)create_array(pmax_T+1,lsize_T);
	for(l=0;l<lsize_T;l++){
		for(i=0;i<pmax_T+1;i++){
			basis_T[i][l] = get_basis_late_T(l,i);
		}
	}
	if(do_polarisation==1){
		pmax_E = get_pmax_late_E();
		basis_E = (double **)create_array(pmax_E+1,lsize_E);
		for(l=0;l<lsize_E;l++){
			for(i=0;i<pmax_E+1;i++){
				basis_E[i][l] = get_basis_late_E(l,i);
			}
		}
	}

	//Initialise M_ij_T matrix
	Mij_T = (double ***)create_3Darray(pmax_T+1,pmax_T+1,npt);
	if(do_polarisation==1){
		Mij_E = (double ***)create_3Darray(pmax_E+1,pmax_E+1,npt);
	}
	for(i=0;i<pmax_T+1;i++){
		for(j=0;j<pmax_T+1;j++){
			for(l=0;l<npt;l++){
				Mij_T[i][j][l] = 0e0;
			}
		}
	}
	if(do_polarisation==1){
		for(i=0;i<pmax_E+1;i++){
			for(j=0;j<pmax_E+1;j++){
				for(l=0;l<npt;l++){
					Mij_E[i][j][l] = 0e0;
				}
			}
		}
	}
	return 0;
}

void calculate_orthoMij_T(int i, int j){
	
	int r,s,l;
	double sum;

	int lmax = eflag_T_lmax;
	int lmin = eflag_T_lmin;
	int lsize = lmax - lmin + 1;
	//printf("npt = %d\t fifxth time\n", npt);
	
	double *vec = create_vector(lsize);
	
	#pragma omp parallel for
	for(l=0;l<lsize;l++){
		vec[l] = basis_T[i][l]*basis_T[j][l];
	}
	
	#pragma omp parallel for private(sum,l,s)
	for(s=0;s<npt;s++){
		sum = 0e0;
		for(l=lmin;l<lmax+1;l++){
			//printf("gl_Pl[%d][%d] = %f\tvec[%d] = %f\n", s, l, gl_Pl[s][l], l-lmin, vec[l-lmin]);
			sum += gl_Pl[s][l]*vec[l-lmin];
			//printf("sum = %f\n", sum);
		}
		Mij_T[i][j][s] = sum;
		//printf("Mij_T[%d][%d][%d] = %f\n", i, j, s, Mij_T[i][j][s]);
	}
	
	return;
}

void calculate_orthoMij_E(int i, int j){
	
	int r,s,l;
	double sum;

	int lmax = eflag_E_lmax;
	int lmin = eflag_E_lmin;
	int lsize = lmax - lmin + 1;
	
	double *vec = create_vector(lsize);
	
	#pragma omp parallel for
	for(l=0;l<lsize;l++){
		vec[l] = basis_E[i][l]*basis_E[j][l];
	}
	
	#pragma omp parallel for private(sum,l,s)
	for(s=0;s<npt;s++){
		sum = 0e0;
		for(l=lmin;l<lmax+1;l++){
			sum += gl_Pl[s][l]*vec[l-lmin];
		}
		Mij_E[i][j][s] = sum;
	}
	
	return;
}

void sync_orthoMij(){
	
	int i,j,k,l,m,n;
	double* Mij_T_send;
	double* Mij_T_recv;
	double* Mij_E_send;
	double* Mij_E_recv;
	int Tsize, Esize;
	int pmax_T,pmax_E;
	pmax_T = get_pmax_late_T();
	if(do_polarisation==1)pmax_E = get_pmax_late_E();
	
	Tsize = (pmax_T+1)*(pmax_T+1)*npt;
	Esize = (pmax_E+1)*(pmax_E+1)*npt;
	Mij_T_send = (double *)malloc( Tsize * sizeof(double) );
	Mij_T_recv = (double *)malloc( Tsize * sizeof(double) );
	if(do_polarisation==1){
		Mij_E_send = (double *)malloc( Esize * sizeof(double) );
		Mij_E_recv = (double *)malloc( Esize * sizeof(double) );
	}

	m=0;
	n=0;
	for(i=0;i<pmax_T+1;i++){
		for(j=0;j<pmax_T+1;j++){
			for(l=0;l<npt;l++){
				Mij_T_send[m] = Mij_T[i][j][l];
				m++;
			}
		}
		}
	if(do_polarisation==1){
		for(i=0;i<pmax_E+1;i++){
			for(j=0;j<pmax_E+1;j++){
				for(l=0;l<npt;l++){
					Mij_E_send[n] = Mij_E[i][j][l];
					n++;
				}
			}
		}
	}
	
	// printf("Tsize=%d\n",Tsize);
	// sleep(5);
	MPI_Allreduce(Mij_T_send,Mij_T_recv,Tsize,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
	if(do_polarisation==1){
		MPI_Allreduce(Mij_E_send,Mij_E_recv,Esize,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
	}
	
	m=0;
	n=0;
	for(i=0;i<pmax_T+1;i++){
		for(j=0;j<pmax_T+1;j++){
			for(l=0;l<npt;l++){
				Mij_T[i][j][l] = Mij_T_recv[m];
				m++;
			}
		}
	}
	if(do_polarisation==1){
		for(i=0;i<pmax_E+1;i++){
			for(j=0;j<pmax_E+1;j++){
				for(l=0;l<npt;l++){
					Mij_E[i][j][l] = Mij_E_recv[n];
					n++;
				}
			}
		}
	}

	free(Mij_T_send);
	free(Mij_T_recv);
	if(do_polarisation==1){
		free(Mij_E_send);
		free(Mij_E_recv);
	}
	
	return;
}

int calculate_orthol(int m, int n, double* result) {

	//printf("npt = %d\t seventh time\n", npt);
	int l,i,r,s;
	int Tm1,Tm2,Tm3,Tm4,Tm5,Tm6;
	int Em1,Em2,Em3,Em4,Em5,Em6;
	int Tn1,Tn2,Tn3,Tn4,Tn5,Tn6;
	int En1,En2,En3,En4,En5,En6;
	double Tsum1,Tsum2,Tsum3,Tsum4,Tsum5,Tsum6,Tsum7,Tsum8,Tsum9,Tsum10,Tsum11,Tsum12,Tsum13,Tsum14;
	double Esum1,Esum2,Esum3,Esum4,Esum5,Esum6,Esum7,Esum8,Esum9,Esum10,Esum11,Esum12,Esum13,Esum14;
	double s1,s2,s3,s4,s5,s6;
	double wgt1,wgt2;
	double sum;

	find_perm_late_TTT(m,&Tm1,&Tm2,&Tm3);
	find_perm_late_TTT(n,&Tn1,&Tn2,&Tn3);

	if(do_polarisation==1){
		find_perm_late_TTE(m,&Tm4,&Tm5,&Em6);
		find_perm_late_TTE(n,&Tn4,&Tn5,&En6);
		find_perm_late_TEE(m,&Tm6,&Em4,&Em5);
		find_perm_late_TEE(n,&Tn6,&En4,&En5);
		find_perm_late_EEE(m,&Em1,&Em2,&Em3);
		find_perm_late_EEE(n,&En1,&En2,&En3);
	}
	
	result[0] = 0e0;
	if(do_polarisation==1){
		result[1] = 0e0;
		result[2] = 0e0;
		result[3] = 0e0;
	}

	#pragma omp parallel
	{
		double r1,r2,r3,r4;
		r1 = 0e0;
		r2 = 0e0;
		r3 = 0e0;
		r4 = 0e0;
		#pragma omp for private(sum,s1,s2,s3,s4,s5,s6,wgt1,wgt2)
		for(i=0;i<npt;i++){

			//printf("Gauss-Legendre Weighting(%d) = %f\n", i, gl_wgt[i]);
			wgt1 = gl_wgt[i]/(48e0*M_PI);
			wgt2 = gl_wgt[i]/(16e0*M_PI);
			
			//printf("Mij_T[%d][%d][%d] = %f\n",Tm1, Tn1, i, Mij_T[Tm1][Tn1][i]);
			s1 = Mij_T[Tm1][Tn1][i]*Mij_T[Tm2][Tn2][i]*Mij_T[Tm3][Tn3][i];
			s2 = Mij_T[Tm2][Tn1][i]*Mij_T[Tm3][Tn2][i]*Mij_T[Tm1][Tn3][i];
			s3 = Mij_T[Tm3][Tn1][i]*Mij_T[Tm1][Tn2][i]*Mij_T[Tm2][Tn3][i];
			s4 = Mij_T[Tm3][Tn1][i]*Mij_T[Tm2][Tn2][i]*Mij_T[Tm1][Tn3][i];
			s5 = Mij_T[Tm2][Tn1][i]*Mij_T[Tm1][Tn2][i]*Mij_T[Tm3][Tn3][i];
			s6 = Mij_T[Tm1][Tn1][i]*Mij_T[Tm3][Tn2][i]*Mij_T[Tm2][Tn3][i];
			//printf("s1 s2 s3 s4 s5 s6 = %f %f %f %f %f %f\n", s1, s2, s3, s4, s5, s6);
			sum = s1+s2+s3+s4+s5+s6;
			r1 += wgt1*sum;
		
			if(do_polarisation==1){
			
				s1 = Mij_T[Tm4][Tn4][i]*Mij_T[Tm5][Tn5][i]*Mij_E[Em6][En6][i];
				s2 = Mij_T[Tm5][Tn4][i]*Mij_T[Tm4][Tn5][i]*Mij_E[Em6][En6][i];
			
				sum = s1+s2;
				r2 += wgt2*sum;
			
				s1 = Mij_T[Tm6][Tn6][i]*Mij_E[Em4][En4][i]*Mij_E[Em5][En5][i];
				s2 = Mij_T[Tm6][Tn6][i]*Mij_E[Em5][En4][i]*Mij_E[Em4][En5][i];
			
				sum = s1+s2;
				r3 += wgt2*sum;
			
				s1 = Mij_E[Em1][En1][i]*Mij_E[Em2][En2][i]*Mij_E[Em3][En3][i];
				s2 = Mij_E[Em2][En1][i]*Mij_E[Em3][En2][i]*Mij_E[Em1][En3][i];
				s3 = Mij_E[Em3][En1][i]*Mij_E[Em1][En2][i]*Mij_E[Em2][En3][i];
				s4 = Mij_E[Em3][En1][i]*Mij_E[Em2][En2][i]*Mij_E[Em1][En3][i];
				s5 = Mij_E[Em2][En1][i]*Mij_E[Em1][En2][i]*Mij_E[Em3][En3][i];
				s6 = Mij_E[Em1][En1][i]*Mij_E[Em3][En2][i]*Mij_E[Em2][En3][i];
			
				sum = s1+s2+s3+s4+s5+s6;
				r4 += wgt1*sum;
			}
		}
		#pragma omp critical
		{
			result[0] += r1;
			if(do_polarisation==1){
				result[1] += r2;
				result[2] += r3;
				result[3] += r4;				
			}
			
		}
	}
	//printf("RESULT[0] = %f\n", result[0]);

	return 0;
}

int calculate_orthol_3D(int m, int n, double* result) {

	int i,j,k,t1,t2,t3;
	double x1,x2,x3,x4,y1,y2,y3,y4;
	int lsize = get_lmax()+1;
	int *lvec = create_ivector(lsize);
	for(i=0;i<lsize;i++){
		lvec[i] = i;
	}
	int lmaxT = eflag_T_lmax;
	int lmaxE = eflag_E_lmax;
	double s1,s2,s3,w1;

	int doTTT, doTTE, doTEE, doEEE;
	
	if(do_polarisation==1){
		result[0] = 0.0;
		result[1] = 0.0;
		result[2] = 0.0;
		result[3] = 0.0;	
		for(i=2;i<lsize;i++){
			s1 = pow(2.0*i+1.0,1.0/3.0);
			if(i<lmaxT){
				doTTT = 0;
				doTTE = 0;
				doTEE = 0;
			}else{
				doTTT = 1;
				doTTE = 1;
				doTEE = 1;
			}
			if(i<lmaxE){
				doEEE = 0;
			}else{
				doEEE = 1;
			}
			for(j=i;j<lsize;j++){
				s2 = pow(2.0*j+1.0,1.0/3.0);
				t1 = i+j;
				t2 = j%2;
				t3 = t1%2;
				if(t1>get_lmax())t1=get_lmax();
				if(j<lmaxT){
					doTTT = 0;
					doTTE = 0;
				}
				if(j<lmaxE){
					doTEE = 0;
					doEEE = 0;
				}
				if(t2==0&&t3==0){
					for(k=j;k<t1+1;k+=2){
						if(k<lmaxT){
							doTTT = 0;
						}
						if(k<lmaxE){
							doTTE = 0;
							doTEE = 0;
							doEEE = 0;
						}
						if(doTTT){
							x1 = plijk_TTT(m,i,j,k);
							y1 = plijk_TTT(n,i,j,k);
						}else{
							x1=0e0;
							y1=0e0;
						}
						if(doTTE){
							x2 = plijk_TTE(m,i,j,k);
							y2 = plijk_TTE(n,i,j,k);
						}else{
							x2=0e0;
							y2=0e0;
						}
						if(doTEE){
							x3 = plijk_TEE(m,i,j,k);
							y3 = plijk_TEE(n,i,j,k);
						}else{
							x3=0e0;
							y3=0e0;
						}
						if(doEEE){
							x4 = plijk_EEE(m,i,j,k);
							y4 = plijk_EEE(n,i,j,k);
						}else{
							x4=0e0;
							y4=0e0;
						}
						s3 = pow(2.0*k+1.0,1.0/3.0);
						w1 = permsix(i,j,k)*calculate_geometric(i,j,k)/(s1*s2*s3);
						result[0] += w1*x1*y1;
						result[1] += w1*x2*y2;
						result[2] += w1*x3*y3;
						result[3] += w1*x4*y4;
					}
				}else if(t2==0&&t3==1){
					for(k=j+1;k<t1+1;k+=2){
						if(k<lmaxT){
							doTTT = 0;
						}
						if(k<lmaxE){
							doTTE = 0;
							doTEE = 0;
							doEEE = 0;
						}
						if(doTTT){
							x1 = plijk_TTT(m,i,j,k);
							y1 = plijk_TTT(n,i,j,k);
						}else{
							x1=0e0;
							y1=0e0;
						}
						if(doTTE){
							x2 = plijk_TTE(m,i,j,k);
							y2 = plijk_TTE(n,i,j,k);
						}else{
							x2=0e0;
							y2=0e0;
						}
						if(doTEE){
							x3 = plijk_TEE(m,i,j,k);
							y3 = plijk_TEE(n,i,j,k);
						}else{
							x3=0e0;
							y3=0e0;
						}
						if(doEEE){
							x4 = plijk_EEE(m,i,j,k);
							y4 = plijk_EEE(n,i,j,k);
						}else{
							x4=0e0;
							y4=0e0;
						}
						s3 = pow(2.0*k+1.0,1.0/3.0);
						w1 = permsix(i,j,k)*calculate_geometric(i,j,k)/(s1*s2*s3);
						result[0] += w1*x1*y1;
						result[1] += w1*x2*y2;
						result[2] += w1*x3*y3;
						result[3] += w1*x4*y4;
					}
				}else if(t2==1&&t3==0){
					for(k=j+1;k<t1+1;k+=2){
						if(k<lmaxT){
							doTTT = 0;
						}
						if(k<lmaxE){
							doTTE = 0;
							doTEE = 0;
							doEEE = 0;
						}
						if(doTTT){
							x1 = plijk_TTT(m,i,j,k);
							y1 = plijk_TTT(n,i,j,k);
						}else{
							x1=0e0;
							y1=0e0;
						}
						if(doTTE){
							x2 = plijk_TTE(m,i,j,k);
							y2 = plijk_TTE(n,i,j,k);
						}else{
							x2=0e0;
							y2=0e0;
						}
						if(doTEE){
							x3 = plijk_TEE(m,i,j,k);
							y3 = plijk_TEE(n,i,j,k);
						}else{
							x3=0e0;
							y3=0e0;
						}
						if(doEEE){
							x4 = plijk_EEE(m,i,j,k);
							y4 = plijk_EEE(n,i,j,k);
						}else{
							x4=0e0;
							y4=0e0;
						}
						s3 = pow(2.0*k+1.0,1.0/3.0);
						w1 = permsix(i,j,k)*calculate_geometric(i,j,k)/(s1*s2*s3);
						result[0] += w1*x1*y1;
						result[1] += w1*x2*y2;
						result[2] += w1*x3*y3;
						result[3] += w1*x4*y4;
					}
				}else if(t2==1&&t3==1){
					for(k=j;k<t1+1;k+=2){
						if(k<lmaxT){
							doTTT = 0;
						}
						if(k<lmaxE){
							doTTE = 0;
							doTEE = 0;
							doEEE = 0;
						}
						if(doTTT){
							x1 = plijk_TTT(m,i,j,k);
							y1 = plijk_TTT(n,i,j,k);
						}else{
							x1=0e0;
							y1=0e0;
						}
						if(doTTE){
							x2 = plijk_TTE(m,i,j,k);
							y2 = plijk_TTE(n,i,j,k);
						}else{
							x2=0e0;
							y2=0e0;
						}
						if(doTEE){
							x3 = plijk_TEE(m,i,j,k);
							y3 = plijk_TEE(n,i,j,k);
						}else{
							x3=0e0;
							y3=0e0;
						}
						if(doEEE){
							x4 = plijk_EEE(m,i,j,k);
							y4 = plijk_EEE(n,i,j,k);
						}else{
							x4=0e0;
							y4=0e0;
						}
						s3 = pow(2.0*k+1.0,1.0/3.0);
						w1 = permsix(i,j,k)*calculate_geometric(i,j,k)/(s1*s2*s3);
						result[0] += w1*x1*y1;
						result[1] += w1*x2*y2;
						result[2] += w1*x3*y3;
						result[3] += w1*x4*y4;
					}
				}
			}
		}
	}else{
		result[0] = 0.0;
		for(i=2;i<lsize;i++){
			s1 = pow(2.0*i+1.0,1.0/3.0);
			for(j=i;j<lsize;j++){
				s2 = pow(2.0*j+1.0,1.0/3.0);
				t1 = i+j;
				t2 = j%2;
				t3 = t1%2;
				if(t1>get_lmax())t1=get_lmax();
				if(t2==0&&t3==0){
					for(k=j;k<t1+1;k+=2){
						x1 = plijk_TTT(m,i,j,k);
						y1 = plijk_TTT(n,i,j,k);
						s3 = pow(2.0*k+1.0,1.0/3.0);
						w1 = permsix(i,j,k)*calculate_geometric(i,j,k)/(s1*s2*s3);
						result[0] += w1*x1*y1;
					}
				}else if(t2==0&&t3==1){
					for(k=j+1;k<t1+1;k+=2){
						x1 = plijk_TTT(m,i,j,k);
						y1 = plijk_TTT(n,i,j,k);
						s3 = pow(2.0*k+1.0,1.0/3.0);
						w1 = permsix(i,j,k)*calculate_geometric(i,j,k)/(s1*s2*s3);
						result[0] += w1*x1*y1;
					}
				}else if(t2==1&&t3==0){
					for(k=j+1;k<t1+1;k+=2){
						x1 = plijk_TTT(m,i,j,k);
						y1 = plijk_TTT(n,i,j,k);
						s3 = pow(2.0*k+1.0,1.0/3.0);
						w1 = permsix(i,j,k)*calculate_geometric(i,j,k)/(s1*s2*s3);
						result[0] += w1*x1*y1;
					}
				}else if(t2==1&&t3==1){
					for(k=j;k<t1+1;k+=2){
						x1 = plijk_TTT(m,i,j,k);
						y2 = plijk_TTT(n,i,j,k);
						s3 = pow(2.0*k+1.0,1.0/3.0);
						w1 = permsix(i,j,k)*calculate_geometric(i,j,k)/(s1*s2*s3);
						result[0] += w1*x1*y1;
					}
				}
			}
		}
	}
	return 0;
}

double plijk_TTT(int r, int l1, int l2, int l3){
	
	int p1,p2,p3;
	find_perm_late_TTT(r,&p1,&p2,&p3);
	double result;
	double b1,b2,b3,b4,b5,b6;
	int i,j,k;
	if(l1<eflag_T_lmin||l1>eflag_T_lmax||l2<eflag_T_lmin||l2>eflag_T_lmax||l3<eflag_T_lmin||l3>eflag_T_lmax){
		result = 0e0;
	}else{
		i = l1-eflag_T_lmin;
		j = l2-eflag_T_lmin;
		k = l3-eflag_T_lmin;
	 	b1 = get_basis_late_T(i,p1)*get_basis_late_T(j,p2)*get_basis_late_T(k,p3);
	 	b2 = get_basis_late_T(i,p2)*get_basis_late_T(j,p3)*get_basis_late_T(k,p1);
	 	b3 = get_basis_late_T(i,p3)*get_basis_late_T(j,p1)*get_basis_late_T(k,p2);
	 	b4 = get_basis_late_T(i,p3)*get_basis_late_T(j,p2)*get_basis_late_T(k,p1);
	 	b5 = get_basis_late_T(i,p2)*get_basis_late_T(j,p1)*get_basis_late_T(k,p3);
	 	b6 = get_basis_late_T(i,p1)*get_basis_late_T(j,p3)*get_basis_late_T(k,p2);
	
	 	result = (b1+b2+b3+b4+b5+b6)/(6e0);
	}
 	
	return result;
}

double plijk_TTE(int r, int l1, int l2, int l3){
	
	int p1,p2,p3;
	find_perm_late_TTE(r,&p1,&p2,&p3);
	double result;
	double b1,b2;
	int i,j,k;
	
	if(l1<eflag_T_lmin||l1>eflag_T_lmax||l2<eflag_T_lmin||l2>eflag_T_lmax||l3<eflag_E_lmin||l3>eflag_E_lmax){
		result = 0e0;
	}else{
		i = l1-eflag_T_lmin;
		j = l2-eflag_T_lmin;
		k = l3-eflag_E_lmin;
	
		double b1,b2;
	 	b1 = get_basis_late_T(i,p1)*get_basis_late_T(j,p2)*get_basis_late_E(k,p3);
	 	b2 = get_basis_late_T(i,p2)*get_basis_late_T(j,p1)*get_basis_late_E(k,p3);
	
	 	result = (b1+b2)/(2e0);
 	}
	return result;
}

double plijk_TEE(int r, int l1, int l2, int l3){
	
	int p1,p2,p3;
	find_perm_late_TEE(r,&p1,&p2,&p3);
	double result;
	double b1,b2;
	int i,j,k;
	
	if(l1<eflag_T_lmin||l1>eflag_T_lmax||l2<eflag_E_lmin||l2>eflag_E_lmax||l3<eflag_E_lmin||l3>eflag_E_lmax){
		result = 0e0;
	}else{
		i = l1-eflag_T_lmin;
		j = l2-eflag_E_lmin;
		k = l3-eflag_E_lmin;
	 	b1 = get_basis_late_T(i,p1)*get_basis_late_E(j,p2)*get_basis_late_E(k,p3);
	 	b2 = get_basis_late_T(i,p1)*get_basis_late_E(j,p3)*get_basis_late_E(k,p2);
	
	 	result = (b1+b2)/(2e0);
 	}
	return result;
}

double plijk_EEE(int r, int l1, int l2, int l3){
	
	int p1,p2,p3;
	find_perm_late_EEE(r,&p1,&p2,&p3);
	double result;
	double b1,b2,b3,b4,b5,b6;
	int i,j,k;
	
	if(l1<eflag_E_lmin||l1>eflag_E_lmax||l2<eflag_E_lmin||l2>eflag_E_lmax||l3<eflag_E_lmin||l3>eflag_E_lmax){
		result = 0e0;
	}else{
		i = l1-eflag_E_lmin;
		j = l2-eflag_E_lmin;
		k = l3-eflag_E_lmin;
	
	 	b1 = get_basis_late_E(i,p1)*get_basis_late_E(j,p2)*get_basis_late_E(k,p3);
	 	b2 = get_basis_late_E(i,p2)*get_basis_late_E(j,p3)*get_basis_late_E(k,p1);
	 	b3 = get_basis_late_E(i,p3)*get_basis_late_E(j,p1)*get_basis_late_E(k,p2);
	 	b4 = get_basis_late_E(i,p3)*get_basis_late_E(j,p2)*get_basis_late_E(k,p1);
	 	b5 = get_basis_late_E(i,p2)*get_basis_late_E(j,p1)*get_basis_late_E(k,p3);
	 	b6 = get_basis_late_E(i,p1)*get_basis_late_E(j,p3)*get_basis_late_E(k,p2);

	 	result = (b1+b2+b3+b4+b5+b6)/(6e0);
	}
 	
	return result;
}
 
double rlijk_TTT(int r,int i,int j,int k){

	int n;
	double result = 0;
	
	for(n=0;n<r+1;n++){
		result += get_lambdal_TTT(r,n)*plijk_TTT(n,i,j,k);
	}
 	
	return result;
}
 
double rlijk_TTE(int r,int i,int j,int k){

	int n;
	double result = 0;
	
	for(n=0;n<r+1;n++){
		result += get_lambdal_TTE(r,n)*plijk_TTE(n,i,j,k);
	}
 	
	return result;
}
 
double rlijk_TEE(int r,int i,int j,int k){

	int n;
	double result = 0;
	
	for(n=0;n<r+1;n++){
		result += get_lambdal_TEE(r,n)*plijk_TEE(n,i,j,k);
	}
 	
	return result;
}
 
double rlijk_EEE(int r,int i,int j,int k){

	int n;
	double result = 0;
	
	for(n=0;n<r+1;n++){
		result += get_lambdal_EEE(r,n)*plijk_EEE(n,i,j,k);
	}
 	
	return result;
}

double bijk_TTT(int i,int j,int k){

	int n;
	double result = 0;
	int s=get_terms_late();
	
	for(n=0;n<s;n++){
		result += get_modes_TTT(n)*plijk_TTT(n,i,j,k);
	}
 	
	return result;
}

double bijk_TTE(int i,int j,int k){

	int n;
	double result = 0;
	int s=get_terms_late();
	
	for(n=0;n<s;n++){
		result += get_modes_TTE(n)*plijk_TTE(n,i,j,k);
	}
 	
	return result;
}

double bijk_TEE(int i,int j,int k){

	int n;
	double result = 0;
	int s=get_terms_late();
	
	for(n=0;n<s;n++){
		result += get_modes_TEE(n)*plijk_TEE(n,i,j,k);
	}
 	
	return result;
}

double bijk_EEE(int i,int j,int k){

	int n;
	double result = 0;
	int s=get_terms_late();
	
	for(n=0;n<s;n++){
		result += get_modes_EEE(n)*plijk_EEE(n,i,j,k);
	}
 	
	return result;
}