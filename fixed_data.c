#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include "global.h"
#include "offload_util.h"

static double ***alpha;
static double **basis_prim;
static double **basis_late_T;
static double **basis_late_E;
static double **basis_tri_prim;
static double **basis_tri_late;
static double ***qtilde_T;
static double ***qtilde_E;
static double ***beta_tri;
static double ***bispectrum;
static double ****trispectrum;
static double ***decompose;
static double *eigen;
static double *eigenR;
static double *eigen_tri;
static double *eigenR_tri;
static double *modes_TTT;
static double *modes_TTE;
static double *modes_TEE;
static double *modes_EEE;
static double *modesR_TTT;
static double *modesR_TTE;
static double *modesR_TEE;
static double *modesR_EEE;
static double *modes_tri;
static double *modesR_tri;
static double **ortho;
static double **orthol_TTT;
static double **orthol_TTE;
static double **orthol_TEE;
static double **orthol_EEE;
static double **ortho_tri;
static double **orthol_tri;
static double **lambda;
static double **lambdal_TTT;
static double **lambdal_TTE;
static double **lambdal_TEE;
static double **lambdal_EEE;
static double **lambda_tri;
static double **lambdal_tri;
static double **gamma_TTT;
static double **gamma_TTE;
static double **gamma_TEE;
static double **gamma_EEE;
static double **gamma_tri;
static double *cls_TT;
static double *cls_TE;
static double *cls_EE;
static double *cls_P;
static double *cls_TP;
static double *beam_TT;
static double *beam_TE;
static double *beam_EE;
static double *noise_TT;
static double *noise_TE;
static double *noise_EE;
static double *t_wgt;
static double **transfer;
static double **fisher_tri;

static int lmax;
 static int pmax_prim;
 static int pmax_late_T;
 static int pmax_late_E;
 static int terms_prim;
 static int terms_late;
// static int **order_prim;
 static int *order_prim_flat;
 static int *order_late_TTT_flat;
 static int *order_late_TTE_flat;
 static int *order_late_TEE_flat;
 static int *order_late_EEE_flat;
static int **order_tri_prim;
static int **order_tri_late;
static int qtilde_T_lsize;
static int qtilde_E_lsize;
 static int qtilde_xsize;
static int *qtilde_T_lvec;
static int *qtilde_E_lvec;
 static double *qtilde_xvec;
static int bt_lsize;
static int bt_xsize;
static int *bt_lvec;
static double *bt_xvec;

int init_lmax(int l){
	lmax = l;
	return 0;
}

int get_lmax(){
	return lmax;
}

int set_terms_prim(){
	/*
	When expanding in primordial basis Q_n(k) = SUM_{ijk}	q_i(k)q_j(k)q_k(k),
	there is a 1-to-1 mapping n <----> {ijk}.
	E.g. n = 0, ijk = 000
		n = 1  ijk = 001 
		n = 2  ijk = 010
		...
		up to n = basis:terms,
		which denotes the total number of basis elements Q(k) basis:terms

	This ordering of terms in the sum is defined by this function, 
	and is determined by the parameter basis:order_prim in parameter.ini file.
	*/
	int r,i,j,k,n,s;
	int* order_raw = malloc(sizeof(int)*3*MAXLINES);
	int* size = malloc(sizeof(int));
	
	terms_prim = alpha_max;	//maximum index number of basis functions
	pmax_prim = 0;
	create_order_prim();	//necessary to calculate order_prim_flat
	int (* restrict order_prim)[3] = (int (*restrict)[3]) order_prim_flat;

	//Different orderings define different mappings, which is set here using the switch		
	switch(eflag_order_prim){
		case 1:
			i=j=k=s=0;
			for (n=0;n<terms_prim;n++){
				order_prim[n][0] = i;
				order_prim[n][1] = j;
				order_prim[n][2] = k;
				if(k-j<2&&j-i<2&&k-i<2){
					s++;
					i=j=0;
					k=s;
				} else if(k-j>1){
					k--;
					j++;
				} else {
					i++;
					j=i;
					k=s-2*i;
				}
				if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
			}
		break;
		case 2:
			i=j=k=s=0;
			for (n=0;n<terms_prim;n++){
				order_prim[n][0] = i;
				order_prim[n][1] = j;
				order_prim[n][2] = k;
				if(i<j){
					j--;
					k++;
				} else if(i>0){
					i--;
					j = (s-i)/2;
					k = s-j-i;
				} else {
					s++;
					i = s/3;
					j = (s-i)/2;
					k = s-i-j;
				}
				if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
			}
		break;
		case 3:
			load_three_int(edist_file, order_raw, size);
			if(terms_prim>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_prim);
				exit;
			} else {
				for (n=0;n<terms_prim;n++){
					order_prim[n][0] = order_raw[3*n];
					order_prim[n][1] = order_raw[3*n+1];
					order_prim[n][2] = order_raw[3*n+2];
					if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
				}
			}
		break;
		case 4:
			terms_prim = alpha_max;
			pmax_prim = 0;
			//create_order_prim();
			load_three_int(edist_file, order_raw, size);
			if(terms_prim>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_prim);
				exit;
			} else {
				for (n=0;n<terms_prim-1;n++){
					order_prim[n][0] = order_raw[3*n];
					order_prim[n][1] = order_raw[3*n+1];
					order_prim[n][2] = order_raw[3*n+2];
					if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
				}
			}
		
			for (n=terms_prim-1;n>1;n--){
				order_prim[n][0] = order_prim[n-1][0];
				order_prim[n][1] = order_prim[n-1][1];
				order_prim[n][2] = order_prim[n-1][2];
			}
			
			order_prim[1][0] = pmax_prim+1;
			order_prim[1][1] = pmax_prim+1;
			order_prim[1][2] = pmax_prim+2;
		break;
		case 5:
			load_three_int(edist_file, order_raw, size);
			if(terms_prim>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_prim);
				exit;
			} else {
				for (n=0;n<terms_prim;n++){
					order_prim[n][0] = order_raw[3*n];
					order_prim[n][1] = order_raw[3*n+1];
					order_prim[n][2] = order_raw[3*n+2];
					if (order_prim[n][2]>pmax_prim)pmax_prim = order_prim[n][2];
				}
			}
		break;
		default:
			printf("invalid mode specified for eigen function ordering\n");
			exit;
		break;
	}

	pmax_prim = 0;
	for (n=0;n<terms_prim;n++){
		if(order_prim[n][2]>pmax_prim) pmax_prim = order_prim[n][2];
	}
	
	return 0;
}

int get_terms_prim(){
//	gives the number of terms present in the expansion of primordial basis
	return terms_prim;
}


int get_pmax_prim(){
// gives the maximum index of the {ijk} in the primordial basis expansion of Q_n (k)
	return pmax_prim;
}

int create_order_prim(){
/*
	need three times as large memory, because we are dealing with 
	a product of three q's:
	q_i * q_j * q_k
*/
	order_prim_flat = (int *)malloc(terms_prim*3*sizeof(int));
	return 0;
}


void find_perm_prim(int r,int* p1,int* p2,int* p3){
//	For a given index r, return the mapped indices {p1,p2,p3}
    int (* restrict order_prim)[3] = (int (*restrict)[3]) order_prim_flat;
	*p1 = order_prim[r][0];
	*p2 = order_prim[r][1];
	*p3 = order_prim[r][2];
	return;
}

int set_terms_late(){
/*
When expanding in late-time basis Q_n(l) = SUM_{ijk}	q_i(l)q_j(l)q_k(l),
there is a 1-to-1 mapping n <----> {ijk}.
E.g. n = 0, ijk = 000
	 n = 1  ijk = 001 
	 n = 2  ijk = 010
	 ...
	 up to n = alphamax,
	  which denotes the total number of basis elements Q(k)
	  alphamax = basis:terms

This ordering of terms in the sum is defined by this function, 
and is determined by the parameter basis:order_late in parameter.ini file.

This procedure is more algebraically tedious than primordial calculation, 
because polarisation must be considered
*/	
	int r,i,j,k,n,s;
	int* order_raw = malloc(sizeof(int)*3*MAXLINES);
	int* size = malloc(sizeof(int));
	
	terms_late = alpha_max; //total number of terms in the (truncated) expansion
	pmax_late_T = 0;
	create_order_late();	//prepare order_late_XXX_flat amounts of memory (X= T or E)

	int (* restrict order_late_TTT)[3] = (int (*restrict)[3]) order_late_TTT_flat;
	int (* restrict order_late_TTE)[3] = (int (*restrict)[3]) order_late_TTE_flat;
	int (* restrict order_late_TEE)[3] = (int (*restrict)[3]) order_late_TEE_flat;
	int (* restrict order_late_EEE)[3] = (int (*restrict)[3]) order_late_EEE_flat;

//	Create the mapping n <---> {ijk} based on the ordering defined 
	switch(eflag_order_late){
		case 1:
			i=j=k=s=0;
			for (n=0;n<terms_late;n++){
				order_late_TTT[n][0] = i;
				order_late_TTT[n][1] = j;
				order_late_TTT[n][2] = k;
				if(k-j<2&&j-i<2&&k-i<2){
					s++;
					i=j=0;
					k=s;
				} else if(k-j>1){
					k--;
					j++;
				} else {
					i++;
					j=i;
					k=s-2*i;
				
				}
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 2:
			i=j=k=s=0;
			for (n=0;n<terms_late;n++){
				order_late_TTT[n][0] = i;
				order_late_TTT[n][1] = j;
				order_late_TTT[n][2] = k;
				if(i<j){
					j--;
					k++;
				} else if(i>0){
					i--;
					j = (s-i)/2;
					k = s-j-i;
				} else {
					s++;
					i = s/3;
					j = (s-i)/2;
					k = s-i-j;
				}
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 3:
			load_three_int(mdist_TTT_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in TTT distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late;n++){
					order_late_TTT[n][0] = order_raw[3*n];
					order_late_TTT[n][1] = order_raw[3*n+1];
					order_late_TTT[n][2] = order_raw[3*n+2];
				}
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 4:
			load_three_int(mdist_TTT_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in TTT distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-1;n++){
					order_late_TTT[n][0] = order_raw[3*n];
					order_late_TTT[n][1] = order_raw[3*n+1];
					order_late_TTT[n][2] = order_raw[3*n+2];
				}
			}
			
			for (n=terms_late-1;n>1;n--){
				order_late_TTT[n][0] = order_late_TTT[n-1][0];
				order_late_TTT[n][1] = order_late_TTT[n-1][1];
				order_late_TTT[n][2] = order_late_TTT[n-1][2];
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 5:
			load_three_int(mdist_TTT_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in TTT distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-3;n++){
					order_late_TTT[n][0] = order_raw[3*n];
					order_late_TTT[n][1] = order_raw[3*n+1];
					order_late_TTT[n][2] = order_raw[3*n+2];
				}
			}
			
			for (n=terms_late-1;n>3;n--){
				order_late_TTT[n][0] = order_late_TTT[n-3][0];
				order_late_TTT[n][1] = order_late_TTT[n-3][1];
				order_late_TTT[n][2] = order_late_TTT[n-3][2];
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 6:
			load_three_int(mdist_TTT_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in TTT distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-4;n++){
					order_late_TTT[n][0] = order_raw[3*n];
					order_late_TTT[n][1] = order_raw[3*n+1];
					order_late_TTT[n][2] = order_raw[3*n+2];
				}
			}
			
			for (n=terms_late-1;n>4;n--){
				order_late_TTT[n][0] = order_late_TTT[n-4][0];
				order_late_TTT[n][1] = order_late_TTT[n-4][1];
				order_late_TTT[n][2] = order_late_TTT[n-4][2];
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 7:
			load_three_int(mdist_TTT_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in TTT distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-1;n++){
					order_late_TTT[n][0] = order_raw[3*n];
					order_late_TTT[n][1] = order_raw[3*n+1];
					order_late_TTT[n][2] = order_raw[3*n+2];
				}
			}
			
			for (n=terms_late-1;n>1;n--){
				order_late_TTT[n][0] = order_late_TTT[n-1][0];
				order_late_TTT[n][1] = order_late_TTT[n-1][1];
				order_late_TTT[n][2] = order_late_TTT[n-1][2];
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 8:
			load_three_int(mdist_TTT_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in TTT distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late;n++){
					order_late_TTT[n][0] = order_raw[3*n];
					order_late_TTT[n][1] = order_raw[3*n+1];
					order_late_TTT[n][2] = order_raw[3*n+2];
				}
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 9:
			load_three_int(mdist_TTT_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in TTT distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-1;n++){
					order_late_TTT[n][0] = order_raw[3*n];
					order_late_TTT[n][1] = order_raw[3*n+1];
					order_late_TTT[n][2] = order_raw[3*n+2];
				}
			}
			
			for (n=terms_late-1;n>1;n--){
				order_late_TTT[n][0] = order_late_TTT[n-1][0];
				order_late_TTT[n][1] = order_late_TTT[n-1][1];
				order_late_TTT[n][2] = order_late_TTT[n-1][2];
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		case 10:
			load_three_int(mdist_TTT_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in TTT distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late-1;n++){
					order_late_TTT[n][0] = order_raw[3*n];
					order_late_TTT[n][1] = order_raw[3*n+1];
					order_late_TTT[n][2] = order_raw[3*n+2];
				}
			}
			
			for (n=terms_late-1;n>=3;n--){
				order_late_TTT[n][0] = order_late_TTT[n-3][0];
				order_late_TTT[n][1] = order_late_TTT[n-3][1];
				order_late_TTT[n][2] = order_late_TTT[n-3][2];
			}
			
			if(do_polarisation==1){
				load_three_int(mdist_TTE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TTE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TTE[n][0] = order_raw[3*n];
						order_late_TTE[n][1] = order_raw[3*n+1];
						order_late_TTE[n][2] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_TEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in TEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_TEE[n][1] = order_raw[3*n];
						order_late_TEE[n][2] = order_raw[3*n+1];
						order_late_TEE[n][0] = order_raw[3*n+2];
					}
				}
				load_three_int(mdist_EEE_file, order_raw, size);
				if(terms_late>*size){
					printf("eigenfunction %d exceeds that in EEE distance file\n",terms_late);
					exit;
				} else {
					for (n=0;n<terms_late;n++){
						order_late_EEE[n][0] = order_raw[3*n];
						order_late_EEE[n][1] = order_raw[3*n+1];
						order_late_EEE[n][2] = order_raw[3*n+2];
					}
				}
			}
		break;
		default:
		printf("invalid mode specified for eigen function ordering\n");
			exit;
		break;
	}
	
	pmax_late_T = 0;
	for (n=0;n<terms_late;n++){
		if(order_late_TTT[n][2]>pmax_late_T) pmax_late_T = order_late_TTT[n][2];
	}
	if(do_polarisation==1){
		pmax_late_E = 0;
		for (n=0;n<terms_late;n++){
			if(order_late_TTE[n][1]>pmax_late_T) pmax_late_T = order_late_TTE[n][1];
			if(order_late_TEE[n][0]>pmax_late_T) pmax_late_T = order_late_TEE[n][0];
			if(order_late_TTE[n][2]>pmax_late_E) pmax_late_E = order_late_TTE[n][2];
			if(order_late_TEE[n][2]>pmax_late_E) pmax_late_E = order_late_TEE[n][2];
			if(order_late_EEE[n][2]>pmax_late_E) pmax_late_E = order_late_EEE[n][2];
		}

	}

	switch(eflag_order_late){
		case 4:
			order_late_TTT[1][0] = pmax_late_T+1;
			order_late_TTT[1][1] = pmax_late_T+1;
			order_late_TTT[1][2] = pmax_late_T+2;
			
			pmax_late_T = pmax_late_T+2;
		break;
		case 5:
			order_late_TTT[1][0] = pmax_late_T+1;
			order_late_TTT[1][1] = pmax_late_T+1;
			order_late_TTT[1][2] = pmax_late_T+1;
			
			order_late_TTT[2][0] = pmax_late_T+2;
			order_late_TTT[2][1] = pmax_late_T+2;
			order_late_TTT[2][2] = pmax_late_T+3;
			
			order_late_TTT[3][0] = pmax_late_T+4;
			order_late_TTT[3][1] = pmax_late_T+4;
			order_late_TTT[3][2] = pmax_late_T+5;
			
			pmax_late_T = pmax_late_T+5;
		break;
		case 6:
			order_late_TTT[1][0] = pmax_late_T+1;
			order_late_TTT[1][1] = pmax_late_T+1;
			order_late_TTT[1][2] = pmax_late_T+2;
			
			order_late_TTT[2][0] = pmax_late_T+3;
			order_late_TTT[2][1] = pmax_late_T+5;
			order_late_TTT[2][2] = pmax_late_T+8;
			
			order_late_TTT[3][0] = pmax_late_T+3;
			order_late_TTT[3][1] = pmax_late_T+6;
			order_late_TTT[3][2] = pmax_late_T+7;
			
			order_late_TTT[4][0] = pmax_late_T+4;
			order_late_TTT[4][1] = pmax_late_T+5;
			order_late_TTT[4][2] = pmax_late_T+6;
			
			pmax_late_T = pmax_late_T+8;
		break;
		case 7:
			order_late_TTT[1][0] = pmax_late_T+1;
			order_late_TTT[1][1] = pmax_late_T+1;
			order_late_TTT[1][2] = pmax_late_T+2;
			
			pmax_late_T = pmax_late_T+2;
		break;
		case 9:
			order_late_TTT[1][0] = pmax_late_T+1;
			order_late_TTT[1][1] = pmax_late_T+1;
			order_late_TTT[1][2] = pmax_late_T+2;
			
			pmax_late_T = pmax_late_T+2;
		break;
		case 10:
			order_late_TTT[0][0] = pmax_late_T+2;
			order_late_TTT[0][1] = pmax_late_T+2;
			order_late_TTT[0][2] = pmax_late_T+4;
			order_late_TTT[1][0] = pmax_late_T+1;
			order_late_TTT[1][1] = pmax_late_T+3;
			order_late_TTT[1][2] = pmax_late_T+4;
			order_late_TTT[2][0] = pmax_late_T+1;
			order_late_TTT[2][1] = pmax_late_T+1;
			order_late_TTT[2][2] = pmax_late_T+5;
			
			pmax_late_T = pmax_late_T+5;
		break;
	}

	/*#pragma offload_transfer target(mic:offload_target) \
	    in(pmax_late_T) \
	    in(order_late_TTT_flat[0:terms_late*3] : ALLOC RETAIN)
	if (do_polarisation)
    {
        #pragma offload_transfer target(mic:offload_target) \
            in(pmax_late_E) \
            in(order_late_TTE_flat[0:terms_late*3] : ALLOC RETAIN) \
            in(order_late_TEE_flat[0:terms_late*3] : ALLOC RETAIN) \
            in(order_late_EEE_flat[0:terms_late*3] : ALLOC RETAIN)
    }*/
	
	return 0;
}
 
int get_terms_late(){
	return terms_late;
}
 
int get_pmax_late_T(){
	return pmax_late_T;
}
 
int get_pmax_late_E(){
	return pmax_late_E;
}

int create_order_late(){
/*
	Allocate memory for three q(l) polynomials in the product,
	and additionally if polarisation is considered	
*/

	order_late_TTT_flat = (int *) malloc(terms_late*3 * sizeof(int));
	if(do_polarisation==1){
        order_late_TTE_flat = (int *) malloc(terms_late*3 * sizeof(int));
        order_late_TEE_flat = (int *) malloc(terms_late*3 * sizeof(int));
        order_late_EEE_flat = (int *) malloc(terms_late*3 * sizeof(int));
	}
	return 0;
}


void find_perm_late_TTT(int r,int* p1,int* p2,int* p3){
	
	int (* restrict order_late_TTT)[3] = (int (*restrict)[3]) order_late_TTT_flat;
	*p1 = order_late_TTT[r][0];
	*p2 = order_late_TTT[r][1];
	*p3 = order_late_TTT[r][2];
	return;
}
 
void find_perm_late_TTE(int r,int* p1,int* p2,int* p3){
	
	int (* restrict order_late_TTE)[3] = (int (*restrict)[3]) order_late_TTE_flat;
	*p1 = order_late_TTE[r][0];
	*p2 = order_late_TTE[r][1];
	*p3 = order_late_TTE[r][2];
	return;
}
 
void find_perm_late_TEE(int r,int* p1,int* p2,int* p3){
	
	int (* restrict order_late_TEE)[3] = (int (*restrict)[3]) order_late_TEE_flat;
	*p1 = order_late_TEE[r][0];
	*p2 = order_late_TEE[r][1];
	*p3 = order_late_TEE[r][2];
	return;
}
 
void find_perm_late_EEE(int r,int* p1,int* p2,int* p3){
	
	int (* restrict order_late_EEE)[3] = (int (*restrict)[3]) order_late_EEE_flat;
	*p1 = order_late_EEE[r][0];
	*p2 = order_late_EEE[r][1];
	*p3 = order_late_EEE[r][2];
	return;
}

int set_terms_tri_prim(){
	
	int r,i,j,k,n,s;
	int* order_raw = malloc(sizeof(int)*4*MAXLINES);
	int* size = malloc(sizeof(int));
	
	terms_prim = alpha_max;
	pmax_late_T = 0;
	create_order_tri_prim();
	
	switch(eflag_order_prim){
		case 3:
			load_four_int(edist_tri_file, order_raw, size);
			if(terms_prim>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_prim);
				exit;
			} else {
				for (n=0;n<terms_prim;n++){
					order_tri_prim[n][0] = order_raw[4*n];
					order_tri_prim[n][1] = order_raw[4*n+1];
					order_tri_prim[n][2] = order_raw[4*n+2];
					order_tri_prim[n][3] = order_raw[4*n+3];
					if (order_tri_prim[n][3]>pmax_prim)pmax_prim = order_tri_prim[n][3];
				}
			}
		break;
		default:
			printf("invalid mode specified for eigen function ordering\n");
			exit;
		break;
	}
	return 0;
}

int create_order_tri_prim(){
	order_tri_prim = (int **)create_iarray(terms_prim,4);
	return 0;
}

void find_perm_tri_prim(int r,int* p1,int* p2,int* p3,int* p4){
	
	*p1 = order_tri_prim[r][0];
	*p2 = order_tri_prim[r][1];
	*p3 = order_tri_prim[r][2];
	*p4 = order_tri_prim[r][3];
	return;
}

int set_terms_tri_late(){
	
	int r,i,j,k,n,s;
	int* order_raw = malloc(sizeof(int)*4*MAXLINES);
	int* size = malloc(sizeof(int));
	
	terms_late = alpha_max;
	pmax_late_T = 0;
	create_order_tri_late();
	
	switch(eflag_order_late){
		case 3:
			load_four_int(edist_tri_file, order_raw, size);
			if(terms_late>*size){
				printf("eigenfunction %d exceeds that in distance file\n",terms_late);
				exit;
			} else {
				for (n=0;n<terms_late;n++){
					order_tri_late[n][0] = order_raw[4*n];
					order_tri_late[n][1] = order_raw[4*n+1];
					order_tri_late[n][2] = order_raw[4*n+2];
					order_tri_late[n][3] = order_raw[4*n+3];
					if (order_tri_late[n][3]>pmax_late_T)pmax_late_T = order_tri_late[n][3];
				}
			}
		break;
		default:
			printf("invalid mode specified for eigen function ordering\n");
			exit;
		break;
	}
	return 0;
}


int create_order_tri_late(){
	order_tri_late = (int **)create_iarray(terms_late,4);
	return 0;
}

void find_perm_tri_late(int r,int* p1,int* p2,int* p3,int* p4){
	
	*p1 = order_tri_late[r][0];
	*p2 = order_tri_late[r][1];
	*p3 = order_tri_late[r][2];
	*p4 = order_tri_late[r][3];
	return;
}

int load_cl(){
	
	cls_TT = create_vector(eflag_T_lmax+1);
	if(do_polarisation==1){
		cls_TE = create_vector(eflag_E_lmax+1);
		cls_EE = create_vector(eflag_E_lmax+1);
	}
	
	if(clflag_uselens==1){
		double* cls_data = malloc( sizeof(double)*MAXLINES*5);
		int* cl_len = malloc( sizeof(int));
		load_txt_dbl(cls_lensed_data_file, 5, cls_data, cl_len);
		
		int cl_size = *cl_len;
		
		double l_raw[cl_size];
		double cl_TT_raw[cl_size];
		double cl_TE_raw[cl_size];
		double cl_BB_raw[cl_size];
		double cl_EE_raw[cl_size];
		double pt,clpt1,clpt2;
		
		int i,j;
		
		j=0;
		
		for (i=0; i<cl_size; i++){
			l_raw[i] = cls_data[j++];
			cl_TT_raw[i] = cls_data[j++];
			cl_EE_raw[i] = cls_data[j++];
			cl_BB_raw[i] = cls_data[j++];
			cl_TE_raw[i] = cls_data[j++];
		}
		
		if (eflag_T_lmax>l_raw[cl_size-1]){
			printf("Cls do not contain enough l's, max data %d max cl: %d\n", eflag_T_lmax, (int)l_raw[cl_size-1]);
			return 1;
			exit;
		}
		
		gsl_spline* spcltt =  gsl_spline_alloc(gsl_interp_cspline, cl_size);
		gsl_interp_accel* acccltt = gsl_interp_accel_alloc();
		
		gsl_spline_init(spcltt,l_raw,cl_TT_raw,cl_size);
		
		for (i=0; i<eflag_T_lmax+1; i++){
			pt = (double)i;
			if(pt>l_raw[0]-5e-1){
				clpt1 = gsl_spline_eval(spcltt,pt,acccltt);
				cls_TT[i] = 2e0*M_PI*clpt1/(pt*(pt+1e0));
			}else{
				cls_TT[i] = 0e0;
			}
		}
		
		gsl_spline_free(spcltt);
		gsl_interp_accel_free(acccltt);
		
		if(do_polarisation==1){
			
			gsl_spline* spclte =  gsl_spline_alloc(gsl_interp_cspline, cl_size);
			gsl_spline* spclee =  gsl_spline_alloc(gsl_interp_cspline, cl_size);
			gsl_interp_accel* accclte = gsl_interp_accel_alloc();
			gsl_interp_accel* accclee = gsl_interp_accel_alloc();
			gsl_spline_init(spclte,l_raw,cl_TE_raw,cl_size);
			gsl_spline_init(spclee,l_raw,cl_EE_raw,cl_size);
			
			for (i=0; i<eflag_E_lmax+1; i++){
				pt = (double)i;
				if(pt>l_raw[0]-5e-1){
					clpt1 = gsl_spline_eval(spclte,pt,accclte);
					clpt2 = gsl_spline_eval(spclee,pt,accclee);
					cls_TE[i] = 2e0*M_PI*clpt1/(pt*(pt+1e0));
					cls_EE[i] = 2e0*M_PI*clpt2/(pt*(pt+1e0));
				}else{
					cls_TE[i] = 0e0;
					cls_EE[i] = 0e0;
				}
			}
			
			gsl_spline_free(spclte);
			gsl_spline_free(spclee);
			gsl_interp_accel_free(accclte);
			gsl_interp_accel_free(accclee);
		}
	}else{
		double* cls_data = malloc( sizeof(double)*MAXLINES*6);
		int* cl_len = malloc( sizeof(int));
		load_txt_dbl(cls_scalar_data_file, 6, cls_data, cl_len);
		
		int cl_size = *cl_len;
		
		double l_raw[cl_size];
		double cl_TT_raw[cl_size];
		double cl_TE_raw[cl_size];
		double cl_EE_raw[cl_size];
		double cl_P_raw[cl_size];
		double cl_PT_raw[cl_size];
		double pt,clpt1,clpt2;
		
		int i,j;
		
		j=0;
		
		for (i=0; i<cl_size; i++){
			l_raw[i] = cls_data[j++];
			cl_TT_raw[i] = cls_data[j++];
			cl_EE_raw[i] = cls_data[j++];
			cl_TE_raw[i] = cls_data[j++];
			cl_P_raw[i] = cls_data[j++];
			cl_PT_raw[i] = cls_data[j++];
		}
		
		if (eflag_T_lmax>l_raw[cl_size-1]){
			printf("Cls do not contain enough l's, max data %d max cl: %d\n", eflag_T_lmax, (int)l_raw[cl_size-1]);
			return 1;
			exit;
		}
		
		gsl_spline* spcltt =  gsl_spline_alloc(gsl_interp_cspline, cl_size);
		gsl_interp_accel* acccltt = gsl_interp_accel_alloc();
		
		gsl_spline_init(spcltt,l_raw,cl_TT_raw,cl_size);
		
		for (i=0; i<eflag_T_lmax+1; i++){
			pt = (double)i;
			if(pt>l_raw[0]-5e-1){
				clpt1 = gsl_spline_eval(spcltt,pt,acccltt);
				cls_TT[i] = 2e0*M_PI*clpt1/(pt*(pt+1e0));
			}else{
				cls_TT[i] = 0e0;
			}
		}
			
		gsl_spline_free(spcltt);
		gsl_interp_accel_free(acccltt);
		
		if(do_polarisation==1){
			
			gsl_spline* spclte =  gsl_spline_alloc(gsl_interp_cspline, cl_size);
			gsl_spline* spclee =  gsl_spline_alloc(gsl_interp_cspline, cl_size);
			gsl_interp_accel* accclte = gsl_interp_accel_alloc();
			gsl_interp_accel* accclee = gsl_interp_accel_alloc();
			gsl_spline_init(spclte,l_raw,cl_TE_raw,cl_size);
			gsl_spline_init(spclee,l_raw,cl_EE_raw,cl_size);
			
			for (i=0; i<eflag_E_lmax+1; i++){
				pt = (double)i;
				if(pt>l_raw[0]-5e-1){
					clpt1 = gsl_spline_eval(spclte,pt,accclte);
					clpt2 = gsl_spline_eval(spclee,pt,accclee);
					cls_TE[i] = 2e0*M_PI*clpt1/(pt*(pt+1e0));
					cls_EE[i] = 2e0*M_PI*clpt2/(pt*(pt+1e0));
				}else{
					cls_TE[i] = 0e0;
					cls_EE[i] = 0e0;
				}
			}
			
			gsl_spline_free(spclte);
			gsl_spline_free(spclee);
			gsl_interp_accel_free(accclte);
			gsl_interp_accel_free(accclee);
		}
	}
	
	return 0;
}

double get_cl_TT(int i){
	double value;
	value = cls_TT[i];
	return value;
}

double get_cl_TE(int i){
	double value;
	value = cls_TE[i];
	return value;
}

double get_cl_EE(int i){
	double value;
	value = cls_EE[i];
	return value;
}

int load_BN(){
	double* pixw_data = malloc(sizeof(double)*MAXLINES*2);
	int* pixw_len = malloc(sizeof(int));
	load_txt_dbl(pixw_file, 2, pixw_data, pixw_len);
	
	int pixw_size = *pixw_len;
	double pixw_raw[pixw_size];
	int i,j,l;
	double pt;
	j=0;
	for (i=0; i<pixw_size; i++){
		l = pixw_data[j++];
		pixw_raw[i] = pixw_data[j++];
	}

	beam_TT = create_vector(eflag_T_lmax+1);
	noise_TT = create_vector(eflag_T_lmax+1);
	
	double* BN_data = malloc( sizeof(double)*MAXLINES*3);
	int* BN_len = malloc( sizeof(int));
	load_txt_dbl(BN_TT_file, 3, BN_data, BN_len);
	
	int BN_size = *BN_len;
	
	double l_raw[BN_size];
	double beam_raw[BN_size];
	double noise_raw[BN_size];
	
	j=0;
	for (i=0; i<BN_size; i++){
		l_raw[i] = BN_data[j++];
		beam_raw[i] = BN_data[j++];
		noise_raw[i] = BN_data[j++];
	}
	if (eflag_T_lmax>l_raw[BN_size-1]){
		printf("BN do not contain enough l's, max data %d max BN: %d\n", eflag_T_lmax, (int)l_raw[BN_size]);
		return 1;
		exit;
	}
	
	gsl_spline* spB =  gsl_spline_alloc (gsl_interp_cspline, BN_size);
	gsl_spline* spN =  gsl_spline_alloc (gsl_interp_cspline, BN_size);
	gsl_interp_accel* accB = gsl_interp_accel_alloc();
	gsl_interp_accel* accN = gsl_interp_accel_alloc();
	
	gsl_spline_init(spB,l_raw,beam_raw,BN_size);
	gsl_spline_init(spN,l_raw,noise_raw,BN_size);
	
	for (i=0; i<eflag_T_lmax+1; i++){
		pt = (double)i;
		beam_TT[i] = 1e0;
		if(pt>l_raw[0]-5e-1){
			beam_TT[i] = pixw_raw[i]*gsl_spline_eval(spB,pt,accB);
		}else{
			beam_TT[i] = 0e0;
		}
		//printf("PIXW_RAW[%d] = %e\t", i, pixw_raw[i]);
		//printf("GSL SPLINE EVAL = %e\t", gsl_spline_eval(spB,pt,accB));
	}
	
	for (i=0; i<eflag_T_lmax+1; i++){
		pt = (double)i;
		noise_TT[i] = 0e0;
		if(pt>l_raw[0]-5e-1){
			noise_TT[i] = gsl_spline_eval(spN,pt,accN);
		}else{
			noise_TT[i] = 0e0;
		}
	}
	
	// 	free(l_raw);
	// 	free(beam_raw);
	// 	free(noise_raw);
	
	gsl_spline_free(spB);
	gsl_spline_free(spN);
	gsl_interp_accel_free(accB);
	gsl_interp_accel_free(accN);
	
	if(do_polarisation==1){
		
		beam_EE = (double *)create_vector(eflag_E_lmax+1);
		noise_EE = (double *)create_vector(eflag_E_lmax+1);
		
		double* BN_EE_data = malloc( sizeof(double)*MAXLINES*3);
		int* BN_EE_len = malloc( sizeof(int));
		load_txt_dbl(BN_EE_file, 3, BN_EE_data, BN_EE_len);
		
		int BN_EE_size = *BN_EE_len;
		
		double l_EE_raw[BN_EE_size];
		double beam_EE_raw[BN_EE_size];
		double noise_EE_raw[BN_EE_size];
		
		j=0;
		for (i=0; i<BN_EE_size; i++){
			l_EE_raw[i] = BN_EE_data[j++];
			beam_EE_raw[i] = BN_EE_data[j++];
			noise_EE_raw[i] = BN_EE_data[j++];
		}
		if (eflag_E_lmax>l_EE_raw[BN_size-1]){
			printf("BN EE do not contain enough l's, max data %d max BN: %d\n", eflag_E_lmax, (int)l_EE_raw[BN_size]);
			return 1;
			exit;
		}
		
		gsl_spline* spB_EE =  gsl_spline_alloc (gsl_interp_cspline, BN_EE_size);
		gsl_spline* spN_EE =  gsl_spline_alloc (gsl_interp_cspline, BN_EE_size);
		gsl_interp_accel* accB_EE = gsl_interp_accel_alloc();
		gsl_interp_accel* accN_EE = gsl_interp_accel_alloc();
		
		gsl_spline_init(spB_EE,l_EE_raw,beam_EE_raw,BN_EE_size);
		gsl_spline_init(spN_EE,l_EE_raw,noise_EE_raw,BN_EE_size);
		
		for (i=0; i<eflag_E_lmax+1; i++){
			pt = (double)i;
			beam_EE[i] = 1e0;
			if(pt>l_EE_raw[0]-5e-1){
				beam_EE[i] = pixw_raw[i]*gsl_spline_eval(spB_EE,pt,accB_EE);
			}else{
				beam_EE[i] = 0e0;
			}
		}
		
		for (i=0; i<eflag_E_lmax+1; i++){
			pt = (double)i;
			noise_EE[i] = 0e0;
			if(pt>l_EE_raw[0]-5e-1){
				noise_EE[i] = gsl_spline_eval(spN_EE,pt,accN_EE);
			}else{
				noise_EE[i] = 0e0;
			}
		}
		
		gsl_spline_free(spB_EE);
		gsl_spline_free(spN_EE);
		gsl_interp_accel_free(accB_EE);
		gsl_interp_accel_free(accN_EE);
		
		beam_TE = (double *)create_vector(eflag_E_lmax+1);
		noise_TE = (double *)create_vector(eflag_E_lmax+1);
		
		for (i=0; i<eflag_E_lmax+1; i++){
			beam_TE[i] = sqrt(beam_TT[i]*beam_EE[i]);
			noise_TE[i] = 0.0;
		}
	}
	return 0;
}

int load_TL(int l_size ,int* l_values){
	
	t_wgt = (double *)create_vector(eflag_T_lmax+1);
	
	double* T_data = malloc( sizeof(double)*MAXLINES*2);
	int* T_len = malloc( sizeof(int));
	load_txt_dbl(transfer_wgt_file, 2, T_data, T_len);
	
	int T_size = *T_len;
	
	double l_raw[T_size];
	double Tl_raw[T_size];
	
	int i,j;
	double pt;
	
	j=0;
	for (i=0; i<T_size; i++){
		l_raw[i] = T_data[j++];
		Tl_raw[i] = T_data[j++];
// 		printf("%d\t%d\t%d\t%e\n", i,j,(int)l_raw[i],Tl_raw[i]);
	}
	if (l_values[l_size-1]>l_raw[T_size-1]){
		printf("T do not contain enough l's, max data %d max T: %d\n", l_values[l_size-1], (int)l_raw[T_size]);
		return 1;
		exit;
	}
	
	gsl_spline* spT =  gsl_spline_alloc (gsl_interp_cspline, T_size);
	gsl_interp_accel* accT = gsl_interp_accel_alloc();
	
	gsl_spline_init(spT,l_raw,Tl_raw,T_size);
	
	for (i=0; i<l_size; i++){
		pt = (double)l_values[i];
		t_wgt[i] = 0.0;
		if(pt!=0){
			t_wgt[i] = gsl_spline_eval(spT,pt,accT);
		}
// 		printf("%d\t%d\t%e\n", i,l_values[i],t_wgt[i]);
	}
	
	gsl_spline_free(spT);
	gsl_interp_accel_free(accT);
	
	return 0;
}


double get_beam_TT(int i){
	double value = beam_TT[i];
	return value;
}

double get_beam_TE(int i){
	double value = beam_TE[i];
	return value;
}

double get_beam_EE(int i){
	double value = beam_EE[i];
	return value;
}

double get_noise_TT(int i){
	double value = noise_TT[i];
	return value;
}

double get_noise_TE(int i){
	double value = noise_TE[i];
	return value;
}

double get_noise_EE(int i){
	double value = noise_EE[i];
	return value;
}

int create_t_wgt(int l_size){
	t_wgt = (double *)create_vector(l_size);
	return 0;
}

double get_t_wgt(int i){
	double value = t_wgt[i];
	return value;
}

int load_lens(){
	
	cls_TP = (double *)create_vector(eflag_T_lmax+1);
	
	double* cls_data = malloc( sizeof(double)*MAXLINES*6);
	int* cl_len = malloc( sizeof(int));
	load_txt_dbl(cls_scalar_data_file, 6, cls_data, cl_len);
	
	int cl_size = *cl_len;
	
	double l_raw[cl_size];
	double cl_TT_raw[cl_size];
	double cl_TE_raw[cl_size];
	double cl_EE_raw[cl_size];
	double cl_P_raw[cl_size];
	double cl_TP_raw[cl_size];
	double pt,clpt;
	
	int i,j;
	
	j=0;
	
	for (i=0; i<cl_size; i++){
		l_raw[i] = cls_data[j++];
		cl_TT_raw[i] = cls_data[j++];
		cl_EE_raw[i] = cls_data[j++];
		cl_TE_raw[i] = cls_data[j++];
		cl_P_raw[i] = cls_data[j++];
		cl_TP_raw[i] = cls_data[j++];
	}
	
	if (eflag_T_lmax>l_raw[cl_size-1]){
		printf("Cls do not contain enough l's, max data %d max cl: %d\n", eflag_T_lmax, (int)l_raw[cl_size-1]);
		return 1;
		exit;
	}
	
	gsl_spline* spcltp =  gsl_spline_alloc(gsl_interp_cspline, cl_size);
	gsl_interp_accel* acccltp = gsl_interp_accel_alloc();
	
	gsl_spline_init(spcltp,l_raw,cl_TP_raw,cl_size);
	
	for (i=0; i<eflag_T_lmax+1; i++){
		pt = (double)i;
		if(pt>l_raw[0]-5e-1){
			clpt = gsl_spline_eval(spcltp,pt,acccltp);
			cls_TP[i] = clpt/(pt*pt*pt);
		}else{
			cls_TP[i] = 0e0;
		}
	}
	
	gsl_spline_free(spcltp);
	gsl_interp_accel_free(acccltp);
	
	return 0;
}

double get_cl_TP(int i){
	double value = cls_TP[i];
	return value;
}
	
int create_alpha(){

	alpha = (double ***)create_3Darray(alpha_max+1,alpha_max+1,alpha_max+1);
	
	return 0;
}

double get_alpha(int i, int j, int k){
	double value = alpha[i][j][k];
	return value;
}
	
int update_alpha(double *results){

	int i,j,k;
	double value;
	
	i = (int)results[0];
	j = (int)results[1];
	k = (int)results[2];
	value = results[3];
	
	//printf("update: %d\t%d\t%d\t%e\n",i,j,k,value);
	
	alpha[i][j][k] = value;
	alpha[j][k][i] = value;
	alpha[k][i][j] = value;
	alpha[k][j][i] = value;
	alpha[j][i][k] = value;
	alpha[i][k][j] = value;
	
	return 0;
}

int output_alpha(){
	
	int number = 1;
	int i,j,k;
	
	int *alpha_size = create_ivector(number);
	alpha_size[0] = alpha_max + 1;
	
	int alpha_total_size = alpha_size[0]*alpha_size[0]*alpha_size[0];
	array_write(&alpha_total_size, alpha_data_file, &alpha[0][0][0]);
	
	return 0;
	
}

int read_alpha(){
	
	int number = 1;
	int i,j,k;
	
	int *alpha_size = create_ivector(number);
	alpha_size[0] = alpha_max + 1;
	
	alpha = (double ***)create_3Darray(alpha_size[0],alpha_size[0],alpha_size[0]);
	int alpha_total_size = alpha_size[0]*alpha_size[0]*alpha_size[0];
	array_read(&alpha_total_size, alpha_data_file, &alpha[0][0][0]);
	
	return 0;
	
}

int create_basis_prim(int size, double max, double *vec){
	/*
		Order_prim does not only define the ordering of indices,
		but it also determines what functions will be used as our basis polynomials.
		They can be either: (1-3) Legendre, (4) Fourier or (5) SinLog. 
		Each of them can have an added custom feature to test for inflationary potentials...

		size = number of intervals in discretization of the function domain
		max = maximum value of the argument
		*vec = array of numbers corresponding to the underlying discretization of arguments
	*/
	int i,j,n;
	double x;

	// prim basis is a function of wave vector k and of basis index i	
	basis_prim = (double **)create_array(size,pmax_prim+1);
	
	if(eflag_order_prim<4){
		basis_functions_legendre(basis_prim, size, pmax_prim, 0e0, max, vec);
	}
	else if(eflag_order_prim==4){
		double **basis_temp = (double **)create_array(size,pmax_prim-1);
	 	basis_functions_fourier(basis_temp, size, pmax_prim-2, 0e0, max, vec);
		
		for (i=0;i<size;i++){
			for (j=0;j<pmax_prim-1;j++){
				basis_prim[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_prim[0][pmax_prim-1] = 0.0;
		basis_prim[0][pmax_prim] = 0.0;
		
		for (i=1;i<size;i++){
			x = vec[i]/max;
			basis_prim[i][pmax_prim-1] = 1.0/x;
			basis_prim[i][pmax_prim] = x*x;
		}
	}
	else if(eflag_order_prim==5){
		double min = 2e0/get_tau0();
		basis_functions_sinlog(basis_prim, size, pmax_prim, min, max, vec);
	}
	return 0;
}

double get_basis_prim(int i, int j){
	double value = basis_prim[i][j];
	return value;
}

int destroy_basis_prim(){
	destroy_array(basis_prim);
	return 0;
}

int create_basis_tri_prim(int size, double max, double *vec){

	int i,j,n;
	double x;
	
	basis_tri_prim = (double **)create_array(size,pmax_prim+1);
	
	if(eflag_order_prim<4){
		basis_functions_tri(basis_tri_prim, size, pmax_prim, 0e0, max, vec);
	}else if(eflag_order_prim==5){
		double min = 2e0/get_tau0();
		basis_functions_sinlog(basis_tri_prim, size, pmax_prim, min, max, vec);
	}
	return 0;
}

double get_basis_tri_prim(int i, int j){
	double value = basis_tri_prim[i][j];
	return value;
}

int destroy_basis_tri_prim(){
	destroy_array(basis_tri_prim);
	return 0;
}


int create_basis_late(int Tsize, int Esize, double Tmin, double Tmax, double Emin, double Emax, double *Tvec, double *Evec){
	/* 
		Late time (CMB) basis includes photons which can have polarization as well.
		Order_late defines the ordering of indices, but also the type of polynomial functions one will use
		to construct a basis (e.g. Fourier, Legendre, Sin-Log etc.)
	*/
	int i,j,l,n;
	double c,x,y,z,w;
	double l1,l2;
	
	basis_late_T = (double **)create_array(Tsize,pmax_late_T+1);
	
	if(do_polarisation==1){
		basis_late_E = (double **)create_array(Esize,pmax_late_E+1);
	}
	
	if(eflag_order_late<4){
		basis_functions_bi(basis_late_T, Tsize, pmax_late_T, Tmin, Tmax, Tvec);
		
		if(do_polarisation==1){
			basis_functions_bi(basis_late_E, Esize, pmax_late_E, Emin, Emax, Evec);
		}
		
	}
	else if(eflag_order_late==4){
		
		double **basis_temp = (double **)create_array(Tsize,pmax_late_T-1);
	 	basis_functions_fourier(basis_temp, Tsize, pmax_late_T-2, Tmin, Tmax, Tvec);
		
		for (i=0;i<Tsize;i++){
			for (j=0;j<pmax_late_T-1;j++){
				basis_late_T[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late_T[0][pmax_late_T-1] = 0.0;
		basis_late_T[0][pmax_late_T] = 0.0;
		
		for (i=1;i<Tsize;i++){
			x = Tvec[i]*(Tvec[i]+1e0)/(Tmax*(Tmax+1e0));
			y = (2e0*Tvec[i]+1e0)/(2e0*Tmax+1e0);
			z = sqrt(x);
			w = pow(y,1e0/6e0);
			basis_late_T[i][pmax_late_T-1] = w/z;
			basis_late_T[i][pmax_late_T] = w*z;
		}
		
		if(do_polarisation==1){
			basis_functions_fourier(basis_late_E, Esize, pmax_late_E, Emin, Emax, Evec);
		}
		
	}
	else if(eflag_order_late==5){
		double **basis_temp = (double **)create_array(Tsize,pmax_late_T-4);
	 	basis_functions_fourier(basis_temp, Tsize, pmax_late_T-5, Tmin, Tmax, Tvec);
		
		for (i=0;i<Tsize;i++){
			for (j=0;j<pmax_late_T-4;j++){
				basis_late_T[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late_T[0][pmax_late_T-4] = t_wgt[0];
		basis_late_T[0][pmax_late_T-3] = 0.0;
		basis_late_T[0][pmax_late_T-2] = 0.0;
		basis_late_T[0][pmax_late_T-1] = 0.0;
		basis_late_T[0][pmax_late_T] = 0.0;
		
		for (i=1;i<Tsize;i++){
			x = Tvec[i]*(Tvec[i]+1.0)/(Tmax*(Tmax+1e0));
			y = (2e0*Tvec[i]+1.0)/(2e0*Tmax+1.0);
			z = sqrt(x);
			w = pow(y,1.0/6e0);
			basis_late_T[i][pmax_late_T-4] = t_wgt[i];
			basis_late_T[i][pmax_late_T-3] = w/z;
			basis_late_T[i][pmax_late_T-2] = w*z;
			basis_late_T[i][pmax_late_T-1] = t_wgt[i]*w/z;
			basis_late_T[i][pmax_late_T] = t_wgt[i]*w*z;
		}
		
		if(do_polarisation==1){
			basis_functions_fourier(basis_late_E, Esize, pmax_late_E, Emin, Emax, Evec);
		}
	}
	else if(eflag_order_late==6){
		double **basis_temp = (double **)create_array(Tsize,pmax_late_T-7);
		
	 	basis_functions_fourier(basis_temp, Tsize, pmax_late_T-8, Tmin, Tmax, Tvec);
		
		
		for (i=0;i<Tsize;i++){
			for (j=0;j<pmax_late_T-7;j++){
				basis_late_T[i][j] = basis_temp[i][j];
			}
		}
		
		destroy_array(basis_temp);
		
		basis_late_T[0][pmax_late_T-7] = 0.0;
		basis_late_T[0][pmax_late_T-6] = 0.0;
		basis_late_T[0][pmax_late_T-5] = 0.0;
		basis_late_T[0][pmax_late_T-4] = 0.0;
		basis_late_T[0][pmax_late_T-3] = 0.0;
		basis_late_T[0][pmax_late_T-2] = 0.0;
		basis_late_T[0][pmax_late_T-1] = 0.0;
		basis_late_T[0][pmax_late_T] = 0.0;

		
		for (i=1;i<Tsize;i++){
			l = i+Tmin;
			x = Tvec[i]*(Tvec[i]+1.0)/(Tmax*(Tmax+1e0));
			y = (2e0*Tvec[i]+1.0)/(2e0*Tmax+1.0);
			z = sqrt(x);
			w = pow(y,1.0/6e0);					
			c = sqrt(get_cl_TT(l)+get_noise_TT(l)/(get_beam_TT(l)*get_beam_TT(l)));
			if(c!=0.0)c = w/c;
			basis_late_T[i][pmax_late_T-7] = w/z;
			basis_late_T[i][pmax_late_T-6] = w*z;
			basis_late_T[i][pmax_late_T-5] = c;
			basis_late_T[i][pmax_late_T-4] = c*x;
			basis_late_T[i][pmax_late_T-3] = c*get_cl_TT(l);
			basis_late_T[i][pmax_late_T-2] = c*get_cl_TP(l);
			basis_late_T[i][pmax_late_T-1] = c*x*get_cl_TT(l);
			basis_late_T[i][pmax_late_T] = c*x*get_cl_TP(l);

		}
		
		if(do_polarisation==1){
			basis_functions_fourier(basis_late_E, Esize, pmax_late_E, Emin, Emax, Evec);
		}
	}
	else if(eflag_order_late==7){
		
		double **basis_temp = (double **)create_array(Tsize,pmax_late_T-1);
	 	basis_functions_bi(basis_temp, Tsize, pmax_late_T-2, Tmin, Tmax, Tvec);
		
		for (i=0;i<Tsize;i++){
			for (j=0;j<pmax_late_T-1;j++){
				basis_late_T[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late_T[0][pmax_late_T-1] = 0.0;
		basis_late_T[0][pmax_late_T] = 0.0;
		
		for (i=1;i<Tsize;i++){
			x = Tvec[i]*(Tvec[i]+1e0)/(Tmax*(Tmax+1e0));
			y = (2e0*Tvec[i]+1e0)/(2e0*Tmax+1e0);
			z = sqrt(x);
			w = pow(y,1e0/6e0);
			basis_late_T[i][pmax_late_T-1] = w/z;
			basis_late_T[i][pmax_late_T] = w*z;
		}
		
		if(do_polarisation==1){
			basis_functions_bi(basis_late_E, Esize, pmax_late_E, Emin, Emax, Evec);
		}
	}
	else if(eflag_order_late==8){
		double min = 2e0;
		basis_functions_sinlog(basis_late_T, Tsize, pmax_late_T, Tmin, Tmax, Tvec);
		
		if(do_polarisation==1){
			basis_functions_sinlog(basis_late_E, Esize, pmax_late_E, Emin, Emax, Evec);
		}
	}
	else if(eflag_order_late==9){
		double **basis_temp = (double **)create_array(Tsize,pmax_late_T-1);
		basis_functions_legendre(basis_temp, Tsize, pmax_late_T-2, Tmin, Tmax, Tvec);
		

		for (i=0;i<Tsize;i++){
			for (j=0;j<pmax_late_T-1;j++){
				basis_late_T[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late_T[0][pmax_late_T-1] = 0.0;
		basis_late_T[0][pmax_late_T] = 0.0;
		
		for (i=0;i<Tsize;i++){
			x = Tvec[i]*(Tvec[i]+1e0)/(Tmax*(Tmax+1e0));
			y = (2e0*Tvec[i]+1e0)/(2e0*Tmax+1e0);
			z = sqrt(x);
			w = pow(y,1e0/6e0);
			basis_late_T[i][pmax_late_T-1] = w/z;
			basis_late_T[i][pmax_late_T] = w*z;
		}
		
		if(do_polarisation==1){
			basis_functions_legendre(basis_late_E, Esize, pmax_late_E, Emin, Emax, Evec);
		}
	}
	else if(eflag_order_late==10){
		
		double **basis_temp = (double **)create_array(Tsize,pmax_late_T-1);
		basis_functions_legendre(basis_temp, Tsize, pmax_late_T-4, Tmin, Tmax, Tvec);
		
		for (i=0;i<Tsize;i++){
			for (j=0;j<pmax_late_T-1;j++){
				basis_late_T[i][j] = basis_temp[i][j];
			}
		}
		destroy_array(basis_temp);
		
		basis_late_T[0][pmax_late_T-1] = 0.0;
		basis_late_T[0][pmax_late_T] = 0.0;
		
		for (i=1;i<Tsize;i++){
			x = Tvec[i]*(Tvec[i]+1e0)/(Tmax*(Tmax+1e0));
			y = (2e0*Tvec[i]+1e0)/(2e0*Tmax+1e0);
			z = sqrt(x);
			w = pow(y,1e0/6e0);
			l1 = Tvec[i]/Tmax;
			l2 = l1*l1;

			basis_late_T[i][pmax_late_T-4] = w/(l1*z);
			basis_late_T[i][pmax_late_T-3] = w/z;
			basis_late_T[i][pmax_late_T-2] = w*l1/z;
			basis_late_T[i][pmax_late_T-1] = w*z;
			basis_late_T[i][pmax_late_T] = w*l2*z;
		}
		
		if(do_polarisation==1){
			basis_functions_legendre(basis_late_E, Esize, pmax_late_E, Emin, Emax, Evec);
		}
	}
	return 0;
}

double* get_basis_late_T_flat(){
    return &basis_late_T[0][0];
}

double get_basis_late_T(int i, int j){
	double value = basis_late_T[i][j];
	return value;
}

int destroy_basis_late_T(){
	destroy_array(basis_late_T);
	return 0;
}

double get_basis_late_E(int i, int j){
	double value = basis_late_E[i][j];
	return value;
}

int destroy_basis_late_E(){
	destroy_array(basis_late_E);
	return 0;
}

int create_basis_tri_late(int size, double min, double max, double *vec){

	int i,j,n;
	double x;
	
	basis_tri_late = (double **)create_array(size,pmax_late_T+1);
	
	if(eflag_order_late<4){
		basis_functions_tri(basis_tri_late, size, pmax_late_T, min, max, vec);
	}else if(eflag_order_late==5){
		if(min<=1e0){
			min = 2e0/get_tau0();
		}else{
			min = min/get_tau0();
		}
		basis_functions_sinlog(basis_tri_late, size, pmax_late_T, min, max, vec);
	}
	return 0;
}

double get_basis_tri_late(int i, int j){
	double value = basis_tri_late[i][j];
	return value;
}

int destroy_basis_tri_late(){
	destroy_array(basis_tri_late);
	return 0;
}

//===============================================================================================================
//-----------------------START OF QTILDE FUNCTIONS---------------------------------------------------------------
//===============================================================================================================
int create_qtilde(int lsize_T, int lsize_E, int x_size){
	/*	Create an array qtilde_X[p][l][x]
		where:	
				p = index in the summation (0 <= p <= pmax)
				l = l value (lsize_X of them)
				x = position (up to xmax discretized in x_size elements)
	*/	
	qtilde_T = (double ***)create_3Darray(pmax_prim+1,lsize_T,x_size);
	if(do_polarisation==1){
		qtilde_E = (double ***)create_3Darray(pmax_prim+1,lsize_E,x_size);
	}
	return 0;
}

double get_qtilde_T(int i, int j, int k){
	double value = qtilde_T[i][j][k];
	return value;
}
double get_qtilde_E(int i, int j, int k){
	double value = qtilde_E[i][j][k];
	return value;
}
int update_qtilde_T(int i, int j, int k, double value){
	qtilde_T[i][j][k] = value;
	return 0;
}
int update_qtilde_E(int i, int j, int k, double value){
	qtilde_E[i][j][k] = value;
	return 0;
}

int output_qtilde(int xsize, double *xvec, int lsize_T, int *lvec_T, int lsize_E, int *lvec_E){

	int i,l,n;
	int two = 2;
	char filename[200];
	
	double **data = create_array(lsize_T+1, xsize+1);
	int *sizes = create_ivector(two);
	sizes[0] = lsize_T+1;
	sizes[1] = xsize+1;
	ivector_write(&two, proj_T_size_file, &sizes[0]);
	
	int data_size = (lsize_T+1) * (xsize+1);
	for (i=0; i<xsize; i++) data[0][i+1] = xvec[i];
	for (l=0; l<lsize_T; l++) data[l+1][0] = (double)lvec_T[l];
	
	if(eflag_order_prim!=4){
		for(n=0;n<pmax_prim+1;n++){
			char suffix[3] = "";
			sprintf(suffix, "%d", n);
			
			filename[0] = '\0';
			strcat(filename, proj_T_data_file);
			strcat(filename, "_");
			strcat(filename, suffix);
			
			data[0][0] = (double)n;
			for(i=0; i<xsize; i++){
				for(l=0; l<lsize_T; l++){
					data[l+1][i+1] = qtilde_T[n][l][i];
					//printf("OUTPUT: qtilde_T = %e\n", qtilde_T[n][l][i]);
				}
			}
			array_write(&data_size, filename, &data[0][0]);
		}
	}
	else if(eflag_order_prim==4){
		for(n=0;n<pmax_prim-1;n++){
			char suffix[3] = "";
			sprintf(suffix, "%d", n);
			
			filename[0] = '\0';
			strcat(filename, proj_T_data_file);
			strcat(filename, "_");
			strcat(filename, suffix);
			
			data[0][0] = (double)n;
			for(i=0; i<xsize; i++){
				for(l=0; l<lsize_T; l++){
					data[l+1][i+1] = qtilde_T[n][l][i];
					//printf("OUTPUT: qtilde_T = %e\n", qtilde_T[n][l][i]);
				}
			}
			array_write(&data_size, filename, &data[0][0]);
		}
		filename[0] = '\0';
		strcat(filename, proj_T_data_file);
		strcat(filename, "_l1");
		
		for(i=0; i<xsize; i++){
			for(l=0; l<lsize_T; l++){
				data[l+1][i+1] = qtilde_T[pmax_prim-1][l][i];
			}
		}
		array_write(&data_size, filename, &data[0][0]);
			
		filename[0] = '\0';
		strcat(filename, proj_T_data_file);
		strcat(filename, "_l2");
		
		for(i=0; i<xsize; i++){
			for(l=0; l<lsize_T; l++){
				data[l+1][i+1] = qtilde_T[pmax_prim][l][i];
			}
		}
		array_write(&data_size, filename, &data[0][0]);
	}
	free(data);
	
	if(do_polarisation==1){			//	Repeat same procedure for E-modes
		double **data2 = create_array(lsize_E+1, xsize+1);
		int *sizes2 = create_ivector(two);
		sizes2[0] = lsize_E+1;
		sizes2[1] = xsize+1;
		ivector_write(&two, proj_E_size_file, &sizes2[0]);
		
		int data2_size = (lsize_E+1) * (xsize+1);
		for (i=0; i<xsize; i++) data2[0][i+1] = xvec[i];
		for (l=0; l<lsize_E; l++) data2[l+1][0] = (double)lvec_E[l];
		
		if(eflag_order_prim!=4){
			for(n=0;n<pmax_prim+1;n++){
				char suffix[3] = "";
				sprintf(suffix, "%d", n);
			
				filename[0] = '\0';
				strcat(filename, proj_E_data_file);
				strcat(filename, "_");
				strcat(filename, suffix);

				data2[0][0] = (double)n;
				for(i=0; i<xsize; i++){
					for(l=0; l<lsize_E; l++){
						data2[l+1][i+1] = qtilde_E[n][l][i];
					}
				}
				printf("Array write: %s\n", filename);
				array_write(&data2_size, filename, &data2[0][0]);
			}
		}else if(eflag_order_prim==4){
			for(n=0;n<pmax_prim-1;n++){
				char suffix[3] = "";
				sprintf(suffix, "%d", n);
			
				filename[0] = '\0';
				strcat(filename, proj_E_data_file);
				strcat(filename, "_");
				strcat(filename, suffix);
			
				data2[0][0] = (double)n;
				for(i=0; i<xsize; i++){
					for(l=0; l<lsize_E; l++){
						data2[l+1][i+1] = qtilde_E[n][l][i];
					}
				}
				array_write(&data2_size, filename, &data2[0][0]);
			}
			filename[0] = '\0';
			strcat(filename, proj_E_data_file);
			strcat(filename, "_l1");
		
			for(i=0; i<xsize; i++){
				for(l=0; l<lsize_E; l++){
					data2[l+1][i+1] = qtilde_E[pmax_prim-1][l][i];
				}
			}
			array_write(&data2_size, filename, &data2[0][0]);
			
			filename[0] = '\0';
			strcat(filename, proj_E_data_file);
			strcat(filename, "_l2");
		
			for(i=0; i<xsize; i++){
				for(l=0; l<lsize_E; l++){
					data2[l+1][i+1] = qtilde_E[pmax_prim][l][i];
				}
			}
			array_write(&data2_size, filename, &data2[0][0]);
		}
		free(data2);
	}
	return 0;
}

int read_qtilde(){
	int number = 2;
	int *sizes = create_ivector(number);
	ivector_read(&number, proj_T_size_file, &sizes[0]);
	qtilde_T_lsize = sizes[0]-1;
	qtilde_xsize = sizes[1]-1;
	qtilde_T_lvec = create_ivector(qtilde_T_lsize);
	qtilde_xvec = create_vector(qtilde_xsize);
	int data_size = (qtilde_T_lsize+1) * (qtilde_xsize+1);
	double **data = create_array(qtilde_T_lsize+1, qtilde_xsize+1);
	char filename[200];
	//printf("Read qtilde 1\n");
	qtilde_E_lsize = 0;
	if(do_polarisation==1){
		ivector_read(&number, proj_E_size_file, &sizes[0]);
		qtilde_E_lsize = sizes[0]-1;
		qtilde_E_lvec = create_ivector(qtilde_E_lsize);
	}
	//printf("Read qtilde 2\n");
	create_qtilde(qtilde_T_lsize, qtilde_E_lsize, qtilde_xsize);
	int n,l,i;
	
	if(eflag_order_prim!=4){
		for(n=0;n<pmax_prim+1;n++){
			char suffix[3] = "";
			sprintf(suffix, "%d", n);
			
			filename[0] = '\0';
			strcat(filename, proj_T_data_file);
			strcat(filename, "_");
			strcat(filename, suffix);
			//printf("Read qtilde 3\n");
			array_read(&data_size, filename, &data[0][0]);
			
			if(n==0){
				for (i=0; i<qtilde_xsize; i++) qtilde_xvec[i] = data[0][i+1];
				for (l=0; l<qtilde_T_lsize; l++) qtilde_T_lvec[l] = (int)data[l+1][0];
				
			}
	
			for(l=0; l<qtilde_T_lsize; l++) qtilde_T[n][l][0] = 0.0;
			
			for(i=0; i<qtilde_xsize; i++){
				for(l=0; l<qtilde_T_lsize; l++){
					qtilde_T[n][l][i] = data[l+1][i+1];
					//printf("READ: qtilde_T = %e\n", qtilde_T[n][l][i]);
				}
			}
			
		}
	}
	else{
		for(n=0;n<pmax_prim-1;n++){
			char suffix[3] = "";
			sprintf(suffix, "%d", n);
			
			filename[0] = '\0';
			strcat(filename, proj_T_data_file);
			strcat(filename, "_");
			strcat(filename, suffix);
			//printf("Filename %s\n",filename);
			fflush;
			array_read(&data_size, filename, &data[0][0]);
			
			if(n==0){
				for (i=0; i<qtilde_xsize; i++) qtilde_xvec[i] = data[0][i+1];
				for (l=0; l<qtilde_T_lsize; l++) qtilde_T_lvec[l] = (int)data[l+1][0];
			}
	
			for(l=0; l<qtilde_T_lsize; l++) qtilde_T[n][l][0] = 0.0;
			
			for(i=0; i<qtilde_xsize; i++){
				for(l=0; l<qtilde_T_lsize; l++){
					qtilde_T[n][l][i] = data[l+1][i+1];
					//printf("READ: qtilde_T = %e\n", qtilde_T[n][l][i]);
				}
			}
		}
			
		filename[0] = '\0';
		strcat(filename, proj_T_data_file);
		strcat(filename, "_l1");
			
		array_read(&data_size, filename, &data[0][0]);
	
		for(l=0; l<qtilde_T_lsize; l++) qtilde_T[pmax_prim-1][l][0] = 0.0;
			
		for(i=0; i<qtilde_xsize; i++){
			for(l=0; l<qtilde_T_lsize; l++){
				qtilde_T[pmax_prim-1][l][i] = data[l+1][i+1];
			}
		}
		filename[0] = '\0';
		strcat(filename, proj_T_data_file);
		strcat(filename, "_l2");
			
		array_read(&data_size, filename, &data[0][0]);
	
		for(l=0; l<qtilde_T_lsize; l++) qtilde_T[pmax_prim][l][0] = 0.0;
			
		for(i=0; i<qtilde_xsize; i++){
			for(l=0; l<qtilde_T_lsize; l++){
				qtilde_T[pmax_prim][l][i] = data[l+1][i+1];
			}
		}
	}
	free(data);
	
	if(do_polarisation==1){		//	Repeat same procedure for E-modes
		int data2_size = (qtilde_E_lsize+1) * (qtilde_xsize+1);
		double **data2 = create_array(qtilde_E_lsize+1, qtilde_xsize+1);
		char filenameE[200];
	
		if(eflag_order_prim!=4){
			for(n=0;n<pmax_prim+1;n++){
				char suffix[3] = "";
				sprintf(suffix, "%d", n);
			
				filenameE[0] = '\0';
				strcat(filenameE, proj_E_data_file);
				strcat(filenameE, "_");
				strcat(filenameE, suffix);
				//printf("Read qtilde 4\n");
				printf("Array read: %s\n", filenameE);
				array_read(&data2_size, filenameE, &data2[0][0]);
			
				if(n==0){
					for (i=0; i<qtilde_xsize; i++) qtilde_xvec[i] = data2[0][i+1];
					for (l=0; l<qtilde_E_lsize; l++) qtilde_E_lvec[l] = (int)data2[l+1][0];
				
				}

				for(l=0; l<qtilde_E_lsize; l++) qtilde_E[n][l][0] = 0.0;
			
				for(i=0; i<qtilde_xsize; i++){
					for(l=0; l<qtilde_E_lsize; l++){
						qtilde_E[n][l][i] = data2[l+1][i+1];
					}
				}
			
			}
			//printf("Read qtilde 5\n");
		}
		else{
			for(n=0;n<pmax_prim-1;n++){
				char suffix[3] = "";
				sprintf(suffix, "%d", n);
			
				filenameE[0] = '\0';
				strcat(filenameE, proj_E_data_file);
				strcat(filenameE, "_");
				strcat(filenameE, suffix);
			
				array_read(&data2_size, filenameE, &data2[0][0]);
			
				if(n==0){
					for (i=0; i<qtilde_xsize; i++) qtilde_xvec[i] = data2[0][i+1];
					for (l=0; l<qtilde_E_lsize; l++) qtilde_E_lvec[l] = (int)data2[l+1][0];
				}
	
				for(l=0; l<qtilde_E_lsize; l++) qtilde_E[n][l][0] = 0.0;
			
				for(i=0; i<qtilde_xsize; i++){
					for(l=0; l<qtilde_E_lsize; l++){
						qtilde_E[n][l][i] = data2[l+1][i+1];
					}
				}
			
			}
			
			filenameE[0] = '\0';
			strcat(filenameE, proj_E_data_file);
			strcat(filenameE, "_l1");
			
			array_read(&data2_size, filenameE, &data2[0][0]);
	
			for(l=0; l<qtilde_E_lsize; l++) qtilde_E[pmax_prim-1][l][0] = 0.0;
			
			for(i=0; i<qtilde_xsize; i++){
				for(l=0; l<qtilde_E_lsize; l++){
					qtilde_E[pmax_prim-1][l][i] = data2[l+1][i+1];
				}
			}
			
			filenameE[0] = '\0';
			strcat(filenameE, proj_E_data_file);
			strcat(filenameE, "_l2");
			
			array_read(&data2_size, filenameE, &data2[0][0]);
	
			for(l=0; l<qtilde_E_lsize; l++) qtilde_E[pmax_prim][l][0] = 0.0;
			
			for(i=0; i<qtilde_xsize; i++){
				for(l=0; l<qtilde_E_lsize; l++){
					qtilde_E[pmax_prim][l][i] = data2[l+1][i+1];
				}
			}
		}
		free(data2);
	}

	/*#pragma offload_transfer target(mic:offload_target) \
	    in(qtilde_xsize) \
	    in(qtilde_xvec[0:qtilde_xsize] : ALLOC RETAIN)*/
	return 0;
}

int destroy_qtilde(){

	destroy_3Darray(qtilde_T);
	if(do_polarisation==1){
		destroy_3Darray(qtilde_E);
	}
	
	return 0;
}

int get_qtilde_T_lsize(){
	return qtilde_T_lsize;
}

int get_qtilde_E_lsize(){
	return qtilde_E_lsize;
}
 
int get_qtilde_xsize(){
	return qtilde_xsize;
}

int get_qtilde_T_lvec(int *vec){
	
	int i;
	for(i=0;i<qtilde_T_lsize;i++){
		vec[i] = qtilde_T_lvec[i];
	}
	return 0;

}

int get_qtilde_E_lvec(int *vec){
	
	int i;
	for(i=0;i<qtilde_E_lsize;i++){
		vec[i] = qtilde_E_lvec[i];
	}
	return 0;

}


int get_qtilde_xvec(double *vec){

	int i;
	for(i=0;i<qtilde_xsize;i++){
		vec[i] = qtilde_xvec[i];
	}
	return 0;
}
//===============================================================================================================
//-----------------------END OF QTILDE FUNCTIONS---------------------------------------------------------------
//===============================================================================================================


//===============================================================================================================
//-------------------------START OF BETA FUNCTIONS---------------------------------------------------------------
//===============================================================================================================
int create_beta_tri(int l_size, int x_size){
	beta_tri = (double ***)create_3Darray(pmax_prim+1,l_size,x_size);
	return 0;
}

double get_beta_tri(int i, int j, int k){
	double value = beta_tri[i][j][k];
	return value;
}
	
int update_beta_tri(int size, double *results){

	int l,a,k;
	
	a = (int)results[0];
	l = (int)results[1];
	
	for (k=2;k<size;k++){
		beta_tri[a][l][k-2] = results[k];
	}
	
// 	printf("test: %d\t%d\n", l,a);
// 	if(a==5)printf("test: %d\t%d\t%e\n", l,a,beta_tri[a][l][100]);
	return 0;
}

int send_beta_tri(int nproc, int rank, int l_size, int x_size){

	int size = l_size*(pmax_prim+1)*x_size;
	int i,j,k,n;
	MPI_Status status;
	
	double ***array = create_3Darray(pmax_prim+1,l_size,x_size);
	
	for (n=1;n<nproc;n++){
		if (rank == 0){
			MPI_Send(&beta_tri[0][0][0], size, MPI_DOUBLE, n, 0, MPI_COMM_WORLD);
		}
		
		if (rank == n){
			MPI_Recv(&array[0][0][0], size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			
			for(i=0;i<pmax_prim+1;i++){
				for(j=0;j<l_size;j++){
					for(k=0;k<x_size;k++){
					beta_tri[i][j][k] = array[i][j][k];
					}
				}
			}
			
		}
	}
	
	destroy_3Darray(array);
	return 0;
}

int destroy_beta_tri(int l_size){

	destroy_3Darray(beta_tri);
	
	return 0;
}

int output_beta_tri(int xsize, double *xvec, int lsize, int *lvec){

	int i,l,n;
	double **data = create_array(lsize+1, xsize);
	int two = 2;
	
	int *sizes = create_ivector(two);
	sizes[0] = lsize+1;
	sizes[1] = xsize;
	ivector_write(&two, proj_tri_size_file, &sizes[0]);
	
	int data_size = (lsize+1) * xsize;
	for (i=1; i<xsize; i++) data[0][i] = xvec[i];
	for (l=0; l<lsize; l++) data[l+1][0] = (double)lvec[l];
// 	for (l=0; l<lsize; l++) printf("%d\n",lvec[l]);
	
	for(n=0;n<pmax_prim+1;n++){
		char suffix[3] = "";
		sprintf(suffix, "%d", n);
		
		char filename[100] = "";
		strcat(filename, proj_tri_data_file);
		strcat(filename, suffix);
		
		data[0][0] = (double)n;
		for(i=1; i<xsize; i++){
			for(l=0; l<lsize; l++){
				data[l+1][i] = beta_tri[n][l][i];
			}
		}
		
// 		for(l=0; l<lsize; l++){
// 			printf("test: %e\n", beta_tri[n][l][xsize-1]);
// 		}
		
		array_write(&data_size, filename, &data[0][0]);
		
	}
	
	return 0;
	
}

int read_beta_tri(){

	int number = 2;
	int *sizes = create_ivector(number);
	ivector_read(&number, proj_tri_size_file, &sizes[0]);
	bt_lsize = sizes[0]-1;
	bt_xsize = sizes[1]-1;
	bt_lvec = create_ivector(bt_lsize);
	bt_xvec = create_vector(bt_xsize);
	int data_size = (bt_lsize+1) * (bt_xsize+1);
	double **data = create_array(bt_lsize+1, bt_xsize+1);
	
	create_beta_tri(bt_lsize, bt_xsize);
	
	int n,l,i;
	for(n=0;n<pmax_prim+1;n++){
		char suffix[3] = "";
		sprintf(suffix, "%d", n);
		
		char filename[100] = "";
		strcat(filename, proj_tri_data_file);
		strcat(filename, suffix);
		
		array_read(&data_size, filename, &data[0][0]);
		
		if(n==0){
			for (i=0; i<bt_xsize; i++) bt_xvec[i] = data[0][i+1];
			for (l=0; l<bt_lsize; l++) bt_lvec[l] = (int)data[l+1][0];
		}

		for(l=0; l<bt_lsize; l++) beta_tri[n][l][0] = 0.0;
		
		for(i=0; i<bt_xsize; i++){
			for(l=0; l<bt_lsize; l++){
				beta_tri[n][l][i] = data[l+1][i+1];
			}
		}
		
	}
	return 0;
	
}

int get_bt_lsize(){
	return bt_lsize;
}

int get_bt_xsize(){
	return bt_xsize;
}

int get_bt_lvec(int *vec){
	
	int i;
	for(i=0;i<bt_lsize;i++){
		vec[i] = bt_lvec[i];
	}
	return 0;

}

int get_bt_xvec(double *vec){

	int i;
	for(i=0;i<bt_xsize;i++){
		vec[i] = bt_xvec[i];
	}
	return 0;
}


int create_bispectrum(int size){

	bispectrum = (double ***)create_3Darray(size,size,size);
	
	return 0;
}

double get_bispectrum(int i, int j, int k){
	double value = bispectrum[i][j][k];
	return value;
}
	
int update_bispectrum(double *results){

	int i,j,k;
	double area;
	
	i = (int)results[0];
	j = (int)results[1];
	k = (int)results[2];
	area = results[3];
	
// 	printf("%d\t%d\t%d\t%e\n",i,j,k,area);
	
	bispectrum[i][j][k] = area;
	bispectrum[j][k][i] = area;
	bispectrum[k][i][j] = area;
	bispectrum[k][j][i] = area;
	bispectrum[j][i][k] = area;
	bispectrum[i][k][j] = area;
	
	return 0;
}

int output_bispectrum(int size, int *lvalues){
	
	int number = 1;
	int i,j,k;
	
	int *l_size = create_ivector(number);
	l_size[0] = size;	
	ivector_write(&number, l_size_file, &l_size[0]);
	
	int *l_values = create_ivector(size);
	for (i=0;i<size;i++){l_values[i] = lvalues[i];}
	ivector_write(&size, l_data_file, &l_values[0]);
	
	int bispectrum_size = size*size*size;
	array_write(&bispectrum_size, bispectrum_file, &bispectrum[0][0][0]);
	
	destroy_ivector(l_values);
	
	return 0;
	
}

int create_trispectrum(int size){

	trispectrum = (double ****)create_4Darray(size,size,size,size);
	
	return 0;
}

double get_trispectrum(int i, int j, int k, int l){
	double value = trispectrum[i][j][k][l];
	return value;
}
	
int update_trispectrum(double *results){

	int i,j,k,l;
	double area;
	
	i = (int)results[0];
	j = (int)results[1];
	k = (int)results[2];
	l = (int)results[3];
	area = results[4];
	
	trispectrum[i][j][k][l] = area;
	trispectrum[j][k][i][l] = area;
	trispectrum[k][i][j][l] = area;
	trispectrum[k][j][i][l] = area;
	trispectrum[j][i][k][l] = area;
	trispectrum[i][k][j][l] = area;
	
	trispectrum[i][j][l][k] = area;
	trispectrum[j][k][l][i] = area;
	trispectrum[k][i][l][j] = area;
	trispectrum[k][j][l][i] = area;
	trispectrum[j][i][l][k] = area;
	trispectrum[i][k][l][j] = area;
	
	trispectrum[i][l][j][k] = area;
	trispectrum[j][l][k][i] = area;
	trispectrum[k][l][i][j] = area;
	trispectrum[k][l][j][i] = area;
	trispectrum[j][l][i][k] = area;
	trispectrum[i][l][k][j] = area;
	
	trispectrum[l][i][j][k] = area;
	trispectrum[l][j][k][i] = area;
	trispectrum[l][k][i][j] = area;
	trispectrum[l][k][j][i] = area;
	trispectrum[l][j][i][k] = area;
	trispectrum[l][i][k][j] = area;
	
	return 0;
}

int output_trispectrum(int size, int *lvalues){
	
	int number = 1;
	int i,j,k,l;
	
	int *l_size = create_ivector(number);
	l_size[0] = size;	
	ivector_write(&number, l_size_file, &l_size[0]);
	
	int *l_values = create_ivector(size);
	for (i=0;i<size;i++){l_values[i] = lvalues[i];}
	ivector_write(&size, l_data_file, &l_values[0]);
	
	int trispectrum_size = size*size*size*size;
	array_write(&trispectrum_size, trispectrum_file, &trispectrum[0][0][0][0]);
	
	destroy_ivector(l_values);
	
	return 0;
	
}

int create_decompose(){
	decompose = (double ***)create_3Darray(alpha_max+1,alpha_max+1,alpha_max+1);
	return 0;
}

double get_decompose(int i, int j, int k){
	double value = decompose[i][j][k];
	return value;
}
	
int update_decompose(double *results){

	int i,j,k;
	double value;
	
	i = (int)results[0];
	j = (int)results[1];
	k = (int)results[2];
	value = results[3];
	
	//printf("update: %d\t%d\t%d\t%e\n",i,j,k,value);
	
	decompose[i][j][k] = value;
	decompose[j][k][i] = value;
	decompose[k][i][j] = value;
	decompose[k][j][i] = value;
	decompose[j][i][k] = value;
	decompose[i][k][j] = value;
	
	return 0;
}

int output_decompose(){
	
	int number = 1;
	int i,j,k;
	
	int *decompose_size = create_ivector(number);
	decompose_size[0] = alpha_max + 1;	
	ivector_write(&number, decompose_size_file, &decompose_size[0]);
	
	int decompose_total_size = decompose_size[0]*decompose_size[0]*decompose_size[0];
	array_write(&decompose_total_size, decompose_data_file, &decompose[0][0][0]);
	
	return 0;
	
}

int read_decompose(char *directory){
	
	int number = 1;
	int i,j,k;
	
	int *decompose_size = create_ivector(number);
	ivector_read(&number, decompose_size_file, &decompose_size[0]);
	
	if (decompose_size[0] != alpha_max + 1){
		printf("alpha max in parameter does not match alpha max in loaded file exiting.\n");
		exit(1);
	}
	
	decompose = (double ***)create_3Darray(decompose_size[0],decompose_size[0],decompose_size[0]);
	int decompose_total_size = decompose_size[0]*decompose_size[0]*decompose_size[0];
	array_read(&decompose_total_size, decompose_data_file, &decompose[0][0][0]);
	
	return 0;
	
}

int create_eigen(){

	int n = terms_prim;
	eigen = (double *)create_vector(n);
	
	return 0;
}

double get_eigen(int i){
	double value = eigen[i];
	return value;
}
	
int update_eigen(double *results){

	int n;
	
	n = (int)results[0];
	
	eigen[n] = results[1];
	
	return 0;
}

int output_eigen(int p1, int p2, int p3){

	char strmodel[10] = "";
	char strp1[10] = "";
	char strp2[10] = "";
	char strp3[10] = "";
	sprintf(strmodel,"_%d",model);
	if(slim.p1_end!=0)sprintf(strp1,"_%d",p1);
	if(slim.p2_end!=0)sprintf(strp2,"_%d",p2);
	if(slim.p3_end!=0)sprintf(strp3,"_%d",p3);
	
	char eigen_output_file[200] = "";
	strcat(eigen_output_file, eigen_data_file);
	strcat(eigen_output_file, strmodel);
	strcat(eigen_output_file, strp1);
	strcat(eigen_output_file, strp2);
	strcat(eigen_output_file, strp3);
	
	int number = 1;
	int i,j,k;
	
	int *eigen_size = create_ivector(number);
	eigen_size[0] = terms_prim;
	
	int eigen_total_size = eigen_size[0];
	vector_write(&eigen_total_size, eigen_output_file, &eigen[0]);
	printf("%s\n", eigen_output_file);
	return 0;
	
}

int read_eigen(int p1, int p2, int p3){

	char strmodel[10] = "";
	char strp1[10] = "";
	char strp2[10] = "";
	char strp3[10] = "";
	sprintf(strmodel,"_%d",model);
	if(slim.p3_end!=0){
		sprintf(strp3,"_%d",p3);
		if(slim.p2_end!=0){
			sprintf(strp2,"_%d",p2);
			if(slim.p1_end!=0){
				sprintf(strp1,"_%d",p1);
			}else{
				sprintf(strp1,"_%d",0);
			}
		}else{
			sprintf(strp2,"_%d",0);
			if(slim.p1_end!=0){
				sprintf(strp1,"_%d",p1);
			}else{
				sprintf(strp1,"_%d",0);
			}
		}
	}else{
		if(slim.p2_end!=0){
			sprintf(strp2,"_%d",p2);
			if(slim.p1_end!=0){
				sprintf(strp1,"_%d",p1);
			}else{
				sprintf(strp1,"_%d",0);
			}
		}else{
			if(slim.p1_end!=0){
				sprintf(strp1,"_%d",p1);
			}
		}
	}
	
	char eigen_output_file[200] = "";
	strcat(eigen_output_file, eigen_data_file);
	strcat(eigen_output_file, strmodel);
	strcat(eigen_output_file, strp1);
	strcat(eigen_output_file, strp2);
	strcat(eigen_output_file, strp3);

	printf("File = %s\n", eigen_output_file);
	
	int number = 1;
	int i,j,k;
	
	int *eigen_size = create_ivector(number);
	eigen_size[0] = terms_prim;
	
	eigen = (double *)create_vector(eigen_size[0]);
	int eigen_total_size = eigen_size[0];
	vector_read(&eigen_total_size, eigen_output_file, &eigen[0]);
	
	return 0;
	
}

int create_eigenR(){

	int n = terms_prim;
	eigenR = (double *)create_vector(n);
	
	return 0;
}

double get_eigenR(int i){
	double value = eigenR[i];
	return value;
}
	
int update_eigenR(double *results){

	int n;
	n = (int)results[0];
	eigenR[n] = results[1];
	return 0;
}

/* -------------------------------------------------------------------
	Create a gamma matrix for each possible polarisation combination
---------------------------------------------------------------------*/
int create_gamma(){

	int n = terms_prim;
	gamma_TTT = (double **)create_array(n,n);
	if(do_polarisation==1){
		gamma_TTE = (double **)create_array(n,n);
		gamma_TEE = (double **)create_array(n,n);
		gamma_EEE = (double **)create_array(n,n);
	}
	return 0;
}

/* ----------------------------------------------------------------
	Return the elements of gamma_XXX matrix
-------------------------------------------------------------------*/
double get_gamma_TTT(int i,int j){
	double value = gamma_TTT[i][j];
	return value;
}
double get_gamma_TTE(int i,int j){
	double value = gamma_TTE[i][j];
	return value;
}
double get_gamma_TEE(int i,int j){
	double value = gamma_TEE[i][j];
	return value;
}
double get_gamma_EEE(int i,int j){
	double value = gamma_EEE[i][j];
	return value;
}

/*	-----------------------------------------------------------------------
	Transfer the values from the array "results" into a static gamma matrix
-------------------------------------------------------------------------*/
int update_gamma(double *results){

/*	First two elements of the "results" always carry the information
	about which element in the gamma matrix we should acces.
	Elements indexed 2 -- 5 carry the value for each polarization mode.

	This function saves temporary results in a static variable gamma matrix.
*/
	int i,j;
	i = results[0];
	j = results[1];
	gamma_TTT[i][j] = results[2];
	if(do_polarisation==1){
		gamma_TTE[i][j] = results[3];
		gamma_TEE[i][j] = results[4];
		gamma_EEE[i][j] = results[5];
	}
	return 0;
}

/*---------------------------------------------------------------------
	Print the elements of the gamma matrix in a 
	dedicated file. This is useful because one can re-use
	the calculations from previously covered models and parameters,
	without the need to re-run the entire set of codes.
-------------------------------------------------------------------------*/
int output_gamma(){
	
	int number = 1;
	int i,j,k;
//	Print gamma_TTT into a gamma_TTT_data_file etc.
	int gamma_total_size = terms_prim*terms_prim;	//matrix size
	array_write(&gamma_total_size, gamma_TTT_data_file, &gamma_TTT[0][0]);
	
	if(do_polarisation==1){
		array_write(&gamma_total_size, gamma_TTE_data_file, &gamma_TTE[0][0]);
		array_write(&gamma_total_size, gamma_TEE_data_file, &gamma_TEE[0][0]);
		array_write(&gamma_total_size, gamma_EEE_data_file, &gamma_EEE[0][0]);
	}
	return 0;
}

/*---------------------------------------------------------------------
	Read the gamma_XXX_data_file that has been computed before, and 
	store its values in a static variable gamma_XXX
-------------------------------------------------------------------------*/
int read_gamma(){
	
	int number = 1;
	int i,j,k;
	
	int gamma_total_size = terms_prim*terms_prim;	//matrix size
	
	gamma_TTT = (double **)create_array(terms_prim,terms_prim);
	array_read(&gamma_total_size, gamma_TTT_data_file, &gamma_TTT[0][0]);
	
	if(do_polarisation==1){
		gamma_TTE = (double **)create_array(terms_prim,terms_prim);
		gamma_TEE = (double **)create_array(terms_prim,terms_prim);
		gamma_EEE = (double **)create_array(terms_prim,terms_prim);
		array_read(&gamma_total_size, gamma_TTE_data_file, &gamma_TTE[0][0]);
		array_read(&gamma_total_size, gamma_TEE_data_file, &gamma_TEE[0][0]);
		array_read(&gamma_total_size, gamma_EEE_data_file, &gamma_EEE[0][0]);
	}
	return 0;
	
}







int create_eigen_tri(){

	int n = terms_prim;
	eigen_tri = (double *)create_vector(n);
	
	return 0;
}

double get_eigen_tri(int i){
	double value = eigen_tri[i];
	return value;
}
	
int update_eigen_tri(double *results){

	int n;
	
	n = (int)results[0];
	
	eigen_tri[n] = results[1];
	
	return 0;
}

int output_eigen_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *eigen_tri_size = create_ivector(number);
	eigen_tri_size[0] = terms_prim;	
	
	int eigen_tri_total_size = eigen_tri_size[0];
	vector_write(&eigen_tri_total_size, eigen_tri_data_file, &eigen_tri[0]);
	
	return 0;
	
}

int read_eigen_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *eigen_tri_size = create_ivector(number);
	eigen_tri_size[0] = terms_prim;	
	
	eigen_tri = (double *)create_vector(eigen_tri_size[0]);
	int eigen_tri_total_size = eigen_tri_size[0];
	vector_read(&eigen_tri_total_size, eigen_tri_data_file, &eigen_tri[0]);
	
	return 0;
	
}

int destroy_eigen_tri(){
	destroy_vector(eigen_tri);
	return 0;
}

int create_eigenR_tri(){

	int n = terms_prim;
	eigenR_tri = (double *)create_vector(n);
	
	return 0;
}

double get_eigenR_tri(int i){
	double value = eigenR_tri[i];
	return value;
}
	
int update_eigenR_tri(double *results){

	int n;
	
	n = (int)results[0];
	
	eigenR_tri[n] = results[1];
	
	return 0;
}

int create_gamma_tri(){
	
	int n = terms_prim;
	gamma_tri = (double **)create_array(n,n);
	
	return 0;
}

double get_gamma_tri(int i,int j){
	double value = gamma_tri[i][j];
	return value;
}

int update_gamma_tri(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	gamma_tri[i][j] = results[2];
	
	return 0;
}

int output_gamma_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *gamma_tri_size = create_ivector(number);
	gamma_tri_size[0] = terms_prim;
	
	int gamma_tri_total_size = gamma_tri_size[0]*gamma_tri_size[0];
	array_write(&gamma_tri_total_size, gamma_tri_data_file, &gamma_tri[0][0]);
	
	return 0;

}

int read_gamma_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *gamma_tri_size = create_ivector(number);
	gamma_tri_size[0] = terms_prim;
	
	double **gamma1_tri = (double **)create_array(gamma_tri_size[0],gamma_tri_size[0]);
	int gamma_tri_total_size = gamma_tri_size[0]*gamma_tri_size[0];
	array_read(&gamma_tri_total_size, gamma_tri_data_file, &gamma1_tri[0][0]);
	
	int size = terms_prim;
	gamma_tri = (double **)create_array(size,size);
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			gamma_tri[i][j] = gamma1_tri[i][j];
		}
	}
	
	return 0;
	
}

int create_ortho(){
// creates a square matrix with number of terms × number of terms
	int n = terms_prim;
	ortho = (double **)create_array(n,n);
	return 0;
}

double get_ortho(int i,int j){
	double value = ortho[i][j];
	return value;
}
	
int update_ortho(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	ortho[i][j] = results[2];
	ortho[j][i] = results[2];
	
	return 0;
}

int output_ortho(){
	
	int number = 1;
	int i,j,k;
	
	int *ortho_size = create_ivector(number);
	ortho_size[0] = terms_prim;
	
	int ortho_total_size = ortho_size[0]*ortho_size[0];
	array_write(&ortho_total_size, ortho_data_file, &ortho[0][0]);
	
	return 0;
	
}

int read_ortho(char *directory){
	
	int number = 1;
	int i,j,k;
	
	int *ortho_size = create_ivector(number);
	ortho_size[0] = terms_prim;
	
	double **ortho1 = (double **)create_array(ortho_size[0],ortho_size[0]);
	int ortho_total_size = ortho_size[0]*ortho_size[0];
	array_read(&ortho_total_size, ortho_data_file, &ortho1[0][0]);
	printf("%s\n", ortho_data_file);
	
	int size = terms_prim;
	ortho = (double **)create_array(size,size);
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			ortho[i][j] = ortho1[i][j];
		}
	}
	
	return 0;
	
}

int create_ortho_tri(){
	
	int n = terms_prim;
	ortho_tri = (double **)create_array(n,n);
	
	return 0;
}

double get_ortho_tri(int i,int j){
	double value = ortho_tri[i][j];
	return value;
}
	
int update_ortho_tri(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	ortho_tri[i][j] = results[2];
	ortho_tri[j][i] = results[2];
	
	return 0;
}

int output_ortho_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *ortho_tri_size = create_ivector(number);
	ortho_tri_size[0] = terms_prim;
	
	int ortho_tri_total_size = ortho_tri_size[0]*ortho_tri_size[0];
	array_write(&ortho_tri_total_size, ortho_tri_data_file, &ortho_tri[0][0]);
	
	return 0;
	
}

int read_ortho_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *ortho_tri_size = create_ivector(number);
	ortho_tri_size[0] = terms_prim;
	
	double **ortho1 = (double **)create_array(ortho_tri_size[0],ortho_tri_size[0]);
	int ortho_tri_total_size = ortho_tri_size[0]*ortho_tri_size[0];
	array_read(&ortho_tri_total_size, ortho_tri_data_file, &ortho1[0][0]);
	
	int size = terms_prim;
	ortho_tri = (double **)create_array(size,size);
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			ortho_tri[i][j] = ortho1[i][j];
		}
	}
	
	return 0;
	
}

int create_orthol(){
	
	int n = terms_late;
	orthol_TTT = (double **)create_array(n,n);
	
	if(do_polarisation==1){
		orthol_TTE = (double **)create_array(n,n);
		orthol_TEE = (double **)create_array(n,n);
		orthol_EEE = (double **)create_array(n,n);
	}
	
	return 0;
}

double get_orthol_TTT(int i,int j){
	double value = orthol_TTT[i][j];
	return value;
}

double get_orthol_TTE(int i,int j){
	double value = orthol_TTE[i][j];
	return value;
}

double get_orthol_TEE(int i,int j){
	double value = orthol_TEE[i][j];
	return value;
}

double get_orthol_EEE(int i,int j){
	double value = orthol_EEE[i][j];
	return value;
}
	
int update_orthol(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	orthol_TTT[i][j] = results[2];
	if(do_polarisation==1){
		orthol_TTE[i][j] = results[3];
		orthol_TEE[i][j] = results[4];
		orthol_EEE[i][j] = results[5];
	}
	
	return 0;
}

int output_orthol(){
	
	int number = 1;
	int i,j,k;
	
	int orthol_total_size = terms_late*terms_late;
	array_write(&orthol_total_size, orthol_TTT_data_file, &orthol_TTT[0][0]);
	
	if(do_polarisation==1){
		array_write(&orthol_total_size, orthol_TTE_data_file, &orthol_TTE[0][0]);
		array_write(&orthol_total_size, orthol_TEE_data_file, &orthol_TEE[0][0]);
		array_write(&orthol_total_size, orthol_EEE_data_file, &orthol_EEE[0][0]);
	}
	
	return 0;
	
}

/*
	If one has previously calculated orthol_XXX factors, they will have been saved
	in orthol_XXX_data_file.
*/
int read_orthol(char *directory){

	int number = 1;
	int i,j,k;
	
	int orthol_total_size = terms_late*terms_late;
	orthol_TTT = (double **)create_array(terms_late,terms_late);
	array_read(&orthol_total_size, orthol_TTT_data_file, &orthol_TTT[0][0]);
	
	if(do_polarisation==1){
		orthol_TTE = (double **)create_array(terms_late,terms_late);
		orthol_TEE = (double **)create_array(terms_late,terms_late);
		orthol_EEE = (double **)create_array(terms_late,terms_late);
		array_read(&orthol_total_size, orthol_TTE_data_file, &orthol_TTE[0][0]);
		array_read(&orthol_total_size, orthol_TEE_data_file, &orthol_TEE[0][0]);
		array_read(&orthol_total_size, orthol_EEE_data_file, &orthol_EEE[0][0]);
	}
	return 0;
}

int create_orthol_tri(){
	
	int n = terms_late;
	orthol_tri = (double **)create_array(n,n);
	
	return 0;
}

double get_orthol_tri(int i,int j){
	double value = orthol_tri[i][j];
	return value;
}

int update_orthol_tri(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	orthol_tri[i][j] = results[2];
	orthol_tri[j][i] = results[2];
	
	return 0;
}

int output_orthol_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *orthol_tri_size = create_ivector(number);
	orthol_tri_size[0] = terms_late;
	
	int orthol_tri_total_size = orthol_tri_size[0]*orthol_tri_size[0];
	array_write(&orthol_tri_total_size, orthol_tri_data_file, &orthol_tri[0][0]);
	
	return 0;

}

int read_orthol_tri(char *directory){
	
	int number = 1;
	int i,j,k;
	
	int *orthol_tri_size = create_ivector(number);
	orthol_tri_size[0] = terms_late;
	
	double **orthol1_tri = (double **)create_array(orthol_tri_size[0],orthol_tri_size[0]);
	int orthol_tri_total_size = orthol_tri_size[0]*orthol_tri_size[0];
	array_read(&orthol_tri_total_size, orthol_tri_data_file, &orthol1_tri[0][0]);
	
	int size = terms_late;
	orthol_tri = (double **)create_array(size,size);
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			orthol_tri[i][j] = orthol1_tri[i][j];
// 			printf("%e\t",orthol[i][j]);
		}
// 		printf("\n");
	}
	
	return 0;
	
}


int create_lambda(){
	
	int n = terms_prim;
	lambda = (double **)create_array(n,n);
	
	return 0;
}

double get_lambda(int i,int j){
	double value = lambda[i][j];
	return value;
}
	
int update_lambda(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	lambda[i][j] = results[2];
	lambda[j][i] = results[2];
	
	return 0;
}

int output_lambda(char *directory){
	
	int number = 1;
	int i,j,k;
	
	int *lambda_size = create_ivector(number);
	lambda_size[0] = terms_prim;
	
	int lambda_total_size = lambda_size[0]*lambda_size[0];
	array_write(&lambda_total_size, lambda_data_file, &lambda[0][0]);
	
	return 0;
	
}

int read_lambda(){
	
	int number = 1;
	int i,j,k;
	
	int *lambda_size = create_ivector(number);
	lambda_size[0] = terms_prim;
	
	double **lambda1 = (double **)create_array(lambda_size[0],lambda_size[0]);
	int lambda_total_size = lambda_size[0]*lambda_size[0];
	array_read(&lambda_total_size, lambda_data_file, &lambda1[0][0]);
	printf("%s\n", lambda_data_file);
	
	int size = terms_prim;
	lambda = (double **)create_array(size,size);
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			lambda[i][j] = lambda1[i][j];
		}
	}
	
	return 0;
	
}

int create_lambda_tri(){
	
	int n = terms_prim;
	lambda_tri = (double **)create_array(n,n);
	
	return 0;
}

double get_lambda_tri(int i,int j){
	double value = lambda_tri[i][j];
	return value;
}
	
int update_lambda_tri(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	lambda_tri[i][j] = results[2];
	lambda_tri[j][i] = results[2];
	
	return 0;
}

int output_lambda_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *lambda_tri_size = create_ivector(number);
	lambda_tri_size[0] = terms_prim;
	
	int lambda_tri_total_size = lambda_tri_size[0]*lambda_tri_size[0];
	array_write(&lambda_tri_total_size, lambda_tri_data_file, &lambda_tri[0][0]);
	
	return 0;
	
}

int read_lambda_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *lambda_tri_size = create_ivector(number);
	lambda_tri_size[0] = terms_prim;
	
	if (lambda_tri_size[0] < terms_prim){
		printf("Alpha max in parameter is greater than alpha max in loaded file. Exiting.\n");
		exit(1);
	}
	
	double **lambda1 = (double **)create_array(lambda_tri_size[0],lambda_tri_size[0]);
	int lambda_tri_total_size = lambda_tri_size[0]*lambda_tri_size[0];
	array_read(&lambda_tri_total_size, lambda_tri_data_file, &lambda1[0][0]);
	
	int size = terms_prim;
	lambda_tri = (double **)create_array(size,size);
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			lambda_tri[i][j] = lambda1[i][j];
		}
	}
	
	return 0;
	
}

int create_lambdal(){
	
	int n = terms_late;
	lambdal_TTT = (double **)create_array(n,n);
	
	if(do_polarisation==1){
		lambdal_TTE = (double **)create_array(n,n);
		lambdal_TEE = (double **)create_array(n,n);
		lambdal_EEE = (double **)create_array(n,n);
	}
	
	return 0;
}

double get_lambdal_TTT(int i,int j){
	double value = lambdal_TTT[i][j];
	return value;
}

double get_lambdal_TTE(int i,int j){
	double value = lambdal_TTE[i][j];
	return value;
}

double get_lambdal_TEE(int i,int j){
	double value = lambdal_TEE[i][j];
	return value;
}

double get_lambdal_EEE(int i,int j){
	double value = lambdal_EEE[i][j];
	return value;
}
	
int update_lambdal(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	lambdal_TTT[i][j] = results[2];
	if(do_polarisation==1){
		lambdal_TTE[i][j] = results[3];
		lambdal_TEE[i][j] = results[4];
		lambdal_EEE[i][j] = results[5];
	}
	
	return 0;
}

int output_lambdal(){
	
	int number = 1;
	int i,j,k;
	
	int lambdal_total_size = terms_late*terms_late;
	array_write(&lambdal_total_size, lambdal_TTT_data_file, &lambdal_TTT[0][0]);
	
	if(do_polarisation==1){
		array_write(&lambdal_total_size, lambdal_TTE_data_file, &lambdal_TTE[0][0]);
		array_write(&lambdal_total_size, lambdal_TEE_data_file, &lambdal_TEE[0][0]);
		array_write(&lambdal_total_size, lambdal_EEE_data_file, &lambdal_EEE[0][0]);
	}
	
	return 0;
	
}

int read_lambdal(char *directory){
	
	int number = 1;
	int i,j,k;
	
	int lambdal_total_size = terms_late*terms_late;
	lambdal_TTT = (double **)create_array(terms_late,terms_late);
	array_read(&lambdal_total_size, lambdal_TTT_data_file, &lambdal_TTT[0][0]);
	
	if(do_polarisation==1){
		lambdal_TTE = (double **)create_array(terms_late,terms_late);
		lambdal_TEE = (double **)create_array(terms_late,terms_late);
		lambdal_EEE = (double **)create_array(terms_late,terms_late);
		array_read(&lambdal_total_size, lambdal_TTE_data_file, &lambdal_TTE[0][0]);
		array_read(&lambdal_total_size, lambdal_TEE_data_file, &lambdal_TEE[0][0]);
		array_read(&lambdal_total_size, lambdal_EEE_data_file, &lambdal_EEE[0][0]);
	}
	
	return 0;
	
}

int create_lambdal_tri(){
	
	int n = terms_late;
	lambdal_tri = (double **)create_array(n,n);
	
	return 0;
}

double get_lambdal_tri(int i,int j){
	double value = lambdal_tri[i][j];
	return value;
}

int update_lambdal_tri(double *results){

	int i,j;
	
	i = results[0];
	j = results[1];
	lambdal_tri[i][j] = results[2];
// 	lambdal_tri[j][i] = results[2];
	
	return 0;
}

int output_lambdal_tri(char *directory){
	
	int number = 1;
	int i,j,k;
	
	int *lambdal_tri_size = create_ivector(number);
	lambdal_tri_size[0] = terms_late;
	
	int lambdal_tri_total_size = lambdal_tri_size[0]*lambdal_tri_size[0];
	array_write(&lambdal_tri_total_size, lambdal_tri_data_file, &lambdal_tri[0][0]);
	
	return 0;

}

int read_lambdal_tri(char *directory){
	
	int number = 1;
	int i,j,k;
	
	int *lambdal_tri_size = create_ivector(number);
	lambdal_tri_size[0] = terms_late;
	
	double **lambdal1_tri = (double **)create_array(lambdal_tri_size[0],lambdal_tri_size[0]);
	int lambdal_tri_total_size = lambdal_tri_size[0]*lambdal_tri_size[0];
	array_read(&lambdal_tri_total_size, lambdal_tri_data_file, &lambdal1_tri[0][0]);
	
	int size = terms_late;
	lambdal_tri = (double **)create_array(size,size);
	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			lambdal_tri[i][j] = lambdal1_tri[i][j];
		}
	}
	
	return 0;
	
}

int create_modesR(){

	int n = terms_late;
	modesR_TTT = (double *)create_vector(n);
	if(do_polarisation==1){
		modesR_TTE = (double *)create_vector(n);
		modesR_TEE = (double *)create_vector(n);
		modesR_EEE = (double *)create_vector(n);
	}
	return 0;
}

double get_modesR_TTT(int i){
	double value = modesR_TTT[i];
	return value;
}

double get_modesR_TTE(int i){
	double value = modesR_TTE[i];
	return value;
}

double get_modesR_TEE(int i){
	double value = modesR_TEE[i];
	return value;
}

double get_modesR_EEE(int i){
	double value = modesR_EEE[i];
	return value;
}
	
int update_modesR(double *results){

	int n;
	
	n = (int)results[0];
	
	modesR_TTT[n] = results[1];

	if(do_polarisation==1){
		modesR_TTE[n] = results[2];
		modesR_TEE[n] = results[3];
		modesR_EEE[n] = results[4];
	}
	
	return 0;
}

int create_modes(){

	int n = terms_late;
	modes_TTT = (double *)create_vector(n);
	if(do_polarisation==1){
		modes_TTE = (double *)create_vector(n);
		modes_TEE = (double *)create_vector(n);
		modes_EEE = (double *)create_vector(n);
	}
	
	return 0;
}

double get_modes_TTT(int i){
	double value = modes_TTT[i];
	return value;
}

double get_modes_TTE(int i){
	double value = modes_TTE[i];
	return value;
}

double get_modes_TEE(int i){
	double value = modes_TEE[i];
	return value;
}

double get_modes_EEE(int i){
	double value = modes_EEE[i];
	return value;
}
	
int update_modes(double *results){

	int n;
	
	n = (int)results[0];
	
	modes_TTT[n] = results[1];

	if(do_polarisation==1){
		modes_TTE[n] = results[2];
		modes_TEE[n] = results[3];
		modes_EEE[n] = results[4];
	}
	
	return 0;
}

int output_modes(int p1, int p2, int p3){

	char strmodel[10] = "";
	char strp1[10] = "";
	char strp2[10] = "";
	char strp3[10] = "";
	sprintf(strmodel,"_%d",model);
	if(slim.p3_end!=0){
		sprintf(strp3,"_%d",p3);
		if(slim.p2_end!=0){
			sprintf(strp2,"_%d",p2);
			if(slim.p1_end!=0){
				sprintf(strp1,"_%d",p1);
			}else{
				sprintf(strp1,"_%d",0);
			}
		}else{
			sprintf(strp2,"_%d",0);
			if(slim.p1_end!=0){
				sprintf(strp1,"_%d",p1);
			}else{
				sprintf(strp1,"_%d",0);
			}
		}
	}else{
		if(slim.p2_end!=0){
			sprintf(strp2,"_%d",p2);
			if(slim.p1_end!=0){
				sprintf(strp1,"_%d",p1);
			}else{
				sprintf(strp1,"_%d",0);
			}
		}else{
			if(slim.p1_end!=0){
				sprintf(strp1,"_%d",p1);
			}
		}
	}
	
	char modes_TTT_output_file[200] = "";
	strcat(modes_TTT_output_file, modes_TTT_data_file);
	strcat(modes_TTT_output_file, strmodel);
	strcat(modes_TTT_output_file, strp1);
	strcat(modes_TTT_output_file, strp2);
	strcat(modes_TTT_output_file, strp3);
	
	int modes_total_size = terms_late;
	vector_write(&modes_total_size, modes_TTT_output_file, &modes_TTT[0]);

	if(do_polarisation==1){
		
		char modes_TTE_output_file[200] = "";
		strcat(modes_TTE_output_file, modes_TTE_data_file);
		strcat(modes_TTE_output_file, strmodel);
		strcat(modes_TTE_output_file, strp1);
		strcat(modes_TTE_output_file, strp2);
		strcat(modes_TTE_output_file, strp3);
		vector_write(&modes_total_size, modes_TTE_output_file, &modes_TTE[0]);
		
		char modes_TEE_output_file[200] = "";
		strcat(modes_TEE_output_file, modes_TEE_data_file);
		strcat(modes_TEE_output_file, strmodel);
		strcat(modes_TEE_output_file, strp1);
		strcat(modes_TEE_output_file, strp2);
		strcat(modes_TEE_output_file, strp3);
		vector_write(&modes_total_size, modes_TEE_output_file, &modes_TEE[0]);
		
		char modes_EEE_output_file[200] = "";
		strcat(modes_EEE_output_file, modes_EEE_data_file);
		strcat(modes_EEE_output_file, strmodel);
		strcat(modes_EEE_output_file, strp1);
		strcat(modes_EEE_output_file, strp2);
		strcat(modes_EEE_output_file, strp3);
		vector_write(&modes_total_size, modes_EEE_output_file, &modes_EEE[0]);
			
		//printf("%s %s %s %s\n", modes_TTT_output_file, modes_TTE_output_file, modes_TEE_output_file, modes_EEE_output_file);

	}
	
	return 0;
	
}

int read_modes(int p1, int p2, int p3){

	char strmodel[10] = "";
	char strp1[10] = "";
	char strp2[10] = "";
	char strp3[10] = "";
	sprintf(strmodel,"_%d",model);
	if(slim.p1_end!=0)sprintf(strp1,"_%d",p1);
	if(slim.p2_end!=0)sprintf(strp2,"_%d",p2);
	if(slim.p3_end!=0)sprintf(strp3,"_%d",p3);
	
	char modes_TTT_output_file[200] = "";
	strcat(modes_TTT_output_file, modes_TTT_data_file);
	strcat(modes_TTT_output_file, strmodel);
	strcat(modes_TTT_output_file, strp1);
	strcat(modes_TTT_output_file, strp2);
	strcat(modes_TTT_output_file, strp3);
	
	int modes_total_size = terms_late;
	modes_TTT = (double *)create_vector(terms_late);
	vector_read(&modes_total_size, modes_TTT_output_file, &modes_TTT[0]);

	if(do_polarisation==1){
		
		char modes_TTE_output_file[200] = "";
		strcat(modes_TTE_output_file, modes_TTE_data_file);
		strcat(modes_TTE_output_file, strmodel);
		strcat(modes_TTE_output_file, strp1);
		strcat(modes_TTE_output_file, strp2);
		strcat(modes_TTE_output_file, strp3);
		modes_TTE = (double *)create_vector(terms_late);
		vector_read(&modes_total_size, modes_TTE_output_file, &modes_TTE[0]);
		
		char modes_TEE_output_file[200] = "";
		strcat(modes_TEE_output_file, modes_TEE_data_file);
		strcat(modes_TEE_output_file, strmodel);
		strcat(modes_TEE_output_file, strp1);
		strcat(modes_TEE_output_file, strp2);
		strcat(modes_TEE_output_file, strp3);
		modes_TEE = (double *)create_vector(terms_late);
		vector_read(&modes_total_size, modes_TEE_output_file, &modes_TEE[0]);
		
		char modes_EEE_output_file[200] = "";
		strcat(modes_EEE_output_file, modes_EEE_data_file);
		strcat(modes_EEE_output_file, strmodel);
		strcat(modes_EEE_output_file, strp1);
		strcat(modes_EEE_output_file, strp2);
		strcat(modes_EEE_output_file, strp3);
		modes_EEE = (double *)create_vector(terms_late);
		vector_read(&modes_total_size, modes_EEE_output_file, &modes_EEE[0]);
		
	}
	
	return 0;
	
}

int destroy_modes(){
	destroy_vector(modes_TTT);
	if(do_polarisation==1){
		destroy_vector(modes_TTE);
		destroy_vector(modes_TEE);
		destroy_vector(modes_EEE);
	}
	return 0;
}

int create_modesR_tri(){

	int n = terms_late;
	modesR_tri = (double *)create_vector(n);
	
	return 0;
}

double get_modesR_tri(int i){
	double value = modesR_tri[i];
	return value;
}
	
int update_modesR_tri(double *results){

	int n;
	
	n = (int)results[0];
	
	modesR_tri[n] = results[1];
	
	return 0;
}

int create_modes_tri(){

	modes_tri = (double *)create_vector(terms_late);
	
	return 0;
}

double get_modes_tri(int i){
	double value = modes_tri[i];
	return value;
}
	
int update_modes_tri(double *results){

	int n;
	
	n = (int)results[0];
	
	modes_tri[n] = results[1];
	
	return 0;
}

int output_modes_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *modes_tri_size = create_ivector(number);
	modes_tri_size[0] = terms_late;
	
	int modes_tri_total_size = modes_tri_size[0];
	vector_write(&modes_tri_total_size, modes_tri_data_file, &modes_tri[0]);
	
	return 0;
	
}

int read_modes_tri(){
	
	int number = 1;
	int i,j,k;
	
	int *modes_tri_size = create_ivector(number);
	modes_tri_size[0] = terms_late;
	
	modes_tri = (double *)create_vector(modes_tri_size[0]);
	int modes_tri_total_size = modes_tri_size[0];
	vector_read(&modes_tri_total_size, modes_tri_data_file, &modes_tri[0]);
	
	return 0;
	
}

int create_transfer(int l_size, int t_size){
	transfer = (double **)create_array(l_size,t_size);
	return 0;
}
	
int update_transfer(int size, double *results){

	int l,a,k;
	
	l = (int)results[0];
	
	for (k=1;k<size;k++){
		transfer[l][k] = results[k];
	}
	
	return 0;
}

double get_transfer(int i, int j){
	double value = transfer[i][j];
	return value;
}

int create_fisher_tri(int xsize){
	fisher_tri = (double **)create_array(xsize,xsize);
	return 0;
}
	
int update_fisher_tri(double *results){

	int i,j;
	
	i = (int)results[0];
	j = (int)results[1];
	fisher_tri[i][j] = results[2];
	fisher_tri[j][i] = results[2];
	
	return 0;
}

double get_fisher_tri(int i, int j){
	double value = fisher_tri[i][j];
	return value;
}

int tri(int i,int j,int k){
	
	if(i>j+k || j>k+i || k>i+j){
		return 0;
	} else {
		return 1;
	}
}

int quad(int i,int j,int k, int l){
	
	if(i>j+k+l || j>k+l+i || k>l+i+j || l>i+j+k){
		return 0;
	} else {
		return 1;
	}
}

double permsix(int i, int j, int k){
	double a;
	if(i==k&&i==j){
		a=1.0;
	}else if(i==j||j==k||i==k){
		a=3.0;
	}else{
		a=6.0;
	}
	return a;
}

void fourfig(int i, char* suffix){
	char text[5];
	suffix[0] = '\0';
	text[0] = '\0';
	
	sprintf(text, "%d", i);
	if(i<10){
		strcat(suffix, "000");
		strcat(suffix, text);
	}else if(i<100){
		strcat(suffix, "00");
		strcat(suffix, text);
	}else if(i<1000){
		strcat(suffix, "0");
		strcat(suffix, text);
	}else{
		strcat(suffix, text);
	}
}

