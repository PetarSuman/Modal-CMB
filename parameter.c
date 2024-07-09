#include "global.h"

// JB ICC tries to compile this as offload for no reason??
#ifndef __MIC__
//#include <iniparser.h>
#include <../iniparser/src/iniparser.h>
#endif

char alpha_data_file[MAXLEN];
char alpha_size_file[MAXLEN];
char bessel_data_file[MAXLEN];
char bessel_size_file[MAXLEN];
char transfer_T_data_file[MAXLEN];
char transfer_E_data_file[MAXLEN];
char transfer_size_file[MAXLEN];
char source_data_file[MAXLEN];
char source_size_file[MAXLEN];
char limber_data_file[MAXLEN];
char limber_size_file[MAXLEN];
char decompose_data_file[MAXLEN];
char decompose_size_file[MAXLEN];
char eigen_data_file[MAXLEN];
char eigen_tri_data_file[MAXLEN];
char ortho_data_file[MAXLEN];
char orthol_TTT_data_file[MAXLEN];
char orthol_TTE_data_file[MAXLEN];
char orthol_TEE_data_file[MAXLEN];
char orthol_EEE_data_file[MAXLEN];
char ortho_tri_data_file[MAXLEN];
char orthol_tri_data_file[MAXLEN];
char lambda_data_file[MAXLEN];
char lambdal_TTT_data_file[MAXLEN];
char lambdal_TTE_data_file[MAXLEN];
char lambdal_TEE_data_file[MAXLEN];
char lambdal_EEE_data_file[MAXLEN];
char lambda_tri_data_file[MAXLEN];
char lambdal_tri_data_file[MAXLEN];
char gamma_TTT_data_file[MAXLEN];
char gamma_TTE_data_file[MAXLEN];
char gamma_TEE_data_file[MAXLEN];
char gamma_EEE_data_file[MAXLEN];
char gamma_tri_data_file[MAXLEN];
char proj_T_data_file[MAXLEN];
char proj_E_data_file[MAXLEN];
char proj_T_size_file[MAXLEN];
char proj_E_size_file[MAXLEN];
char proj_tri_data_file[MAXLEN];
char proj_tri_size_file[MAXLEN];
char modes_TTT_data_file[MAXLEN];
char modes_TTE_data_file[MAXLEN];
char modes_TEE_data_file[MAXLEN];
char modes_EEE_data_file[MAXLEN];
char modes_tri_data_file[MAXLEN];
char flat_data_file[MAXLEN];
char flat_size_file[MAXLEN];
char l_data_file[MAXLEN];
char l_size_file[MAXLEN];
char bispectrum_file[MAXLEN];
char trispectrum_file[MAXLEN];
char interpolated_l_data_file[MAXLEN];
char interpolated_l_size_file[MAXLEN];
char interpolated_bispectrum_file[MAXLEN];
char pixw_file[MAXLEN];
char BN_TT_file[MAXLEN];
char BN_EE_file[MAXLEN];
char transfer_wgt_file[MAXLEN];
char cls_scalar_data_file[MAXLEN];
char cls_lensed_data_file[MAXLEN];
char restart_file[MAXLEN];
char bload_file[MAXLEN];
char edist_file[MAXLEN];
char edist_tri_file[MAXLEN];
char mdist_TTT_file[MAXLEN];
char mdist_TTE_file[MAXLEN];
char mdist_TEE_file[MAXLEN];
char mdist_EEE_file[MAXLEN];
char mdist_tri_file[MAXLEN];

double nscalar;
double deltaphi;
double kpivot;
int l_flat;
int bessel_max;
int bessel_points;
int multiplier;
int section_count;
double one_accuracy;
double tilt;
int initial_grid_size;
int depth;
double tri_accuracy;
int alpha_max;
int alpha_points;
int step;
int l1,l2,l3,l4;
int model;
int model2;

shape_limits slim;

// Flags
int do_polarisation;
int use_l_file;
char lini_file[MAXLEN];
int eflag_order_prim;
int eflag_order_late;
int eflag_T_lmax;
int eflag_T_lmin;
int eflag_E_lmax;
int eflag_E_lmin;
double fskyT;
double fskyE;
int rflag_do3D;
int clflag_uselens;
int tflag_bessel;
int tflag_transfer;
int tflag_flat;
int tflag_limber;
int bflag_lset;
int bflag_load;
int dflag_out;
int dflag_dmax;
int gflag_pca;
int oflag_load;

#ifndef __MIC__
int initilise(char *inifile){
	/*
	Read all the relevant parameters and files from the given directory
	*/
	dictionary *d;
	d = iniparser_load(inifile);
	if (d==NULL) {
		printf("cannot parse file [%s]", inifile);
		return -1 ;
	}
	
	char *string = malloc(MAXLEN*(sizeof(char)));
	char suffix1[5],suffix2[5],suffix3[5],suffix4[5],suffix5[5],suffix6[5],suffix7[5];
	char directory[MAXLEN];

	//	Load all relevant model parameters and survey-related files
	string = iniparser_getstring(d, "model:directory", "/nfs/st01/hpc-gr-epss/ps792/");
	strcpy(directory,string);
	
	eflag_T_lmax = iniparser_getint(d, "model:T_lmax", 1000);
	eflag_E_lmax = iniparser_getint(d, "model:E_lmax", 1000);
	eflag_T_lmin = iniparser_getint(d, "model:T_lmin", 2);
	eflag_E_lmin = iniparser_getint(d, "model:E_lmin", 2);
	
	fskyT = iniparser_getdouble(d, "model:fskyT", 1e0);
	fskyE = iniparser_getdouble(d, "model:fskyE", 1e0);
	
	nscalar = iniparser_getdouble(d, "model:nscalar", 9.60891e-1);
	deltaphi = iniparser_getdouble(d, "model:delta_phi", 1.5714e-8);
	kpivot = iniparser_getdouble(d, "model:kpivot", 5e-2);
	
	string = iniparser_getstring(d, "model:pixw_file", "/nfs/st01/hpc-gr-epss/ps792/PlanckV3/pixwgt_2048.txt");
	if(string[0] == '/'){
		strcpy(pixw_file,string);
	}else{
		pixw_file[0] = '\0';
		strcat(pixw_file, directory);
		strcat(pixw_file, string);
	}
	
	string = iniparser_getstring(d, "model:BN_TT_file", "/nfs/st01/hpc-gr-epss/ps792/DX12/BN_DX12_v3_smica_T.txt");
	if(string[0] == '/'){
		strcpy(BN_TT_file,string);
	}else{
		BN_TT_file[0] = '\0';
		strcat(BN_TT_file, directory);
		strcat(BN_TT_file, string);
	}
	
	string = iniparser_getstring(d, "model:BN_EE_file", "/nfs/st01/hpc-gr-epss/ps792/DX12/BN_DX12_v3_smica_E.txt");
	if(string[0] == '/'){
		strcpy(BN_EE_file,string);
	}else{
		BN_EE_file[0] = '\0';
		strcat(BN_EE_file, directory);
		strcat(BN_EE_file, string);
	}
	
	use_l_file = iniparser_getint(d, "model:use_l_file", 0);	
	string = iniparser_getstring(d, "model:lini_file", "lini.txt");
	if(string[0] == '/'){
		strcpy(lini_file,string);
	}else{
		lini_file[0] = '\0';
		strcat(lini_file, directory);
		strcat(lini_file, string);
	}
	 
	do_polarisation = iniparser_getint(d, "model:do_polarisation", 0);

	 
	slim.p1_bgn = iniparser_getint(d, "model:shape1_bgn", 0);
	slim.p1_stp = iniparser_getint(d, "model:shape1_stp", 1);
	slim.p1_end = iniparser_getint(d, "model:shape1_end", 0);
	slim.p2_bgn = iniparser_getint(d, "model:shape2_bgn", 0);
	slim.p2_stp = iniparser_getint(d, "model:shape2_stp", 1);
	slim.p2_end = iniparser_getint(d, "model:shape2_end", 0);
	slim.p3_bgn = iniparser_getint(d, "model:shape3_bgn", 0);
	slim.p3_stp = iniparser_getint(d, "model:shape3_stp", 1);
	slim.p3_end = iniparser_getint(d, "model:shape3_end", 0);

	//	Basis parameters
	alpha_max = iniparser_getint(d, "basis:terms", 101);
	eflag_order_prim = iniparser_getint(d, "basis:order_prim", 3);
	eflag_order_late = iniparser_getint(d, "basis:order_late", 3);
	rflag_do3D = iniparser_getint(d, "basis:do3D", 1);
	clflag_uselens = iniparser_getint(d, "basis:uselens", 1);
	alpha_points = iniparser_getint(d, "basis:alpha_points", 100);

	
	string = iniparser_getstring(d, "model:transfer_wgt_file", "wmap7_Tl.txt");
	if(string[0] == '/'){
		strcpy(transfer_wgt_file,string);
	}else{
		transfer_wgt_file[0] = '\0';
		strcat(transfer_wgt_file, directory);
		strcat(transfer_wgt_file, string);
	}
	
	string = iniparser_getstring(d, "model:cls_lensed_file", "cls_planck_DDX9_lensed.txt");
	if(string[0] == '/'){
		strcpy(cls_lensed_data_file,string);
	}else{
		cls_lensed_data_file[0] = '\0';
		strcat(cls_lensed_data_file, directory);
		strcat(cls_lensed_data_file, string);
	}
	
	string = iniparser_getstring(d, "model:cls_scalar_file", "cls_planck_DDX9_scalar.txt");
	if(string[0] == '/'){
		strcpy(cls_scalar_data_file,string);
	}else{
		cls_scalar_data_file[0] = '\0';
		strcat(cls_scalar_data_file, directory);
		strcat(cls_scalar_data_file, string);
	}
	
	suffix1[0] = '\0';
	sprintf(suffix1, "%d", eflag_T_lmax);
	suffix2[0] = '\0';
	sprintf(suffix2, "%d", eflag_E_lmax);
	suffix3[0] = '\0';
	sprintf(suffix3, "%d", alpha_max);
	suffix4[0] = '\0';
	sprintf(suffix4, "%d", eflag_order_prim);
	suffix5[0] = '\0';
	sprintf(suffix5, "%d", eflag_order_late);
	suffix6[0] = '\0';
	sprintf(suffix6, "%d", eflag_T_lmin);
	suffix7[0] = '\0';
	sprintf(suffix7, "%d", eflag_E_lmin);
	
	string = iniparser_getstring(d, "model:order_file_prim", "order_dist.txt");
	if(string[0] == '/'){
		strcpy(edist_file,string);
	}else{
		edist_file[0] = '\0';
		strcat(edist_file, directory);
		strcat(edist_file, string);
	}
	
	string = iniparser_getstring(d, "model:order_TTT_file_late", "order_dist.txt");
	if(string[0] == '/'){
		strcpy(mdist_TTT_file,string);
	}else{
		mdist_TTT_file[0] = '\0';
		strcat(mdist_TTT_file, directory);
		strcat(mdist_TTT_file, string);
	}
	
	string = iniparser_getstring(d, "model:order_TTE_file_late", "order_dist_XXY.txt");
	if(string[0] == '/'){
		strcpy(mdist_TTE_file,string);
	}else{
		mdist_TTE_file[0] = '\0';
		strcat(mdist_TTE_file, directory);
		strcat(mdist_TTE_file, string);
	}
	
	
	string = iniparser_getstring(d, "model:order_TEE_file_late", "order_dist_XXY.txt");
	if(string[0] == '/'){
		strcpy(mdist_TEE_file,string);
	}else{
		mdist_TEE_file[0] = '\0';
		strcat(mdist_TEE_file, directory);
		strcat(mdist_TEE_file, string);
	}
	
	
	string = iniparser_getstring(d, "model:order_EEE_file_late", "order_dist.txt");
	if(string[0] == '/'){
		strcpy(mdist_EEE_file,string);
	}else{
		mdist_EEE_file[0] = '\0';
		strcat(mdist_EEE_file, directory);
		strcat(mdist_EEE_file, string);
	}
	
	string = iniparser_getstring(d, "model:order_file_prim_tri", "order_dist_tri.txt");
	if(string[0] == '/'){
		strcpy(edist_tri_file,string);
	}else{
		edist_tri_file[0] = '\0';
		strcat(edist_tri_file, directory);
		strcat(edist_tri_file, string);
	}
	
	string = iniparser_getstring(d, "model:order_file_late_tri", "order_dist_tri.txt");
	if(string[0] == '/'){
		strcpy(mdist_tri_file,string);
	}else{
		mdist_tri_file[0] = '\0';
		strcat(mdist_tri_file, directory);
		strcat(mdist_tri_file, string);
	}
	
	string = iniparser_getstring(d, "basis:eigen_file", "eigen");
	if(string[0] == '/'){
		strcpy(eigen_data_file,string);
	}else{
		eigen_data_file[0] = '\0';
		strcat(eigen_data_file, directory);
		strcat(eigen_data_file, string);
		strcat(eigen_data_file, "_");
		strcat(eigen_data_file, suffix1);
		strcat(eigen_data_file, "_");
		strcat(eigen_data_file, suffix4);
		strcat(eigen_data_file, "_");
		strcat(eigen_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:eigen_file_tri", "eigen_tri");
	if(string[0] == '/'){
		strcpy(eigen_tri_data_file,string);
	}else{
		eigen_tri_data_file[0] = '\0';
		strcat(eigen_tri_data_file, directory);
		strcat(eigen_tri_data_file, string);
		strcat(eigen_tri_data_file, "_");
		strcat(eigen_tri_data_file, suffix1);
		strcat(eigen_tri_data_file, "_");
		strcat(eigen_tri_data_file, suffix4);
		strcat(eigen_tri_data_file, "_");
		strcat(eigen_tri_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:modes_file", "modes");
	if(string[0] == '/'){
		strcpy(modes_TTT_data_file,string);
	}else{
		modes_TTT_data_file[0] = '\0';
		strcat(modes_TTT_data_file, directory);
		strcat(modes_TTT_data_file, string);
		strcat(modes_TTT_data_file, "_TTT_");
		strcat(modes_TTT_data_file, suffix6);
		strcat(modes_TTT_data_file, "-");
		strcat(modes_TTT_data_file, suffix1);
		strcat(modes_TTT_data_file, "_");
		strcat(modes_TTT_data_file, suffix5);
		strcat(modes_TTT_data_file, "_");
		strcat(modes_TTT_data_file, suffix3);
		modes_TTE_data_file[0] = '\0';
		strcat(modes_TTE_data_file, directory);
		strcat(modes_TTE_data_file, string);
		strcat(modes_TTE_data_file, "_TTE_");
		strcat(modes_TTE_data_file, suffix6);
		strcat(modes_TTE_data_file, "-");
		strcat(modes_TTE_data_file, suffix1);
		strcat(modes_TTE_data_file, "_");
		strcat(modes_TTE_data_file, suffix7);
		strcat(modes_TTE_data_file, "-");
		strcat(modes_TTE_data_file, suffix2);
		strcat(modes_TTE_data_file, "_");
		strcat(modes_TTE_data_file, suffix5);
		strcat(modes_TTE_data_file, "_");
		strcat(modes_TTE_data_file, suffix3);
		modes_TEE_data_file[0] = '\0';
		strcat(modes_TEE_data_file, directory);
		strcat(modes_TEE_data_file, string);
		strcat(modes_TEE_data_file, "_TEE_");
		strcat(modes_TEE_data_file, suffix6);
		strcat(modes_TEE_data_file, "-");
		strcat(modes_TEE_data_file, suffix1);
		strcat(modes_TEE_data_file, "_");
		strcat(modes_TEE_data_file, suffix7);
		strcat(modes_TEE_data_file, "-");
		strcat(modes_TEE_data_file, suffix2);
		strcat(modes_TEE_data_file, "_");
		strcat(modes_TEE_data_file, suffix5);
		strcat(modes_TEE_data_file, "_");
		strcat(modes_TEE_data_file, suffix3);
		modes_EEE_data_file[0] = '\0';
		strcat(modes_EEE_data_file, directory);
		strcat(modes_EEE_data_file, string);
		strcat(modes_EEE_data_file, "_EEE_");
		strcat(modes_EEE_data_file, suffix7);
		strcat(modes_EEE_data_file, "-");
		strcat(modes_EEE_data_file, suffix2);
		strcat(modes_EEE_data_file, "_");
		strcat(modes_EEE_data_file, suffix5);
		strcat(modes_EEE_data_file, "_");
		strcat(modes_EEE_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:modes_file_tri", "modes_tri");
	if(string[0] == '/'){
		strcpy(modes_tri_data_file,string);
	}else{
		modes_tri_data_file[0] = '\0';
		strcat(modes_tri_data_file, directory);
		strcat(modes_tri_data_file, string);
		strcat(modes_tri_data_file, "_");
		strcat(modes_tri_data_file, suffix1);
		strcat(modes_tri_data_file, "_");
		strcat(modes_tri_data_file, suffix5);
		strcat(modes_tri_data_file, "_");
		strcat(modes_tri_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:ortho_file", "ortho");
	if(string[0] == '/'){
		strcpy(ortho_data_file,string);
	}else{
		ortho_data_file[0] = '\0';
		strcat(ortho_data_file, directory);
		strcat(ortho_data_file, string);
		strcat(ortho_data_file, "_");
		strcat(ortho_data_file, suffix1);
		strcat(ortho_data_file, "_");
		strcat(ortho_data_file, suffix4);
		strcat(ortho_data_file, "_");
		strcat(ortho_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:ortho_file_tri", "ortho_tri");
	if(string[0] == 'cd/'){
		strcpy(ortho_tri_data_file,string);
	}
	else{
		ortho_tri_data_file[0] = '\0';
		strcat(ortho_tri_data_file, directory);
		strcat(ortho_tri_data_file, string);
		strcat(ortho_tri_data_file, "_");
		strcat(ortho_tri_data_file, suffix1);
		strcat(ortho_tri_data_file, "_");
		strcat(ortho_tri_data_file, suffix4);
		strcat(ortho_tri_data_file, "_");
		strcat(ortho_tri_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:orthol_file", "orthol");
	if(string[0] == '/'){
		strcpy(orthol_TTT_data_file,string);
	}
	else{
		orthol_TTT_data_file[0] = '\0';
		strcat(orthol_TTT_data_file, directory);
		strcat(orthol_TTT_data_file, string);
		strcat(orthol_TTT_data_file, "_TTT_");
		strcat(orthol_TTT_data_file, suffix6);
		strcat(orthol_TTT_data_file, "-");
		strcat(orthol_TTT_data_file, suffix1);
		strcat(orthol_TTT_data_file, "_");
		strcat(orthol_TTT_data_file, suffix5);
		strcat(orthol_TTT_data_file, "_");
		strcat(orthol_TTT_data_file, suffix3);
		orthol_TTE_data_file[0] = '\0';
		strcat(orthol_TTE_data_file, directory);
		strcat(orthol_TTE_data_file, string);
		strcat(orthol_TTE_data_file, "_TTE_");
		strcat(orthol_TTE_data_file, suffix6);
		strcat(orthol_TTE_data_file, "-");
		strcat(orthol_TTE_data_file, suffix1);
		strcat(orthol_TTE_data_file, "_");
		strcat(orthol_TTE_data_file, suffix7);
		strcat(orthol_TTE_data_file, "-");
		strcat(orthol_TTE_data_file, suffix2);
		strcat(orthol_TTE_data_file, "_");
		strcat(orthol_TTE_data_file, suffix5);
		strcat(orthol_TTE_data_file, "_");
		strcat(orthol_TTE_data_file, suffix3);
		orthol_TEE_data_file[0] = '\0';
		strcat(orthol_TEE_data_file, directory);
		strcat(orthol_TEE_data_file, string);
		strcat(orthol_TEE_data_file, "_TEE_");
		strcat(orthol_TEE_data_file, suffix6);
		strcat(orthol_TEE_data_file, "-");
		strcat(orthol_TEE_data_file, suffix1);
		strcat(orthol_TEE_data_file, "_");
		strcat(orthol_TEE_data_file, suffix7);
		strcat(orthol_TEE_data_file, "-");
		strcat(orthol_TEE_data_file, suffix2);
		strcat(orthol_TEE_data_file, "_");
		strcat(orthol_TEE_data_file, suffix5);
		strcat(orthol_TEE_data_file, "_");
		strcat(orthol_TEE_data_file, suffix3);
		orthol_EEE_data_file[0] = '\0';
		strcat(orthol_EEE_data_file, directory);
		strcat(orthol_EEE_data_file, string);
		strcat(orthol_EEE_data_file, "_EEE_");
		strcat(orthol_EEE_data_file, suffix7);
		strcat(orthol_EEE_data_file, "-");
		strcat(orthol_EEE_data_file, suffix2);
		strcat(orthol_EEE_data_file, "_");
		strcat(orthol_EEE_data_file, suffix5);
		strcat(orthol_EEE_data_file, "_");
		strcat(orthol_EEE_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:orthol_file_tri", "orthol_tri");
	if(string[0] == '/'){
		strcpy(orthol_tri_data_file,string);
	}else{
		orthol_tri_data_file[0] = '\0';
		strcat(orthol_tri_data_file, directory);
		strcat(orthol_tri_data_file, string);
		strcat(orthol_tri_data_file, "_");
		strcat(orthol_tri_data_file, suffix1);
		strcat(orthol_tri_data_file, "_");
		strcat(orthol_tri_data_file, suffix5);
		strcat(orthol_tri_data_file, "_");
		strcat(orthol_tri_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:lambda_file", "lambda");
	if(string[0] == '/'){
		strcpy(lambda_data_file,string);
	}else{
		lambda_data_file[0] = '\0';
		strcat(lambda_data_file, directory);
		strcat(lambda_data_file, string);
		strcat(lambda_data_file, "_");
		strcat(lambda_data_file, suffix1);
		strcat(lambda_data_file, "_");
		strcat(lambda_data_file, suffix4);
		strcat(lambda_data_file, "_");
		strcat(lambda_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:lambda_file_tri", "lambda_tri");
	if(string[0] == '/'){
		strcpy(lambda_tri_data_file,string);
	}else{
		lambda_tri_data_file[0] = '\0';
		strcat(lambda_tri_data_file, directory);
		strcat(lambda_tri_data_file, string);
		strcat(lambda_tri_data_file, "_");
		strcat(lambda_tri_data_file, suffix1);
		strcat(lambda_tri_data_file, "_");
		strcat(lambda_tri_data_file, suffix4);
		strcat(lambda_tri_data_file, "_");
		strcat(lambda_tri_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:lambdal_file", "lambdal");
	if(string[0] == '/'){
		strcpy(lambdal_TTT_data_file,string);
	}else{
		lambdal_TTT_data_file[0] = '\0';
		strcat(lambdal_TTT_data_file, directory);
		strcat(lambdal_TTT_data_file, string);
		strcat(lambdal_TTT_data_file, "_TTT_");
		strcat(lambdal_TTT_data_file, suffix6);
		strcat(lambdal_TTT_data_file, "-");
		strcat(lambdal_TTT_data_file, suffix1);
		strcat(lambdal_TTT_data_file, "_");
		strcat(lambdal_TTT_data_file, suffix5);
		strcat(lambdal_TTT_data_file, "_");
		strcat(lambdal_TTT_data_file, suffix3);
		lambdal_TTE_data_file[0] = '\0';
		strcat(lambdal_TTE_data_file, directory);
		strcat(lambdal_TTE_data_file, string);
		strcat(lambdal_TTE_data_file, "_TTE_");
		strcat(lambdal_TTE_data_file, suffix6);
		strcat(lambdal_TTE_data_file, "-");
		strcat(lambdal_TTE_data_file, suffix1);
		strcat(lambdal_TTE_data_file, "_");
		strcat(lambdal_TTE_data_file, suffix7);
		strcat(lambdal_TTE_data_file, "-");
		strcat(lambdal_TTE_data_file, suffix2);
		strcat(lambdal_TTE_data_file, "_");
		strcat(lambdal_TTE_data_file, suffix5);
		strcat(lambdal_TTE_data_file, "_");
		strcat(lambdal_TTE_data_file, suffix3);
		lambdal_TEE_data_file[0] = '\0';
		strcat(lambdal_TEE_data_file, directory);
		strcat(lambdal_TEE_data_file, string);
		strcat(lambdal_TEE_data_file, "_TEE_");
		strcat(lambdal_TEE_data_file, suffix6);
		strcat(lambdal_TEE_data_file, "-");
		strcat(lambdal_TEE_data_file, suffix1);
		strcat(lambdal_TEE_data_file, "_");
		strcat(lambdal_TEE_data_file, suffix7);
		strcat(lambdal_TEE_data_file, "-");
		strcat(lambdal_TEE_data_file, suffix2);
		strcat(lambdal_TEE_data_file, "_");
		strcat(lambdal_TEE_data_file, suffix5);
		strcat(lambdal_TEE_data_file, "_");
		strcat(lambdal_TEE_data_file, suffix3);
		lambdal_EEE_data_file[0] = '\0';
		strcat(lambdal_EEE_data_file, directory);
		strcat(lambdal_EEE_data_file, string);
		strcat(lambdal_EEE_data_file, "_EEE_");
		strcat(lambdal_EEE_data_file, suffix7);
		strcat(lambdal_EEE_data_file, "-");
		strcat(lambdal_EEE_data_file, suffix2);
		strcat(lambdal_EEE_data_file, "_");
		strcat(lambdal_EEE_data_file, suffix5);
		strcat(lambdal_EEE_data_file, "_");
		strcat(lambdal_EEE_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:lambdal_file_tri", "lambdal_tri");
	if(string[0] == '/'){
		strcpy(lambdal_tri_data_file,string);
	}else{
		lambdal_tri_data_file[0] = '\0';
		strcat(lambdal_tri_data_file, directory);
		strcat(lambdal_tri_data_file, string);
		strcat(lambdal_tri_data_file, "_");
		strcat(lambdal_tri_data_file, suffix1);
		strcat(lambdal_tri_data_file, "_");
		strcat(lambdal_tri_data_file, suffix5);
		strcat(lambdal_tri_data_file, "_");
		strcat(lambdal_tri_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:gamma_file", "gamma");
	if(string[0] == '/'){
		strcpy(gamma_TTT_data_file,string);
		strcat(gamma_TTT_data_file, "_TTT");
		strcpy(gamma_TTE_data_file,string);
		strcat(gamma_TTE_data_file, "_TTE");
		strcpy(gamma_TEE_data_file,string);
		strcat(gamma_TEE_data_file, "_TEE");
		strcpy(gamma_EEE_data_file,string);
		strcat(gamma_EEE_data_file, "_EEE");
	}else{
		gamma_TTT_data_file[0] = '\0';
		strcat(gamma_TTT_data_file, directory);
		strcat(gamma_TTT_data_file, string);
		strcat(gamma_TTT_data_file, "_TTT_");
		strcat(gamma_TTT_data_file, suffix6);
		strcat(gamma_TTT_data_file, "-");
		strcat(gamma_TTT_data_file, suffix1);
		strcat(gamma_TTT_data_file, "_");
		strcat(gamma_TTT_data_file, suffix4);
		strcat(gamma_TTT_data_file, "_");
		strcat(gamma_TTT_data_file, suffix5);
		strcat(gamma_TTT_data_file, "_");
		strcat(gamma_TTT_data_file, suffix3);
		gamma_TTE_data_file[0] = '\0';
		strcat(gamma_TTE_data_file, directory);
		strcat(gamma_TTE_data_file, string);
		strcat(gamma_TTE_data_file, "_TTE_");
		strcat(gamma_TTE_data_file, suffix6);
		strcat(gamma_TTE_data_file, "-");
		strcat(gamma_TTE_data_file, suffix1);
		strcat(gamma_TTE_data_file, "_");
		strcat(gamma_TTE_data_file, suffix7);
		strcat(gamma_TTE_data_file, "-");
		strcat(gamma_TTE_data_file, suffix2);
		strcat(gamma_TTE_data_file, "_");
		strcat(gamma_TTE_data_file, suffix4);
		strcat(gamma_TTE_data_file, "_");
		strcat(gamma_TTE_data_file, suffix5);
		strcat(gamma_TTE_data_file, "_");
		strcat(gamma_TTE_data_file, suffix3);
		gamma_TEE_data_file[0] = '\0';
		strcat(gamma_TEE_data_file, directory);
		strcat(gamma_TEE_data_file, string);
		strcat(gamma_TEE_data_file, "_TEE_");
		strcat(gamma_TEE_data_file, suffix6);
		strcat(gamma_TEE_data_file, "-");
		strcat(gamma_TEE_data_file, suffix1);
		strcat(gamma_TEE_data_file, "_");
		strcat(gamma_TEE_data_file, suffix7);
		strcat(gamma_TEE_data_file, "-");
		strcat(gamma_TEE_data_file, suffix2);
		strcat(gamma_TEE_data_file, "_");
		strcat(gamma_TEE_data_file, suffix4);
		strcat(gamma_TEE_data_file, "_");
		strcat(gamma_TEE_data_file, suffix5);
		strcat(gamma_TEE_data_file, "_");
		strcat(gamma_TEE_data_file, suffix3);
		gamma_EEE_data_file[0] = '\0';
		strcat(gamma_EEE_data_file, directory);
		strcat(gamma_EEE_data_file, string);
		strcat(gamma_EEE_data_file, "_EEE_");
		strcat(gamma_EEE_data_file, suffix7);
		strcat(gamma_EEE_data_file, "-");
		strcat(gamma_EEE_data_file, suffix2);
		strcat(gamma_EEE_data_file, "_");
		strcat(gamma_EEE_data_file, suffix4);
		strcat(gamma_EEE_data_file, "_");
		strcat(gamma_EEE_data_file, suffix5);
		strcat(gamma_EEE_data_file, "_");
		strcat(gamma_EEE_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:gamma_tri_file", "gamma_tri");
	if(string[0] == '/'){
		strcpy(gamma_tri_data_file,string);
	}else{
		gamma_tri_data_file[0] = '\0';
		strcat(gamma_tri_data_file, directory);
		strcat(gamma_tri_data_file, string);
		strcat(gamma_tri_data_file, "_");
		strcat(gamma_tri_data_file, suffix1);
		strcat(gamma_tri_data_file, "_");
		strcat(gamma_tri_data_file, suffix4);
		strcat(gamma_tri_data_file, "_");
		strcat(gamma_tri_data_file, suffix5);
		strcat(gamma_tri_data_file, "_");
		strcat(gamma_tri_data_file, suffix3);
	}
	
	string = iniparser_getstring(d, "basis:proj_file", "proj");
	if(string[0] == '/'){
		strcpy(proj_T_data_file,string);
		strcat(proj_T_data_file, "_T");
		strcpy(proj_E_data_file,string);
		strcat(proj_E_data_file, "_E");
		strcat(proj_T_size_file, string);
		strcat(proj_T_size_file, "_T_size");
		strcat(proj_E_size_file, string);
		strcat(proj_E_size_file, "_E_size");
	}else{
		proj_T_data_file[0] = '\0';
		strcat(proj_T_data_file, directory);
		strcat(proj_T_data_file, string);
		strcat(proj_T_data_file, "_T_");
		strcat(proj_T_data_file, suffix1);
		strcat(proj_T_data_file, "_");
		strcat(proj_T_data_file, suffix4);
		proj_E_data_file[0] = '\0';
		strcat(proj_E_data_file, directory);
		strcat(proj_E_data_file, string);
		strcat(proj_E_data_file, "_E_");
		strcat(proj_E_data_file, suffix2);
		strcat(proj_E_data_file, "_");
		strcat(proj_E_data_file, suffix4);
		proj_T_size_file[0] = '\0';
		strcat(proj_T_size_file, directory);
		strcat(proj_T_size_file, string);
		strcat(proj_T_size_file, "_T_");
		strcat(proj_T_size_file, suffix1);
		strcat(proj_T_size_file, "_size");
		proj_E_size_file[0] = '\0';
		strcat(proj_E_size_file, directory);
		strcat(proj_E_size_file, string);
		strcat(proj_E_size_file, "_E_");
		strcat(proj_E_size_file, suffix2);
		strcat(proj_E_size_file, "_size");
	}
	
	string = iniparser_getstring(d, "basis:proj_tri_file", "proj_tri");
	if(string[0] == '/'){
		strcpy(proj_tri_data_file,string);
		proj_tri_size_file[0] = '\0';
		strcat(proj_tri_size_file, string);
		strcat(proj_tri_size_file, "_size");
	}else{
		proj_tri_data_file[0] = '\0';
		strcat(proj_tri_data_file, directory);
		strcat(proj_tri_data_file, string);
		strcat(proj_tri_data_file, "_");
		strcat(proj_tri_data_file, suffix1);
		strcat(proj_tri_data_file, "_");
		strcat(proj_tri_data_file, suffix4);
		proj_tri_size_file[0] = '\0';
		strcat(proj_tri_size_file, directory);
		strcat(proj_tri_size_file, string);
		strcat(proj_tri_size_file, "_");
		strcat(proj_tri_size_file, suffix1);
		strcat(proj_tri_size_file, "_size");
	}
	
	//	Transfer parameters
	tflag_bessel = iniparser_getint(d, "transfer:flag_bessel", 0);
	tflag_transfer = iniparser_getint(d, "transfer:flag_transfer", 0);
	tflag_flat = iniparser_getint(d, "transfer:flag_flat", 0);
	tflag_limber = iniparser_getint(d, "transfer:flag_limber", 0);
	
	bessel_max = iniparser_getint(d, "transfer:bessel_max", 10000);
	bessel_points = iniparser_getint(d, "transfer:bessel_points", 8000);
	
	string = iniparser_getstring(d, "transfer:source_file", "source");
	if(string[0] == '/'){
		strcpy(source_data_file,string);
		source_size_file[0] = '\0';
		strcat(source_size_file, string);
		strcat(source_size_file, "_size");
	}else{
		source_data_file[0] = '\0';
		strcat(source_data_file, directory);
		strcat(source_data_file, string);
		source_size_file[0] = '\0';
		strcat(source_size_file, directory);
		strcat(source_size_file, string);
		strcat(source_size_file, "_size");
	}
	
	string = iniparser_getstring(d, "transfer:bessel_file", "bessel");
	if(string[0] == '/'){
		strcpy(bessel_data_file,string);
		bessel_size_file[0] = '\0';
		strcat(bessel_size_file, string);
		strcat(bessel_size_file, "_size");
	}else{
		bessel_data_file[0] = '\0';
		strcat(bessel_data_file, directory);
		strcat(bessel_data_file, string);
		bessel_size_file[0] = '\0';
		strcat(bessel_size_file, directory);
		strcat(bessel_size_file, string);
		strcat(bessel_size_file, "_size");
	}
	
	string = iniparser_getstring(d, "transfer:transfer_file", "transfer");
	if(string[0] == '/'){
		strcpy(transfer_T_data_file,string);
		strcat(transfer_T_data_file, "_T");
		transfer_size_file[0] = '\0';
		strcat(transfer_size_file, string);
		strcat(transfer_size_file, "_size");
	}else{
		transfer_T_data_file[0] = '\0';
		strcat(transfer_T_data_file, directory);
		strcat(transfer_T_data_file, string);
		strcat(transfer_T_data_file, "_T");
		transfer_size_file[0] = '\0';
		strcat(transfer_size_file, directory);
		strcat(transfer_size_file, string);
		strcat(transfer_size_file, "_size");
	}
	if(do_polarisation==1){
		if(string[0] == '/'){
			strcpy(transfer_E_data_file,string);
			strcat(transfer_E_data_file, "_E");
		}else{
			transfer_E_data_file[0] = '\0';
			strcat(transfer_E_data_file, directory);
			strcat(transfer_E_data_file, string);
			strcat(transfer_E_data_file, "_E");
		}
	}
	
	string = iniparser_getstring(d, "transfer:flat_file", "flat");
	if(string[0] == '/'){
		strcpy(flat_data_file,string);
		flat_size_file[0] = '\0';
		strcat(flat_size_file, string);
		strcat(flat_size_file, "_size");
	}else{
		flat_data_file[0] = '\0';
		strcat(flat_data_file, directory);
		strcat(flat_data_file, string);
		flat_size_file[0] = '\0';
		strcat(flat_size_file, directory);
		strcat(flat_size_file, string);
		strcat(flat_size_file, "_size");
	}
	
	string = iniparser_getstring(d, "transfer:limber_file", "limber");
	if(string[0] == '/'){
		strcpy(limber_data_file,string);
		limber_size_file[0] = '\0';
		strcat(limber_size_file, string);
		strcat(limber_size_file, "_size");
	}else{
		limber_data_file[0] = '\0';
		strcat(limber_data_file, directory);
		strcat(limber_data_file, string);
		limber_size_file[0] = '\0';
		strcat(limber_size_file, directory);
		strcat(limber_size_file, string);
		strcat(limber_size_file, "_size");
	}
	
	//	Bispectrum
	bflag_lset = iniparser_getint(d, "bispectrum:flag_lset", 0);
	bflag_load = iniparser_getint(d, "bispectrum:flag_load", 0);
	multiplier = iniparser_getint(d, "bispectrum:multiplier", 20);
	section_count = iniparser_getint(d, "bispectrum:section_count", 10000);
	one_accuracy = iniparser_getdouble(d, "bispectrum:one_accuracy", 1e-6);
	initial_grid_size = iniparser_getint(d, "bispectrum:initial_grid_size", 100);
	depth = iniparser_getint(d, "bispectrum:depth", 3);
	tri_accuracy = iniparser_getdouble(d, "bispectrum:tri_accuracy", 5e-2);
	l_flat = iniparser_getint(d, "bispectrum:l_flat", 150);
	
	string = iniparser_getstring(d, "bispectrum:l_file", "l_values");
	if(string[0] == '/'){
		strcpy(l_data_file,string);
		l_size_file[0] = '\0';
		strcat(l_size_file, string);
		strcat(l_size_file, "_size");
	}else{
		l_data_file[0] = '\0';
		strcat(l_data_file, directory);
		strcat(l_data_file, string);
		l_size_file[0] = '\0';
		strcat(l_size_file, directory);
		strcat(l_size_file, string);
		strcat(l_size_file, "_size");
	}
	
	string = iniparser_getstring(d, "bispectrum:l_int_file", "l_values_int");
	if(string[0] == '/'){
		strcpy(interpolated_l_data_file,string);
		interpolated_l_size_file[0] = '\0';
		strcat(interpolated_l_size_file, string);
		strcat(interpolated_l_size_file, "_size");
	}else{
		interpolated_l_data_file[0] = '\0';
		strcat(interpolated_l_data_file, directory);
		strcat(interpolated_l_data_file, string);
		interpolated_l_size_file[0] = '\0';
		strcat(interpolated_l_size_file, directory);
		strcat(interpolated_l_size_file, string);
		strcat(interpolated_l_size_file, "_size");
	}
	
	string = iniparser_getstring(d, "bispectrum:bispectrum_file", "bispectrum");
	if(string[0] == '/'){
		strcpy(bispectrum_file,string);
	}else{
		bispectrum_file[0] = '\0';
		strcat(bispectrum_file, directory);
		strcat(bispectrum_file, string);
	}
	
	string = iniparser_getstring(d, "bispectrum:bispectrum_int_file", "bispectrum_int");
	if(string[0] == '/'){
		strcpy(interpolated_bispectrum_file,string);
	}else{
		interpolated_bispectrum_file[0] = '\0';
		strcat(interpolated_bispectrum_file, directory);
		strcat(interpolated_bispectrum_file, string);
	}
	
	string = iniparser_getstring(d, "bispectrum:trispectrum_file", "trispectrum");
	if(string[0] == '/'){
		strcpy(trispectrum_file,string);
	}else{
		trispectrum_file[0] = '\0';
		strcat(trispectrum_file, directory);
		strcat(trispectrum_file, string);
	}
	
	string = iniparser_getstring(d, "bispectrum:bispectrum_file", "bispectrum_load");
	if(string[0] == '/'){
		strcpy(bload_file,string);
	}else{
		bload_file[0] = '\0';
		strcat(bload_file, directory);
		strcat(bload_file, string);
	}
	
	//	Old parameters
	step = iniparser_getint(d, "old:step", 20);
	dflag_out = iniparser_getint(d, "old:dflag_out", 1);
	dflag_dmax = iniparser_getint(d, "old:dflag_dmax", 10);
	gflag_pca = iniparser_getint(d, "old:gflag_pca", 1);
	tilt = iniparser_getdouble(d, "old:tilt", 0e0);
	
	string = iniparser_getstring(d, "old:restart_file", "restart");
	if(string[0] == '/'){
		strcpy(restart_file,string);
	}else{
		restart_file[0] = '\0';
		strcat(restart_file, directory);
		strcat(restart_file, string);
	}
	
	string = iniparser_getstring(d, "old:alpha_file", "alpha");
	if(string[0] == '/'){
		strcpy(alpha_data_file,string);
		alpha_size_file[0] = '\0';
		strcat(alpha_size_file, string);
		strcat(alpha_size_file, "_size");
	}else{
		alpha_data_file[0] = '\0';
		strcat(alpha_data_file, directory);
		strcat(alpha_data_file, string);
		alpha_size_file[0] = '\0';
		strcat(alpha_size_file, directory);
		strcat(alpha_size_file, string);
		strcat(alpha_size_file, "_size");
	}
	
	string = iniparser_getstring(d, "old:decompose_file", "decompose");
	if(string[0] == '/'){
		strcpy(decompose_data_file,string);
		decompose_size_file[0] = '\0';
		strcat(decompose_size_file, string);
		strcat(decompose_size_file, "_size");
	}else{
		decompose_data_file[0] = '\0';
		strcat(decompose_data_file, directory);
		strcat(decompose_data_file, string);
		decompose_size_file[0] = '\0';
		strcat(decompose_size_file, directory);
		strcat(decompose_size_file, string);
		strcat(decompose_size_file, "_size");
	}
	
	return 0;
}
#endif
