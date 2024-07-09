
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "global.h"

static double kmax;
static double xmax;
static double f_max;
static double ftau;
static double tau0;
static double tmax;
static int tlsize;
static int flsize;
static int blsize;
static int lmlsize;
static int ksize;
static int xsize;
static int fsize;
static int tausize;
static double **transfer_T;
static double **transfer_E;
static double **transfer_flat;
static double **bessel;
static double **limber;

int load_txt_int(char* filename, int columns, int* values, int* size){
	FILE* fp;
	char line[MAXLEN];
	int n;
	int i = -1;
	// int j = -1;
	char** cptr = (char**) malloc(sizeof(char*));

	if ( !(fp = fopen(filename, "r")) ) {
		perror(filename);
		return 0; 
	}
	// printf("beginning\n");
	
	while (fgets(line, MAXLEN, fp)) {
		if (*line == '#') continue;
		// j++;
		// printf("Line j=%d\n",j);
		
		values[++i] = (int)strtod( line, cptr);

		// printf("here1 i=%d\n",i);
		for(n=1;n<columns;n++){
			if (cptr!=NULL ){
				values[++i] = (int)strtod( *cptr, cptr);
				// printf("here2 i=%d\n",i);
			}
		}
		
		if ( i>(columns*(MAXLINES-1)) ) break;
	}
	*size = i/columns+1;
	if(*size<0)*size=0;
	return 1;
}

int load_txt_dbl(char* filename, int columns, double* values, int* size){
	FILE* fp;
	char line[MAXLEN];
	int n;
	long int i = -1;
	// long int j = -1;
	char** cptr = (char**) malloc(sizeof(char*));

	if ( !(fp = fopen(filename, "r")) ) {
		perror(filename);
		return 0;
	}
	while (fgets(line, MAXLEN, fp)) {
		if (*line == '#') continue;
		
		values[++i] = strtod( line, cptr);
		for(n=1;n<columns;n++){
			if (cptr!=NULL ){
				values[++i] = strtod( *cptr, cptr);
			}
			//printf("VALUES[%d] = %e\n",i, values[i]);
		}
		
		if ( i>(columns*(MAXLINES-1)) ) break;
	}
	*size = (i/columns)+1;
	if(*size<0)*size=0;
// 	*size = i;
	return 1;
}

int load_txt_dbl_long(char* filename, int columns, double* values, long int *size){
	FILE* fp;
	char line[MAXLEN];
	int n;
	long int i = -1;
	char** cptr = (char**) malloc(sizeof(char*));

	if ( !(fp = fopen(filename, "r")) ) {
		perror(filename);
		return 0;
	}
	
	while (fgets(line, MAXLEN, fp)) {
		if (*line == '#') continue;
		
		values[++i] = strtod( line, cptr);

		for(n=1;n<columns;n++){
			if (cptr!=NULL ) values[++i] = strtod( *cptr, cptr);
		}
		
		if ( i>(columns*(MAXLINES-1)) ) break;
	}
	*size = i/columns+1;
// 	*size = i;
	return 1;
}
// 
int load_one(char* filename, int* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    values[++i] = (int)d1;
    if ( i>(MAXLINES-1) ) break;
  }
  *size = i+1;

  return 1;
}

int load_one_double(char* filename, double* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    values[++i] = d1;
    if ( i>(MAXLINES-1) ) break;
  }
  *size = i+1;

  return 1;
}

int load_one_double_long(char* filename, double* values, long int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1;
  long int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    values[++i] = d1;
    if ( i>(MAXLINES-1) ) break;
  }
  *size = i+1;

  return 1;
}

int load_two(char* filename, double* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = d2 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
  
    if (*line == '#') continue;	
    d1 = d2 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    values[++i] = d1;
    values[++i] = d2;
    if ( i>(2*MAXLINES-2) ) break;
  }
  *size = i/2 + 1;
  
  return 1;
}
int load_three_int(char* filename, int* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2,d3;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d3 = strtod( *cptr, cptr);
    values[++i] = (int)d1;
    values[++i] = (int)d2;
    values[++i] = (int)d3;
    if ( i>(3*MAXLINES-3) ) break;
  }
  *size = i/3+1;

  return 1;
}

int load_three(char* filename, double* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2,d3;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d3 = strtod( *cptr, cptr);
    values[++i] = d1;
    values[++i] = d2;
    values[++i] = d3;
    if ( i>(3*MAXLINES-3) ) break;
  }
  *size = i/3+1;

  return 1;
}

int load_four_int(char* filename, int* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2,d3,d4;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d3 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d4 = strtod( *cptr, cptr);
    values[++i] = (int)d1;
    values[++i] = (int)d2;
    values[++i] = (int)d3;
    values[++i] = (int)d4;
    if ( i>(4*MAXLINES-4) ) break;
  }
  *size = i/4+1;

  return 1;
}
int load_four(char* filename, double* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2,d3,d4;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d3 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d4 = strtod( *cptr, cptr);
    values[++i] = d1;
    values[++i] = d2;
    values[++i] = d3;
    values[++i] = d4;
    if ( i>(4*MAXLINES-4) ) break;
  }
  *size = i/4+1;

  return 1;
}
int load_five_int(char* filename, int* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2,d3,d4,d5;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d3 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d4 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d5 = strtod( *cptr, cptr);
    values[++i] = (int)d1;
    values[++i] = (int)d2;
    values[++i] = (int)d3;
    values[++i] = (int)d4;
    values[++i] = (int)d5;
    if ( i>(5*MAXLINES-5) ) break;
  }
  *size = i/5+1;

  return 1;
}
int load_six(char* filename, double* values, int* size){
  FILE* fp;
  char line[MAXLEN];
  double d1,d2,d3,d4,d5,d6;
  int i = -1;
  char** cptr = (char**) malloc(sizeof(char*));

  d1 = 0;

  if ( !(fp = fopen(filename, "r")) ) {
    perror(filename);
    return 0;
  }

  while (fgets(line, MAXLEN, fp)) {
    if (*line == '#') continue;	
    d1 = 0;
    d1 = strtod( line, cptr);
    if (cptr!=NULL ) d2 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d3 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d4 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d5 = strtod( *cptr, cptr);
    if (cptr!=NULL ) d6 = strtod( *cptr, cptr);
    values[++i] = d1;
    values[++i] = d2;
    values[++i] = d3;
    values[++i] = d4;
    values[++i] = d5;
    values[++i] = d6;
    if ( i>(6*MAXLINES-6) ) break;
  }
  *size = i/6+1;

  return 1;
}

int load_transfer(){
	
	int number = 2;
	int *t_sizes = create_ivector(number);
	ivector_read(&number, transfer_size_file, &t_sizes[0]);
	tlsize = t_sizes[0];
	ksize = t_sizes[1];
	
	transfer_T = (double **)create_array(tlsize,ksize);
	int t_total = tlsize*ksize;
	
	array_read(&t_total, transfer_T_data_file, &transfer_T[0][0]);
	
	if(do_polarisation==1){
		transfer_E = (double **)create_array(tlsize,ksize);
		array_read(&t_total, transfer_E_data_file, &transfer_E[0][0]);
	}
	
	return 0;
}

int load_flat(){
	
	int number = 2;
	int *f_sizes = create_ivector(number);
	ivector_read(&number, flat_size_file, &f_sizes[0]);
	flsize = f_sizes[0];
	fsize = f_sizes[1];
	
	transfer_flat = (double **)create_array(flsize,fsize);
	int f_total = flsize*fsize;
	
	array_read(&f_total, flat_data_file, &transfer_flat[0][0]);
	
	f_max = transfer_flat[0][fsize-1];
	ftau = transfer_flat[0][0];
	
	return 0;
}

int load_bessel(){
	
	int number = 2;
	int *b_sizes = create_ivector(number);
	ivector_read(&number, bessel_size_file, &b_sizes[0]);
	blsize = b_sizes[0];
	xsize = b_sizes[1];
	
	bessel = (double **)create_array(blsize,xsize);
	int b_total = blsize*xsize;
	
	array_read(&b_total, bessel_data_file, &bessel[0][0]);
	tau0 = bessel[0][0];
	bessel[0][0] = 0.0;
	xmax = bessel[0][xsize-1];
	
	kmax = eflag_T_lmax / tau0;
	if(eflag_E_lmax>eflag_T_lmax)kmax = eflag_E_lmax / tau0;

	return 0;
}

int load_limber(){
	
	int number = 2;
	int *lm_sizes = create_ivector(number);
	ivector_read(&number, limber_size_file, &lm_sizes[0]);
	lmlsize = lm_sizes[0];
	tausize = lm_sizes[1];
	
	limber = (double **)create_array(lmlsize,tausize);
	int lm_total = lmlsize*tausize;
	
	array_read(&lm_total, limber_data_file, &limber[0][0]);

	tmax = limber[0][tausize-1];

	return 0;
}

double get_kmax(){
//	kmax calculated in load_bessel() via l_max	
	return kmax;
}

double get_xmax(){
	return xmax;
}

double get_fmax(){
	return f_max;
}
	
double get_ftau(){
	return ftau;
}

double get_tau0(){
	return tau0;
}

double get_tmax(){
	return tmax;
}

int get_tlsize(){
	return tlsize;
}

int get_flsize(){
	return flsize;
}

int get_blsize(){
	return blsize;
}

int get_lmlsize(){
	return lmlsize;
}

int get_ksize(){
	return ksize;
}

int get_xsize(){
	return xsize;
}

int get_fsize(){
	return fsize;
}

int get_tausize(){
	return tausize;
}

int get_kvec(double *vec){

	int i;
	vec[0] = 0;
	for(i=1;i<ksize;i++){
		vec[i] = transfer_T[0][i];
	}
	
	return 0;
}

int get_xvec(double *vec){
	
	int i;
	vec[0] = 0;
	for(i=1;i<xsize;i++){
		vec[i] = bessel[0][i];
	}
	
	return 0;

}

int get_fkvec(double *vec){

	int i;
	vec[0] = 0;
	for(i=1;i<fsize;i++){
		vec[i] = transfer_flat[0][i];
	}
	
	return 0;
}

int get_tauvec(double *vec){

	int i;
	vec[0] = 0;
	for(i=1;i<tausize;i++){
		vec[i] = limber[0][i];
	}
	
	return 0;
}

int get_tvec_T(int l, double *vec){
	
	int i,ltest,n=0;		
	for ( i = 1; i < tlsize; i++ ) {
		ltest = (int)transfer_T[i][0];
		if (ltest == l) {
			n = i;
			break;
		}
	}
			
	if ( n == 0){
		printf("Warning: Some l values not in transfer table - %d\n", l);
		exit(1);
	}
	
	vec[0] = 0;
	for(i=1;i<ksize;i++){
		vec[i] = transfer_T[n][i];
	}
	
	return 0;
}

int get_tvec_E(int l, double *vec){
	
	int i,ltest,n=0;		
	for ( i = 1; i < tlsize; i++ ) {
		ltest = (int)transfer_E[i][0];
		if (ltest == l) {
			n = i;
			break;
		}
	}
			
	if ( n == 0){
		printf("Warning: Some l values not in transfer table - %d\n", l);
		exit(1);
	}
	
	vec[0] = 0;
	for(i=1;i<ksize;i++){
		vec[i] = transfer_E[n][i];
	}
	
	return 0;
}

int get_fvec(int l, double *vec){
	
	int i,ltest,n=0;		
	for ( i = 1; i < flsize; i++ ) {
		ltest = (int)transfer_flat[i][0];
		if (ltest == l) {
			n = i;
			break;
		}
	}
	
	if ( n == 0){
		printf("Warning: Some l values not in transfer table - %d\n", l);
		exit(1);
	}
	
	vec[0] = 0;
	for(i=1;i<fsize;i++){
		vec[i] = transfer_flat[n][i];
	}
	
	return 0;
}

int get_bvec(int l, double *vec){
	
	int i,ltest,n=0;		
	for ( i = 1; i < blsize; i++ ) {
		ltest = (int)bessel[i][0];
		if (ltest == l) {
			n = i;
			break;
		}
	}
			
	if ( n == 0){
		printf("Warning: Some l values not in bessel table - %d\n", l);
		exit(1);
	}
	
	vec[0] = 0;
	for(i=1;i<xsize;i++){
		vec[i] = bessel[n][i];
	}
	
	return 0;

}

int get_lmvec(int l, double *vec){
	
	int i,ltest,n=0;		
	for ( i = 1; i < lmlsize; i++ ) {
		ltest = (int)limber[i][0];
		if (ltest == l) {
			n = i;
			break;
		}
	}
	
	if ( n == 0){
		printf("Warning: Some l values not in limber table - %d\n", l);
		exit(1);
	}
	
	vec[0] = 0;
	for(i=1;i<tausize;i++){
		vec[i] = limber[n][i];
	}
	
	return 0;

}