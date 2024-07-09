#include <stdio.h>
#include <stdlib.h>
#include "global.h"

/* Write out array */
void array_write(int *n, char filename[100], double grid[(*n)]){
/*
  This function writes the array "grid" of lenght "n" in a file "filename"
*/

   FILE *fp;
   fp = fopen( filename, "w" );
   if (fp == NULL) {
   		printf("**** Error opening file: %s. ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fwrite(grid, (size_t) sizeof(double), (size_t) (*n), fp);
   fclose(fp);
   return;
}

//void array_write_long(long int *n, char filename[100], double grid[(*n)]){
void array_write_long(long int *n, char filename[100], double* grid){
  
	FILE *fp;
	fp = fopen( filename, "w" );
	//printf("Write size: %zd\n", (size_t) (*n));
	
	if (fp == NULL) {
		printf("**** Error opening file: %s. ****\n", filename);
		printf("**** Program Aborted ****\n");
		exit(1);
	}
	
	size_t output = 0;
	
	int m = 100000000;
	int loop = 1 + (*n / m);
	
	int i;
	double* ptr = grid;
	
	for(i=0;i<loop;i++){
		if(i==loop-1) m = *n - m*(loop-1);
		output += fwrite(ptr, (size_t) sizeof(double), (size_t) (m), fp);
		ptr += m;
	}
	//printf("Write output: %zd\n", output);
	fclose(fp);

	return;
}

/* Read in array */
void array_read(int *n, char filename[100], double grid[(*n)]){

   FILE *fp;
   fp = fopen( filename, "r");
   if (fp == NULL) {
   		printf("**** Error opening file: %s ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fread(grid, (size_t) sizeof(double), (size_t) (*n), fp);
   fclose(fp);

   return;
}

/* Read in array */
//void array_read_long(long int *n, char filename[100], double grid[(*n)]){
void array_read_long(long int *n, char filename[100], double* grid){

  FILE *fp;
	fp = fopen( filename, "r");
  //printf("Read size: %zd\n", (size_t) (*n));
	if (fp == NULL) {
 		printf("**** Error opening file: %s ****\n", filename);
 		printf("**** Program Aborted ****\n");
		exit(1);
  }
  size_t output;
  output = fread(grid, (size_t) sizeof(double), (size_t) (*n), fp);
  //printf("Read output: %zd\n", output);
  fclose(fp);

  return;
}

/* Read in nth column of an array on disk into ptr */
int array_readcol(int n, char filename[100], int nrow, int ncol, double *ptr){
  FILE *fp;
  int i = 0;
  long int skip, offset;
  
  skip = (ncol-1)*sizeof(double);

  if ( (fp = fopen(filename, "r")) != NULL )
    {
      offset = n*sizeof(double);
      fseek( fp, offset, SEEK_SET);
      for ( i=0; i<nrow; i++ )
	{
	  fread (ptr, sizeof(double), 1, fp);
	  ptr++;
	  if ( fseek(fp, skip, SEEK_CUR) != 0  )
	    {
	      printf("**** Fseek error ****\n");
	      exit(1);
	    }
	}
      fclose(fp);
    } else {
      printf("**** Error opening file: %s. ****\n", filename);
      printf("**** Program Aborted ****\n");
      exit(1);
    }
  
  return 0;
}

/* Read in nth row of an array on disk into ptr */
int array_readrow(int n, char filename[100], int nrow, int ncol, double *ptr){
  FILE *fp;
  int i = 0;
  int ret;
  long int offset;
  
  offset = n*ncol*sizeof(double);
  
  if ( (fp = fopen(filename, "r")) != NULL )
    {
      fseek( fp, offset, SEEK_SET);
      ret = fread( ptr, sizeof(double), ncol, fp);
      fclose(fp);
    } else {
      printf("**** Error opening file: %s. ****\n", filename);
      printf("**** Program Aborted ****\n");
      exit(1);
    }
  
  return 0;
}

/* Read in points structure from disk - */
int array_readpoints(int l0, int l1, int l2, int l3, char filename[100], int nrow, int ncol, one_points *points){
  double *drow = create_vector(ncol);
  int i;

  array_readrow( l0, filename, nrow, ncol, drow );
  for ( i=0; i<ncol; i++)
    points[i].x = drow[i];

  array_readrow( l1, filename, nrow, ncol, drow );
  for ( i=0; i<ncol; i++)
    points[i].y1 = drow[i];

  array_readrow( l2, filename, nrow, ncol, drow );
  for ( i=0; i<ncol; i++)
    points[i].y2 = drow[i];

  array_readrow( l3, filename, nrow, ncol, drow );
  for ( i=0; i<ncol; i++)
    points[i].y3 = drow[i];

  return 0;
}

/* Double vector */
double *create_vector(int length){
  return (double*) malloc( (size_t) length*sizeof(double));
}

void destroy_vector(double *vector){
  free(vector);
}

void vector_write(int *n, char filename[100], double grid[(*n)]){

   FILE *fp;
   fp = fopen( filename, "w" );
   if (fp == NULL) {
   		printf("**** Error opening file: %s. ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fwrite(grid, (size_t) sizeof(double), (size_t) (*n), fp);
   fclose(fp);
   
   return;
}

void vector_read(int *n, char filename[100], double grid[(*n)]){

   FILE *fp;
   fp = fopen( filename, "r");
   if (fp == NULL) {
   		printf("**** Error opening file: %s ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fread(grid, (size_t) sizeof(double), (size_t) (*n), fp);
   fclose(fp);

   return;
}
/* Double integer vector */
int *create_ivector(int length){
   return (int*) malloc( (size_t) length*sizeof(int));
}

void destroy_ivector(int *vector){
   free(vector);
}

void ivector_write(int *n, char filename[100], int grid[(*n)]){

   FILE *fp;
   fp = fopen( filename, "w" );
   if (fp == NULL) {
   		printf("**** Error opening file: %s. ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fwrite(grid, (size_t) sizeof(int), (size_t) (*n), fp);
   fclose(fp);
   
   return;
}

void ivector_read(int *n, char filename[100], int grid[(*n)]){

   FILE *fp;
   fp = fopen( filename, "r");
   if (fp == NULL) {
   		printf("**** Error opening file: %s ****\n", filename);
   		printf("**** Program Aborted ****\n");
		exit(1);
   }
   fread(grid, (size_t) sizeof(int), (size_t) (*n), fp);
   fclose(fp);

   return;
}

/* Double 2-D array */

int **create_iarray(int dim_x, int dim_y){
   int i; int **array;
   array=(int**) malloc(dim_x * sizeof(int*) );
   array[0] = (int *) malloc(dim_x*dim_y * sizeof(int) );
   for(i=1;i< dim_x;i++)    /* loop over rows */
     array[i]=array[0] + i * dim_y;
   return array;
}

// int **create_iarray_long(long int dim_x, long int dim_y){
//    int i;
//    int **array;
//    array=(int**) malloc(dim_x * sizeof(int*) );
//    array[0] = (int *) malloc(dim_x*dim_y * sizeof(int) );
//    for(i=1;i< dim_x;i++)    /* loop over rows */
//      array[i]=array[0] + i * dim_y;
//    return array;
// }
void destroy_iarray(int **array){
   free(array[0]);
   free(array);
}

double **create_array(int dim_x, int dim_y){
   int i; double **array;
   
   long int size1;
   long int size2;

   size1 = dim_x;
   size2 = size1*dim_y;

   array = (double**) malloc(size1 * sizeof(double*));
   array[0] = (double*) malloc(size2 * sizeof(double));

   for(i=1;i< dim_x;i++)    /* loop over rows */
     array[i]=array[0] + i * dim_y;

   return array;
}

void destroy_array(double **array){
   free(array[0]);
   free(array);
}

/* Double 3-D array */
double ***create_3Darray(int dim_x, int dim_y, int dim_z){
   int i, j;
   double ***array;
   long int size1;
   long int size2;
   long int size3;

   size1 = dim_x;
   size2 = size1*dim_y;
   size3 = size2*dim_z;

   array = (double***) malloc(size1 * sizeof(double**));
   array[0] = (double **) malloc(size2 * sizeof(double*));
   array[0][0] = (double *) malloc(size3 * sizeof(double));

   for(j=1; j<dim_y; j++)
     array[0][j] = array[0][j-1] + dim_z;

   for(i=1; i<dim_x; i++){
     array[i] = array[i-1] + dim_y;
     array[i][0] = array[i-1][dim_y-1] + dim_z;
     for(j=1; j<dim_y; j++)
       array[i][j] = array[i][j-1] + dim_z;
   }
   return array;
}

void destroy_3Darray(double ***array){
   free(array[0][0]);
   free(array[0]);
   free(array);
}

/* Double 3-D array */
double ****create_4Darray(int dim_w, int dim_x, int dim_y, int dim_z){
   int i, j, k;
   double ****array;
   
   long int size1;
   long int size2;
   long int size3;
   long int size4;

   size1 = dim_w;
   size2 = size1*dim_x;
   size3 = size2*dim_y;
   size4 = size3*dim_z;

   array = (double****) malloc(size1 * sizeof(double***));
   array[0] = (double***) malloc(size2 * sizeof(double**));
   array[0][0] = (double **) malloc(size3 * sizeof(double*));
   array[0][0][0] = (double *) malloc(size4 * sizeof(double));

   for(j=1; j<dim_y; j++){
   	array[0][0][j] = array[0][0][j-1] + dim_z;
   }

   for(i=1; i<dim_x; i++){
   	array[0][i] = array[0][i-1] + dim_y;
   	array[0][i][0] = array[0][i-1][dim_y-1] + dim_z;
   	for(j=1; j<dim_y; j++){
   		array[0][i][j] = array[0][i][j-1] + dim_z;
   	}
   }
// 	printf("4d 1\n");
   for(i=1; i<dim_w; i++){
   	array[i] = array[i-1] + dim_x;
   	array[i][0] = array[i-1][dim_x-1] + dim_y;
   	array[i][0][0] = array[i-1][dim_x-1][dim_y-1] + dim_z;
   	for(k=1; k<dim_y; k++){
   		array[i][0][k] = array[i][0][k-1] + dim_z;
   	}
   	for(j=1; j<dim_x; j++){
   		array[i][j] = array[i][j-1] + dim_y;
   		array[i][j][0] = array[i][j-1][dim_y-1] + dim_z;
   		for(k=1; k<dim_y; k++){
   			array[i][j][k] = array[i][j][k-1] + dim_z;
   		}
   	}
   }

   return array;
}

void destroy_4Darray(double ****array){
   free(array[0][0][0]);
   free(array[0][0]);
   free(array[0]);
   free(array);
}


/* Double 3-D array */
double ***create_3Darray_long(long int dim_x, long int dim_y, long int dim_z){
   int i, j;
   double ***array;

   array=(double***) malloc(dim_x * sizeof(double**) );
   array[0] = (double **) malloc(dim_x*dim_y * sizeof(double*) );
   array[0][0] = (double *) malloc(dim_x*dim_y*dim_z*sizeof(double));

   for(j=1; j<dim_y; j++)
     array[0][j] = array[0][j-1] + dim_z;
   for(i=1; i<dim_x; i++){
     array[i] = array[i-1] + dim_y;
     array[i][0] = array[i-1][dim_y-1] + dim_z;
     for(j=1; j<dim_y; j++)
       array[i][j] = array[i][j-1] + dim_z;
   }
   return array;
}

void destroy_3Darray_long(double ***array){
   free(array[0][0]);
   free(array[0]);
   free(array);
}
double ****create_4Darray_long(long int dim_w, long int dim_x, long int dim_y, long int dim_z){
   int i, j, k;
   double ****array;

   array=(double****) malloc(dim_w * sizeof(double***) );
   array[0]=(double***) malloc(dim_w*dim_x * sizeof(double**) );
   array[0][0] = (double **) malloc(dim_w*dim_x*dim_y * sizeof(double*) );
   array[0][0][0] = (double *) malloc(dim_w*dim_x*dim_y*dim_z*sizeof(double));

   for(j=1; j<dim_y; j++){
   	array[0][0][j] = array[0][0][j-1] + dim_z;
   }

   for(i=1; i<dim_x; i++){
   	array[0][i] = array[0][i-1] + dim_y;
   	array[0][i][0] = array[0][i-1][dim_y-1] + dim_z;
   	for(j=1; j<dim_y; j++){
   		array[0][i][j] = array[0][i][j-1] + dim_z;
   	}
   }
// 	printf("4d 1\n");
   for(i=1; i<dim_w; i++){
   	array[i] = array[i-1] + dim_x;
   	array[i][0] = array[i-1][dim_x-1] + dim_y;
   	array[i][0][0] = array[i-1][dim_x-1][dim_y-1] + dim_z;
   	for(k=1; k<dim_y; k++){
   		array[i][0][k] = array[i][0][k-1] + dim_z;
   	}
   	for(j=1; j<dim_x; j++){
   		array[i][j] = array[i][j-1] + dim_y;
   		array[i][j][0] = array[i][j-1][dim_y-1] + dim_z;
   		for(k=1; k<dim_y; k++){
   			array[i][j][k] = array[i][j][k-1] + dim_z;
   		}
   	}
   }

   return array;
}