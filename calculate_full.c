#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <gsl/gsl_spline.h>
#include "global.h"

static int b_size;
static int t_size;

static double* b_x;
static double* t_x;

static gsl_interp_accel* b_acc1;
static gsl_interp_accel* b_acc2;
static gsl_interp_accel* b_acc3;
static gsl_interp_accel* t_acc1;
static gsl_interp_accel* t_acc2;
static gsl_interp_accel* t_acc3;

static gsl_spline* b_sp1;
static gsl_spline* b_sp2;
static gsl_spline* b_sp3;
static gsl_spline* t_sp1;
static gsl_spline* t_sp2;
static gsl_spline* t_sp3;

int calculate_area(int l1, int l2, int l3, double **array, int *cells_old, int *cells_new, cell *cell_list_old, cell *cell_list_new, int *level, double *estimate);
double calculate_point(double **array, double x, double y);
double calculate_integral_b(double k1, double k2, double k3);
double calculate_integral_t(double k1, double k2, double k3);

double calculate_full(int l1, int l2, int l3){
	int i,j,k,m,n,l,count;
	
	double h1;
	double h2;
	double h3;
	double total_area = 0;
	int orientation;
	int step = pow(2,depth);

	int tri_size = initial_grid_size * step + 1;
	int cells_old = initial_grid_size * initial_grid_size;
	int row_end;
	int flag;
	double x;
	double y;
	double area;
	double point;
	double result;
	double estimate;
	
	double **k_triangle = create_array(tri_size, tri_size);
	
	for (i = 0; i < tri_size; i++) {
		for (j = 0; j < tri_size; j++) {
			k_triangle[i][j] = -1;
		}
	}
	
// 	cell cell_list[cells_old];
	cell* cell_list_old = (cell*)malloc(sizeof(cell) * cells_old);

	b_size = get_xsize();
	t_size = get_ksize();
	
	b_x = create_vector(b_size);
	t_x = create_vector(t_size);
	
	double b_y1[b_size];
	double b_y2[b_size];
	double b_y3[b_size];
	
	double t_y1[t_size];
	double t_y2[t_size];
	double t_y3[t_size];
	
	get_xvec(b_x);
	get_bvec(l1,b_y1);
	get_bvec(l2,b_y2);
	get_bvec(l3,b_y3);
	
	get_kvec(t_x);
	get_tvec(l1,t_y1);
	get_tvec(l2,t_y2);
	get_tvec(l3,t_y3);
	
// 	for (i = 0; i < t_size; i++){
//  		printf("%d\t%e\n", i, t_x[i]);
// 	}
	
	for (i = 0; i < cells_old; i++) {
		cell_list_old[i].i = 0;
		cell_list_old[i].j = 0;
		cell_list_old[i].orientation = 0;
		cell_list_old[i].flag = 0;
		cell_list_old[i].area = 0;
	}
	
	b_acc1 = gsl_interp_accel_alloc();
	b_acc2 = gsl_interp_accel_alloc();
	b_acc3 = gsl_interp_accel_alloc();
	t_acc1 = gsl_interp_accel_alloc();
	t_acc2 = gsl_interp_accel_alloc();
	t_acc3 = gsl_interp_accel_alloc();
	
	b_sp1 =  gsl_spline_alloc (gsl_interp_cspline, b_size);
	b_sp2 =  gsl_spline_alloc (gsl_interp_cspline, b_size);
	b_sp3 =  gsl_spline_alloc (gsl_interp_cspline, b_size);
	t_sp1 =  gsl_spline_alloc (gsl_interp_cspline, t_size);
	t_sp2 =  gsl_spline_alloc (gsl_interp_cspline, t_size);
	t_sp3 =  gsl_spline_alloc (gsl_interp_cspline, t_size);
	
	
	gsl_spline_init(b_sp1,b_x,b_y1,b_size);
	gsl_spline_init(b_sp2,b_x,b_y2,b_size);
	gsl_spline_init(b_sp3,b_x,b_y3,b_size);
	gsl_spline_init(t_sp1,t_x,t_y1,t_size);
	gsl_spline_init(t_sp2,t_x,t_y2,t_size);
	gsl_spline_init(t_sp3,t_x,t_y3,t_size);

/*
	free(b_y1);
	free(b_y2);
	free(b_y3);

	free(t_y1);
	free(t_y2);
	free(t_y3);*/
	
	for ( i = 0; i < tri_size; i = i + step) {
		x = (double)i / (tri_size - 1);
		for ( j = 0; j < tri_size - i; j = j + step ) {
			y = (double)(i+2*j - tri_size + 1) / (tri_size - 1);
			if ( k_triangle[i][j] == -1 ) {
				k_triangle[i][j] = calculate_point(k_triangle, x, y);
			}
		}
	}
	count = 0;
	row_end = 2*initial_grid_size;
	estimate = 0;
	int cells_new = cells_old;
	for ( m = 1; m < initial_grid_size + 1; m++) {
		orientation = 1;
		i = m*step;
		j = 0;
		for ( n = 1; n < row_end; n++ ) {
			switch (orientation) {
				case 1:
					h1 = k_triangle[i][j];
					h2 = k_triangle[i-step][j];
					h3 = k_triangle[i-step][j+step];
				break;
			
				case 2:
					h1 = k_triangle[i][j];
					h2 = k_triangle[i+step][j];
					h3 = k_triangle[i+step][j-step];
				break;
					
				default:
					h1 = 0;
					h2 = 0;
					h3 = 0;
				break;
			}
			
			area = (h1 + h2 + h3) / (3 * initial_grid_size * initial_grid_size);
			estimate += area;
			
			cell_list_old[count].i = i;
			cell_list_old[count].j = j;
			cell_list_old[count].orientation = orientation;
			if (area == 0) {
				cell_list_old[count].flag = 0;
			} else {
				cell_list_old[count].flag = 1;
				cells_new += 3;
			}
			cell_list_old[count].area = area;
			
			count++;
			
			if (orientation == 1) {
				i = i - step;
				j = j + step;
				orientation = 2;
			} else {
				i = i + step;
				orientation = 1;
			}
		}
		row_end = row_end - 2;	
	}
	
	cell *cell_list_new = (cell*)malloc(sizeof(cell) * cells_new);
	flag = 1;
	if(depth==0)flag=0;
	
	for(i=0;flag;i++){
		flag = calculate_area(l1,l2,l3,k_triangle,&cells_old,&cells_new,cell_list_old,cell_list_new,&i,&estimate);
		cell_list_old = (cell*)realloc(cell_list_old,sizeof(cell) * cells_old);
		if (cell_list_old == NULL)printf ("Oh no!\n");

		for(m=0;m<cells_old;m++){
			cell_list_old[m].i = cell_list_new[m].i;
			cell_list_old[m].j = cell_list_new[m].j;
			cell_list_old[m].orientation = cell_list_new[m].orientation;
			cell_list_old[m].flag = cell_list_new[m].flag;
			cell_list_old[m].area = cell_list_new[m].area;
		}
		cell_list_new = (cell*)realloc(cell_list_new,sizeof(cell) * cells_new);
		if (cell_list_new == NULL)printf ("Oh no!\n");
	}
	
	gsl_spline_free(b_sp1);
	gsl_spline_free(b_sp2);
	gsl_spline_free(b_sp3);
	
	gsl_spline_free(t_sp1);
	gsl_spline_free(t_sp2);
	gsl_spline_free(t_sp3);
	
	gsl_interp_accel_free(b_acc1);
	gsl_interp_accel_free(b_acc2);
	gsl_interp_accel_free(b_acc3);
	
	gsl_interp_accel_free(t_acc1);
	gsl_interp_accel_free(t_acc2);
	gsl_interp_accel_free(t_acc3);
	
	total_area = M_2_PI * M_2_PI * M_2_PI * estimate;
	
	return total_area;

}


int calculate_area(int l1, int l2, int l3, double **array, int *cells_old, int *cells_new, cell *cell_list_old, cell *cell_list_new, int *level, double *estimate){

// 	printf("here1\n");
	double total_area = 0;
	int myrank;
	
// 	cell *cell_list_new = (cell*)malloc(sizeof(cell) * *cells_new);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
// 	printf("here2\n");
	int i1;
	int j1;
	int orientation1;
	int k;
	int m;
	int r,s;

	int flag;

	double test;
	double area1 = 0;
	double area2 = 0;
	double area3 = 0;
	double area4 = 0;
	double new_area = 0;

	int step1 = pow(2,depth - *level);
	int step2 = step1 / 2;
	int cell_step1 = step1*step1;
	int cell_step2 = cell_step1 / 4;
	int count = 0;
	int flag_count = 0;
	
	double h1;
	double h2;
	double h3;
	double h4;
	double h5;
	double h6;

	double x1;
	double x2;
	double y;
	double y1;
	double y2;

	double z1,z2,z3;

// 	printf("here3\n");
	int* lvec = create_ivector(3);
	int tri_size = initial_grid_size * pow(2,depth) + 1;
	
	double max_error = fabs(*estimate * tri_accuracy / (initial_grid_size * initial_grid_size * pow(4,*level)));
	
// 	printf("here4\t%d\t%d\t%d\t%d\t%d\t%d\t%e\n",l1,l2,l3,*cells_old,*cells_new,*level,*estimate);
	for (m = 0; m < *cells_old; m ++) {
// **1**
		i1 = cell_list_old[m].i;
		j1 = cell_list_old[m].j;
		
		if (cell_list_old[m].flag == 1) {

// 			printf("here4.1\t%d\n",m);
// 			printf("here4.1\t%d\t%d\t%d\t%d\t%e\n",i1,j1,cell_list[m].orientation,cell_list[m].flag,cell_list[m].area);
			
			switch (cell_list_old[m].orientation) {
				case 1:
// 					printf("here4.2\t%d\n",m);
					x1 = (double)(i1-step2) / (tri_size - 1);
					x2 = (double)(i1-2*step2) / (tri_size - 1);
					y = (double)(i1+2*j1 - tri_size + 1) / (tri_size - 1);
					y1 = (double)(i1+2*j1-step2 - tri_size + 1) / (tri_size - 1);
					y2 = (double)(i1+2*j1+step2 - tri_size + 1) / (tri_size - 1);
// **2**				
// 					printf("here4.3\t%d\n",m);
					if (array[i1-step2][j1] == -1) { array[i1-step2][j1] = calculate_point(array, x1, y1);}
					if (array[i1-step2][j1+step2] == -1) { array[i1-step2][j1+step2] = calculate_point(array, x1, y2);}
					if (array[i1-2*step2][j1+step2] == -1) { array[i1-2*step2][j1+step2] = calculate_point(array, x2, y);}
				
// 					printf("here4.4\t%d\n",m);
					h1 = array[i1][j1];
					h2 = array[i1-step2][j1];
					h3 = array[i1-step2][j1+step2];
					h4 = array[i1-step1][j1];
					h5 = array[i1-2*step2][j1+step2];
					h6 = array[i1-step1][j1+step1];
// 					printf("here4.5\t%d\n",m);
				break;
	
				case 2:
// 					printf("here4.2\t%d\n",m);
					x1 = (double)(i1+step2) / (tri_size - 1);
					x2 = (double)(i1+2*step2) / (tri_size - 1);
					y = (double)(i1+2*j1 - tri_size + 1) / (tri_size - 1);
					y1 = (double)(i1+2*j1-step2 - tri_size + 1) / (tri_size - 1);
					y2 = (double)(i1+2*j1+step2 - tri_size + 1) / (tri_size - 1);
					
// 					printf("here4.3\t%d\t%d\t%d\t%e\n",m,i1+step2,j1-step2,array[i1+step2][j1-step2]);
					
					r = i1+step2;
					s = j1-step2;
// 					printf("here4.4\t%d\t%d\t%d\t%e\n",m,r,s,z1);
					
					if (array[i1+step2][j1-step2] == -1) {
						z1 = calculate_point(array, x1, y1);
// 						printf("here4.4\t%d\t%d\t%d\t%e\n",m,r,s,z1);
						array[r][s] = z1;
// 						printf("here4.5\t%d\t%e\n",m,array[r][s]);
					}
					if (array[i1+step2][j1] == -1) {
						z2 = calculate_point(array, x1, y2);
						array[i1+step2][j1] = z2;
					}
					if (array[i1+2*step2][j1-step2] == -1) {
						z3 = calculate_point(array, x2, y);
						array[i1+2*step2][j1-step2] = z3;
					}
					
// 					if (array[i1+step2][j1-step2] == -1) { array[i1+step2][j1-step2] = calculate_point(array, x1, y1);}
// 					printf("here4.3.1\t%d\n",m);
// 					if (array[i1+step2][j1] == -1) { array[i1+step2][j1] = calculate_point(array, x1, y2);}
// 					printf("here4.3.2\t%d\n",m);
// 					if (array[i1+2*step2][j1-step2] == -1) { array[i1+2*step2][j1-step2] = calculate_point(array, x2, y);}
// 					printf("here4.3.3\t%d\n",m);
				
// 					printf("here4.4\t%d\n",m);
					h1 = array[i1][j1];
					h2 = array[i1+step2][j1];
					h3 = array[i1+step2][j1-step2];
					h4 = array[i1+step1][j1];
					h5 = array[i1+2*step2][j1-step2];
					h6 = array[i1+step1][j1-step1];
// 					printf("here4.5\t%d\n",m);
				break;
			
				default:
					h1 = 0;
					h2 = 0;
					h3 = 0;
					h4 = 0;
					h5 = 0;
					h6 = 0;
				break;
			}
// 			printf("here4.6\t%d\n",m);
// **3**				
			area1 = (h1 + h2 + h3) / (3 * initial_grid_size * initial_grid_size * pow(4,*level+1));
			area2 = (h2 + h3 + h5) / (3 * initial_grid_size * initial_grid_size * pow(4,*level+1));
			area3 = (h2 + h4 + h5) / (3 * initial_grid_size * initial_grid_size * pow(4,*level+1));
			area4 = (h3 + h5 + h6) / (3 * initial_grid_size * initial_grid_size * pow(4,*level+1));
					
// 			printf("here4.7\t%d\n",m);
			test = fabs(cell_list_old[m].area - area1 - area2 - area3 - area4);
// **4**				
// 			printf("here4.8\t%d\n",m);
			if (test > max_error) { flag = 1; } else { flag = 0; }
// **5**				
// 			printf("here4.9\t%d\n",m);
			cell_list_new[count].i = i1;
			cell_list_new[count].j = j1;
			cell_list_new[count].orientation = cell_list_old[m].orientation;
			cell_list_new[count].flag = flag;
			cell_list_new[count].area = area1;
			
			cell_list_new[count + 1].flag = flag;
			cell_list_new[count + 1].area = area2;
		
			cell_list_new[count + 2].flag = flag;
			cell_list_new[count + 2].area = area3;
			cell_list_new[count + 2].orientation = cell_list_old[m].orientation;
			
			cell_list_new[count + 3].flag = flag;
			cell_list_new[count + 3].area = area4;
			cell_list_new[count + 3].orientation = cell_list_old[m].orientation;
			
// 			printf("here4.10\t%d\n",m);
			if (cell_list_old[m].orientation == 1) {
				cell_list_new[count + 1].i = i1 - 2*step2;
				cell_list_new[count + 1].j = j1 + step2;
				cell_list_new[count + 1].orientation = 2;
		
				cell_list_new[count + 2].i = i1 - step2;
				cell_list_new[count + 2].j = j1;
			
				cell_list_new[count + 3].i = i1 - step2;
				cell_list_new[count + 3].j = j1 + step2;
// 			printf("here4.11\t%d\n",m);
			} else {
				cell_list_new[count + 1].i = i1 + 2*step2;
				cell_list_new[count + 1].j = j1 - step2;
				cell_list_new[count + 1].orientation = 1;
			
				cell_list_new[count + 2].i = i1 + step2;
				cell_list_new[count + 2].j = j1;
			
				cell_list_new[count + 3].i = i1 + step2;
				cell_list_new[count + 3].j = j1 - step2;
// 			printf("here4.12\t%d\n",m);
			}
			count += 4;
		} else {
// 			printf("here4.2\t%d\n",m);
			cell_list_new[count].i = i1;
			cell_list_new[count].j = j1;
			cell_list_new[count].orientation = cell_list_old[m].orientation;
			cell_list_new[count].flag = cell_list_old[m].flag;
			cell_list_new[count].area = cell_list_old[m].area;
			count++;
		}
	}
// **6**
// 	printf("here5\n");
	for (m = 0; m < count; m++) {
		total_area += cell_list_new[m].area;
		flag_count += cell_list_new[m].flag;
	}

// 	printf("here6\n");
	test = fabs((*estimate - total_area) / total_area);
	*cells_new = count + 3 * flag_count;
	
	printf("[%d]==> (finished computing level %d) %d %d %d \n",myrank,*level,l1,l2,l3);
	
	*estimate = total_area;
	*cells_old = count;

// 	free(cell_list);
// 	cell_list = (cell*)realloc(cell_list,sizeof(cell) * *cells_new);
// 	for(m=0;m<*cells_old;m++){
// 		cell_list[m].i = cell_list_new[m].i;
// 		cell_list[m].j = cell_list_new[m].j;
// 		cell_list[m].orientation = cell_list_new[m].orientation;
// 		cell_list[m].flag = cell_list_new[m].flag;
// 		cell_list[m].area = cell_list_new[m].area;
// 		printf("here6.1\t%d\n",m);
// 	}
// 	free(cell_list_new);
	
// 	printf("here7\n");
	if (*level < depth && test > 0.005) {
		flag = 1;
	} else {
		flag = 0;
	}
	
// 	printf("here8\t%d\t%d\n",*cells_old,*cells_new);
	return flag;
}

double calculate_point(double **array, double x, double y){

// 	MPI sync
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if ( myrank == 0 ) sync_tasks(3,4);
	
	double result;
	double tolerance = 1/ (initial_grid_size * pow(2,depth + 1));
	
	double part1;
	double part2;
	double part3;
	
	double k1 = (1-x);
	double k2 = (1+x+y)/2;
	double k3 = (1+x-y)/2;
// **1**
	if (k1 > 1) {k1 = 1;}
	if (k2 > 1) {k2 = 1;}
	if (k3 > 1) {k3 = 1;}
// **2**

	if (k1 < tolerance || k2 < tolerance ||  k3 < tolerance ) {
		result = 0;
	} else {
		part1 = slice(k1, k2, k3);
		if (part1 == 0) {
			part2 = 0;
			part3 = 0;
		} else {
			part2 = calculate_integral_t(k1, k2, k3);
			if (part2 == 0) {
				part3 = 0;
			} else {
				part3 = calculate_integral_b(k1, k2, k3);
			}
		}
		result = part1 * part2 * part3;
	}

// **4**
	if ( 1 - k1 < tolerance || 1 - k2 < tolerance ||  1 - k3 < tolerance ) {
		result = 2 * result;
	}
	
// 	printf("point\t%e\n",result);
	return  result;
}

double calculate_integral_b(double k1, double k2, double k3){

	int n;
	double z;
	double z1;
	double z2;
	double z3;
	
	double factor;
	double result = 0;
	double part1;
	double part2;
	double part3;
	
	double min = b_x[0];
	double max = b_x[b_size-1];
	
	double y[b_size];
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, b_size);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	
	y[0] = 0;
	
	for (n = 1; n < b_size; n++) {
		
		z = b_x[n];
		
		z1 = z*k1;
		z2 = z*k2;
		z3 = z*k3;

		factor = z * z;
		
		part1 = gsl_spline_eval(b_sp1,z1,b_acc1);
		part2 = gsl_spline_eval(b_sp2,z2,b_acc2);
		part3 = gsl_spline_eval(b_sp3,z3,b_acc3);
		
		y[n] = factor * part1 * part2 * part3;
		
	}
	
	gsl_spline_init(sp,b_x,y,b_size);
	result = gsl_spline_eval_integ(sp,min,max,acc);
	
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
	return result;
}

double calculate_integral_t(double k1, double k2, double k3){

	int n;
	double z;
	double z1;
	double z2;
	double z3;
	
	double factor;
	double result = 0;
	double part1;
	double part2;
	double part3;
	
	double min = t_x[0];
	double max = t_x[t_size-1];
	
	double y[t_size];
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, t_size);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	
	y[0] = 0;
	
	for (n = 1; n < t_size; n++) {
		
		z = t_x[n];
		
		z1 = z*k1;
		z2 = z*k2;
		z3 = z*k3;

		factor = scale(2.0*z/3.0) / z;
		
		part1 = gsl_spline_eval(t_sp1,z1,t_acc1);
		part2 = gsl_spline_eval(t_sp2,z2,t_acc2);
		part3 = gsl_spline_eval(t_sp3,z3,t_acc3);
		
		y[n] = factor * part1 * part2 * part3;
		
	}
	
	gsl_spline_init(sp,t_x,y,t_size);
	result = gsl_spline_eval_integ(sp,min,max,acc);
	
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
	return result;
}
