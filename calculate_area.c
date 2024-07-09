#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_spline.h>
#include "global.h"

double calculate_area(int l1, int l2, int l3, double **array, int cells_old, int cells_new, cell *cell_list, int level, double estimate, int b_size, double *b_x, gsl_spline* b_sp1, gsl_interp_accel* b_acc1, gsl_spline* b_sp2, gsl_interp_accel* b_acc2, gsl_spline* b_sp3, gsl_interp_accel* b_acc3, int t_size, double *t_x, gsl_spline* t_sp1, gsl_interp_accel* t_acc1, gsl_spline* t_sp2, gsl_interp_accel* t_acc2, gsl_spline* t_sp3, gsl_interp_accel* t_acc3){

	double total_area = 0;
	int myrank;
	
	cell cell_list_new[cells_new];
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	printf("[%d]==> (finished computing level %d) %d %d %d \n",myrank,level,l1,l2,l3);
	
	int i1;
	int j1;
	int orientation1;
	int k;
	int m;

	int flag;

	double test;
	double area1 = 0;
	double area2 = 0;
	double area3 = 0;
	double area4 = 0;
	double new_area = 0;

	int step1 = pow(2,depth - level);
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

	int* lvec = create_ivector(3);
	int tri_size = initial_grid_size * pow(2,depth) + 1;
	
	double max_error = fabs(estimate * tri_accuracy / (initial_grid_size * initial_grid_size * pow(4,level)));
	
	printf("here1\n");
	for (m = 0; m < cells_old; m ++) {
// **1**	
		if (cell_list[m].flag == 1) {

			switch (cell_list[m].orientation) {
				case 1:
					x1 = (double)(cell_list[m].i-step2) / (tri_size - 1);
					x2 = (double)(cell_list[m].i-2*step2) / (tri_size - 1);
					y = (double)(cell_list[m].i+2*cell_list[m].j - tri_size + 1) / (tri_size - 1);
					y1 = (double)(cell_list[m].i+2*cell_list[m].j-step2 - tri_size + 1) / (tri_size - 1);
					y2 = (double)(cell_list[m].i+2*cell_list[m].j+step2 - tri_size + 1) / (tri_size - 1);
// **2**						
					if (array[cell_list[m].i-step2][cell_list[m].j] == -1) { array[cell_list[m].i-step2][cell_list[m].j] = calculate_point(array, x1, y1, b_size, b_x, b_sp1, b_acc1, b_sp2, b_acc2, b_sp3, b_acc3, t_size, t_x, t_sp1, t_acc1, t_sp2, t_acc2, t_sp3, t_acc3);}
					if (array[cell_list[m].i-step2][cell_list[m].j+step2] == -1) { array[cell_list[m].i-step2][cell_list[m].j+step2] = calculate_point(array, x1, y2, b_size, b_x, b_sp1, b_acc1, b_sp2, b_acc2, b_sp3, b_acc3, t_size, t_x, t_sp1, t_acc1, t_sp2, t_acc2, t_sp3, t_acc3);}
					if (array[cell_list[m].i-2*step2][cell_list[m].j+step2] == -1) { array[cell_list[m].i-2*step2][cell_list[m].j+step2] = calculate_point(array, x2, y, b_size, b_x, b_sp1, b_acc1, b_sp2, b_acc2, b_sp3, b_acc3, t_size, t_x, t_sp1, t_acc1, t_sp2, t_acc2, t_sp3, t_acc3);}
				
					h1 = array[cell_list[m].i][cell_list[m].j];
					h2 = array[cell_list[m].i-step2][cell_list[m].j];
					h3 = array[cell_list[m].i-step2][cell_list[m].j+step2];
					h4 = array[cell_list[m].i-step1][cell_list[m].j];
					h5 = array[cell_list[m].i-2*step2][cell_list[m].j+step2];
					h6 = array[cell_list[m].i-step1][cell_list[m].j+step1];
				break;
	
				case 2:
					x1 = (double)(cell_list[m].i+step2) / (tri_size - 1);
					x2 = (double)(cell_list[m].i+2*step2) / (tri_size - 1);
					y = (double)(cell_list[m].i+2*cell_list[m].j - tri_size + 1) / (tri_size - 1);
					y1 = (double)(cell_list[m].i+2*cell_list[m].j-step2 - tri_size + 1) / (tri_size - 1);
					y2 = (double)(cell_list[m].i+2*cell_list[m].j+step2 - tri_size + 1) / (tri_size - 1);
					
					if (array[cell_list[m].i+step2][cell_list[m].j-step2] == -1) { array[cell_list[m].i+step2][cell_list[m].j-step2] = calculate_point(array, x1, y1, b_size, b_x, b_sp1, b_acc1, b_sp2, b_acc2, b_sp3, b_acc3, t_size, t_x, t_sp1, t_acc1, t_sp2, t_acc2, t_sp3, t_acc3);}
					if (array[cell_list[m].i+step2][cell_list[m].j] == -1) { array[cell_list[m].i+step2][cell_list[m].j] = calculate_point(array, x1, y2, b_size, b_x, b_sp1, b_acc1, b_sp2, b_acc2, b_sp3, b_acc3, t_size, t_x, t_sp1, t_acc1, t_sp2, t_acc2, t_sp3, t_acc3);}
					if (array[cell_list[m].i+2*step2][cell_list[m].j-step2] == -1) { array[cell_list[m].i+2*step2][cell_list[m].j-step2] = calculate_point(array, x2, y, b_size, b_x, b_sp1, b_acc1, b_sp2, b_acc2, b_sp3, b_acc3, t_size, t_x, t_sp1, t_acc1, t_sp2, t_acc2, t_sp3, t_acc3);}
				
					h1 = array[cell_list[m].i][cell_list[m].j];
					h2 = array[cell_list[m].i+step2][cell_list[m].j];
					h3 = array[cell_list[m].i+step2][cell_list[m].j-step2];
					h4 = array[cell_list[m].i+step1][cell_list[m].j];
					h5 = array[cell_list[m].i+2*step2][cell_list[m].j-step2];
					h6 = array[cell_list[m].i+step1][cell_list[m].j-step1];
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
// **3**				
			area1 = (h1 + h2 + h3) / (3 * initial_grid_size * initial_grid_size * pow(4,level+1));
			area2 = (h2 + h3 + h5) / (3 * initial_grid_size * initial_grid_size * pow(4,level+1));
			area3 = (h2 + h4 + h5) / (3 * initial_grid_size * initial_grid_size * pow(4,level+1));
			area4 = (h3 + h5 + h6) / (3 * initial_grid_size * initial_grid_size * pow(4,level+1));
						
			test = fabs(cell_list[m].area - area1 - area2 - area3 - area4);
// **4**						
			if (test > max_error) { flag = 1; } else { flag = 0; }
// **5**				
			cell_list_new[count].i = cell_list[m].i;
			cell_list_new[count].j = cell_list[m].j;
			cell_list_new[count].orientation = cell_list[m].orientation;
			cell_list_new[count].flag = flag;
			cell_list_new[count].area = area1;
			
			cell_list_new[count + 1].flag = flag;
			cell_list_new[count + 1].area = area2;
		
			cell_list_new[count + 2].flag = flag;
			cell_list_new[count + 2].area = area3;
			cell_list_new[count + 2].orientation = cell_list[m].orientation;
			
			cell_list_new[count + 3].flag = flag;
			cell_list_new[count + 3].area = area4;
			cell_list_new[count + 3].orientation = cell_list[m].orientation;
			
			if (cell_list[m].orientation == 1) {
				cell_list_new[count + 1].i = cell_list[m].i - 2*step2;
				cell_list_new[count + 1].j = cell_list[m].j + step2;
				cell_list_new[count + 1].orientation = 2;
		
				cell_list_new[count + 2].i = cell_list[m].i - step2;
				cell_list_new[count + 2].j = cell_list[m].j;
			
				cell_list_new[count + 3].i = cell_list[m].i - step2;
				cell_list_new[count + 3].j = cell_list[m].j + step2;
			} else {
				cell_list_new[count + 1].i = cell_list[m].i + 2*step2;
				cell_list_new[count + 1].j = cell_list[m].j - step2;
				cell_list_new[count + 1].orientation = 1;
			
				cell_list_new[count + 2].i = cell_list[m].i + step2;
				cell_list_new[count + 2].j = cell_list[m].j;
			
				cell_list_new[count + 3].i = cell_list[m].i + step2;
				cell_list_new[count + 3].j = cell_list[m].j - step2;
			}
			count += 4;
		} else {
			cell_list_new[count].i = cell_list[m].i;
			cell_list_new[count].j = cell_list[m].j;
			cell_list_new[count].orientation = cell_list[m].orientation;
			cell_list_new[count].flag = cell_list[m].flag;
			cell_list_new[count].area = cell_list[m].area;
			count++;
		}
	}
// **6**
	printf("here2\n");
	for (m = 0; m < count; m++) {
		total_area += cell_list_new[m].area;
		flag_count += cell_list_new[m].flag;
	}

	printf("here3\n");
	test = fabs((estimate - total_area) / total_area);
	cells_new = count + 3 * flag_count;
// **7**
	printf("test %d\t%e\n",cells_new,test);
	
	level = level+1;

	if (level < depth && test > 0.005) {
	printf("here4\n");
		total_area = calculate_area(l1, l2, l3, array, count, cells_new, cell_list_new, level, total_area, b_size, b_x, b_sp1, b_acc1, b_sp2, b_acc2, b_sp3, b_acc3, t_size, t_x, t_sp1, t_acc1, t_sp2, t_acc2, t_sp3, t_acc3);
	}
	
	return total_area;
}
