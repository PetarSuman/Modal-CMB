#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
// #include <nagmk21.h>
#include "global.h"


/*
1. set up timing variables
2. read in command line arguments for directories and l values
3. variable definitions
4. zero triangle and cell list to stop odd behaviour
5. read in Bessel and transfer data from file
6. create sparse grid of points over the triangle
7. calculate rough estimate of integral for input to area calculations
8. start adaptive algorithm for triangles
9. print result
*/

int main( int argc, char *argv[] ){

	double pi = 3.141592653589793;
	
// **1**

	double time1, time2, time3, time4, time5, duration;


	if (argc != 3) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile shape \n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);

	load_bessel();
	load_transfer();
	load_flat();
	
	model = atoi(argv[2]);
	
// ckpt init
	ckpt_register();

// ini file
	int* lvalues = malloc( sizeof(int)*MAXLINES);
	int* lceil = malloc( sizeof(int));
	int l_size;
	int  lcol = 0;
	int restart = 0;
	int i,j,k,n;
	double total_area = 0;
	double factor = 0;
	int ltotal;
	int l1,l2,l3;

// mpi vars
	int rank, nproc, lnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	printf("[%d] Process[%d]  \n", rank, getpid());


// load l's from ini file or cli 

	printf("[%d] Ini file: %s\n", rank, inifile);

	if(use_l_file==1){
		if ( load_one(inifile, lvalues, lceil ) ) {			
			printf("[%d] Ini load complete\n", rank);
		} else {
			printf("[%d] Ini file error\n", rank);
			MPI_Finalize();
			exit(1);
		}
		l_size = *lceil;
	}else{
		l_size = eflag_lmax+1;
		for(i=0;i<l_size;i++){
			lvalues[i] = i;
		}
	}
	
	
// 	int k_size = get_ksize();
// 	double k_vec[k_size];
// 	double t_vec[k_size];
// 	get_kvec(k_vec);
// 	for ( i = 0; i < k_size; i++ ) {
// 			printf("%e\t",k_vec[i]);
// 		for ( j = 0; j < l_size; j++ ) {
// 			l1 = lvalues[j];
// 			get_tvec(l1, t_vec);
// 			printf("%e\t",t_vec[i]);
// 		}
// 		printf("\n");
// 	}		
	
	if (rank == 0){
		create_bispectrum(l_size);
	}
	
	int n1,n2,n3;
	int nmax = 50;
	
	if(bflag_lset==0) {
		ltotal = (l_size)*(l_size+1)*(l_size+2)/6;
	} else if(bflag_lset==1) {
		ltotal = l_size;
	} else if(bflag_lset==2) {
		ltotal = l_size;
	} else if(bflag_lset==3) {
		n=0;
		n1=0;
		n2=nmax;
		n3=nmax;
		for(i=0;n1<=(2*nmax)/3;i++){
			for(j=0;n3>=n2;j++){
				n2 += 1;
				n3 -= 1;
				n++;
			}
			n1++;
			n2 = nmax-n1;
			if(n2<n1)n2 = n1;
			n3 = 2*nmax-n2-n1;
		}
		ltotal = n;
	} else if(bflag_lset==4) {
		ltotal = 1;
	}
	
	int ltriples[ltotal][6];
	
	if(bflag_lset==0){
		
		n=0;
		for (i=0;i<l_size;i++){
			for (j=i;j<l_size;j++){
				for (k=j;k<l_size;k++){
					ltriples[n][0] = lvalues[i];
					ltriples[n][1] = lvalues[j];
					ltriples[n][2] = lvalues[k];
					ltriples[n][3] = i;
					ltriples[n][4] = j;
					ltriples[n][5] = k;
					n++;
				}
			}
		}
	} else if(bflag_lset==1) {
		
		for (k=0;k<l_size;k++){
			ltriples[k][0] = lvalues[k];
			ltriples[k][1] = lvalues[k];
			ltriples[k][2] = lvalues[k];
			ltriples[k][3] = k;
			ltriples[k][4] = k;
			ltriples[k][5] = k;
		}
		
	} else if(bflag_lset==2) {
		
		for (k=0;k<l_size;k++){
			ltriples[k][0] = lvalues[0];
			ltriples[k][1] = lvalues[k];
			ltriples[k][2] = lvalues[k];
			ltriples[k][3] = 0;
			ltriples[k][4] = k;
			ltriples[k][5] = k;
		}
		
	} else if(bflag_lset==3) {
		
		n=0;
		n1=0;
		n2=nmax;
		n3=nmax;
		for(i=0;n1<=(2*nmax)/3;i++){
			for(j=0;n3>=n2;j++){
				
				ltriples[n][0] = lvalues[n1];
				ltriples[n][1] = lvalues[n2];
				ltriples[n][2] = lvalues[n3];
				ltriples[n][3] = n1;
				ltriples[n][4] = n2;
				ltriples[n][5] = n3;
// 				printf("[%d]\t%d\t%d\t%d\t%d\n",rank,n,lvalues[n1],lvalues[n2],lvalues[n3]);
				
				n2 += 1;
				n3 -= 1;
				n++;
			}
			
			n1++;
			n2 = nmax-n1;
			if(n2<n1)n2 = n1;
			n3 = 2*nmax-n2-n1;
			
		}
		
	} else if(bflag_lset==4) {
		ltriples[0][0] = 2;
		ltriples[0][1] = 250;
		ltriples[0][2] = 250;
		ltriples[0][3] = 0;
		ltriples[0][4] = 5;
		ltriples[0][5] = 5;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	int mode;
	
	double *results =  malloc( sizeof(double)*4);
	results[0] = 0.0;
	results[1] = 0.0;
	results[2] = 0.0;
	results[3] = 0.0;
	
// 	double *test123 =  malloc( sizeof(double)*get_tlsize());
// 	for (i=0;i<l_size;i++)get_tvec(lvalues[i],test123);
// 	printf("%d\n",l_size);

	if(bflag_load==1){
		double* data_raw = malloc( sizeof(double)*MAXLINES*4);
		int* r_size = malloc( sizeof(int));
		
		load_four(bload_file,data_raw,r_size);
		
		int c_size = *r_size;
		double completed[c_size][4];
		
		j = 0;
		for(i=0;i<c_size;i++){
			completed[i][2] = data_raw[j++];
			completed[i][1] = data_raw[j++];
			completed[i][0] = data_raw[j++];
			completed[i][3] = data_raw[j++];
// 			printf("[%d]\t%d\t%d\t%d\t%e\n",rank,(int)completed[i][0],(int)completed[i][1],(int)completed[i][2],completed[i][3]);
		}
		
		int check;
		int l4,l5,l6;
		for(i=0;i<ltotal;i++){
			
			l1 = ltriples[i][0];
			l2 = ltriples[i][1];
			l3 = ltriples[i][2];
			check = 0;
			
			for(j=0;j<c_size;j++){
				l4 = (int)completed[j][0];
				l5 = (int)completed[j][1];
				l6 = (int)completed[j][2];
				if( l4==l1 && l5==l2 && l6==l3 ){
					check = 1;
					ltriples[i][0] = 0;
					ltriples[i][1] = 0;
					ltriples[i][2] = 0;
					if(rank==0){
						results[0] = ltriples[i][3];
						results[1] = ltriples[i][4];
						results[2] = ltriples[i][5];
						results[3] = completed[j][3];
						update_bispectrum(results);
						printf("[%d]\t%d\t%d\t%d\t%e\n",rank,l1,l2,l3,results[3]);
					}
				}
			}
			
			if (check == 0){
				printf("[%d]\t%d\t%d\t%d\t%e\n",rank,l1,l2,l3,0.0);
				ltriples[i][0] = 0;
				ltriples[i][1] = 0;
				ltriples[i][2] = 0;
				if(rank==0){
					results[0] = ltriples[i][3];
					results[1] = ltriples[i][4];
					results[2] = ltriples[i][5];
					results[3] = 0.0;
					update_bispectrum(results);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	

	// master record total number of tasks	
	if ( rank == 0 ) record_tasks(ltotal);
	lnext = 0;
	// mpi loop over the lvalues 
	while ( (lnext = get_next_task(3,4,rank,results)) < ltotal ) {
		
		time1 = csecond();
		
		l1 = ltriples[lnext][0];
		l2 = ltriples[lnext][1];
		l3 = ltriples[lnext][2];
		
		int test = l1+l2+l3 - 2*((l1+l2+l3)/2);
		
		if (l1 == 0 || l2 == 0 || l3 == 0 || test || l3 > l1 + l2){
			// undefined
			mode = 0;
			continue;
		} else if (l1 < l_flat || l2 < l_flat || l3 < l_flat){
			// do full calculation
			mode = 1;
		} else if (l1 >= l_flat && l2 >= l_flat && l3 >= l_flat){
			// use flat sky approximation
			mode = 2;
		}
		
		if (mode == 0) {
			total_area = 0.0;
			
		} else if (mode == 1){
			
			total_area = calculate_full(l1, l2, l3);
			
		} else if (mode == 2) {
			
			total_area = calculate_flat(l1, l2, l3);
			
		}
		
// 		factor = 6.0*pow(2.0*M_PI,6);
		
		results[0] = (double)ltriples[lnext][3];
		results[1] = (double)ltriples[lnext][4];
		results[2] = (double)ltriples[lnext][5];
// 		results[3] = string(l1,l2,l3);
		results[3] = total_area;
// 		results[3] = total_area*factor;
		
		sync_tasks(3,4);
		
		time2 = csecond();
		if(mode!=0)printf("[%d]\t%d\t%d\t%d\t%e\t%e\n", rank, l1, l2, l3, results[3], time2-time1);
		
	} // end of MPI loop

	if (rank == 0){
		signal_end_tasks(3,4);
	}

	double x1,x2,x3,y1,y2,y3,z1,z2,z3;
	double bi;

	if (rank == 0 && bflag_lset == 3){
		for(n=0;n<ltotal;n++){
			i = ltriples[n][3];
			j = ltriples[n][4];
			k = ltriples[n][5];
			z1 = (double)lvalues[i]*(lvalues[i]+1);
			z2 = (double)lvalues[j]*(lvalues[j]+1);
			z3 = (double)lvalues[k]*(lvalues[k]+1);
			x1 = (double)(-i+j+k)/(i+j+k);
			x2 = (double)(i-j+k)/(i+j+k);
			x3 = (double)(i+j-k)/(i+j+k);
			y1 = 2.0*(j-k)/(i+j+k);
			y2 = 2.0*(k-i)/(i+j+k);
			y3 = 2.0*(i-j)/(i+j+k);
			factor = (648.0/3.1415927)*(z1*z2*z3)/(z1+z2+z3);
			bi = get_bispectrum(i,j,k);
			total_area = factor*bi;
			printf("%e\t%e\t%e\n",x1,y1,total_area);
			printf("%e\t%e\t%e\n",x1,-y1,total_area);
			printf("%e\t%e\t%e\n",x2,y2,total_area);
			printf("%e\t%e\t%e\n",x2,-y2,total_area);
			printf("%e\t%e\t%e\n",x3,y3,total_area);
			printf("%e\t%e\t%e\n",x3,-y3,total_area);
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank == 0){
		
		output_bispectrum(l_size, lvalues);
		
		double ***check = create_3Darray(l_size,l_size,l_size);
		
		int bi_size = l_size * l_size * l_size;
		
		array_read(&bi_size, bispectrum_file, &check[0][0][0]);
		for (i=0;i<l_size;i++){
			printf("%d\t%d\t%d\t%e\n", lvalues[i], lvalues[i], lvalues[i], check[i][i][i]);
		}
		
	}

	printf("[%d] Done.\n", rank);
	
	MPI_Finalize();

	// End of code
	return 0;
}
