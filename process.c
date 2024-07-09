#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include "global.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "wigxjpf.h"

struct record
	{
		int l1,l2,l3;
		float bisp;
	};

int main( int argc, char *argv[] ){
	
// **1**

	double time1, time2, time3, time4, time5, duration;


 	if (argc != 2) {
//	if (argc != 4) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	
// mpi vars
	int rank, nproc, onext, mnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	int lmax = 2000;
	int lmin = 2;
	
	double factor;
	int l1,l2,l3;
	double bisp, bisp2, ppp;
	double x1,x2,x3,y1,y2,y3,s1,s2,s3,b1,b2,b3;
	double ***bispectrum;
	double *bi_vec;
	
	int loops;
	int auxloop;
	int start_loop;
	int end_loop;
	
	bispectrum = create_3Darray(lmax,lmax,lmax);
	
	FILE *ptr_myfile;
	long int size = 334831501;
	struct record my_record[size];


    int i, j, k, n, r;
    int i0, j0, k0;
    int i1, j1, k1;
    int i2, j2, k2;
	int test;
	
	i0 = lmin;
	j0 = lmin;
	k0 = lmin;
	
	check_ijk_min_XXX(lmin,&i0,&j0,&k0);
	
    // Step 1) Count the number of iterations in the original loop.
    long int count = 0;
	i = i0;
	j = j0;
	k = k0;

	for (count = 1; count < lmax*lmax*lmax+1; count++){
		test = get_ijk_next_XXX(lmax,&i,&j,&k);
		if(!test) break;
    }
	count++;
	
	long int bi_size = count;
	bi_vec = create_vector(bi_size);

	ptr_myfile=fopen("/fast/space/projects/planck/jf334.private/TSS_Bispectrum/Blll_TSS.unf","rb");
	if (!ptr_myfile){
		printf("Unable to open file!");
		return 1;
	}
	
	//     fseek(ptr_myfile, 0L, SEEK_END);     /* Position to end of file */
	//     lFileLen = ftell(ptr_myfile);        /* Get file length */
	//     rewind(ptr_myfile);                  /* Back to start of file */
	//
	// s1 = 0e0;
	// s2 = 0e0;
	// s3 = 0e0;
	
	fread(&my_record, size*sizeof(struct record), 1, ptr_myfile);
	// for(i=0;i<20;i++){
	// 	l1 = my_record[i].l1;
	// 	l2 = my_record[i].l2;
	// 	l3 = my_record[i].l3;
	// 	bisp = my_record[i].bisp;
	// 	// printf("%d\t%d\t%d\t%e\n",l1,l2,l3,bisp);
	// }
	fclose(ptr_myfile);
	
	// printf("\n\n\n");

	// loops=Mij_total/nproc;
	// auxloop = fmod(Mij_total,nproc);
	// start_loop = rank*loops;
	// end_loop = (rank+1)*loops-1;
	//
	// if (auxloop != 0){
	// 	if (rank < auxloop){
	// 		start_loop = start_loop + rank;
	// 		end_loop = end_loop + rank + 1;
	// 	}else{
	// 		start_loop = start_loop + auxloop;
	// 		end_loop = end_loop + auxloop;
	// 	}
	// }

	
    // wig_table_init(2*lmax,3);
    // wig_temp_init(2*lmax);
	
	int L,L1,L2,L3;
	double twopi = sqrt(2e0/3.1415927);
	double thrd = 1e0/3e0;
	double sxth = 1e0/6e0;
	double p1,p2,p3,p4,p5,p6,p7,p8;
	
	s1 = 0e0;
	s2 = 0e0;
	s3 = 0e0;
		
	// #pragma omp parallel for private(l1,l2,l3,bisp,L,L1,L2,L3,p1,p2,p3,p4,p5,p6,factor)
	for(i=0;i<size;i++){
		l1 = my_record[i].l1;
		l2 = my_record[i].l2;
		l3 = my_record[i].l3;
		bisp = my_record[i].bisp;
		L = l1+l2+l3;
		L1 = L-2*l1;
		L2 = L-2*l2;
		L3 = L-2*l3;
		p1 = sqrt((2e0*l1+1e0)*(2e0*l2+1e0)*(2e0*l2+1e0)/4e0*3.1415927e0);
		p2 = pow(-1,(L/2)%2)*twopi*pow(L+1e0,-5e-1);
		p3 = pow(L+thrd,5e-1)*pow(L+sxth,-2.5e-1);
		p4 = pow(L1+sxth,2.5e-1)*pow(L1+thrd,-5e-1);
		p5 = pow(L2+sxth,2.5e-1)*pow(L2+thrd,-5e-1);
		p6 = pow(L3+sxth,2.5e-1)*pow(L3+thrd,-5e-1);
		// factor = sqrt((2e0*l1+1e0)*(2e0*l2+1e0)*(2e0*l2+1e0)/4e0*3.1415927e0)*wig3jj(2*l1, 2*l2, 2*l3, 0, 0, 0);
		factor = p1*p2*p3*p4*p5*p6;
		
		x1 = (l3*l3+l2*l2-l1*l1)/(2e0*l3*l2);
		x2 = (l3*l3+l1*l1-l2*l2)/(2e0*l3*l1);
		x3 = (l1*l1+l2*l2-l3*l3)/(2e0*l1*l2);
		y1 = (x1*x1 - 5.5e-1) / (l2*(l2+1e0)*l3*(l3+1e0));
		y2 = (x2*x2 - 5.5e-1) / (l1*(l1+1e0)*l3*(l3+1e0));
		y3 = (x3*x3 - 5.5e-1) / (l2*(l2+1e0)*l1*(l1+1e0));
		// y1 = 1e0 / (l2*(l2+1e0)*l3*(l3+1e0));
		// y2 = 1e0 / (l1*(l1+1e0)*l3*(l3+1e0));
		// y3 = 1e0 / (l2*(l2+1e0)*l1*(l1+1e0));
	    if(l1==40){
	        b1 = 1e0;
	    }else if (l1 < 150){
	        b1 = 4e1*sin((l1-4e1)/4e1)/(l1-4e1);
	    }else{
	        b1 = 4e1/(l1+1e2);
	    }
	    if(l2==40){
	        b2 = 1e0;
	    }else if (l2 < 150){
	        b2 = 4e1*sin((l2-4e1)/4e1)/(l2-4e1);
	    }else{
	        b2 = 4e1/(l2+1e2);
	    }
	    if(l3==40){
	        b3 = 1e0;
	    }else if (l3 < 150){
	        b3 = 4e1*sin((l3-4e1)/4e1)/(l3-4e1);
	    }else{
	        b3 = 4e1/(l3+1e2);
	    }
		b1 = 1e0;
		b2 = 1e0;
		b3 = 1e0;
		
		bisp2 = (b1*b2*y3 + b1*y2*b3 + y1*b2*b3);
		
		ppp = (l1*(l1+1e0)*l2*(l2+1e0)*l3*(l3+1e0));
		
		s1 = s1 + bisp*bisp*ppp;
		s2 = s2 + factor*factor*bisp2*bisp2*ppp;
		s3 = s3 + factor*bisp2*bisp*ppp;
		

		// if(i<20)printf("%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\n",l1,l2,l3,p2,p3,p4,p5,p6,p2*p3*p4*p5*p6);
		// bispectrum[l1][l2][l3] = bisp/factor;
		
		if(l1+l2+l3==1002)printf("%d\t%d\t%d\t%e\t%e\n",l1,l2,l3,bisp/factor,bisp2);
	}
		
	printf("%e\t%e\t%e\t%e\n",s1,s2,s3,s3/sqrt(s1*s2));
	
	// while(fread(&my_record1000, size*sizeof(struct record), size, ptr_myfile)==1){
	// 	l1 = my_record.l1;
	// 	l2 = my_record.l2;
	// 	l3 = my_record.l3;
	// 	bisp = my_record.bisp;
	// 	factor = sqrt((2e0*l1+1e0)*(2e0*l2+1e0)*(2e0*l2+1e0)/4e0*3.1415927e0)*wig3jj(2*l1, 2*l2, 2*l3, 0, 0, 0);
	// 	bispectrum[l1][l2][l3] = bisp/factor;
	//
	//
	//
	// 	// x1 = (l3*l3+l2*l2-l1*l1)/(2e0*l3*l2);
	// 	// x2 = (l3*l3+l1*l1-l2*l2)/(2e0*l3*l1);
	// 	// x3 = (l1*l1+l2*l2-l3*l3)/(2e0*l1*l2);
	// 	// y1 = (x1*x1 - 8e-1) / (l2*(l2+1e0)*l3*(l3+1e0));
	// 	// y2 = (x2*x2 - 8e-1) / (l1*(l1+1e0)*l3*(l3+1e0));
	// 	// y3 = (x3*x3 - 8e-1) / (l2*(l2+1e0)*l1*(l1+1e0));
	// 	//
	// 	// bisp2 = y1+y2+y3;
	// 	//
	// 	// ppp = l1*(l1+1e0)*l2*(l2+1e0)*l3*(l3+1e0);
	// 	//
	// 	// s1 = s1 + ppp*bisp*bisp;
	// 	// s2 = s2 + ppp*bisp*bisp2;
	// 	// s3 = s3 + ppp*bisp2*bisp2;
	// 	//
	// 	// if(l1%10==0 && l2%10==0 && l3%10==0)printf("%d\t%d\t%d\t%e\t%e\n",l1,l2,l3,bisp,bisp2);
	// }
	// fclose(ptr_myfile);
	
	i0 = lmin;
	j0 = lmin;
	k0 = lmin;
	
	check_ijk_min_XXX(lmin,&i0,&j0,&k0);
	
    // Step 1) Count the number of iterations in the original loop.
    count = 0;
	i = i0;
	j = j0;
	k = k0;
	
	for (count = 1; count < lmax*lmax*lmax+1; count++)
	{
		bi_vec[count] = bispectrum[i][j][k];
		test = get_ijk_next_XXX(lmax,&i,&j,&k);
		if(!test) break;
    }
	count++;
	
	ptr_myfile=fopen("/fast/space/projects/planck/jf334.private/TSS_Bispectrum/blll_vec.unf","wb");
	fwrite(bi_vec, (size_t) sizeof(double), (size_t) (bi_size), ptr_myfile);
    fclose(ptr_myfile);
	
	// printf("%e\t%e\t%e\n",s1,s2,s3);
	
	//     x = wig3jj(2*10 , 2*12 , 2*14 , 2*2, 2*0 , 2*(-2));
	// printf("Wigner test (-0.0214552): \t%e\n",x);
	
	
	
	// double* bisp_data = malloc( sizeof(double)*MAXLINES*4);
	// int* bisp_len = malloc( sizeof(int));
	// load_txt_dbl("/home/cosmos/users/jf334/PlanckV3/F90/output", 4, bisp_data, bisp_len);
	//
	// int l1,l2,l3;
	// double bisp;
	// int bisp_size = *bisp_len;
	// double x1,x2,x3,y1,y2,y3;
	// j=0;
	//
	// for (i=0; i<bisp_size; i++){
	// 	l1 = bisp_data[j++];
	// 	l2 = bisp_data[j++];
	// 	l3 = bisp_data[j++];
	// 	// y1 = (1e0*l3+l2-l1)/2e0;
	// 	// y2 = (1e0*l1+l3-l2)/2e0;
	// 	// y3 = (1e0*l1+l2-l3)/2e0;
	// 	// x2 = (1e0-l2/2000e0);
	// 	// y2 = (l3-l1)/2000e0;
	// 	// x3 = (1e0-l3/2000e0);
	// 	// y3 = (l1-l2)/2000e0;
	// 	bisp = bisp_data[j++];
	//     factor = sqrt((2e0*l1+1e0)*(2e0*l2+1e0)*(2e0*l2+1e0)/4e0*3.1415927e0)*wig3jj(2*l1, 2*l2, 2*l3, 0, 0, 0);
	//
	// 	x1 = (l3*l3+l2*l2-l1*l1)/(2e0*l3*l2);
	// 	x2 = (l3*l3+l1*l1-l2*l2)/(2e0*l3*l1);
	// 	x3 = (l1*l1+l2*l2-l3*l3)/(2e0*l1*l2);
	// 	y1 = (x1*x1 - 8e-1) / (l2*(l2+1e0)*l3*(l3+1e0));
	// 	y2 = (x2*x2 - 8e-1) / (l1*(l1+1e0)*l3*(l3+1e0));
	// 	y3 = (x3*x3 - 8e-1) / (l2*(l2+1e0)*l1*(l1+1e0));
	//
	// 	printf("%d\t%d\t%d\t%e\t%e\n",l1,l2,l3,bisp/factor,y1+y2+y3);
	//
	//
	// 	// printf("%e\t%e\t%e\n",(double)l2,y1,bisp/factor);
	// 	// printf("%e\t%e\t%e\n",(double)l3,y1,bisp/factor);
	// 	// printf("%e\t%e\t%e\n",(double)l1,y2,bisp/factor);
	// 	// printf("%e\t%e\t%e\n",(double)l3,y2,bisp/factor);
	// 	// printf("%e\t%e\t%e\n",(double)l1,y3,bisp/factor);
	// 	// printf("%e\t%e\t%e\n",(double)l2,y3,bisp/factor);
	// 	// printf("%e\t%e\t%e\n",x1,-y1,bisp/factor);
	// 	// printf("%e\t%e\t%e\n",x2,-y2,bisp/factor);
	// 	// printf("%e\t%e\t%e\n",x3,-y3,bisp/factor);
	// }
	
	MPI_Finalize();
}