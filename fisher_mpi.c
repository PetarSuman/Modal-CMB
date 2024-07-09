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

int main( int argc, char *argv[] ){

// set up timing
	double time1, time2, time3, time4, time5, duration;
	
	if (argc != 4) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s dir1 dir2 model\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	int mode = atoi(argv[3]);
	
	time1 = csecond();

// mpi vars
	int rank, nproc, lnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);


// read in data
	
	char directory1[100] = "";
	strcat(directory1, argv[1]);
	strcat(directory1, "/");

	char l_number_file1[100] = "";
	strcat(l_number_file1, directory1);
// 	strcat(l_number_file1, interpolated_l_size_file);
	strcat(l_number_file1, l_size_file);

	char l_file1[100] = "";
	strcat(l_file1, directory1);
// 	strcat(l_file1, interpolated_l_data_file);
	strcat(l_file1, l_data_file);
	
	char reduced_bi_file1[100] = "";
	strcat(reduced_bi_file1, directory1);
// 	strcat(reduced_bi_file1, interpolated_bispectrum_file);
// 	strcat(reduced_bi_file1, bispectrum_file);
	strcat(reduced_bi_file1, "bispectrum_smpi0");
	
	char directory2[100] = "";
	strcat(directory2, argv[2]);
	strcat(directory2, "/");

	char l_number_file2[100] = "";
	strcat(l_number_file2, directory2);
	strcat(l_number_file2, interpolated_l_size_file);
// 	strcat(l_number_file2, l_size_file);

	char l_file2[100] = "";
	strcat(l_file2, directory2);
	strcat(l_file2, interpolated_l_data_file);
// 	strcat(l_file2, l_data_file);
	
	char reduced_bi_file2[100] = "";
	strcat(reduced_bi_file2, directory2);
	strcat(reduced_bi_file2, interpolated_bispectrum_file);
// 	strcat(reduced_bi_file2, bispectrum_file);
// 	strcat(reduced_bi_file2, "bispectrum_int");
// 	strcat(reduced_bi_file2, atoi(mode));
	
// 	sprintf(reduced_bi_file2, "%sbispectrum_rec%d", directory2,mode);
// 	sprintf(reduced_bi_file2, "%sbispectrum_rec", directory2);
	
	int i,j,k,r,s,t,m,n;
	int number = 1;
	printf("%s\n", l_number_file1);
	printf("%s\n", l_file1);
	printf("%s\n", reduced_bi_file1);
	printf("%s\n", l_number_file2);
	printf("%s\n", l_file2);
	printf("%s\n", reduced_bi_file2);
	
// 	double* data1 = malloc( sizeof(double)*MAXLINES*4);
// 	double* data2 = malloc( sizeof(double)*MAXLINES*4);
// 	int* size1 = malloc( sizeof(int));
// 	int* size2 = malloc( sizeof(int));
// 	load_four("/home/cosmos/users/jf334/Live/bmpi_local", data1, size1);
// 	load_four("/home/cosmos/users/jf334/Live/smpi_local", data2, size2);
// 	int l_size = 499;
// 	
// 	long int l_size_long = l_size;
// 	double ***bispectrum1 = create_3Darray_long(l_size_long,l_size_long,l_size_long);
// 	double ***bispectrum2 = create_3Darray_long(l_size_long,l_size_long,l_size_long);
// 	int *l_values = create_ivector(l_size);
// 	
// 	for(i=0;i<l_size;i++) l_values[i] = i+2;
// 	
// 	for(i=0;i<*size1;i++){
// 		
// 		r = (int)data1[4*i]-2;
// 		s = (int)data1[4*i+1]-2;
// 		t = (int)data1[4*i+2]-2;
// 		n = 4*i+3;
// 		bispectrum1[r][s][t] = data1[n];
// 		
// 		r = (int)data2[4*i]-2;
// 		s = (int)data2[4*i+1]-2;
// 		t = (int)data2[4*i+2]-2;
// 		n = 4*i+3;
// 		bispectrum2[r][s][t] = data2[n];
// 	}
// 	
// 	for(i=0;i<l_size;i++){
// 		for(j=i;j<l_size;j++){
// 			for(k=j;k<l_size;k++){
// 				bispectrum1[j][k][i] = bispectrum1[i][j][k];
// 				bispectrum1[k][i][j] = bispectrum1[i][j][k];
// 				bispectrum1[k][j][i] = bispectrum1[i][j][k];
// 				bispectrum1[j][i][k] = bispectrum1[i][j][k];
// 				bispectrum1[i][k][j] = bispectrum1[i][j][k];
// 				
// 				bispectrum2[j][k][i] = bispectrum2[i][j][k];
// 				bispectrum2[k][i][j] = bispectrum2[i][j][k];
// 				bispectrum2[k][j][i] = bispectrum2[i][j][k];
// 				bispectrum2[j][i][k] = bispectrum2[i][j][k];
// 				bispectrum2[i][k][j] = bispectrum2[i][j][k];
// 			}
// 		}
// 	}
	
// 	int l_size = 501;
// 	int *l_values = create_ivector(l_size);

// 	
	int *l_data_size = create_ivector(number);
	ivector_read(&number, l_number_file1, &l_data_size[0]);

	int l_size = l_data_size[0];

	ivector_read(&number, l_number_file2, &l_data_size[0]);

	int l_size2 = l_data_size[0];

// 	if (l_size!=l_size2){
// 		printf("L file sizes do not match\t%d\t%d\n",l_size,l_size2);
// 		return 1;
// 		exit;
// 	}
// 	
	int *l_values = create_ivector(l_size);
	ivector_read(&l_size, l_file1, &l_values[0]);

	int *l_values2 = create_ivector(l_size);
	ivector_read(&l_size, l_file2, &l_values2[0]);

	for(i=0;i<l_size;i++){
		printf("%d\t%d\t%d\n",i,l_values[i],l_values2[i]);
	}

	for (i=0;i<l_size;i++){
		if (l_values[i]!=l_values2[i]){
			printf("L files do not match, row: %d l1: %d l2: %d\n", i, l_values[i], l_values2[i]);
			return 1;
			exit;
		}
	}
	
// 	double* cls_raw = malloc( sizeof(double)*MAXLINES*2);
// 	int* cl_len = malloc( sizeof(int));
// 
// 	char cls_file[100] = "";
// 	strcat(cls_file, directory1);
// 	strcat(cls_file, "cls.dat");
// 	load_two(cls_file, cls_raw, cl_len);
// 	
// 	int cl_size = *cl_len;
// 	
// 	//printf("cls_len %d\n", *cl_len);
// 	
// 	double cls[cl_size+1];
// 	double cls_l[cl_size+1];
// 	
// 	
// 	j=0;
// 	cls[0] = 0;
// 	cls_l[0] = 0;
// 	for (i=1; i<cl_size+1; i++){
// 		cls_l[i] = cls_raw[j++];
// 		cls[i] = cls_raw[j++] * cls_l[i] * (cls_l[i]+1.0);
// 	}
// 	
// 	if (l_values[l_size-1]>cls_l[cl_size-1]){
// 		printf("Cls do not contain enough l's, max data %d max cl: %d\n", l_values[l_size-1], (int)cls_l[cl_size-1]);
// 		return 1;
// 		exit;
// 	}
// 	
// 	double cls_all[l_size];
// 	
// 	int spl_clsize = cl_size + 4;
// 	int wrk_clsize = 6 * cl_size + 16;
// 	double spl_clk[spl_clsize];
// 	double spl_clc[spl_clsize];
// 	double wrk_cl[wrk_clsize];
// 	int ifail;
// 	double pt, clpt;
// 	
// 	e01baf_(&cl_size,cls_l,cls,spl_clk,spl_clc,&spl_clsize,wrk_cl,&wrk_clsize,&ifail);
// 	
// 	for (i=0; i<l_size; i++){
// 		pt = (double)l_values[i];
// 		e02bbf_(&spl_clsize, spl_clk, spl_clc, &pt, &clpt, &ifail);
// 		cls_all[i] = clpt / (l_values[i]*(l_values[i]+1));
// 	}
	
	long int l_size_long = l_size;
	double ***bispectrum1 = create_3Darray_long(l_size_long,l_size_long,l_size_long);
	double ***bispectrum2 = create_3Darray_long(l_size_long,l_size_long,l_size_long);
	long int bi_size = l_size_long*l_size_long*l_size_long;
	array_read_long(&bi_size, reduced_bi_file1, &bispectrum1[0][0][0]);
	array_read_long(&bi_size, reduced_bi_file2, &bispectrum2[0][0][0]);
	
	create_cl(l_size);
	create_beam(l_size);
	create_noise(l_size);
	create_lens(l_size);
	load_cl(directory1, l_size, l_values);
	load_BN(directory1, l_size, l_values);
	load_lens(directory1, l_size, l_values);

	double factor;
	double b1,b2,b3;
// 	printf("here1\n");
// 	for(i=0;i<l_size;i++){
// 		for(j=i;j<l_size;j++){
// 			for(k=j;k<l_size;k++){
// // 				printf("%d\t%d\t%d\n",i,j,k);
// 				factor = calculate_ISW(l_size,l_values,i,j,k);
// // 				printf("%d\t%d\t%d\t%e\n",i,j,k,factor);
// 				bispectrum1[i][j][k] = factor;
// 				bispectrum1[j][k][i] = factor;
// 				bispectrum1[k][i][j] = factor;
// 				bispectrum1[k][j][i] = factor;
// 				bispectrum1[j][i][k] = factor;
// 				bispectrum1[i][k][j] = factor;
// 			}
// 		}
// 	}

// 	for(i=0;i<l_size;i++){
// 		b1=get_beam(i);
// 		for(j=i;j<l_size;j++){
// 			b2=get_beam(j);
// 			for(k=j;k<l_size;k++){
// 				b3=get_beam(k);
// 				factor = bispectrum2[i][j][k] / (b1*b2*b3);
// 				bispectrum2[i][j][k] = factor;
// 				bispectrum2[j][k][i] = factor;
// 				bispectrum2[k][i][j] = factor;
// 				bispectrum2[k][j][i] = factor;
// 				bispectrum2[j][i][k] = factor;
// 				bispectrum2[i][k][j] = factor;
// 			}
// 		}
// 	}
// 	for(i=0;i<l_size;i++){
// 		b1=get_beam(i);
// 		for(j=i;j<l_size;j++){
// 			b2=get_beam(j);
// 			for(k=j;k<l_size;k++){
// 				b3=get_beam(k);
// 				factor = bispectrum1[i][j][k] * (b1*b2*b3);
// 				bispectrum1[i][j][k] = factor;
// 				bispectrum1[j][k][i] = factor;
// 				bispectrum1[k][i][j] = factor;
// 				bispectrum1[k][j][i] = factor;
// 				bispectrum1[j][i][k] = factor;
// 				bispectrum1[i][k][j] = factor;
// 			}
// 		}
// 	}
	
// 	double fisher=0;
	double fisher1=0;
	double fisher2=0;
	double fisher3=0;
	/*
	double fisher4=0;
	double fisher5=0;
	double fisher6=0;
	double fisher7=0;
	double fisher8=0;
	*/
	double term;
	double cl1,cl2,cl3;
	int l1,l2,l3;
	int test;
	double geometric;
	int l_min = 2;
	int l_max = 500;

	for(i=0;i<l_size;i+=2)printf("%d\t%e\t%e\n", l_values[i], bispectrum1[2][i][i], bispectrum2[2][i][i]);
	
	for(i=0;i<l_size;i++){
		l1 = l_values[i];
// 		cl1 = get_cl(i);
		cl1 = get_beam(i)*get_beam(i)*get_cl(i)+get_noise(i);
		time1 = csecond();
		test = 0;
		if (l1<l_min||l1>l_max) test = 1;
		
		if(test)continue;
		
		for(j=0;j<l_size;j++){
			l2 = l_values[j];
// 			cl2 = get_cl(j);
			cl2 = get_beam(j)*get_beam(j)*get_cl(j)+get_noise(j);
			test = 0;
			if (l2<l_min||l2>l_max) test = 1;
			
			if(test)continue;
			
			for(k=0;k<l_size;k++){
				l3 = l_values[k];
// 				cl3 = get_cl(k);
				cl3 = get_beam(k)*get_beam(k)*get_cl(k)+get_noise(k);
				
				test = (l1+l2+l3) - 2*( (l1+l2+l3)/2 );
				if (l1>l2+l3||l2>l1+l3||l3>l1+l2) test = 1;
				if (cl1==0||cl2==0||cl3==0) test = 1;
				
				if (l3<l_min||l3>l_max) test = 1;
				
				if(test)continue;
				
				geometric = calculate_geometric(l1, l2, l3) * pow((get_beam(i)*get_beam(j)*get_beam(k)),2);
// 				geometric = calculate_geometric(l1, l2, l3);
				term = geometric * (bispectrum1[i][j][k] * bispectrum1[i][j][k]) / (cl1 * cl2 * cl3);
				fisher1 += term;
				term = geometric * (bispectrum2[i][j][k] * bispectrum2[i][j][k]) / (cl1 * cl2 * cl3);
				fisher2 += term;
				term = geometric * (bispectrum1[i][j][k] * bispectrum2[i][j][k]) / (cl1 * cl2 * cl3);
				fisher3 += term;
			}
			
		}
		
		time2 = csecond();
		//printf("Done l1 = %d in %e\n", l1, time2-time1);
		
	}
	
	printf("Fisher coeffecient =\t%e\t%e\t%e\t%e\n",fisher3/sqrt(fisher1*fisher2),fisher1,fisher2,fisher3);
// 	printf("Correlation =\t%e\n", fisher3/sqrt(fisher1*fisher2));
// 	printf("Fisher coeffecients\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", fisher1, fisher2, fisher3, fisher4, fisher5, fisher6, fisher7, fisher8);

// 	int gstep, gmax, gnum;
// 	gstep = 4;
// 	gmax = 500;
// 	gnum = gmax/gstep+1;
// 	int u,v,w;
// 	double x,y;

// 	for(i=0;i<gnum;i++){
// 		for(j=0;j<i+1;j++){
// 			l1 = gstep*i;
// 			l2 = gmax+gstep*(j-i);
// 			l3 = gmax-gstep*j;
// 			u = l1+l2+l3;
// 			v = l2+l3-l1;
// 			w = 2*(l2-l3);
// 			x = (double)v/(double)u;
// 			y = (double)w/(double)u;
// 			printf("%e\t%e\t%e\n",x,y,bispectrum1[l1][l2][l3]);
// 		}
// 	}
	
// 	for(i=0;i<l_size;i+=10){
// 		l1 = l_values[i];
// 		if (l1==0) l1=2;
// 		for(j=i;j<l_size;j+=10){
// 			l2 = l_values[j];
// 			if (l2==0) l2=2;
// 			for(k=j;k<l_size;k+=10){
// 				l3 = l_values[k];
// 				if (l3==0) l3=2;
// 
// 				test = 0;
// 				if (l1>l2+l3||l2>l1+l3||l3>l1+l2) test = 1;
// 				if (cl1==0||cl2==0||cl3==0) test = 1;
// 
// 				if (l3<l_min||l3>l_max) test = 1;
// 
// 				if(test)continue;
// 				printf("%d\t%d\t%d\t%e\n",l1,l2,l3,bispectrum1[l1][l2][l3]);
// 
// 			}
// 
// 		}
// 	}
	
	MPI_Finalize();
	
	return 0;	
}

