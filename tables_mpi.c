#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <mpi.h>
#include <time.h>
#include "global.h"

int main( int argc, char *argv[] ) {
	
	double time1, time2, time3, time4, time5, duration;
	
// mpi vars
	int rank, nproc, lnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);


	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s inifile\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
	}

// ini file
	
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	printf("Loaded l.ini\n");
	
	int i,j,k,m,n,l;
	
	int lsize;
	int* lvalues = malloc( sizeof(int)*5000);
	int* lceil = malloc( sizeof(int));

	if(use_l_file==1){
		if ( load_one(inifile, lvalues, lceil ) ) {			
			printf("[%d] Ini load complete\n", rank);
		} else {
			printf("[%d] Ini file error\n", rank);
			MPI_Finalize();
			exit(1);
		}
		lsize = *lceil;
	}else{
		lsize = eflag_T_lmax-1;
		for(i=2;i<eflag_T_lmax+1;i++){
			lvalues[i-2] = i;
			// printf("%d\t%d\n",i,lvalues[i-2]);
		}
	}
	
// Load source data

	
	int two = 2;
	double z;
	
	int *source_sizes = create_ivector(two);
	ivector_read(&two, source_size_file, &source_sizes[0]);
	int tau_size_raw = source_sizes[0];
	int k_size_raw = source_sizes[1];
	
	double **s_data = create_array( tau_size_raw, k_size_raw);
	int s_data_size = tau_size_raw * k_size_raw;
	
// 	if(rank<nproc/2){
	printf("[%d] reading source\n", rank);
	array_read(&s_data_size, source_data_file, &s_data[0][0]);
	printf("[%d] read source\n", rank);
// 	}
	
	
	double tau_vec_raw[tau_size_raw];
	double k_vec_raw[k_size_raw];
	tau_vec_raw[0] = 0;
	k_vec_raw[0] = 0;
	
	for (i = 1; i < tau_size_raw; i++){
		tau_vec_raw[i] = s_data[i][0];
		s_data[i][0] = 0;
	}
	
	for (i = 1; i < k_size_raw; i++){
		k_vec_raw[i] = s_data[0][i];
		s_data[0][i] = 0;
	}
	
	double tau0R = s_data[0][0];
	s_data[0][0] = 0;
 
	double tau_min = tau_vec_raw[0];
	double tau0 = tau_vec_raw[tau_size_raw-1];
	double kmin = k_vec_raw[0];
	double kmax = k_vec_raw[k_size_raw-1];
	
	double step = 2e-5;
	double pt = 0.0;
	n=0;
	for (i = 0; pt<kmax;i++){
		pt+=step;
		if(kmin>pt)n=0;
		n++;
	}
	int k_size = n-1;
	double k_vec[k_size];
	pt = 0.0;
	n = 0;
	for (i=0; i<k_size;i++){
		pt+=step;
		if(kmin>pt)n=0;
		k_vec[n] = pt;
		n++;
	}
	
	double **src_tmp = create_array(k_size, tau_size_raw);
	
	double srck[k_size_raw];
	gsl_spline* sp_k_raw =  gsl_spline_alloc (gsl_interp_cspline, k_size_raw);
	gsl_interp_accel* acc_k_raw = gsl_interp_accel_alloc();
	
	for (i = 0; i < tau_size_raw; i++) {
		for (j = 0; j < k_size_raw; j++) {
			srck[j] = 5.0*s_data[i][j]/3.0;
// 			printf("%d\t%e\n",j,k_vec[j]);
		}
		gsl_spline_init(sp_k_raw,k_vec_raw,srck,k_size_raw);
		for (j = 0; j < k_size; j++) {
			pt = k_vec[j];
			src_tmp[j][i] = gsl_spline_eval(sp_k_raw,pt,acc_k_raw);
		}
	}
	
	gsl_spline_free(sp_k_raw);
	gsl_interp_accel_free(acc_k_raw);

	printf("[%d] K:%d\t%e\n", rank, k_size, kmax);
	
	
	step = 2e-5*tau0;
	pt = 0.0;
	n=0;
	for (i = 0; pt<tau0;i++){
		pt+=step;
		if(tau_min>pt)n=0;
		n++;
	}
	int tau_size = n-1;
	double tau_vec[tau_size];
	pt = 0.0;
	n = 0;
	for (i=0; i<tau_size;i++){
		pt+=step;
		if(tau_min>pt)n=0;
		tau_vec[n] = pt;
		n++;
	}
	
	double **src = create_array(k_size, tau_size);
	
	double srct[tau_size_raw];
	gsl_spline* sp_tau_raw =  gsl_spline_alloc (gsl_interp_cspline, tau_size_raw);
	gsl_interp_accel* acc_tau_raw = gsl_interp_accel_alloc();
	
	for (i = 0; i < k_size; i++) {
		for (j = 0; j < tau_size_raw; j++) {
			srct[j] = src_tmp[i][j];
		}
		gsl_spline_init(sp_tau_raw,tau_vec_raw,srct,tau_size_raw);
		for (j = 0; j < tau_size; j++) {
			pt = tau_vec[j];
			src[i][j] = gsl_spline_eval(sp_tau_raw,pt,acc_tau_raw);
		}
	}
	
	gsl_spline_free(sp_tau_raw);
	gsl_interp_accel_free(acc_tau_raw);

	destroy_array(src_tmp);
	
	printf("[%d] Tau:%d\t%e\t%e\n", rank, tau_size, tau0, tau0R);
  
// **** Create Bessel ****
	if (tflag_bessel==1 && rank==0){
	printf("[%d] Calculating bessel\n",rank);

		int b_size = tau_size;
// 		int b_size = k_size;
		double **b_data = create_array(lsize, b_size);
		
		b_data[0][0] = tau0;

		for ( i = 1; i < b_size; i++ ) {
			b_data[0][i] = tau_vec[i];
// 			b_data[0][i] = tau0*k_vec[i];
		}

		for ( j = 1; j < lsize; j++ ) {
			time1 = csecond();
			l = lvalues[j - 1];
			b_data[j][0] = (double)l;
			for ( i = 1; i < b_size; i++ ) {
				z = b_data[0][i];
// 				b_data[j][i] = gsl_sf_bessel_jl(l,z);
				b_data[j][i] = spherical_bessel(l,z);
			}
			printf("[%d] Bessel: Done l: %d in %e\n",rank, l, csecond()-time1);
		}

		
		int b_data_size = lsize * b_size;
		array_write(&b_data_size, bessel_data_file, &b_data[0][0]);
		
		int *bessel_sizes = create_ivector(two);
		bessel_sizes[0] = lsize;
		bessel_sizes[1] = b_size;
		ivector_write(&two, bessel_size_file, &bessel_sizes[0]);
		
		destroy_array(b_data);
		destroy_ivector(bessel_sizes);
		
		printf("[%d] Bessel table successfully created\n",rank);
	}
	
// **** Create Transfer ****

	
	if (tflag_transfer==1){
	printf("[%d] Calculating transfer\n", rank);
	
		double x[tau_size];
		double y[tau_size];
		for (i = 0; i < tau_size; i++) {
			x[i] = tau_vec[i];
		}
		
		double *results =  malloc( sizeof(double)*(1+k_size));
		for(i=0;i<1+k_size;i++) results[i] = 0.0;
		
		if ( rank == 0 ) record_tasks(lsize-1);
		MPI_Barrier(MPI_COMM_WORLD);

// 		if(rank<nproc/2){
	
		if ( rank == 0 ) create_transfer(lsize, k_size);
		
		time1 = csecond();

		int next;
		
		gsl_spline* sp_tau =  gsl_spline_alloc (gsl_interp_cspline, tau_size);
		gsl_interp_accel* acc_tau = gsl_interp_accel_alloc();
		
		while ( (next = get_next_task(17,1+k_size,rank,results)) < lsize-1 ) {
			time1 = csecond();
			l = lvalues[next];
			results[0] = next+1;
			for ( i = 1; i < k_size; i++ ) {

				for (j = 0; j < tau_size; j++) {
					pt = k_vec[i]*(tau0 - tau_vec[j]);
					y[j] = src[i][j]*spherical_bessel(l,pt);
				}
				gsl_spline_init(sp_tau,x,y,tau_size);
				results[i] = gsl_spline_eval_integ(sp_tau,x[0],x[tau_size-1],acc_tau);
			}
			printf("[%d] Transfer: Done l: %d in %e\n", rank, l, csecond()-time1);

		}
		
		gsl_spline_free(sp_tau);
		gsl_interp_accel_free(acc_tau);
		
		if(rank==0){
			signal_end_tasks(17,1+k_size);
			
			double **t_data = create_array(lsize,k_size);
			
			t_data[0][0] = 0.0;
			
			for (i=1; i<k_size; i++) t_data[0][i] = k_vec[i];
			for (l=1; l<lsize; l++) t_data[l][0] = lvalues[l-1];
			
			for (l=1; l<lsize; l++){
				for (i=1; i<k_size; i++){
// 					if(k_vec[i]>0.5*l/tau0){
						t_data[l][i] = get_transfer(l,i);
// 					}else{
// 						t_data[l][i] = 0.0;
// 					}
// 					printf("%d\t%d\t%e\n",l,i,t_data[l][i]);
				}
			}
			
			int t_data_size = lsize * k_size;
			array_write(&t_data_size, transfer_T_data_file, &t_data[0][0]);
				
			int *transfer_sizes = create_ivector(two);
			transfer_sizes[0] = lsize;
			transfer_sizes[1] = k_size;
			ivector_write(&two, transfer_size_file, &transfer_sizes[0]);
		}
// 		}
		MPI_Barrier(MPI_COMM_WORLD);
		printf("[%d] Transfer table successfully created\n",rank);
	} else if(tflag_transfer==2){
		printf("[%d] Reformatting transfer\n",rank);
	
		double* rawdata = malloc(sizeof(double)*3*MAXLINES);
		int *size = malloc(sizeof(int));
		char filename[200];
		char suffix[10];
		double lwgt;
		
		suffix[0] = '\0';
		fourfig(lvalues[0],suffix);
		
		// printf("%s\n",suffix);
		filename[0] = '\0';
		// printf("%s\n",filename);
		strcat(filename, "/home/cosmos/users/jf334/CAMB-ML/Output/transfer_planck_unlensed_2015_");
		// strcat(filename, "/nfs/local-cosmos2/projects/planck/jf334.private/CAMB/transfer_planck_unlensed_");
		// printf("%s\n",filename);
		strcat(filename, suffix);
		// printf("%s\n",filename);
		strcat(filename, ".dat");
		load_txt_dbl(filename, 3, rawdata, size);
		k_size = *size+1;
		printf("%s\t%d\t%d\n",filename,lsize,k_size);
			
		double **t_data = create_array(lsize+1,k_size);
		double **e_data = create_array(lsize+1,k_size);
		
		for (l=1; l<lsize+1; l++) t_data[l][0] = lvalues[l-1];
		if(do_polarisation==1)for (l=1; l<lsize+1; l++) e_data[l][0] = lvalues[l-1];
		
		for (n=0;n<*size;n++){
			t_data[0][n+1] = rawdata[3*n];
			t_data[1][n+1] = rawdata[3*n+1];
			if(do_polarisation==1){
				e_data[0][n+1] = rawdata[3*n];
				e_data[1][n+1] = sqrt((lvalues[0]+2e0)*(lvalues[0]+1e0)*lvalues[0]*(lvalues[0]-1e0))*rawdata[3*n+2];
			}
		}

		for(i=2;i<lsize+1;i++){
			lwgt = sqrt((lvalues[i-1]+2e0)*(lvalues[i-1]+1e0)*lvalues[i-1]*(lvalues[i-1]-1e0));
			filename[0] = '\0';
			strcat(filename, "/home/cosmos/users/jf334/CAMB-ML/Output/transfer_planck_unlensed_2015_");
			// strcat(filename, "/nfs/local-cosmos2/projects/planck/jf334.private/CAMB/transfer_planck_unlensed_");
			fourfig(t_data[i][0],suffix);
			// sprintf(suffix,"%04d",t_data[i][0]);
			strcat(filename, suffix);
			strcat(filename, ".dat");
			load_txt_dbl(filename, 3, rawdata, size);
			printf("%s\t%d\t%d\n",filename,lsize,k_size);
		
			for (n=0;n<*size;n++){
				t_data[i][n+1] = rawdata[3*n+1];
				if(do_polarisation==1)e_data[i][n+1] = lwgt*rawdata[3*n+2];
			}
		}
			
		int t_data_size = (lsize+1) * k_size;
		array_write(&t_data_size, transfer_T_data_file, &t_data[0][0]);
		if(do_polarisation==1)array_write(&t_data_size, transfer_E_data_file, &e_data[0][0]);
			
		int *transfer_sizes = create_ivector(two);
		transfer_sizes[0] = (lsize+1);
		transfer_sizes[1] = k_size;
		ivector_write(&two, transfer_size_file, &transfer_sizes[0]);

// 		for (i=0; i<k_size; i++){
// 			printf("%e\t%e\t%e\t%e\t%e\n",t_data[0][i],t_data[1][i],t_data[19][i],t_data[199][i],t_data[1999][i]);
// // 			for (l=200; l<221; l++){
// // 				printf("%e\t",t_data[l][i]);
// // 			}
// // 			printf("\n");
// 		}
// 		printf("\n");
// 		printf("\n");
// 		for (i=0; i<k_size; i++){
// 			printf("%e\t%e\t%e\t%e\t%e\n",e_data[0][i],e_data[1][i],e_data[19][i],e_data[199][i],e_data[1999][i]);
// // 			for (l=200; l<221; l++){
// // 				printf("%e\t",e_data[l][i]);
// // 			}
// // 			printf("\n");
// 		}
		printf("[%d] Transfer table successfully created\n",rank);
	}
// **** Create Flat ****

	if (tflag_flat==1 && rank==0){
		
		double s2[k_size], f2[tau_size];
// 		double ***s2_data = create_3Darray(lsize, k_size, tau_size);
		double **s2_data = create_array(k_size, tau_size);
		
		double ti;
		double t2;
		double kz;
		double k2;
		double pt;
		double cs;
		double l2;
		
// 		for (i=0; i<tau_size; i++) {
// 			time1 = csecond();
// 			
// 			ti = 1.0 / (tau0-tau_vec[i]);
// 			t2 = ti*ti;
// 			
// 			for (j=0; j<k_size; j++) {
// 				s2[j] = src[j][i];
// 			}
// 			
// 			gsl_spline* sp_k =  gsl_spline_alloc (gsl_interp_cspline, k_size);
// 			gsl_interp_accel* acc_k = gsl_interp_accel_alloc();
// 			
// 			gsl_spline_init(sp_k,k_vec,s2,k_size);
// 			
// 			for (j=0; j<k_size; j++) {
// 				kz = k_vec[j];
// 				k2 = kz*kz;
// 				cs = M_PI_2*cos(kz*tau_vec[i])*t2;
// 				for (k=1; k<lsize; k++) {
// 					l2 = lvalues[k]*lvalues[k]*t2;
// 					pt = sqrt(k2+l2);
// 					if(pt>=kmin && pt<=kmax){
// 						s2_data[k][j][i] = gsl_spline_eval(sp_k,pt,acc_k) * cs;
// 					} else {
// 						s2_data[k][j][i] = 0.0;
// 					}
// 				}
// 			}
// 			gsl_spline_free(sp_k);
// 			gsl_interp_accel_free(acc_k);
// 			printf("Flat part 1: Done tau: %e in %e\n", tau_vec[i], csecond()-time1);
// 		}
// 		
// 		int f_size = k_size;
// 		double **f_data = create_array(lsize, f_size);
// 		
// 		for (k=1; k<lsize; k++) {
// 			time1 = csecond();
// 			
// 			for (j=1; j<k_size; j++) {
// 				
// 				gsl_spline* sp_k =  gsl_spline_alloc (gsl_interp_cspline, tau_size);
// 				gsl_interp_accel* acc_k = gsl_interp_accel_alloc();
// 				
// 				for (i=0; i<tau_size; i++) f2[i] = s2_data[k][j][i];
// 				
// 				gsl_spline_init(sp_k,tau_vec,f2,tau_size);
// 				f_data[k][j] = gsl_spline_eval_integ(sp_k,tau_min,tau0,acc_k);
// 				
// 				gsl_spline_free(sp_k);
// 				gsl_interp_accel_free(acc_k);
// 				
// 			}
// 			printf("Flat part 2: Done l: %e in %e\n", lvalues[k], csecond()-time1);
// 		}

		int f_size = k_size;
		double **f_data = create_array(lsize, f_size);

		for (l=1; l<lsize; l++) {
		
			time1 = csecond();
			for (i=0; i<tau_size; i++) {
			
				ti = 1.0 / (tau0-tau_vec[i]);
				t2 = ti*ti;
				l2 = lvalues[l]*lvalues[l]*t2;
	
				for (k=0; k<k_size; k++) {
					s2[k] = src[k][i];
// 					printf("%d\t%d\t%e\n",k,i,src[k][i]);
				}
					
				gsl_spline* sp_k =  gsl_spline_alloc (gsl_interp_cspline, k_size);
				gsl_interp_accel* acc_k = gsl_interp_accel_alloc();
				
				gsl_spline_init(sp_k,k_vec,s2,k_size);
				
				for (k=0; k<k_size; k++) {
					kz = k_vec[k];
					k2 = kz*kz;
					pt = sqrt(k2+l2);
					if(pt>=kmin && pt<=kmax){
						cs = M_PI_2*cos(kz*tau_vec[i])*t2;
						s2_data[k][i] = gsl_spline_eval(sp_k,pt,acc_k) * cs;
					} else {
						s2_data[k][i] = 0.0;
					}
				}
				gsl_spline_free(sp_k);
				gsl_interp_accel_free(acc_k);
// 				printf("Flat part 1: Done tau: %e in %e\n", tau_vec[i], csecond()-time1);
			}
			
			for (k=1; k<k_size; k++) {

				gsl_spline* sp_k =  gsl_spline_alloc (gsl_interp_cspline, tau_size);
				gsl_interp_accel* acc_k = gsl_interp_accel_alloc();

				for (i=0; i<tau_size; i++) f2[i] = s2_data[k][i];

				gsl_spline_init(sp_k,tau_vec,f2,tau_size);
				f_data[l][k] = gsl_spline_eval_integ(sp_k,tau_min,tau0,acc_k);

				gsl_spline_free(sp_k);
				gsl_interp_accel_free(acc_k);

			}
			printf("Flat Done l: %d in %e\n", lvalues[l], csecond()-time1);
		}
		
		
		
		f_data[0][0] = tau0-tau0R;
// 		f_data[0][0] = tau0R;
		
		for (i=1; i<f_size; i++) f_data[0][i] = k_vec[i];
		for (l=1; l<lsize; l++) f_data[l][0] = lvalues[l];
		
		int f_data_size = lsize * f_size;
		array_write(&f_data_size, flat_data_file, &f_data[0][0]);
		
		int *flat_sizes = create_ivector(two);
		flat_sizes[0] = lsize;
		flat_sizes[1] = f_size;	
		ivector_write(&two, flat_size_file, &flat_sizes[0]);
		
		destroy_array(f_data);
		destroy_ivector(flat_sizes);
		
		printf("Flat table successfully created\n");
	}
	
// **** Create limber ****
/*
	if (tflag_limber==1){
		double **lm_data = create_array(lsize, tau2_size);
		
		lm_data[0][0] = 0;
		
		for(i=1;i<lsize;i++) lm_data[i][0] = (double)lvalues[i-1];
		for(i=1;i<tau2_size;i++) lm_data[0][i] = tau2_vec[i];
		
		for(j=0;j<tau2_size;j++){
			//x[j] = (7.0/6.0)*(tau0 - tau2_vec[j]);
			x[j] = tau0 - tau2_vec[j];
			if (x[j]>tau0) x[j] = tau0;
			if (x[j]<0.0) x[j] = 0;
		}
		
		for(i=1;i<lsize;i++){
			time1 = csecond();
			l = lvalues[i - 1];
			
			y[0] = kmax;
			for(j=1;j<tau2_size;j++){
				y[j] = (double)l / tau2_vec[j];
				//y[j] = (6.0/7.0)*((double)l / tau2_vec[j]);
				if (y[j]>kmax) y[j] = kmax;
			}
			
			e02def_(&tau2_size,&nxest,&nyest,x,y,lamda,mu,(double *)coefficents,f,wrk,iwrk,&ifail);
			
			for(j=1;j<tau2_size;j++){
				if (y[j] == kmax){
					lm_data[i][j] = 0;
				} else {
					lm_data[i][j] = f[j];
				}
			}
			printf("Limber: Done l: %d in %e\n", l, csecond()-time1);
		}
		
		int out_data_size = lsize * tau2_size;
		array_write(&out_data_size, limber_data_file, &lm_data[0][0]);	
		
		int *limber_sizes = create_ivector(two);
		limber_sizes[0] = lsize;
		limber_sizes[1] = tau2_size;	
		ivector_write(&two, limber_size_file, &limber_sizes[0]);
		
		destroy_array(lm_data, lsize);
		destroy_ivector(limber_sizes);
		
		printf("Limber table successfully created\n");
	}
*/
	MPI_Finalize();
	
	return 0;
}







