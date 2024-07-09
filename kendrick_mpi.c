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

double integrand(int model1, int model2, int flag, double k1, double k2, double k3);
double weight(double k1, double k2, double k3);
double M(double l1, double l2, double l3);

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


	if (argc != 4) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s data model1 model2\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	char directory[100] = "";
	strcat(directory, argv[1]);
	strcat(directory, "/");
	
// 	int i,j,k,l,m,n;

// mpi vars
	int rank, nproc, lnext;

// mpi init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	
	int models =11;
	char filename[100] = "/nfs/local-cosmos2/projects/planck/jf334.private/bispectrum/smpiout";
// 	char filename[100] = "/home/cosmos/ccc-cam/jf334/Matlab/Data/Bispectrum3D_Equi_SMPI_kcut";
	
	FILE* fp;
	char line[MAXLEN];
	char** endptr = (char**) malloc(sizeof(char*));
	double bispectrum[models];
	int error;
	
	int i,j,k,t1,t2,t3,n;
	double x,y,z;
	int lmax = 2000;
	double result = 0.0;
	double s1,s2,s3;
	int terms = get_terms_late();
	
	if ( !(fp = fopen(filename, "r")) ) {
		perror(filename);
	}
	
	n=0;
	time1 = csecond();
	
	while (fgets(line, MAXLEN, fp)) {
		error=0;
		if (strlen(line) < 6+8*models+1) continue;
		
		l1 = (int)strtod( line, endptr);
// 		if (endptr==NULL || l1<2 || l1>lmax || (l1!=2 && l1%10!=0) ) continue;
		if (endptr==NULL || l1<2 || l1>lmax || (l1!=2 && l1%10!=0) ) continue;
		l2 = (int)strtod( *endptr, endptr);
		if (endptr==NULL || l2<2 || l2>lmax || (l2!=2 && l2%10!=0) ) continue;
		l3 = (int)strtod( *endptr, endptr);
		if (endptr==NULL || l3<2 || l3>lmax || (l3!=2 && l3%10!=0) ) continue;
		
		for(i=0;i<models;i++){
			bispectrum[i] = strtod( *endptr, endptr);
			if (endptr==NULL || fabs(bispectrum[i]>1.0)){
				error = 1;
				break;
			}
		}
		if(error==1) continue;
// 		
		n++;
// 		for(i=0;i<models;i++){
// 			x = bispectrum[i];
// 		}
		printf("%s",line);
// 		if(n%1000==0){
// 			duration = csecond() - time1;
// 			printf("done %d in %e\n",n,duration);
// 			time1 = csecond();
// 			break;
// 		}
// 		printf("\n");
	}

	

// 	printf("[%d] Process[%d]  \n", rank, getpid());
/*	
	load_bessel(directory);
	load_transfer(directory);
	
	int flag = 2;
	
	int model1 = atoi(argv[2]);
	int model2 = atoi(argv[3]);
	
 	int l_max = 800;
 	int resolution;
 	if(flag==0){
 		resolution = l_max;
	} else {
		resolution = 400;
	}
	
 	double xmax = get_xmax();
	double kmin = 2.0/xmax;
	double kmax = (double)l_max/xmax;
	double lvec[resolution+1];
	double kvec[resolution+1];
	for(i=0;i<resolution+1;i++){
		lvec[i] = (double)i;
		kvec[i] = kmin + i*(kmax-kmin)/(double)resolution;
	}
	
	time1 = csecond();
	
	double integral1=0;
	double integral2=0;
	double integral3=0;
	double ***points1 = create_3Darray(2,2,2);
	double ***points2 = create_3Darray(2,2,2);
	double ***points3 = create_3Darray(2,2,2);
	double k1,k2,k3;
	double l1,l2,l3;
 	double cube_size;
 	
 	if(flag==0){
		cube_size = lvec[1]*lvec[1]*lvec[1];
	} else {
		cube_size = (kvec[1]-kvec[0])*(kvec[1]-kvec[0])*(kvec[1]-kvec[0]);
	}
	

	for(i=0;i<resolution;i++){
		for(j=0;j<resolution;j++){
			for(k=0;k<resolution;k++){
				if(i<=j+k+1 && j<=k+i+1 && k<=i+j+1){
					
					for(l=0;l<2;l++){
						for(m=0;m<2;m++){
							for(n=0;n<2;n++){
								l1 = lvec[i+l];
								l2 = lvec[j+m];
								l3 = lvec[k+n];
								k1 = kvec[i+l];
								k2 = kvec[j+m];
								k3 = kvec[k+n];
 								
 								if(flag==0){
									points1[l][m][n] = M(l1,l2,l3)*integrand(model1,model2,flag,k1,k2,k3);
								} else {
									points1[l][m][n] = integrand(model1,model1,flag,k1,k2,k3);
									points2[l][m][n] = integrand(model2,model2,flag,k1,k2,k3);
									points3[l][m][n] = integrand(model1,model2,flag,k1,k2,k3);
								}
							}
						}
					}
					
					integral1 += calculate_volume(i,j,k,points1) * cube_size;
					integral2 += calculate_volume(i,j,k,points2) * cube_size;
					integral3 += calculate_volume(i,j,k,points3) * cube_size;
					
				}
				
			}
		}
	}
	
	time2 = csecond();
	printf("Volume:\t%e\t%e\t%e\t%e\t%e\n", integral3/sqrt(integral2*integral1), integral1, integral2, integral3, time2-time1);
	
// 	printf("[%d] Done.\n", rank);
*/
	MPI_Finalize();
	
	// End of code
	return 0;
}

double integrand(int model1, int model2, int flag, double k1, double k2, double k3){

	double result;
	
	if(k1+k2>=k3 && k2+k3>=k1 && k3+k1>=k2){
		
		double k12 = k1*k1;
		double k22 = k2*k2;
		double k32 = k3*k3;
		double k13 = k12*k1;
		double k23 = k22*k2;
		double k33 = k32*k3;
		double p222 = k12*k22*k32;
		double p333 = k13*k23*k33;
		double p1 = k1+k2+k3;
// 		double w = weight(k1,k2,k3);
		double w = 1.0/p1;
		double s1;
		double s2;
		
		switch (model1) {
			case 1: s1 = dbi(k1, k2, k3); break;
			case 2: s1 = equilateral(k1, k2, k3); break;
			case 3: s1 = smooth2(k1, k2, k3); break;
			case 4: s1 = feature(k1, k2, k3); break;
			case 5: s1 = ghost(k1, k2, k3); break;
			case 6: s1 = local(k1, k2, k3); break;
			case 7: s1 = single(k1, k2, k3); break;
			case 8: s1 = warm(k1, k2, k3); break;
			case 9: s1 = smooth1(k1, k2, k3); break;
			case 10: s1 = nonlocal(k1, k2, k3); break;
			case 11: s1 = maldacena(k1, k2, k3); break;
			case 12: s1 = string(k1, k2, k3); break;
			case 13: s1 = smooth3(k1, k2, k3); break;
			case 14: s1 = constant(k1, k2, k3); break;
			case 15: s1 = chen(k1, k2, k3); break;
			case 16: s1 = orthogonal(k1, k2, k3); break;
			case 17: s1 = mode1(k1, k2, k3); break;
		}
		
		switch (model2) {
			case 1: s2 = dbi(k1, k2, k3); break;
			case 2: s2 = equilateral(k1, k2, k3); break;
			case 3: s2 = smooth2(k1, k2, k3); break;
			case 4: s2 = feature(k1, k2, k3); break;
			case 5: s2 = ghost(k1, k2, k3); break;
			case 6: s2 = local(k1, k2, k3); break;
			case 7: s2 = single(k1, k2, k3); break;
			case 8: s2 = warm(k1, k2, k3); break;
			case 9: s2 = smooth1(k1, k2, k3); break;
			case 10: s2 = nonlocal(k1, k2, k3); break;
			case 11: s2 = maldacena(k1, k2, k3); break;
			case 12: s2 = string(k1, k2, k3); break;
			case 13: s2 = smooth3(k1, k2, k3); break;
			case 14: s2 = constant(k1, k2, k3); break;
			case 15: s2 = chen(k1, k2, k3); break;
			case 16: s2 = orthogonal(k1, k2, k3); break;
			case 17: s2 = mode1(k1, k2, k3); break;
		}
		
		if(flag==0){
			result = p333*s1*s2;
		} else if(flag==1){
			result = p222*fabs(s1);
		} else {
			result = p222*s1*p222*s2*weight(k1,k2,k3);
// 			result = p222*s1*p222*s2/p1;
		}
		
// 		if(fabs(result)>1e2) printf("%e\t%e\t%e\t%e\n",k1,k2,k3,result);
		
	} else {
		result=0;
	}
	
	return result;
 
}

double weight(double k1, double k2, double k3){

	double result;
	
 	double one = 1.0/get_xmax();

	double k = k1+k2+k3;
	double k1t = k2+k3-k1;
	double k2t = k3+k1-k2;
	double k3t = k1+k2-k3;
	
	double p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13;
	
	p1  = k1*(k1+one)*k2*(k2+one)*k3*(k3+one);
	p2  = k + 3*one/2.0;
	p3  = k + one/3.0;
	p4  = sqrt(k1t+one/6.0)*sqrt(k2t+one/6.0)*sqrt(k3t+one/6.0);
	
	p5  = (2*k1+one)*(2*k2+one)*(2+k3+one);
	p6  = k + one;
	p7  = k + 3*one;
	p8  = k1t + one/3.0;
	p9  = k2t + one/3.0;
	p10 = k3t + one/3.0;
	p11 = sqrt(k+one/6.0);
	
	p12 = p1*p2*p2*p3*p4;
	p13 = k*p5*p6*p6*p7*p7*p8*p9*p10*p11;
	result = 0;
	if(p13!=0)result = p12/p13;

	return result;
}

double M(double l1, double l2, double l3){

	double lt = l1+l2+l3;
	double l1t = l2+l3-l1;
	double l2t = l3+l1-l2;
	double l3t = l1+l2-l3;
	
	double p1 = l1*l2*l3;
	double p2 = 2*M_PI*M_PI*M_PI*sqrt(lt*l1t*l2t*l3t);
	double p3 = 0;
	if(p2>0) p3 = 1/p2;
	
	double result = p1*p3;
 
	return result;

}

