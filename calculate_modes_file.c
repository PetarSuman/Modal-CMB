#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include "global.h"

// Declare external functions

void calculate_modes_file(char* filename, int models, double **modearray){
	
	FILE* fp;
	char line[MAXLEN];
	char** endptr = (char**) malloc(sizeof(char*));
	double bispectrum[models];
	int error;
	double time1, duration;
	
	int i,j,k,t1,t2,t3,n;
	double x,y,z;
	int lmax = get_lmax();
	double result = 0.0;
	double s1,s2,s3;
	int terms = get_terms_late();
	
	if ( !(fp = fopen(filename, "r")) ) {
		perror(filename);
	}

	for(i=0;i<models;i++){
		for(j=0;j<terms;j++){
			modearray[i][j] = 0.0;
		}
	}
	
	n=0;
	time1 = csecond();
	
	while (fgets(line, MAXLEN, fp)) {
		error=0;
// 		printf("read line\n");
		if (strlen(line) < 6+8*models+1) continue;
// 		printf("long enought\n");
		
		l1 = (int)strtod( line, endptr);
// 		printf("%d\n",l1);
		if (endptr==NULL || l1<2 || l1>lmax ) continue;
		l2 = (int)strtod( *endptr, endptr);
// 		printf("%d\n",l2);
		if (endptr==NULL || l2<2 || l2>lmax ) continue;
		l3 = (int)strtod( *endptr, endptr);
// 		printf("%d\n",l3);
		if (endptr==NULL || l3<2 || l3>lmax ) continue;
		
		for(i=0;i<models;i++){
			bispectrum[i] = strtod( *endptr, endptr);
			if (endptr==NULL || fabs(bispectrum[i]>1.0)){
				error = 1;
				break;
			}
		}
		if(error==1) continue;
		
		n++;
		
		s1 = pow(2.0*l1+1.0,1.0/3.0)*(get_cl(l1)+get_noise(l1)/(get_beam(l1)*get_beam(l1)));
		s2 = pow(2.0*l2+1.0,1.0/3.0)*(get_cl(l2)+get_noise(l2)/(get_beam(l2)*get_beam(l2)));
		s3 = pow(2.0*l3+1.0,1.0/3.0)*(get_cl(l3)+get_noise(l3)/(get_beam(l3)*get_beam(l3)));
		z = permsix(l1,l2,l3)*calculate_geometric(l1,l2,l3)/sqrt(s1 * s2 * s3);
// 		printf("%d\t%d\t%d",l1,l2,l3);
		for(i=0;i<models;i++){
			x = bispectrum[i];
// 			printf("\t%e",x);
			for(j=0;j<terms;j++){
				y = plijk(j,l1,l2,l3);
				modearray[i][j] += x*y*z;
			}
		}
		if(n%10000==0){
			duration = csecond() - time1;
			printf("done %d in %e\n",n,duration);
			time1 = csecond();
		}
// 		printf("\n");
	}

	fclose(fp);
}