#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "global.h"

double cell1(double a, double b, double c, double d, double e, double f, double g, double h);
double cell2(double a, double b, double c, double d, double e, double f, double g, double h);
double cell3(double a, double b, double c, double d, double e, double f, double g, double h);
double cell4(double a, double b, double c, double d, double e, double f, double g, double h);
double cell5(double a, double b, double c, double d, double e, double f, double g, double h);

double calculate_weight(int min, int max, int i, int j, int k) {

	int l,m,n,sum;
	int i0,j0,k0;
	int grid[3][3][3];
	double pt[3][3][3];
	double weight = 0e0;
	
	sum = 0;
	for(l=0; l<3; l++) {
		i0 = i+l-1;
		if(i0<min || i0>max){
			for(m=0; m<3; m++) {
				for(n=0; n<3; n++) {
					grid[l][m][n] = -8;
				}
			}
		}else{
			for(m=0; m<3; m++) {
				j0 = j+m-1;
				if(j0<min || j0>max){
					for(n=0; n<3; n++) {
						grid[l][m][n] = -8;
					}
				}else{
					for(n=0; n<3; n++) {
						k0 = k+n-1;
						if(k0<min || k0>max){
							grid[l][m][n] = -8;
						}else if(i0>(j0+k0) || j0>(i0+k0) || k0>(i0+j0)){
							grid[l][m][n] = 0;
						} else {
							grid[l][m][n] = 1;
							sum++;
						}
					}
				}
			}
		}
	}
	
	if(sum==27){
		weight = 1e0;
	}else{
		
		int s1,s2,s3;
		int a,b,c,d,e,f,g,h;

		for(l=0; l<3; l++) {
			for(m=0; m<3; m++) {
				for(n=0; n<3; n++) {
					pt[l][m][n] = 0e0;
				}
			}
		}
		pt[1][1][1]=1e0;
		
		for(s1=0; s1<2; s1++) {
			for(s2=0; s2<2; s2++) {
				for(s3=0; s3<2; s3++) {
					
					sum=0;
					for(l=0; l<2; l++) {
						for(m=0; m<2; m++) {
							for(n=0; n<2; n++) {
								sum += grid[l+s1][m+s2][n+s3];
							}
						}
					}
					a = pt[s1][s2][s3];
					b = pt[s1+1][s2][s3];
					c = pt[s1][s2+1][s3];
					d = pt[s1+1][s2+1][s3];
					e = pt[s1][s2][s3+1];
					f = pt[s1+1][s2][s3+1];
					g = pt[s1][s2+1][s3+1];
					h = pt[s1+1][s2+1][s3+1];
					
					if (sum==4) {
		
						if (grid[1+s1][0+s2][1+s3] == 1) weight += cell2(f,e,h,g,b,a,d,c);
						if (grid[1+s1][1+s2][0+s3] == 1) weight += cell2(d,c,b,a,h,g,f,e);
						if (grid[0+s1][1+s2][1+s3] == 1) weight += cell2(g,h,e,f,c,d,a,b);
		
					} else if (sum==5) {
		
						weight += cell1(a,b,c,d,e,f,g,h);
		
					} else if (sum==6) {
		
						if (grid[0+s1][0+s2][1+s3] == 1) weight += cell3(g,h,e,f,c,d,a,b);
						if (grid[1+s1][0+s2][0+s3] == 1) weight += cell3(f,h,b,d,e,g,a,c);
						if (grid[0+s1][1+s2][0+s3] == 1) weight += cell3(d,h,c,g,b,f,a,e);
		
					} else if (sum==7) {
		
						if (grid[0+s1][1+s2][0+s3] == 0) weight += cell4(f,e,h,g,b,a,d,c);
						if (grid[0+s1][0+s2][1+s3] == 0) weight += cell4(d,c,b,a,h,g,f,e);
						if (grid[1+s1][0+s2][0+s3] == 0) weight += cell4(g,h,e,f,c,d,a,b);
		
					} else if (sum==8) {
	
						weight += 1e0/8e0;
		
					}
					
				}
			}
		}
	}
	
	if(i!=k){
		if(i==j || j==k){
			weight *= 3e0;
		}else{
			weight *= 6e0;
		}
	}
 	return weight;

}

double calculate_volume(int i, int j, int k, double ***points) {

	double a = points[0][0][0];
	double b = points[1][0][0];
	double c = points[0][1][0];
	double d = points[1][1][0];
	double e = points[0][0][1];
	double f = points[1][0][1];
	double g = points[0][1][1];
	double h = points[1][1][1];

	int l,m,n,sum=0;
	int grid[2][2][2];
	double volume;
	
	for(l=0; l<2; l++) {
		for(m=0; m<2; m++) {
			for(n=0; n<2; n++) {
				if ( (i+l)>(j+k+m+n) || (j+m)>(i+k+l+n) || (k+n)>(i+j+l+m) ){
					grid[l][m][n] = 0;
				} else {
					grid[l][m][n] = 1;
					sum++;
				}
			}
		}
	}
	
	if (sum==4) {
		
		if (grid[1][0][1] == 1) volume = cell2(f,e,h,g,b,a,d,c);
		if (grid[1][1][0] == 1) volume = cell2(d,c,b,a,h,g,f,e);
		if (grid[0][1][1] == 1) volume = cell2(g,h,e,f,c,d,a,b);
		
	} else if (sum==5) {
		
		volume = cell1(a,b,c,d,e,f,g,h);
		
	} else if (sum==6) {
		
		if (grid[0][0][1] == 1) volume = cell3(g,h,e,f,c,d,a,b);
		if (grid[1][0][0] == 1) volume = cell3(f,h,b,d,e,g,a,c);
		if (grid[0][1][0] == 1) volume = cell3(d,h,c,g,b,f,a,e);
		
	} else if (sum==7) {
		
		if (grid[0][1][0] == 0) volume = cell4(f,e,h,g,b,a,d,c);
		if (grid[0][0][1] == 0) volume = cell4(d,c,b,a,h,g,f,e);
		if (grid[1][0][0] == 0) volume = cell4(g,h,e,f,c,d,a,b);
		
	} else if (sum==8) {
	
		volume = cell5(a,b,c,d,e,f,g,h);
		
	} else {
		volume = 0;
	}

 	return volume;

}

double calculate_volume_tri(int i, int j, int k, int l, double ****points) {

	double p1,p2,p3,p4,p5,p6;

	int a,b,c,d,sum=0;
	double volume = 0.0;
	
	for(a=0; a<2; a++) {
		for(b=0; b<2; b++) {
			for(c=0; c<2; c++) {
				for(d=0; d<2; d++) {
					if (quad(i+a,j+b,k+c,l+d)){
						sum++;
					}
				}
			}
		}
	}
// 	printf("%d\t%d\t%d\t%d\t%d\n",i,j,k,l,sum);
	if (sum==5) {
		p1=points[0][0][0][1];
		p2=points[0][0][0][0]+points[1][0][0][1]+points[0][0][1][1]+points[0][1][0][1];
		p3=points[1][0][0][0]+points[0][0][1][0]+points[0][1][0][0]+points[1][0][1][1]+points[1][1][0][1]+points[0][1][1][1];
		p4=points[1][0][1][0]+points[1][1][0][0]+points[0][1][1][0]+points[1][1][1][1];
		p5=points[1][1][1][0];
		volume = (p1+7.0*p2+41.0*p3+191.0*p4+641.0*p5) / 40320.0;
	} else if (sum==11) {
		p1=points[0][0][0][1];
		p2=points[0][0][0][0]+points[1][0][0][1]+points[0][0][1][1]+points[0][1][0][1];
		p3=points[1][0][0][0]+points[0][0][1][0]+points[0][1][0][0]+points[1][0][1][1]+points[1][1][0][1]+points[0][1][1][1];
		p4=points[1][0][1][0]+points[1][1][0][0]+points[0][1][1][0]+points[1][1][1][1];
		p5=points[1][1][1][0];
		volume = (55.0*p1+151.0*p2+315.0*p3+479.0*p4+575.0*p5) / 10080.0;
	} else if (sum==12) {
		p1=points[0][0][0][0]+points[0][0][0][1]+points[0][0][1][0]+points[0][1][0][0]+points[1][0][0][0];
		p2=points[0][0][1][1]+points[0][1][0][1]+points[1][0][0][1]+points[0][1][1][0]+points[1][0][1][0]+points[1][1][0][0];
		p3=points[0][1][1][1]+points[1][0][1][1]+points[1][1][0][1]+points[1][1][1][0];
		p4=points[1][1][1][1];
		volume = (439.0*p1+531.0*p2+599.0*p3+623.0*p4) / 10080.0;
	} else if (sum==14) {
		p1=points[0][0][1][0]+points[0][0][0][1];
		p2=points[0][0][0][0]+points[0][0][1][1];
		p3=points[1][0][0][1]+points[1][0][1][0]+points[0][1][0][1]+points[0][1][1][0];
		p4=points[1][0][0][0]+points[0][1][0][0]+points[0][1][1][1]+points[1][0][1][1];
		p5=points[1][1][0][1]+points[1][1][1][0];
		p6=points[1][1][1][1]+points[1][1][0][0];
		volume = (919.0*p1 + 1069.0*p2 + 1161.0*p3 + 1219.0*p4 + 1239.0*p5 + 1253.0*p6)/20160.0;
	} else if (sum==15) {
		p1=points[0][0][0][1];
		p2=points[0][0][0][0]+points[1][0][0][1]+points[0][0][1][1]+points[0][1][0][1];
		p3=points[1][0][0][0]+points[0][0][1][0]+points[0][1][0][0]+points[1][0][1][1]+points[1][1][0][1]+points[0][1][1][1];
		p4=points[1][0][1][0]+points[1][1][0][0]+points[0][1][1][0]+points[1][1][1][1];
		p5=points[1][1][1][0];
		volume = (1879.0*p1+2329.0*p2+2479.0*p3+2513.0*p4+2519.0*p5) / 40320.0;
	} else if (sum==16) {
		p1=points[0][0][0][0];
		p2=points[0][0][0][1]+points[0][0][1][0]+points[0][1][0][0]+points[1][0][0][0];
		p3=points[0][0][1][1]+points[0][1][0][1]+points[1][0][0][1]+points[0][1][1][0]+points[1][0][1][0]+points[1][1][0][0];
		p4=points[0][1][1][1]+points[1][0][1][1]+points[1][1][0][1]+points[1][1][1][0];
		p5=points[1][1][1][1];
		volume = (p1+p2+p3+p4+p5) / 16.0;
	} else {
		volume = 0;
	}

 	return volume;

}

double cell1(double a, double b, double c, double d, double e, double f, double g, double h){

// | a   |  |   f |
// |   d |  | g h |

	double result = (11*(a+b+c+e)+17*(d+g+f)+25*h)/240.0;
	return result;

}

double cell2(double a, double b, double c, double d, double e, double f, double g, double h){

// | a b |  | e   |
// | c   |  |     |

	double result = (47*a+19*(b+c+e)+5*(d+f+g)+h)/720.0;
	return result;

}
double cell3(double a, double b, double c, double d, double e, double f, double g, double h){

// | a b |  |   f |
// | c d |  | g   |

	double result = (40*(b+c)+35*(a+d)+26*(f+g)+19*(e+h))/360.0;
	return result;

}
double cell4(double a, double b, double c, double d, double e, double f, double g, double h){

// | a b |  | e f |
// | c d |  | g   |

	double result = (89*a+85*(b+c+e)+71*(d+f+g)+43*h)/720.0;
	return result;

}
double cell5(double a, double b, double c, double d, double e, double f, double g, double h){

// | a b |  | e f |
// | c d |  | g h |

	double result = (a+b+c+d+e+f+g+h)/8.0;
	return result;
	
}
