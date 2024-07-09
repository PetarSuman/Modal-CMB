#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <mkl.h>
#include <omp.h>

#include "global.h"
//#include "offload_util.h"     we deleted _OFFLOADABLE
//#include "local_barrier.h"

 void decompose_tetrapyd_XXX(tetrapyd_limits* block, int rank, int nranks, int nthreads, int lmin, int lmax) 
{
	
	printf("[%d]Start: %d\t%d\t%d\t%d\n",rank,nranks,nthreads,lmin,lmax);
    int i, j, k, n, r;
    int i0, j0, k0;
    int i1, j1, k1;
    int i2, j2, k2;
	int mpi_i0, mpi_j0, mpi_k0;
	
	i0 = lmin;
	j0 = lmin;
	k0 = lmin;
	
	// printf("[%d]lmin: %d\t%d\t%d\n",rank,i0,j0,k0);
	
	check_ijk_min_XXX(lmin,&i0,&j0,&k0);
	
	// printf("[%d]lmin check: %d\t%d\t%d\n",rank,i0,j0,k0);
	
    // Step 1) Count the number of iterations in the original loop.
    long int count;
	i = i0;
	j = j0;
	k = k0;
	int test;
	long int countmax;
	countmax = lmax;
	countmax *= lmax;
	countmax *= lmax;
	countmax += 1;
	for (count = 1; count < countmax; count++)
	{

		// printf("[%d]test bef: %d\t%d\t%d\t%d\n",rank,test,i,j,k);
		test = get_ijk_next_XXX(lmax,&i,&j,&k);
		// printf("[%d]test aft: %d\t%d\t%d\t%d\n",rank,test,i,j,k);
		if(!test) 
            break;
    }
	count++;
	
	// printf("[%d]count: %d\t%d\t%d\n",rank,count,lmin,lmax);
	
    // Step 2) Divvy up the total iterations between the ranks and threads
    // number of iterations for this task
	
	long int work_per_rank[nranks];
	divide_tasks(count,nranks,work_per_rank);
	long int work_per_thread[nthreads];
	long int tcount = work_per_rank[rank];
	divide_tasks(tcount,nthreads,work_per_thread);

	i = i0;
	j = j0;
	k = k0;
	
	// Find starting ijk for this rank
	long int t;
    for (r = 0; r < rank; r++)
    {
	    for (t = 0; t < work_per_rank[r]; t++)
	    {
    		test = get_ijk_next_XXX(lmax,&i,&j,&k);
		}
    }

    // Find starting ijk for each thread on this rank
    block[0].i_bgn = i;
    block[0].j_bgn = j;
    block[0].k_bgn = k;
    block[0].loops = work_per_thread[0];
	for(n = 1; n < nthreads; n++)
	{
	    for (t = 0; t < work_per_thread[n-1]; t++)
	    {
			test = get_ijk_next_XXX(lmax,&i,&j,&k);
		}
		block[n].i_bgn = i;
	    block[n].j_bgn = j;
	    block[n].k_bgn = k;
	    block[n].loops = work_per_thread[n];
	}

	#ifdef DEBUG
    for (n = 0; n < nthreads; n++)
    {
        printf("[%d] (%d, %d, %d) %d\n", n, block[n].i_bgn, block[n].j_bgn, block[n].k_bgn, block[n].loops);
        fflush(0);
    }
	#endif
}

 int get_ijk_next_XXX(int lmax, int *i1, int *j1, int *k1)
{
    int i,j,k;
	int test = 1;
	i=*i1;
	j=*j1;
	k=*k1;
	
	if(i%2 + j == lmax){
		i++;
		j=i;
		k=j+i%2;
	}else{
		if(k+2 > lmax || k>=i+j){
			j++;
			k=j+i%2;
		}else{
			k=k+2;
		}
	}
	
	if(i+j+k+1>=3*lmax){
		test = 0;
		if(lmax%2==0){
			i=lmax;
			j=lmax;
			k=lmax;
		}else{
			i=lmax-1;
			j=lmax;
			k=lmax;
		}
	}
	
    *i1 = i;
    *j1 = j;
    *k1 = k;
    return test;
}

 void check_ijk_min_XXX(int lmin, int *i, int *j, int *k)
{
	int n,i0,j0,k0;
	
	// Check points are within allowed values
	i0 = *i < lmin ? lmin : *i;
	j0 = *j < lmin ? lmin : *j;
	k0 = *k < lmin ? lmin : *k;
	
	// Check ordering of points to make sure that i<=j<=k
	if(i0>j0){
		n=j0;
		j0=i0;
		i0=n;
	}
	if(j0>k0){
		n=j0;
		j0=k0;
		k0=n;
		if(i0>j0){
			n=j0;
			j0=i0;
			i0=n;
		}
	}
	
	// Check triangle and parity conditions are met by points and move forward one if not 
	if(i0+j0<k0){
		k0 = i0+j0;
	}else if((i0+j0+k0)%2==1){
		k0=k0+1;
	}
	
	*i = i0;
	*j = j0;
	*k = k0; 
}

 void check_ijk_max_XXX(int lmax, int *i, int *j, int *k)
{
	int n,i0,j0,k0;
	
	// Check end points are within allowed values
	i0 = *i > lmax ? lmax : *i;
	j0 = *j > lmax ? lmax : *j;
	k0 = *k > lmax ? lmax : *k;
	
	// Check ordering of points to make sure that i<=j<=k
	if(i0>j0){
		n=j0;
		j0=i0;
		i0=n;
	}
	if(j0>k0){
		n=j0;
		j0=k0;
		k0=n;
		if(i0>j0){
			n=j0;
			j0=i0;
			i0=n;
		}
	}
	
	// Check triangle and parity conditions are met by end points and move back one if not
	if(i0+j0<k0){
		k0 = i0+j0;
	}else if((i0+j0+k0)%2==1){
		if(k0>j0){
			k0=k0-1;
		}else if(j0>i0){
			j0=j0-1;
		}else{
			i0=i0-1;
		}
	}
	
	*i = i0;
	*j = j0;
	*k = k0; 
}

 void decompose_tetrapyd_XXY(tetrapyd_limits* block, int rank, int nranks, int nthreads, int lminX, int lminY, int lmaxX, int lmaxY) 
{
    int i, j, k, n, r;
    int i0, j0, k0;
    int i1, j1, k1;
    int i2, j2, k2;
	int mpi_i0, mpi_j0, mpi_k0;
	
	i0 = lminX;
	j0 = lminX;
	k0 = lminY;
	
	check_ijk_min_XXY(lminX,lminY,&i0,&j0,&k0);
	
    // Step 1) Count the number of iterations in the original loop.
    long int count = 0;
	i = i0;
	j = j0;
	k = k0;
	int test;
	for (count = 1; count < lmaxX*lmaxX*lmaxY+1; count++)
	{
		test = get_ijk_next_XXY(lminX,lminY,lmaxX,lmaxY,&i,&j,&k);
		if(!test) 
            break;
    }
	count++;
	
    // Step 2) Divvy up the total iterations between the ranks and threads
    // number of iterations for this task
	
	long int work_per_rank[nranks];
	divide_tasks(count,nranks,work_per_rank);
	long int work_per_thread[nthreads];
	long int tcount = work_per_rank[rank];
	divide_tasks(tcount,nthreads,work_per_thread);

	i = i0;
	j = j0;
	k = k0;
	
	// Find starting ijk for this rank
	long int t;
    for (r = 0; r < rank; r++)
    {
	    for (t = 0; t < work_per_rank[r]; t++)
	    {
    		test = get_ijk_next_XXY(lminX,lminY,lmaxX,lmaxY,&i,&j,&k);
		}
    }

    // Find starting ijk for each thread on this rank
    block[0].i_bgn = i;
    block[0].j_bgn = j;
    block[0].k_bgn = k;
    block[0].loops = work_per_thread[0];
	for(n = 1; n < nthreads; n++)
	{
	    for (t = 0; t < work_per_thread[n-1]; t++)
	    {
			test = get_ijk_next_XXY(lminX,lminY,lmaxX,lmaxY,&i,&j,&k);
		}
		block[n].i_bgn = i;
	    block[n].j_bgn = j;
	    block[n].k_bgn = k;
	    block[n].loops = work_per_thread[n];
	}

	#ifdef DEBUG
    for (n = 0; n < nthreads; n++)
    {
        printf("[%d] (%d, %d, %d) %d\n", n, block[n].i_bgn, block[n].j_bgn, block[n].k_bgn, block[n].loops);
        fflush(0);
    }
	#endif
}

 int get_ijk_next_XXY(int lminX, int lminY, int lmaxX, int lmaxY, int *i1, int *j1, int *k1)
{
    int i,j,k;
	int test = 1;
	i=*i1;
	j=*j1;
	k=*k1;
	
	if(j == lmaxX && k+2 > lmaxY){
		i++;
		j=i;
		k=(i+j)%2+lminY;
	}else{
		if(k+2 > min(i+j,lmaxY)){
			j++;
			k=max((j-i),lminY);
			k = k + (i+j+k)%2;
		}else{
			k=k+2;
		}
	}
	
	if(i+j+k+1>=2*lmaxX+lmaxY){
		test = 0;
		if(lmaxY%2==0){
			i=lmaxX;
			j=lmaxX;
			k=lmaxY;
		}else{
			i=lmaxX;
			j=lmaxX;
			k=lmaxY-1;
		}
	}
	
    *i1 = i;
    *j1 = j;
    *k1 = k;
    return test;
}

 void check_ijk_min_XXY(int lminX, int lminY, int *i, int *j, int *k)
{
	int n,i0,j0,k0;
	
	// Check points are within allowed values
	i0 = *i < lminX ? lminX : *i;
	j0 = *j < lminX ? lminX : *j;
	k0 = *k < lminY ? lminY : *k;
	
	// Check ordering of points to make sure that i<=j
	if(i0>j0){
		n=j0;
		j0=i0;
		i0=n;
	}
	
	// Check triangle and parity conditions are met by points and move forward one if not 
	if(i0+j0<k0){
		if(i0+j0<lminY){
			if(lminY>2*lminX){
				i0=lminY/2;
				j0=i0 + lminY%2;
				k0=lminY;
			}else{
				i0=lminX;
				j0=lminX;
				k0=lminY+lminY%2;
			}
		}else{
			k0 = i0+j0;
		}
	}else if(i0+k0<j0){
		j0 = i0+k0;
	}else if((i0+j0+k0)%2==1){
		k0=k0+1;
	}
	
	*i = i0;
	*j = j0;
	*k = k0; 
}

 void check_ijk_max_XXY(int lmaxX, int lmaxY, int *i, int *j, int *k)
{
	int n,i0,j0,k0;
	
	// Check end points are within allowed values
	i0 = *i > lmaxX ? lmaxX : *i;
	j0 = *j > lmaxX ? lmaxX : *j;
	k0 = *k > lmaxY ? lmaxY : *k;
	
	// Check ordering of points to make sure that i<=j
	if(i0>j0){
		n=j0;
		j0=i0;
		i0=n;
	}
	
	// Check triangle and parity conditions are met by end points and move back one if not
	if(i0+j0<k0){
		k0 = i0+j0;
	}else if(i0+k0<j0){
		j0=i0+k0;
	}else if((i0+j0+k0)%2==1){
		k0=k0-1;
	}
	
	*i = i0;
	*j = j0;
	*k = k0; 
}

void decompose_tetrapyd_prim(tetrapyd_limits* block, int rank, int nranks, int nthreads, int min, int max) 
{
    int i, j, k, n, r;
    int i0, j0, k0;
    int i1, j1, k1;
    int i2, j2, k2;
	int mpi_i0, mpi_j0, mpi_k0;
	
	i0 = min;
	j0 = min;
	k0 = min;
	
	check_ijk_min_prim(min,&i0,&j0,&k0);
	
    // Step 1) Count the number of iterations in the original loop.
    long int count = 0;
	i = i0;
	j = j0;
	k = k0;
	int test;
	for (count = 1; count < (max+1)*(max+1)*(max+1); count++)
	{
		test = get_ijk_next_prim(max,&i,&j,&k);
		if(!test)
		    break;
    }
	count++;

    // Step 2) Divvy up the total iterations between the ranks and threads
    // number of iterations for this task
	
	long int work_per_rank[nranks];
	divide_tasks(count,nranks,work_per_rank);
	long int work_per_thread[nthreads];
	long int tcount = work_per_rank[rank];
	divide_tasks(tcount,nthreads,work_per_thread);
	
	i = i0;
	j = j0;
	k = k0;
	
	// Find starting ijk for this rank
	long int t;
    for (r = 0; r < rank; r++)
    {
	    for (t = 0; t < work_per_rank[r]; t++)
	    {
    		test = get_ijk_next_prim(max,&i,&j,&k);
		}
    }

    // Find starting ijk for each thread on this rank
    block[0].i_bgn = i;
    block[0].j_bgn = j;
    block[0].k_bgn = k;
    block[0].loops = work_per_thread[0];
	for(n = 1; n < nthreads; n++)
	{
	    for (t = 0; t < work_per_thread[n-1]; t++)
	    {
			test = get_ijk_next_prim(max,&i,&j,&k);
		}
		block[n].i_bgn = i;
	    block[n].j_bgn = j;
	    block[n].k_bgn = k;
	    block[n].loops = work_per_thread[n];
	}
}

int get_ijk_next_prim(int max, int *i1, int *j1, int *k1)
{
    int i,j,k;
	int test = 1;
	i=*i1;
	j=*j1;
	k=*k1;
	
	if(j == max){
		i++;
		j=i;
		k=j;
	}else{
		if(k == max || k>=(i+j+2)){
			j++;
			k=j;
		}else{
			k++;
		}
	}
	
	if(i>=max){
		test = 0;
		i=max;
		j=max;
		k=max;
	}
	
    *i1 = i;
    *j1 = j;
    *k1 = k;
    return test;
}

void check_ijk_min_prim(int  min, int *i, int *j, int *k)
{
	int n,i0,j0,k0;
	
	// Check points are within allowed values
	i0 = *i < min ? min : *i;
	j0 = *j < min ? min : *j;
	k0 = *k < min ? min : *k;
	
	// Check ordering of points to make sure that i<=j<=k
	if(i0>j0){
		n=j0;
		j0=i0;
		i0=n;
	}
	if(j0>k0){
		n=j0;
		j0=k0;
		k0=n;
		if(i0>j0){
			n=j0;
			j0=i0;
			i0=n;
		}
	}
	
	// Check triangle and parity conditions are met by points and move forward if not 
	if(i0+j0+2<k0){
		k0 = i0+j0+2;
	}
	
	*i = i0;
	*j = j0;
	*k = k0; 
	return;
}

void check_ijk_max_prim(int max, int *i, int *j, int *k)
{
	int n,i0,j0,k0;
	
	// Check end points are within allowed values
	i0 = *i > max ? max : *i;
	j0 = *j > max ? max : *j;
	k0 = *k > max ? max : *k;
	
	// Check ordering of points to make sure that i<=j<=k
	if(i0>j0){
		n=j0;
		j0=i0;
		i0=n;
	}
	if(j0>k0){
		n=j0;
		j0=k0;
		k0=n;
		if(i0>j0){
			n=j0;
			j0=i0;
			i0=n;
		}
	}
	
	// Check triangle and parity conditions are met by end points and move back if not
	if(i0+j0+2<k0){
		k0 = i0+j0+2;
		if(k0>max)k0=max;
	}
	
	*i = i0;
	*j = j0;
	*k = k0; 
	return;
}


void divide_tasks(long int tasks, int nranks, long int* tasks_per_rank)
{
	int size = tasks / nranks;
	int extra = tasks % nranks;
	int i;
	for(i = 0; i < extra; i++)
	{
		tasks_per_rank[i] = size+1;
	}
	for(i = extra; i < nranks; i++)
	{
		tasks_per_rank[i] = size;
	}
}
