#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <omp.h>
#include "global.h"
#include "calculate_gamma_3d.h"
//#include "offload_util.h"

/*
	After computing (1) primordial basis, (2) late-time projection of primordial basis and 
	(3) late time basis, we can relate the projected early time modes with the actual late
	time modes via the transformation matrix \Gamma:

		\Gamma_{np} = \sum_r \gamma^{-1}_{nr} * \langle Q_r, \tilde{Q}_p \rangle _l

	This function (gamma_mpi.c) computes and stores the values of \Gamma in a file,
	which can be used later in other calculations.
*/

void gamma_3d_host(double * restrict gamma_flat, int gamma_flat_size, int rank, int numranks, double workshare);		//_OFFLOADABLE

void gamma_3d_host(double * restrict gamma_flat, int gamma_size, int rank, int numranks, double workshare){				//_OFFLOADABLE
    int i, j;
    double (*restrict gamma)[gamma_size] = (double (*restrict)[gamma_size]) gamma_flat;

    int xsize_spl = get_qtilde_xsize();
    int xsize_pad = (xsize_spl + 7) & ~7;	
    double* xvec_spl = (double*) _mm_malloc(xsize_pad * sizeof(double) + 64, 64);
    get_qtilde_xvec(xvec_spl);

    // SJP: Pre-compute the xdiff and inverse xdiff for the inline spline calculation.
    double* xdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
    double* ixdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
    for (i = 0; i < xsize_spl-1; i++)
    {
        xdiff_spl[i] = xvec_spl[i+1] - xvec_spl[i];
        ixdiff_spl[i] = 1.0 / xdiff_spl[i];
    }

    // jb zero the gamma array
    int m, n;
    for (m=0; m<gamma_size; m++)
        for (n=0; n<gamma_size; n++)
            gamma[m][n] = 0.0;

    // JB thread data for CO algorithm
    thread_g3d_data_t* thread_data;
    corebarrier_t* core_barriers;

    // JB start omp parallel region here
    #pragma omp parallel default(none) private(n,m,i,j) \
    shared(rank, numranks, gamma, gamma_size, xsize_spl, xsize_pad, xvec_spl, xdiff_spl, ixdiff_spl, thread_data, core_barriers, workshare)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int maxblocksize = MAXBLOCKSIZE;

        // JB allocate the thread data for all threads
        #pragma omp barrier
        #pragma omp master
        {
            thread_data = (thread_g3d_data_t*) _mm_malloc(nthreads * sizeof(thread_g3d_data_t),64);
            core_barriers = (corebarrier_t*) _mm_malloc(nthreads * sizeof(corebarrier_t), 64);
        }
        #pragma omp barrier

        // allocate private thread data
        int gamma_size_pad = (gamma_size+7) & ~7;
        thread_data[tid].mvec = (double*) _mm_malloc(gamma_size_pad * gamma_size_pad * sizeof(double), 64);
        thread_data[tid].intgrlvec = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);
        thread_data[tid].plijkz = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);
        thread_data[tid].yvec = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);
        thread_data[tid].xsize = xsize_spl;
        thread_data[tid].xvec = xvec_spl;
        thread_data[tid].xdiff = xdiff_spl;
        thread_data[tid].ixdiff = ixdiff_spl;

        // core barriers (though not useful on the host)
        core_barriers[tid].userbarrier_arrive=0;
        core_barriers[tid].userbarrier_depart=0;
        thread_data[tid].bar = (barrier_t*) _mm_malloc(sizeof(barrier_t), 64);
        thread_data[tid].bar->me = &(core_barriers[tid]);
        thread_data[tid].bar->usersense = 1;
        thread_data[tid].bar->mycoretid = 0;
        thread_data[tid].bar->coreval = 0x00000001;    // assuming 1 threads per core

        // zero the memory
        for (m=0; m<gamma_size; m++)
            for (n=0; n<gamma_size; n++)
                thread_data[tid].mvec[m*gamma_size_pad+n]=0.0;

        #pragma omp master
        {
            decompose_gamma_TTT_3D_mpi(thread_data, rank, numranks, workshare);
        }
        #pragma omp barrier

        double t0, t1;

        // SJP: For timing accuracy, ensure that all threads start and end sync'd
        #pragma omp barrier
        t0 = omp_get_wtime();

        // main calculation
        calculate_gamma_TTT_3D_nested(thread_data,maxblocksize);

        #pragma omp barrier
        t1 = omp_get_wtime();
        #pragma omp master
        printf("[%d]\tgamma_3D \t%d\t%e\n", rank, n, t1-t0);

        // SJP: Just do this as a critical for now.  Swap for an array reduction later.
        #pragma omp barrier
        t0 = omp_get_wtime();
        #pragma omp critical
        for (m = 0; m < gamma_size; m++)
        {
            for (n = 0; n < gamma_size; n++)
            {
                gamma[m][n] += thread_data[tid].mvec[m*gamma_size_pad + n];
            }
        }
        #pragma omp barrier
        t1 = omp_get_wtime();
        #pragma omp master
        printf("[%d]\treduction\t%d\t%e\n", rank, n, t1-t0);

        /* JB free spline data */

        #pragma omp master
        {
            _mm_free(core_barriers);
        }
        _mm_free(thread_data[tid].mvec);
        _mm_free(thread_data[tid].intgrlvec);
        _mm_free(thread_data[tid].plijkz);
        _mm_free(thread_data[tid].yvec);
        _mm_free(thread_data[tid].bar);
    }
    _mm_free(xvec_spl);
    _mm_free(xdiff_spl);
    _mm_free(ixdiff_spl);
    _mm_free(thread_data);
}

int main( int argc, char *argv[] ){

	//	Set up timing
	double time1, time2, time3, time4, time5, duration;
	
	if (argc != 2) {
		printf ( "**** Incorrect number of arguments	****\n" );
		printf( "Usage is >:%s dir\n", argv[0] );
		printf ( "**** Program terminated ****\n" );
		exit (1);
	}
	
	time1 = csecond();

	//	MPI Variables and Initialization
	int rank, nproc, lnext;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);


	//	Read in parameters and Initialize
	char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);
	
	set_terms_prim();
	set_terms_late();

	int terms = get_terms_prim();
	int i,j,k,r,s,t,l,m,n;
	
	read_qtilde();

	double Tmax = (double)eflag_T_lmax;
	double Tmin = (double)eflag_T_lmin;
	int Tsize = Tmax-Tmin+1;  
	double Tvec[Tsize];
	for(i=0;i<Tsize;i++)Tvec[i] = (double)i+Tmin;
	
	double Emax = (double)eflag_E_lmax;
	double Emin = (double)eflag_E_lmin;
	int Esize = Emax-Emin+1;
	double Evec[Esize];
	for(i=0;i<Esize;i++)Evec[i] = (double)i+Emin;
	
	create_basis_late(Tsize, Esize, Tmin, Tmax, Emin, Emax, Tvec, Evec);
	
	load_cl();
	load_BN();
	load_lens();
		
	double sum1,sum2,sum3,sum4;
	double x1,x2,x3,x4;
	int psize, psize_T, psize_E;
	psize = get_pmax_prim()+1;
	psize_T = get_pmax_late_T()+1;
	if(do_polarisation==1) psize_E = get_pmax_late_E()+1;
	
	int gamma_size = terms;
	int i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4;
	
	// Helpful printing statements:
	if(rank==0){
		printf("lmax T: %d\n", (int)Tmax);
		printf("lmin T: %d\n", (int)Tmin);
		printf("fsky T: %e\n", fskyT);
		if(do_polarisation==1)printf("lmax E: %d\n", (int)Emax);
		if(do_polarisation==1)printf("lmin E: %d\n", (int)Emin);
		if(do_polarisation==1)printf("fsky E: %e\n", fskyE);
		printf("prim pmax: %d\n", get_pmax_prim());
		printf("terms: %d\n", get_terms_prim());
		printf("late T pmax: %d\n", get_pmax_late_T());
		if(do_polarisation==1)printf("late E pmax: %d\n", get_pmax_late_E());
		printf("xsize: %d\n",get_qtilde_xsize());
		if(do_polarisation==1){
			for(n=0;n<gamma_size;n++){
				find_perm_prim(n,&i,&j,&k);
				find_perm_late_TTT(n,&i1,&j1,&k1);
				find_perm_late_TTE(n,&i2,&j2,&k2);
				find_perm_late_TEE(n,&i3,&j3,&k3);
				find_perm_late_EEE(n,&i4,&j4,&k4);
				//printf("%d\t(%d,%d,%d)\t(%d,%d,%d)\t(%d,%d,%d)\t(%d,%d,%d)\t(%d,%d,%d)\n",n,i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4);
			}
		}
		else{
			for(n=0;n<gamma_size;n++){
				find_perm_prim(n,&i,&j,&k);
				find_perm_late_TTT(n,&i1,&j1,&k1);
				//printf("%d\tprimordial: (%d,%d,%d)\tlate-time: (%d,%d,%d)\n",n,i,j,k,i1,j1,k1);
			}
		}
	}
	
	int loops;
	int auxloop ;
	int start_loop;
	int end_loop;

	MPI_Barrier(MPI_COMM_WORLD);

	// Read in the data from files with previously calculated values	
	read_orthol();
	read_lambdal();

	double **orthoinv_TTT = (double **)create_array(gamma_size,gamma_size);
	double **orthoinv_TTE;
	double **orthoinv_TEE;
	double **orthoinv_EEE;
	 
	if(do_polarisation==1){
		orthoinv_TTE = (double **)create_array(gamma_size,gamma_size);
		orthoinv_TEE = (double **)create_array(gamma_size,gamma_size);
		orthoinv_EEE = (double **)create_array(gamma_size,gamma_size);
	}
	if(do_polarisation==1){
		for(i=0;i<gamma_size;i++){
			for(j=0;j<gamma_size;j++){
				orthoinv_TTT[i][j]=0.0;
				orthoinv_TTE[i][j]=0.0;
				orthoinv_TEE[i][j]=0.0;
				orthoinv_EEE[i][j]=0.0;
				for(k=0;k<gamma_size;k++){
					orthoinv_TTT[i][j] += get_lambdal_TTT(k,i)*get_lambdal_TTT(k,j);
					orthoinv_TTE[i][j] += get_lambdal_TTE(k,i)*get_lambdal_TTE(k,j);
					orthoinv_TEE[i][j] += get_lambdal_TEE(k,i)*get_lambdal_TEE(k,j);
					orthoinv_EEE[i][j] += get_lambdal_EEE(k,i)*get_lambdal_EEE(k,j);
				}
			}
		}
	}
	else{
		for(i=0;i<gamma_size;i++){
			for(j=0;j<gamma_size;j++){
				orthoinv_TTT[i][j]=0.0;
				for(k=0;k<gamma_size;k++){
					orthoinv_TTT[i][j] += get_lambdal_TTT(k,i)*get_lambdal_TTT(k,j);
				}
			}
		}
	}
	create_gamma();
	int gamma_total = terms*terms;

    double *gamma_TTT_flat = (double*) _mm_malloc(gamma_size*gamma_size * sizeof(double), 64);
    double *gamma_TTE_flat;
    double *gamma_TEE_flat;
    double *gamma_EEE_flat;
	 
	if(do_polarisation==1){
		gamma_TTE_flat = (double*) _mm_malloc(gamma_size*gamma_size * sizeof(double), 64);
		gamma_TEE_flat = (double*) _mm_malloc(gamma_size*gamma_size * sizeof(double), 64);
		gamma_EEE_flat = (double*) _mm_malloc(gamma_size*gamma_size * sizeof(double), 64);
	}

	// alias flattened arrays as multi-dim
    double (*restrict gamma_TTT)[gamma_size] = (double (*restrict)[gamma_size]) gamma_TTT_flat;
    double (*restrict gamma_TTE)[gamma_size] = (double (*restrict)[gamma_size]) gamma_TTE_flat;
    double (*restrict gamma_TEE)[gamma_size] = (double (*restrict)[gamma_size]) gamma_TEE_flat;
    double (*restrict gamma_EEE)[gamma_size] = (double (*restrict)[gamma_size]) gamma_EEE_flat;
	
	if(rflag_do3D == 0){
		
		init_gamma_3Dint();
		
		loops=gamma_size/nproc;
		auxloop = fmod(gamma_size,nproc);
		start_loop = rank*loops;
		end_loop = (rank+1)*loops-1;

		if (auxloop != 0){
			if (rank < auxloop){
				start_loop = start_loop + rank;
				end_loop = end_loop + rank + 1;
			}else{
				start_loop = start_loop + auxloop;
				end_loop = end_loop + auxloop;
			}
		}
		
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				gamma_TTT[i][j] = 0e0;
				if (do_polarisation)
                {
                    gamma_TTE[i][j] = 0e0;
                    gamma_TEE[i][j] = 0e0;
                    gamma_EEE[i][j] = 0e0;
                }
			}
		}
 		
        // main gamma 3d calculation.
        gamma_3d_host(&(gamma_TTT[0][0]), gamma_size, rank, nproc, 0.0);
		
		// polarisation currently disabled for 3d
        if (do_polarisation)
        {
            perror("WARNING: Polarisation currently disabled for 3D gamma calculations.");
            perror("WARNING: Skipping polarisation calculation...");
        }
 		/*
		printf("3D mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
		MPI_Barrier(MPI_COMM_WORLD);
		
		double *mvec = (double *)create_vector(gamma_size);
		
		time1 = csecond();
		for(n=start_loop;n<end_loop+1;n++){
			calculate_gamma_TTT_3D(n,mvec);
			for (r=0;r<gamma_size;r++)gamma_TTT[r][n] = mvec[r];
			calculate_gamma_TTE_3D(n,mvec);
			for (r=0;r<gamma_size;r++)gamma_TTE[r][n] = mvec[r];
			calculate_gamma_TEE_3D(n,mvec);
			for (r=0;r<gamma_size;r++)gamma_TEE[r][n] = mvec[r];
			calculate_gamma_EEE_3D(n,mvec);
			for (r=0;r<gamma_size;r++)gamma_EEE[r][n] = mvec[r];
			
			duration = csecond() - time1;
			printf("[%d]\t3D QQ\t%d\t%e\n",rank,n, duration);
			time1 = csecond();
			
		}
		printf("[%d] finished\n", rank);
		*/
		
	}
	else if(eflag_order_late==6){
		
		if(do_polarisation==1){
			init_gamma_glint();
			int Mij_total;
			if(do_polarisation==1){
				Mij_total = psize*(psize_T+psize_E);
			}else{
				Mij_total = psize*psize_T;
			}
			
			int Mij_pairs[Mij_total][3];
			n=0;
			for (i=0;i<psize;i++) {
				for (j=0;j<psize_T;j++) {
					Mij_pairs[n][0] = 0;
					Mij_pairs[n][1] = i;
					Mij_pairs[n][2] = j;
					n++;
				}
				if(do_polarisation==1){
					for (j=0;j<psize_E;j++) {
						Mij_pairs[n][0] = 1;
						Mij_pairs[n][1] = i;
						Mij_pairs[n][2] = j;
						n++;
					}
				}
			}
			
			loops=Mij_total/nproc;
			auxloop = fmod(Mij_total,nproc);
			start_loop = rank*loops;
			end_loop = (rank+1)*loops-1;
			
			if (auxloop != 0){
				if (rank < auxloop){
					start_loop = start_loop + rank;
					end_loop = end_loop + rank + 1;
				}
				else{
					start_loop = start_loop + auxloop;
					end_loop = end_loop + auxloop;
				}
			}
			
			time1 = csecond();
			for(n=start_loop;n<end_loop+1;n++){
				
				i = Mij_pairs[n][1];
				j = Mij_pairs[n][2];
				if(Mij_pairs[n][0]==0){
					calculate_gammaMij_T(i,j);
				}else{
					calculate_gammaMij_E(i,j);
				}
				duration = csecond() - time1;
				time1 = csecond();
			}
			
			sync_gammaMij();
			
			int gamma_pairs[gamma_total][2];
			n=0;
			for (i=0;i<gamma_size;i++) {
				for (j=0;j<gamma_size;j++) {
					gamma_pairs[n][0] = i;
					gamma_pairs[n][1] = j;
					n++;
				}
			}
			
			loops=gamma_total/nproc;
			auxloop = fmod(gamma_total,nproc);
			start_loop = rank*loops;
			end_loop = (rank+1)*loops-1;
			
			if (auxloop != 0){
				if (rank < auxloop){
					start_loop = start_loop + rank;
					end_loop = end_loop + rank + 1;
				}else{
					start_loop = start_loop + auxloop;
					end_loop = end_loop + auxloop;
				}
			}

			for (i=0;i<gamma_size;i++) {
				for (j=0;j<gamma_size;j++) {
					gamma_TTE[i][j] = 0e0;
					gamma_TEE[i][j] = 0e0;
					gamma_EEE[i][j] = 0e0;
				}
			}
			
			double* result;
			
			result = (double *)malloc( 4 * sizeof(double) );
			
			time1 = csecond();
			if(do_polarisation==1){
				for(n=start_loop;n<end_loop+1;n++){
					
					i = gamma_pairs[n][0];
					j = gamma_pairs[n][1];
					
					calculate_gamma(i,j,result);
					gamma_TTE[i][j] = result[1];
					gamma_TEE[i][j] = result[2];
					gamma_EEE[i][j] = result[3];
					
					duration = csecond() - time1;
					//printf("[%d]\t%d\t%d\t%e\t%e\t%e\t%e\n", rank, i, j, gamma_TTE[i][j], gamma_TEE[i][j], gamma_EEE[i][j], duration);
					time1 = csecond();
					
				} // end of MPI loop
			}
			printf("[%d] finished\n", rank);
		}
			
		init_gamma_3Dint();


		/*
		loops=gamma_size/nproc;
		auxloop = fmod(gamma_size,nproc);
		start_loop = rank*loops;
		end_loop = (rank+1)*loops-1;
			
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				gamma_TTT[i][j] = 0e0;
			}
		}
		
		if (auxloop != 0){
			if (rank < auxloop){
				start_loop = start_loop + rank;
				end_loop = end_loop + rank + 1;
			}else{
				start_loop = start_loop + auxloop;
				end_loop = end_loop + auxloop;
			}
		}
 		
		printf("3D mpirank %d auxloop %d start %d end %d\n",rank,auxloop,start_loop,end_loop);
		MPI_Barrier(MPI_COMM_WORLD);
		
		double *mvec = (double *)create_vector(gamma_size);
		
		time1 = csecond();
		for(n=start_loop;n<end_loop+1;n++){
			calculate_gamma_TTT_3D(n,mvec);
			for (r=0;r<gamma_size;r++)gamma_TTT[r][i] += mvec[r];
			
			duration = csecond() - time1;
			printf("[%d]\t3D QQ\t%d\t%d\t%d\t%e\n",rank,n, i, j, duration);
			time1 = csecond();
			
		}
		*/
	}
	else{
		
		double t0 = MPI_Wtime();
		init_gamma_glint();
		double t1 = MPI_Wtime();

		int Mij_total;

		if(do_polarisation==1){
			Mij_total = psize*(psize_T+psize_E);
		}
		else{
			Mij_total = psize*psize_T;
		}
		
		int Mij_pairs[Mij_total][3];
		n=0;
		for (i=0;i<psize;i++){
			for (j=0;j<psize_T;j++){
				Mij_pairs[n][0] = 0;
				Mij_pairs[n][1] = i;
				Mij_pairs[n][2] = j;
				n++;
			}
			if(do_polarisation==1){
				for (j=0;j<psize_E;j++) {
					Mij_pairs[n][0] = 1;
					Mij_pairs[n][1] = i;
					Mij_pairs[n][2] = j;
					n++;
				}
			}
		}
		
		loops=Mij_total/nproc;
		auxloop = fmod(Mij_total,nproc);
		start_loop = rank*loops;
		end_loop = (rank+1)*loops-1;
		
		if (auxloop != 0){
			if (rank < auxloop){
				start_loop = start_loop + rank;
				end_loop = end_loop + rank + 1;
			}
			else{
				start_loop = start_loop + auxloop;
				end_loop = end_loop + auxloop;
			}
		}
		
		for (i=0;i<gamma_size;i++){
			for (j=0;j<gamma_size;j++){
				gamma_TTT[i][j] = 0e0;
				if(do_polarisation==1){
					gamma_TTE[i][j] = 0e0;
					gamma_TEE[i][j] = 0e0;
                    gamma_EEE[i][j] = 0e0;
                }
            }
        }
		
		time1 = csecond();
		t0 = MPI_Wtime();

        /*#pragma offload target(mic:offload_target)
        {
            precompute_gammaMij_T();

            if (do_polarisation)
                precompute_gammaMij_E();
        }*/

		precompute_gammaMij_T();
		if (do_polarisation==1) precompute_gammaMij_E();

		MPI_Barrier(MPI_COMM_WORLD);
		printf("[%d] Finished Mij\n",rank);
		t1 = MPI_Wtime();

	    t0 = MPI_Wtime();

		printf("[%d] Syncronised Mij\n",rank);
		t1 = MPI_Wtime();

		MPI_Barrier(MPI_COMM_WORLD);
		
		loops=gamma_total/nproc;
		auxloop = fmod(gamma_total,nproc);
		start_loop = rank*loops;
		end_loop = (rank+1)*loops-1;
		
		if (auxloop != 0){
			if (rank < auxloop){
				start_loop = start_loop + rank;
				end_loop = end_loop + rank + 1;
			}else{
				start_loop = start_loop + auxloop;
				end_loop = end_loop + auxloop;
			}
		}

        time1 = csecond();
        t0 = MPI_Wtime();

//============================================================================================================================
        #pragma offload_transfer target(mic:offload_target) \
                in(gamma_TTT_flat[0:gamma_total] : ALLOC RETAIN)

        if (do_polarisation)
        {
            #pragma offload_transfer target(mic:offload_target) \
                    in(gamma_TTE_flat[0:gamma_total] : ALLOC RETAIN) \
                    in(gamma_TEE_flat[0:gamma_total] : ALLOC RETAIN) \
                    in(gamma_EEE_flat[0:gamma_total] : ALLOC RETAIN)
        }
        #pragma offload target(mic:offload_target) in(start_loop, end_loop, gamma_size, do_polarisation, deltaphi) \
            nocopy(gamma_TTT_flat, gamma_TTE_flat, gamma_TEE_flat, gamma_EEE_flat)
        {

            #pragma omp parallel default(none) shared(start_loop,end_loop,\
                    gamma_TTT_flat,\
                    gamma_TTE_flat, gamma_TEE_flat,\
                    gamma_EEE_flat, do_polarisation, \
                    gamma_size) \
                    private(n,i,j)
            {
                double result[4];

                #pragma omp for
                for(n=start_loop;n<end_loop+1;n++)
                {
                    i = n / gamma_size;
                    j = n % gamma_size;

                    calculate_gamma_2d(i,j,&result[0]);

                    gamma_TTT_flat[n] = result[0];
					//printf("Gamma_TTT_flat = %e\n", gamma_TTT_flat[n]);
                    if (do_polarisation)
                    {
                        gamma_TTE_flat[n] = result[1];
                        gamma_TEE_flat[n] = result[2];
                        gamma_EEE_flat[n] = result[3];
                    }
                }
            }
        }
        #pragma offload_transfer target(mic:offload_target) \
                out(gamma_TTT_flat[0:gamma_total] : REUSE FREE)

        if (do_polarisation)
        {
            #pragma offload_transfer target(mic:offload_target) \
                    out(gamma_TTE_flat[0:gamma_total] : REUSE FREE) \
                    out(gamma_TEE_flat[0:gamma_total] : REUSE FREE) \
                    out(gamma_EEE_flat[0:gamma_total] : REUSE FREE)
        }
		printf("[%d] finished\n", rank);
		t1 = MPI_Wtime();
		//if (rank==0)
		   // printf("[%d] Mij compute time = %f\n", rank, t1-t0);
	}

//============================================================================================================================

	//sync_gammaMij();  --- this seems to be equivalent to what is happening with send and recv below

	double* gamma_TTT_send;
	double* gamma_TTT_recv;
	double* gamma_TTE_send;
	double* gamma_TTE_recv;
	double* gamma_TEE_send;
	double* gamma_TEE_recv;
	double* gamma_EEE_send;
	double* gamma_EEE_recv;
	
	gamma_TTT_send = (double *)malloc( gamma_total * sizeof(double) );
	gamma_TTT_recv = (double *)malloc( gamma_total * sizeof(double) );
	if(do_polarisation==1){
		gamma_TTE_send = (double *)malloc( gamma_total * sizeof(double) );
		gamma_TTE_recv = (double *)malloc( gamma_total * sizeof(double) );
		gamma_TEE_send = (double *)malloc( gamma_total * sizeof(double) );
		gamma_TEE_recv = (double *)malloc( gamma_total * sizeof(double) );
		gamma_EEE_send = (double *)malloc( gamma_total * sizeof(double) );
		gamma_EEE_recv = (double *)malloc( gamma_total * sizeof(double) );
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	n=0;
	for (i=0;i<gamma_size;i++) {
		for (j=0;j<gamma_size;j++) {
			gamma_TTT_send[n] = gamma_TTT[i][j];
			n++;
		}
	}
	if(do_polarisation==1){
		n=0;
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				gamma_TTE_send[n] = gamma_TTE[i][j];
				n++;
			}
		}
		n=0;
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				gamma_TEE_send[n] = gamma_TEE[i][j];
				n++;
			}
		}
		n=0;
		for (i=0;i<gamma_size;i++) {
			for (j=0;j<gamma_size;j++) {
				gamma_EEE_send[n] = gamma_EEE[i][j];
				n++;
			}
		}
	}
	
	MPI_Reduce(gamma_TTT_send,gamma_TTT_recv,gamma_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
	if(do_polarisation==1){
		MPI_Reduce(gamma_TTE_send,gamma_TTE_recv,gamma_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(gamma_TEE_send,gamma_TEE_recv,gamma_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(gamma_EEE_send,gamma_EEE_recv,gamma_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		n=0;
		for (i=0;i<gamma_size;i++){
			for (j=0;j<gamma_size;j++){
				
				gamma_TTT[i][j] = gamma_TTT_recv[n];
// 				if(i!=j)gamma_TTT[j][i] = gamma_TTT[i][j];
				
				if(do_polarisation==1){
					gamma_TTE[i][j] = gamma_TTE_recv[n];
// 					if(i!=j)gamma_TTE[j][i] = gamma_TTE[i][j];
					gamma_TEE[i][j] = gamma_TEE_recv[n];
// 					if(i!=j)gamma_TEE[j][i] = gamma_TEE[i][j];
					gamma_EEE[i][j] = gamma_EEE_recv[n];
// 					if(i!=j)gamma_EEE[j][i] = gamma_EEE[i][j];
				}
				
				n++;
			}
		}
		printf("\nOutput gamma calc\n");
        for (i=0;i<gamma_size;i++){
             for (j=0;j<gamma_size;j++) {
                 printf("%e\n",gamma_TTT[i][j]);
             }
             //printf("\n");
        }
			
		double **gamma_TTT_temp;
		double **gamma_TTE_temp;
		double **gamma_TEE_temp;
		double **gamma_EEE_temp;
	 
		if(do_polarisation==1){
			gamma_TTT_temp = (double **)create_array(gamma_size,gamma_size);
			gamma_TTE_temp = (double **)create_array(gamma_size,gamma_size);
			gamma_TEE_temp = (double **)create_array(gamma_size,gamma_size);
			gamma_EEE_temp = (double **)create_array(gamma_size,gamma_size);
			for (i=0;i<gamma_size;i++) {
				for (j=0;j<gamma_size;j++) {
					gamma_TTT_temp[i][j] = 0e0;
					gamma_TTE_temp[i][j] = 0e0;
					gamma_TEE_temp[i][j] = 0e0;
					gamma_EEE_temp[i][j] = 0e0;
					for (n=0;n<gamma_size;n++) {
						gamma_TTT_temp[i][j] +=orthoinv_TTT[i][n]*gamma_TTT[n][j];
						gamma_TTE_temp[i][j] +=orthoinv_TTE[i][n]*gamma_TTE[n][j];
						gamma_TEE_temp[i][j] +=orthoinv_TEE[i][n]*gamma_TEE[n][j];
						gamma_EEE_temp[i][j] +=orthoinv_EEE[i][n]*gamma_EEE[n][j];
						printf("Gamma_TTT_temp = %e\n",gamma_TTT_temp[i][j]);
						printf("Gamma_TTE_temp = %e\n",gamma_TTE_temp[i][j]);
						printf("Gamma_TEE_temp = %e\n",gamma_TEE_temp[i][j]);
						printf("Gamma_EEE_temp = %e\n",gamma_EEE_temp[i][j]);
					}
				}
			}
		}else{
			gamma_TTT_temp = (double **)create_array(gamma_size,gamma_size);
			for (i=0;i<gamma_size;i++){
				for (j=0;j<gamma_size;j++){
					gamma_TTT_temp[i][j] = 0e0;
					for (n=0;n<gamma_size;n++) {
						gamma_TTT_temp[i][j] +=orthoinv_TTT[i][n]*gamma_TTT[n][j];
						//printf("Gamma_TTT_temp = %e\n",gamma_TTT_temp[i][j]);
					}
				}
			}
		}
		
		 //printf("\n-------------------------------------------------------------\n");
		 for (j=0;j<gamma_size;j++) {
		 	//printf("Gamma_TTT_temp = %e\n",gamma_TTT_temp[j][0]);
		 }
		 //printf("\n");
		
			
		double *results_g;
		if(do_polarisation==1){
			results_g =  malloc( sizeof(double)*6);
			results_g[0] = 0.0;
			results_g[1] = 0.0;
			results_g[2] = 0.0;
			results_g[3] = 0.0;
			results_g[4] = 0.0;
			results_g[5] = 0.0;
		}else{
			results_g =  malloc( sizeof(double)*3);
			results_g[0] = 0.0;
			results_g[1] = 0.0;
			results_g[2] = 0.0;
		}
		
		if(do_polarisation==1){
			for (i=0;i<gamma_size;i++) {
				for (j=0;j<gamma_size;j++) {
					results_g[0] = (double)i;
					results_g[1] = (double)j;
					results_g[2] = gamma_TTT_temp[i][j];
					results_g[3] = gamma_TTE_temp[i][j];
					results_g[4] = gamma_TEE_temp[i][j];
					results_g[5] = gamma_EEE_temp[i][j];
					update_gamma(results_g);
				}
			}
		}else{
			for (i=0;i<gamma_size;i++) {
				for (j=0;j<gamma_size;j++) {
					results_g[0] = (double)i;
					results_g[1] = (double)j;
					results_g[2] = gamma_TTT_temp[i][j];
					update_gamma(results_g);
				}
			}
		}
		output_gamma();
		printf("[%d] Written Gamma\n", rank);
		
		printf("\n");

		for (n=0;n<gamma_size;n++) {
			for (i=0;i<gamma_size;i++) {
				gamma_TTT[i][n] = 0.0;
				for (j=i;j<gamma_size;j++) {
					gamma_TTT[i][n] += gamma_TTT_temp[j][n]*get_orthol_TTT(j,i);
				}
				printf("%e\t",gamma_TTT[i][n]);
			}
			//printf("\n");
		}

		/*
		// #ifdef DEBUG
		//         printf("\n");
		//         for (i=0; i<gamma_size; i++)
		//         {
		//             for (j=0;j<gamma_size;j++) {
		//                 printf("%e ",gamma_TTT[i][j]);
		//             }
		//             printf("\n");
		//         }
		//         if (do_polarisation)
		//         {
		//             printf("\n");
		//             for (i=0; i<gamma_size; i++)
		//             {
		//                 for (j=0;j<gamma_size;j++) {
		//                     printf("%e ",gamma_TTE[i][j]);
		//                 }
		//                 printf("\n");
		//             }
		//             printf("\n");
		//             for (i=0; i<gamma_size; i++)
		//             {
		//                 for (j=0;j<gamma_size;j++) {
		//                     printf("%e ",gamma_TEE[i][j]);
		//                 }
		//                 printf("\n");
		//             }
		//             printf("\n");
		//             for (i=0; i<gamma_size; i++)
		//             {
		//                 for (j=0;j<gamma_size;j++) {
		//                     printf("%e ",gamma_EEE[i][j]);
		//                 }
		//                 printf("\n");
		//             }
		//         }
		//         #else
        // printf("\n");
        // for (i=0; i<gamma_size; i++)
        // {
        //     printf("%e\n",gamma_TTT[i][0]);
        // }
        // printf("\n");
        // #endif
		
		// for (i=0;i<gamma_size;i++) {
		// 	for (j=0;j<gamma_size;j++) {
		// 		x1=x2=x3=0.0;
		// 		for (n=0;n<gamma_size;n++) {
		// 			x1 += gamma_TTT[n][i]*gamma_TTT[n][i];
		// 			x2 += gamma_TTT[n][j]*gamma_TTT[n][j];
		// 			x3 += gamma_TTT[n][i]*gamma_TTT[n][j];
		// 		}
		// 		x4 = x3 / sqrt(x1*x2);
		// 		//if(j==0)printf("%e\t",x4);
		// 		printf("%e\t",x4);
		// 	}
		// 	printf("\n");
		// }
		*/
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

