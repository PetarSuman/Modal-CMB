#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <mkl.h>
#include <omp.h>

#include "global.h"
#include "offload_util.h"
#include "local_barrier.h"
#include "calculate_gamma_3d.h"

static double *lweight3D;
static double ***qtilde3D_T;
static double ***qtilde3D_E;

static double calculate_sijk_TTT(int n, int xsize, double* restrict x, double* restrict y, int l1, int l2, int l3, double* restrict xdiff, double* restrict ixdiff);
static double calculate_sijk_TTE(int n, int l1, int l2, int l3);
static double calculate_sijk_TEE(int n, int l1, int l2, int l3);
static double calculate_sijk_EEE(int n, int l1, int l2, int l3);

_OFFLOADABLE static const double sixth = (1.0/6.0);

void init_gamma_3Dint()
{
    double s1,s2,s3,s4,s5,s6,s7,s8,s9;	
    int i,j,l;

    int xsize = get_qtilde_xsize();
    int lsize_T = get_qtilde_T_lsize();
    int *lvec_T = create_ivector(lsize_T);
    get_qtilde_T_lvec(lvec_T);

    lweight3D = (double *)create_vector(lsize_T);

    for(i = 0; i < lsize_T; i++){
        l = lvec_T[i];
        if(l > 1){
            s1 = pow(2e0*l+1e0,1e0/3e0);
            lweight3D[i] =  1e0/s1;
        } else {
            lweight3D[i] = 0e0;
        }
    }

    int pmax = get_pmax_prim();
    qtilde3D_T = (double ***)create_3Darray(lsize_T,pmax+1,xsize);
    for(l=0;l<lsize_T;l++){
        if(l>1){
            s1 = pow(2e0*l+1e0,1e0/6e0);
            s2 = get_cl_TT(l)+get_noise_TT(l)/(get_beam_TT(l)*get_beam_TT(l));
            s3 = 1e0 / sqrt(s2);
            s4 =  s1*s3;
        }else{
            s4 = 0e0;
        }
        for(i = 0; i < xsize; i++){
            for(j = 0; j < pmax+1; j++){
                qtilde3D_T[l][j][i] = s4*get_qtilde_T(j,l,i);
            }
        }
    }

    if(do_polarisation==1){
        int lsize_E = get_qtilde_E_lsize();
        int *lvec_E = create_ivector(lsize_E);
        get_qtilde_E_lvec(lvec_E);

        qtilde3D_E = (double ***)create_3Darray(lsize_E,pmax+1,xsize);
        for(l = 0; l < lsize_E; l++){
            if(l>1){
                s1 = pow(2e0*l+1e0,1e0/6e0);
                s2 =get_beam_TT(l)*get_beam_TT(l)*get_cl_TT(l)+get_noise_TT(l);
                s3 =get_beam_TT(l)*get_beam_EE(l)*get_cl_TE(l)+get_noise_TE(l);
                s4 =get_beam_EE(l)*get_beam_EE(l)*get_cl_EE(l)+get_noise_EE(l);
                s5 = get_beam_EE(l)/sqrt(s2*(s2*s4-s3*s3));
                s4 =  s1*s3;
            }else{
                s4 = 0e0;
            }
            for(i = 0; i < xsize; i++){
                for(j = 0; j < pmax+1; j++){
                    qtilde3D_E[l][j][i] = s5*(s2*get_qtilde_E(j,l,i)-s3*get_qtilde_T(j,l,i));
                }
            }
        }
    }
    destroy_qtilde();

    return;
}

void decompose_gamma_TTT_3D_mpi(thread_g3d_data_t* thread_data, int rank, int numranks, double workshare) 
{
    int i, j, k;
    //int lmax = get_lmax();
    int lsize_T = get_qtilde_T_lsize();

    // Step 1) Count the number of iterations in the original loop.
    int count = 0;
    for (i = 2; i < lsize_T; i++) {
        for (j = i; j < lsize_T; j++) {
            for (k = j + (i % 2); k < min(i+j,(lsize_T-1))+1; k += 2) {
                count++;
            }
        }
    }

    // Step 2) Divvy up the total iterations between the MPI ranks
    // number of iterations for this task
    int blk_size_raw = count / numranks;
    int left_over = count - (blk_size_raw * numranks);

    // Step 3) Find the start and end points of the iterations for this MPI rank
    int r;
    int i0, j0, k0;
    int i1, j1, k1;
    int dummy;

    i0 = 2;
    j0 = i0;
    k0 = j0;
    // iterate through the preceding ranks up to and excluding this rank
    for (r = 0; r < rank; r++)
    {
        int rankwork = blk_size_raw;
        if (r < left_over) rankwork += 1;

        get_next_ijk_bound(i0, j0, k0, &i1, &j1, &k1, lsize_T, lsize_T, lsize_T, &dummy, rankwork);
        i0 = i1;
        j0 = j1;
        k0 = k1;
    }

    // get the bounds for this MPI rank
    // give raw amount of work to every worker in every task
    int work = 0;
    work += blk_size_raw;

    if (rank < left_over)
        work += 1;

    // Get bounds for the host
    int host_work = (int) ( (1.0 - workshare)*work );
    get_next_ijk_bound(i0, j0, k0, &i1, &j1, &k1, lsize_T, lsize_T, lsize_T, &dummy, host_work);

    // if on device then get the bounds and work count for the device
    #ifdef __MIC__
    i0 = i1;
    j0 = j1;
    k0 = k1;
    int mic_work = work - host_work;
    get_next_ijk_bound(i0, j0, k0, &i1, &j1, &k1, lsize_T, lsize_T, lsize_T, &dummy, mic_work);
    work = mic_work;

    #else
    work = host_work;
    #endif

    // Step 4) Divvy up the iterations for this MPI rank between the OMP threads
    int nthreads = omp_get_num_threads();

    // hack. host has no nested parallelism
    #ifdef __MIC__
    const int ntasks = MIC_NCORES;
    #else
    const int ntasks = nthreads;
    #endif
    int nworkers = nthreads / ntasks;

    int thrd_work = 0;
    int rank_count = work;

    int rank_blk_size_raw = rank_count / ntasks;
    int rank_left_over = rank_count - (rank_blk_size_raw * ntasks);

    // give raw amount of work to every worker in every task
    thrd_work += rank_blk_size_raw;

    // go through all the threads and assign them their upper and lower bounds
    int final_i = i1;
    int final_j = j1;
    int final_k = k1;

    int t;
    for (t = 0; t < nthreads; t+=nworkers)
    {
        thrd_work = rank_blk_size_raw;
        int taskid = t / nworkers;

        if (taskid < rank_left_over) thrd_work += 1;
        get_next_ijk_bound(i0, j0, k0, &i1, &j1, &k1, final_i, final_j, final_k, &dummy, thrd_work);

        #ifdef DEBUG
        printf("[%d] tid=%d: -- start = (%d, %d, %d), end = (%d, %d, %d), work = %d\n",rank, t, i0, j0, k0, i1, j1, k1, thrd_work);
        #endif

        int wkr;
        for (wkr = 0; wkr < nworkers; wkr++)
        {
            thread_data[t+wkr].i0 = i0;
            thread_data[t+wkr].j0 = j0;
            thread_data[t+wkr].k0 = k0;
            thread_data[t+wkr].i1 = i1;
            thread_data[t+wkr].j1 = j1;
            thread_data[t+wkr].k1 = k1;
        }
        i0 = i1;
        j0 = j1;
        k0 = k1;
    }
}

void get_next_ijk_bound(int i0, int j0, int k0, int *i1, int *j1, int *k1, \
        int final_i, int final_j, int final_k, int* blocksize, int maxblocksize)
{
    int lsize_T = get_qtilde_T_lsize();
    int i,j,k;
    int _blocksize=0;

    int iout, jout, kout;
    iout = final_i;
    jout = final_j;
    kout = final_k;

    int start_i = i0;
    int end_i = final_i + 1;
    for (i=start_i; i<end_i; i++)
    {
        int start_j = i == i0 ? j0 : i;
        int end_j = i == final_i ? final_j+1 : lsize_T;
        for (j=start_j; j<end_j; j++)
        {
            int start_k = j+(i%2);
            int end_k = min(i+j,(lsize_T-1))+1;
            if ((i==i0) && (j==j0))           
            {
                start_k = k0;
            }
            if ((i==final_i) && (j==final_j))
            {
                end_k = final_k;
            }

            for (k = start_k; k < end_k; k += 2) 
            {
                if (_blocksize == maxblocksize)
                {
                    iout = i;
                    jout = j;
                    kout = k;
                    goto end;
                }
                _blocksize++;
            }
        }
    }
end:
    // cap the output at the final value
    if ( (iout>=final_i) && (jout>=final_j) && (kout>final_k) )
    {
        iout = final_i;
        jout = final_j;
        kout = final_k;
    }

    // output the last i,j,k values and the block size
    *i1 = iout;
    *j1 = jout;
    *k1 = kout;
    *blocksize = _blocksize;
    return;
}

void calculate_gamma_TTT_3D_nested(thread_g3d_data_t* thread_data, int maxblocksize)
{
    int i,j,k,n,m,s,t1,t2,t3;
    int terms = get_terms_prim();
    int pmax_late_T = get_pmax_late_T();
    int terms_pad = (terms + 7) & ~7;
    int lsize_T = get_qtilde_T_lsize();

    // get the thread, task and worker ids
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    // hack. host has no nested parallelism
    #ifdef __MIC__
    const int ntasks = MIC_NCORES;
    #else
    const int ntasks = nthreads;
    #endif
    int nworkers = nthreads / ntasks;
    int taskid = tid / nworkers;
    int workerid = tid % nworkers;

    int xsize = thread_data[tid].xsize;
    double* xvec = thread_data[tid].xvec;
    double* yvec = thread_data[tid].yvec;
    double (*restrict mvec)[terms_pad] = (double (*restrict)[terms_pad]) thread_data[tid].mvec;
    double* xdiff = thread_data[tid].xdiff;
    double* ixdiff = thread_data[tid].ixdiff;
    double (*restrict x)[terms_pad] = (double (*restrict)[terms_pad]) thread_data[tid].intgrlvec;
    double (*restrict plijkz)[terms_pad] = (double (*restrict )[terms_pad]) thread_data[tid].plijkz;
    double xmax = xvec[xsize-1];

    double* restrict basis_late_flat = get_basis_late_T_flat();
    double (*restrict basis_late)[pmax_late_T+1] = (double (*restrict)[pmax_late_T+1]) basis_late_flat;

    // store the ultimate loop bounds for this task
    const int initial_i = thread_data[tid].i0;
    const int initial_j = thread_data[tid].j0;
    const int initial_k = thread_data[tid].k0;
    const int final_i = thread_data[tid].i1;
    const int final_j = thread_data[tid].j1;
    const int final_k = thread_data[tid].k1;

    // i0,i1,j0,j1...etc are the loop bounds for a subblock within the task
    // initialise subblock lower bound to lower task bound
    int i0 = initial_i;
    int j0 = initial_j;
    int k0 = initial_k;
    int i1;
    int j1;
    int k1;
    int blocksize;

    // get upper bounds of subblock
    get_next_ijk_bound(i0,j0,k0,&i1,&j1,&k1,final_i,final_j,final_k,&blocksize,maxblocksize);
    int hasWork = blocksize > 0;

    // main loop
    while (hasWork)
    {
        double t0 = omp_get_wtime();

        int start_i = i0;
        int end_i = i1+1;
        int l;

        // do blocksize iterations of loop nest
        l=0;
        for (i = start_i; i < end_i; i++)
        {
            int start_j = i == i0 ? j0 : i;
            int end_j = i == i1 ? j1+1 : lsize_T;
            for (j = start_j; j < end_j; j++)
            {
                int start_k = j+(i%2);
                int end_k = min(i+j,(lsize_T-1))+1;
                if ((i==i0) && (j==j0))
                {
                    start_k = k0;
                }
                if ((i==i1) && (j==j1))
                {
                    end_k = k1;
                }

                for (k = start_k; k < end_k; k += 2) 
                {
                    int start_n = workerid;
                    int end_n = terms;
                    for (n = start_n; n < end_n; n+=nworkers)
                    { 
                        x[l][n] = calculate_sijk_TTT(n, xsize, xvec, yvec, i, j, k, xdiff, ixdiff);
                    }
                    l++;
                }
            }
        }

        l=0;
        for (i = start_i; i < end_i; i++)
        {
            int start_j = i == i0 ? j0 : i;
            int end_j = i == i1 ? j1+1 : lsize_T;
            for (j = start_j; j < end_j; j++)
            {
                int start_k = j+(i%2);
                int end_k = min(i+j,(lsize_T-1))+1;
                if ((i==i0) && (j==j0))
                {
                    start_k = k0;
                }
                if ((i==i1) && (j==j1))
                {
                    end_k = k1;
                }

                for (k = start_k; k < end_k; k += 2) 
                {
                    double s1 = lweight3D[i];
                    double s2 = lweight3D[j];
                    double s3 = lweight3D[k];
                    double z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;

                    int start_m = workerid;
                    int end_m = terms;
                    for (m = start_m; m < end_m; m+=nworkers)
                    {
                        plijkz[l][m] = plijk_TTT(m,i,j,k) * z;
                    }
                    l++;
                }
            }
        }
        userCoreBarrier(thread_data[tid].bar);

        int start_m = workerid;
        int end_m = terms;

        // one worker matrix multiply x and plijkz and store in mvec; other workers wait.
        if (workerid == 0)
        {
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, terms, terms, blocksize, 1.0, &plijkz[0][0], terms_pad, &x[0][0], terms_pad, 1.0, &mvec[0][0], terms_pad);
        }
        userCoreBarrier(thread_data[tid].bar);

        // set new lower block bounds to old upper block bounds
        // and get next upper bounds
        i0 = i1;
        j0 = j1;
        k0 = k1;
        get_next_ijk_bound(i0,j0,k0,&i1,&j1,&k1,final_i,final_j,final_k,&blocksize,maxblocksize);
        hasWork = blocksize > 0;    // stopping condition

        userCoreBarrier(thread_data[tid].bar);
    }
    userCoreBarrier(thread_data[tid].bar);

    int start_m = workerid;
    int end_m = terms;
    for (m = start_m; m < end_m; m+=nworkers)
    {
        //#pragma vector aligned
        for (n = 0; n < terms; n++)
        {
            mvec[m][n] *= 6.0*deltaphi*deltaphi;
        }
    }
}


void calculate_gamma_TTT_3D(int n, double *mvec) {

    int i,j,k,m,t1,t2,t3;
    double x,y,z;
    double b1,b2,b3,b4;
    double q1,q2,q3,q4;

    int lsize_T = get_qtilde_T_lsize();
    int terms = get_terms_prim();
    double s1,s2,s3;

    for(m=0;m<terms;m++) mvec[m] = 0e0;

    for(i=2;i<lsize_T;i+=2){
        s1 = lweight3D[i];
        for(j=i;j<lsize_T;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_T)t1=lsize_T;
            for(k=j;k<t1;k+=2){
                //x = calculate_sijk_TTT(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TTT(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
        for(j=i+1;j<lsize_T;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_T)t1=lsize_T;
            for(k=j;k<t1;k+=2){
                //x = calculate_sijk_TTT(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TTT(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }

    for(i=3;i<lsize_T;i+=2){
        s1 = lweight3D[i];
        for(j=i;j<lsize_T;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_T)t1=lsize_T;
            for(k=j+1;k<t1;k+=2){
                //x = calculate_sijk_TTT(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TTT(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
        for(j=i+1;j<lsize_T;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_T)t1=lsize_T;
            for(k=j+1;k<t1;k+=2){
                //x = calculate_sijk_TTT(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TTT(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }
    for(m=0;m<terms;m++) mvec[m] *= 6.0*deltaphi*deltaphi;

    return;
}

void calculate_gamma_TTE_3D(int n, double *mvec) {

    int i,j,k,m,t1,t2,t3;
    double x,y,z;
    double b1,b2,b3,b4;
    double q1,q2,q3,q4;

    int lsize_T = get_qtilde_T_lsize();
    int lsize_E = get_qtilde_E_lsize();
    int terms = get_terms_prim();
    double s1,s2,s3;

    for(m=0;m<terms;m++) mvec[m] = 0e0;

    for(i=2;i<lsize_T;i+=2){
        s1 = lweight3D[i];
        for(j=i;j<lsize_T;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j;k<t1;k+=2){
                x = calculate_sijk_TTE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TTE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
        for(j=i+1;j<lsize_T;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j;k<t1;k+=2){
                x = calculate_sijk_TTE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TTE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }

    for(i=3;i<lsize_T;i+=2){
        s1 = lweight3D[i];
        for(j=i;j<lsize_T;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j+1;k<t1;k+=2){
                x = calculate_sijk_TTE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TTE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
        for(j=i+1;j<lsize_T;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j+1;k<t1;k+=2){
                x = calculate_sijk_TTE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TTE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }
    for(m=0;m<terms;m++) mvec[m] *= 6.0*deltaphi*deltaphi;

    return;
}

void calculate_gamma_TEE_3D(int n, double *mvec) {

    int i,j,k,m,t1,t2,t3;
    double x,y,z;
    double b1,b2,b3,b4;
    double q1,q2,q3,q4;

    int lsize_T = get_qtilde_T_lsize();
    int lsize_E = get_qtilde_E_lsize();
    int terms = get_terms_prim();
    double s1,s2,s3;

    for(m=0;m<terms;m++) mvec[m] = 0e0;

    for(i=2;i<lsize_T;i+=2){
        s1 = lweight3D[i];
        for(j=i;j<lsize_E;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j;k<t1;k+=2){
                x = calculate_sijk_TEE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TEE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
        for(j=i+1;j<lsize_E;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j;k<t1;k+=2){
                x = calculate_sijk_TEE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TEE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }

    for(i=3;i<lsize_T;i+=2){
        s1 = lweight3D[i];
        for(j=i;j<lsize_E;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j+1;k<t1;k+=2){
                x = calculate_sijk_TEE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TEE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
        for(j=i+1;j<lsize_E;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j+1;k<t1;k+=2){
                x = calculate_sijk_TEE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_TEE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }
    for(m=0;m<terms;m++) mvec[m] *= 6.0*deltaphi*deltaphi;

    return;
}

void calculate_gamma_EEE_3D(int n, double *mvec) {

    int i,j,k,m,t1,t2,t3;
    double x,y,z;
    double b1,b2,b3,b4;
    double q1,q2,q3,q4;

    int lsize_E = get_qtilde_E_lsize();
    int terms = get_terms_prim();
    double s1,s2,s3;

    for(m=0;m<terms;m++) mvec[m] = 0e0;

    for(i=2;i<lsize_E;i+=2){
        s1 = lweight3D[i];
        for(j=i;j<lsize_E;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j;k<t1;k+=2){
                x = calculate_sijk_EEE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_EEE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
        for(j=i+1;j<lsize_E;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j;k<t1;k+=2){
                x = calculate_sijk_EEE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_EEE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }

    for(i=3;i<lsize_E;i+=2){
        s1 = lweight3D[i];
        for(j=i;j<lsize_E;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j+1;k<t1;k+=2){
                x = calculate_sijk_EEE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_EEE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
        for(j=i+1;j<lsize_E;j+=2){
            s2 = lweight3D[j];
            t1 = i+j;
            if(t1>lsize_E)t1=lsize_E;
            for(k=j+1;k<t1;k+=2){
                x = calculate_sijk_EEE(n,i,j,k);
                s3 = lweight3D[k];
                z = permsix(i,j,k)*calculate_geometric(i,j,k)*s1*s2*s3;
                for(m=0;m<terms;m++){
                    y = plijk_EEE(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }
    for(m=0;m<terms;m++) mvec[m] *= 6.0*deltaphi*deltaphi;

    return;
}

_OFFLOADABLE
inline double Delta(double* xdiff, int i) {
	return xdiff[i];
}
_OFFLOADABLE
inline double delta(double* restrict xdiff, double* restrict ixdiff, double* restrict y, int i) {
	return (y[i+1] - y[i]) * ixdiff[i];
}
_OFFLOADABLE
inline double deriv(double* restrict xdiff, double* restrict ixdiff, double* restrict y, int i) {
        //return (Delta(xdiff,i-1)*delta(xdiff, ixdiff,y,i) + Delta(xdiff,i)*delta(xdiff,ixdiff,y,i-1)) / (Delta(xdiff,i-1) + Delta(xdiff,i));
	return (y[i+1] - y[i]) * ixdiff[i];
}

#define TRAPEZIUM
//#define HERMITE
_OFFLOADABLE double integrate(int size, double* restrict xdiff, double* restrict ixdiff, double* restrict y)
{

#ifdef HERMITE
	__assume_aligned(xdiff, 64);
	__assume_aligned(y, 64);
	__assume_aligned(ixdiff, 64);

	int i;
	double r = 0;
	#pragma loop count (213)
	for (i = 1; i < size-2; i++) {
		r += 0.5 * xdiff[i] * (y[i] + y[i+1] + xdiff[i] * (deriv(xdiff,ixdiff,y,i) - deriv(xdiff,ixdiff,y,i+1)) * sixth);
	}
	return r;
#endif

#ifdef TRAPEZIUM
	int i;
	double r = 0;
	for (i = 0; i < size-1; i++) {
		r += (xdiff[i]) * (y[i] + y[i+1]) * 0.5;
	}
	return r;
#endif

}

double calculate_sijk_TTT(int n, int xsize, double* restrict x, double* restrict y, int l1, int l2, int l3, double* restrict xdiff, double* restrict ixdiff)
{
    int i,a1,a2,a3;
    find_perm_prim(n,&a1,&a2,&a3);

    double p1,p2,p3,p4,p5,p6,p7,p8,p9;
    double t1,t2,t3,t4,t5,t6;

    for (i = 0; i < xsize; i++)
    {
        p1 = qtilde3D_T[l1][a1][i];
        p2 = qtilde3D_T[l1][a2][i];
        p3 = qtilde3D_T[l1][a3][i];
        p4 = qtilde3D_T[l2][a1][i];
        p5 = qtilde3D_T[l2][a2][i];
        p6 = qtilde3D_T[l2][a3][i];
        p7 = qtilde3D_T[l3][a1][i];
        p8 = qtilde3D_T[l3][a2][i];
        p9 = qtilde3D_T[l3][a3][i];

        t1 = p1*p5*p9;
        t2 = p2*p6*p7;
        t3 = p3*p4*p8;
        t4 = p3*p5*p7;
        t5 = p2*p4*p9;
        t6 = p1*p6*p8;

        y[i] = x[i] * x[i] * (t1+t2+t3+t4+t5+t6) * sixth;
    }

    // set spurious large values in first 10 entries to zero
    for (i = 0; i < xsize; i++)
        if( i<10 && fabs(y[i])>1e-18 ) y[i]=0;

    /*
    gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, xsize);
    gsl_interp_accel* acc = gsl_interp_accel_alloc();

    gsl_spline_init(sp,xvec,yvec,xsize);
    result = gsl_spline_eval_integ(sp,xmin,xmax,acc);

    gsl_spline_free(sp);
    gsl_interp_accel_free(acc);
    // 	free(yvec);
    */

    return integrate(xsize, xdiff, ixdiff, y);
}

double calculate_sijk_TTE(int n, int l1, int l2, int l3){

    /*
    double result;
    int i,a1,a2,a3;
    find_perm_prim(n,&a1,&a2,&a3);

    double p1,p2,p3,p4,p5,p6,p7,p8,p9;
    double t1,t2,t3,t4,t5,t6;

    int xsize = get_qtilde_xsize();
    double yvec[xsize];

    for (i=0;i<xsize;i++){

        p1 = qtilde3D_T[l1][a1][i];
        p2 = qtilde3D_T[l1][a2][i];
        p3 = qtilde3D_T[l1][a3][i];
        p4 = qtilde3D_T[l2][a1][i];
        p5 = qtilde3D_T[l2][a2][i];
        p6 = qtilde3D_T[l2][a3][i];
        p7 = qtilde3D_E[l3][a1][i];
        p8 = qtilde3D_E[l3][a2][i];
        p9 = qtilde3D_E[l3][a3][i];

        t1 = p1*p5*p9;
        t2 = p2*p6*p7;
        t3 = p3*p4*p8;
        t4 = p3*p5*p7;
        t5 = p2*p4*p9;
        t6 = p1*p6*p8;

        yvec[i] = x[i] * x[i] * (t1+t2+t3+t4+t5+t6)/6e0;
        if( i<10 && fabs(yvec[i])>1e-18 )yvec[i]=0;
    }

    gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, xsize);
    gsl_interp_accel* acc = gsl_interp_accel_alloc();

    gsl_spline_init(sp,xvec,yvec,xsize);
    result = gsl_spline_eval_integ(sp,xmin,xmax,acc);

    gsl_spline_free(sp);
    gsl_interp_accel_free(acc);
    // 	free(yvec);
    return result;
    */
    return -1;
}

double calculate_sijk_TEE(int n, int l1, int l2, int l3) {

    /*
    double result;
    int i,a1,a2,a3;
    find_perm_prim(n,&a1,&a2,&a3);

    double p1,p2,p3,p4,p5,p6,p7,p8,p9;
    double t1,t2,t3,t4,t5,t6;

    int xsize = get_qtilde_xsize();
    double yvec[xsize];

    for (i=0;i<xsize;i++){

        p1 = qtilde3D_T[l1][a1][i];
        p2 = qtilde3D_T[l1][a2][i];
        p3 = qtilde3D_T[l1][a3][i];
        p4 = qtilde3D_E[l2][a1][i];
        p5 = qtilde3D_E[l2][a2][i];
        p6 = qtilde3D_E[l2][a3][i];
        p7 = qtilde3D_E[l3][a1][i];
        p8 = qtilde3D_E[l3][a2][i];
        p9 = qtilde3D_E[l3][a3][i];

        t1 = p1*p5*p9;
        t2 = p2*p6*p7;
        t3 = p3*p4*p8;
        t4 = p3*p5*p7;
        t5 = p2*p4*p9;
        t6 = p1*p6*p8;

        yvec[i] = xvec[i] * xvec[i] * (t1+t2+t3+t4+t5+t6)/6e0;
        if( i<10 && fabs(yvec[i])>1e-18 )yvec[i]=0;
    }

    gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, xsize);
    gsl_interp_accel* acc = gsl_interp_accel_alloc();

    gsl_spline_init(sp,xvec,yvec,xsize);
    result = gsl_spline_eval_integ(sp,xmin,xmax,acc);

    gsl_spline_free(sp);
    gsl_interp_accel_free(acc);
    // 	free(yvec);
    return result;
    */
    return -1;
}

double calculate_sijk_EEE(int n, int l1, int l2, int l3) {

    /*
    double result;
    int i,a1,a2,a3;
    find_perm_prim(n,&a1,&a2,&a3);

    double p1,p2,p3,p4,p5,p6,p7,p8,p9;
    double t1,t2,t3,t4,t5,t6;

    int xsize = get_qtilde_xsize();
    double yvec[xsize];

    for (i=0;i<xsize;i++){

        p1 = qtilde3D_E[l1][a1][i];
        p2 = qtilde3D_E[l1][a2][i];
        p3 = qtilde3D_E[l1][a3][i];
        p4 = qtilde3D_E[l2][a1][i];
        p5 = qtilde3D_E[l2][a2][i];
        p6 = qtilde3D_E[l2][a3][i];
        p7 = qtilde3D_E[l3][a1][i];
        p8 = qtilde3D_E[l3][a2][i];
        p9 = qtilde3D_E[l3][a3][i];

        t1 = p1*p5*p9;
        t2 = p2*p6*p7;
        t3 = p3*p4*p8;
        t4 = p3*p5*p7;
        t5 = p2*p4*p9;
        t6 = p1*p6*p8;

        yvec[i] = xvec[i] * xvec[i] * (t1+t2+t3+t4+t5+t6)/6e0;
        if( i<10 && fabs(yvec[i])>1e-18 )yvec[i]=0;
    }

    gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, xsize);
    gsl_interp_accel* acc = gsl_interp_accel_alloc();

    gsl_spline_init(sp,xvec,yvec,xsize);
    result = gsl_spline_eval_integ(sp,xmin,xmax,acc);

    gsl_spline_free(sp);
    gsl_interp_accel_free(acc);
    // 	free(yvec);
    return result;
    */
    return -1;
}
