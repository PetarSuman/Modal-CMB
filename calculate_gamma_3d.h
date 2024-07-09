#ifndef H_GAMMA_3D
#define H_GAMMA_3D

#define I_CUTOFF 3

#include "local_barrier.h"

#define MAXBLOCKSIZE 64
/**
 * thread private data in the multi-threaded gamma_3d algorithm 
 */
typedef struct
{
    int i0;
    int j0;
    int k0;
    int i1;
    int j1;
    int k1;
    int xsize;
    double* mvec;
    double* xvec;
    double* yvec;
    double* xdiff;
    double* ixdiff;
    double* intgrlvec;  // array of x in calculate_gamma_3d 
    double* plijkz;     // array of plijk*z in calculate_gamma_3d
    barrier_t* bar;     // ptr to struct for core barrier data
} thread_g3d_data_t;

/**
 * Precompute look-up arrays for calculate_gamma_3D_* functions.
 */
void init_gamma_3Dint();

/**
 * Calculates the domain decomposition of the l1l2l3 loop for the TTT case by assigning every 
 * thread on every MPI rank unique upper and lower loop coordinates in its entry to the
 * array `thread_data`.
 *
 * @param thread_data -- Pointer to an array of type thread_g3d_data_t containing the data 
 *                          for each thread. The i0,j0,k0,i1,j1,k1 entries are updated in
 *                          this function. 
 * @param rank -- MPI rank of the caller..
 * @param numranks -- The number of MPI ranks to decompose the iteration space.
 * @param workshare -- Ratio of work done in offload to work done on the host. Allowed values
 *                      are in the range [0,1].
 */
void decompose_gamma_TTT_3D_mpi(thread_g3d_data_t* thread_data, int rank, int numranks, double workshare); 

/**
 * Function that moves through the l1l2l3 iteration space in blocks. Given a starting
 * coordinate in iteration space (i0,j0,k0) and a desire to move maxblocksize number of iterations
 * forwards and a boundary coordinate of (final_i, final_j, final_k), return the next coordinates
 * (i1,j1,k1) in iteration space and the blocksize. Where blocksize is the actual number of
 * iterations moved forwards. In the case where the iterations go beyond the boundary coordinate 
 * (i1,j1,k1) is set to the boundary coordinate blocksize is set the the trunkcated number of 
 * iterations (<=maxblocksize).
 *
 * @param i0 -- lower i coordinate
 * @param j0 -- lower j coordinate
 * @param k0 -- lower k coordinate
 * @param[out] i1 -- Upper i coordinate
 * @param[out] j1 -- Upper j coordinate
 * @param[out] k1 -- Upper k coordinate
 * @param final_i -- Boundary i coordinate
 * @param final_j -- Boundary j coordinate
 * @param final_k -- Boundary k coordinate
 * @param blocksize -- Number of iterations between (i0,j0,k0) and (i1,j1,k1)
 * @param[out] maxblocksize -- Max number of iterations to try to move foward from (i0,j0,k0)
 */
void get_next_ijk_bound(int i0, int j0, int k0, int *i1, int *j1, int *k1, \
        int final_i, int final_j, int final_k, int* blocksize, int maxblocksize);


/**
 * Calculate the Gamma matrix for the TTT case using the optimised multi-threaded 3D method.
 *
 * @param thread_data -- Pointer to an array of type thread_g3d_data_t containing the data 
 *                          for each thread needed to execute this function in parallel.
 * @param maxblocksize -- Optimisation parameter used for tuning cache blocking. Typical value=64.
 */
void calculate_gamma_TTT_3D_nested(thread_g3d_data_t* thread_data, int maxblocksize);

/**
 * Calculate the Gamma matrix for the TTT case using the old method. DEPRECATED.
 */
void calculate_gamma_TTT_3D(int n, double *mvec);

/**
 * Calculate the Gamma matrix for the TTE case using the old method. DEPRECATED.
 */
void calculate_gamma_TTE_3D(int n, double *mvec);

/**
 * Calculate the Gamma matrix for the TEE case using the old method. DEPRECATED.
 */
void calculate_gamma_TEE_3D(int n, double *mvec);

/**
 * Calculate the Gamma matrix for the EEE case using the old method. DEPRECATED.
 */
void calculate_gamma_EEE_3D(int n, double *mvec);

#endif
