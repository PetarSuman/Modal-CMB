#ifndef H_TETRAPYD
#define H_TETRAPYD

#include "offload_util.h"

/**
 * tetrapyd_tools.h
 * 
 * Functions for decomposing and looping over l1,l2,l3 space subject to
 * triangle and parity rules.
 *
 * Example usage:
 * int nthreads = omp_get_max_threads();
 * tetrapyd_limits* block;
 * block = (tetrapyd_limits*) malloc(nthreads * sizeof(tetrapyd_limits));
 * decompose_tetrapyd_XXX(block, rank, nproc, nthreads, 0, 1D_size);
 *
 * // iterare through region is:
 * #pragma omp
 * 
 * t = omp_get_thread_num();
 * i = block[t].i_bgn;
 * j = block[t].j_bgn;
 * k = block[t].k_bgn;
 * for(n=0;n<block[t].loops;n++){
 * ... * calculation * ...
 * 	get_ijk_next_XXX(max,&i,&j,&k);
 *
 */

typedef struct
{
    int i_bgn;
    int j_bgn;
    int k_bgn;
    int loops;
} tetrapyd_limits;

/** 
 * Decompose region defined by triangle and parity rules into chunks for 3D integration
 * usage for initilisation. Finds the start point and loop count for this MPI process and
 * subdivides that further, returning the start points and loop counts for each of this rank's
 * OpenMP threads.
 *
 * @param[out] thread_limits -- Array of tetrapyd_limits nthreads long, each entry contains the
 *                          start ijk and loop count for each thread/
 * @param rank -- MPI rank of this process.
 * @param numranks -- Number of MPI ranks in this MPI world.
 * @param nthreads -- Number of OpenMP threads in this rank.
 * @param lmin -- Minimum l value.
 * @param lmax -- Maximum l value.
 */
_OFFLOADABLE void decompose_tetrapyd_XXX(tetrapyd_limits* thread_limits, int rank, int numranks, int nthreads, int lmin, int lmax);

/**
 * Given values for i,j,k return their next values in the sequence subject to the tetrapyd rules
 * and return whether or not it was successful.
 *
 * @param lmax -- Maximum l value
 * @param[inout] i/j/k -- Current i/j/k goes in, next valid i/j/k comes out.
 * @return 0 if next ijk is out-of-bounds, else 1.
 */
_OFFLOADABLE int get_ijk_next_XXX(int lmax, int *i, int *j, int *k);
_OFFLOADABLE void check_ijk_max_XXX(int lmax, int *i, int *j, int *k);
_OFFLOADABLE void check_ijk_min_XXX(int lmin, int *i, int *j, int *k);

/**
 * Decompose region defined by triangle rule (with 2 steps over boundary for correct interpolation) into chunks for 3D integration.
 * Finds the start point and loop count for this MPI process and
 * subdivides that further, returning the start points and loop counts for each of this rank's
 * OpenMP threads.
 * @param thread_limits[out] -- Array of tetrapyd_limits nthreads long, each entry contains the
 *                          start ijk and loop count for each thread/
 * @param rank -- MPI rank of this process.
 * @param numranks -- Number of MPI ranks in this MPI world.
 * @param nthreads -- Number of OpenMP threads in this rank.
 * @param min -- Minimum l value.
 * @param max -- Maximum l value.
**/

_OFFLOADABLE void decompose_tetrapyd_XXY(tetrapyd_limits* thread_limits, int rank, int numranks, int nthreads, int lminX, int lminY, int lmaxX, int lmaxY);
_OFFLOADABLE int get_ijk_next_XXY(int lminX, int lminY, int lmaxX, int lmaxY, int *i, int *j, int *k);
_OFFLOADABLE void check_ijk_max_XXY(int lmaxX, int lmaxY, int *i, int *j, int *k);
_OFFLOADABLE void check_ijk_min_XXY(int lminX, int lminY, int *i, int *j, int *k);


void decompose_tetrapyd_prim(tetrapyd_limits* thread_limits, int rank, int numranks, int nthreads, int min, int max);

/**
 * Given values for i,j,k return their next values in the sequence subject to the tetrapyd rules
 * and with 2 steps over the boundary for correct interpolation.
 * return whether or not it was successful.
 *
 * @param lmax -- Maximum l value
 * @param[inout] i/j/k -- Current i/j/k goes in, next valid i/j/k comes out.
 * @return 0 if next ijk is out-of-bounds, else 1.
 */
int get_ijk_next_prim(int max, int *i, int *j, int *k);
void check_ijk_max_prim(int max, int *i, int *j, int *k);
void check_ijk_min_prim(int min, int *i, int *j, int *k);

/**
 * Given an amount of work 'tasks' and some parallel workers 'nranks', 
 * evenly divide the work between the workers and return an array of 
 * the amount of work for each worker.
 *
 * @param tasks -- The amount of work or iterations to be done.
 * @param nranks -- The number of parallel workers (MPI or threads).
 * @param[out] tasks_per_rank -- The amount of work for each worker.
 */
_OFFLOADABLE void divide_tasks(long int tasks, int nranks, long int* tasks_per_rank);

#endif
