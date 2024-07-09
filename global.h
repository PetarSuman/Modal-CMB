#ifndef H_GLOBAL
#define H_GLOBAL

// Declare user defined structures

#define MAXLEN 200   // max length of char arrays
#define MAXLINES 10000000 // max lines in input files

#include <gsl/gsl_spline.h>
#include "offload_util.h"
#include "tetrapyd_tools.h"

// useful helper macros
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

typedef struct {
	int i;
	int j;
	int orientation;
	int flag;
	double area;
} cell;

typedef struct {
	double x;
	double y1;
	double y2;
	double y3;
} one_points;

typedef struct {
	double a1;
	double a2;
	double a3;
} shape_params;

typedef struct {
	int p1_bgn;
	int p1_stp;
	int p1_end;
	int p2_bgn;
	int p2_stp;
	int p2_end;
	int p3_bgn;
	int p3_stp;
	int p3_end;
} shape_limits;

/*
typedef struct
{
    int i_bgn;
    int j_bgn;
    int k_bgn;
    int loops;

} tetrapyd_limits;

void decompose_tetrapyd_XXX(tetrapyd_limits* thread_limits, int rank, int numranks, int lmin, int lmax);
int get_ijk_next_XXX(int lmax, int *i, int *j, int *k);
void check_ijk_max_XXX(int lmax, int *i, int *j, int *k);
void check_ijk_min_XXX(int lmin, int *i, int *j, int *k);


void decompose_tetrapyd_prim(tetrapyd_limits* thread_limits, int rank, int numranks, int min, int max);
int get_ijk_next_prim(int max, int *i, int *j, int *k);
void check_ijk_max_prim(int max, int *i, int *j, int *k);
void check_ijk_min_prim(int min, int *i, int *j, int *k);

void divide_tasks(long int tasks, int ranks, int* ranksize);
*/

// Declare external variables
extern char alpha_data_file[MAXLEN];
extern char alpha_size_file[MAXLEN];
extern char bessel_data_file[MAXLEN];
extern char bessel_size_file[MAXLEN];
extern char transfer_T_data_file[MAXLEN];
extern char transfer_E_data_file[MAXLEN];
extern char transfer_size_file[MAXLEN];
extern char source_data_file[MAXLEN];
extern char source_size_file[MAXLEN];
extern char limber_data_file[MAXLEN];
extern char limber_size_file[MAXLEN];
extern char decompose_data_file[MAXLEN];
extern char decompose_size_file[MAXLEN];
extern char eigen_data_file[MAXLEN];
extern char eigen_tri_data_file[MAXLEN];
extern char ortho_data_file[MAXLEN];
extern char orthol_TTT_data_file[MAXLEN];
extern char orthol_TTE_data_file[MAXLEN];
extern char orthol_TEE_data_file[MAXLEN];
extern char orthol_EEE_data_file[MAXLEN];
extern char ortho_tri_data_file[MAXLEN];
extern char orthol_tri_data_file[MAXLEN];
extern char lambda_data_file[MAXLEN];
extern char lambdal_TTT_data_file[MAXLEN];
extern char lambdal_TTE_data_file[MAXLEN];
extern char lambdal_TEE_data_file[MAXLEN];
extern char lambdal_EEE_data_file[MAXLEN];
extern char lambda_tri_data_file[MAXLEN];
extern char lambdal_tri_data_file[MAXLEN];
extern char gamma_TTT_data_file[MAXLEN];
extern char gamma_TTE_data_file[MAXLEN];
extern char gamma_TEE_data_file[MAXLEN];
extern char gamma_EEE_data_file[MAXLEN];
extern char gamma_tri_data_file[MAXLEN];
extern char proj_T_data_file[MAXLEN];
extern char proj_E_data_file[MAXLEN];
extern char proj_T_size_file[MAXLEN];
extern char proj_E_size_file[MAXLEN];
extern char proj_tri_data_file[MAXLEN];
extern char proj_tri_size_file[MAXLEN];
extern char modes_TTT_data_file[MAXLEN];
extern char modes_TTE_data_file[MAXLEN];
extern char modes_TEE_data_file[MAXLEN];
extern char modes_EEE_data_file[MAXLEN];
extern char modes_tri_data_file[MAXLEN];
extern char flat_data_file[MAXLEN];
extern char flat_size_file[MAXLEN];
extern char l_data_file[MAXLEN];
extern char l_size_file[MAXLEN];
extern char bispectrum_file[MAXLEN];
extern char trispectrum_file[MAXLEN];
extern char interpolated_l_data_file[MAXLEN];
extern char interpolated_l_size_file[MAXLEN];
extern char interpolated_bispectrum_file[MAXLEN];
extern char pixw_file[MAXLEN];
extern char BN_TT_file[MAXLEN];
extern char BN_EE_file[MAXLEN];
extern char transfer_wgt_file[MAXLEN];
extern char cls_scalar_data_file[MAXLEN];
extern char cls_lensed_data_file[MAXLEN];
extern char restart_file[MAXLEN];
extern char bload_file[MAXLEN];
extern char edist_file[MAXLEN];
extern char edist_tri_file[MAXLEN];
extern char mdist_TTT_file[MAXLEN];
extern char mdist_TTE_file[MAXLEN];
extern char mdist_TEE_file[MAXLEN];
extern char mdist_EEE_file[MAXLEN];
extern char mdist_tri_file[MAXLEN];

extern double nscalar;
extern _OFFLOADABLE double deltaphi;
extern double kpivot;
extern int l_flat;
extern int bessel_max;
extern int bessel_points;
extern int multiplier;
extern int section_count;
extern double one_accuracy;
extern double tilt;
extern int initial_grid_size;
extern int depth;
extern double tri_accuracy;
extern int alpha_max;
extern int alpha_points;
extern int step;
extern int l1,l2,l3,l4;
extern int model;
extern int model2;
extern shape_limits slim;

// Flags
extern _OFFLOADABLE int do_polarisation;
extern int use_l_file;
extern char lini_file[MAXLEN];
extern int eflag_order_prim;
extern int eflag_order_late;
extern _OFFLOADABLE int eflag_T_lmax;
extern _OFFLOADABLE int eflag_T_lmin;
extern _OFFLOADABLE int eflag_E_lmax;
extern _OFFLOADABLE int eflag_E_lmin;
extern double fskyT;
extern double fskyE;
extern int rflag_do3D;
extern int clflag_uselens;
extern int tflag_bessel;
extern int tflag_transfer;
extern int tflag_flat;
extern int tflag_limber;
extern int bflag_lset;
extern int bflag_load;
extern int dflag_out;
extern int dflag_dmax;
extern int gflag_pca;
extern int oflag_load;

// named

extern int initilise(char *inifile);

extern double spherical_bessel(int l , double z);

extern double calculate_flat(int l1, int l2, int l3);
extern double calculate_full(int l1, int l2, int l3);
extern double calculate_limber1(int l1, int l2, int l3);
extern double calculate_limber2(int l1, int l2, int l3);
extern double calculate_limber3(int l1, int l2, int l3);
extern void calculate_ortho(int min, int max, tetrapyd_limits* block, double* ortho);
// extern double calculate_ortho(int r,int s);
extern double calculate_ortho_mc(int r,int s);
extern double calculate_ortho_tri(int r,int s);
extern double calculate_check_tri(int r,int s);
// extern double calculate_orthol(int r, int s, int size, int* vec);
extern double init_orthol_glint();
extern void calculate_orthoMij_T(int i, int j);
extern void calculate_orthoMij_E(int i, int j);
extern int calculate_orthol(int r, int s, double* result);
extern int calculate_orthol_3D(int r, int s, double* result);
extern double calculate_orthol_tri_old(int r,int s, int size, int* vec);
extern double calculate_delta(int r,int s, int size);
extern double calculate_deltal(int r, int s, int size, int* vec);
extern double integrate_template2(int min, int max, shape_params params, tetrapyd_limits* block, double* full_res); 
extern void calculate_eigen(int min, int max, shape_params params, tetrapyd_limits* block, double* alpha);
extern double calculate_eigen_mc(int r);
extern double calculate_eigen_tri(int r);
extern double calculate_eigen_tri_mc(int r);
extern double check_eigen_tri(int num, double max);
extern void modify_bispectrum(int size, int* vec, double*** bispectrum);
extern double calculate_modes(int r, int size, double*** bispectrum);
extern void calculate_modes_file(char *filename, int models, double **modearray);
extern double calculate_modes_tri(int r, int size, int* vec, double**** trispectrum);
extern double init_ISW_glint();
extern double calculate_ISW(int i);
extern double calculate_ISW_norm();
extern void calculate_ISW_3D(int min, int max, tetrapyd_limits* block, double* alpha);
extern double calculate_ISW_norm_3D();
extern double calculate_ISW_bispectrum(int i,int j,int k);
extern double init_PS_glint();
extern double calculate_PS(int i);
extern double calculate_PS_norm();
extern void calculate_PS_3D(int min, int max, tetrapyd_limits* block, double* alpha);
extern double calculate_PS_norm_3D();
extern void init_3D(int mode);
extern void calculate_mode_3D(int min, int max, tetrapyd_limits* block, double* alpha);
extern double calculate_3D_bispectrum(int l1, int l2, int l3);
extern double*** reconstruct_b(int size, int* vec);
extern double*** reconstruct_bn(int max, int size, int* vec);
extern int get_size();
extern double check(shape_params params);
extern double correlation_prim(shape_params params1, shape_params params2);
extern double check_tri(double kmax);
extern double checkmode(int num, shape_params params);
extern double check_eigen_tri(int num, double max);
extern double check_norm_tri(double kmax);
extern void find_perm(int r,int* p1,int* p2,int* p3);
extern double pijk(int r,int i,int j,int k);
extern double plijk_TTT(int r,int i,int j,int k);
extern double plijk_TTE(int r,int i,int j,int k);
extern double plijk_TEE(int r,int i,int j,int k);
extern double plijk_EEE(int r,int i,int j,int k);
extern double rijk(int r,int i,int j,int k);
extern double rlijk_TTT(int r,int i,int j,int k);
extern double rlijk_TTE(int r,int i,int j,int k);
extern double rlijk_TEE(int r,int i,int j,int k);
extern double rlijk_EEE(int r,int i,int j,int k);
extern double sijk(int i,int j,int k);
extern double bijk_TTT(int i,int j,int k);
extern double bijk_TTE(int i,int j,int k);
extern double bijk_TEE(int i,int j,int k);
extern double bijk_EEE(int i,int j,int k);

extern double triQ(int r,int i,int j,int k,int l);
extern double triR(int r,int i,int j,int k,int l);
extern double trilR(int r,int i,int j,int k,int l);
extern double triS(int i,int j,int k,int l);

extern double shape(double k1, double k2, double k3, shape_params params);
extern double shape2(double k1, double k2, double k3, shape_params params);
extern double shape3(double k1, double k2, double k3, shape_params params);
extern double scale(double k, shape_params params);
extern double slice(double k1, double k2, double k3, shape_params params);
extern double smooth1_in(double k1, double k2, double k3, shape_params params);
extern double smooth1(double k1, double k2, double k3, shape_params params);
extern double smooth2_in(double k1, double k2, double k3, shape_params params);
extern double smooth2(double k1, double k2, double k3, shape_params params);
extern double smooth3_in(double k1, double k2, double k3, shape_params params);
extern double smooth3(double k1, double k2, double k3, shape_params params);

extern double cosmo_collider(double k1, double k2, double k3, shape_params params);
extern double equicollider(double k1, double k2, double k3, shape_params params);
extern double multispeed(double k1, double k2, double k3, shape_params params);
extern double quasilocalfield(double k1, double k2, double k3, shape_params params);

extern double dbi(double k1, double k2, double k3, shape_params params);
extern double equilateral(double k1, double k2, double k3, shape_params params);
extern double NBD(double k1, double k2, double k3, shape_params params);
extern double NBD_cut(double k1, double k2, double k3, shape_params params);
extern double ghost(double k1, double k2, double k3, shape_params params);
extern double local(double k1, double k2, double k3, shape_params params);
extern double mode1(double k1, double k2, double k3, shape_params params);
extern double single(double k1, double k2, double k3, shape_params params);
extern double nicola(double k1, double k2, double k3, shape_params params);
extern double warm(double k1, double k2, double k3, shape_params params);
extern double warm_cut(double k1, double k2, double k3, shape_params params);

extern double nonlocal(double k1, double k2, double k3, shape_params params);
extern double feature(double k1, double k2, double k3, shape_params params);
extern double feature_scale(double k1, double k2, double k3, shape_params params);
extern double feature_slice(double k1, double k2, double k3, shape_params params);
extern double sinlog(double k1, double k2, double k3, shape_params params);
extern double sinlogequi(double k1, double k2, double k3, shape_params params);
extern double sinlogflat(double k1, double k2, double k3, shape_params params);
extern double sinlogtest(double k1, double k2, double k3, shape_params params);
extern double sinloglocal(double k1, double k2, double k3, shape_params params);
extern double sinlog_scale(double k1, double k2, double k3, shape_params params);
extern double sinlog_slice(double k1, double k2, double k3, shape_params params);
extern double maldacena(double k1, double k2, double k3, shape_params params);
extern double string(double k1, double k2, double k3, shape_params params);
extern double constant(double k1, double k2, double k3, shape_params params);
extern double chen(double k1, double k2, double k3, shape_params params);
extern double QSF(double k1, double k2, double k3, shape_params params);
extern double adshead(double k1, double k2, double k3, shape_params params);
extern double orthogonal(double k1, double k2, double k3, shape_params params);
extern double flat(double k1, double k2, double k3, shape_params params);
extern double neil1(double k1, double k2, double k3, shape_params params);
extern double neil2(double k1, double k2, double k3, shape_params params);
extern double NBD_ivan_cos1(double k1, double k2, double k3, shape_params params);
extern double NBD_ivan_cos2(double k1, double k2, double k3, shape_params params);
extern double NBD_ivan_sin1(double k1, double k2, double k3, shape_params params);
extern double NBD_ivan_sin2(double k1, double k2, double k3, shape_params params);
extern double NBD_xingang(double k1, double k2, double k3, shape_params params);
extern double NBD_DBI(double k1, double k2, double k3, shape_params params);
extern double NBD_sin(double k1, double k2, double k3, shape_params params);
extern double sinfit(double k1, double k2, double k3, shape_params params);
extern double sinflat(double k1, double k2, double k3, shape_params params);
extern double NBD_sinlog(double k1, double k2, double k3, shape_params params);
extern double EFT1(double k1, double k2, double k3, shape_params params);
extern double EFT2(double k1, double k2, double k3, shape_params params);
extern double featureENV(double k1, double k2, double k3, shape_params params);
extern double featureENV2(double k1, double k2, double k3, shape_params params);
extern double featureENV3(double k1, double k2, double k3, shape_params params);
extern double featureENV4(double k1, double k2, double k3, shape_params params);
extern double featureENV5(double k1, double k2, double k3, shape_params params);
extern double featureENV6(double k1, double k2, double k3, shape_params params);
extern double featureENV7(double k1, double k2, double k3, shape_params params);
extern double nicola529(double k1, double k2, double k3, shape_params params);
extern double nicolanum(double k1, double k2, double k3, shape_params params);
extern double ribero(double k1, double k2, double k3, shape_params params);
extern double ksin(double k1, double k2, double k3, shape_params params);
extern double expsin(double k1, double k2, double k3, shape_params params);
extern double expcos(double k1, double k2, double k3, shape_params params);
extern double equiexpsin(double k1, double k2, double k3, shape_params params);
extern double equiexpcos(double k1, double k2, double k3, shape_params params);
extern double flatexpsin(double k1, double k2, double k3, shape_params params);
extern double flatexpcos(double k1, double k2, double k3, shape_params params);
extern double cs_feature(double k1, double k2, double k3, shape_params params);
extern double freqdrift(double k1, double k2, double k3, shape_params params);
extern void loadnicola();

extern double shape_tri(double k1, double k2, double k3, double k4);
extern double shape_tri3(double k1, double k2, double k3, double k4);
extern double const_tri(double k1, double k2, double k3, double k4);
extern double equi_tri(double k1, double k2, double k3, double k4);
extern double localg_tri(double k1, double k2, double k3, double k4);
extern double localt_tri(double k1, double k2, double k3, double k4);

extern int basis_functions_legendre(double **array, int size, int pmax, double min, double max, double *values);
extern int basis_functions_chebyshev(double **array, int size, int pmax, double min, double max, double *values);
extern int basis_functions_fourier(double **array, int size, int pmax, double min, double max, double *values);
extern int basis_functions_sinlog(double **array, int size, int pmax, double min, double max, double *values);
extern double sinlog_pt(double kpt, int j);
extern int basis_functions_bi(double **array, int size, int pmax, double min, double max, double *values);
extern int basis_functions_tri(double **array, int size, int pmax, double min, double max, double *values);
extern double calculate_alpha(int rank, double *kvalues, int i, int j, int k);
extern int cubic_interpolation(double ***in, double ***out, int *cell);
extern double calculate_geometric(int l1, int l2, int l3);
extern double calculate_geometric_tri(int l1, int l2, int l3, int l4);
extern double calculate_kint_T(int alpha, int l, double p);
extern double calculate_kint_E(int alpha, int l, double p);
extern double calculate_xint_TTT(int l1, int l2, int l3, int n, int xsize, double *xvec);
extern double calculate_xint_TTE(int l1, int l2, int l3, int n, int xsize, double *xvec);
extern double calculate_xint_TEE(int l1, int l2, int l3, int n, int xsize, double *xvec);
extern double calculate_xint_EEE(int l1, int l2, int l3, int n, int xsize, double *xvec);
extern double calculate_kint_tri(int l, int alpha, double p);
extern double calculate_xint_tri(int l1, int l2, int l3, int l4, int n, int xsize, double *xvec);

extern double calculate_volume(int i, int j, int k, double ***points);
extern double calculate_weight(int min, int max, int i, int j, int k);
extern double calculate_volume_tri(int i, int j, int k, int l, double ****points);

extern double calculate_fisher_tri(int i, int j);
extern double calculate_gamma_tri(int m, int n);
extern double calculate_orthol_tri(int m, int n);

extern void init_gamma_glint();
extern void init_gamma_3Dint();
extern _OFFLOADABLE void precompute_gammaMij_T(); // jb 
extern _OFFLOADABLE void precompute_gammaMij_E(); // jb
extern _OFFLOADABLE void calculate_gammaMij_T(int i, int j);
extern _OFFLOADABLE void calculate_gammaMij_E(int i, int j);
extern void sync_gammaMij();
extern void calculate_gamma(int m, int n, double* result);
extern _OFFLOADABLE void calculate_gamma_2d(int m, int n, double *result); //jb
extern _OFFLOADABLE double integrate_cubic_spline(double *x, double *y, int npts); //jb
extern void calculate_gamma_TTT_3D(int m, double *mvec);
extern void calculate_gamma_TTE_3D(int m, double *mvec);
extern void calculate_gamma_TEE_3D(int m, double *mvec);
extern void calculate_gamma_EEE_3D(int m, double *mvec);
extern _OFFLOADABLE double gamma_2d_pt_TTT(int m, int n, int i); // jb x4
extern _OFFLOADABLE double gamma_2d_pt_TTE(int m, int n, int i);
extern _OFFLOADABLE double gamma_2d_pt_TEE(int m, int n, int i);
extern _OFFLOADABLE double gamma_2d_pt_EEE(int m, int n, int i);

extern double init_gamma_tri_glint();
extern double calculate_gamma_tri(int m, int n);
extern void calculate_gamma_tri_4D(int m, int i, double *mvec);

extern double init_norm_glint();
extern double calculate_norm(int m, int n);

extern double calculate_correlation(double kmin, double kmax, double nu1, double nu2);
extern double calculate_localnorm(double kmin, double kmax, double nu1, double nu2);
extern void init_rng(int seed);
extern double get_rng_gauss();

// gl integration
extern int asy(double *nodes, double*weights, unsigned long int n);

// fixed data

extern int init_lmax(int l);
extern int get_lmax();

extern int set_terms_prim();
extern int get_terms_prim();
extern _OFFLOADABLE int get_pmax_prim();
extern int create_order_prim();
extern int create_order_tri_prim();
extern _OFFLOADABLE void find_perm_prim(int r,int* p1,int* p2,int* p3);

extern int set_terms_late();
extern _OFFLOADABLE int get_terms_late();
//extern _OFFLOADABLE int get_pmax_late();
extern _OFFLOADABLE int get_pmax_late_T();
extern _OFFLOADABLE int get_pmax_late_E();
extern int create_order_late();
extern int create_order_tri_late();
extern void find_perm_late(int r,int* p1,int* p2,int* p3);
extern _OFFLOADABLE void find_perm_late_TTT(int r,int* p1,int* p2,int* p3);
extern _OFFLOADABLE void find_perm_late_TTE(int r,int* p1,int* p2,int* p3);
extern _OFFLOADABLE void find_perm_late_TEE(int r,int* p1,int* p2,int* p3);
extern _OFFLOADABLE void find_perm_late_EEE(int r,int* p1,int* p2,int* p3);

extern int set_terms_tri();
extern int get_terms_tri();
extern int get_qmax();
extern int create_order_tri();
extern void find_perm_tri(int r,int* p1,int* p2,int* p3,int* p4);
extern void find_perm_tri_prim(int r,int* p1,int* p2,int* p3,int* p4);

extern int load_cl();
extern double get_cl_TT(int i);
extern double get_cl_TE(int i);
extern double get_cl_EE(int i);
extern double get_cl_TP(int i);

extern int load_BN();
extern int load_Tl(int l_size ,int* l_values);
extern int load_lens();
extern double get_beam_TT(int i);
extern double get_beam_TE(int i);
extern double get_beam_EE(int i);
extern double get_noise_TT(int i);
extern double get_noise_TE(int i);
extern double get_noise_EE(int i);
extern int create_t_wgt(int l_size);
extern double get_t_wgt(int i);

extern int create_alpha();
extern double get_alpha(int i, int j, int k);
extern int update_alpha(double *results);
extern int output_alpha();
extern int read_alpha();

extern int create_basis_prim(int ksize, double kmax, double *k);
extern double get_basis_prim(int i, int j);
extern int destroy_basis_prim();

extern int create_basis_late(int Tsize, int Esize, double Tmin, double Tmax, double Emin, double Emax, double *Tvec, double *Evec);
extern double get_basis_late_T(int i, int j);
extern double get_basis_late_E(int i, int j);
extern double* get_basis_late_T_flat();
extern int destroy_basis_late();

extern int create_basis_tri_prim(int ksize, double kmax, double *k);
extern double get_basis_tri_prim(int i, int j);
extern int destroy_basis_tri_prim();

extern int create_basis_tri_late(int ksize, double kmin, double kmax, double *k);
extern double get_basis_tri_late(int i, int j);
extern int destroy_basis_tri_late();

extern int create_qtilde(int lsize_T, int lsize_E, int x_size);
extern double get_qtilde_T(int i, int j, int k);
extern int update_qtilde_T(int i, int j, int k, double value);
extern double get_qtilde_E(int i, int j, int k);
extern int update_qtilde_E(int i, int j, int k, double value);
extern int output_qtilde(int xsize, double *xvec, int lsize_T, int *lvec_T, int lsize_E, int *lvec_E);
extern int read_qtilde();
extern int destroy_qtilde();

extern int get_qtilde_T_lsize();
extern int get_qtilde_E_lsize();
extern _OFFLOADABLE int get_qtilde_xsize();
extern int get_qtilde_T_lvec(int *vec);
extern int get_qtilde_E_lvec(int *vec);
extern _OFFLOADABLE int get_qtilde_xvec(double *vec);

extern int create_beta_tri(int l_size, int x_size);
extern double get_beta_tri(int i, int j, int k);
extern int update_beta_tri(int size, double *results);
extern int send_beta_tri(int nproc, int rank, int l_size, int x_size);
extern int output_beta_tri(int xsize, double *xvec, int lsize, int *lvec);
extern int read_beta_tri();
extern int destroy_beta_tri();

extern int get_bt_lsize();
extern int get_bt_xsize();
extern int get_bt_lvec(int *vec);
extern int get_bt_xvec(double *vec);

extern int create_bispectrum(int size);
extern double get_bispectrum(int i, int j, int k);
extern int update_bispectrum(double *results);
extern int output_bispectrum(int size, int *lvalues);

extern int create_trispectrum(int size);
extern double get_trispectrum(int i, int j, int k, int l);
extern int update_trispectrum(double *results);
extern int output_trispectrum(int size, int *lvalues);

extern int create_decompose();
extern double get_decompose(int i, int j, int k);
extern int update_decompose(double *results);
extern int output_decompose();
extern int read_decompose();

extern int create_eigen();
extern double get_eigen(int i);
extern int update_eigen(double *results);
extern int output_eigen(int p1, int p2, int p3);
extern int read_eigen(int p1, int p2, int p3);
extern int destroy_eigen_tri();

extern int create_eigen_tri();
extern double get_eigen_tri(int i);
extern int update_eigen_tri(double *results);
extern int output_eigen_tri();
extern int read_eigen_tri();

extern int create_eigenR();
extern double get_eigenR(int i);
extern int update_eigenR(double *results);

extern int create_eigenR_tri();
extern double get_eigenR_tri(int i);
extern int update_eigenR_tri(double *results);

extern int create_ortho();
extern double get_ortho(int i,int j);
extern int update_ortho(double *results);
extern int output_ortho();
extern int read_ortho();

extern int create_orthol();
extern double get_orthol_TTT(int i,int j);
extern double get_orthol_TTE(int i,int j);
extern double get_orthol_TEE(int i,int j);
extern double get_orthol_EEE(int i,int j);
extern int update_orthol(double *results);
extern int output_orthol();
extern int read_orthol();

extern int create_ortho_tri();
extern double get_ortho_tri(int i,int j);
extern int update_ortho_tri(double *results);
extern int output_ortho_tri();
extern int read_ortho_tri();

extern int create_orthol_tri();
extern double get_orthol_tri(int i,int j);
extern int update_orthol_tri(double *results);
extern int output_orthol_tri();
extern int read_orthol_tri();

extern int create_lambda();
extern double get_lambda(int i,int j);
extern int update_lambda(double *results);
extern int output_lambda();
extern int read_lambda();

extern int create_lambdal();
extern double get_lambdal_TTT(int i,int j);
extern double get_lambdal_TTE(int i,int j);
extern double get_lambdal_TEE(int i,int j);
extern double get_lambdal_EEE(int i,int j);
extern int update_lambdal(double *results);
extern int output_lambdal();
extern int read_lambdal();

extern int create_lambda_tri();
extern double get_lambda_tri(int i,int j);
extern int update_lambda_tri(double *results);
extern int output_lambda_tri();
extern int read_lambda_tri();

extern int create_lambdal_tri();
extern double get_lambdal_tri(int i,int j);
extern int update_lambdal_tri(double *results);
extern int output_lambdal_tri();
extern int read_lambdal_tri();

extern int create_gamma();
extern double get_gamma_TTT(int i,int j);
extern double get_gamma_TTE(int i,int j);
extern double get_gamma_TEE(int i,int j);
extern double get_gamma_EEE(int i,int j);
extern int update_gamma(double *results);
extern int output_gamma();
extern int read_gamma();

extern int create_gamma_tri();
extern double get_gamma_tri(int i,int j);
extern int update_gamma_tri(double *results);
extern int output_gamma_tri();
extern int read_gamma_tri();

extern int create_modes();
extern double get_modes_TTT(int i);
extern double get_modes_TTE(int i);
extern double get_modes_TEE(int i);
extern double get_modes_EEE(int i);
extern int update_modes(double *results);
extern int output_modes(int p1, int p2, int p3);
extern int read_modes(int p1, int p2, int p3);
extern int destroy_modes();

extern int create_modes_tri();
extern double get_modes_tri(int i);
extern int update_modes_tri(double *results);
extern int output_modes_tri();
extern int read_modes_tri();

extern int create_modesR();
extern double get_modesR_TTT(int i);
extern double get_modesR_TTE(int i);
extern double get_modesR_TEE(int i);
extern double get_modesR_EEE(int i);
extern int update_modesR(double *results);

extern int create_modesR_tri();
extern double get_modesR_tri(int i);
extern int update_modesR_tri(double *results);

extern int create_transfer(int l_size, int t_size);
extern int update_transfer(int size, double *results);
extern double get_transfer(int i, int j);

extern int create_fisher_tri(int xsize);
extern int update_fisher_tri(double *results);
extern double get_fisher_tri(int i, int j);

extern int tri(int i,int j,int k);
extern int quad(int i,int j,int k,int l);
extern double permsix(int i, int j, int k);
extern void fourfig(int i, char* suffix);

 // arrays
extern void ivector_write(int *n, char filename[MAXLEN], int grid[(*n)]);
extern void ivector_read(int *n, char filename[MAXLEN], int grid[(*n)]);
extern int *create_ivector(int length);
extern void destroy_ivector(int *vector);
extern void vector_write(int *n, char filename[MAXLEN], double grid[(*n)]);
extern void vector_read(int *n, char filename[MAXLEN], double grid[(*n)]);
extern double *create_vector(int length);
extern void destroy_vector(double *vector);
extern int **create_iarray(int dim_x, int dim_y);
extern void destroy_iarray(int **array);
extern double **create_array(int dim_x, int dim_y);
extern void destroy_array(double **array);
extern double ***create_3Darray(int dim_x, int dim_y, int dim_z);
extern double ***create_3Darray_long(long int dim_x, long int dim_y, long int dim_z);
extern void destroy_3Darray(double ***array);
extern void destroy_3Darray_long(double ***array);
extern double ****create_4Darray(int dim_w, int dim_x, int dim_y, int dim_z);
extern void destroy_4Darray(double ****array);
extern double ****create_4Darray_long(long int dim_w, long int dim_x, long int dim_y, long int dim_z);
extern void array_write(int *n, char filename[MAXLEN], double grid[(*n)]);
extern void array_write_long(long int *n, char filename[MAXLEN], double grid[(*n)]);
extern void array_read(int *n, char filename[MAXLEN], double grid[(*n)]);
extern void array_read_long(long int *n, char filename[MAXLEN], double grid[(*n)]);
extern int array_readcol(int n, char filename[MAXLEN], int nrow, int ncol, double *ptr);
extern int array_readrow(int n, char filename[MAXLEN], int nrow, int ncol, double *ptr);
extern int array_readpoints(int l0, int l1, int l2, int l3, char filename[MAXLEN], int nrow, int ncol, one_points *points);

// timing
extern double csecond();

// checkpointing
extern void ckpt_handle_sigurg(int sig);
extern int ckpt_check();
extern void ckpt_block();
extern void ckpt_unblock();
extern void ckpt_reset();
extern void ckpt_register();
extern int ckpt_write(int* lvec, double** array, int n, int rank);
extern int ckpt_load(int* lvec, double** array, int n, int rank);

// io
extern int load_txt_int(char* filename, int columns, int* values, int* size);
extern int load_txt_dbl(char* filename, int columns, double* values, int* size);
extern int load_txt_dbl_long(char* filename, int columns, double* values, long int* size);
extern int load_one(char* filename, int* values, int* size);
extern int load_one_double(char* filename, double* values, int* size);
extern int load_one_double_long(char* filename, double* values, long int* size);
extern int load_two(char* filename, double* values, int* size);
extern int load_three(char* filename, double* values, int* size);
extern int load_three_int(char* filename, int* values, int* size);
extern int load_four(char* filename, double* values, int* size);
extern int load_four_int(char* filename, int* values, int* size);
extern int load_six(char* filename, double* values, int* size);
extern int load_transfer();
extern int load_flat();
extern int load_bessel();
extern int load_limber();
extern double get_kmax();
extern double get_fmax();
extern double get_ftau();
extern double get_tau0();
extern double get_xmax();
extern double get_tmax();
extern int get_tlsize();
extern int get_flsize();
extern int get_blsize();
extern int get_lmlsize();
extern int get_ksize();
extern int get_fsize();
extern int get_xsize();
extern int get_tausize();
extern int get_kvec(double *vec);
extern int get_fkvec(double *vec);
extern int get_xvec(double *vec);
extern int get_tauvec(double *vec);
extern int get_tvec_T(int l, double *vec);
extern int get_tvec_E(int l, double *vec);
extern int get_fvec(int l, double *vec);
extern int get_bvec(int l, double *vec);
extern int get_lmvec(int l, double *vec);

// MPI
extern int sync_tasks(int mode, int size);
extern int get_next_task(int mode, int size, int rank, double *results);
extern int signal_end_tasks(int mode, int size);
extern int record_tasks(int ntasks);

#endif
