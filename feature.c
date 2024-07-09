#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

typedef struct {
        double a1;
        double a2;
        double a3;
} shape_params;

// Python script to call (no need to change)
#define CMD "python ./SRFT_call.py"
#define MAX_LINE_LEN 500

// Prototypes
int cs_feature_prior_check(shape_params params);
void cs_feature_init(shape_params params, int *n_total, double **k, double **dpp);
double cs_feature(double k1, double k2, double k3, shape_params params);

// Checking that the parameters are within the prior bounds
int cs_feature_prior_check(shape_params params) {

  // Get the parameters
  double log10u = params.a1;
  double log10s = params.a2;
  double tau0 = params.a3;
  double aux, logbeta, minsr;
  aux = pow(10., log10u)/(exp(0.5)+pow(10., log10u));
  logbeta = log(2.*pow((pow(10., log10s)/aux), 2));
  minsr = log10(0.03);

  if (
      // hard limits
      log10u < -4 || log10u > 0 ||
      log10s < -2 || log10s > 0 ||
      tau0 > 0 ||
      // slow-roll limit
      (log10u < minsr && log10s < minsr) ||
      // beta limit
      logbeta < 0 || logbeta > 14) {
    return 0;
  };
  return 1;
}


// Initialisation: get the vector of samples of deltaP/P
void cs_feature_init(shape_params params, int *n_total, double **k, double **dpp) {

  // Get the parameters
  double log10u = params.a1;
  double log10s = params.a2;
  double tau0 = params.a3;
  
  // Prepare the command
  char cmd_args[MAX_LINE_LEN];
  sprintf(cmd_args, "%s %g %g %g", CMD, log10u, log10s, tau0);

  // Prepare the vars and arrays to track and store the output
  int i = 0, status;
  char line[MAX_LINE_LEN];
  int n_guess = 100; // Initial guess for the output size.
                     // If too small, arrays are doubled on the fly.
  double *tmp_k, *tmp_dpp, k_i, dpp_i;
  *k   = (double *)malloc(n_guess*sizeof(double));
  *dpp = (double *)malloc(n_guess*sizeof(double));
  // Launch the process
  FILE *process;
  printf("# Calling command: %s\n", cmd_args);
  process = popen(cmd_args, "r");
  if (process == NULL) {
  printf("ERROR: failed to set up the environment when running '%s'\n", cmd_args);
  }

  // Read the output from the buffer
  while (fgets(line, sizeof(line)-1, process) != NULL) {
      sscanf(line, "%lf %lf", &k_i, &dpp_i);
      if ((i+1) > n_guess) {
        n_guess *= 2;
        tmp_k   = (double *)realloc(*k,   n_guess*sizeof(double));
        tmp_dpp = (double *)realloc(*dpp, n_guess*sizeof(double));
        if (tmp_k == NULL || tmp_dpp == NULL) {
            printf("ERROR: allocation failed!");
        };
        *k   = tmp_k;
        *dpp = tmp_dpp;
      };
      (*k  )[i] = k_i;
      (*dpp)[i] = dpp_i;
      i++;
  };
  *n_total = i;

  // Close the process (or check if it failed)
  status = pclose(process);
  if (status != 0) {
    printf("ERROR: Something failed along the process.\n");
    printf("Run it by hand to test it: '%s'",cmd_args);
  };

  // Uncomment for tests.
  // Output must be the same as the python script when called directly
  /*
  int j;
  for (j=0; j<i;j++) {printf("%d : %15.8e %15.8e\n", j, (*k)[j], (*dpp)[j]);};
  */
  
};


// Get B(k1,k2,k3). First call initialises the shape 
double cs_feature(double k1, double k2, double k3, shape_params params) {

  // Static set of parameters and spline helpers
  static shape_params params_current = {.a1 = 0, .a2 = 0, .a3 = 0};
  static gsl_interp_accel *acc = NULL;
  static gsl_spline *spline = NULL;
  static double k_min_spline, k_max_spline;
  
  // If the parameters changed or this is the 1st call, initialise
  if (params.a1 != params_current.a1 ||
      params.a2 != params_current.a2 ||
      params.a3 != params_current.a3) {
    printf("# Params changed or 1st call -- Initialising...\n");
    // Check the prior constraints, and fail if they fall outside the prior
    if (cs_feature_prior_check(params) == 0) {
      printf("ERROR: Outside prior bounds!\n");
    };
    // Release the memory used by the interpolation, if 1st call
    if (params_current.a1 != 0 &&
        params_current.a2 != 0 &&
        params_current.a3 != 0) {
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
    };
    // Define the new interpolating function
    int n_samples;
    double *k, *dpp;
    cs_feature_init(params, &n_samples, &k, &dpp);
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc(gsl_interp_cspline, n_samples);
    gsl_spline_init(spline, k, dpp, n_samples);
    // Finally, save the parameters so we know the spline is initialised
    params_current = params;
    k_min_spline = k[0];
    k_max_spline = k[n_samples-1];
    free(k);
    free(dpp);
    printf("# Interpolating spline created!\n");
  };
  
  // Evaluate the shape at the given k's
  double kt, kto2, ki2, ki12, ki123;
  kt   = k1+k2+k3;
  kto2 = kt/2.;  
  if (kto2 <= k_min_spline || kto2 >= k_max_spline) {
    return 0;
  };
  ki2   = pow(k1,2)+pow(k2,2)+pow(k3,2);
  ki12  = k1*k2 + k1*k3 + k2*k3;
  ki123 = k1*k2*k3;
  double c0, c1, c2;
  c0 = pow(kt,-2) * (-1.)   * ki12 +
       pow(kt,-1) * ( 1./4.) * (pow(k1,4)+pow(k2,4)+pow(k3,4))/ki123 +
       pow(kt,-1) * (-3./2.) * (k1*k2/k3+k1*k3/k2+k2*k3/k1) +
       1.         * (-5./4.) +
       pow(kt, 1) * ( 1./4.) * ki12/ki123;
  c1 = pow(kt,-1) * ( 1./2.) * ki12 +
       1.         * (-1./8.) * ((ki2-pow(k1,2))/k1 + (ki2-pow(k2,2))/k2 + (ki2-pow(k3,2))/k3);
  c2 = ki2        * (1./16.);

  return c0 * gsl_spline_eval(spline, kto2, acc) +
         c1 * gsl_spline_eval_deriv(spline, kto2, acc) +
         c2 * gsl_spline_eval_deriv2(spline, kto2, acc);    
};


// Your call goes here
int main() {

 shape_params params;
 params.a1 = -3.; // log10u
 params.a2 = -0.8; // log10s
 params.a3 = -1200; // tau0

 // Check the prior constraints, and fail if they fall outside the prior
 if (cs_feature_prior_check(params) == 0) {
   printf("ERROR: Outside prior bounds!\n");
 };
 /* Actually, this function is already called inside the cs_feature,
 but you probably prefer to make the check here, since if that fails
 once cs_feature has been called, the only effect is that a more or
 less flat B (modulo numerical errors) is returned, which may be
 problematic later.
 Better to check it first and avoid those points! */

 // Examples 
 double S, k1, k2, k3, k_min=0.0001, k_max=0.1;
 int i, N=1200;
 printf("# Equilateral direction\n");
 for (i=1; i<=N; i++) {
    k1 = (k_min+(double)i/(double)N*(k_max-k_min))/3.;
    k2 = (k_min+(double)i/(double)N*(k_max-k_min))/3.;
    k3 = (k_min+(double)i/(double)N*(k_max-k_min))/3.;
    S  = cs_feature(k1, k2, k3, params);
    printf("%15.8e %15.8e %15.8e %15.8e\n", k1, k2, k3, S);
  };
 printf("# Folded direction\n");
 for (i=1; i<=N; i++) {
    k1 = (k_min+(double)i/(double)N*(k_max-k_min))/2.;
    k2 = (k_min+(double)i/(double)N*(k_max-k_min))/4.;
    k3 = (k_min+(double)i/(double)N*(k_max-k_min))/4.;
    S  = cs_feature(k1, k2, k3, params);
    printf("%15.8e %15.8e %15.8e %15.8e\n", k1, k2, k3, S);
  };
 printf("# Squeezed direction\n");
 for (i=1; i<=N; i++) {
    k1 = ((double)i/(double)N*(k_max-k_min))/2.;
    k2 = ((double)i/(double)N*(k_max-k_min))/2.;
    k3 = k_min;
    S  = cs_feature(k1, k2, k3, params);
    printf("%15.8e %15.8e %15.8e %15.8e\n", k1, k2, k3, S);
  };
 
}
