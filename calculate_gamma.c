#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "global.h"
#include "offload_util.h"
//#include "mkl.h" //comment out, otherwise I get compilation errors saying that Cblas has already been declared

_OFFLOADABLE static int npt;
_OFFLOADABLE static int lsize_T;
_OFFLOADABLE static int lsize_E;
_OFFLOADABLE static int lmax;
_OFFLOADABLE static int npt_pad;
_OFFLOADABLE static int lsize_T_pad;
_OFFLOADABLE static int lsize_E_pad;
_OFFLOADABLE static int lmaxPlus1_pad;

static double *TTweight;
static double *gl_pts;
//static double **gl_Pl;
_OFFLOADABLE static double *gl_Pl_flat;
_OFFLOADABLE static double *gl_wgt;
//static double **basis_T;
_OFFLOADABLE static double *basis_T_flat;
//static double ***qtilde_T;
_OFFLOADABLE static double *qtilde_T_flat;
static double *EEweight;
//static double **basis_E;
_OFFLOADABLE static double *basis_E_flat;
//static double ***qtilde_E;
_OFFLOADABLE static double *qtilde_E_flat;
static double *lweight3D;
//static int xsize;
static double *xvec;
static double xmin, xmax;
static double ***qtilde3D_T;
static double ***qtilde3D_E;
//static double ****Mij_T;
//static double ****Mij_E;
_OFFLOADABLE static double *Mij_T_flat;
_OFFLOADABLE static double *Mij_E_flat;


void init_gamma_glint(){

    double s1,s2,s3,s4,s5,s6,s7,s8,s9;	
    int i,j,k,l;

    int lmax_T = eflag_T_lmax;
    int lmax_E = eflag_E_lmax;
    int lmin_T = eflag_T_lmin;
    int lmin_E = eflag_E_lmin;
    lmax = lmax_T;
    if(do_polarisation==1){
        if(lmax_E>lmax_T)lmax = lmax_E;
    }

    lsize_T = lmax_T-lmin_T+1;
    int lvec_T[lsize_T];
    for(i=0;i<lsize_T;i++)lvec_T[i] = i+lmin_T;

    int xsize = get_qtilde_xsize();

    npt = 3*lmax/2 + 1;

    lmaxPlus1_pad = ((lmax+1)+7) & ~7;
    lsize_T_pad = (lsize_T + 7) & ~7;
    npt_pad = (npt + 7) & ~7;

    TTweight = (double *)create_vector(lsize_T);

    for(i=0;i<lsize_T;i++){
        l = lvec_T[i];
        if(l > 1){
            s1 = pow(2e0*l+1e0,5e0/6e0);
            s2 = get_beam_TT(l)*get_beam_TT(l)*get_cl_TT(l)+get_noise_TT(l);
            s3 = get_beam_TT(l) / sqrt(s2);
            TTweight[i] =  s1*s3;
        } else {
            TTweight[i] = 0e0;
        }
        //printf("L = %d\ts1 = %e\ts2 = %e\ts3 = %e\tsqrt(s2) = %e\n",l,s1,s2,s3, sqrt(s2));
        //printf("1 / sqrt(s2) = %e\t",1/sqrt(s2));
    }
    gsl_integration_glfixed_table *table =  gsl_integration_glfixed_table_alloc(npt);

    gl_pts = (double *)create_vector(npt);
    gl_wgt = (double *) _mm_malloc(npt * sizeof(double), 64);

    asy(gl_pts,gl_wgt,npt);

    //gl_Pl = (double **)create_array(npt,lmax+1);
    gl_Pl_flat = (double*) _mm_malloc(npt*(lmaxPlus1_pad)*sizeof(double), 64);
    double (*restrict gl_Pl)[lmaxPlus1_pad] = (double (*restrict)[lmaxPlus1_pad]) gl_Pl_flat;

    for(i=0;i<npt;i++){
        for(l=0;l<lmax+1;l++){
            gl_Pl[i][l] = gsl_sf_legendre_Pl(l,gl_pts[i]);
        }
    }

    int pmax = get_pmax_prim();
    int pmax_T = get_pmax_late_T();
    int pmax_E;
    if(do_polarisation==1)
        pmax_E = get_pmax_late_E();

    //basis_T = (double **)create_array(pmax_T+1,lsize_T);
    basis_T_flat = (double*) _mm_malloc((pmax_T+1)*lsize_T_pad*sizeof(double), 64);
    double (*restrict basis_T)[lsize_T_pad] = (double (*restrict)[lsize_T_pad]) basis_T_flat;
    for(l=0;l<lsize_T;l++){
        for(i=0;i<pmax_T+1;i++){
            basis_T[i][l] = get_basis_late_T(l,i);
        }
    }

    //qtilde_T = (double ***)create_3Darray(xsize,pmax+1,lsize_T);
    qtilde_T_flat = (double *) _mm_malloc(xsize*(pmax+1)*lsize_T_pad*sizeof(double), 64);
    double (*restrict qtilde_T)[pmax+1][lsize_T_pad] = (double (*restrict)[pmax+1][lsize_T_pad]) qtilde_T_flat;
    for(l=0;l<lsize_T;l++){
        for(i=0;i<xsize;i++){
            for(j=0;j<pmax+1;j++){
                qtilde_T[i][j][l] = TTweight[l]*get_qtilde_T(j,lvec_T[l],i);
                //printf("TTweight[l = %d] = %e\tqtilde_T[%d][%d][%d] = %e\n", l, TTweight[l], i,j,l,qtilde_T[i][j][l]);
            }
        }
    }

    if(do_polarisation==1){


        lsize_E = lmax_E-lmin_E+1;
        int lvec_E[lsize_E];
        for(i=0;i<lsize_E;i++)lvec_E[i] = i+lmin_E;

        EEweight = (double *)create_vector(lsize_E);

        for(i=0;i<lsize_E;i++){
            l = lvec_E[i];
            if(l > 1){
                s1 = pow(2e0*l+1e0,5e0/6e0);
                s2 =get_beam_TT(l)*get_beam_TT(l)*get_cl_TT(l)+get_noise_TT(l);
                s3 =get_beam_TT(l)*get_beam_EE(l)*get_cl_TE(l)+get_noise_TE(l);
                s4 =get_beam_EE(l)*get_beam_EE(l)*get_cl_EE(l)+get_noise_EE(l);
                s5 = get_beam_EE(l)/sqrt(s2*(s2*s4-s3*s3));
                EEweight[i] =  s1*s5;

				//printf("PS-E:\t%d\t%e\t%e\t%e\n",l,s2,s3,s4);
            } else {
                EEweight[i] = 0e0;
            }
        }

        //basis_E = (double **)create_array(pmax_E+1,lsize_E);
        lsize_E_pad = (lsize_E + 7) & ~7;
        basis_E_flat = (double*) _mm_malloc((pmax_E+1)*lsize_E_pad*sizeof(double), 64);
        double (*restrict basis_E)[lsize_E_pad] = (double (*restrict)[lsize_E_pad]) basis_E_flat;
        for(l=0;l<lsize_E;l++){
            for(i=0;i<pmax_E+1;i++){
                basis_E[i][l] = get_basis_late_E(l,i);
            }
        }


        pmax = get_pmax_prim();
        //qtilde_E = (double ***)create_3Darray(xsize,pmax+1,lsize_E);
        qtilde_E_flat = (double*) _mm_malloc(xsize*(pmax+1)*lsize_E_pad*sizeof(double),64);
        double (*restrict qtilde_E)[pmax+1][lsize_E_pad] = (double (*restrict)[pmax+1][lsize_E_pad]) qtilde_E_flat;
        for(k=0;k<lsize_E;k++){
            l = lvec_E[k];
            s1 = get_beam_TT(l)*get_beam_TT(l)*get_cl_TT(l)+get_noise_TT(l);
            s2 = sqrt(fskyE/fskyT)*get_beam_TE(l)*get_beam_TE(l)*get_cl_TE(l)+get_noise_TE(l);
            s3 = EEweight[k];
            for(i=0;i<xsize;i++){
                for(j=0;j<pmax+1;j++){
                    qtilde_E[i][j][k] = s3*(s1*get_qtilde_E(j,l,i)-s2*get_qtilde_T(j,l,i));
                }
            }
        }

    }

    //Mij_T = (double ****)create_4Darray(pmax+1,pmax_T+1,xsize,npt);
    Mij_T_flat = (double *) _mm_malloc((pmax+1)*(pmax_T+1)*xsize*npt_pad*sizeof(double), 64);
    if(do_polarisation==1){
        //Mij_E = (double ****)create_4Darray(pmax+1,pmax_E+1,xsize,npt);
        Mij_E_flat = (double *) _mm_malloc((pmax+1)*(pmax_E+1)*xsize*npt_pad*sizeof(double), 64);
    }


    double (*restrict Mij_T)[pmax_T+1][xsize][npt_pad] = (double (*restrict)[pmax_T+1][xsize][npt_pad]) Mij_T_flat;
    for(i=0;i<pmax+1;i++){
        for(j=0;j<pmax_T+1;j++){
            for(k=0;k<xsize;k++){
                for(l=0;l<npt;l++){
                    Mij_T[i][j][k][l] = 0e0;
                }
            }
        }
        if(do_polarisation==1){
            double (*restrict Mij_E)[pmax_E+1][xsize][npt_pad] = (double (*restrict)[pmax_E+1][xsize][npt_pad]) Mij_E_flat;
            for(j=0;j<pmax_E+1;j++){
                for(k=0;k<xsize;k++){
                    for(l=0;l<npt;l++){
                        Mij_E[i][j][k][l] = 0e0;
                    }
                }
            }
        }
    }

/*  
    Petar:
    In the calculations of Mij_T, the code was returning zeros,
    in paticular the term proportional to qtilde was zero when calling 
    the function precompute_gammaMij_T(). It seems reasonable NOT to 
    destroy the qtilde in this function call, since we will need it later too.
    So let's comment it out: (Alex said actually let's leave it...)        
*/
    destroy_qtilde();

    /*
    // Comment this out to make it work - Alex

    #pragma offload_transfer target(mic:offload_target) \
        in(eflag_T_lmin, eflag_E_lmin) \
        in(npt, npt_pad, lsize_T, lmax, lsize_T_pad, lmaxPlus1_pad) \
        in(Mij_T_flat[0:(pmax+1)*(pmax_T+1)*xsize*npt_pad] : ALLOC RETAIN) \
        in(qtilde_T_flat[0:xsize*(pmax+1)*lsize_T_pad] : ALLOC RETAIN) \
        in(gl_Pl_flat[0:npt*lmaxPlus1_pad] : ALLOC RETAIN) \
        in(gl_wgt[0:npt] : ALLOC RETAIN) \
        in(basis_T_flat[0:(pmax_T+1)*lsize_T_pad] : ALLOC RETAIN)

    if (do_polarisation)
    {
        #pragma offload_transfer target(mic:offload_target) \
            in(lsize_E, lsize_E_pad) \
            in(Mij_E_flat[0:(pmax+1)*(pmax_E+1)*xsize*npt_pad] : ALLOC RETAIN) \
            in(qtilde_E_flat[0:xsize*(pmax+1)*lsize_E_pad] : ALLOC RETAIN) \
            in(basis_E_flat[0:(pmax_E+1)*lsize_E_pad] : ALLOC RETAIN)
    }
    */

    //printf("Mij_FLAT = %e\t%e\t%e\t\n", Mij_T_flat[0], Mij_T_flat[1],Mij_T_flat[2]);
    return;
}


_OFFLOADABLE
void precompute_gammaMij_T(){
    int xsize = get_qtilde_xsize();
    int pmax = get_pmax_prim();
    int pmax_T = get_pmax_late_T();
    double (*restrict Mij_T)[pmax_T+1][xsize][npt_pad] = (double (*restrict)[pmax_T+1][xsize][npt_pad]) Mij_T_flat;
    double (*restrict qtilde_T)[pmax+1][lsize_T_pad] = (double (*restrict)[pmax+1][lsize_T_pad]) qtilde_T_flat;
    double (*restrict gl_Pl)[lmaxPlus1_pad] = (double (*restrict)[lmaxPlus1_pad]) gl_Pl_flat;
    double (*restrict basis_T)[lsize_T_pad] = (double (*restrict)[lsize_T_pad]) basis_T_flat;

    #pragma omp parallel
    {
        double* vec = (double*) _mm_malloc(lsize_T * sizeof(double), 64);

        int jj; 
        const int JBF=8;

        int p,q,i,j,l;
        #pragma omp for
        for (jj = 0; jj < npt; jj += JBF)
        {
            int jmax = (jj + JBF > npt) ? npt : jj + JBF;
            for (p = 0; p < pmax+1; p++) {
                for (q = 0; q < pmax_T+1; q++) {
                    for (i = 0; i < xsize; i++) {
                        #pragma vector aligned
                        for (l = 0; l < lsize_T; l++) {
                            vec[l] = qtilde_T[i][p][l]*basis_T[q][l];
                            //printf("qtilde_T[i = %d][p = %d][l = %d] = %e\t", i,p,l,qtilde_T[i][p][l]); 
                            //printf("basis_T[q = %d][l = %d] = %e\n", q,l,basis_T[q][l]);
                        }
                        //printf("\n\n");
                        for (j = jj; j < jmax; j++) {
                            double sum = 0;
                            #pragma vector unaligned
                            for ( l = 0; l < lsize_T; l++){
                                sum += gl_Pl[j][l+eflag_T_lmin]*vec[l];
                                //printf("vec [l = %d] = %e\t", l, vec[l]);  //this seems to be 0
                                //printf("gl_Pl[j = %d][l + eflag_T_lmin = %d] = %e\n", j, l+eflag_T_lmin, gl_Pl[j][l+eflag_T_lmin]);
                            }
                            Mij_T[p][q][i][j] = sum;
                            //printf("Mij_T = %e\t", sum);
                        }
                    }
                }
            }
        }
        _mm_free(vec);
    }
    return;
}

_OFFLOADABLE
void precompute_gammaMij_E(){
    int xsize = get_qtilde_xsize();
    int pmax = get_pmax_prim();
    int pmax_E = get_pmax_late_E();
    double (*restrict Mij_E)[pmax_E+1][xsize][npt_pad] = (double (*restrict)[pmax_E+1][xsize][npt_pad]) Mij_E_flat;
    double (*restrict qtilde_E)[pmax+1][lsize_E_pad] = (double (*restrict)[pmax+1][lsize_E_pad]) qtilde_E_flat;
    double (*restrict gl_Pl)[lmaxPlus1_pad] = (double (*restrict)[lmaxPlus1_pad]) gl_Pl_flat;
    double (*restrict basis_E)[lsize_E_pad] = (double (*restrict)[lsize_E_pad]) basis_E_flat;

    #pragma omp parallel
    {
        double* vec = (double*) _mm_malloc(lsize_E_pad * sizeof(double), 64);

        int jj; 
        const int JBF=8;

        int p,q,i,j,l;
        #pragma omp for
        for (jj = 0; jj < npt; jj += JBF)
        {
            int jmax = (jj + JBF > npt) ? npt : jj + JBF;
            for (p = 0; p < pmax+1; p++) {
                for (q = 0; q < pmax_E+1; q++) {
                    for (i = 0; i < xsize; i++) {

                        #pragma vector aligned
                        for (l = 0; l < lsize_E; l++) {
                            vec[l] = qtilde_E[i][p][l]*basis_E[q][l];
                        }

                        for (j = jj; j < jmax; j++) {
                            double sum = 0;
                            #pragma vector unaligned
                            for ( l = 0; l < lsize_E; l++) {
                                sum += gl_Pl[j][l+eflag_E_lmin]*vec[l];
                            }
                            Mij_E[p][q][i][j] = sum;
                        }
                    }
                }
            }
        }
        _mm_free(vec);
    }
    return;
}


void calculate_gammaMij_T(int i, int j){

    double sum;
    int xsize = get_qtilde_xsize();
    int pmax = get_pmax_prim();
    int pmax_T = get_pmax_late_T();
    double (*restrict Mij_T)[pmax_T+1][xsize][npt_pad] = (double (*restrict)[pmax_T+1][xsize][npt_pad]) Mij_T_flat;
    double (*restrict qtilde_T)[pmax+1][lsize_T_pad] = (double (*restrict)[pmax+1][lsize_T_pad]) qtilde_T_flat;
    double (*restrict gl_Pl)[lmaxPlus1_pad] = (double (*restrict)[lmaxPlus1_pad]) gl_Pl_flat;
    double (*restrict basis_T)[lsize_T_pad] = (double (*restrict)[lsize_T_pad]) basis_T_flat;

    #pragma omp parallel private(sum)
    {
        int r,s,l;
        //double *vec = create_vector(lsize_T);
        double* vec = (double*) _mm_malloc(lsize_E_pad * sizeof(double), 64);
        #pragma omp for
        for(r=0;r<xsize;r++){
            for(l=0;l<lsize_T;l++){
                vec[l] = qtilde_T[r][i][l]*basis_T[j][l];
            }
            for(s=0;s<npt;s++){
                sum = 0e0;
                for(l=0;l<lsize_T;l++){
                    sum += gl_Pl[s][l+eflag_T_lmin]*vec[l];
                }
                Mij_T[i][j][r][s] = sum;
            }
        }
        _mm_free(vec);
    }
    return;
}

void calculate_gammaMij_E(int i, int j){

    double sum;
    int xsize = get_qtilde_xsize();

    int pmax = get_pmax_prim();
    int pmax_E = get_pmax_late_E();
    double (*restrict Mij_E)[pmax_E+1][xsize][npt_pad] = (double (*restrict)[pmax_E+1][xsize][npt_pad]) Mij_E_flat;
    double (*restrict qtilde_E)[pmax+1][lsize_E_pad] = (double (*restrict)[pmax+1][lsize_E_pad]) qtilde_E_flat;
    double (*restrict gl_Pl)[lmaxPlus1_pad] = (double (*restrict)[lmaxPlus1_pad]) gl_Pl_flat;
    double (*restrict basis_E)[lsize_E_pad] = (double (*restrict)[lsize_E_pad]) basis_E_flat;

    #pragma omp parallel private(sum)
    {
        //double *vec = create_vector(lsize_E);
        double* vec = (double*) _mm_malloc(lsize_E_pad * sizeof(double), 64);
        int r,s,l;

        #pragma omp for
        for(r=0;r<xsize;r++){
            for(l=0;l<lsize_E;l++){
                vec[l] = qtilde_E[r][i][l]*basis_E[j][l];
            }
            for(s=0;s<npt;s++){
                sum = 0e0;
                for(l=0;l<lsize_E;l++){
                    sum += gl_Pl[s][l+eflag_E_lmin]*vec[l];
                }
                Mij_E[i][j][r][s] = sum;
            }
        }
        _mm_free(vec);
    }

    return;
}

void sync_gammaMij(){

    int i,j,k,l,m,n;
    double* Mij_T_send;
    double* Mij_T_recv;
    double* Mij_E_send;
    double* Mij_E_recv;
    int Tsize, Esize;
    int pmax,pmax_T,pmax_E;
    pmax = get_pmax_prim();
    pmax_T = get_pmax_late_T();
    if(do_polarisation==1)pmax_E = get_pmax_late_E();
    int xsize = get_qtilde_xsize();

    Tsize = (pmax_T+1)*xsize*npt;
    Esize = (pmax_E+1)*xsize*npt;
    Mij_T_send = (double *)malloc( Tsize * sizeof(double) );
    Mij_T_recv = (double *)malloc( Tsize * sizeof(double) );
    if(do_polarisation==1){
        Mij_E_send = (double *)malloc( Esize * sizeof(double) );
        Mij_E_recv = (double *)malloc( Esize * sizeof(double) );
    }
    double (*restrict Mij_T)[pmax_T+1][xsize][npt_pad] = (double (*restrict)[pmax_T+1][xsize][npt_pad]) Mij_T_flat;
    double (*restrict Mij_E)[pmax_E+1][xsize][npt_pad] = (double (*restrict)[pmax_E+1][xsize][npt_pad]) Mij_E_flat;

    for(i=0;i<pmax+1;i++){
        m=0;
        n=0;
        for(j=0;j<pmax_T+1;j++){
            for(k=0;k<xsize;k++){
                for(l=0;l<npt;l++){
                    Mij_T_send[m] = Mij_T[i][j][k][l];
                    m++;
                }
            }
        }
        if(do_polarisation==1){
            for(j=0;j<pmax_E+1;j++){
                for(k=0;k<xsize;k++){
                    for(l=0;l<npt;l++){
                        Mij_E_send[n] = Mij_E[i][j][k][l];
                        n++;
                    }
                }
            }
        }

        MPI_Allreduce(Mij_T_send,Mij_T_recv,Tsize,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
        if(do_polarisation==1){
            MPI_Allreduce(Mij_E_send,Mij_E_recv,Esize,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
        }

        m=0;
        n=0;
        for(j=0;j<pmax_T+1;j++){
            for(k=0;k<xsize;k++){
                for(l=0;l<npt;l++){
                    Mij_T[i][j][k][l] = Mij_T_recv[m];
                    m++;
                }
            }
        }

        if(do_polarisation==1){
            for(j=0;j<pmax_E+1;j++){
                for(k=0;k<xsize;k++){
                    for(l=0;l<npt;l++){
                        Mij_E[i][j][k][l] = Mij_E_recv[n];
                        n++;
                    }
                }
            }
        }
    }

    free(Mij_T_send);
    free(Mij_T_recv);
    // 	destroy_3Darray(qtilde3D_T);
    if(do_polarisation==1){
        free(Mij_E_send);
        free(Mij_E_recv);
        // 		destroy_3Darray(qtilde3D_E);
    }
    return;
}

_OFFLOADABLE
double gamma_2d_pt_TTT(int m, int n, int i){
    int l,r,s,j;
    int p1,p2,p3;
    int qt1,qt2,qt3,qt4,qt5,qt6;
    int pmax_T = get_pmax_late_T();
    int xsize = get_qtilde_xsize();

    double (*restrict Mij_T)[pmax_T+1][xsize][npt_pad] = (double (*restrict)[pmax_T+1][xsize][npt_pad]) Mij_T_flat;

    find_perm_prim(n,&p1,&p2,&p3);
    find_perm_late_TTT(m,&qt1,&qt2,&qt3);

    double result = 0.0;

    #pragma vector aligned
    for (j = 0; j < npt; j++){
        double sum1 = Mij_T[p1][qt1][i][j];
        double sum2 = Mij_T[p1][qt2][i][j];
        double sum3 = Mij_T[p1][qt3][i][j];
        double sum4 = Mij_T[p2][qt1][i][j];
        double sum5 = Mij_T[p2][qt2][i][j];
        double sum6 = Mij_T[p2][qt3][i][j];
        double sum7 = Mij_T[p3][qt1][i][j];
        double sum8 = Mij_T[p3][qt2][i][j];
        double sum9 = Mij_T[p3][qt3][i][j];
        //printf("sum1 --- sum9 = %e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9);

        double s1  = sum1*sum5*sum9;
        double s2  = sum1*sum6*sum8;
        double s3  = sum2*sum4*sum9;
        double s4  = sum2*sum6*sum7;
        double s5  = sum3*sum4*sum8;
        double s6  = sum3*sum5*sum7;

        //printf("s1 --- s6 = %e\t%e\t%e\t%e\t%e\t%e\n", s1,s2,s3,s4,s5,s6);


        double sum = s1+s2+s3+s4+s5+s6;
        result += gl_wgt[j]*sum;
    }
    return result;

    //printf("GAMMA_2D_PT_TTT = %e\n", result);
}

_OFFLOADABLE
double gamma_2d_pt_TTE(int m, int n, int i){
    int l,r,s,j;
    int p1,p2,p3;
    int qt1,qt2,qt3,qt4,qt5,qt6;
    int qe1,qe2,qe3,qe4,qe5,qe6;
    int pmax_T = get_pmax_late_T();
    int pmax_E = get_pmax_late_E();
    int xsize = get_qtilde_xsize();

    double (*restrict Mij_T)[pmax_T+1][xsize][npt_pad] = (double (*restrict)[pmax_T+1][xsize][npt_pad]) Mij_T_flat;
    double (*restrict Mij_E)[pmax_E+1][xsize][npt_pad] = (double (*restrict)[pmax_E+1][xsize][npt_pad]) Mij_E_flat;

    find_perm_prim(n,&p1,&p2,&p3);
    find_perm_late_TTE(m,&qt1,&qt2,&qe3);

    double result = 0.0;

    #pragma vector aligned
    for (j = 0; j < npt; j++){
        double sum1 = Mij_T[p1][qt1][i][j];
        double sum2 = Mij_T[p1][qt2][i][j];
        double sum3 = Mij_E[p1][qe3][i][j];
        double sum4 = Mij_T[p2][qt1][i][j];
        double sum5 = Mij_T[p2][qt2][i][j];
        double sum6 = Mij_E[p2][qe3][i][j];
        double sum7 = Mij_T[p3][qt1][i][j];
        double sum8 = Mij_T[p3][qt2][i][j];
        double sum9 = Mij_E[p3][qe3][i][j];

        double s1  = sum1*sum5*sum9;
        double s2  = sum1*sum6*sum8;
        double s3  = sum2*sum4*sum9;
        double s4  = sum2*sum6*sum7;
        double s5  = sum3*sum4*sum8;
        double s6  = sum3*sum5*sum7;

        double sum = s1+s2+s3+s4+s5+s6;
        result += gl_wgt[j]*sum;
    }
    return result;
}

_OFFLOADABLE
double gamma_2d_pt_TEE(int m, int n, int i) 
{
    int l,r,s,j;
    int p1,p2,p3;
    int qt1,qt2,qt3,qt4,qt5,qt6;
    int qe1,qe2,qe3,qe4,qe5,qe6;
    int pmax_T = get_pmax_late_T();
    int pmax_E = get_pmax_late_E();
    int xsize = get_qtilde_xsize();

    double (*restrict Mij_T)[pmax_T+1][xsize][npt_pad] = (double (*restrict)[pmax_T+1][xsize][npt_pad]) Mij_T_flat;
    double (*restrict Mij_E)[pmax_E+1][xsize][npt_pad] = (double (*restrict)[pmax_E+1][xsize][npt_pad]) Mij_E_flat;

    find_perm_prim(n,&p1,&p2,&p3);
    find_perm_late_TEE(m,&qt1,&qe2,&qe3);

    double result = 0.0;

    #pragma vector aligned
    for (j = 0; j < npt; j++) 
    {
        double sum1 = Mij_T[p1][qt1][i][j];
        double sum2 = Mij_E[p1][qe2][i][j];
        double sum3 = Mij_E[p1][qe3][i][j];
        double sum4 = Mij_T[p2][qt1][i][j];
        double sum5 = Mij_E[p2][qe2][i][j];
        double sum6 = Mij_E[p2][qe3][i][j];
        double sum7 = Mij_T[p3][qt1][i][j];
        double sum8 = Mij_E[p3][qe2][i][j];
        double sum9 = Mij_E[p3][qe3][i][j];

        double s1  = sum1*sum5*sum9;
        double s2  = sum1*sum6*sum8;
        double s3  = sum2*sum4*sum9;
        double s4  = sum2*sum6*sum7;
        double s5  = sum3*sum4*sum8;
        double s6  = sum3*sum5*sum7;

        double sum = s1+s2+s3+s4+s5+s6;
        result += gl_wgt[j]*sum;
    }
    return result;
}

_OFFLOADABLE
double gamma_2d_pt_EEE(int m, int n, int i) 
{
    int l,r,s,j;
    int p1,p2,p3;
    int qt1,qt2,qt3,qt4,qt5,qt6;
    int qe1,qe2,qe3,qe4,qe5,qe6;
    int pmax_T = get_pmax_late_T();
    int pmax_E = get_pmax_late_E();
    int xsize = get_qtilde_xsize();

    double (*restrict Mij_E)[pmax_E+1][xsize][npt_pad] = (double (*restrict)[pmax_E+1][xsize][npt_pad]) Mij_E_flat;

    find_perm_prim(n,&p1,&p2,&p3);
    find_perm_late_EEE(m,&qe1,&qe2,&qe3);

    double result = 0.0;

    #pragma vector aligned
    for (j = 0; j < npt; j++) 
    {
        double sum1 = Mij_E[p1][qe1][i][j];
        double sum2 = Mij_E[p1][qe2][i][j];
        double sum3 = Mij_E[p1][qe3][i][j];
        double sum4 = Mij_E[p2][qe1][i][j];
        double sum5 = Mij_E[p2][qe2][i][j];
        double sum6 = Mij_E[p2][qe3][i][j];
        double sum7 = Mij_E[p3][qe1][i][j];
        double sum8 = Mij_E[p3][qe2][i][j];
        double sum9 = Mij_E[p3][qe3][i][j];

        double s1  = sum1*sum5*sum9;
        double s2  = sum1*sum6*sum8;
        double s3  = sum2*sum4*sum9;
        double s4  = sum2*sum6*sum7;
        double s5  = sum3*sum4*sum8;
        double s6  = sum3*sum5*sum7;

        double sum = s1+s2+s3+s4+s5+s6;
        result += gl_wgt[j]*sum;
    }
    return result;
}

void gamma_pt(int m, int n, int i, double* ypt){

    int l,r,s,j;
    int p1,p2,p3;
    int qt1,qt2,qt3,qt4,qt5,qt6;
    int qe1,qe2,qe3,qe4,qe5,qe6;
    double s1,s2,s3,s4,s5,s6;
    double wgt;
    double sum, factor;
    int pmax_T = get_pmax_late_T();
    int pmax_E = get_pmax_late_E();
    int xsize = get_qtilde_xsize();
    double (*restrict Mij_T)[pmax_T+1][xsize][npt_pad] = (double (*restrict)[pmax_T+1][xsize][npt_pad]) Mij_T_flat;
    double (*restrict Mij_E)[pmax_E+1][xsize][npt_pad] = (double (*restrict)[pmax_E+1][xsize][npt_pad]) Mij_E_flat;

    find_perm_prim(n,&p1,&p2,&p3);

    find_perm_late_TTT(m,&qt1,&qt2,&qt3);
    if(do_polarisation==1){
        find_perm_late_TTE(m,&qt4,&qt5,&qe6);
        find_perm_late_TEE(m,&qt6,&qe4,&qe5);
        find_perm_late_EEE(m,&qe1,&qe2,&qe3);
    }

    ypt[0] = 0e0;
    if(do_polarisation==1){
        ypt[1] = 0e0;
        ypt[2] = 0e0;
        ypt[3] = 0e0;
    }

    for(j=0;j<npt;j++){

        wgt = gl_wgt[j];

        s1 = Mij_T[p1][qt1][i][j]*Mij_T[p2][qt2][i][j]*Mij_T[p3][qt3][i][j];
        s2 = Mij_T[p2][qt1][i][j]*Mij_T[p3][qt2][i][j]*Mij_T[p1][qt3][i][j];
        s3 = Mij_T[p3][qt1][i][j]*Mij_T[p1][qt2][i][j]*Mij_T[p2][qt3][i][j];
        s4 = Mij_T[p3][qt1][i][j]*Mij_T[p2][qt2][i][j]*Mij_T[p1][qt3][i][j];
        s5 = Mij_T[p2][qt1][i][j]*Mij_T[p1][qt2][i][j]*Mij_T[p3][qt3][i][j];
        s6 = Mij_T[p1][qt1][i][j]*Mij_T[p3][qt2][i][j]*Mij_T[p2][qt3][i][j];

        sum = s1+s2+s3+s4+s5+s6;
        ypt[0] += wgt*sum;

        if(do_polarisation==1){

            s1 = Mij_T[p1][qt4][i][j]*Mij_T[p2][qt5][i][j]*Mij_E[p3][qe6][i][j];
            s2 = Mij_T[p2][qt4][i][j]*Mij_T[p3][qt5][i][j]*Mij_E[p1][qe6][i][j];
            s3 = Mij_T[p3][qt4][i][j]*Mij_T[p1][qt5][i][j]*Mij_E[p2][qe6][i][j];
            s4 = Mij_T[p3][qt4][i][j]*Mij_T[p2][qt5][i][j]*Mij_E[p1][qe6][i][j];
            s5 = Mij_T[p2][qt4][i][j]*Mij_T[p1][qt5][i][j]*Mij_E[p3][qe6][i][j];
            s6 = Mij_T[p1][qt4][i][j]*Mij_T[p3][qt5][i][j]*Mij_E[p2][qe6][i][j];

            sum = s1+s2+s3+s4+s5+s6;
            ypt[1] += wgt*sum;

            s1 = Mij_T[p1][qt6][i][j]*Mij_E[p2][qe4][i][j]*Mij_E[p3][qe5][i][j];
            s2 = Mij_T[p2][qt6][i][j]*Mij_E[p3][qe4][i][j]*Mij_E[p1][qe5][i][j];
            s3 = Mij_T[p3][qt6][i][j]*Mij_E[p1][qe4][i][j]*Mij_E[p2][qe5][i][j];
            s4 = Mij_T[p3][qt6][i][j]*Mij_E[p2][qe4][i][j]*Mij_E[p1][qe5][i][j];
            s5 = Mij_T[p2][qt6][i][j]*Mij_E[p1][qe4][i][j]*Mij_E[p3][qe5][i][j];
            s6 = Mij_T[p1][qt6][i][j]*Mij_E[p3][qe4][i][j]*Mij_E[p2][qe5][i][j];

            sum = s1+s2+s3+s4+s5+s6;
            ypt[2] += wgt*sum;

            s1 = Mij_E[p1][qe1][i][j]*Mij_E[p2][qe2][i][j]*Mij_E[p3][qe3][i][j];
            s2 = Mij_E[p2][qe1][i][j]*Mij_E[p3][qe2][i][j]*Mij_E[p1][qe3][i][j];
            s3 = Mij_E[p3][qe1][i][j]*Mij_E[p1][qe2][i][j]*Mij_E[p2][qe3][i][j];
            s4 = Mij_E[p3][qe1][i][j]*Mij_E[p2][qe2][i][j]*Mij_E[p1][qe3][i][j];
            s5 = Mij_E[p2][qe1][i][j]*Mij_E[p1][qe2][i][j]*Mij_E[p3][qe3][i][j];
            s6 = Mij_E[p1][qe1][i][j]*Mij_E[p3][qe2][i][j]*Mij_E[p2][qe3][i][j];

            sum = s1+s2+s3+s4+s5+s6;
            ypt[3] += wgt*sum;
        }

    }

    return;
}

_OFFLOADABLE
void calculate_gamma_2d(int m, int n, double *result){

    //printf("CALCULATE_GAMMA_2D CALLED\n");

    int xsize = get_qtilde_xsize();
    double *xvec = (double*) alloca(xsize*sizeof(double));
    get_qtilde_xvec(xvec);

    double *yvec_TTT;
    double *yvec_TTE;
    double *yvec_TEE;
    double *yvec_EEE;

    yvec_TTT = (double*) alloca(xsize*sizeof(double));

    if (do_polarisation)
    {
        yvec_TTE = (double*) alloca(xsize*sizeof(double));
        yvec_TEE = (double*) alloca(xsize*sizeof(double));
        yvec_EEE = (double*) alloca(xsize*sizeof(double));
    }

    double xmin = xvec[0];
    double xmax = xvec[xsize-1];

    int i;
    for (i = 0; i < xsize; i++)
    {
        yvec_TTT[i] = xvec[i]*xvec[i]*gamma_2d_pt_TTT(m,n,i);
        //printf("(1) yvec_TTT[%d] = %e\n", i, yvec_TTT[i]);
        if( i<10 && fabs(yvec_TTT[i]) > 1e-18 )
            yvec_TTT[i]=0;
        //printf("(2) yvec_TTT[%d] = %e\n", i, yvec_TTT[i]);
    }

    if (do_polarisation){
        for (i = 0; i < xsize; i++)
        {
            yvec_TTE[i] = xvec[i]*xvec[i]*gamma_2d_pt_TTE(m,n,i);
            if( i<10 && fabs(yvec_TTE[i]) > 1e-18 )
                yvec_TTE[i]=0;
        }
        for (i = 0; i < xsize; i++)
        {
            yvec_TEE[i] = xvec[i]*xvec[i]*gamma_2d_pt_TEE(m,n,i);
            if( i<10 && fabs(yvec_TEE[i]) > 1e-18 )
                yvec_TEE[i]=0;
        }
        for (i = 0; i < xsize; i++)
        {
            yvec_EEE[i] = xvec[i]*xvec[i]*gamma_2d_pt_EEE(m,n,i);
            if( i<10 && fabs(yvec_EEE[i]) > 1e-18 )
                yvec_EEE[i]=0;
        }
    }

    result[0] = integrate_cubic_spline(xvec, yvec_TTT, xsize);
    result[0] *= deltaphi*deltaphi/(8.0*M_PI);

    //printf("(3) RESULT = %e\n", result[0]);

    if (do_polarisation){
        result[1] = integrate_cubic_spline(xvec, yvec_TTE, xsize);
        result[1] *= deltaphi*deltaphi/(8.0*M_PI);
        result[2] = integrate_cubic_spline(xvec, yvec_TEE, xsize);
        result[2] *= deltaphi*deltaphi/(8.0*M_PI);
        result[3] = integrate_cubic_spline(xvec, yvec_EEE, xsize);
        result[3] *= deltaphi*deltaphi/(8.0*M_PI);
    }
}

_OFFLOADABLE
double integrate_cubic_spline(double *x, double *y, int npts)
{
    double *spl[5];	 // cubic spline coefficients at each point
    int i;

    // sanity check
    if(npts < 1) return 0.0;
    if(npts == 1) return 0.0;

    // allocate spl
    spl[0] = malloc(npts * sizeof *spl[0]); //x 
    spl[1] = malloc(npts * sizeof *spl[0]); //y
    spl[2] = malloc(npts * sizeof *spl[0]); //a
    spl[3] = malloc(npts * sizeof *spl[0]); //b
    spl[4] = malloc(npts * sizeof *spl[0]); //c

    // create the cubic polynml coefficients at each point
    for (i=0; i<npts; i++)
    {   
        spl[0][i] = x[i];
        spl[1][i] = y[i];
    }   

    spl[4][0] = spl[0][1] - spl[0][0];
    spl[3][1] = (spl[1][1] - spl[1][0])/spl[4][0];

    for(i=1; i<npts-1; i++)
    {   
        spl[4][i] = x[i+1] - spl[0][i];
        spl[2][i] = 2.0*(spl[4][i-1] + spl[4][i]);
        spl[3][i+1] = (spl[1][i+1] - spl[1][i])/spl[4][i]; 
        spl[3][i] = spl[3][i+1] - spl[3][i];
    }   

    spl[2][0] = -spl[4][0];
    spl[2][npts-1] = -spl[4][npts-2];
    spl[3][0] = 0.0;
    spl[3][npts-1] = 0.0;

    if(npts>2)
    {
        spl[3][0] = spl[3][2]/(spl[0][3]-spl[0][1]) - spl[3][1]/(spl[0][2]-spl[0][0]);
        spl[3][npts-1] = spl[3][npts-2]/(spl[0][npts-1]-spl[0][npts-3]) - spl[3][npts-3]/(spl[0][npts-2]-spl[0][npts-4]);
        spl[3][0] = spl[3][0]*spl[4][0]*spl[4][0] / (spl[0][3]-spl[0][0]);
        spl[3][npts-1] = -spl[3][npts-1]*spl[4][npts-2]*spl[4][npts-2] / (spl[0][npts-1]-spl[0][npts-4]);
    }

    for (i=1; i<npts; i++)
    {
        double h = spl[4][i-1]/spl[2][i-1];
        spl[2][i] = spl[2][i] - h*spl[4][i-1];
        spl[3][i] = spl[3][i] - h*spl[3][i-1];
    }

    spl[3][npts-1] = spl[3][npts-1]/spl[2][npts-1];
    for (i=npts-2; i >= 0; i--)
    {
        spl[3][i] = (spl[3][i] - spl[4][i]*spl[3][i+1])/spl[2][i];
    }

    spl[2][npts-1] = (spl[1][npts-1] - spl[1][npts-2])/spl[4][npts-2] + spl[4][npts-2]*(spl[3][npts-2] + 2.0*spl[3][npts-1]);

    for (i=0; i<npts-1; i++)
    {
        spl[2][i] = (spl[1][i+1] - spl[1][i])/spl[4][i] - spl[4][i]*(spl[3][i+1] + 2.0*spl[3][i]);
        spl[4][i] = (spl[3][i+1] - spl[3][i])/spl[4][i];
        spl[3][i] = 3.0*spl[3][i];
    }

    spl[3][npts-1] = 3.0*spl[3][npts-1];
    spl[4][npts-1] = spl[4][npts-2];

    /* integrate with the cubic spline */
    double intgrl = 0.0;

    for(i=0; i<npts-1; i++)
    {
        double delX = x[i+1] - x[i];
        intgrl += delX*(spl[1][i] + delX*(spl[2][i]/2.0 + delX*(spl[3][i]/3.0 + delX*spl[4][i]/4.0)) );
    }   

    // free the cubic spline
    free(spl[0]);
    free(spl[1]);
    free(spl[2]);
    free(spl[3]);
    free(spl[4]);

    return intgrl;
}

void calculate_gamma(int m, int n, double *result) {

    int i;

    int xsize = get_qtilde_xsize();
    double *xvec = create_vector(xsize);
    get_qtilde_xvec(xvec);
    double *yvec_TTT;
    double *yvec_TTE;
    double *yvec_TEE;
    double *yvec_EEE;

    if(do_polarisation==1){
        yvec_TTT = create_vector(xsize);
        yvec_TTE = create_vector(xsize);
        yvec_TEE = create_vector(xsize);
        yvec_EEE = create_vector(xsize);
    }else{
        yvec_TTT = create_vector(xsize);
    }

    double xmin = xvec[0];
    double xmax = xvec[xsize-1];

    #pragma omp parallel
    {
        double *ypt;
        if(do_polarisation==1){
            ypt = (double*)calloc(4,sizeof(double));
        }else{
            ypt = (double*)calloc(1,sizeof(double));
        }
        #pragma omp for
        for(i=0;i<xsize;i++){
            gamma_pt(m,n,i,ypt);
            yvec_TTT[i] = xvec[i]*xvec[i]*ypt[0];
            if(i<10&&fabs(yvec_TTT[i])>1e-18) yvec_TTT[i]=0e0;
            if(do_polarisation==1){
                yvec_TTE[i] = xvec[i]*xvec[i]*ypt[1];
                yvec_TEE[i] = xvec[i]*xvec[i]*ypt[2];
                yvec_EEE[i] = xvec[i]*xvec[i]*ypt[3];
                if(i<10&&fabs(yvec_TTE[i])>1e-18) yvec_TTE[i]=0e0;
                if(i<10&&fabs(yvec_TEE[i])>1e-18) yvec_TEE[i]=0e0;
                if(i<10&&fabs(yvec_EEE[i])>1e-18) yvec_EEE[i]=0e0;
            }
        }
        free(ypt);
    }

    gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, xsize);
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline_init(sp,xvec,yvec_TTT,xsize);
    result[0] = gsl_spline_eval_integ(sp,xmin,xmax,acc);
    if(do_polarisation==1){
        gsl_spline_init(sp,xvec,yvec_TTE,xsize);
        result[1] = gsl_spline_eval_integ(sp,xmin,xmax,acc);
        gsl_spline_init(sp,xvec,yvec_TEE,xsize);
        result[2] = gsl_spline_eval_integ(sp,xmin,xmax,acc);
        gsl_spline_init(sp,xvec,yvec_EEE,xsize);
        result[3] = gsl_spline_eval_integ(sp,xmin,xmax,acc);
    }
    gsl_spline_free(sp);
    gsl_interp_accel_free(acc);

    result[0] *= deltaphi*deltaphi/(8e0*M_PI);
    if(do_polarisation==1){
        result[1] *= deltaphi*deltaphi/(8e0*M_PI);
        result[2] *= deltaphi*deltaphi/(8e0*M_PI);
        result[3] *= deltaphi*deltaphi/(8e0*M_PI);
    }

    free(xvec);
    free(yvec_TTT);
    if(do_polarisation==1){
        free(yvec_TTE);
        free(yvec_TEE);
        free(yvec_EEE);
    }

    return;
}

