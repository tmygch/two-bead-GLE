//* library

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>

/*
 * Public/reproducibility version.
 *
 * The original research code used Numerical Recipes routines (ran2, gasdev,
 * and fit).  They are not distributed with this source.  Random numbers are
 * generated here with the C++ standard library, and the unweighted linear
 * regression used below is implemented directly.  Consequently, trajectories
 * will not be bitwise identical to those generated with the original RNG,
 * although the stochastic model is unchanged.
 */

//* macros
//** general setting

#define PI     3.141592653589793
#define PI2    (2.0 * PI)

//** function like macros

#define SQ(a)         ((a)*(a))
#define CU(a)         ((a)*(a)*(a))

#define MAX(x,y)      (((x) > (y) ? (x) : (y)))
#define MIN(x,y)      (((x) < (y) ? (x) : (y)))

//** system parameters

#define N       6             // # of segments
#ifndef MMM1
#define MMM1    0             // override at compile, e.g. -DMMM1=1
#endif
#define TN1     MMM1            // tagged segment
#define TN2     5               // tagged segment
#define TN1P    (TN1+1)         // tagged segment (index for theory)
#define TN2P    (TN2+1)         // tagged segment (index for theory)

static_assert(TN1 >= 0 && TN1 < N, "TN1 must be a valid bead index");
static_assert(TN2 >= 0 && TN2 < N, "TN2 must be a valid bead index");
static_assert(TN1 != TN2, "TN1 and TN2 must be different beads");
#define SAMPLES 30
#define SAMPLES_EQL 5
#define TMSD_SAMPLES 8

#define M_END 100                // Delta = 5 x 10^4
#define T_END (M_END*3000)          // t = 1 x 10^8
#define DELTA 5                   // Δ = DELTA * TINT
#define DT    0.01
#define TINT  (DT*1.0)
#define MEAS_STEP ((int)(T_END/TINT))

#define BIN_NUM 100
#define BIN_MIN 0.0
#define BIN_MAX 0.5
#define BIN_WID ((BIN_MAX-BIN_MIN)/(double)BIN_NUM)  // bin width

#define ETA0   1.0              // viscosity
#define D      1.0              // diffusion coeff. of segment
#define ZETA   1.0              // friction coeff.
#define D2     (D * 2.0)        // diffusion coeff. for pair coordinate 
#define ZETA2  (ZETA/2.0)       // friction coeff for pair coordinate
#define B      1.0              // bond length
#define RDS    1.0              // radius of segment
#define kT     1.0              // Boltzmann factor
#define K      1.0              // spring const.

//** print settings

#define PRN 1000

//* global variables

FILE   *fp[10];

int    hist[BIN_NUM+1];
double t_ser[MEAS_STEP+DELTA+1][3];
double f_ser[MEAS_STEP+DELTA+1][3];
double e_msd[MEAS_STEP+DELTA+1];
double avdr1 = 0, avdr2 = 0;
double t_msd_tensor[3][3] = {0};
double c_fn1[MEAS_STEP/DELTA+1], w_fn1[MEAS_STEP/DELTA+1];
double c_fn2[MEAS_STEP/DELTA+1], w_fn2[MEAS_STEP/DELTA+1];
double t_msd[(int)(M_END/TINT)+DELTA+1];
double w_msd[(int)(M_END/TINT)+DELTA+1];
double v_cor[(int)(M_END/TINT)+DELTA+1];
double w_cor[(int)(M_END/TINT)+DELTA+1];
double f_cor[(int)(M_END/TINT)+DELTA+1];
double g_cor[(int)(M_END/TINT)+DELTA+1];
double mu[(int)(M_END/TINT)+DELTA+1];
double log_t_msd[M_END+1];
double tamsd[TMSD_SAMPLES+1][(int)(M_END/TINT)+DELTA+1];
double K2 = 0.6;

//* random number generator

static std::mt19937_64 rng;
static std::normal_distribution<double> normal01(0.0, 1.0);
static std::uniform_real_distribution<double> uniform01(0.0, 1.0);

double uniform_random_open () {

  double u;
  do {
    u = uniform01(rng);
  } while (u <= 0.0);
  return u;
}

//* functions
//** initial setup

void init_array () {

  for ( int i = 0; i < BIN_NUM; i++ )  hist[i] = 0;
  for ( int i = 0; i < MEAS_STEP; i++ )  e_msd[i] = 0;
  for ( int i = 0; i <= MEAS_STEP/DELTA; i++ )  w_fn1[i] = 0;
  for ( int i = 0; i <= MEAS_STEP/DELTA; i++ )  w_fn2[i] = 0;
  for ( int i = 0; i < M_END/TINT; i++ )  t_msd[i] = 0;
  for ( int i = 0; i < M_END/TINT; i++ )  w_msd[i] = 0;
  for ( int i = 0; i <= 2; i++ ) t_ser[0][i] = 0;

  if (TN1 == 0) K2= 0.375;
}

void init_print_setting () {

  fp[0] = fopen("time_series.dat", "w");
  fp[1] = fopen("t_msd.dat", "w");
  fp[2] = fopen("transport_coefficient.dat", "w");
  fp[3] = fopen("mmt_trans_coef.dat", "w");
  fp[4] = fopen("setup.dat", "w");

  fprintf( fp[4], "TN1 = %d | D= %10.8f | K2= %f\n", TN1, D, K2);
  printf("TN1 = %d | D= %10.8f | K2= %f\n", TN1, D, K2);
}

void init ( long& seed ) {

  ///// seed for random number generator
  seed = 19750927;
  //seed = 20160325;
  rng.seed(static_cast<std::mt19937_64::result_type>(seed));
  normal01.reset();

  init_array();
  init_print_setting();
}

//** vector operation

double square_vec( double x[] ){

  return (SQ(x[0]) + SQ(x[1]) + SQ(x[2]));
}

void normalize_vec ( double x[] ){

  double norm = square_vec(x);
  for ( int i = 0 ; i <= 2; i++ ) x[i] /= sqrt(norm);
}

double inner_product ( double u[], double x[] ){

  double ip = 0;
  for ( int i = 0 ; i <= 2; i++ ) ip += u[i]*x[i];

  return (ip);
}

void outer_product ( double o_prdct[3], double u[], double x[] ){

  o_prdct[0] = u[1]*x[2] - u[2]*x[1];
  o_prdct[1] = u[2]*x[0] - u[0]*x[2];
  o_prdct[2] = u[0]*x[1] - u[1]*x[0];
}

void set_vector ( double x[], double y[] ){

  for ( int i = 0 ; i <= 2; i++ ) x[i] = y[i];
}
void set_vector_diff ( double dr[], double x[], double y[] ){

  for ( int i = 0 ; i <= 2; i++ ) dr[i] = x[i] - y[i];
}

//** noise generation

double noise_exp_dist(long* a) {

  (void)a;
  return -log(uniform_random_open());
}

void noise_on_unitsphere ( double dw[], long& a ){

  (void)a;
  for ( int i = 0 ; i <= 2; i++ ) dw[i] = uniform01(rng) - 0.5;

  normalize_vec(dw);
}

void gauss_noise ( double dw[], long& a ){

  (void)a;
  for ( int i = 0 ; i <= 2; i++ )  dw[i] = normal01(rng);
}

void gauss_noise_hi ( double dw[], long& a ){

  (void)a;
  for ( int n = 0 ; n <= 3*N-1; n++ ) dw[n] = normal01(rng);
}

void thermal_noise_vec ( double noise[], double dt, long& a ){

  double dw[3];  gauss_noise(dw, a);
  for ( int i = 0; i <= 2; i++ )  noise[i] = sqrt(2 * D * dt) * dw[i] * ZETA;
}

//** center of mass

void center_of_mass ( double x[][3], double com[] ) {

  for ( int i = 0 ; i <= 2; i++ ) {

    com[i] = 0;
    for ( int n = 0 ; n < N; n++ )  com[i] += x[n][i];
    com[i] /= (double)N;
  }
}
//** initial configuration

void set_com_origin ( double x[][3], double com[3] ) {

  center_of_mass (x, com);
  for ( int n = 0 ; n < N; n++ ) {

    for ( int i = 0 ; i <= 2; i++ )  x[n][i] -= com[i];
  }
  center_of_mass (x, com);      // set com[i] = 0
}

void set_center_particle_origin ( double x[][3] ) {

  double tx[3];
  set_vector (tx, x[TN1]);
  for ( int n = 0 ; n < N; n++ ) {

    for ( int i = 0 ; i <= 2; i++ )  x[n][i] -= tx[i];
  }
}

void init_config ( double x[][3], double com[3], long& a ){

  x[0][0] = x[0][1] = x[0][2] = 0;
  for ( int n = 1 ; n < N; n++ ) {

  double dw[3]; gauss_noise (dw, a);

    for ( int i = 0 ; i <= 2; i++ )  x[n][i] = x[n-1][i] + B * dw[i];
  }
  set_center_particle_origin (x);
}

//** t_ser, f_ser,  and e_msd

int time_series (double x[][3], double com[], double t, int n_meas){

  int l  = (t+TINT/20.0) / TINT;
  double dx[3];
  if (l > n_meas) { 

    for ( int i = 0 ; i <= 2; i++ ){

      t_ser[l][i]  = x[TN1][i] - x[TN2][i];
      f_ser[l][i]  = -K2 * t_ser[l][i];
    }

    set_vector_diff(dx, x[TN1], x[TN2]);
    e_msd[l] += square_vec(dx);
    n_meas = l;
  }
  
  return (n_meas);
}

//** radius of gyration

double rad_gyration (double x[][3]){

  double rod = 0;
  for ( int m = 0; m <= N-1; m++ ) {
    for ( int n = 0; n <= N-1; n++ ) {
      for ( int i = 0 ; i <= 2; i++ ) {
        rod += SQ(x[m][i]-x[n][i]);
      }
    }
  }
  return (rod / (2*N*N));
}

//** toy6 dynamics

double L[6][6] =
  {{ 1, -1,  0,  0,  0,  0},
   {-1,  3, -1,  0, -1,  0},
   { 0, -1,  4, -1, -1, -1},
   { 0,  0, -1,  1,  0,  0},
   { 0, -1, -1,  0,  2,  0},
   { 0,  0, -1,  0,  0,  1}};

void set_force_rouse (double x[][3], double force[N][3], long& a ) {

  int k, l, i;
  for (int k = 0; k < N; k++) {
    for (int i = 0; i < 3; i++) force[k][i] = 0.0;
    for (int l = 0; l < N; l++) {
      for (int i = 0; i < 3; i++)  force[k][i] += -K * L[k][l] * x[l][i] * DT;
    }
  }
}

//** diffusion tensor

double k_delta (int i, int j) {

  int delta = 0;
  if ( i == j)  delta = 1;
  return ( (double)delta );
}

void set_diffusion_tensor ( double x[][3], double dc[3*N][3*N] ) {

  double dr[3], dr_sq; 
  for ( int i = 0 ; i < 3*N; i++ ) {
    for ( int j = 0 ; j < 3*N; j++ ) {

      ////// Rouse model
      dc[i][j] = D * k_delta (i, j);
      /*
      ////// Rotne-Prager-Yamakawa tensor
      // if i & j are coordinates from the same particles
      if ( i/3 == j/3 )  dc[i][j] = D * k_delta (i%3, j%3);

      // if i & j are coordinates from different particles
      else {                

        set_vector_diff (dr, x[i/3], x[j/3]);
        dr_sq = square_vec(dr);
        
        if ( dr_sq > SQ(2.0*RDS) ) {
        
          dc[i][j]  = k_delta (i%3, j%3) / 3.0 - dr[i%3] * dr[j%3] / dr_sq;
          dc[i][j] *= 2.0 * SQ(RDS) / dr_sq; 
          dc[i][j] += k_delta (i%3, j%3) + dr[i%3] * dr[j%3] / dr_sq;
          dc[i][j] *= kT / (8.0 * PI * ETA0 * sqrt (dr_sq));
        }
        else {

          dc[i][j]  = (1.0 - 9.0/32.0 * sqrt(dr_sq) / RDS) * k_delta (i%3, j%3);
          dc[i][j] += 3.0/32.0 * dr[i%3] * dr[j%3] / (RDS * sqrt(dr_sq));
          dc[i][j] *= kT / (6.0 * PI * ETA0 * RDS);
        }
       }
      */
    }
  }
}

//** thermal noise

void set_thermal_noise_hi ( double noise_hi[N][3], double dt, long& a ){

  double dw[3*N];  gauss_noise_hi(dw, a);
  for ( int i = 0 ; i <= 3*N-1; i++ ) noise_hi[i/3][i%3] = sqrt(2.0 * D* dt) * dw[i];
}

//** hydrodynamic forces

void set_force_hi ( double force_hi[N][3], double force[N][3], double dt ) {

  for ( int i = 0 ; i <= 3*N-1; i++ ) {

    force_hi[i/3][i%3] = D * force[i/3][i%3] / kT;
  }
}

//** time evolution and calculation of ensemble averaged msd

void euler_method (double x[][3], double com[], double& t, long& a, int& n_meas, int smpl){

  // update diffusion tensor & sigma_ij
  double force[N][3];  set_force_rouse (x, force, a);
  // hydrodynamic force & noise
  double force_hi[N][3]; set_force_hi (force_hi, force, DT);
  double noise_hi[N][3]; set_thermal_noise_hi (noise_hi, DT, a);
  for ( int n = 0 ; n <= N-1; n++ ) {

    for ( int i = 0 ; i <= 2; i++ )  x[n][i] += force_hi[n][i] + noise_hi[n][i];
  }

  ///// data analysis
  n_meas = time_series (x, com, t, n_meas);
}

void time_evolve ( double x[][3], double com[], long& a, int& n_meas, double& t,
                   FILE *fp, int smpl ) {

  if ( smpl == 2 && t <= 1000 && (int)(t*DT) % 10000 == 0) {

    for ( int n = 0 ; n < N; n++ )  fprintf ( fp, "%d %f %f %f\n", n, x[n][0], x[n][1], x[n][2]);
    fprintf (fp, "\n\n");
  }

  euler_method (x, com, t, a, n_meas, smpl);
}

//** time averaged velocity autocorrelation

void v_correlation (){

  for ( int m = 0; m < M_END/TINT; m+=(m+4)/4 ) {

    v_cor[m] = 0;
    f_cor[m] = 0;
    
    int sample; double vt1, vt2, ft1, ft2;
    for ( sample = 1; sample < (T_END - M_END)/TINT; sample++ ) {

      for ( int i = 0 ; i <= 2; i++ ) {

        vt1 = (t_ser[sample+m][i] - t_ser[sample+m-1][i]) / DT;
        vt2 = (t_ser[sample  ][i] - t_ser[sample  -1][i]) / DT;
        ft1 =  f_ser[sample+m][i];
        ft2 =  f_ser[sample  ][i];
        
        v_cor[m] += vt1 * vt2;
        f_cor[m] += ft1 * vt2; // + ft2 * vt1 - ft1 * ft2 / ZETA;
      }
    }
    v_cor[m] /= (double)(sample-1);
    f_cor[m] /= (double)(sample-1);
  }
}

void v_correlation_full (){

  for ( int m = 0; m < M_END/TINT; m++ ) {

    v_cor[m] = 0;
    f_cor[m] = 0;
    
    int sample; double vt1, vt2, ft1, ft2;
    for ( sample = 1; sample < (T_END - M_END)/TINT; sample++ ) {

      for ( int i = 0; i <= 2; i++ ) {

        vt1 = (t_ser[sample+m][i] - t_ser[sample+m-1][i]) / DT;
        vt2 = (t_ser[sample  ][i] - t_ser[sample  -1][i]) / DT;
        ft1 =  f_ser[sample+m][i];
        ft2 =  f_ser[sample  ][i];
        
        v_cor[m] += vt1 * vt2;
        f_cor[m] += ft1 * vt2; // + ft2 * vt1 - ft1 * ft2 / ZETA;
      }
    }
    v_cor[m] /= (double)(sample-1);
    f_cor[m] /= (double)(sample-1);
  }
}


//** calculation of time averaged msd

void time_av_msd ( ) {

  for ( int m = 1; m < M_END/TINT; m+=(m+4)/4 ) {

    t_msd[m] = 0;

    int sample = 0;
    for ( sample = 0; sample < (T_END - M_END)/TINT; sample++ ) {

      for ( int i = 0 ; i <= 2; i++ ) {
        t_msd[m] += SQ( t_ser[sample+m][i] - t_ser[sample][i] );
      }
    }
    t_msd[m] /= (double)sample;
  }
}

//** tmsd tensor

void tmsd_tensor ( ) {

  double dr[3];
  for ( int m = 0 ; m < MEAS_STEP-DELTA; m++ ) {

    set_vector_diff (dr, t_ser[m+DELTA], t_ser[m]);
    for ( int i = 0 ; i <= 2; i++ ) {
      for ( int j = 0 ; j <= 2; j++ ) {

        t_msd_tensor[i][j] += dr[i] * dr[j] / (double)(MEAS_STEP-DELTA);
      }
    }
  }
}

void correlation_func_avdr ( ) {

  tmsd_tensor ( );

  avdr1 = 0; avdr2 = 0;
  for ( int i = 0 ; i <= 2; i++ ) {
    
    avdr1 += t_msd_tensor[i][i];
    
    for ( int j = 0 ; j <= 2; j++ ) {
      avdr2 += t_msd_tensor[i][j] * t_msd_tensor[j][i];
    }
  }
}

void normalize_t_msd_tensor ( ) {

  for ( int i = 0 ; i <= 2; i++ ) {
    for ( int j = 0 ; j <= 2; j++ ) {
      t_msd_tensor[i][j] /= (double)SAMPLES_EQL; 
    }
  }

  avdr1 /= (double)SAMPLES_EQL;
  avdr2 /= SQ((double)SAMPLES_EQL);
}

//** calculation of correlation function 1

void correlation_func1 ( ) {

  double dr[3], integral = 0;
  for ( int m = 0 ; m <= MEAS_STEP - DELTA; m += DELTA ) {

    for ( int i = 0 ; i < DELTA; i++ ) {
      set_vector_diff (dr, t_ser[m+i+DELTA], t_ser[m+i]);
      integral += square_vec(dr);
    }

    c_fn1[m/DELTA] = SQ( integral / (double)(m+DELTA) - avdr1);
  }
}

//** calculation of correlation function 2

void tensor_product (double dr1dr2[3][3], double dr1[3], double dr2[3]) {

  for ( int i = 0 ; i <= 2; i++ ) {
    for ( int j = 0 ; j <= 2; j++ ) dr1dr2[i][j] = dr1[i] * dr2[j];
  }
}
void set_tensor_diff (double dr1dr2_c[3][3], double dr1dr2[3][3], double t_msd_tensor[3][3]) {

  for ( int i = 0 ; i <= 2; i++ ) {
    for ( int j = 0 ; j <= 2; j++ ) {

      dr1dr2_c[i][j] = dr1dr2[i][j] - t_msd_tensor[i][j];
    }
  }
}
double make_contraction (double drdr1[3][3], double drdr2[3][3]) {

  double contraction = 0;
  for ( int i = 0 ; i <= 2; i++ ) {
    for ( int j = 0 ; j <= 2; j++ ) {
      contraction += drdr1[i][j] * drdr2[j][i];
    }
  }
  return(contraction);
}

void correlation_func2 ( ) {

  double dr1[3], dr2[3], dr1dr1[3][3], dr1dr1_c[3][3], dr2dr2[3][3], dr2dr2_c[3][3], integral = 0;
  for ( int m1 = 0 ; m1 <= MEAS_STEP - DELTA; m1 += DELTA ) {
    for ( int i = 0 ; i < DELTA; i++ ) {

      set_vector_diff (dr1, t_ser[m1+i+DELTA], t_ser[m1+i]);
      tensor_product (dr1dr1, dr1, dr1);
      set_tensor_diff (dr1dr1_c, dr1dr1, t_msd_tensor);
      
      for ( int m2 = 0; m2 <= m1; m2 += DELTA ) {
        for ( int j = 0 ; j < DELTA; j++ ) {
          if (m2+j < m1+i) {
            set_vector_diff (dr2, t_ser[m2+j+DELTA], t_ser[m2+j]);
            tensor_product (dr2dr2, dr2, dr2);
            set_tensor_diff (dr2dr2_c, dr2dr2, t_msd_tensor);
            integral += make_contraction(dr1dr1_c, dr2dr2_c);
          }
          else if (m2+j == m1+i) {

            set_vector_diff (dr2, t_ser[m2+j+DELTA], t_ser[m2+j]);
            tensor_product (dr2dr2, dr2, dr2);
            set_tensor_diff (dr2dr2_c, dr2dr2, t_msd_tensor);
            integral += make_contraction(dr1dr1_c, dr2dr2_c) / 2.0;
          }
        }
      }
    }
    c_fn2[m1/DELTA] = 2.0 / SQ((double)m1+DELTA) * integral;
  }
}

//** histgram

void  histgram ( double value ) {

  if ( value < BIN_MAX && value > BIN_MIN )
    hist[(int)(( value - BIN_MIN ) / BIN_WID)] ++;
}

//** transport coefficient

/*
 * Unweighted straight-line fit y = a + b x.
 * The arrays are treated as 1-based here because the original analysis code
 * stores the fitted data in elements 1,...,ndata.  This is the only fit mode
 * used by this program.  The argument list is kept compatible with the old
 * call site; sig is unused when mwt == 0.
 */
void fit (const double x[], const double y[], int ndata, const double sig[], int mwt,
          double *a, double *b, double *siga, double *sigb, double *chi2, double *q) {

  (void)sig;
  if (mwt != 0) {
    fprintf(stderr, "Error: this public version implements only unweighted fitting (mwt=0).\n");
    exit(EXIT_FAILURE);
  }
  if (ndata < 3) {
    fprintf(stderr, "Error: at least three points are required for linear fitting.\n");
    exit(EXIT_FAILURE);
  }

  double sx = 0.0, sy = 0.0;
  for (int i = 1; i <= ndata; ++i) {
    sx += x[i];
    sy += y[i];
  }
  const double xbar = sx / static_cast<double>(ndata);
  const double ybar = sy / static_cast<double>(ndata);

  double sxx = 0.0, sxy = 0.0;
  for (int i = 1; i <= ndata; ++i) {
    const double dx = x[i] - xbar;
    sxx += dx * dx;
    sxy += dx * (y[i] - ybar);
  }
  if (sxx == 0.0) {
    fprintf(stderr, "Error: all x values are identical in linear fitting.\n");
    exit(EXIT_FAILURE);
  }

  *b = sxy / sxx;
  *a = ybar - (*b) * xbar;

  *chi2 = 0.0;
  for (int i = 1; i <= ndata; ++i) {
    const double r = y[i] - (*a + *b * x[i]);
    *chi2 += r * r;
  }

  const double sigma = sqrt(*chi2 / static_cast<double>(ndata - 2));
  *sigb = sigma / sqrt(sxx);
  *siga = sigma * sqrt(1.0 / static_cast<double>(ndata) + xbar * xbar / sxx);
  *q = 1.0;  // no weighted chi-square goodness-of-fit probability is used here
}

double mmt1=0.0, mmt2=0.0;

void transport_coefficient ( FILE *fp ) {

  double a, b, siga, sigb, chi2, q;
  double log_m[M_END], normal_m[M_END], sig[M_END];

  ///// log scale fit
  // for ( int m = 1; m < M_END; m++ ) {

  //   log_m[m] = log((double)m * TINT);
  //   log_t_msd[m] = log (t_msd[m]);
  // }

  // fit(log_m, log_t_msd, M_END-1, sig, 0, &a, &b, &siga, &sigb, &chi2, &q);

  // fprintf( fp, "%f\n", exp(a));
  // histgram (exp(a));

  ///// normal scale fit

  for ( int m = 1; m < M_END; m++ ) {

    normal_m[m] = (double)m * TINT;
  }

  fit(normal_m, t_msd, M_END-1, sig, 0, &a, &b, &siga, &sigb, &chi2, &q);

  fprintf(fp, "%f %f %f %f\n", b, sigb, a, siga);

  histgram (b);

  mmt1 += b / SAMPLES;
  mmt2 += b*b / SAMPLES;
}

//** print
//*** print time series (fp)

void print_t_ser ( FILE *fp , int n, int i, int x, double t ) {

  fprintf( fp, "%d %f %d\n", i, t, x);
}

//*** log average

void log_average ( double& av_x, double& av_y, int& pj, int j, double data[] ) {

  int num  = 0;

  av_x = 0;  av_y = 0;
  for ( int i = pj ; i <= j; i++ ) {

    num ++;
    av_x += i;
    av_y += data[i];
  }
  
  av_x /= (double) num;
  av_y /= (double) num;

  pj = j+1;
}

//*** print ensemble averaged msd (tp)

void print_e_msd ( FILE *tp , int smpl ) {

  if ( smpl % MAX(SAMPLES/PRN, 1) == 0 ) {

    tp = fopen("e_msd.dat", "w");

    double av_x, av_y;
    for ( int j = 1, pj = 1; j < T_END/TINT; j += (j+4)/4 ) {

      log_average (av_x, av_y, pj, j, e_msd);
      fprintf( tp, "%f %f\n", av_x*TINT, av_y / (double)smpl );
    }
    fclose(tp);
  }
}

//*** print time averaged msd (fp)

void print_t_msd ( FILE *fp,  FILE *tp, int smpl ) {

  if ( smpl <= SAMPLES ) {

    for ( int m = 1; m < M_END/TINT; m+=(m+4)/4 ) {

      w_msd[m] += t_msd[m];
      if ( smpl <= TMSD_SAMPLES ) {

        tamsd[smpl-1][m] = t_msd[m];
        //fprintf( fp, "%f %f\n", (double)m*TINT, t_msd[m]);
      }
    }

    if ( smpl == TMSD_SAMPLES) {

      for ( int m = 1; m < M_END/TINT; m+=(m+4)/4 ) {

        fprintf( fp, "%f ", (double)m*TINT);

        for ( int i = 0; i < smpl; i++ ) {

          fprintf( fp, "%f ", tamsd[i][m]);
        }
        fprintf( fp, "\n");
      }
      fclose (fp);
    }

  if ( smpl % MAX(SAMPLES/PRN, 1) == 0 ) {

      tp = fopen("et_msd.dat", "w");

      for ( int m = 1; m < M_END/TINT; m+=(m+4)/4 ) {

        fprintf( tp, "%f %f\n", (double)m*TINT, w_msd[m]/smpl);
      }

      fclose(tp);
    }
  }
}

//*** print time averaged vcor (fp)

void print_v_cor ( FILE *tp, FILE *tp2, int smpl ) {

  for ( int m = 0; m < M_END/TINT; m++ ) {
      
    w_cor[m] += v_cor[m] / D2;   // psi(t)
    g_cor[m] += f_cor[m] / kT; // <∇U(t)v(0)> + <∇U(0)v(t)> - <∇U(t)∇U(0)>/ZETA
  }
  w_cor[0] -= 6.0 / DT;  // remove delta part
  // simple interporation: y = (y2-y1)/(x2-x1) * (x-x1) + y1
  //w_cor[0] = w_cor[1] - (w_cor[2]-w_cor[1]);
  
  double denom = smpl * 3;  // 3 <= x, y, z
  mu[0] = -w_cor[0] / denom + g_cor[0]/ denom;
  //printf("denominator=%f\n", 1 + DT/2.0 * w_cor[0] / denom);
  for ( int m = 1; m < M_END/TINT; m++ ) {
    double sum = 0; 
    for (int k=1; k< m; k++)  sum += mu[k] * w_cor[m-k]/denom;

    //printf("%d: %f | %f | m=%d; %f \n", smpl, sum, mu[0], m-1, mu[m-1]);
    mu[m] = -w_cor[m]/denom - DT/2.0 * mu[0] * w_cor[m]/denom - sum * DT + g_cor[m]/denom;
    mu[m] /= (1 + DT/2.0 * w_cor[0] / denom);
  }
    
  if ( smpl % MAX(SAMPLES/PRN, 1) == 0 ) {

    tp = fopen("et_vcor_coarse.dat", "w");
    tp2 = fopen("et_vcor.dat", "w");

    for ( int m = 0; m < M_END/TINT; m+=(m+5)/5 ) {

      // output negative value, because the vcor is negative 
      fprintf( tp, "%f %f %f %f\n", (double)m*TINT, -w_cor[m]/denom, g_cor[m]/denom, mu[m]);
    }

    double av_x1, av_y1, av_x2, av_y2, av_x3, av_y3;
    for ( int m = 1, pm1=1, pm2=1, pm3=1; m < M_END/TINT; m+=(m+30)/30 ) {

      log_average (av_x1, av_y1, pm1, m, w_cor);
      log_average (av_x2, av_y2, pm2, m, g_cor);
      log_average (av_x3, av_y3, pm3, m, mu);
      fprintf( tp2, "%f %f %f %f\n", av_x1*TINT, -av_y1/denom, av_y2/denom, av_y3);
    }
    fclose(tp); fclose(tp2);
  }
}

//*** print fn1 (fp)

void print_et_fn1 ( FILE *tp, int smpl ) {

  if (smpl % 10 == 0) printf ("fn1: %5d | %16.13f %16.13f\n", smpl, avdr1*100,  6.0 * D / N * DELTA * DT * 100);

  for ( int m = 0; m < MEAS_STEP/DELTA; m++ ) {

    w_fn1[m] += c_fn1[m] / SQ(avdr1);
  }

  if ( smpl % MAX(SAMPLES/PRN, 1) == 0 ) {

    tp = fopen("et_fn1.dat", "w");

    double av_x, av_y;
    for ( int m = 0, pm = 0; m < MEAS_STEP/DELTA; m+=(m+2)/2 ) {

      log_average (av_x, av_y, pm, m, w_fn1);
      fprintf( tp, "%15.13f %15.13f\n", (av_x+1)*TINT*DELTA, av_y / (double)smpl);
      //fprintf( tp, "%15.13f %15.13f\n", (m+1)*TINT*DELTA, w_fn1[m]/(double)smpl);
    }

    fclose(tp);
  }
}

//*** print fn2 (fp)

void print_et_fn2 ( FILE *tp, int smpl ) {

  if (smpl % 10 == 0) printf ("fn2: %5d | %16.13f %16.13f\n", smpl, avdr2*100, 1200.0 * SQ(DELTA*TINT*D/N));
  for ( int m = 0; m < MEAS_STEP/DELTA; m++ ) {

    w_fn2[m] += c_fn2[m] / avdr2;
  }

  if ( smpl % MAX(SAMPLES/PRN, 1) == 0 ) {

    tp = fopen("et_fn2.dat", "w");

    double av_x, av_y;
    for ( int m = 0, pm = 0; m < MEAS_STEP/DELTA; m+=(m+2)/2 ) {

      log_average (av_x, av_y, pm, m, w_fn2);
      fprintf( tp, "%13.11f %13.11f\n", (av_x+1)*TINT*DELTA, av_y / (double)smpl);
      //fprintf( tp, "%15.13f %15.13f\n", (m+1)*TINT*DELTA, w_fn2[m]/(double)smpl);
    }

    fclose(tp);
  }
}

//*** print histgram (tp)

void print_histgram ( FILE *tp, FILE *tp2, int smpl ) {

  if ( smpl % 1000 == 0 ) {

    printf( "** %5.1f percent finish **\n", (double)smpl / (double)SAMPLES * 100.0);

    double all = 0, norm = 0;
    for ( int i = 0; i < BIN_NUM; i++ )  all  += (double)hist[i] * BIN_WID;
    for ( int i = 0; i < BIN_NUM; i++ )
      norm += ((double)i+0.5) * BIN_WID * (double)hist[i] * BIN_WID / all;

    tp = fopen("h.transport_coefficient.dat", "w");
    tp2 = fopen("nh.transport_coefficient.dat", "w");

    for ( int i = 0; i < BIN_NUM; i++ ) {

      fprintf( tp,  "%f %f\n", ((double)i+0.5) * BIN_WID, (double)hist[i] / all);
      fprintf( tp2, "%f %f\n",
               ((double)i+0.5) * BIN_WID / norm, norm * (double)hist[i] / all);
    }

    fclose(tp); fclose(tp2);
  }
}

//* main function

int main ( void ) {
  
  long a;  FILE *tp[10];

  ////// initial setup
  int n_meas = 0;  double x[N][3], com[3], av_pos = 0.0;

  ///// equilibration & calculation of avdr1
  init (a); init_config (x, com, a); 
  for ( int smpl = 1; smpl <= SAMPLES_EQL; smpl++ ) {

    ///// time evolution
    n_meas = 0; set_center_particle_origin (x);
    for ( double t = DT; t <= T_END+DT*DELTA; t+= DT ) {
      time_evolve (x, com, a, n_meas, t, fp[0], smpl);
    }
    correlation_func_avdr ();
    if (smpl % 1 == 0) printf ("equilibrate: %3d (%15.10f %15.10f) (%15.10f %15.10f)\n",
                               smpl, avdr1/smpl, 6.0 * D / N * DELTA * DT,
                               avdr2/SQ(smpl), 12.0 * SQ(DELTA*TINT*D/N));
  }
  normalize_t_msd_tensor();

  ///// ensembles
  int smpl_rod = 0; double rod = 0;
  init_array (); 
  for ( int smpl = 1; smpl <= SAMPLES; smpl++ ) {

    printf("K2=%f\n", K2);

    ///// time evolution
    n_meas = 0; set_center_particle_origin (x);
    //init_config (x, com, a);
    int t_step = 0;
    for ( double t = DT; t <= T_END+DT*DELTA; t+= DT ) {
      t_step ++;
      //if ((int)(t/DT)%10000 == 0) printf ("%5.1f % fin. (smpl = %d)\n", t/T_END * 100, smpl);
      time_evolve (x, com, a, n_meas, t, fp[0], smpl);
      if (t_step % 100 == 0) {rod += rad_gyration(x); smpl_rod ++;}
    }
    //fprintf( fp[0], "\n\n");
    
    av_pos += x[TN1][0] - x[TN2][0];
    time_av_msd();
    // v_correlation();
    v_correlation_full();
    // correlation_func1 ( );
    // correlation_func2 ( );
    //transport_coefficient ( fp[2] );

    print_e_msd ( tp[0], smpl );
    print_t_msd ( fp[1], tp[3], smpl );
    print_v_cor ( tp[4], tp[5], smpl );
    // print_et_fn1 ( tp[6], smpl );
    // print_et_fn2 ( tp[5], smpl );
    print_histgram ( tp[1], tp[2], smpl );

    printf("(TN1,TN2)=(%3d, %3d); sample=%5d; rad of gyraion=%8.5f\n", TN1, TN2, smpl, rod / (double)smpl_rod);
  }

  double var = mmt2-mmt1*mmt1;
  fprintf( fp[3], "%f %f %f %f %f\n", (double)T_END*TINT, mmt1, mmt2, sqrt(var) / mmt1, rod/(double)smpl_rod);
}
