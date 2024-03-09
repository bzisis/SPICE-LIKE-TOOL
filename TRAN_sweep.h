#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include "csparse.h"
#include <gsl/gsl_linalg.h>
#include <math.h>
//#include "netlist_parser.h"

gsl_matrix *left_hand_BE;   // (G + ((1/h)*C))
gsl_vector *right_hand_BE;  // (e(tk) + ((1/h)*C*x(tk-1)))

gsl_matrix *alt_left_hand_BE;

gsl_vector *BE_x_curr;
gsl_vector *BE_x_prev;

gsl_permutation *permutation_BE;

cs *sparse_left_hand_BE;
cs *sparse_right_hand_BE;

css *tran_S;
csn *tran_N;
double *tran_y;

double *SPARSE_BE_x_curr;
double *SPARSE_BE_x_prev;

int k_tran_variable;

int k_pulse_variable;

gsl_matrix *tran_M_iterative;
gsl_matrix *tran_M_inverse_iterative;

gsl_vector *tran_r_iterative;

gsl_vector *tran_p_iterative;
gsl_vector *tran_p_scale_beta;
gsl_vector *tran_p_scale_alpha;

gsl_vector *tran_q_iterative;
gsl_vector *tran_q_scale_alpha;

gsl_vector *tran_z_iterative;

gsl_vector *tran_Ax_iterative;
gsl_vector *tran_bminusAx_iterative;

gsl_permutation *tran_permut_iterative;

gsl_vector *tran_r_tilde_iterative;
gsl_vector *tran_p_tilde_iterative;
gsl_vector *tran_p_tilde_scale_beta;
gsl_vector *tran_q_tilde_iterative;
gsl_vector *tran_z_tilde_iterative;

double *tran_sparse_y_iterative;
double *tran_sparse_r_iterative;
double *tran_sparse_r_tilde_iterative;
double *tran_sparse_p_iterative;
double *tran_sparse_p_tilde_iterative;
double *tran_sparse_q_iterative;
double *tran_sparse_q_tilde_iterative;
double *tran_sparse_z_iterative;
double *tran_sparse_z_tilde_iterative;
double *tran_sparse_m_iterative;
double *tran_sparse_right_hand_BE;

cs *tran_m_cmprsd_iterative;
cs *tran_m_cc_iterative;

///////////////////////////////////////////////////////////////////////////////////////////////////

gsl_matrix *left_hand_TR;   // (G + ((2 / h) * C))
gsl_vector *right_hand_TR;  // e(tk) + e(tk-1) - (((G - (2 / h) * C) * x(tk-1)))

gsl_matrix *alt_left_hand_TR;

gsl_vector *TR_x_curr;
gsl_vector *TR_x_prev;

gsl_vector *tran_MNA_e_prev;

gsl_matrix *right_hand_TR_sub;

gsl_matrix *alt_tran_MNA_G;

gsl_permutation *permutation_TR;

cs *sparse_left_hand_TR;
cs *sparse_right_hand_TR;

css *tran_S_TR;
csn *tran_N_TR;
double *tran_y_TR;

double *SPARSE_TR_x_curr;
double *SPARSE_TR_x_prev;

double *tran_SPARSE_e_prev;

double *tran_sparse_right_hand_TR;

cs *alt_tran_sparse_right_hand_TR;

void init_BE_hands();
void init_TR_hands();

void choose_tran_method();

void trapezoidal_method();
void SPARSE_trapezoidal_method();

void backwardeuler_method();
void SPARSE_backwardeuler_method();

double EXP_lookup(double, double, double, double, double, double, double, double);
double SIN_lookup(double, double, double, double, double, double, double, double);
double PULSE_lookup(double, double, double, double, double, double, double, double);

void TRAN_sweep_voltagesource(double);
void SPARSE_TRAN_sweep_voltagesource(double);

void TRAN_sweep_currentsource(int, int, double);
void SPARSE_TRAN_sweep_currentsource(int, int, double);

void tran_LU_sweep(int, int);
void tran_Cholesky_sweep(int, int);
void tran_BiCG_sweep(int, int);
void tran_CG_sweep(int, int);

void TR_tran_LU_sweep(int, int);
void TR_tran_Cholesky_sweep(int, int);
void TR_tran_BiCG_sweep(int, int);
void TR_tran_CG_sweep(int, int);

void tran_sparse_inverse_matrix();
void tran_sparse_preconditioner_solve();
void tran_sparse_transpose_preconditioner();

void tran_preconditioner_solve();
void tran_transpose_preconditioner_solve();

void tran_SPARSE_LU_sweep(int, int);
void tran_SPARSE_Cholesky_sweep(int, int);
void tran_SPARSE_BiCG_sweep(int, int);
void tran_SPARSE_CG_sweep(int, int);

void TR_tran_SPARSE_LU_sweep(int, int);
void TR_tran_SPARSE_Cholesky_sweep(int, int);
void TR_tran_SPARSE_BiCG_sweep(int, int);
void TR_tran_SPARSE_CG_sweep(int, int);

void tran_inverse_matrix();

void free_BE_hands();

void free_TR_hands();