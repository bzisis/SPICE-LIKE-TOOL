#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include "csparse.h"
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <complex.h>
#include <cs.h>

gsl_matrix_complex *AC_left_hand;
gsl_vector_complex *AC_right_hand;
gsl_vector_complex *AC_x;

gsl_permutation *AC_permutation;

gsl_vector_complex *AC_r_iterative;
gsl_vector_complex *AC_r_herm_iterative;

gsl_vector_complex *AC_p_iterative;
gsl_vector_complex *AC_p_herm_beta;
gsl_vector_complex *AC_p_herm_iterative;
gsl_vector_complex *AC_p_scale_iterative;

gsl_vector_complex *AC_q_iterative;
gsl_vector_complex *AC_q_herm_iterative;
gsl_vector_complex *AC_q_herm_alpha;

gsl_vector_complex *AC_z_iterative;
gsl_vector_complex *AC_z_herm_iterative;

gsl_vector_complex *AC_Ax_iterative;
gsl_vector_complex *AC_bminusAx_iterative;

gsl_matrix_complex *AC_M_iterative;
gsl_matrix_complex *AC_M_inverse_iterative;

int AC_rowsA;
int AC_columnsA;

int k_LH;
int k_RH;

void choose_AC_method();

void init_AC_matrices();

void AC_LU_sweep();

void AC_BiCG_sweep();

void AC_BiCG_solver();

void AC_voltagesource_LH(int, int);
void AC_resistance_LH(int, int, double);
void AC_capacitance_LH(int, int, double, double);
void AC_inductance_LH(int, int, double, double);

void AC_inverse_matrix();