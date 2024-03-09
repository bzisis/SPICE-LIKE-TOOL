#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include "csparse.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
//#include <suitesparse/SuiteSparse_config.h>
//#include <suitesparse/cs.h>

double *sparse_x_iterative;

double *sparse_y_iterative;

double *sparse_r_iterative;
double *sparse_r_tilde_iterative;

double *sparse_q_iterative;
double *sparse_q_tilde_iterative;

double *sparse_p_iterative;
double *sparse_p_tilde_iterative;

double *sparse_m_iterative;

double *sparse_z_iterative;
double *sparse_z_tilde_iterative;

cs *sparse_m_cmprsd_iterative;

cs *sparse_m_cc_iterative;

gsl_matrix *sparse_M_iterative;

gsl_matrix *sparse_M_inverse_iterative;

gsl_permutation *sparse_permut_iterative;

void choose_sparse_iterative_decomp(char *);

void sparse_BiCG_decomp();

void sparse_BiCG_DC_op(char *);

void sparse_free_BiCG();

void sparse_CG_decomp();

void sparse_CG_DC_op(char *);

void sparse_free_CG();

void sparse_preconditioner_solve();

void sparse_transpose_preconditioner_solve();

double dot_product(double *, double *, int);

void sparse_inverse_matrix();

void jacobiPreconditioner(const cs *, double *);