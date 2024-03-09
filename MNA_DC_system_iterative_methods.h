#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

gsl_matrix *A_iterative;    //  ((n - 1) + m2) x ((n - 1) + m2)

gsl_vector *X_iterative;    //  ((n - 1) + m2) x 1

gsl_vector *B_iterative;    //  ((n - 1) + m2) x 1

gsl_vector *r_iterative;    //  ((n - 1) + m2) x 1  ->  (B - Ax)

gsl_vector *p_iterative;    //  ((n - 1) + m2) x 1

gsl_vector *p_scale_beta;

gsl_vector *p_scale_alpha;

gsl_vector *q_iterative;    //  ((n - 1) + m2) x 1  ->  (A * p)

gsl_vector *q_scale_alpha;

gsl_vector *z_iterative;    //  ((n - 1) + m2) x 1  ->  M * z = r => z = (M ^ (-1)) * r

gsl_matrix *M_iterative;    //  M = diag(A)

gsl_matrix *M_inverse_iterative;    // (M ^ (-1))

gsl_permutation *permut_iterative;

gsl_vector *p_tilde_iterative;

gsl_vector *z_tilde_iterative;

gsl_vector *r_tilde_iterative;

gsl_vector *q_tilde_iterative;

gsl_vector *p_scale_beta_tilde;

gsl_vector *Ax_iterative;

gsl_vector *bminusAx_iterative;

void choose_iterative_decomp(char *);

void CG_decomp();

void BiCG_decomp();

void preconditioner_solve();

void transpose_preconditioner_solve();

void inverse_matrix();

void free_CG();

void free_BiCG();

void CG_DC_op(char *);

void BiCG_DC_op(char *);