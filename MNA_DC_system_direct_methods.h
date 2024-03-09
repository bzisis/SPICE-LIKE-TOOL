#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <gsl/gsl_linalg.h>

//  global vectors    //
gsl_matrix *A;

gsl_vector *X;

gsl_vector *B;

gsl_permutation *P;

void LU_decomp();

void cholesky_decomp();

void LU_DC_op(char *);

void cholesky_DC_op(char *);

void choose_decomp(char *);

void dc_sweep();

int return_newpos_node_name(char *);