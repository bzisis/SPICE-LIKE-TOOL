#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include "csparse.h"
#include <gsl/gsl_linalg.h>
//#include <suitesparse/SuiteSparse_config.h>
//#include <suitesparse/cs.h>

css *sparse_S;

csn *sparse_N;

double *sparse_y;

double *sparse_b_plot;
double *sparse_b_plot_temp;

void sparse_LU_decomp();

void sparse_cholesky_decomp();

void choose_sparse_decomp(char *);

void sparse_cholesky_DC_op(char *);

void sparse_LU_DC_op(char *);

void sparse_dc_sweep();