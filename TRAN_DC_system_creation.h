#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include "csparse.h"
#include <gsl/gsl_linalg.h>

gsl_matrix *tran_MNA_G;

gsl_matrix *tran_MNA_C;

gsl_vector *tran_MNA_X0;

gsl_vector *tran_MNA_e;

//cs *tran_SPARSE_tilde_G;

//cs *tran_SPARSE_tilde_C;

double *tran_SPARSE_X0;

double *tran_SPARSE_e;

int tran_MNA_rows;
int tran_MNA_columns;

int tran_SPARSE_rows;
int tran_SPARSE_columns;

void init_tran_MNA_matrices();

void init_tran_SPARSE_matrices();

void choose_tran_solution_method();

void free_tran_MNA_matrices();

void free_tran_SPARSE_matrices();