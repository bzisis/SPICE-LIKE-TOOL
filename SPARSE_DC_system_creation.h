#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include "csparse.h"
//#include <suitesparse/cs.h>
//#include <suitesparse/cs.h>

//   global variables //

int non_zeros;

cs *sparse_A;

cs *sparse_C;

double *sparse_b;

int sparse_rowsA, sparse_columnsA;

int sparse_k_variable;

int sparse_counter;

cs *tran_sparse_C;

cs *tran_sparse_tilde_C;

cs *tran_SPARSE_tilde_G;

int non_zeros_tran_C;

int tran_sparse_counter;

//   function declarations                     //

void init_global_sparse_vars();

void count_non_zeros();

void sparse_resistance(int, int);

void sparse_voltagesource(int, int);

void TRAN_sparse_capacitance(int, int);

void TRAN_sparse_inductance();

void triplet_format_method();

void triplet_format_resistance(int, int, double);

void triplet_format_voltagesource(int, int, double);

void triplet_format_currentsource(int, int, double);

void triplet_format_inductance(int, int, double);

void sparse_triplet_format_capacitance(int, int, double);

void destroy_sparse_matrices();