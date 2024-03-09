#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>

//   global variables //

double **array_A;

double *array_b;

double *array_x;

int rowsA, columnsA;

int rowsb;

int rowsx;

int k_variable;

double **array_C;

//   function declarations                     //

void init_global_arrays();

void global_arrays_allocation();

void print_global_arrays(double **, int, int);

void print_global_vectors(double *, int);

//  MNA Algorithm                              //

void MNA_resistance(int, int, double);

void MNA_currentsource(int, int, double);

void MNA_voltagesource(int, int, double);

void MNA_inductance(int, int, double);

void MNA_algorithm();

void destroy_global_arrays(double **, int);

void destroy_global_vectors(double *);
/////////////////////////////////////////////////
void TRAN_inductance(int, int, double);

void TRAN_capacitance(int, int, double);