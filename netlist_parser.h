/**************************************************/
/*         Circuit Simulation Algorithms          */
/**************************************************/ 
/* Authors: Georgia-Despoina Tzanetou             */
/*          Zisis Balatsos                        */
/*                                                */
/*       University of Thessaly Volos, 2023       */
/*                                                */
/**************************************************/

///////////////////////////////////////////////////////////////////////////////////////////////////
//   Library of the parsing tool, containing all necessary "C" libraries, structs, global        //
//   variables and function declarations for completing a successful parse of a specified        //
//   netlist file of the following form:                                                         //
//                                                                                               //
//                  ComponentTypeComponentName_PositiveNode_NegativeNode_Value                   //
//                                                                                               //
//   For example:   V1 5 0 2                                                                     //
//                                                                                               //
//   where: V = ComponentType                                                                    //
//          1 = ComponentName                                                                    //
//          5 = PositiveNode                                                                     //
//          0 = NegativeNode                                                                     //
//          2 = Value                                                                            //
//                                                                                               //
//   The main idea behind this parsing algorithm follows the next steps:                         //
//   1) Open a .txt netlist file, where each line represents a circuit component and is in the   //
//      form of the above example (each component is represented in one line).                   //
//   2) It is taken for granted that the .txt netlist file is of the mentioned form in order for //
//      the parsing to be work correctly.                                                        //
//   3) Representations:                                                                         //
//                      (i) V -> voltage source                                                  //
//                     (ii) I -> current source                                                  //
//                    (iii) R -> resistance                                                      //
//                     (iv) L -> inductance                                                      //
//                      (v) C -> capacitance                                                     //
//                     (vi) * -> comment in the beginning of a line (ignored by the parser)      //
//   (*) This tool uses non-case sensitive parsing.                                              //
//                                                                                               //
//   4) Scan each line of the file from the beginning until the end of it and save all th needed //
//      information of the component as mentioned above, in a singly linked list node. Every new //
//      component that is read, creates a new node that is connected with the next pointer of a  //
//      previous insertion. The head of the list is a global variable and always points to the   //
//      first insertion of a component to the list.                                              //
//   5) The final node points to NULL, indicating the last component of the netlist file.        //
//                                                                                               //
//   This netlist only uses one token (" ") to implement the lexical analysis.                   //
///////////////////////////////////////////////////////////////////////////////////////////////////

//   Useful "C" Libraries                      //

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <stdint.h>
#include <complex.h>
#include <math.h>
#include <gsl/gsl_complex.h>

//  NETLIST_NUM defined as 1, is a macro indicating that argv[1] is the user's correct input and  //
//  is used as the user's netlist file input. In order to read the file from terminal,            //
//  NETLIST_ARGC must exclusively be equal to 2, otherwise the program stops and an error message //
//  is returned.                                                                                  //

#define NETLIST_NUM 1
#define NETLIST_ARGC 2

#define HASHTABLE_DEPTH 15

//   structs          //

// Node structure to represent each element in the hash table
typedef struct hashtable_t{
    char* old_name[HASHTABLE_DEPTH];
    int new_name_value[HASHTABLE_DEPTH];
}hashtable;

hashtable *ht;

//   main component struct for storing all required information                                    //
struct pos_node_t{
    char *pos_node_name;
    int new_pos_node_name;
    double pos_node_value;
};

struct neg_node_t{
    char *neg_node_name;
    int new_neg_node_name;
    double neg_node_value;
};

typedef struct pos_node_t positive_node;

typedef struct neg_node_t negative_node;

struct exp_t{
    double exp_i1;
    double exp_i2;
    double exp_td1;
    double exp_tc1;
    double exp_td2;
    double exp_tc2;
};

typedef struct exp_t exp_spec;

struct sin_t{
    double sin_i1;
    double sin_ia;
    double sin_fr;
    double sin_td;
    double sin_df;
    double sin_ph;
};

typedef struct sin_t sin_spec;

struct pulse_t{
    double pulse_i1;
    double pulse_i2;
    double pulse_td;
    double pulse_tr;
    double pulse_tf;
    double pulse_pw;
    double pulse_per;
};

typedef struct pulse_t pulse_spec;

struct pwl_t{
    double pwl_start;
    double pwl_end;
};

typedef struct pwl_t* pwl_spec;

struct complex_num_t{
    double real_part;
    double imag_part;
    double complex complex_num;
};

typedef struct complex_num_t complex_num;

struct component_t{
    char comp_type;
    char *comp_name;
    positive_node pos_node;
    negative_node neg_node;
    double value;
    struct component_t *nxt_comp;   // pointer to the next singly linked list node
    int vsource_pos;
    int exp_flag;
    int pulse_flag;
    int sin_flag;
    int pwl_flag;
    int dc_flag;
    exp_spec exp_field;
    sin_spec sin_field;
    pulse_spec pulse_field;
    pwl_spec pwl_field;
    int comp_pwl_counter;
    double AC_mag;
    double AC_phase;
    int AC_component;
    //complex_num AC_complex_number;
    gsl_complex AC_complex_number;
};


//   typedefs         //

typedef struct component_t component;

typedef struct component_t* ptr_comp;

//   global variables //

ptr_comp list_head; // global variable used as a head of the singly linked list, for always pointing to the first node of the list.

int G2_components;

int G1_components;

unsigned long component_counter;    // counts the number of circuit components of each netlist (R, L, C, I, V, r, l, c, v, i).

char **oldname_node_table;

int curr_pos_oldname_table; // number of nodes except node0, which is the ground

int ht_size;
int ht_depth;

int spd_flag;

struct dcsweep_t{
    int dcsweep_flag;
    char *dcsweep_component;
    double dcsweep_startval;
    double dcsweep_endval;
    double dcsweep_step;
};

typedef struct dcsweep_t dcsweep;

dcsweep DCsweep;

int plot_flag;

char *plot_node;

char *print_node;

int iter_flag;

double ITOL;

int sparse_flag;

int tran_flag;

double tran_time_step;

double tran_fin_time;

int trapezoidal_flag;
int backwardeuler_flag;

int tran_m_parameter;

int AC_flag;

double AC_start_freq;
double AC_end_freq;

int AC_LIN_flag;
int AC_LOG_flag;

int AC_points;

//   function declarations                     //

//   initialization of the global variables performed before the start of the parsing.           //
void init_global_vars();

//   initialization of each field of every new component before its insertion to the singly      //
//   linked list.                                                                                //
void init_components(ptr_comp);

//   main parsing function, which reads each line of the netlist file using the space-token and  //
//   stores component's information to the corresponding struct's field                          //
void add_components(FILE *, const char *);

//   printing of all the components and their information, used for debugging reasons            //
void print_components();
void print_updated_components();

//   destroying the singly linked list after all the necessary frees                             //
void destroy_components();

unsigned long hash(char *);

void ht_dimensions();

void ht_creation();

void add_ht_oldname();

void add_oldname_pos();

void print_ht();

void free_ht();

char *return_new_filename(char *);

char *return_solution_filename(char *filename);

void choose_method(char *);