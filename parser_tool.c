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
//   Main function of the parsing tool consisting of the following steps:                        //
//   1) Open and read Netlist File (.txt form) from terminal as the second argument after the    //
//      running the executable. If .txt file is not opened or read successfully, print specific  //
//      error message.                                                                           //
//                                                                                               //
//   2) Store all the essential information of each circuit component: i) Component's type.      //
//                                                                     ii) Component's name.     //
//                                                                     iii) Component's positive //
//                                                                          and negative nodes   //
//                                                                          names and values.    // 
//                                                                     iv) Component's value.    //
//                                                                                               //
//   3) Print all stored components and their information for debugging reasons.                 //
//                                                                                               //
//   4) Clear the singly linked list that was created for connecting each component of the given //
//      netlist file, along with all the allocated space for its field of the component's        //
//      struct.                                                                                  //
//                                                                                               //
//   5) Close the opened .txt Netlist file and terminate parsing tool.                           //
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "netlist_parser.h"
#include "MNA_DC_system_creation.h"
#include "MNA_DC_system_direct_methods.h"
#include "MNA_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_creation.h"

int main(int argc, char *argv[]){
    FILE *fptrnetlist;
    char *filename;
    char *newfilename;
    char *solutionfilename;

    if(argc == NETLIST_ARGC){
        fptrnetlist = fopen(argv[NETLIST_NUM], "r");
        filename = strdup(argv[NETLIST_NUM]);
        
        newfilename = return_new_filename(filename);
        solutionfilename = return_solution_filename(filename);
    }
    else{
        printf("\nFormat should be: ./parser_tool MY_NETLIST_FILE.txt\n");
    }

    if(fptrnetlist == NULL){
        printf("\nNetlist file could not be opened.\n");
        
        exit(0);
    }

    add_components(fptrnetlist, " ");

    add_ht_oldname();

    add_oldname_pos();

    //print_ht();
    //print_components();
    //print_updated_components();

    //global_arrays_allocation();

    //print_global_arrays(array_A, rowsA, columnsA);
    //print_global_vectors(array_b, rowsb);
    //print_global_vectors(array_x, rowsx);

/***************************************************/
    choose_method(solutionfilename);
/***************************************************/

    destroy_components();

    destroy_sparse_matrices();

    fclose(fptrnetlist);
    
    free(filename);
    //free(solutionfilename);

    return(0);
}
