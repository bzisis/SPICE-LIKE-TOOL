#include "MNA_DC_system_creation.h"
#include "netlist_parser.h"

void init_global_arrays(){
    array_A = NULL;

    array_b = NULL;

    array_x = NULL;

    rowsA = curr_pos_oldname_table + G2_components; // curr_pos_oldname_table + G2_components = n = (n - 1) + m2
    columnsA = curr_pos_oldname_table +  G2_components;

    rowsb = curr_pos_oldname_table + G2_components;
    //columnsb = 1;


    rowsx = curr_pos_oldname_table + G2_components;
    //columnsx = 1;

    k_variable = 0;

    array_C = NULL;
}

void global_arrays_allocation(){
    int i;

    init_global_arrays();

    array_A = (double **)calloc(rowsA, sizeof(double *));
    if(!array_A){
        printf("\nMemory allocation for global arrays failed.\n");

        exit(4);
    }

    for(i = 0; i < rowsA; i++){
        array_A[i] = (double *)calloc(columnsA, sizeof(double));
        if(!array_A[i]){
            printf("\nMemory allocation for global arrays failed.\n");

            exit(5);
        }
    }

    array_b = (double *)calloc(rowsb, sizeof(double));
    if(!array_b){
        printf("\nMemory allocation for global arrays failed.\n");

        exit(4);
    }

    array_x = (double *)calloc(rowsx, sizeof(double));
    if(!array_x){
        printf("\nMemory allocation for global arrays failed.\n");

        exit(4);
    }

    array_C = (double **)calloc(rowsA, sizeof(double *));
    if(!array_C){
        printf("\nMemory allocation for global arrays failed.\n");

        exit(4);
    }
    for(i = 0; i < rowsA; i++){
        array_C[i] = (double *)calloc(columnsA, sizeof(double));
        if(!array_C[i]){
            printf("\nMemory allocation for global arrays failed.\n");

            exit(5);
        }
    }

}

void print_global_arrays(double **array, int rows, int columns){
    int i, j;

    printf("\n");
    if(array == array_A){
        printf("Array A is:\n");
        
        for(i = 0; i < rows; i++){
            for(j = 0; j < columns; j++){
                printf("%lf\t", array[i][j]); 
            }
            printf("\n");
        }
    }
}

void print_global_vectors(double *array, int rows){
    int i;

    printf("\n");
    if(array == array_b){
        printf("Array b is:\n");
    }

    else if(array == array_x){
        printf("Array x is:\n");
    }

    for(i = 0; i < rows; i++){
        printf("%lf\t", array[i]); 
        printf("\n");
    }
}

void MNA_resistance(int MNA_pos_node, int MNA_neg_node, double MNA_value){
    int row_position, column_position;

    double MNA_g;

    MNA_g = 1 / (MNA_value);

    row_position = MNA_pos_node - 1;
    column_position = MNA_neg_node - 1;

    if(MNA_pos_node == 0){
        array_A[column_position][column_position] = array_A[column_position][column_position] + MNA_g;
        
        return;
    }

    if(MNA_neg_node == 0){
        array_A[row_position][row_position] = array_A[row_position][row_position] + MNA_g;

        return;
    }
    // if((MNA_pos_node && MNA_neg_node ) != 0)
    array_A[row_position][row_position] = array_A[row_position][row_position] + MNA_g;

    array_A[row_position][column_position] = array_A[row_position][column_position] - MNA_g;

    array_A[column_position][row_position] = array_A[column_position][row_position] - MNA_g;

    array_A[column_position][column_position] = array_A[column_position][column_position] + MNA_g;
}

void MNA_currentsource(int MNA_pos_node, int MNA_neg_node, double MNA_value){
    int pos_position, neg_position;

    pos_position = MNA_pos_node - 1;
    neg_position = MNA_neg_node - 1;

    if(MNA_pos_node == 0){
        array_b[neg_position] = array_b[neg_position] + MNA_value;

        return;
    }

    if(MNA_neg_node == 0){
        array_b[pos_position] = array_b[pos_position] - MNA_value;

        return;
    }

    array_b[pos_position] = array_b[pos_position] - MNA_value;

    array_b[neg_position] = array_b[neg_position] + MNA_value;

}

void MNA_voltagesource(int MNA_pos_node, int MNA_neg_node, double MNA_value){
    int pos_position, neg_position;

    int k_expression;

    pos_position = MNA_pos_node - 1;
    neg_position = MNA_neg_node - 1;

    k_expression = curr_pos_oldname_table + k_variable;

    if(MNA_pos_node == 0){
        array_A[k_expression][neg_position] = array_A[k_expression][neg_position] - 1;

        array_A[neg_position][k_expression] = array_A[neg_position][k_expression] - 1;

        array_b[k_expression] = array_b[k_expression] + MNA_value;

        k_variable++;

        return;
    }
    
    if(MNA_neg_node == 0){
        array_A[pos_position][k_expression] = array_A[pos_position][k_expression] + 1;

        array_A[k_expression][pos_position] = array_A[k_expression][pos_position] + 1;

        array_b[k_expression] = array_b[k_expression] + MNA_value;

        k_variable++;

        return;
    }

    array_A[pos_position][k_expression] = array_A[pos_position][k_expression] + 1;
    array_A[k_expression][pos_position] = array_A[k_expression][pos_position] + 1;

    array_A[k_expression][neg_position] = array_A[k_expression][neg_position] - 1;
    array_A[neg_position][k_expression] = array_A[neg_position][k_expression] - 1;

    array_b[k_expression] = array_b[k_expression] + MNA_value;

    k_variable++;
}
//  MNA_inductance not needed at this point because we manipulate inductance as a voltage source //
//  with zero MNA_value.                                                                         //

void MNA_inductance(int MNA_pos_node, int MNA_neg_node, double MNA_value){
    int pos_position, neg_position;

    int k_expression;

    pos_position = MNA_pos_node - 1;
    neg_position = MNA_neg_node - 1;

    k_expression = curr_pos_oldname_table + k_variable;

    if(MNA_pos_node == 0){
        array_A[k_expression][neg_position] = array_A[k_expression][neg_position] - 1;

        array_A[neg_position][k_expression] = array_A[neg_position][k_expression] - 1;

        array_b[k_expression] = array_b[k_expression] + 0;

        array_C[k_expression][k_expression] = -MNA_value;

        k_variable++;

        return;
    }
    
    if(MNA_neg_node == 0){
        array_A[pos_position][k_expression] = array_A[pos_position][k_expression] + 1;

        array_A[k_expression][pos_position] = array_A[k_expression][pos_position] + 1;

        array_b[k_expression] = array_b[k_expression] + 0;

        array_C[k_expression][k_expression] = -MNA_value;

        k_variable++;

        return;
    }

    array_A[pos_position][k_expression] = array_A[pos_position][k_expression] + 1;
    array_A[k_expression][pos_position] = array_A[k_expression][pos_position] + 1;

    array_A[k_expression][neg_position] = array_A[k_expression][neg_position] - 1;
    array_A[neg_position][k_expression] = array_A[neg_position][k_expression] - 1;

    array_b[k_expression] = array_b[k_expression] + 0;

    array_C[k_expression][k_expression] = -MNA_value;

    k_variable++;
}

void MNA_algorithm(){
    ptr_comp MNA_ptr;

    MNA_ptr = list_head;

    while(MNA_ptr != NULL){
        switch(MNA_ptr -> comp_type){
            case 'R':
            case 'r':
                MNA_resistance(MNA_ptr -> pos_node.new_pos_node_name, MNA_ptr -> neg_node.new_neg_node_name, MNA_ptr -> value);
                break;

            case 'I':
            case 'i':
                MNA_currentsource(MNA_ptr -> pos_node.new_pos_node_name, MNA_ptr -> neg_node.new_neg_node_name, MNA_ptr -> value);
                break;

            case 'V':
            case 'v':
                MNA_voltagesource(MNA_ptr -> pos_node.new_pos_node_name, MNA_ptr -> neg_node.new_neg_node_name, MNA_ptr -> value);
                break;
            
            case 'L':
            case 'l':
                //MNA_voltagesource(MNA_ptr -> pos_node.new_pos_node_name, MNA_ptr -> neg_node.new_neg_node_name, 0);
                MNA_inductance(MNA_ptr -> pos_node.new_pos_node_name, MNA_ptr -> neg_node.new_neg_node_name, MNA_ptr -> value);
                break;
            
            case 'C':
            case 'c':
                TRAN_capacitance(MNA_ptr -> pos_node.new_pos_node_name, MNA_ptr -> neg_node.new_neg_node_name, MNA_ptr -> value);
                break;
                
            default:
                break;
        }

        MNA_ptr = MNA_ptr -> nxt_comp;
    }
}

void TRAN_capacitance(int MNA_pos_node, int MNA_neg_node, double MNA_value){
    int row_position, column_position;

    row_position = MNA_pos_node - 1;
    column_position = MNA_neg_node - 1;

    if(MNA_pos_node == 0){
        array_C[column_position][column_position] = array_C[column_position][column_position] + MNA_value;
        
        return;
    }

    if(MNA_neg_node == 0){
        array_C[row_position][row_position] = array_C[row_position][row_position] + MNA_value;

        return;
    }
    // if((MNA_pos_node && MNA_neg_node ) != 0)
    array_C[row_position][row_position] = array_C[row_position][row_position] + MNA_value;

    array_C[row_position][column_position] = array_C[row_position][column_position] - MNA_value;

    array_C[column_position][row_position] = array_C[column_position][row_position] - MNA_value;

    array_C[column_position][column_position] = array_C[column_position][column_position] + MNA_value;

}

void destroy_global_arrays(double **array, int rows){
    int i;

    for(i = 0; i < rows; i++){
        free(array[i]);
    }
    free(array);
}

void destroy_global_vectors(double *array){
    free(array);
}