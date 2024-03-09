#include "SPARSE_DC_system_creation.h"
#include "netlist_parser.h"
#include "SPARSE_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_direct_methods.h"

void init_global_sparse_vars(){
    non_zeros = 0;

    sparse_rowsA = curr_pos_oldname_table + G2_components;  //sparse_A->m = sparse_rowsA;

    sparse_columnsA = curr_pos_oldname_table + G2_components;   //sparse_A->n = sparse_columnsA;

    sparse_k_variable = 0;

    sparse_counter = 0;

    sparse_b = (double *)calloc(sparse_rowsA, sizeof(double));

    if(sparse_b == NULL){
        printf("\nMemory allocation for sparse matrix B failed.\n");

        exit(8);
    }

    non_zeros_tran_C = 0;

    tran_sparse_counter = 0;
}

void count_non_zeros(){
    init_global_sparse_vars();

    ptr_comp list_ptr;

    list_ptr = list_head;

    while(list_ptr != NULL){
        switch(list_ptr->comp_type){
        case 'R':
        case 'r':
            sparse_resistance(list_ptr->pos_node.new_pos_node_name, list_ptr->neg_node.new_neg_node_name);
            break;

        case 'L':
        case 'l':
            sparse_voltagesource(list_ptr->pos_node.new_pos_node_name, list_ptr->neg_node.new_neg_node_name);
            TRAN_sparse_inductance();
            break;
        
        case 'V':
        case 'v':
            sparse_voltagesource(list_ptr->pos_node.new_pos_node_name, list_ptr->neg_node.new_neg_node_name);
            break;

        case 'C':
        case 'c':
            TRAN_sparse_capacitance(list_ptr->pos_node.new_pos_node_name, list_ptr->neg_node.new_neg_node_name);
        
        default:
            break;
        }


        list_ptr = list_ptr -> nxt_comp;
    }
}

void sparse_resistance(int posnode, int negnode){
    if((posnode == 0) && (negnode == 0)){
        return;
    }

    if(((posnode != 0) && (negnode == 0)) || ((posnode == 0) && (negnode != 0))){
        non_zeros = non_zeros + 1;

        return;
    }

    if((posnode != 0) && (negnode != 0)){
        non_zeros = non_zeros + 4;

        return;
    }

}

void sparse_voltagesource(int posnode, int negnode){
    if((posnode == 0) && (negnode == 0)){
        return;
    }

    if(((posnode != 0) && (negnode == 0)) || ((posnode == 0) && (negnode != 0))){
        non_zeros = non_zeros + 2;

        return;
    }

    if((posnode != 0) && (negnode != 0)){
        non_zeros = non_zeros + 4;

        return;
    }
}

void TRAN_sparse_capacitance(int posnode, int negnode){
    if((posnode == 0) && (negnode == 0)){
        return;
    }

    if(((posnode != 0) && (negnode == 0)) || ((posnode == 0) && (negnode != 0))){
        non_zeros_tran_C = non_zeros_tran_C + 1;

        return;
    }

    if((posnode != 0) && (negnode != 0)){
        non_zeros_tran_C = non_zeros_tran_C + 4;

        return;
    }    
}

void TRAN_sparse_inductance(){
    non_zeros_tran_C = non_zeros_tran_C + 1;
}

void triplet_format_method(){
    count_non_zeros();

    ptr_comp triplet_format_ptr;

    triplet_format_ptr = list_head;

    //sparse_A->nz = non_zeros;

    sparse_A = cs_spalloc(sparse_rowsA, sparse_columnsA, non_zeros, 1, 1);
    sparse_A->nz = non_zeros;

    tran_sparse_C = cs_spalloc(sparse_rowsA, sparse_columnsA, non_zeros_tran_C, 1, 1);
    tran_sparse_C->nz = non_zeros_tran_C;

    //tran_SPARSE_tilde_G = cs_spalloc(sparse_rowsA, sparse_columnsA, non_zeros, 1, 1);
    //tran_SPARSE_tilde_G -> nz = non_zeros_tran_C;

    while(triplet_format_ptr != NULL){
        switch(triplet_format_ptr -> comp_type){
            case 'R':
            case 'r':
                triplet_format_resistance(triplet_format_ptr->pos_node.new_pos_node_name, 
                                                triplet_format_ptr->neg_node.new_neg_node_name, triplet_format_ptr->value);
                break;
            
            case 'V':
            case 'v':
                triplet_format_voltagesource(triplet_format_ptr->pos_node.new_pos_node_name, 
                                                triplet_format_ptr->neg_node.new_neg_node_name, triplet_format_ptr->value);
                break;

            case 'L':
            case 'l':
                triplet_format_inductance(triplet_format_ptr->pos_node.new_pos_node_name, 
                                                triplet_format_ptr->neg_node.new_neg_node_name, triplet_format_ptr->value);
                break;
        
            case 'I':
            case 'i':
                triplet_format_currentsource(triplet_format_ptr->pos_node.new_pos_node_name, 
                                                triplet_format_ptr->neg_node.new_neg_node_name, triplet_format_ptr->value);
                break;

            case 'C':
            case 'c':
                sparse_triplet_format_capacitance(triplet_format_ptr->pos_node.new_pos_node_name, 
                                                triplet_format_ptr->neg_node.new_neg_node_name, triplet_format_ptr->value);
                break;

            default:
                break;
        }


        triplet_format_ptr = triplet_format_ptr -> nxt_comp;
    }

    sparse_C = cs_compress(sparse_A);
    cs_dupl(sparse_C);

    tran_sparse_tilde_C = cs_compress(tran_sparse_C);
    cs_dupl(tran_sparse_tilde_C);

    tran_SPARSE_tilde_G = cs_compress(sparse_A);
    cs_dupl(tran_SPARSE_tilde_G);

    //cs_print(tran_SPARSE_tilde_G, "tilde_G", 1);
    //cs_print(sparse_C, "sparse_C", 1);
    //cs_print(tran_sparse_tilde_C, "tilde_C", 1);
}

void triplet_format_resistance(int posnode, int negnode, double value){
    if((posnode == 0) && (negnode == 0)){
        return;
    }

    if(((posnode != 0) && (negnode == 0))){
        sparse_A->i[sparse_counter] = posnode - 1;
        sparse_A->p[sparse_counter] = posnode - 1;
        sparse_A->x[sparse_counter] = (1 / value);

        sparse_counter = sparse_counter + 1;

        return;
    }

    if((posnode == 0) && (negnode != 0)){
        sparse_A->i[sparse_counter] = negnode - 1;
        sparse_A->p[sparse_counter] = negnode - 1;
        sparse_A->x[sparse_counter] = (1 / value);

        sparse_counter = sparse_counter + 1;

        return;        
    }

    if((posnode != 0) && (negnode != 0)){
        sparse_A->i[sparse_counter] = posnode - 1;
        sparse_A->p[sparse_counter] = posnode - 1;

        sparse_A->x[sparse_counter] = (1 / value);

        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = posnode - 1;
        sparse_A->p[sparse_counter] = negnode - 1;

        sparse_A->x[sparse_counter] = - (1 / value);

        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = negnode - 1;
        sparse_A->p[sparse_counter] = posnode - 1;

        sparse_A->x[sparse_counter] = - (1 / value);

        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = negnode - 1;
        sparse_A->p[sparse_counter] = negnode - 1;

        sparse_A->x[sparse_counter] = (1 / value);

        sparse_counter = sparse_counter + 1;

        return;
    }
}

void triplet_format_voltagesource(int posnode, int negnode, double value){
    if((posnode == 0) && (negnode == 0)){
        return;
    }

    if((posnode != 0) && (negnode == 0)){
        sparse_A->i[sparse_counter] = posnode - 1;
        sparse_A->p[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->x[sparse_counter] = 1;
        
        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->p[sparse_counter] = posnode - 1;
        sparse_A->x[sparse_counter] = 1;

        sparse_b[curr_pos_oldname_table + sparse_k_variable] = sparse_b[curr_pos_oldname_table + sparse_k_variable] + value;

        sparse_counter = sparse_counter + 1;
        sparse_k_variable = sparse_k_variable + 1;

        return;
    }

    if((posnode == 0) && (negnode != 0)){
        sparse_A->i[sparse_counter] = negnode - 1;
        sparse_A->p[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->x[sparse_counter] = -1;
        
        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->p[sparse_counter] = negnode - 1;
        sparse_A->x[sparse_counter] = -1;

        sparse_b[curr_pos_oldname_table + sparse_k_variable] = sparse_b[curr_pos_oldname_table + sparse_k_variable] + value;

        sparse_counter = sparse_counter + 1;
        sparse_k_variable = sparse_k_variable + 1;

        return;
    }

    if((posnode != 0) && (negnode != 0)){
        sparse_A->i[sparse_counter] = posnode - 1;
        sparse_A->p[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->x[sparse_counter] = 1;
        
        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->p[sparse_counter] = posnode - 1;
        sparse_A->x[sparse_counter] = 1;

        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = negnode - 1;
        sparse_A->p[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->x[sparse_counter] = -1;
        
        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->p[sparse_counter] = negnode - 1;
        sparse_A->x[sparse_counter] = -1;

        sparse_b[curr_pos_oldname_table + sparse_k_variable] = sparse_b[curr_pos_oldname_table + sparse_k_variable] + value;

        sparse_counter = sparse_counter + 1;
        sparse_k_variable = sparse_k_variable + 1;

        return; 
    }
}

void triplet_format_currentsource(int posnode, int negnode, double value){
    if((posnode != 0) && (negnode == 0)){
        sparse_b[posnode - 1] = sparse_b[posnode - 1] - value;

        return;
    }

    if((posnode == 0) && (negnode != 0)){
        sparse_b[negnode - 1] = sparse_b[negnode - 1] + value;

        return;
    }

    if((posnode != 0) && (negnode != 0)){
        sparse_b[posnode - 1] = sparse_b[posnode - 1] - value;
        sparse_b[negnode - 1] = sparse_b[negnode - 1] + value;

        return;
    }
}

void triplet_format_inductance(int posnode, int negnode, double value){
    tran_sparse_C->i[tran_sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
    tran_sparse_C->p[tran_sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
    tran_sparse_C->x[tran_sparse_counter] = -value;

    tran_sparse_counter = tran_sparse_counter + 1;

    if((posnode == 0) && (negnode == 0)){
        return;
    }

    if((posnode != 0) && (negnode == 0)){
        sparse_A->i[sparse_counter] = posnode - 1;
        sparse_A->p[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->x[sparse_counter] = 1;
        
        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->p[sparse_counter] = posnode - 1;
        sparse_A->x[sparse_counter] = 1;

        sparse_b[curr_pos_oldname_table + sparse_k_variable] = sparse_b[curr_pos_oldname_table + sparse_k_variable] + 0;

        sparse_counter = sparse_counter + 1;
        sparse_k_variable = sparse_k_variable + 1;

        return;
    }

    if((posnode == 0) && (negnode != 0)){
        sparse_A->i[sparse_counter] = negnode - 1;
        sparse_A->p[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->x[sparse_counter] = -1;
        
        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->p[sparse_counter] = negnode - 1;
        sparse_A->x[sparse_counter] = -1;

        sparse_b[curr_pos_oldname_table + sparse_k_variable] = sparse_b[curr_pos_oldname_table + sparse_k_variable] + 0;

        sparse_counter = sparse_counter + 1;
        sparse_k_variable = sparse_k_variable + 1;

        return;
    }

    if((posnode != 0) && (negnode != 0)){
        sparse_A->i[sparse_counter] = posnode - 1;
        sparse_A->p[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->x[sparse_counter] = 1;
        
        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->p[sparse_counter] = posnode - 1;
        sparse_A->x[sparse_counter] = 1;

        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = negnode - 1;
        sparse_A->p[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->x[sparse_counter] = -1;
        
        sparse_counter = sparse_counter + 1;

        sparse_A->i[sparse_counter] = curr_pos_oldname_table + sparse_k_variable;
        sparse_A->p[sparse_counter] = negnode - 1;
        sparse_A->x[sparse_counter] = -1;

        sparse_b[curr_pos_oldname_table + sparse_k_variable] = sparse_b[curr_pos_oldname_table + sparse_k_variable] + 0;

        sparse_counter = sparse_counter + 1;
        sparse_k_variable = sparse_k_variable + 1;

        return; 
    }  
}

void sparse_triplet_format_capacitance(int posnode, int negnode, double value){
    if((posnode == 0) && (negnode == 0)){
        return;
    }

    if(((posnode != 0) && (negnode == 0))){
        tran_sparse_C->i[tran_sparse_counter] = posnode - 1;
        tran_sparse_C->p[tran_sparse_counter] = posnode - 1;
        tran_sparse_C->x[tran_sparse_counter] = value;

        tran_sparse_counter = tran_sparse_counter + 1;

        return;
    }

    if((posnode == 0) && (negnode != 0)){
        tran_sparse_C->i[tran_sparse_counter] = negnode - 1;
        tran_sparse_C->p[tran_sparse_counter] = negnode - 1;
        tran_sparse_C->x[tran_sparse_counter] = (1 / value);

        tran_sparse_counter = tran_sparse_counter + 1;

        return;        
    }

    if((posnode != 0) && (negnode != 0)){
        tran_sparse_C->i[tran_sparse_counter] = posnode - 1;
        tran_sparse_C->p[tran_sparse_counter] = posnode - 1;

        tran_sparse_C->x[tran_sparse_counter] = (1 / value);

        tran_sparse_counter = tran_sparse_counter + 1;

        tran_sparse_C->i[tran_sparse_counter] = posnode - 1;
        tran_sparse_C->p[tran_sparse_counter] = negnode - 1;

        tran_sparse_C->x[tran_sparse_counter] = - (1 / value);

        tran_sparse_counter = tran_sparse_counter + 1;

        tran_sparse_C->i[tran_sparse_counter] = negnode - 1;
        tran_sparse_C->p[tran_sparse_counter] = posnode - 1;

        tran_sparse_C->x[tran_sparse_counter] = - (1 / value);

        tran_sparse_counter = tran_sparse_counter + 1;

        tran_sparse_C->i[tran_sparse_counter] = negnode - 1;
        tran_sparse_C->p[tran_sparse_counter] = negnode - 1;

        tran_sparse_C->x[tran_sparse_counter] = (1 / value);

        tran_sparse_counter = tran_sparse_counter + 1;

        return;
    }   
}

void destroy_sparse_matrices(){
    cs_spfree(tran_SPARSE_tilde_G);
    cs_spfree(tran_sparse_tilde_C);
    cs_spfree(tran_sparse_C);
    
    cs_spfree(sparse_A);
    cs_spfree(sparse_C);

    free(sparse_b);

    free(sparse_b_plot);
}