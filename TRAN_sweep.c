#include "TRAN_sweep.h"
#include "TRAN_DC_system_creation.h"
#include "SPARSE_DC_system_creation.h"
#include "MNA_DC_system_creation.h"
#include "netlist_parser.h"
#include "MNA_DC_system_direct_methods.h"
#include "MNA_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_direct_methods.h"


double PWL_lookup(double, double, int, ptr_comp);

void choose_tran_method(){
    if((trapezoidal_flag == 1) || ((trapezoidal_flag == 0) && (backwardeuler_flag == 0) && (tran_flag == 1))){
        if(sparse_flag == 0){
            trapezoidal_method();
        }

        else{
            SPARSE_trapezoidal_method();
        }
    }

    if(backwardeuler_flag == 1){
        if(sparse_flag == 0){
            backwardeuler_method();
        }

        else{
            SPARSE_backwardeuler_method();
        }
    }
}

void trapezoidal_method(){
    init_TR_hands();

    int k, i;
    k_pulse_variable = 0;

    double tk;
    double h;

    int first_iter = 1;

    double lookup_result = 0;

    double tran_node_result = 0;

    char *tran_plot_ptr = NULL;
    int tran_plot_node_name = 0;

    double tran_plot_x = 0;
    double tran_plot_y = 0;

    FILE *tran_fileptr = NULL;
    FILE *tran_gnuplotpipe = NULL;
    char *tran_gnuplotcommands = {"plot 'tranexample.tmp' with linespoints"};

    tran_fileptr = fopen("tranexample.tmp", "w");
    tran_gnuplotpipe = popen("gnuplot -persistent", "w");

    ptr_comp TR_ptr = NULL;

    h = (2 / tran_time_step);

    // set left_hand
    gsl_matrix_memcpy(alt_tran_MNA_G, tran_MNA_G);

    gsl_matrix_scale(tran_MNA_C, h);    // (1 / h) * Ctilde

    gsl_matrix_add(tran_MNA_G, tran_MNA_C); // (Gtilde + ((1 / h) * Ctilde))
    gsl_matrix_memcpy(left_hand_TR, tran_MNA_G);

    gsl_matrix_sub(alt_tran_MNA_G, tran_MNA_C); // (Gtilde - ((2/h) * Ctilde))
    gsl_matrix_memcpy(right_hand_TR_sub, alt_tran_MNA_G);

    if(plot_node != NULL){
        tran_plot_ptr = strdup(plot_node);
    }
    
    if(print_node != NULL){
        tran_plot_ptr = strdup(print_node);
    }

    if(tran_plot_ptr != NULL){
        tran_plot_node_name = return_newpos_node_name(tran_plot_ptr);
        
        free(tran_plot_ptr);
    }

    for(k = 1; k <= tran_m_parameter; k++){
        k_tran_variable = 0;

        tk = k * tran_time_step;

        TR_ptr = list_head;

        for(i = 0; i < tran_MNA_rows; i++){
            gsl_vector_set(tran_MNA_e, i, 0);
        }

        while(TR_ptr != NULL){
            if((TR_ptr->comp_type == 'V') || (TR_ptr->comp_type == 'v')){
                if(TR_ptr -> dc_flag == 1){
                    TRAN_sweep_voltagesource(TR_ptr -> value);
                }

                if(TR_ptr -> exp_flag == 1){
                    lookup_result = EXP_lookup(tk, TR_ptr->exp_field.exp_i1, TR_ptr->exp_field.exp_i2, TR_ptr->exp_field.exp_td1, 
                                                TR_ptr->exp_field.exp_td2, TR_ptr->exp_field.exp_tc1, TR_ptr->exp_field.exp_tc2, 
                                                tran_fin_time);

                    TRAN_sweep_voltagesource(lookup_result);
                }

                if(TR_ptr -> pwl_flag == 1){
                    lookup_result = PWL_lookup(tk, tran_fin_time, TR_ptr->comp_pwl_counter, TR_ptr);

                    TRAN_sweep_voltagesource(lookup_result);
                }

                if(TR_ptr -> sin_flag == 1){
                    lookup_result = SIN_lookup(tk, TR_ptr->sin_field.sin_i1, TR_ptr->sin_field.sin_ia, TR_ptr->sin_field.sin_fr,
                                                TR_ptr->sin_field.sin_td, TR_ptr->sin_field.sin_df, TR_ptr->sin_field.sin_ph,
                                                tran_fin_time);
                    
                    TRAN_sweep_voltagesource(lookup_result);
                }

                if(TR_ptr -> pulse_flag == 1){
                    lookup_result = PULSE_lookup(tk, TR_ptr->pulse_field.pulse_i1, TR_ptr->pulse_field.pulse_i2,
                                                    TR_ptr->pulse_field.pulse_td, TR_ptr->pulse_field.pulse_tr,
                                                    TR_ptr->pulse_field.pulse_tf, TR_ptr->pulse_field.pulse_pw, 
                                                    TR_ptr->pulse_field.pulse_per);

                    TRAN_sweep_voltagesource(lookup_result);
                }
            }

            if((TR_ptr->comp_type == 'I') || (TR_ptr->comp_type == 'i')){
                if(TR_ptr -> dc_flag == 1){
                    TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, TR_ptr->value);
                }

                 if(TR_ptr -> exp_flag == 1){
                    lookup_result = EXP_lookup(tk, TR_ptr->exp_field.exp_i1, TR_ptr->exp_field.exp_i2, TR_ptr->exp_field.exp_td1, 
                                                TR_ptr->exp_field.exp_td2, TR_ptr->exp_field.exp_tc1, TR_ptr->exp_field.exp_tc2, 
                                                tran_fin_time);

                    TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(TR_ptr -> pwl_flag == 1){
                    lookup_result = PWL_lookup(tk, tran_fin_time, TR_ptr->comp_pwl_counter, TR_ptr);

                    TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(TR_ptr -> sin_flag == 1){
                    lookup_result = SIN_lookup(tk, TR_ptr->sin_field.sin_i1, TR_ptr->sin_field.sin_ia, TR_ptr->sin_field.sin_fr,
                                                TR_ptr->sin_field.sin_td, TR_ptr->sin_field.sin_df, TR_ptr->sin_field.sin_ph,
                                                tran_fin_time);
                    
                    TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(TR_ptr -> pulse_flag == 1){
                    lookup_result = PULSE_lookup(tk, TR_ptr->pulse_field.pulse_i1, TR_ptr->pulse_field.pulse_i2,
                                                    TR_ptr->pulse_field.pulse_td, TR_ptr->pulse_field.pulse_tr,
                                                    TR_ptr->pulse_field.pulse_tf, TR_ptr->pulse_field.pulse_pw, 
                                                    TR_ptr->pulse_field.pulse_per);

                    TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, lookup_result);
                }
            }

            if((TR_ptr->comp_type == 'L') || (TR_ptr->comp_type == 'l')){
                TRAN_sweep_voltagesource(0);
            }

            TR_ptr = TR_ptr -> nxt_comp;
        }

        if((spd_flag == 0) && (iter_flag == 0)){
            // LU
            TR_tran_LU_sweep(first_iter, k);

            first_iter = 0;
        }

        if((spd_flag == 1) && (iter_flag == 0)){
            // Cholesky
            TR_tran_Cholesky_sweep(first_iter, k);

            first_iter = 0;
        }

        if((spd_flag == 0) && (iter_flag == 1)){
            // BiCG
            TR_tran_BiCG_sweep(first_iter, k);

            first_iter = 0;
        }

        if((spd_flag == 1) && (iter_flag == 1)){
            // CG
            TR_tran_CG_sweep(first_iter, k);

            first_iter = 0;
        }

        if((tran_plot_node_name != -1) && (tran_plot_node_name != 0)){
            tran_node_result = gsl_vector_get(TR_x_curr, (tran_plot_node_name - 1));
        }

        tran_plot_x = tk;

        tran_plot_y = tran_node_result;

        fprintf(tran_fileptr, "%f %f\n", tran_plot_x, tran_plot_y);
    }

    fprintf(tran_gnuplotpipe, "%s\n", tran_gnuplotcommands);

    fprintf(tran_gnuplotpipe, "exit\n");

    fclose(tran_fileptr);
    pclose(tran_gnuplotpipe);
}

void TR_tran_LU_sweep(int iter, int k_value){
    int i, s;

    // first iteration (k = 1) -> x(0) , e(0) = 0
    // left_hand_TR * TR_x_curr = tran_MNA_e + (tran_MNA_G - ((2 / h) * tran_MNA_C)) * BE_x_prev)
    if(iter == 1){
        gsl_blas_dgemv(CblasNoTrans, -1.0, right_hand_TR_sub, tran_MNA_X0, 0.0, right_hand_TR);
        gsl_vector_add(right_hand_TR, tran_MNA_e);
    }
    
    // remaining iterations (k = 2, 3, ...) -> x(tk)
    if((iter == 0) && (k_value > 1)){
        gsl_blas_dgemv(CblasNoTrans, -1.0, right_hand_TR_sub, TR_x_prev, 0.0, right_hand_TR);
        gsl_vector_add(right_hand_TR, tran_MNA_e);
        gsl_vector_add(right_hand_TR, tran_MNA_e_prev);
    }

    for(i = 0; i < tran_MNA_rows; i++){
        gsl_vector_set(TR_x_curr, i, 0);
    }

    permutation_TR = gsl_permutation_alloc(tran_MNA_rows);

    alt_left_hand_TR = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    gsl_matrix_memcpy(alt_left_hand_TR, left_hand_TR);

    gsl_linalg_LU_decomp(alt_left_hand_TR, permutation_TR, &s);

    gsl_linalg_LU_solve(alt_left_hand_TR, permutation_TR, right_hand_TR, TR_x_curr);

    gsl_vector_memcpy(TR_x_prev, TR_x_curr);

    gsl_vector_memcpy(tran_MNA_e_prev, tran_MNA_e);

    gsl_permutation_free(permutation_TR);
    gsl_matrix_free(alt_left_hand_TR);
}

void TR_tran_Cholesky_sweep(int iter, int k_value){
    int i;

    // first iteration (k = 1) -> x(0) , e(0) = 0
    // left_hand_TR * TR_x_curr = tran_MNA_e + (tran_MNA_G - ((2 / h) * tran_MNA_C)) * BE_x_prev)
    if(iter == 1){
        gsl_blas_dgemv(CblasNoTrans, -1.0, right_hand_TR_sub, tran_MNA_X0, 0.0, right_hand_TR);
        gsl_vector_add(right_hand_TR, tran_MNA_e);
    }
    
    // remaining iterations (k = 2, 3, ...) -> x(tk)
    if((iter == 0) && (k_value > 1)){
        gsl_blas_dgemv(CblasNoTrans, -1.0, right_hand_TR_sub, TR_x_prev, 0.0, right_hand_TR);
        gsl_vector_add(right_hand_TR, tran_MNA_e);
        gsl_vector_add(right_hand_TR, tran_MNA_e_prev);
    }

    for(i = 0; i < tran_MNA_rows; i++){
        gsl_vector_set(TR_x_curr, i, 0);
    }

    alt_left_hand_TR = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    gsl_matrix_memcpy(alt_left_hand_TR, left_hand_TR);

    gsl_linalg_cholesky_decomp(alt_left_hand_TR);

    gsl_linalg_cholesky_solve(alt_left_hand_TR, right_hand_TR, TR_x_curr);

    gsl_vector_memcpy(TR_x_prev, TR_x_curr);

    gsl_vector_memcpy(tran_MNA_e_prev, tran_MNA_e);

    gsl_matrix_free(alt_left_hand_TR);
}

void TR_tran_BiCG_sweep(int iter, int k_value){
    int i, j;

    int tran_iter = 0;

    double tran_alpha = 0;
    double tran_beta = 0;

    double tran_norm_r = 0;
    double tran_norm_b = 0;

    double tran_rho = 0;
    double tran_rho1 = 0;

    double tran_omega = 0;

    tran_r_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_r_tilde_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_p_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_p_tilde_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_p_scale_beta = gsl_vector_alloc(tran_MNA_rows);
    tran_p_tilde_scale_beta = gsl_vector_alloc(tran_MNA_rows);

    tran_q_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_q_tilde_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_z_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_z_tilde_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_Ax_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_bminusAx_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_M_iterative = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    tran_M_inverse_iterative = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);

    tran_permut_iterative = gsl_permutation_alloc(tran_MNA_rows);

    for(i = 0; i < tran_MNA_rows; i++){
        for(j = 0; j < tran_MNA_columns; j++){
            gsl_matrix_set(tran_M_iterative, i, j, 0);

            gsl_matrix_set(tran_M_inverse_iterative, i, j, 0);
        }

        gsl_vector_set(tran_r_iterative, i, 0);
        gsl_vector_set(tran_r_tilde_iterative, i, 0);

        gsl_vector_set(tran_p_iterative, i , 0);
        gsl_vector_set(tran_p_tilde_iterative, i, 0);
        gsl_vector_set(tran_p_scale_beta, i, 0);
        gsl_vector_set(tran_p_tilde_scale_beta, i, 0);

        gsl_vector_set(tran_q_iterative, i, 0);
        gsl_vector_set(tran_q_tilde_iterative, i, 0);

        gsl_vector_set(tran_z_iterative, i, 0);
        gsl_vector_set(tran_z_tilde_iterative, i, 0);

        gsl_vector_set(tran_Ax_iterative, i, 0);
        gsl_vector_set(tran_bminusAx_iterative, i, 0);
    }

    if(iter == 1){
        // b = right_hand_BE
        gsl_blas_dgemv(CblasNoTrans, -1.0, right_hand_TR_sub, tran_MNA_X0, 0.0, right_hand_TR);
        gsl_vector_add(right_hand_TR, tran_MNA_e);
    }

    if((iter == 0) && (k_value > 1)){
        gsl_blas_dgemv(CblasNoTrans, -1.0, right_hand_TR_sub, TR_x_prev, 0.0, right_hand_TR);
        gsl_vector_add(right_hand_TR, tran_MNA_e);
        gsl_vector_add(right_hand_TR, tran_MNA_e_prev);
    }

    // initial guess x(0)
    for(i = 0; i < tran_MNA_rows; i++){
        gsl_vector_set(TR_x_curr, i, 0);
    }

    // A = alt_left_hand_BE
    alt_left_hand_TR = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    gsl_matrix_memcpy(alt_left_hand_TR, left_hand_TR);

    // set tran_M_iterative
    for(i = 0; i < tran_MNA_rows; i++){
        for(j = 0; j < tran_MNA_columns; j++){
            if(i == j){
                if((gsl_matrix_get(alt_left_hand_TR, i, j)) != 0){
                    gsl_matrix_set(tran_M_iterative, i, j, (gsl_matrix_get(alt_left_hand_TR, i, j)));
                }
                else{
                    gsl_matrix_set(tran_M_iterative, i, j, 1);
                }
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1.0, alt_left_hand_TR, TR_x_curr, 0.0, tran_Ax_iterative);

    gsl_vector_memcpy(tran_bminusAx_iterative, right_hand_TR);

    gsl_vector_sub(tran_bminusAx_iterative, tran_Ax_iterative);

    gsl_vector_memcpy(tran_r_iterative, tran_bminusAx_iterative);

    gsl_vector_memcpy(tran_r_tilde_iterative, tran_bminusAx_iterative);

    tran_norm_r = gsl_blas_dnrm2(tran_r_iterative);
    tran_norm_b = gsl_blas_dnrm2(right_hand_TR);

    if(tran_norm_b == 0){
        tran_norm_b = 1;
    }

    tran_inverse_matrix();

    while(((tran_norm_r / tran_norm_b) > ITOL) && (tran_iter < (curr_pos_oldname_table + G2_components))){
        tran_iter = tran_iter + 1;

        tran_preconditioner_solve();
        tran_transpose_preconditioner_solve();

        gsl_blas_ddot(tran_r_tilde_iterative, tran_z_iterative, &tran_rho);

        if((fabs(tran_rho)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        if(tran_iter == 1){
            gsl_vector_memcpy(tran_p_iterative, tran_z_iterative);

            gsl_vector_memcpy(tran_p_tilde_iterative, tran_z_tilde_iterative);
        }
        else{
            tran_beta = (tran_rho / tran_rho1);

            gsl_vector_memcpy(tran_p_scale_beta, tran_p_iterative);
            gsl_vector_scale(tran_p_scale_beta, tran_beta);
            gsl_vector_add(tran_p_scale_beta, tran_z_iterative);
            gsl_vector_memcpy(tran_p_iterative, tran_p_scale_beta);
            
            gsl_vector_memcpy(tran_p_tilde_scale_beta, tran_p_tilde_iterative);
            gsl_vector_scale(tran_p_tilde_scale_beta, tran_beta);
            gsl_vector_add(tran_p_tilde_scale_beta, tran_z_tilde_iterative);
            gsl_vector_memcpy(tran_p_tilde_iterative, tran_p_tilde_scale_beta);
        }

        tran_rho1 = tran_rho;

        gsl_blas_dgemv(CblasNoTrans, 1.0, alt_left_hand_TR, tran_p_iterative, 0.0, tran_q_iterative);
        gsl_blas_dgemv(CblasTrans, 1.0, alt_left_hand_TR, tran_p_tilde_iterative, 0.0, tran_q_tilde_iterative);

        gsl_blas_ddot(tran_p_tilde_iterative, tran_q_iterative, &tran_omega);

        if((fabs(tran_omega)) < 1e-14){
            exit(10);
        }

        tran_alpha = (tran_rho / tran_omega);

        gsl_blas_daxpy(tran_alpha, tran_p_iterative, TR_x_curr);

        gsl_blas_daxpy(-tran_alpha, tran_q_iterative, tran_r_iterative);

        gsl_blas_daxpy(-tran_alpha, tran_q_tilde_iterative, tran_r_tilde_iterative);

        tran_norm_r = gsl_blas_dnrm2(tran_r_iterative);
    }

    gsl_vector_memcpy(TR_x_prev, TR_x_curr);

    gsl_vector_memcpy(tran_MNA_e_prev, tran_MNA_e);

    gsl_matrix_free(alt_left_hand_TR);
        
    gsl_matrix_free(tran_M_iterative);
    gsl_matrix_free(tran_M_inverse_iterative);

    gsl_vector_free(tran_r_iterative);
    gsl_vector_free(tran_r_tilde_iterative);
    gsl_vector_free(tran_p_iterative);
    gsl_vector_free(tran_p_tilde_iterative);
    gsl_vector_free(tran_p_scale_beta);
    gsl_vector_free(tran_p_tilde_scale_beta);
    gsl_vector_free(tran_q_iterative);
    gsl_vector_free(tran_q_tilde_iterative);
    gsl_vector_free(tran_z_iterative);
    gsl_vector_free(tran_z_tilde_iterative);
    gsl_vector_free(tran_Ax_iterative);
    gsl_vector_free(tran_bminusAx_iterative);

    gsl_permutation_free(tran_permut_iterative);
}

void TR_tran_CG_sweep(int iter, int k_value){
    int i, j;

    int tran_iter = 0;

    double tran_alpha = 0;
    double tran_beta = 0;

    double tran_norm_r = 0;
    double tran_norm_b = 0;

    double tran_rho = 0;
    double tran_rho1 = 0;

    double tran_ptrans_q = 0;

    tran_r_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_p_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_p_scale_beta = gsl_vector_alloc(tran_MNA_rows);

    tran_q_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_z_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_Ax_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_bminusAx_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_M_iterative = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    tran_M_inverse_iterative = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);

    tran_permut_iterative = gsl_permutation_alloc(tran_MNA_rows);

    for(i = 0; i < tran_MNA_rows; i++){
        for(j = 0; j < tran_MNA_columns; j++){
            gsl_matrix_set(tran_M_iterative, i, j, 0);

            gsl_matrix_set(tran_M_inverse_iterative, i, j, 0);
        }

        gsl_vector_set(tran_r_iterative, i, 0);

        gsl_vector_set(tran_p_iterative, i , 0);
        gsl_vector_set(tran_p_scale_beta, i, 0);

        gsl_vector_set(tran_q_iterative, i, 0);

        gsl_vector_set(tran_z_iterative, i, 0);

        gsl_vector_set(tran_Ax_iterative, i, 0);
        gsl_vector_set(tran_bminusAx_iterative, i, 0);
    }

    if(iter == 1){
        // b = right_hand_BE
        gsl_blas_dgemv(CblasNoTrans, -1.0, right_hand_TR_sub, tran_MNA_X0, 0.0, right_hand_TR);
        gsl_vector_add(right_hand_TR, tran_MNA_e);
    }

    if((iter == 0) && (k_value > 1)){
        // b = right_hand_BE
        gsl_blas_dgemv(CblasNoTrans, -1.0, right_hand_TR_sub, TR_x_prev, 0.0, right_hand_TR);
        gsl_vector_add(right_hand_TR, tran_MNA_e);
        gsl_vector_add(right_hand_TR, tran_MNA_e_prev);
    }

    // initial guess x(0)
    for(i = 0; i < tran_MNA_rows; i++){
        gsl_vector_set(TR_x_curr, i, 0);
    }
    
    // A = alt_left_hand_BE
    alt_left_hand_TR = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    gsl_matrix_memcpy(alt_left_hand_TR, left_hand_TR);

    // set tran_M_iterative
    for(i = 0; i < tran_MNA_rows; i++){
        for(j = 0; j < tran_MNA_columns; j++){
            if(i == j){
                if((gsl_matrix_get(alt_left_hand_TR, i, j)) != 0){
                    gsl_matrix_set(tran_M_iterative, i, j, (gsl_matrix_get(alt_left_hand_TR, i, j)));
                }
                else{
                    gsl_matrix_set(tran_M_iterative, i, j, 1);
                }
            }
        }
    }

    // Ax = tran_Ax_iterative = (left_hand_BE * BE_x_curr)
    gsl_blas_dgemv(CblasNoTrans, 1.0, alt_left_hand_TR, TR_x_curr, 0.0, tran_Ax_iterative);

    // b - Ax = tran_bminusAx_iterative - tran_Ax_iterative
    gsl_vector_memcpy(tran_bminusAx_iterative, right_hand_TR);
    gsl_vector_sub(tran_bminusAx_iterative, tran_Ax_iterative);
        
    // r = tran_r_iterative = b - Ax
    gsl_vector_memcpy(tran_r_iterative, tran_bminusAx_iterative);

    tran_norm_r = gsl_blas_dnrm2(tran_r_iterative);
    tran_norm_b = gsl_blas_dnrm2(right_hand_TR);

    if(tran_norm_b == 0){    // if right_hand_BE = 0
        tran_norm_b = 1;
    }

    // inverse tran_M_iterative
    tran_inverse_matrix();

    while(((tran_norm_r / tran_norm_b) > ITOL) && (tran_iter < (curr_pos_oldname_table + G2_components))){
        tran_iter = tran_iter + 1;

        // tran_M_inverse_iterative * tran_z_iterative = tran_r_iterative
        tran_preconditioner_solve();

        // rho = tran_r_iterative * tran_z_iterative
        gsl_blas_ddot(tran_r_iterative, tran_z_iterative, &tran_rho);

        if(tran_iter == 1){
            // tran_p_iterative = tran_z_iterative
            gsl_vector_memcpy(tran_p_iterative, tran_z_iterative);
        }

        else{
            tran_beta = (tran_rho / tran_rho1);

            // p = tran_p_iterative = z + (b * p) 
            gsl_vector_memcpy(tran_p_scale_beta, tran_p_iterative);
            gsl_vector_scale(tran_p_scale_beta, tran_beta);
            gsl_vector_add(tran_p_scale_beta, tran_z_iterative);    
            gsl_vector_memcpy(tran_p_iterative, tran_p_scale_beta);
        }

        tran_rho1 = tran_rho;

        // q = tran_q_iterative = A * p = left_hand_BE * tran_p_iterative
        gsl_blas_dgemv(CblasNoTrans, 1.0, left_hand_TR, tran_p_iterative, 0.0, tran_q_iterative);

        gsl_blas_ddot(tran_p_iterative, tran_q_iterative, &tran_ptrans_q);
        tran_alpha = (tran_rho / tran_ptrans_q);

        // x = BE_x_curr = x + (alpha * p)
        gsl_blas_daxpy(tran_alpha, tran_p_iterative, TR_x_curr);

        // r = tran_r_iterative = r - (alpha * q)
        gsl_blas_daxpy(-tran_alpha, tran_q_iterative, tran_r_iterative);

        // calculate next step norm
        tran_norm_r = gsl_blas_dnrm2(tran_r_iterative);
    }

    gsl_vector_memcpy(TR_x_prev, TR_x_curr);

    gsl_vector_memcpy(tran_MNA_e_prev, tran_MNA_e);

    gsl_matrix_free(alt_left_hand_TR);
        
    gsl_matrix_free(tran_M_iterative);
    gsl_matrix_free(tran_M_inverse_iterative);

    gsl_vector_free(tran_r_iterative);
    gsl_vector_free(tran_p_iterative);
    gsl_vector_free(tran_p_scale_beta);
    gsl_vector_free(tran_q_iterative);
    gsl_vector_free(tran_z_iterative);
    gsl_vector_free(tran_Ax_iterative);
    gsl_vector_free(tran_bminusAx_iterative);

    gsl_permutation_free(tran_permut_iterative);
}

void SPARSE_trapezoidal_method(){
    init_TR_hands();

    int k, i;
    k_pulse_variable = 0;

    double tk;
    double h;

    int first_iter = 1;

    double lookup_result = 0;

    double tran_node_result = 0;

    char *tran_plot_ptr = NULL;
    int tran_plot_node_name = 0;

    double tran_plot_x = 0;
    double tran_plot_y = 0;

    FILE *tran_fileptr = NULL;
    FILE *tran_gnuplotpipe = NULL;
    char *tran_gnuplotcommands = {"plot 'tranexample.tmp' with linespoints"};

    tran_fileptr = fopen("tranexample.tmp", "w");
    tran_gnuplotpipe = popen("gnuplot -persistent", "w");

    ptr_comp TR_ptr = NULL;

    // init_TR_hands()

    h = (2 / tran_time_step);

    // set left_hand
    // sparse_left_hand_BE = sparse_Gtilde + (2/h) * sparse_Ctilde
    sparse_left_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, h);

    if(plot_node != NULL){
        tran_plot_ptr = strdup(plot_node);
    }
    
    if(print_node != NULL){
        tran_plot_ptr = strdup(print_node);
    }

    if(tran_plot_ptr != NULL){
        tran_plot_node_name = return_newpos_node_name(tran_plot_ptr);
        free(tran_plot_ptr);
    }

    for(k = 1; k <= tran_m_parameter; k++){
        k_tran_variable = 0;

        tk = k * tran_time_step;

        TR_ptr = list_head;

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_SPARSE_e[i] = 0;
        }

        // fill up SPARSE_e(tk)
        while(TR_ptr != NULL){
            if((TR_ptr->comp_type == 'V') || (TR_ptr->comp_type == 'v')){
                if(TR_ptr -> dc_flag == 1){
                    SPARSE_TRAN_sweep_voltagesource(TR_ptr -> value);
                }

                if(TR_ptr -> exp_flag == 1){
                    lookup_result = EXP_lookup(tk, TR_ptr->exp_field.exp_i1, TR_ptr->exp_field.exp_i2, TR_ptr->exp_field.exp_td1, 
                                                TR_ptr->exp_field.exp_td2, TR_ptr->exp_field.exp_tc1, TR_ptr->exp_field.exp_tc2, 
                                                tran_fin_time);

                    SPARSE_TRAN_sweep_voltagesource(lookup_result);
                }

                if(TR_ptr -> pwl_flag == 1){
                    lookup_result = PWL_lookup(tk, tran_fin_time, TR_ptr->comp_pwl_counter, TR_ptr);

                    SPARSE_TRAN_sweep_voltagesource(lookup_result);
                }

                if(TR_ptr -> sin_flag == 1){
                    lookup_result = SIN_lookup(tk, TR_ptr->sin_field.sin_i1, TR_ptr->sin_field.sin_ia, TR_ptr->sin_field.sin_fr,
                                                TR_ptr->sin_field.sin_td, TR_ptr->sin_field.sin_df, TR_ptr->sin_field.sin_ph,
                                                tran_fin_time);
                    
                    SPARSE_TRAN_sweep_voltagesource(lookup_result);
                }

                if(TR_ptr -> pulse_flag == 1){
                    lookup_result = PULSE_lookup(tk, TR_ptr->pulse_field.pulse_i1, TR_ptr->pulse_field.pulse_i2,
                                                    TR_ptr->pulse_field.pulse_td, TR_ptr->pulse_field.pulse_tr,
                                                    TR_ptr->pulse_field.pulse_tf, TR_ptr->pulse_field.pulse_pw, 
                                                    TR_ptr->pulse_field.pulse_per);

                    SPARSE_TRAN_sweep_voltagesource(lookup_result);
                }
            }

            if((TR_ptr->comp_type == 'I') || (TR_ptr->comp_type == 'i')){
                if(TR_ptr -> dc_flag == 1){
                    SPARSE_TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, TR_ptr->value);
                }

                 if(TR_ptr -> exp_flag == 1){
                    lookup_result = EXP_lookup(tk, TR_ptr->exp_field.exp_i1, TR_ptr->exp_field.exp_i2, TR_ptr->exp_field.exp_td1, 
                                                TR_ptr->exp_field.exp_td2, TR_ptr->exp_field.exp_tc1, TR_ptr->exp_field.exp_tc2, 
                                                tran_fin_time);

                    SPARSE_TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(TR_ptr -> pwl_flag == 1){
                    lookup_result = PWL_lookup(tk, tran_fin_time, TR_ptr->comp_pwl_counter, TR_ptr);

                    SPARSE_TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(TR_ptr -> sin_flag == 1){
                    lookup_result = SIN_lookup(tk, TR_ptr->sin_field.sin_i1, TR_ptr->sin_field.sin_ia, TR_ptr->sin_field.sin_fr,
                                                TR_ptr->sin_field.sin_td, TR_ptr->sin_field.sin_df, TR_ptr->sin_field.sin_ph,
                                                tran_fin_time);
                    
                    SPARSE_TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(TR_ptr -> pulse_flag == 1){
                    lookup_result = PULSE_lookup(tk, TR_ptr->pulse_field.pulse_i1, TR_ptr->pulse_field.pulse_i2,
                                                    TR_ptr->pulse_field.pulse_td, TR_ptr->pulse_field.pulse_tr,
                                                    TR_ptr->pulse_field.pulse_tf, TR_ptr->pulse_field.pulse_pw, 
                                                    TR_ptr->pulse_field.pulse_per);

                    SPARSE_TRAN_sweep_currentsource(TR_ptr->pos_node.new_pos_node_name, TR_ptr->neg_node.new_neg_node_name, lookup_result);
                }
            }

            if((TR_ptr->comp_type == 'L') || (TR_ptr->comp_type == 'l')){
                SPARSE_TRAN_sweep_voltagesource(0);
            }

            TR_ptr = TR_ptr -> nxt_comp;
        }

        // solve Ax = b for each step
        if((spd_flag == 0) && (iter_flag == 0)){
            // SPARSE LU
            TR_tran_SPARSE_LU_sweep(first_iter, k);
        
            first_iter = 0;
        }

        if((spd_flag == 1) && (iter_flag == 0)){
            // SPARSE Cholesky
            TR_tran_SPARSE_Cholesky_sweep(first_iter, k);

            first_iter = 0;
        }

        if((spd_flag == 0) && (iter_flag == 1)){
            // SPARSE BiCG
            TR_tran_SPARSE_BiCG_sweep(first_iter, k);

            first_iter = 0;
        }

        if((spd_flag == 1) && (iter_flag == 1)){
            // SPARSE CG
            TR_tran_SPARSE_CG_sweep(first_iter, k);

            first_iter = 0;
        }

        if((tran_plot_node_name != -1) && (tran_plot_node_name != 0)){
            tran_node_result = SPARSE_TR_x_curr[tran_plot_node_name - 1];
        }

        tran_plot_x = tk;

        tran_plot_y = tran_node_result;

        fprintf(tran_fileptr, "%f %f\n", tran_plot_x, tran_plot_y);
    }
    
    fprintf(tran_gnuplotpipe, "%s\n", tran_gnuplotcommands);

    fprintf(tran_gnuplotpipe, "exit\n");

    fclose(tran_fileptr);
    pclose(tran_gnuplotpipe);

    cs_spfree(sparse_left_hand_TR);
}

void TR_tran_SPARSE_LU_sweep(int iter, int k_value){
    int i;

    if(iter == 1){
        SPARSE_TR_x_curr = calloc(tran_SPARSE_rows, sizeof(double));
        SPARSE_TR_x_prev = calloc(tran_SPARSE_rows, sizeof(double));
        
        // sparse_right_hand_TR = Gtilde - ((2/h) * Ctilde)
        sparse_right_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, (-(2 / tran_time_step)));

        // SPARSE_TR_x_curr = (Gtilde - ((2/h) * Ctilde)) * X(0)
        cs_gaxpy(sparse_right_hand_TR, tran_SPARSE_X0, SPARSE_TR_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = tran_SPARSE_e[i] - SPARSE_TR_x_curr[i];
        }
        /*// SPARSE_TR_x_curr = -((2/h) * Ctilde) 
        cs_gaxpy(tran_sparse_tilde_C, -1, SPARSE_TR_x_curr);
        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = SPARSE_TR_x_curr[i] * (2 / tran_time_step);
        }
        
        // SPARSE_TR_x_curr = Gtilde - ((2/h) * Ctilde)
        cs_gaxpy(tran_SPARSE_tilde_G, 1, SPARSE_TR_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = SPARSE_TR_x_curr[i] * tran_SPARSE_X0[i];
        }

        // SPARSE_TR_x_curr = e(tk) - (Gtilde - ((2/h) * Ctilde))
        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = tran_SPARSE_e[i] - SPARSE_TR_x_curr[i];
        }*/
    }

    if((iter == 0) && (k_value > 1)){
        free(SPARSE_TR_x_curr);
        SPARSE_TR_x_curr = calloc(tran_SPARSE_rows, sizeof(double));

        /*cs_gaxpy(tran_sparse_tilde_C, SPARSE_TR_x_prev, SPARSE_TR_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = SPARSE_TR_x_curr[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = SPARSE_TR_x_curr[i] + tran_SPARSE_e[i];
        }*/

        // sparse_right_hand_TR = Gtilde - ((2/h) * Ctilde)
        sparse_right_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, (-(2 / tran_time_step)));

        // SPARSE_TR_x_curr = (Gtilde - ((2/h) * Ctilde)) * X(0)
        cs_gaxpy(sparse_right_hand_TR, SPARSE_TR_x_prev, SPARSE_TR_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = tran_SPARSE_e[i] - SPARSE_TR_x_curr[i];
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = SPARSE_TR_x_curr[i] + tran_SPARSE_e_prev[i];
        }
    }

    tran_S_TR = cs_sqr(2, sparse_left_hand_TR, 0);
    tran_N_TR = cs_lu(sparse_left_hand_TR, tran_S_TR, 1);

    tran_y_TR = cs_malloc(tran_SPARSE_rows, sizeof(double));

    cs_ipvec(tran_N_TR->pinv, SPARSE_TR_x_curr, tran_y_TR, tran_SPARSE_rows);
    cs_lsolve(tran_N_TR->L, tran_y_TR);
    cs_usolve(tran_N_TR->U, tran_y_TR);
    cs_ipvec(tran_S_TR->q, tran_y_TR, SPARSE_TR_x_curr, tran_SPARSE_rows);

    for(i = 0; i < tran_SPARSE_rows; i++){
        SPARSE_TR_x_prev[i] = SPARSE_TR_x_curr[i];

        tran_SPARSE_e_prev[i] = tran_SPARSE_e[i];
    }

    cs_free(tran_y_TR);
    cs_sfree(tran_S_TR);
    cs_nfree(tran_N_TR);

    cs_spfree(sparse_right_hand_TR);
}

void TR_tran_SPARSE_Cholesky_sweep(int iter, int k_value){
    int i;

    if(iter == 1){
        SPARSE_TR_x_curr = calloc(tran_SPARSE_rows, sizeof(double));
        SPARSE_TR_x_prev = calloc(tran_SPARSE_rows, sizeof(double));
        
        // sparse_right_hand_TR = Gtilde - ((2/h) * Ctilde)
        sparse_right_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, (-(2 / tran_time_step)));

        // SPARSE_TR_x_curr = (Gtilde - ((2/h) * Ctilde)) * X(0)
        cs_gaxpy(sparse_right_hand_TR, tran_SPARSE_X0, SPARSE_TR_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = tran_SPARSE_e[i] - SPARSE_TR_x_curr[i];
        }
    }

    if((iter == 0) && (k_value > 1)){
        free(SPARSE_TR_x_curr);
        SPARSE_TR_x_curr = calloc(tran_SPARSE_rows, sizeof(double));

        // sparse_right_hand_TR = Gtilde - ((2/h) * Ctilde)
        sparse_right_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, (-(2 / tran_time_step)));

        // SPARSE_TR_x_curr = (Gtilde - ((2/h) * Ctilde)) * X(0)
        cs_gaxpy(sparse_right_hand_TR, SPARSE_TR_x_prev, SPARSE_TR_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = tran_SPARSE_e[i] - SPARSE_TR_x_curr[i];
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_TR_x_curr[i] = SPARSE_TR_x_curr[i] + tran_SPARSE_e_prev[i];
        }
    }

    //tran_S_TR = cs_sqr(2, sparse_left_hand_TR, 0);
    //tran_N_TR = cs_lu(sparse_left_hand_TR, tran_S_TR, 1);

    tran_S_TR = cs_schol(1, sparse_left_hand_TR);
    tran_N_TR = cs_chol(sparse_left_hand_TR, tran_S_TR);

    tran_y_TR = cs_malloc(tran_SPARSE_rows, sizeof(double));

    cs_ipvec(tran_S_TR->pinv, SPARSE_TR_x_curr, tran_y_TR, tran_SPARSE_rows);
    cs_lsolve(tran_N_TR->L, tran_y_TR);
    cs_ltsolve(tran_N_TR->L, tran_y_TR);
    cs_pvec(tran_S_TR->q, tran_y_TR, SPARSE_TR_x_curr, tran_SPARSE_rows);

    for(i = 0; i < tran_SPARSE_rows; i++){
        SPARSE_TR_x_prev[i] = SPARSE_TR_x_curr[i];

        tran_SPARSE_e_prev[i] = tran_SPARSE_e[i];
    }

    cs_free(tran_y_TR);
    cs_sfree(tran_S_TR);
    cs_nfree(tran_N_TR);

    cs_spfree(sparse_right_hand_TR);
}

void TR_tran_SPARSE_BiCG_sweep(int iter, int k_value){
    int i, j, p;

    int tran_iter = 0;

    double tran_alpha = 0;
    double tran_beta = 0;

    double tran_norm_r = 0;
    double tran_norm_b = 0;

    double tran_rho = 0;
    double tran_rho1 = 0;

    double tran_omega = 0;

    tran_sparse_y_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_y_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_y_iterative failed.\n");

        exit(8);
    }

    tran_sparse_r_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_r_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_r_iterative failed.\n");

        exit(8);
    }

    tran_sparse_r_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_r_tilde_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_r_tilde_iterative failed.\n");

        exit(8);
    }
    
    tran_sparse_p_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_p_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_p_iterative failed.\n");

        exit(8);
    }

    tran_sparse_p_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_p_tilde_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_p_tilde_iterative failed.\n");

        exit(8);
    }

    tran_sparse_q_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_q_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_q_iterative failed.\n");

        exit(8);
    }

    tran_sparse_q_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_q_tilde_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_q_tilde_iterative failed.\n");

        exit(8);
    }

    tran_sparse_z_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_z_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_z_iterative failed.\n");

        exit(8);
    }

    tran_sparse_z_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_z_tilde_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_z_tilde_iterative failed.\n");

        exit(8);
    }

    tran_sparse_m_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_m_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_m_iterative failed.\n");

        exit(8);
    }

    tran_sparse_right_hand_TR = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_right_hand_TR == NULL){
        printf("\nMemory allocation for tran_sparse_right_hand_BE failed.\n");

        exit(8);
    }


    for(j = 0; j < tran_SPARSE_rows; j++){
        for(p = sparse_left_hand_TR -> p[j]; p < sparse_left_hand_TR -> p[ j+ 1]; p++){
            if(sparse_left_hand_TR -> i[p] == j){
                tran_sparse_m_iterative[j] = sparse_left_hand_TR -> x[p];
            }
        }
    }

    tran_M_iterative = gsl_matrix_alloc(tran_SPARSE_rows, tran_SPARSE_rows); // check columns
    tran_M_inverse_iterative = gsl_matrix_alloc(tran_SPARSE_rows, tran_SPARSE_rows);

    for(i = 0; i < tran_SPARSE_rows; i++){
        for(j = 0; j < tran_SPARSE_columns; j++){
            gsl_matrix_set(tran_M_iterative, i, j, 0);
            gsl_matrix_set(tran_M_inverse_iterative, i, j, 0);
        }
    }

    tran_permut_iterative = gsl_permutation_alloc(tran_SPARSE_rows);

    if(iter == 1){
        SPARSE_TR_x_curr = calloc(tran_SPARSE_rows, sizeof(double));
        SPARSE_TR_x_prev = calloc(tran_SPARSE_rows, sizeof(double));

        // alt_tran_sparse_right_hand_TR = Gtilde - ((2/h) * Ctilde)
        alt_tran_sparse_right_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, (-(2 / tran_time_step)));

        // tran_sparse_right_hand_TR = alt_tran_sparse_right_hand_TR * x(0)
        cs_gaxpy(alt_tran_sparse_right_hand_TR, tran_SPARSE_X0, tran_sparse_right_hand_TR);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_TR[i] = tran_SPARSE_e[i] - tran_sparse_right_hand_TR[i];
        }
    }

    if((iter == 0) && (k_value > 1)){
        free(SPARSE_TR_x_curr);
        SPARSE_TR_x_curr = calloc(tran_SPARSE_rows, sizeof(double));

        alt_tran_sparse_right_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, (-(2 / tran_time_step)));

        // tran_sparse_right_hand_TR = alt_tran_sparse_right_hand_TR * x(0)
        cs_gaxpy(alt_tran_sparse_right_hand_TR, SPARSE_TR_x_prev, tran_sparse_right_hand_TR);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_TR[i] = tran_SPARSE_e[i] - tran_sparse_right_hand_TR[i];
        }
        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_TR[i] = tran_sparse_right_hand_TR[i] + tran_SPARSE_e_prev[i];
        }
    }

    // y = A * x
    for (j = 0 ; j < tran_SPARSE_rows ; j++){
        for (p = sparse_left_hand_TR -> p[j] ; p < sparse_left_hand_TR -> p[j+1] ; p++){
            tran_sparse_y_iterative[sparse_left_hand_TR -> i[p]] = tran_sparse_y_iterative[sparse_left_hand_TR -> i[p]] + sparse_left_hand_TR -> x[p] * SPARSE_TR_x_curr[j];
        }
    }

    // r1 = b - y
    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_sparse_r_iterative[i] = tran_sparse_right_hand_TR[i] - tran_sparse_y_iterative[i];
    }

    // tran_sparse_r_tilde_iterative = tran_sparse_r_iterative
    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_sparse_r_tilde_iterative[i] = tran_sparse_r_iterative[i];
    }

    tran_norm_r = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_r_iterative, 1);
    tran_norm_b = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_right_hand_TR, 1);

    if(tran_norm_b == 0){
        tran_norm_b = 1;
    }

    tran_sparse_inverse_matrix();

    tran_m_cmprsd_iterative = cs_spalloc(tran_SPARSE_rows, tran_SPARSE_rows, tran_SPARSE_rows, 1, 1);
    tran_m_cmprsd_iterative->nz = tran_SPARSE_rows;

    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_m_cmprsd_iterative -> i[i] = i;
        tran_m_cmprsd_iterative -> p[i] = i;
        tran_m_cmprsd_iterative -> x[i] = gsl_matrix_get(tran_M_inverse_iterative, i, i);
    }

    tran_m_cc_iterative = cs_compress(tran_m_cmprsd_iterative);
    cs_dupl(tran_m_cc_iterative);

    while(((tran_norm_r / tran_norm_b) > ITOL) && (tran_iter < tran_SPARSE_rows)){
        tran_iter = tran_iter + 1;

        tran_sparse_preconditioner_solve();
        tran_sparse_transpose_preconditioner();


        tran_rho = dot_product(tran_sparse_r_tilde_iterative, tran_sparse_z_iterative, tran_SPARSE_rows);

        if((fabs(tran_rho)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        if(tran_iter == 1){
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_sparse_p_iterative[i] = tran_sparse_z_iterative[i];

                tran_sparse_p_tilde_iterative[i] = tran_sparse_z_tilde_iterative[i];
            }
        }

        else{
            tran_beta = (tran_rho / tran_rho1);
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_sparse_p_iterative[i] = (tran_beta * tran_sparse_p_iterative[i]) + tran_sparse_z_iterative[i];

                tran_sparse_p_tilde_iterative[i] = (tran_beta * tran_sparse_p_tilde_iterative[i]) + tran_sparse_z_tilde_iterative[i];
            }

        }

        tran_rho1 = tran_rho;

        for(j = 0; j < tran_SPARSE_rows; j++){
            for (p = sparse_left_hand_TR -> p[j]; p < sparse_left_hand_TR -> p[j + 1] ; p++){
                tran_sparse_q_iterative[sparse_left_hand_TR -> i[p]] = tran_sparse_q_iterative[sparse_left_hand_TR -> i[p]]
                                                        + (sparse_left_hand_TR -> x[p] * tran_sparse_p_iterative[j]);
            }
        }

        for(j = 0; j < tran_SPARSE_rows; j++){
            tran_sparse_q_tilde_iterative[j] = 0;

            for (p = sparse_left_hand_TR -> p[j]; p < sparse_left_hand_TR -> p[j + 1] ; p++){
                tran_sparse_q_tilde_iterative[j] = tran_sparse_q_tilde_iterative[j]
                                                        + (sparse_left_hand_TR -> x[p] * tran_sparse_p_tilde_iterative[sparse_left_hand_TR->i[p]]);
            }
        }

        tran_omega = dot_product(tran_sparse_p_tilde_iterative, tran_sparse_q_iterative, tran_SPARSE_rows);

        if((fabs(tran_omega)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        tran_alpha = (tran_rho / tran_omega);

        cblas_daxpy(tran_SPARSE_rows, tran_alpha, tran_sparse_p_iterative, 1, SPARSE_TR_x_curr, 1);

        cblas_daxpy(tran_SPARSE_rows, -tran_alpha, tran_sparse_q_iterative, 1, tran_sparse_r_iterative, 1);

        cblas_daxpy(tran_SPARSE_rows, -tran_alpha, tran_sparse_q_tilde_iterative, 1, tran_sparse_r_tilde_iterative, 1);

        tran_norm_r = cblas_dnrm2(sparse_rowsA, tran_sparse_r_iterative, 1);

        // reset q_iterative vector to zero after each step
        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_q_iterative[i] = 0;
        }

    }


    for(i = 0; i < tran_SPARSE_rows; i++){
        SPARSE_TR_x_prev[i] = SPARSE_TR_x_curr[i];

        tran_SPARSE_e_prev[i] = tran_SPARSE_e[i];
    }

    gsl_permutation_free(tran_permut_iterative);

    gsl_matrix_free(tran_M_iterative);
    gsl_matrix_free(tran_M_inverse_iterative);

    cs_spfree(tran_m_cmprsd_iterative);
    cs_spfree(tran_m_cc_iterative);

    cs_spfree(alt_tran_sparse_right_hand_TR);

    free(tran_sparse_y_iterative);
    free(tran_sparse_right_hand_TR);
    free(tran_sparse_m_iterative);
    free(tran_sparse_z_tilde_iterative);
    free(tran_sparse_z_iterative);
    free(tran_sparse_q_tilde_iterative);
    free(tran_sparse_q_iterative);
    free(tran_sparse_p_iterative);
    free(tran_sparse_p_tilde_iterative);
    free(tran_sparse_r_iterative);
    free(tran_sparse_r_tilde_iterative);
}

void TR_tran_SPARSE_CG_sweep(int iter, int k_value){
    int i, j, p;

    int tran_iter = 0;

    double tran_norm_r = 0;
    double tran_norm_b = 0;

    double tran_rho = 0;
    double tran_rho1 = 0;
    double tran_beta = 0;
    double tran_alpha = 0;

    tran_sparse_y_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_y_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_y_iterative failed.\n");

        exit(8);
    }

    tran_sparse_r_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_r_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_r_iterative failed.\n");

        exit(8);
    }

    tran_sparse_p_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_p_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_p_iterative failed.\n");

        exit(8);
    }

    tran_sparse_z_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_z_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_z_iterative failed.\n");

        exit(8);
    }

    tran_sparse_q_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_q_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_q_iterative failed.\n");

        exit(8);
    }

    tran_sparse_m_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_m_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_m_iterative failed.\n");

        exit(8);
    }
    
    for(j = 0; j < tran_SPARSE_rows; j++){
        for(p = sparse_left_hand_TR -> p[j]; p < sparse_left_hand_TR -> p[ j+ 1]; p++){
            if(sparse_left_hand_TR -> i[p] == j){
                tran_sparse_m_iterative[j] = sparse_left_hand_TR -> x[p];
            }
        }
    }

    tran_sparse_right_hand_TR = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_right_hand_TR == NULL){
        printf("\nMemory allocation for tran_sparse_right_hand_BE failed.\n");

        exit(8);
    }

    tran_M_iterative = gsl_matrix_alloc(tran_SPARSE_rows, tran_SPARSE_rows); // check columns
    tran_M_inverse_iterative = gsl_matrix_alloc(tran_SPARSE_rows, tran_SPARSE_rows);

    for(i = 0; i < tran_SPARSE_rows; i++){
        for(j = 0; j < tran_SPARSE_columns; j++){
            gsl_matrix_set(tran_M_iterative, i, j, 0);
            gsl_matrix_set(tran_M_inverse_iterative, i, j, 0);
        }
    }

    tran_permut_iterative = gsl_permutation_alloc(tran_SPARSE_rows);

    if(iter == 1){
        SPARSE_TR_x_curr = calloc(tran_SPARSE_rows, sizeof(double));
        SPARSE_TR_x_prev = calloc(tran_SPARSE_rows, sizeof(double));

        // alt_tran_sparse_right_hand_TR = Gtilde - ((2/h) * Ctilde)
        alt_tran_sparse_right_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, (-(2 / tran_time_step)));

        // tran_sparse_right_hand_TR = alt_tran_sparse_right_hand_TR * x(0)
        cs_gaxpy(alt_tran_sparse_right_hand_TR, tran_SPARSE_X0, tran_sparse_right_hand_TR);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_TR[i] = tran_SPARSE_e[i] - tran_sparse_right_hand_TR[i];
        }

        /*cs_gaxpy(tran_sparse_tilde_C, tran_SPARSE_X0, tran_sparse_right_hand_BE);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] + tran_SPARSE_e[i];
        }*/
    }

    if((iter == 0) && (k_value > 1)){
        free(SPARSE_TR_x_curr);
        SPARSE_TR_x_curr = calloc(tran_SPARSE_rows, sizeof(double));

        alt_tran_sparse_right_hand_TR = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, (-(2 / tran_time_step)));

        // tran_sparse_right_hand_TR = alt_tran_sparse_right_hand_TR * x(0)
        cs_gaxpy(alt_tran_sparse_right_hand_TR, SPARSE_TR_x_prev, tran_sparse_right_hand_TR);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_TR[i] = tran_SPARSE_e[i] - tran_sparse_right_hand_TR[i];
        }
        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_TR[i] = tran_sparse_right_hand_TR[i] + tran_SPARSE_e_prev[i];
        }

        /*cs_gaxpy(tran_sparse_tilde_C, SPARSE_BE_x_prev, tran_sparse_right_hand_BE);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] + tran_SPARSE_e[i];
        }*/
    }

     // y = A * x
    for (j = 0 ; j < tran_SPARSE_rows ; j++){
        for (p = sparse_left_hand_TR -> p[j] ; p < sparse_left_hand_TR -> p[j+1] ; p++){
            tran_sparse_y_iterative[sparse_left_hand_TR -> i[p]] = tran_sparse_y_iterative[sparse_left_hand_TR -> i[p]] + sparse_left_hand_TR -> x[p] * SPARSE_TR_x_curr[j];
        }
    }

    // r1 = b - y
    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_sparse_r_iterative[i] = tran_sparse_right_hand_TR[i] - tran_sparse_y_iterative[i];
    }

    tran_norm_r = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_r_iterative, 1);
    tran_norm_b = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_right_hand_TR, 1);

    if(tran_norm_b == 0){
        tran_norm_b = 1;
    }

    tran_sparse_inverse_matrix();

    tran_m_cmprsd_iterative = cs_spalloc(tran_SPARSE_rows, tran_SPARSE_rows, tran_SPARSE_rows, 1, 1);
    tran_m_cmprsd_iterative->nz = tran_SPARSE_rows;

    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_m_cmprsd_iterative -> i[i] = i;
        tran_m_cmprsd_iterative -> p[i] = i;
        tran_m_cmprsd_iterative -> x[i] = gsl_matrix_get(tran_M_inverse_iterative, i, i);
    }

    tran_m_cc_iterative = cs_compress(tran_m_cmprsd_iterative);
    cs_dupl(tran_m_cc_iterative);

    while(((tran_norm_r / tran_norm_b) > ITOL) && (tran_iter < tran_SPARSE_rows)){
        tran_iter = tran_iter + 1;

        tran_sparse_preconditioner_solve();

        tran_rho = dot_product(tran_sparse_r_iterative, tran_sparse_z_iterative, tran_SPARSE_rows);

        if(tran_iter == 1){
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_sparse_p_iterative[i] = tran_sparse_z_iterative[i];
            }
        }

        else{
            tran_beta = (tran_rho / tran_rho1);

            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_sparse_p_iterative[i] = (tran_beta * tran_sparse_p_iterative[i]) + tran_sparse_z_iterative[i];
            }
        }

        tran_rho1 = tran_rho;

        for(j = 0; j < tran_SPARSE_rows; j++){
            for (p = sparse_left_hand_TR -> p[j]; p < sparse_left_hand_TR -> p[j+1] ; p++){
                tran_sparse_q_iterative[sparse_left_hand_TR -> i[p]] = tran_sparse_q_iterative[sparse_left_hand_TR -> i[p]]
                                                        + (sparse_left_hand_TR -> x[p] * tran_sparse_p_iterative[j]);
            }
        }

        tran_alpha = (tran_rho / dot_product(tran_sparse_p_iterative, tran_sparse_q_iterative, tran_SPARSE_rows));

        cblas_daxpy(tran_SPARSE_rows, tran_alpha, tran_sparse_p_iterative, 1, SPARSE_TR_x_curr, 1);

        cblas_daxpy(tran_SPARSE_rows, -tran_alpha, tran_sparse_q_iterative, 1, tran_sparse_r_iterative, 1);

        tran_norm_r = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_r_iterative, 1);

        // reset tran_sparse_q_iterative vector to zero after each step
        for(i = 0; i < sparse_rowsA; i++){
            tran_sparse_q_iterative[i] = 0;
        }
    }

    for(i = 0; i < tran_SPARSE_rows; i++){
        SPARSE_TR_x_prev[i] = SPARSE_TR_x_curr[i];

        tran_SPARSE_e_prev[i] = tran_SPARSE_e[i];
    }

    gsl_permutation_free(tran_permut_iterative);

    gsl_matrix_free(tran_M_iterative);
    gsl_matrix_free(tran_M_inverse_iterative);

    cs_spfree(tran_m_cmprsd_iterative);
    cs_spfree(tran_m_cc_iterative);

    cs_spfree(alt_tran_sparse_right_hand_TR);

    free(tran_sparse_right_hand_TR);
    free(tran_sparse_y_iterative);
    free(tran_sparse_r_iterative);
    free(tran_sparse_p_iterative);
    free(tran_sparse_m_iterative);
    free(tran_sparse_z_iterative);
    free(tran_sparse_q_iterative);
}

void SPARSE_backwardeuler_method(){
    int k, i;
    k_pulse_variable = 0;

    double tk;
    double h;

    int first_iter = 1;

    double lookup_result = 0;

    double tran_node_result = 0;

    char *tran_plot_ptr = NULL;
    int tran_plot_node_name = 0;

    long double tran_plot_x = 0;
    double tran_plot_y = 0;

    FILE *tran_fileptr = NULL;
    FILE *tran_gnuplotpipe = NULL;
    char *tran_gnuplotcommands = {"plot 'tranexample.tmp' with linespoints"};

    tran_fileptr = fopen("tranexample.tmp", "w");
    tran_gnuplotpipe = popen("gnuplot -persistent", "w");

    ptr_comp BE_ptr = NULL;

    init_BE_hands();

    h = (1 / tran_time_step);

    // set left_hand
    // sparse_left_hand_BE = sparse_Gtilde + (1/h) * sparse_Ctilde
    sparse_left_hand_BE = cs_add(tran_SPARSE_tilde_G, tran_sparse_tilde_C, 1, h);

    if(plot_node != NULL){
        tran_plot_ptr = strdup(plot_node);
    }
    
    if(print_node != NULL){
        tran_plot_ptr = strdup(print_node);
    }

    if(tran_plot_ptr != NULL){
        tran_plot_node_name = return_newpos_node_name(tran_plot_ptr);
        free(tran_plot_ptr);
    }

    for(k = 1; k <= tran_m_parameter; k++){
        k_tran_variable = 0;

        tk = k * tran_time_step;

        BE_ptr = list_head;

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_SPARSE_e[i] = 0;
        }

        // fill up SPARSE_e(tk)
        while(BE_ptr != NULL){
            if((BE_ptr->comp_type == 'V') || (BE_ptr->comp_type == 'v')){
                if(BE_ptr -> dc_flag == 1){
                    SPARSE_TRAN_sweep_voltagesource(BE_ptr -> value);
                }

                if(BE_ptr -> exp_flag == 1){
                    lookup_result = EXP_lookup(tk, BE_ptr->exp_field.exp_i1, BE_ptr->exp_field.exp_i2, BE_ptr->exp_field.exp_td1, 
                                                BE_ptr->exp_field.exp_td2, BE_ptr->exp_field.exp_tc1, BE_ptr->exp_field.exp_tc2, 
                                                tran_fin_time);

                    SPARSE_TRAN_sweep_voltagesource(lookup_result);
                }

                if(BE_ptr -> pwl_flag == 1){
                    lookup_result = PWL_lookup(tk, tran_fin_time, BE_ptr->comp_pwl_counter, BE_ptr);

                    SPARSE_TRAN_sweep_voltagesource(lookup_result);
                }

                if(BE_ptr -> sin_flag == 1){
                    lookup_result = SIN_lookup(tk, BE_ptr->sin_field.sin_i1, BE_ptr->sin_field.sin_ia, BE_ptr->sin_field.sin_fr,
                                                BE_ptr->sin_field.sin_td, BE_ptr->sin_field.sin_df, BE_ptr->sin_field.sin_ph,
                                                tran_fin_time);
                    
                    SPARSE_TRAN_sweep_voltagesource(lookup_result);
                }

                if(BE_ptr -> pulse_flag == 1){
                    lookup_result = PULSE_lookup(tk, BE_ptr->pulse_field.pulse_i1, BE_ptr->pulse_field.pulse_i2,
                                                    BE_ptr->pulse_field.pulse_td, BE_ptr->pulse_field.pulse_tr,
                                                    BE_ptr->pulse_field.pulse_tf, BE_ptr->pulse_field.pulse_pw, 
                                                    BE_ptr->pulse_field.pulse_per);

                    SPARSE_TRAN_sweep_voltagesource(lookup_result);
                }
            }

            if((BE_ptr->comp_type == 'I') || (BE_ptr->comp_type == 'i')){
                if(BE_ptr -> dc_flag == 1){
                    SPARSE_TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, BE_ptr->value);
                }

                 if(BE_ptr -> exp_flag == 1){
                    lookup_result = EXP_lookup(tk, BE_ptr->exp_field.exp_i1, BE_ptr->exp_field.exp_i2, BE_ptr->exp_field.exp_td1, 
                                                BE_ptr->exp_field.exp_td2, BE_ptr->exp_field.exp_tc1, BE_ptr->exp_field.exp_tc2, 
                                                tran_fin_time);

                    SPARSE_TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(BE_ptr -> pwl_flag == 1){
                    lookup_result = PWL_lookup(tk, tran_fin_time, BE_ptr->comp_pwl_counter, BE_ptr);

                    SPARSE_TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(BE_ptr -> sin_flag == 1){
                    lookup_result = SIN_lookup(tk, BE_ptr->sin_field.sin_i1, BE_ptr->sin_field.sin_ia, BE_ptr->sin_field.sin_fr,
                                                BE_ptr->sin_field.sin_td, BE_ptr->sin_field.sin_df, BE_ptr->sin_field.sin_ph,
                                                tran_fin_time);
                    
                    SPARSE_TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(BE_ptr -> pulse_flag == 1){
                    lookup_result = PULSE_lookup(tk, BE_ptr->pulse_field.pulse_i1, BE_ptr->pulse_field.pulse_i2,
                                                    BE_ptr->pulse_field.pulse_td, BE_ptr->pulse_field.pulse_tr,
                                                    BE_ptr->pulse_field.pulse_tf, BE_ptr->pulse_field.pulse_pw, 
                                                    BE_ptr->pulse_field.pulse_per);

                    SPARSE_TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, lookup_result);
                }
            }

            if((BE_ptr->comp_type == 'L') || (BE_ptr->comp_type == 'l')){
                SPARSE_TRAN_sweep_voltagesource(0);
            }

            BE_ptr = BE_ptr -> nxt_comp;
        }

        // solve Ax = b for each step
        if((spd_flag == 0) && (iter_flag == 0)){
            // SPARSE LU
            tran_SPARSE_LU_sweep(first_iter, k);
        
            first_iter = 0;
        }

        if((spd_flag == 1) && (iter_flag == 0)){
            // SPARSE Cholesky
            tran_SPARSE_Cholesky_sweep(first_iter, k);

            first_iter = 0;
        }

        if((spd_flag == 0) && (iter_flag == 1)){
            // SPARSE BiCG
            tran_SPARSE_BiCG_sweep(first_iter, k);

            first_iter = 0;
        }

        if((spd_flag == 1) && (iter_flag == 1)){
            // SPARSE CG
            tran_SPARSE_CG_sweep(first_iter, k);

            first_iter = 0;
        }

        if((tran_plot_node_name != -1) && (tran_plot_node_name != 0)){
            tran_node_result = SPARSE_BE_x_curr[tran_plot_node_name - 1];
        }

        tran_plot_x = tk;

        tran_plot_y = tran_node_result;

        fprintf(tran_fileptr, "%LF %lf\n", tran_plot_x, tran_plot_y);
    }

    fprintf(tran_gnuplotpipe, "%s\n", tran_gnuplotcommands);

    fprintf(tran_gnuplotpipe, "exit\n");

    fclose(tran_fileptr);
    pclose(tran_gnuplotpipe);

    cs_spfree(sparse_left_hand_BE);
}

void backwardeuler_method(){
    int k, i;
    k_pulse_variable = 0;

    double tk;
    double h;

    int first_iter = 1;

    double lookup_result = 0;

    double tran_node_result = 0;

    char *tran_plot_ptr = NULL;
    int tran_plot_node_name = 0;

    double tran_plot_x = 0;
    double tran_plot_y = 0;

    FILE *tran_fileptr = NULL;
    FILE *tran_gnuplotpipe = NULL;
    char *tran_gnuplotcommands = {"plot 'tranexample.tmp' with linespoints"};

    tran_fileptr = fopen("tranexample.tmp", "w");
    tran_gnuplotpipe = popen("gnuplot -persistent", "w");

    ptr_comp BE_ptr = NULL;

    init_BE_hands();

    h = (1 / tran_time_step);

    // set left_hand
    gsl_matrix_scale(tran_MNA_C, h);    // (1 / h) * Ctilde

    gsl_matrix_add(tran_MNA_G, tran_MNA_C); // (Gtilde + ((1 / h) * Ctilde))

    gsl_matrix_memcpy(left_hand_BE, tran_MNA_G);

    if(plot_node != NULL){
        tran_plot_ptr = strdup(plot_node);
    }
    
    if(print_node != NULL){
        tran_plot_ptr = strdup(print_node);
    }

    if(tran_plot_ptr != NULL){
        tran_plot_node_name = return_newpos_node_name(tran_plot_ptr);
        free(tran_plot_ptr);
    }

    for(k = 1; k <= tran_m_parameter; k++){
        k_tran_variable = 0;
        //k_pulse_variable = 0;

        //first_iter = 1;

        tk = k * tran_time_step;

        BE_ptr = list_head;

        for(i = 0; i < tran_MNA_rows; i++){
            gsl_vector_set(tran_MNA_e, i, 0);
        }

        // fill up e(tk)
        while(BE_ptr != NULL){
            if((BE_ptr->comp_type == 'V') || (BE_ptr->comp_type == 'v')){
                if(BE_ptr -> dc_flag == 1){
                    TRAN_sweep_voltagesource(BE_ptr -> value);
                }

                if(BE_ptr -> exp_flag == 1){
                    lookup_result = EXP_lookup(tk, BE_ptr->exp_field.exp_i1, BE_ptr->exp_field.exp_i2, BE_ptr->exp_field.exp_td1, 
                                                BE_ptr->exp_field.exp_td2, BE_ptr->exp_field.exp_tc1, BE_ptr->exp_field.exp_tc2, 
                                                tran_fin_time);

                    TRAN_sweep_voltagesource(lookup_result);
                }

                if(BE_ptr -> pwl_flag == 1){
                    lookup_result = PWL_lookup(tk, tran_fin_time, BE_ptr->comp_pwl_counter, BE_ptr);

                    TRAN_sweep_voltagesource(lookup_result);
                }

                if(BE_ptr -> sin_flag == 1){
                    lookup_result = SIN_lookup(tk, BE_ptr->sin_field.sin_i1, BE_ptr->sin_field.sin_ia, BE_ptr->sin_field.sin_fr,
                                                BE_ptr->sin_field.sin_td, BE_ptr->sin_field.sin_df, BE_ptr->sin_field.sin_ph,
                                                tran_fin_time);
                    
                    TRAN_sweep_voltagesource(lookup_result);
                }

                if(BE_ptr -> pulse_flag == 1){
                    lookup_result = PULSE_lookup(tk, BE_ptr->pulse_field.pulse_i1, BE_ptr->pulse_field.pulse_i2,
                                                    BE_ptr->pulse_field.pulse_td, BE_ptr->pulse_field.pulse_tr,
                                                    BE_ptr->pulse_field.pulse_tf, BE_ptr->pulse_field.pulse_pw, 
                                                    BE_ptr->pulse_field.pulse_per);

                    TRAN_sweep_voltagesource(lookup_result);
                }
            }

            if((BE_ptr->comp_type == 'I') || (BE_ptr->comp_type == 'i')){
                if(BE_ptr -> dc_flag == 1){
                    TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, BE_ptr->value);
                }

                 if(BE_ptr -> exp_flag == 1){
                    lookup_result = EXP_lookup(tk, BE_ptr->exp_field.exp_i1, BE_ptr->exp_field.exp_i2, BE_ptr->exp_field.exp_td1, 
                                                BE_ptr->exp_field.exp_td2, BE_ptr->exp_field.exp_tc1, BE_ptr->exp_field.exp_tc2, 
                                                tran_fin_time);

                    TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(BE_ptr -> pwl_flag == 1){
                    lookup_result = PWL_lookup(tk, tran_fin_time, BE_ptr->comp_pwl_counter, BE_ptr);

                    TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(BE_ptr -> sin_flag == 1){
                    lookup_result = SIN_lookup(tk, BE_ptr->sin_field.sin_i1, BE_ptr->sin_field.sin_ia, BE_ptr->sin_field.sin_fr,
                                                BE_ptr->sin_field.sin_td, BE_ptr->sin_field.sin_df, BE_ptr->sin_field.sin_ph,
                                                tran_fin_time);
                    
                    TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, lookup_result);
                }

                if(BE_ptr -> pulse_flag == 1){
                    lookup_result = PULSE_lookup(tk, BE_ptr->pulse_field.pulse_i1, BE_ptr->pulse_field.pulse_i2,
                                                    BE_ptr->pulse_field.pulse_td, BE_ptr->pulse_field.pulse_tr,
                                                    BE_ptr->pulse_field.pulse_tf, BE_ptr->pulse_field.pulse_pw, 
                                                    BE_ptr->pulse_field.pulse_per);

                    TRAN_sweep_currentsource(BE_ptr->pos_node.new_pos_node_name, BE_ptr->neg_node.new_neg_node_name, lookup_result);
                }
            }

            if((BE_ptr->comp_type == 'L') || (BE_ptr->comp_type == 'l')){
                TRAN_sweep_voltagesource(0);
            }

            BE_ptr = BE_ptr -> nxt_comp;
        }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // solve Ax = b for each step
        /*printf("\ne(tk)\n:");
        for(i = 0; i < tran_MNA_rows; i++){
            printf("%lf\n",gsl_vector_get(tran_MNA_e, i));
        }*/


        if((spd_flag == 0) && (iter_flag == 0)){
            // LU
            tran_LU_sweep(first_iter, k);

                first_iter = 0;
            }

            if((spd_flag == 1) && (iter_flag == 0)){
                // Cholesky
                tran_Cholesky_sweep(first_iter, k);

                first_iter = 0;
            }

            if((spd_flag == 0) && (iter_flag == 1)){
                // BiCG
                tran_BiCG_sweep(first_iter, k);

                first_iter = 0;
            }

            if((spd_flag == 1) && (iter_flag == 1)){
                // CG
                tran_CG_sweep(first_iter, k);

                first_iter = 0;
            }
        

        if((tran_plot_node_name != -1) && (tran_plot_node_name != 0)){
            tran_node_result = gsl_vector_get(BE_x_curr, (tran_plot_node_name - 1));
        }

        tran_plot_x = tk;

        tran_plot_y = tran_node_result;

        fprintf(tran_fileptr, "%f %f\n", tran_plot_x, tran_plot_y);
    }

    fprintf(tran_gnuplotpipe, "%s\n", tran_gnuplotcommands);

    fprintf(tran_gnuplotpipe, "exit\n");

    fclose(tran_fileptr);
    pclose(tran_gnuplotpipe);
}

void tran_SPARSE_LU_sweep(int iter, int k_value){
    int i;

    if(iter == 1){
        SPARSE_BE_x_curr = calloc(tran_SPARSE_rows, sizeof(double));
        SPARSE_BE_x_prev = calloc(tran_SPARSE_rows, sizeof(double));
        
        cs_gaxpy(tran_sparse_tilde_C, tran_SPARSE_X0, SPARSE_BE_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_BE_x_curr[i] = SPARSE_BE_x_curr[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_BE_x_curr[i] = SPARSE_BE_x_curr[i] + tran_SPARSE_e[i];
        }
    }

    if((iter == 0) && (k_value > 1)){
        free(SPARSE_BE_x_curr);
        SPARSE_BE_x_curr = calloc(tran_SPARSE_rows, sizeof(double));

        cs_gaxpy(tran_sparse_tilde_C, SPARSE_BE_x_prev, SPARSE_BE_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_BE_x_curr[i] = SPARSE_BE_x_curr[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_BE_x_curr[i] = SPARSE_BE_x_curr[i] + tran_SPARSE_e[i];
        }
    }

    tran_S = cs_sqr(2, sparse_left_hand_BE, 0);
    tran_N = cs_lu(sparse_left_hand_BE, tran_S, 1);

    tran_y = cs_malloc(tran_SPARSE_rows, sizeof(double));

    cs_ipvec(tran_N->pinv, SPARSE_BE_x_curr, tran_y, tran_SPARSE_rows);
    cs_lsolve(tran_N->L, tran_y);
    cs_usolve(tran_N->U, tran_y);
    cs_ipvec(tran_S->q, tran_y, SPARSE_BE_x_curr, tran_SPARSE_rows);

    for(i = 0; i < tran_SPARSE_rows; i++){
        SPARSE_BE_x_prev[i] = SPARSE_BE_x_curr[i];
    }

    cs_free(tran_y);
    cs_sfree(tran_S);
    cs_nfree(tran_N);
}

void tran_SPARSE_Cholesky_sweep(int iter, int k_value){
    int i;

    if(iter == 1){
        SPARSE_BE_x_curr = calloc(tran_SPARSE_rows, sizeof(double));
        SPARSE_BE_x_prev = calloc(tran_SPARSE_rows, sizeof(double));
        
        cs_gaxpy(tran_sparse_tilde_C, tran_SPARSE_X0, SPARSE_BE_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_BE_x_curr[i] = SPARSE_BE_x_curr[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_BE_x_curr[i] = SPARSE_BE_x_curr[i] + tran_SPARSE_e[i];
        }
    }

    if((iter == 0) && (k_value > 1)){
        free(SPARSE_BE_x_curr);
        SPARSE_BE_x_curr = calloc(tran_SPARSE_rows, sizeof(double));

        cs_gaxpy(tran_sparse_tilde_C, SPARSE_BE_x_prev, SPARSE_BE_x_curr);

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_BE_x_curr[i] = SPARSE_BE_x_curr[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            SPARSE_BE_x_curr[i] = SPARSE_BE_x_curr[i] + tran_SPARSE_e[i];
        }
    }

    tran_S = cs_schol(1, sparse_left_hand_BE);
    tran_N = cs_chol(sparse_left_hand_BE, tran_S);

    tran_y = cs_malloc(tran_SPARSE_rows, sizeof(double));

    cs_ipvec(tran_S->pinv, SPARSE_BE_x_curr, tran_y, tran_SPARSE_rows);
    cs_lsolve(tran_N->L, tran_y);
    cs_ltsolve(tran_N->L, tran_y);
    cs_pvec(tran_S->q, tran_y, SPARSE_BE_x_curr, tran_SPARSE_rows);

    for(i = 0; i < tran_SPARSE_rows; i++){
        SPARSE_BE_x_prev[i] = SPARSE_BE_x_curr[i];
    }

    cs_free(tran_y);
    cs_sfree(tran_S);
    cs_nfree(tran_N);
}

void tran_SPARSE_BiCG_sweep(int iter, int k_value){
    int i, j, p;

    int tran_iter = 0;

    double tran_alpha = 0;
    double tran_beta = 0;

    double tran_norm_r = 0;
    double tran_norm_b = 0;

    double tran_rho = 0;
    double tran_rho1 = 0;

    double tran_omega = 0;

    tran_sparse_y_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_y_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_y_iterative failed.\n");

        exit(8);
    }

    tran_sparse_r_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_r_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_r_iterative failed.\n");

        exit(8);
    }

    tran_sparse_r_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_r_tilde_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_r_tilde_iterative failed.\n");

        exit(8);
    }
    
    tran_sparse_p_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_p_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_p_iterative failed.\n");

        exit(8);
    }

    tran_sparse_p_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_p_tilde_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_p_tilde_iterative failed.\n");

        exit(8);
    }

    tran_sparse_q_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_q_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_q_iterative failed.\n");

        exit(8);
    }

    tran_sparse_q_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_q_tilde_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_q_tilde_iterative failed.\n");

        exit(8);
    }

    tran_sparse_z_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_z_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_z_iterative failed.\n");

        exit(8);
    }

    tran_sparse_z_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_z_tilde_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_z_tilde_iterative failed.\n");

        exit(8);
    }

    tran_sparse_m_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_m_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_m_iterative failed.\n");

        exit(8);
    }

    tran_sparse_right_hand_BE = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_right_hand_BE == NULL){
        printf("\nMemory allocation for tran_sparse_right_hand_BE failed.\n");

        exit(8);
    }


    for(j = 0; j < tran_SPARSE_rows; j++){
        for(p = sparse_left_hand_BE -> p[j]; p < sparse_left_hand_BE -> p[ j+ 1]; p++){
            if(sparse_left_hand_BE -> i[p] == j){
                tran_sparse_m_iterative[j] = sparse_left_hand_BE -> x[p];
            }
        }
    }

    tran_M_iterative = gsl_matrix_alloc(tran_SPARSE_rows, tran_SPARSE_rows); // check columns
    tran_M_inverse_iterative = gsl_matrix_alloc(tran_SPARSE_rows, tran_SPARSE_rows);

    for(i = 0; i < tran_SPARSE_rows; i++){
        for(j = 0; j < tran_SPARSE_columns; j++){
            gsl_matrix_set(tran_M_iterative, i, j, 0);
            gsl_matrix_set(tran_M_inverse_iterative, i, j, 0);
        }
    }

    tran_permut_iterative = gsl_permutation_alloc(tran_SPARSE_rows);

    if(iter == 1){
        SPARSE_BE_x_curr = calloc(tran_SPARSE_rows, sizeof(double));
        SPARSE_BE_x_prev = calloc(tran_SPARSE_rows, sizeof(double));

        cs_gaxpy(tran_sparse_tilde_C, tran_SPARSE_X0, tran_sparse_right_hand_BE);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] + tran_SPARSE_e[i];
        }
    }

    if((iter == 0) && (k_value > 1)){
        free(SPARSE_BE_x_curr);
        SPARSE_BE_x_curr = calloc(tran_SPARSE_rows, sizeof(double));

        cs_gaxpy(tran_sparse_tilde_C, SPARSE_BE_x_prev, tran_sparse_right_hand_BE);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] + tran_SPARSE_e[i];
        }
    }

    // y = A * x
    for (j = 0 ; j < tran_SPARSE_rows ; j++){
        for (p = sparse_left_hand_BE -> p[j] ; p < sparse_left_hand_BE -> p[j+1] ; p++){
            tran_sparse_y_iterative[sparse_left_hand_BE -> i[p]] = tran_sparse_y_iterative[sparse_left_hand_BE -> i[p]] + sparse_left_hand_BE -> x[p] * SPARSE_BE_x_curr[j];
        }
    }

    // r1 = b - y
    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_sparse_r_iterative[i] = tran_sparse_right_hand_BE[i] - tran_sparse_y_iterative[i];
    }

    // tran_sparse_r_tilde_iterative = tran_sparse_r_iterative
    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_sparse_r_tilde_iterative[i] = tran_sparse_r_iterative[i];
    }

    tran_norm_r = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_r_iterative, 1);
    tran_norm_b = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_right_hand_BE, 1);

    if(tran_norm_b == 0){
        tran_norm_b = 1;
    }

    tran_sparse_inverse_matrix();

    tran_m_cmprsd_iterative = cs_spalloc(tran_SPARSE_rows, tran_SPARSE_rows, tran_SPARSE_rows, 1, 1);
    tran_m_cmprsd_iterative->nz = tran_SPARSE_rows;

    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_m_cmprsd_iterative -> i[i] = i;
        tran_m_cmprsd_iterative -> p[i] = i;
        tran_m_cmprsd_iterative -> x[i] = gsl_matrix_get(tran_M_inverse_iterative, i, i);
    }

    tran_m_cc_iterative = cs_compress(tran_m_cmprsd_iterative);
    cs_dupl(tran_m_cc_iterative);

    while(((tran_norm_r / tran_norm_b) > ITOL) && (tran_iter < tran_SPARSE_rows)){
        tran_iter = tran_iter + 1;

        tran_sparse_preconditioner_solve();
        tran_sparse_transpose_preconditioner();


        tran_rho = dot_product(tran_sparse_r_tilde_iterative, tran_sparse_z_iterative, tran_SPARSE_rows);

        if((fabs(tran_rho)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        if(tran_iter == 1){
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_sparse_p_iterative[i] = tran_sparse_z_iterative[i];

                tran_sparse_p_tilde_iterative[i] = tran_sparse_z_tilde_iterative[i];
            }
        }

        else{
            tran_beta = (tran_rho / tran_rho1);
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_sparse_p_iterative[i] = (tran_beta * tran_sparse_p_iterative[i]) + tran_sparse_z_iterative[i];

                tran_sparse_p_tilde_iterative[i] = (tran_beta * tran_sparse_p_tilde_iterative[i]) + tran_sparse_z_tilde_iterative[i];
            }

        }

        tran_rho1 = tran_rho;

        for(j = 0; j < tran_SPARSE_rows; j++){
            for (p = sparse_left_hand_BE -> p[j]; p < sparse_left_hand_BE -> p[j + 1] ; p++){
                tran_sparse_q_iterative[sparse_left_hand_BE -> i[p]] = tran_sparse_q_iterative[sparse_left_hand_BE -> i[p]]
                                                        + (sparse_left_hand_BE -> x[p] * tran_sparse_p_iterative[j]);
            }
        }

        for(j = 0; j < tran_SPARSE_rows; j++){
            tran_sparse_q_tilde_iterative[j] = 0;

            for (p = sparse_left_hand_BE -> p[j]; p < sparse_left_hand_BE -> p[j + 1] ; p++){
                tran_sparse_q_tilde_iterative[j] = tran_sparse_q_tilde_iterative[j]
                                                        + (sparse_left_hand_BE -> x[p] * tran_sparse_p_tilde_iterative[sparse_left_hand_BE->i[p]]);
            }
        }

        tran_omega = dot_product(tran_sparse_p_tilde_iterative, tran_sparse_q_iterative, tran_SPARSE_rows);

        if((fabs(tran_omega)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        tran_alpha = (tran_rho / tran_omega);

        cblas_daxpy(tran_SPARSE_rows, tran_alpha, tran_sparse_p_iterative, 1, SPARSE_BE_x_curr, 1);

        cblas_daxpy(tran_SPARSE_rows, -tran_alpha, tran_sparse_q_iterative, 1, tran_sparse_r_iterative, 1);

        cblas_daxpy(tran_SPARSE_rows, -tran_alpha, tran_sparse_q_tilde_iterative, 1, tran_sparse_r_tilde_iterative, 1);

        tran_norm_r = cblas_dnrm2(sparse_rowsA, tran_sparse_r_iterative, 1);

        // reset q_iterative vector to zero after each step
        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_q_iterative[i] = 0;
        }

    }


    for(i = 0; i < tran_SPARSE_rows; i++){
        SPARSE_BE_x_prev[i] = SPARSE_BE_x_curr[i];
    }

    gsl_permutation_free(tran_permut_iterative);

    gsl_matrix_free(tran_M_iterative);
    gsl_matrix_free(tran_M_inverse_iterative);

    cs_spfree(tran_m_cmprsd_iterative);
    cs_spfree(tran_m_cc_iterative);

    free(tran_sparse_y_iterative);
    free(tran_sparse_right_hand_BE);
    free(tran_sparse_m_iterative);
    free(tran_sparse_z_tilde_iterative);
    free(tran_sparse_z_iterative);
    free(tran_sparse_q_tilde_iterative);
    free(tran_sparse_q_iterative);
    free(tran_sparse_p_iterative);
    free(tran_sparse_p_tilde_iterative);
    free(tran_sparse_r_iterative);
    free(tran_sparse_r_tilde_iterative);
}

void tran_sparse_inverse_matrix(){
    int i, j;
    int signum;

    if(backwardeuler_flag == 1){
        for(i = 0; i < tran_SPARSE_rows; i++){
            for(j = 0; j < tran_SPARSE_rows; j++){
                if(i == j){
                    if(sparse_left_hand_BE->x[i] != 0){
                        gsl_matrix_set(tran_M_iterative, i, j, tran_sparse_m_iterative[i]);
                    }
                    else{
                        gsl_matrix_set(tran_M_iterative, i, j, 1);
                    }
                }
            }
        }
    }

    else if(trapezoidal_flag == 1){
        for(i = 0; i < tran_SPARSE_rows; i++){
            for(j = 0; j < tran_SPARSE_rows; j++){
                if(i == j){
                    if(sparse_left_hand_TR->x[i] != 0){
                        gsl_matrix_set(tran_M_iterative, i, j, tran_sparse_m_iterative[i]);
                    }
                    else{
                        gsl_matrix_set(tran_M_iterative, i, j, 1);
                    }
                }
            }
        }
    }

    gsl_linalg_LU_decomp(tran_M_iterative, tran_permut_iterative, &signum);

    gsl_linalg_LU_invert(tran_M_iterative, tran_permut_iterative, tran_M_inverse_iterative);
}

void tran_sparse_preconditioner_solve(){
    int p, j;

    //cs_gaxpy(sparse_m_cc_iterative, sparse_r_iterative, sparse_z_iterative);
    for(j = 0; j < tran_SPARSE_rows; j++){
        for(p = tran_m_cc_iterative -> p[j] ; p < tran_m_cc_iterative -> p[j+1] ; p++){
            tran_sparse_z_iterative[tran_m_cc_iterative -> i[p]] = tran_m_cc_iterative -> x[p] * tran_sparse_r_iterative[j];
        }
    }
}

void tran_sparse_transpose_preconditioner(){
    int p, j;

    for(j = 0; j < tran_SPARSE_rows; j++){
        tran_sparse_z_tilde_iterative[j] = 0;

        for(p = tran_m_cc_iterative -> p[j] ; p < tran_m_cc_iterative -> p[j + 1] ; p++){
            tran_sparse_z_tilde_iterative[j] = tran_m_cc_iterative -> x[p] * tran_sparse_r_iterative[tran_m_cc_iterative->i[p]];
            //sparse_z_tilde_iterative[j] = sparse_m_cc_iterative -> x[p] * sparse_r_iterative[sparse_m_cc_iterative->i[p]];
        }
    }

}

void tran_SPARSE_CG_sweep(int iter, int k_value){
    int i, j, p;

    int tran_iter = 0;

    double tran_norm_r = 0;
    double tran_norm_b = 0;

    double tran_rho = 0;
    double tran_rho1 = 0;
    double tran_beta = 0;
    double tran_alpha = 0;

    tran_sparse_y_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_y_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_y_iterative failed.\n");

        exit(8);
    }

    tran_sparse_r_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_r_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_r_iterative failed.\n");

        exit(8);
    }

    tran_sparse_p_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_p_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_p_iterative failed.\n");

        exit(8);
    }

    tran_sparse_z_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_z_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_z_iterative failed.\n");

        exit(8);
    }

    tran_sparse_q_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_q_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_q_iterative failed.\n");

        exit(8);
    }

    tran_sparse_m_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_m_iterative == NULL){
        printf("\nMemory allocation for tran_sparse_m_iterative failed.\n");

        exit(8);
    }
    
    for(j = 0; j < tran_SPARSE_rows; j++){
        for(p = sparse_left_hand_BE -> p[j]; p < sparse_left_hand_BE -> p[ j+ 1]; p++){
            if(sparse_left_hand_BE -> i[p] == j){
                tran_sparse_m_iterative[j] = sparse_left_hand_BE -> x[p];
            }
        }
    }

    tran_sparse_right_hand_BE = (double *)calloc(sparse_rowsA, sizeof(double));
    if(tran_sparse_right_hand_BE == NULL){
        printf("\nMemory allocation for tran_sparse_right_hand_BE failed.\n");

        exit(8);
    }

    tran_M_iterative = gsl_matrix_alloc(tran_SPARSE_rows, tran_SPARSE_rows); // check columns
    tran_M_inverse_iterative = gsl_matrix_alloc(tran_SPARSE_rows, tran_SPARSE_rows);

    for(i = 0; i < tran_SPARSE_rows; i++){
        for(j = 0; j < tran_SPARSE_columns; j++){
            gsl_matrix_set(tran_M_iterative, i, j, 0);
            gsl_matrix_set(tran_M_inverse_iterative, i, j, 0);
        }
    }

    tran_permut_iterative = gsl_permutation_alloc(tran_SPARSE_rows);

    if(iter == 1){
        SPARSE_BE_x_curr = calloc(tran_SPARSE_rows, sizeof(double));
        SPARSE_BE_x_prev = calloc(tran_SPARSE_rows, sizeof(double));

        cs_gaxpy(tran_sparse_tilde_C, tran_SPARSE_X0, tran_sparse_right_hand_BE);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] + tran_SPARSE_e[i];
        }
    }

    if((iter == 0) && (k_value > 1)){
        free(SPARSE_BE_x_curr);
        SPARSE_BE_x_curr = calloc(tran_SPARSE_rows, sizeof(double));

        cs_gaxpy(tran_sparse_tilde_C, SPARSE_BE_x_prev, tran_sparse_right_hand_BE);

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] * (1 / tran_time_step);
        }

        for(i = 0; i < tran_SPARSE_rows; i++){
            tran_sparse_right_hand_BE[i] = tran_sparse_right_hand_BE[i] + tran_SPARSE_e[i];
        }
    }

     // y = A * x
    for (j = 0 ; j < tran_SPARSE_rows ; j++){
        for (p = sparse_left_hand_BE -> p[j] ; p < sparse_left_hand_BE -> p[j+1] ; p++){
            tran_sparse_y_iterative[sparse_left_hand_BE -> i[p]] = tran_sparse_y_iterative[sparse_left_hand_BE -> i[p]] + sparse_left_hand_BE -> x[p] * SPARSE_BE_x_curr[j];
        }
    }

    // r1 = b - y
    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_sparse_r_iterative[i] = tran_sparse_right_hand_BE[i] - tran_sparse_y_iterative[i];
    }

    tran_norm_r = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_r_iterative, 1);
    tran_norm_b = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_right_hand_BE, 1);

    if(tran_norm_b == 0){
        tran_norm_b = 1;
    }

    tran_sparse_inverse_matrix();

    tran_m_cmprsd_iterative = cs_spalloc(tran_SPARSE_rows, tran_SPARSE_rows, tran_SPARSE_rows, 1, 1);
    tran_m_cmprsd_iterative->nz = tran_SPARSE_rows;

    for(i = 0; i < tran_SPARSE_rows; i++){
        tran_m_cmprsd_iterative -> i[i] = i;
        tran_m_cmprsd_iterative -> p[i] = i;
        tran_m_cmprsd_iterative -> x[i] = gsl_matrix_get(tran_M_inverse_iterative, i, i);
    }

    tran_m_cc_iterative = cs_compress(tran_m_cmprsd_iterative);
    cs_dupl(tran_m_cc_iterative);

    while(((tran_norm_r / tran_norm_b) > ITOL) && (tran_iter < tran_SPARSE_rows)){
        tran_iter = tran_iter + 1;

        tran_sparse_preconditioner_solve();

        tran_rho = dot_product(tran_sparse_r_iterative, tran_sparse_z_iterative, tran_SPARSE_rows);

        if(tran_iter == 1){
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_sparse_p_iterative[i] = tran_sparse_z_iterative[i];
            }
        }

        else{
            tran_beta = (tran_rho / tran_rho1);

            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_sparse_p_iterative[i] = (tran_beta * tran_sparse_p_iterative[i]) + tran_sparse_z_iterative[i];
            }
        }

        tran_rho1 = tran_rho;

        for(j = 0; j < tran_SPARSE_rows; j++){
            for (p = sparse_left_hand_BE -> p[j]; p < sparse_left_hand_BE -> p[j+1] ; p++){
                tran_sparse_q_iterative[sparse_left_hand_BE -> i[p]] = tran_sparse_q_iterative[sparse_left_hand_BE -> i[p]]
                                                        + (sparse_left_hand_BE -> x[p] * tran_sparse_p_iterative[j]);
            }
        }

        tran_alpha = (tran_rho / dot_product(tran_sparse_p_iterative, tran_sparse_q_iterative, tran_SPARSE_rows));

        cblas_daxpy(tran_SPARSE_rows, tran_alpha, tran_sparse_p_iterative, 1, SPARSE_BE_x_curr, 1);

        cblas_daxpy(tran_SPARSE_rows, -tran_alpha, tran_sparse_q_iterative, 1, tran_sparse_r_iterative, 1);

        tran_norm_r = cblas_dnrm2(tran_SPARSE_rows, tran_sparse_r_iterative, 1);

        // reset tran_sparse_q_iterative vector to zero after each step
        for(i = 0; i < sparse_rowsA; i++){
            tran_sparse_q_iterative[i] = 0;
        }
    }

    for(i = 0; i < tran_SPARSE_rows; i++){
        SPARSE_BE_x_prev[i] = SPARSE_BE_x_curr[i];
    }

    gsl_permutation_free(tran_permut_iterative);

    gsl_matrix_free(tran_M_iterative);
    gsl_matrix_free(tran_M_inverse_iterative);

    cs_spfree(tran_m_cmprsd_iterative);
    cs_spfree(tran_m_cc_iterative);

    free(tran_sparse_right_hand_BE);
    free(tran_sparse_y_iterative);
    free(tran_sparse_r_iterative);
    free(tran_sparse_p_iterative);
    free(tran_sparse_m_iterative);
    free(tran_sparse_z_iterative);
    free(tran_sparse_q_iterative);
}

void tran_LU_sweep(int iter, int k_value){
    int i, s;

    // first iteration (k = 1) -> x(0)
    // right_hand_BE * BE_x_curr = tran_MNA_e + (((1 / h) * tran_MNA_C) * BE_x_prev)
    if(iter == 1){
        gsl_blas_dgemv(CblasNoTrans, 1.0, tran_MNA_C, tran_MNA_X0, 0.0, right_hand_BE);
        gsl_vector_add(right_hand_BE, tran_MNA_e);
    }
    
    // remaining iterations (k = 2, 3, ...) -> x(tk)
    if((iter == 0) && (k_value > 1)){
        gsl_blas_dgemv(CblasNoTrans, 1.0, tran_MNA_C, BE_x_prev, 0.0, right_hand_BE);
        gsl_vector_add(right_hand_BE, tran_MNA_e);
    }

    for(i = 0; i < tran_MNA_rows; i++){
        gsl_vector_set(BE_x_curr, i, 0);
    }

    permutation_BE = gsl_permutation_alloc(tran_MNA_rows);

    alt_left_hand_BE = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    gsl_matrix_memcpy(alt_left_hand_BE, left_hand_BE);

    gsl_linalg_LU_decomp(alt_left_hand_BE, permutation_BE, &s);

    gsl_linalg_LU_solve(alt_left_hand_BE, permutation_BE, right_hand_BE, BE_x_curr);

    /*printf("\nX_curr:\n");
    for(i = 0; i < tran_MNA_rows; i++){
        printf("%lf\n", gsl_vector_get(BE_x_curr, i));
    }*/

    gsl_vector_memcpy(BE_x_prev, BE_x_curr);

    gsl_permutation_free(permutation_BE);
    gsl_matrix_free(alt_left_hand_BE);
}

void tran_Cholesky_sweep(int iter, int k_value){
    int i;

    // first iteration (k = 1) -> x(0)
    if(iter == 1){
        gsl_blas_dgemv(CblasNoTrans, 1.0, tran_MNA_C, tran_MNA_X0, 0.0, right_hand_BE);
        gsl_vector_add(right_hand_BE, tran_MNA_e);
    }

    if((iter == 0) && (k_value > 1)){
        gsl_blas_dgemv(CblasNoTrans, 1.0, tran_MNA_C, BE_x_prev, 0.0, right_hand_BE);
        gsl_vector_add(right_hand_BE, tran_MNA_e);
    }

    for(i = 0; i < tran_MNA_rows; i++){
        gsl_vector_set(BE_x_curr, i, 0);
    }

    alt_left_hand_BE = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    gsl_matrix_memcpy(alt_left_hand_BE, left_hand_BE);

    gsl_linalg_cholesky_decomp(alt_left_hand_BE);
    gsl_linalg_cholesky_solve(alt_left_hand_BE, right_hand_BE, BE_x_curr);

    gsl_vector_memcpy(BE_x_prev, BE_x_curr);

    gsl_matrix_free(alt_left_hand_BE);
}

void tran_BiCG_sweep(int iter, int k_value){
    int i, j;

    int tran_iter = 0;

    double tran_alpha = 0;
    double tran_beta = 0;

    double tran_norm_r = 0;
    double tran_norm_b = 0;

    double tran_rho = 0;
    double tran_rho1 = 0;

    double tran_omega = 0;

    tran_r_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_r_tilde_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_p_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_p_tilde_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_p_scale_beta = gsl_vector_alloc(tran_MNA_rows);
    tran_p_tilde_scale_beta = gsl_vector_alloc(tran_MNA_rows);

    tran_q_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_q_tilde_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_z_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_z_tilde_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_Ax_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_bminusAx_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_M_iterative = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    tran_M_inverse_iterative = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);

    tran_permut_iterative = gsl_permutation_alloc(tran_MNA_rows);

    for(i = 0; i < tran_MNA_rows; i++){
        for(j = 0; j < tran_MNA_columns; j++){
            gsl_matrix_set(tran_M_iterative, i, j, 0);

            gsl_matrix_set(tran_M_inverse_iterative, i, j, 0);
        }

        gsl_vector_set(tran_r_iterative, i, 0);
        gsl_vector_set(tran_r_tilde_iterative, i, 0);

        gsl_vector_set(tran_p_iterative, i , 0);
        gsl_vector_set(tran_p_tilde_iterative, i, 0);
        gsl_vector_set(tran_p_scale_beta, i, 0);
        gsl_vector_set(tran_p_tilde_scale_beta, i, 0);

        gsl_vector_set(tran_q_iterative, i, 0);
        gsl_vector_set(tran_q_tilde_iterative, i, 0);

        gsl_vector_set(tran_z_iterative, i, 0);
        gsl_vector_set(tran_z_tilde_iterative, i, 0);

        gsl_vector_set(tran_Ax_iterative, i, 0);
        gsl_vector_set(tran_bminusAx_iterative, i, 0);
    }

    if(iter == 1){
        // b = right_hand_BE
        gsl_blas_dgemv(CblasNoTrans, 1.0, tran_MNA_C, tran_MNA_X0, 0.0, right_hand_BE);
        gsl_vector_add(right_hand_BE, tran_MNA_e);
    }

    if((iter == 0) && (k_value > 1)){
        gsl_blas_dgemv(CblasNoTrans, 1.0, tran_MNA_C, BE_x_prev, 0.0, right_hand_BE);
        gsl_vector_add(right_hand_BE, tran_MNA_e);
    }

    // initial guess x(0)
    for(i = 0; i < tran_MNA_rows; i++){
        gsl_vector_set(BE_x_curr, i, 0);
    }

    // A = alt_left_hand_BE
    alt_left_hand_BE = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    gsl_matrix_memcpy(alt_left_hand_BE, left_hand_BE);

    // set tran_M_iterative
    for(i = 0; i < tran_MNA_rows; i++){
        for(j = 0; j < tran_MNA_columns; j++){
            if(i == j){
                if((gsl_matrix_get(alt_left_hand_BE, i, j)) != 0){
                    gsl_matrix_set(tran_M_iterative, i, j, (gsl_matrix_get(alt_left_hand_BE, i, j)));
                }
                else{
                    gsl_matrix_set(tran_M_iterative, i, j, 1);
                }
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1.0, alt_left_hand_BE, BE_x_curr, 0.0, tran_Ax_iterative);

    gsl_vector_memcpy(tran_bminusAx_iterative, right_hand_BE);

    gsl_vector_sub(tran_bminusAx_iterative, tran_Ax_iterative);

    gsl_vector_memcpy(tran_r_iterative, tran_bminusAx_iterative);

    gsl_vector_memcpy(tran_r_tilde_iterative, tran_bminusAx_iterative);

    tran_norm_r = gsl_blas_dnrm2(tran_r_iterative);
    tran_norm_b = gsl_blas_dnrm2(right_hand_BE);

    if(tran_norm_b == 0){
        tran_norm_b = 1;
    }

    tran_inverse_matrix();

    while(((tran_norm_r / tran_norm_b) > ITOL) && (tran_iter < (curr_pos_oldname_table + G2_components))){
        tran_iter = tran_iter + 1;

        tran_preconditioner_solve();
        tran_transpose_preconditioner_solve();

        gsl_blas_ddot(tran_r_tilde_iterative, tran_z_iterative, &tran_rho);

        if((fabs(tran_rho)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        if(tran_iter == 1){
            gsl_vector_memcpy(tran_p_iterative, tran_z_iterative);

            gsl_vector_memcpy(tran_p_tilde_iterative, tran_z_tilde_iterative);
        }
        else{
            tran_beta = (tran_rho / tran_rho1);

            gsl_vector_memcpy(tran_p_scale_beta, tran_p_iterative);
            gsl_vector_scale(tran_p_scale_beta, tran_beta);
            gsl_vector_add(tran_p_scale_beta, tran_z_iterative);
            gsl_vector_memcpy(tran_p_iterative, tran_p_scale_beta);
            
            gsl_vector_memcpy(tran_p_tilde_scale_beta, tran_p_tilde_iterative);
            gsl_vector_scale(tran_p_tilde_scale_beta, tran_beta);
            gsl_vector_add(tran_p_tilde_scale_beta, tran_z_tilde_iterative);
            gsl_vector_memcpy(tran_p_tilde_iterative, tran_p_tilde_scale_beta);
        }

        tran_rho1 = tran_rho;

        gsl_blas_dgemv(CblasNoTrans, 1.0, alt_left_hand_BE, tran_p_iterative, 0.0, tran_q_iterative);
        gsl_blas_dgemv(CblasTrans, 1.0, alt_left_hand_BE, tran_p_tilde_iterative, 0.0, tran_q_tilde_iterative);

        gsl_blas_ddot(tran_p_tilde_iterative, tran_q_iterative, &tran_omega);

        if((fabs(tran_omega)) < 1e-14){
            exit(10);
        }

        tran_alpha = (tran_rho / tran_omega);

        gsl_blas_daxpy(tran_alpha, tran_p_iterative, BE_x_curr);

        gsl_blas_daxpy(-tran_alpha, tran_q_iterative, tran_r_iterative);

        gsl_blas_daxpy(-tran_alpha, tran_q_tilde_iterative, tran_r_tilde_iterative);

        tran_norm_r = gsl_blas_dnrm2(tran_r_iterative);
    }

    gsl_vector_memcpy(BE_x_prev, BE_x_curr);

    gsl_matrix_free(alt_left_hand_BE);
        
    gsl_matrix_free(tran_M_iterative);
    gsl_matrix_free(tran_M_inverse_iterative);

    gsl_vector_free(tran_r_iterative);
    gsl_vector_free(tran_r_tilde_iterative);
    gsl_vector_free(tran_p_iterative);
    gsl_vector_free(tran_p_tilde_iterative);
    gsl_vector_free(tran_p_scale_beta);
    gsl_vector_free(tran_p_tilde_scale_beta);
    gsl_vector_free(tran_q_iterative);
    gsl_vector_free(tran_q_tilde_iterative);
    gsl_vector_free(tran_z_iterative);
    gsl_vector_free(tran_z_tilde_iterative);
    gsl_vector_free(tran_Ax_iterative);
    gsl_vector_free(tran_bminusAx_iterative);

    gsl_permutation_free(tran_permut_iterative);
}

void tran_transpose_preconditioner_solve(){
    gsl_blas_dgemv(CblasTrans, 1.0, tran_M_inverse_iterative, tran_r_tilde_iterative, 0.0, tran_z_tilde_iterative);
}

void tran_CG_sweep(int iter, int k_value){
    int i, j;

    int tran_iter = 0;

    double tran_alpha = 0;
    double tran_beta = 0;

    double tran_norm_r = 0;
    double tran_norm_b = 0;

    double tran_rho = 0;
    double tran_rho1 = 0;

    double tran_ptrans_q = 0;

    tran_r_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_p_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_p_scale_beta = gsl_vector_alloc(tran_MNA_rows);

    tran_q_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_z_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_Ax_iterative = gsl_vector_alloc(tran_MNA_rows);
    tran_bminusAx_iterative = gsl_vector_alloc(tran_MNA_rows);

    tran_M_iterative = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    tran_M_inverse_iterative = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);

    tran_permut_iterative = gsl_permutation_alloc(tran_MNA_rows);

    for(i = 0; i < tran_MNA_rows; i++){
        for(j = 0; j < tran_MNA_columns; j++){
            gsl_matrix_set(tran_M_iterative, i, j, 0);

            gsl_matrix_set(tran_M_inverse_iterative, i, j, 0);
        }

        gsl_vector_set(tran_r_iterative, i, 0);

        gsl_vector_set(tran_p_iterative, i , 0);
        gsl_vector_set(tran_p_scale_beta, i, 0);

        gsl_vector_set(tran_q_iterative, i, 0);

        gsl_vector_set(tran_z_iterative, i, 0);

        gsl_vector_set(tran_Ax_iterative, i, 0);
        gsl_vector_set(tran_bminusAx_iterative, i, 0);
    }

    if(iter == 1){
        // b = right_hand_BE
        gsl_blas_dgemv(CblasNoTrans, 1.0, tran_MNA_C, tran_MNA_X0, 0.0, right_hand_BE);
        gsl_vector_add(right_hand_BE, tran_MNA_e);
    }

    if((iter == 0) && (k_value > 1)){
        // b = right_hand_BE
        gsl_blas_dgemv(CblasNoTrans, 1.0, tran_MNA_C, BE_x_prev, 0.0, right_hand_BE);
        gsl_vector_add(right_hand_BE, tran_MNA_e);
    }

    // initial guess x(0)
    for(i = 0; i < tran_MNA_rows; i++){
        gsl_vector_set(BE_x_curr, i, 0);
    }
    
    // A = alt_left_hand_BE
    alt_left_hand_BE = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
    gsl_matrix_memcpy(alt_left_hand_BE, left_hand_BE);

    // set tran_M_iterative
    for(i = 0; i < tran_MNA_rows; i++){
        for(j = 0; j < tran_MNA_columns; j++){
            if(i == j){
                if((gsl_matrix_get(alt_left_hand_BE, i, j)) != 0){
                    gsl_matrix_set(tran_M_iterative, i, j, (gsl_matrix_get(alt_left_hand_BE, i, j)));
                }
                else{
                    gsl_matrix_set(tran_M_iterative, i, j, 1);
                }
            }
        }
    }

    // Ax = tran_Ax_iterative = (left_hand_BE * BE_x_curr)
    gsl_blas_dgemv(CblasNoTrans, 1.0, alt_left_hand_BE, BE_x_curr, 0.0, tran_Ax_iterative);

    // b - Ax = tran_bminusAx_iterative - tran_Ax_iterative
    gsl_vector_memcpy(tran_bminusAx_iterative, right_hand_BE);
    gsl_vector_sub(tran_bminusAx_iterative, tran_Ax_iterative);
        
    // r = tran_r_iterative = b - Ax
    gsl_vector_memcpy(tran_r_iterative, tran_bminusAx_iterative);

    tran_norm_r = gsl_blas_dnrm2(tran_r_iterative);
    tran_norm_b = gsl_blas_dnrm2(right_hand_BE);

    if(tran_norm_b == 0){    // if right_hand_BE = 0
        tran_norm_b = 1;
    }

    // inverse tran_M_iterative
    tran_inverse_matrix();

    while(((tran_norm_r / tran_norm_b) > ITOL) && (tran_iter < (curr_pos_oldname_table + G2_components))){
        tran_iter = tran_iter + 1;

        // tran_M_inverse_iterative * tran_z_iterative = tran_r_iterative
        tran_preconditioner_solve();

        // rho = tran_r_iterative * tran_z_iterative
        gsl_blas_ddot(tran_r_iterative, tran_z_iterative, &tran_rho);

        if(tran_iter == 1){
            // tran_p_iterative = tran_z_iterative
            gsl_vector_memcpy(tran_p_iterative, tran_z_iterative);
        }

        else{
            tran_beta = (tran_rho / tran_rho1);

            // p = tran_p_iterative = z + (b * p) 
            gsl_vector_memcpy(tran_p_scale_beta, tran_p_iterative);
            gsl_vector_scale(tran_p_scale_beta, tran_beta);
            gsl_vector_add(tran_p_scale_beta, tran_z_iterative);    
            gsl_vector_memcpy(tran_p_iterative, tran_p_scale_beta);
        }

        tran_rho1 = tran_rho;

        // q = tran_q_iterative = A * p = left_hand_BE * tran_p_iterative
        gsl_blas_dgemv(CblasNoTrans, 1.0, left_hand_BE, tran_p_iterative, 0.0, tran_q_iterative);

        gsl_blas_ddot(tran_p_iterative, tran_q_iterative, &tran_ptrans_q);
        tran_alpha = (tran_rho / tran_ptrans_q);

        // x = BE_x_curr = x + (alpha * p)
        gsl_blas_daxpy(tran_alpha, tran_p_iterative, BE_x_curr);

        // r = tran_r_iterative = r - (alpha * q)
        gsl_blas_daxpy(-tran_alpha, tran_q_iterative, tran_r_iterative);

        // calculate next step norm
        tran_norm_r = gsl_blas_dnrm2(tran_r_iterative);
    }

    gsl_vector_memcpy(BE_x_prev, BE_x_curr);

    gsl_matrix_free(alt_left_hand_BE);
        
    gsl_matrix_free(tran_M_iterative);
    gsl_matrix_free(tran_M_inverse_iterative);

    gsl_vector_free(tran_r_iterative);
    gsl_vector_free(tran_p_iterative);
    gsl_vector_free(tran_p_scale_beta);
    gsl_vector_free(tran_q_iterative);
    gsl_vector_free(tran_z_iterative);
    gsl_vector_free(tran_Ax_iterative);
    gsl_vector_free(tran_bminusAx_iterative);

    gsl_permutation_free(tran_permut_iterative);
}

void tran_inverse_matrix(){
    int s;

    gsl_linalg_LU_decomp(tran_M_iterative, tran_permut_iterative, &s);

    gsl_linalg_LU_invert(tran_M_iterative, tran_permut_iterative, tran_M_inverse_iterative);
}

void tran_preconditioner_solve(){
    gsl_blas_dgemv(CblasNoTrans, 1.0, tran_M_inverse_iterative, tran_r_iterative, 0.0, tran_z_iterative);
}

void SPARSE_TRAN_sweep_voltagesource(double TRAN_value){
    int k_tran_expression;

    k_tran_expression = curr_pos_oldname_table + k_tran_variable;

    tran_SPARSE_e[k_tran_expression] = tran_SPARSE_e[k_tran_expression] + TRAN_value;

    k_tran_variable = k_tran_variable + 1;
}

void TRAN_sweep_voltagesource(double TRAN_value){
    int k_tran_expression;
    double tran_sweep_result;
    double updated_value;

    k_tran_expression = curr_pos_oldname_table + k_tran_variable;

    //tran_MNA_e[k_tran_expression] = tran_MNA_e[k_tran_expression] + TRAN_value;

    tran_sweep_result = gsl_vector_get(tran_MNA_e, k_tran_expression);

    updated_value = tran_sweep_result + TRAN_value;

    gsl_vector_set(tran_MNA_e, k_tran_expression, updated_value);

    k_tran_variable = k_tran_variable + 1;
}

void SPARSE_TRAN_sweep_currentsource(int TRAN_pos_node, int TRAN_neg_node, double TRAN_value){
    int pos_position, neg_position;

    pos_position = TRAN_pos_node - 1;
    neg_position = TRAN_neg_node - 1;

    if(TRAN_pos_node == 0){
        tran_SPARSE_e[neg_position] = tran_SPARSE_e[neg_position] + TRAN_value;

        return;
    }

    if(TRAN_neg_node == 0){
        tran_SPARSE_e[pos_position] = tran_SPARSE_e[pos_position] - TRAN_value;

        return;
    }

    tran_SPARSE_e[pos_position] = tran_SPARSE_e[pos_position] - TRAN_value;

    tran_SPARSE_e[neg_position] = tran_SPARSE_e[neg_position] + TRAN_value;
}

void TRAN_sweep_currentsource(int TRAN_pos_node, int TRAN_neg_node, double TRAN_value){
    int pos_position, neg_position;
    double tran_sweep_result;
    double updated_value;
    double alt_tran_sweep_result;
    double alt_updated_value;

    pos_position = TRAN_pos_node - 1;
    neg_position = TRAN_neg_node - 1;

    if(TRAN_pos_node == 0){
       // array_b[neg_position] = array_b[neg_position] + MNA_value;
        tran_sweep_result = gsl_vector_get(tran_MNA_e, neg_position);
        updated_value = tran_sweep_result + TRAN_value;

        gsl_vector_set(tran_MNA_e, neg_position, updated_value);

        return;
    }

    if(TRAN_neg_node == 0){
        //array_b[pos_position] = array_b[pos_position] - MNA_value;
        tran_sweep_result = gsl_vector_get(tran_MNA_e, pos_position);
        updated_value = tran_sweep_result - TRAN_value;

        gsl_vector_set(tran_MNA_e, pos_position, updated_value);

        return;
    }

    //array_b[pos_position] = array_b[pos_position] - MNA_value;
    tran_sweep_result = gsl_vector_get(tran_MNA_e, pos_position);
    updated_value = tran_sweep_result - TRAN_value;

    gsl_vector_set(tran_MNA_e, pos_position, updated_value);

    //array_b[neg_position] = array_b[neg_position] + MNA_value;
    alt_tran_sweep_result = gsl_vector_get(tran_MNA_e, neg_position);
    alt_updated_value = alt_tran_sweep_result + TRAN_value;

    gsl_vector_set(tran_MNA_e, neg_position, alt_updated_value);
}

void init_TR_hands(){
    if(sparse_flag == 1){
        //sparse_left_hand_BE = cs_spalloc(tran_SPARSE_rows, tran_SPARSE_columns, non_zeros, 1, 1);
        //sparse_left_hand_BE->nz = non_zeros;

        //sparse_right_hand_BE = cs_spalloc(tran_SPARSE_rows, tran_SPARSE_columns, non_zeros_tran_C, 1, 1);
        //sparse_right_hand_BE->nz = non_zeros_tran_C;
        tran_SPARSE_e_prev = NULL;

        tran_SPARSE_e_prev = (double *)calloc(tran_SPARSE_rows, sizeof(double));
        
        if(tran_SPARSE_e_prev == NULL){
            printf("\nMemory allocation for TRANSIENT SPARSE MATRIX failed.\n");

            exit(11);
        }
    }

    else{
        left_hand_TR = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
        right_hand_TR = gsl_vector_alloc(tran_MNA_rows);

        TR_x_curr = gsl_vector_alloc(tran_MNA_rows);
        TR_x_prev = gsl_vector_alloc(tran_MNA_rows);

        tran_MNA_e_prev = gsl_vector_alloc(tran_MNA_rows);

        right_hand_TR_sub = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);

        alt_tran_MNA_G = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);

        //permutation_BE = gsl_permutation_alloc(tran_MNA_rows);
    }
}

void init_BE_hands(){
    if(sparse_flag == 1){
        //sparse_left_hand_BE = cs_spalloc(tran_SPARSE_rows, tran_SPARSE_columns, non_zeros, 1, 1);
        //sparse_left_hand_BE->nz = non_zeros;

        //sparse_right_hand_BE = cs_spalloc(tran_SPARSE_rows, tran_SPARSE_columns, non_zeros_tran_C, 1, 1);
        //sparse_right_hand_BE->nz = non_zeros_tran_C;
    }

    else{
        left_hand_BE = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);
        right_hand_BE = gsl_vector_alloc(tran_MNA_rows);

        BE_x_curr = gsl_vector_alloc(tran_MNA_rows);
        BE_x_prev = gsl_vector_alloc(tran_MNA_rows);

        //permutation_BE = gsl_permutation_alloc(tran_MNA_rows);
    }

}

double EXP_lookup(double tcurr, double i1, double i2, double td1, double td2, double tc1, double tc2, double tfinal){
    double exp_result = 0;

    double exponential_part = 0;
    double exponential_result = 0;

    double alt_exponential_part = 0;
    double alt_exponential_result = 0;

    if((tcurr >= 0) && (tcurr <= td1)){
        exp_result = i1;

        return(exp_result);
    }

    if((tcurr >= td1) && (tcurr <= td2)){
        exponential_part = -((tcurr - td1) / tc1);

        exponential_result = exp(exponential_part);

        exp_result = i1 + ((i2 - i1) * (1 - exponential_result));

        return(exp_result);
    }

    if((tcurr >= td2) && (tcurr <= tfinal)){
        exponential_part = -((tcurr - td2) / tc2);

        exponential_result = exp(exponential_part);

        alt_exponential_part = -((tcurr - td1) / tc1);

        alt_exponential_result = exp(alt_exponential_part);

        exp_result = i1 + ((i2 - i1) * (exponential_result - alt_exponential_result));

        return(exp_result);
    }

    return(-1);
}

double SIN_lookup(double tcurr, double i1, double ia, double fr, double td, double df, double ph, double tfinal){
    double sin_result = 0;

    double exponential_part = 0;
    double exponential_result = 0;

    double sine_part = 0;
    double sine_result = 0;

    double alt_sine_part = 0;
    //double alt_sine_result = 0;

    if((tcurr >= 0) && (tcurr <= td)){
        sine_part = ((2 * M_PI * ph) / (360));

        sine_result = sin(sine_part);

        sin_result = i1 + (ia * sine_result);

        return(sin_result);
    }

    if((tcurr >= td) && (tcurr <= tfinal)){
        sine_part = (2 * M_PI * fr * (tcurr - td));

        alt_sine_part = ((2 * M_PI * ph) / (360));

        sine_result = sin(sine_part + alt_sine_part);

        exponential_part = -((tcurr - td) * df);

        exponential_result = exp(exponential_part);

        sin_result = i1 + (ia * sine_result * exponential_result);

        return(sin_result);
    }

    return(-1);
}

double PULSE_lookup(double tcurr, double i1, double i2, double td, double tr, double tf, double pw, double per){
    double pulse_result = 0;

    double pulse_part = 0;

    if((tcurr >= 0) && (tcurr <= td)){
        pulse_result = i1;
        
        return(pulse_result);
    }

    if((tcurr >= (td + (k_pulse_variable * per))) && (tcurr <= (td + (tr + (k_pulse_variable * per))))){
        pulse_part = ((i2 - i1) / (tr));

        pulse_result = i1 + (pulse_part * (tcurr - (td + (k_pulse_variable * per))));

        //k_pulse_variable++;

        return(pulse_result);
    }

    if((tcurr >= (td + (tr + (k_pulse_variable * per)))) && (tcurr <= (td + tr + (pw + (k_pulse_variable * per))))){
        pulse_result = i2;

        //k_pulse_variable++;

        return(pulse_result);
    }

    if((tcurr >= (td + tr + (pw + (k_pulse_variable * per)))) && (tcurr <= (td + tr + pw + (tf + (k_pulse_variable * per))))){
        pulse_part = ((i1 - i2) / (tf));

        pulse_result = i2 + (pulse_part * (tcurr - td - tr - (pw + (k_pulse_variable * per))));

        //k_pulse_variable++;
        return(pulse_result);
    }

    if((tcurr >= (td + tr + pw + (tf + (k_pulse_variable * per)))) && (tcurr <= (td + (per + (k_pulse_variable * per))))){
        pulse_result = i1;

        if(tcurr == td + (per + (k_pulse_variable * per))){
            //pulse_result = i1;

            k_pulse_variable++;
        }

        return(pulse_result);
    }

    return(-1);
}

double PWL_lookup(double tcurr, double tfinal, int pwl_cntr, ptr_comp pwl_ptr){
    double pwl_result = 0;

    double pwl_part = 0;
    double alt_pwl_part = 0;
    
    int i;

    for(i = 0; i <= pwl_cntr; i++){
        if(tcurr == (pwl_ptr->pwl_field[i].pwl_start)){
            pwl_result = pwl_ptr->pwl_field[i].pwl_end;

            return(pwl_result);
        }

        if((tcurr >= 0) && (tcurr <= (pwl_ptr->pwl_field[0].pwl_start))){
            if(pwl_ptr->pwl_field[0].pwl_start > 0){
                pwl_result = pwl_ptr->pwl_field[0].pwl_end;

                return(pwl_result);
            }
        }

        if((tcurr >= (pwl_ptr->pwl_field[pwl_cntr].pwl_start)) && (tcurr <= tfinal)){
            if((pwl_ptr->pwl_field[pwl_cntr].pwl_start) < tfinal){
                pwl_result = pwl_ptr->pwl_field[pwl_cntr].pwl_end;

                return(pwl_result);
            }
        }

        if((tcurr >= (pwl_ptr->pwl_field[i].pwl_start)) && (tcurr <= (pwl_ptr->pwl_field[i + 1].pwl_start))){
            pwl_part = ((pwl_ptr->pwl_field[i + 1].pwl_end) - (pwl_ptr->pwl_field[i].pwl_end)) / 
                        ((pwl_ptr->pwl_field[i + 1].pwl_start) - (pwl_ptr->pwl_field[i].pwl_start));

            alt_pwl_part = (tcurr - (pwl_ptr->pwl_field[i].pwl_start));

            pwl_result = (pwl_ptr->pwl_field[i].pwl_end) + (pwl_part * alt_pwl_part);

            return(pwl_result);
        }
    }

    return(-1);
}

void free_BE_hands(){
    gsl_matrix_free(left_hand_BE);
    gsl_vector_free(right_hand_BE);
    gsl_vector_free(BE_x_curr);
    gsl_vector_free(BE_x_prev);

    free(SPARSE_BE_x_curr);
    free(SPARSE_BE_x_prev);
}

void free_TR_hands(){
    gsl_matrix_free(left_hand_TR);
    gsl_vector_free(right_hand_TR);
    gsl_matrix_free(right_hand_TR_sub);
    gsl_vector_free(TR_x_curr);
    gsl_vector_free(TR_x_prev);
    gsl_vector_free(tran_MNA_e_prev);
    gsl_matrix_free(alt_tran_MNA_G);

    free(SPARSE_TR_x_curr);
    free(SPARSE_TR_x_prev);
    free(tran_SPARSE_e_prev);
}