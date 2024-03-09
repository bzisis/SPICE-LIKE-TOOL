#include "AC_TRAN_sweep.h"
#include "TRAN_sweep.h"
#include "TRAN_DC_system_creation.h"
#include "SPARSE_DC_system_creation.h"
#include "MNA_DC_system_creation.h"
#include "netlist_parser.h"
#include "MNA_DC_system_direct_methods.h"
#include "MNA_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_direct_methods.h"

void AC_voltagesource_RH(gsl_complex, ptr_comp);
void AC_currentsource_RH(int, int, gsl_complex, ptr_comp);

void init_AC_matrices(){
    AC_rowsA = curr_pos_oldname_table + G2_components;
    AC_columnsA = curr_pos_oldname_table + G2_components;

    //AC_left_hand = gsl_matrix_complex_calloc(AC_rowsA, AC_columnsA);

    //AC_x = gsl_vector_complex_calloc(AC_rowsA);

    AC_right_hand = gsl_vector_complex_alloc(AC_rowsA);

    //AC_permutation = gsl_permutation_alloc(AC_rowsA);
}

void choose_AC_method(){
    init_AC_matrices();

    if(sparse_flag == 0){
        if((spd_flag == 1) && (iter_flag == 1)){
            exit(13);
        }
        if((spd_flag == 0) && (iter_flag == 1)){
            // BiCG
            AC_BiCG_sweep();
        }
        else{
            // LU
            AC_LU_sweep();
        }
    }

    else if(sparse_flag == 1){
        if((spd_flag == 1) && (iter_flag == 1)){
            exit(13);
        }
        if((spd_flag == 0) && (iter_flag == 1)){
            // SPARSE BiCG
            // SPARSE_AC_BiCG_sweep();
        }
        else{
            // SPARSE LU
            // SPARSE_AC_LU_sweep();
        }
    }
}

void AC_BiCG_sweep(){
    double w;

    int s;

    k_RH = 0;
    //k_LH = 0;

    int AC_step;

    ptr_comp AC_sweep_ptr;

    AC_sweep_ptr = list_head;

    char *AC_mag_plot_ptr = NULL;
    char *AC_phase_plot_ptr = NULL;

    int AC_mag_plot_node_name = 0;
    int AC_phase_plot_node_name = 0;

    double AC_mag_plot_x = 0;
    double AC_mag_plot_y = 0;

    double AC_phase_plot_x = 0;
    double AC_phase_plot_y = 0;

    gsl_complex AC_mag_node_result;
    gsl_complex AC_phase_node_result;

    double AC_mag_result = 0;
    double AC_phase_result = 0;

    FILE *AC_mag_fileptr = NULL;
    FILE *AC_mag_gnuplotpipe = NULL;
    char *AC_mag_gnuplotcommands = {"plot 'AC_mag.tmp' with linespoints"};

    FILE *AC_phase_fileptr = NULL;
    FILE *AC_phase_gnuplotpipe = NULL;
    char *AC_phase_gnuplotcommands = {"plot 'AC_phase.tmp' with linespoints"};

    AC_mag_fileptr = fopen("AC_mag.tmp", "w");
    AC_mag_gnuplotpipe = popen("gnuplot -persistent", "w");

    AC_phase_fileptr = fopen("AC_phase.tmp", "w");
    AC_phase_gnuplotpipe = popen("gnuplot -persistent", "w");

    if(plot_node != NULL){
        AC_mag_plot_ptr = strdup(plot_node);
        AC_phase_plot_ptr = strdup(plot_node);
    }
    
    if(print_node != NULL){
        AC_mag_plot_ptr = strdup(print_node);
        AC_phase_plot_ptr = strdup(print_node);
    }

    if((AC_mag_plot_ptr != NULL) && (AC_phase_plot_ptr != NULL)){
        AC_mag_plot_node_name = return_newpos_node_name(AC_mag_plot_ptr);

        AC_phase_plot_node_name = return_newpos_node_name(AC_phase_plot_ptr);
        
        free(AC_mag_plot_ptr);
        free(AC_phase_plot_ptr);
    }

    // fill up AC_right_hand
    while(AC_sweep_ptr != NULL){
        switch(AC_sweep_ptr -> comp_type){
        case 'V':
        case 'v':
            AC_voltagesource_RH(AC_sweep_ptr -> AC_complex_number, AC_sweep_ptr);
            break;

        case 'I':
        case 'i':
            AC_currentsource_RH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                AC_sweep_ptr -> AC_complex_number, AC_sweep_ptr);
            break;

        case 'L':
        case 'l':
            k_RH = k_RH + 1;
            break;
        
        default:
            break;
        }

        AC_sweep_ptr = AC_sweep_ptr -> nxt_comp;
    }

    if(AC_LIN_flag == 1){
        AC_step = ((AC_end_freq - AC_start_freq) / AC_points);

        for(w = AC_start_freq; w <= AC_end_freq; w = (w + AC_step)){
            // reset AC_sweep_ptr
            AC_sweep_ptr = list_head;

            k_LH = 0;

            AC_left_hand = gsl_matrix_complex_calloc(AC_rowsA, AC_columnsA);

            AC_mag_node_result = gsl_complex_rect(0, 0);
            AC_phase_node_result = gsl_complex_rect(0, 0);

            while(AC_sweep_ptr != NULL){
                switch(AC_sweep_ptr -> comp_type){
                    case 'R':
                    case 'r':
                        AC_resistance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                         AC_sweep_ptr -> value);
                        break;

                    case 'V':
                    case 'v':
                        AC_voltagesource_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name);
                        break;

                    case 'L':
                    case 'l':
                        AC_inductance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                         AC_sweep_ptr -> value, w);
                        break;
                    
                    case 'C':
                    case 'c':
                        AC_capacitance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                          AC_sweep_ptr -> value, w);
                        break;

                    default:
                        break;
                }
                
                AC_sweep_ptr = AC_sweep_ptr -> nxt_comp;
            }
            
            /*printf("\nLeft_hand:\n");
            for(int p = 0; p < AC_rowsA; p++){
                for(int t = 0; t < AC_columnsA; t++){
                    gsl_complex wer;

                    wer = gsl_matrix_complex_get(AC_left_hand, p, t);
                    printf("%.4lf+%.4lfI\t",GSL_REAL(wer),GSL_IMAG(wer));
                }
                printf("\n");
            }*/

            AC_x = gsl_vector_complex_calloc(AC_rowsA);

            //AC_permutation = gsl_permutation_alloc(AC_rowsA);

            AC_BiCG_solver();

            /*gsl_linalg_complex_LU_decomp(AC_left_hand, AC_permutation ,&s);
            gsl_linalg_complex_LU_solve(AC_left_hand, AC_permutation, AC_right_hand, AC_x);*/

            if((AC_mag_plot_node_name != -1) && (AC_mag_plot_node_name != 0)){
                AC_mag_node_result = gsl_vector_complex_get(AC_x, (AC_mag_plot_node_name - 1));
            }

            if((AC_phase_plot_node_name != -1) && (AC_phase_plot_node_name != 0)){
                AC_phase_node_result = gsl_vector_complex_get(AC_x, (AC_phase_plot_node_name - 1));
            }

            AC_mag_result = gsl_complex_abs(AC_mag_node_result);
            // phase in degrees
            //AC_phase_result = gsl_complex_arg(AC_phase_node_result) * (180 / M_PI);
            // phase in rads
            AC_phase_result = gsl_complex_arg(AC_phase_node_result);

            AC_mag_plot_x = w;
            AC_mag_plot_y = AC_mag_result;

            AC_phase_plot_x = w;
            AC_phase_plot_y = AC_phase_result;
            
            fprintf(AC_mag_fileptr, "%f %f\n", AC_mag_plot_x, AC_mag_plot_y);
            fprintf(AC_phase_fileptr, "%f %f\n", AC_phase_plot_x, AC_phase_plot_y);

            gsl_permutation_free(AC_permutation);
            gsl_matrix_complex_free(AC_left_hand);
            gsl_vector_complex_free(AC_x);
        }

        fprintf(AC_mag_gnuplotpipe, "%s\n", AC_mag_gnuplotcommands);
        fprintf(AC_mag_gnuplotpipe, "exit\n");
        fclose(AC_mag_fileptr);
        pclose(AC_mag_gnuplotpipe);

        fprintf(AC_phase_gnuplotpipe, "%s\n", AC_phase_gnuplotcommands);
        fprintf(AC_phase_gnuplotpipe, "exit\n");
        fclose(AC_phase_fileptr);
        pclose(AC_phase_gnuplotpipe);

        gsl_vector_complex_free(AC_right_hand);
    }

    else if(AC_LOG_flag == 1){

    }
}

void AC_BiCG_solver(){
    int i, j;
    
    int AC_iter;

    gsl_complex AC_alpha = GSL_COMPLEX_ZERO;
    gsl_complex AC_beta = GSL_COMPLEX_ZERO;

    gsl_complex AC_conj_alpha = GSL_COMPLEX_ZERO;
    gsl_complex AC_conj_beta = GSL_COMPLEX_ZERO;

    double AC_norm_r = 0;
    double AC_norm_b = 0;

    gsl_complex AC_rho = GSL_COMPLEX_ZERO;
    gsl_complex AC_rho1 = GSL_COMPLEX_ZERO;
    
    gsl_complex AC_omega = GSL_COMPLEX_ZERO;

    AC_r_iterative = gsl_vector_complex_calloc(AC_rowsA);
    AC_r_herm_iterative = gsl_vector_complex_calloc(AC_rowsA);

    AC_p_iterative = gsl_vector_complex_calloc(AC_rowsA);
    AC_p_herm_iterative = gsl_vector_complex_calloc(AC_rowsA);
    AC_p_herm_beta = gsl_vector_complex_calloc(AC_rowsA);
    AC_p_scale_iterative = gsl_vector_complex_calloc(AC_rowsA);

    AC_q_iterative = gsl_vector_complex_calloc(AC_rowsA);
    AC_q_herm_iterative = gsl_vector_complex_calloc(AC_rowsA);
    AC_q_herm_alpha = gsl_vector_complex_calloc(AC_rowsA);

    AC_z_iterative = gsl_vector_complex_calloc(AC_rowsA);
    AC_z_herm_iterative = gsl_vector_complex_calloc(AC_rowsA);

    AC_Ax_iterative = gsl_vector_complex_calloc(AC_rowsA);
    AC_bminusAx_iterative = gsl_vector_complex_calloc(AC_rowsA);

    AC_M_iterative = gsl_matrix_complex_calloc(AC_rowsA, AC_columnsA);
    AC_M_inverse_iterative = gsl_matrix_complex_calloc(AC_rowsA, AC_columnsA);

    AC_permutation = gsl_permutation_calloc(AC_rowsA);

    gsl_complex AC_temp_complex_vector;

    /*// 0 + 0I
    gsl_complex_float complex_float_zero;
    complex_float_zero.dat[0] = 0;
    complex_float_zero.dat[1] = 0;

    // 1 + 0I
    gsl_complex_float complex_float_one;
    complex_float_one.dat[0] = 1;
    complex_float_one.dat[1] = 0;*/

    for(i = 0; i < AC_rowsA; i++){
        for(j = 0; j < AC_columnsA; j++){
            if(i == j){
                AC_temp_complex_vector = gsl_matrix_complex_get(AC_left_hand, i, j);

                if((AC_temp_complex_vector.dat[0] != 0) && (AC_temp_complex_vector.dat[1] != 0)){
                    gsl_matrix_complex_set(AC_M_iterative, i, j, AC_temp_complex_vector);
                }
                else{
                    gsl_matrix_complex_set(AC_M_iterative, i, j, GSL_COMPLEX_ONE);
                }
            }
        }
    }

    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, AC_left_hand, AC_x, GSL_COMPLEX_ZERO, AC_Ax_iterative);

    gsl_vector_complex_memcpy(AC_bminusAx_iterative, AC_right_hand);

    gsl_vector_complex_sub(AC_bminusAx_iterative, AC_Ax_iterative);

    gsl_vector_complex_memcpy(AC_r_iterative, AC_bminusAx_iterative);

    gsl_vector_complex_memcpy(AC_r_herm_iterative, AC_bminusAx_iterative);

    AC_norm_b = gsl_blas_dznrm2(AC_right_hand);
    AC_norm_r = gsl_blas_dznrm2(AC_r_iterative);

    if(AC_norm_b == 0){
        AC_norm_b = 1;
    }

    AC_iter = 0;

    AC_inverse_matrix();

    while(((AC_norm_r / AC_norm_b) > ITOL) && (AC_iter < (curr_pos_oldname_table + G2_components))){
        AC_iter = AC_iter + 1;

        gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, AC_M_inverse_iterative, AC_r_iterative, GSL_COMPLEX_ZERO, AC_z_iterative);

        gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, AC_M_inverse_iterative, AC_r_herm_iterative, GSL_COMPLEX_ZERO, AC_z_herm_iterative);

        gsl_blas_zdotc(AC_r_herm_iterative, AC_z_iterative, &AC_rho);

        if(gsl_complex_abs(AC_rho) < 1e-14){
            exit(15);
        }

        if(AC_iter == 1){
            gsl_vector_complex_memcpy(AC_p_iterative, AC_z_iterative);
            gsl_vector_complex_memcpy(AC_p_herm_iterative, AC_z_herm_iterative);
        }

        else{
            AC_beta = gsl_complex_div(AC_rho, AC_rho1);
            AC_conj_beta = gsl_complex_conjugate(AC_beta);

            gsl_vector_complex_memcpy(AC_p_herm_beta, AC_p_iterative);
            gsl_vector_complex_scale(AC_p_herm_beta, AC_beta);
            gsl_vector_complex_add(AC_p_herm_beta, AC_z_iterative);
            gsl_vector_complex_memcpy(AC_p_iterative, AC_p_herm_beta);

            gsl_vector_complex_memcpy(AC_p_scale_iterative, AC_p_herm_iterative);
            gsl_vector_complex_scale(AC_p_scale_iterative, AC_conj_beta);
            gsl_vector_complex_add(AC_p_scale_iterative, AC_z_herm_iterative);
            gsl_vector_complex_memcpy(AC_p_herm_iterative, AC_p_scale_iterative);
        }

        AC_rho1.dat[0] = AC_rho.dat[0];
        AC_rho1.dat[1] = AC_rho.dat[1];

        gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, AC_left_hand, AC_p_iterative, GSL_COMPLEX_ZERO, AC_q_iterative);

        gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, AC_left_hand, AC_p_herm_iterative, GSL_COMPLEX_ZERO, AC_q_herm_iterative);

        gsl_blas_zdotc(AC_p_herm_iterative, AC_q_iterative, &AC_omega);

        if(gsl_complex_abs(AC_omega) < 1e-14){
            exit(15);
        }

        AC_alpha = gsl_complex_div(AC_rho, AC_omega);
        AC_conj_alpha = gsl_complex_conjugate(AC_alpha);

        gsl_blas_zaxpy(AC_alpha, AC_p_iterative, AC_x);

        gsl_blas_zaxpy((gsl_complex_sub(GSL_COMPLEX_ZERO, AC_alpha)), AC_q_iterative, AC_r_iterative);

        gsl_blas_zaxpy((gsl_complex_sub(GSL_COMPLEX_ZERO, AC_conj_alpha)), AC_q_herm_iterative, AC_r_herm_iterative);

        AC_norm_r = gsl_blas_dznrm2(AC_r_iterative);
    }
}

void AC_inverse_matrix(){
    int signum;

    gsl_linalg_complex_LU_decomp(AC_M_iterative, AC_permutation, &signum);

    gsl_linalg_complex_LU_invert(AC_M_iterative, AC_permutation, AC_M_inverse_iterative);
}

void AC_LU_sweep(){
    double w;

    int s;

    k_RH = 0;
    //k_LH = 0;

    int AC_step;

    ptr_comp AC_sweep_ptr;

    AC_sweep_ptr = list_head;

    char *AC_mag_plot_ptr = NULL;
    char *AC_phase_plot_ptr = NULL;

    int AC_mag_plot_node_name = 0;
    int AC_phase_plot_node_name = 0;

    double AC_mag_plot_x = 0;
    double AC_mag_plot_y = 0;

    double AC_phase_plot_x = 0;
    double AC_phase_plot_y = 0;

    gsl_complex AC_mag_node_result;
    gsl_complex AC_phase_node_result;

    double AC_mag_result = 0;
    double AC_phase_result = 0;

    FILE *AC_mag_fileptr = NULL;
    FILE *AC_mag_gnuplotpipe = NULL;
    char *AC_mag_gnuplotcommands = {"set logscale x \nplot 'AC_mag.tmp' with linespoints"};

    FILE *AC_phase_fileptr = NULL;
    FILE *AC_phase_gnuplotpipe = NULL;
    char *AC_phase_gnuplotcommands = {"set logscale x\nplot 'AC_phase.tmp' with linespoints"};

    AC_mag_fileptr = fopen("AC_mag.tmp", "w");
    AC_mag_gnuplotpipe = popen("gnuplot -persistent", "w");

    AC_phase_fileptr = fopen("AC_phase.tmp", "w");
    AC_phase_gnuplotpipe = popen("gnuplot -persistent", "w");

    if(plot_node != NULL){
        AC_mag_plot_ptr = strdup(plot_node);
        AC_phase_plot_ptr = strdup(plot_node);
    }
    
    if(print_node != NULL){
        AC_mag_plot_ptr = strdup(print_node);
        AC_phase_plot_ptr = strdup(print_node);
    }

    if((AC_mag_plot_ptr != NULL) && (AC_phase_plot_ptr != NULL)){
        AC_mag_plot_node_name = return_newpos_node_name(AC_mag_plot_ptr);

        AC_phase_plot_node_name = return_newpos_node_name(AC_phase_plot_ptr);
        
        free(AC_mag_plot_ptr);
        free(AC_phase_plot_ptr);
    }

    // fill up AC_right_hand
    while(AC_sweep_ptr != NULL){
        switch(AC_sweep_ptr -> comp_type){
        case 'V':
        case 'v':
            AC_voltagesource_RH(AC_sweep_ptr -> AC_complex_number, AC_sweep_ptr);
            break;

        case 'I':
        case 'i':
            AC_currentsource_RH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                AC_sweep_ptr -> AC_complex_number, AC_sweep_ptr);
            break;

        case 'L':
        case 'l':
            k_RH = k_RH + 1;
            break;
        
        default:
            break;
        }

        AC_sweep_ptr = AC_sweep_ptr -> nxt_comp;
    }

    /*printf("\nRight_hand:\n");
    for(int y = 0; y < AC_rowsA; y++){
        gsl_complex weee = gsl_vector_complex_get(AC_right_hand, y);

        printf("%lf+%lfI\n",GSL_REAL(weee), GSL_IMAG(weee));
    }*/

    if(AC_LIN_flag == 1){
        AC_step = ((AC_end_freq - AC_start_freq) / AC_points);

        for(w = AC_start_freq; w <= AC_end_freq; w = (w + AC_step)){
            // reset AC_sweep_ptr
            AC_sweep_ptr = list_head;

            k_LH = 0;

            AC_left_hand = gsl_matrix_complex_calloc(AC_rowsA, AC_columnsA);

            AC_mag_node_result = gsl_complex_rect(0, 0);
            AC_phase_node_result = gsl_complex_rect(0, 0);

            while(AC_sweep_ptr != NULL){
                switch(AC_sweep_ptr -> comp_type){
                    case 'R':
                    case 'r':
                        AC_resistance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                         AC_sweep_ptr -> value);
                        break;

                    case 'V':
                    case 'v':
                        AC_voltagesource_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name);
                        break;

                    case 'L':
                    case 'l':
                        AC_inductance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                         AC_sweep_ptr -> value, w);
                        break;
                    
                    case 'C':
                    case 'c':
                        AC_capacitance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                          AC_sweep_ptr -> value, w);
                        break;

                    default:
                        break;
                }
                
                AC_sweep_ptr = AC_sweep_ptr -> nxt_comp;
            }
            
            /*printf("\nLeft_hand:\n");
            for(int p = 0; p < AC_rowsA; p++){
                for(int t = 0; t < AC_columnsA; t++){
                    gsl_complex wer;

                    wer = gsl_matrix_complex_get(AC_left_hand, p, t);
                    printf("%.4lf+%.4lfI\t",GSL_REAL(wer),GSL_IMAG(wer));
                }
                printf("\n");
            }*/

            AC_x = gsl_vector_complex_calloc(AC_rowsA);

            AC_permutation = gsl_permutation_alloc(AC_rowsA);

            gsl_linalg_complex_LU_decomp(AC_left_hand, AC_permutation ,&s);

            gsl_linalg_complex_LU_solve(AC_left_hand, AC_permutation, AC_right_hand, AC_x);

            if((AC_mag_plot_node_name != -1) && (AC_mag_plot_node_name != 0)){
                AC_mag_node_result = gsl_vector_complex_get(AC_x, (AC_mag_plot_node_name - 1));
            }

            if((AC_phase_plot_node_name != -1) && (AC_phase_plot_node_name != 0)){
                AC_phase_node_result = gsl_vector_complex_get(AC_x, (AC_phase_plot_node_name - 1));
            }

            AC_mag_result = gsl_complex_abs(AC_mag_node_result);
            // phase in degrees
            //AC_phase_result = gsl_complex_arg(AC_phase_node_result) * (180 / M_PI);
            // phase in rads
            AC_phase_result = gsl_complex_arg(AC_phase_node_result);

            AC_mag_plot_x = w;
            AC_mag_plot_y = AC_mag_result;

            AC_phase_plot_x = w;
            AC_phase_plot_y = AC_phase_result;
            
            fprintf(AC_mag_fileptr, "%f %f\n", AC_mag_plot_x, AC_mag_plot_y);
            fprintf(AC_phase_fileptr, "%f %f\n", AC_phase_plot_x, AC_phase_plot_y);

            gsl_permutation_free(AC_permutation);
            gsl_matrix_complex_free(AC_left_hand);
            gsl_vector_complex_free(AC_x);
        }

        fprintf(AC_mag_gnuplotpipe, "%s\n", AC_mag_gnuplotcommands);
        fprintf(AC_mag_gnuplotpipe, "exit\n");
        fclose(AC_mag_fileptr);
        pclose(AC_mag_gnuplotpipe);

        fprintf(AC_phase_gnuplotpipe, "%s\n", AC_phase_gnuplotcommands);
        fprintf(AC_phase_gnuplotpipe, "exit\n");
        fclose(AC_phase_fileptr);
        pclose(AC_phase_gnuplotpipe);

        gsl_vector_complex_free(AC_right_hand);
    }

    else if(AC_LOG_flag == 1){
        double num_decades = log10(AC_end_freq / AC_start_freq);
        int total_points = (int)(num_decades * AC_points) + 1; // Total number of points including start and end frequencies

        double *freq = (double *)malloc(total_points * sizeof(double));
        if (freq == NULL)
        {
            printf("Memory allocation failed.\n");
            exit(12);
        }

        int i, j;
        double ratio = pow(10.0, 1.0 / AC_points);
        double curr_freq = AC_start_freq;

        for (i = 0; i < num_decades; i++)
        {
            for (j = 0; j < AC_points; j++)
            {
                freq[i * AC_points + j] = curr_freq;
                //printf("%lf\n", freq[i * AC_points + j]);
                curr_freq *= ratio;
            }
        }

        // Include the end frequency
        freq[total_points - 1] = AC_end_freq;
        //printf("%lf\n", AC_end_freq);

        for(i = 0; i < total_points; i++){
            w = freq[i];
            printf("%lf\n",w);

            AC_sweep_ptr = list_head;

            k_LH = 0;

            AC_left_hand = gsl_matrix_complex_calloc(AC_rowsA, AC_columnsA);

            AC_mag_node_result = gsl_complex_rect(0, 0);
            AC_phase_node_result = gsl_complex_rect(0, 0);

            while(AC_sweep_ptr != NULL){
                switch(AC_sweep_ptr -> comp_type){
                    case 'R':
                    case 'r':
                        AC_resistance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                         AC_sweep_ptr -> value);
                        break;

                    case 'V':
                    case 'v':
                        AC_voltagesource_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name);
                        break;

                    case 'L':
                    case 'l':
                        AC_inductance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                         AC_sweep_ptr -> value, w);
                        break;
                    
                    case 'C':
                    case 'c':
                        AC_capacitance_LH(AC_sweep_ptr -> pos_node.new_pos_node_name, AC_sweep_ptr -> neg_node.new_neg_node_name, 
                                          AC_sweep_ptr -> value, w);
                        break;

                    default:
                        break;
                }
                
                AC_sweep_ptr = AC_sweep_ptr -> nxt_comp;
            }
            
            /*printf("\nLeft_hand:\n");
            for(int p = 0; p < AC_rowsA; p++){
                for(int t = 0; t < AC_columnsA; t++){
                    gsl_complex wer;

                    wer = gsl_matrix_complex_get(AC_left_hand, p, t);
                    printf("%.4lf+%.4lfI\t",GSL_REAL(wer),GSL_IMAG(wer));
                }
                printf("\n");
            }*/

            AC_x = gsl_vector_complex_calloc(AC_rowsA);

            AC_permutation = gsl_permutation_alloc(AC_rowsA);

            gsl_linalg_complex_LU_decomp(AC_left_hand, AC_permutation ,&s);

            gsl_linalg_complex_LU_solve(AC_left_hand, AC_permutation, AC_right_hand, AC_x);

            if((AC_mag_plot_node_name != -1) && (AC_mag_plot_node_name != 0)){
                AC_mag_node_result = gsl_vector_complex_get(AC_x, (AC_mag_plot_node_name - 1));
            }

            if((AC_phase_plot_node_name != -1) && (AC_phase_plot_node_name != 0)){
                AC_phase_node_result = gsl_vector_complex_get(AC_x, (AC_phase_plot_node_name - 1));
            }

            AC_mag_result = gsl_complex_abs(AC_mag_node_result);
            // phase in degrees
            //AC_phase_result = gsl_complex_arg(AC_phase_node_result) * (180 / M_PI);
            // phase in rads
            AC_phase_result = gsl_complex_arg(AC_phase_node_result);

            AC_mag_plot_x = w;
            AC_mag_plot_y = 20 * log10(AC_mag_result);

            AC_phase_plot_x = w;
            AC_phase_plot_y = AC_phase_result;
            
            fprintf(AC_mag_fileptr, "%f %f\n", AC_mag_plot_x, AC_mag_plot_y);
            fprintf(AC_phase_fileptr, "%f %f\n", AC_phase_plot_x, AC_phase_plot_y);

            gsl_permutation_free(AC_permutation);
            gsl_matrix_complex_free(AC_left_hand);
            gsl_vector_complex_free(AC_x);
        }

        fprintf(AC_mag_gnuplotpipe, "%s\n", AC_mag_gnuplotcommands);
        fprintf(AC_mag_gnuplotpipe, "exit\n");
        fclose(AC_mag_fileptr);
        pclose(AC_mag_gnuplotpipe);

        fprintf(AC_phase_gnuplotpipe, "%s\n", AC_phase_gnuplotcommands);
        fprintf(AC_phase_gnuplotpipe, "exit\n");
        fclose(AC_phase_fileptr);
        pclose(AC_phase_gnuplotpipe);

        gsl_vector_complex_free(AC_right_hand);    
    }
}

void AC_voltagesource_RH(gsl_complex phasor, ptr_comp AC_ptr){
    if((AC_ptr -> AC_component) == 0){
        k_RH++;

        return;
    }

    int k_expression;

    gsl_complex temp_complex_dest;
    gsl_complex temp_complex_src;

    k_expression = curr_pos_oldname_table + k_RH;

    //array_b[k_expression] = array_b[k_expression] + MNA_value;
    temp_complex_src = gsl_vector_complex_get(AC_right_hand, k_expression);

    temp_complex_dest = gsl_complex_add(temp_complex_src, phasor);

    gsl_vector_complex_set(AC_right_hand, k_expression, temp_complex_dest);

    k_RH++;

    return;
}

void AC_currentsource_RH(int posnode, int negnode, gsl_complex phasor, ptr_comp AC_ptr){
    if((AC_ptr -> AC_component) == 0){
        return;
    }

    int pos_position, neg_position;

    pos_position = posnode - 1;
    neg_position = negnode - 1;

    gsl_complex temp_complex_dest;
    gsl_complex temp_complex_src;

    if(posnode == 0){
        //array_b[neg_position] = array_b[neg_position] + MNA_value;
        temp_complex_src = gsl_vector_complex_get(AC_right_hand, neg_position);

        temp_complex_dest = gsl_complex_add(temp_complex_src, phasor);

        gsl_vector_complex_set(AC_right_hand, neg_position, temp_complex_dest);
        return;
    }

    if(negnode == 0){
        //array_b[pos_position] = array_b[pos_position] - MNA_value;
        temp_complex_src = gsl_vector_complex_get(AC_right_hand, pos_position);

        temp_complex_dest = gsl_complex_sub(temp_complex_src, phasor);

        gsl_vector_complex_set(AC_right_hand, pos_position, temp_complex_dest);

        return;
    }

    //array_b[pos_position] = array_b[pos_position] - MNA_value;
    //array_b[neg_position] = array_b[neg_position] + MNA_value;
    temp_complex_src = gsl_vector_complex_get(AC_right_hand, neg_position);
    temp_complex_dest = gsl_complex_add(temp_complex_src, phasor);
    gsl_vector_complex_set(AC_right_hand, neg_position, temp_complex_dest);

    temp_complex_src = gsl_vector_complex_get(AC_right_hand, pos_position);
    temp_complex_dest = gsl_complex_sub(temp_complex_src, phasor);
    gsl_vector_complex_set(AC_right_hand, pos_position, temp_complex_dest);
}

void AC_resistance_LH(int posnode, int negnode, double value){
    int row_position, column_position;

    double AC_g;

    gsl_complex temp_complex_dest;
    gsl_complex temp_complex_src;

    AC_g = (1 / value);

    row_position = posnode - 1;
    column_position = negnode - 1;

    if(posnode == 0){
        //array_A[column_position][column_position] = array_A[column_position][column_position] + AC_g;
        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, column_position, column_position);

        temp_complex_dest = gsl_complex_add_real(temp_complex_src, AC_g);

        gsl_matrix_complex_set(AC_left_hand, column_position, column_position, temp_complex_dest);
        return;
    }

    if(negnode == 0){
        //array_A[row_position][row_position] = array_A[row_position][row_position] + AC_g;
        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, row_position, row_position);

        temp_complex_dest = gsl_complex_add_real(temp_complex_src, AC_g);

        gsl_matrix_complex_set(AC_left_hand, row_position, row_position, temp_complex_dest);

        return;
    }
    
    // if((MNA_pos_node && MNA_neg_node ) != 0)

    /*array_A[row_position][row_position] = array_A[row_position][row_position] + AC_g;
    array_A[row_position][column_position] = array_A[row_position][column_position] - AC_g;
    array_A[column_position][row_position] = array_A[column_position][row_position] - AC_g;
    array_A[column_position][column_position] = array_A[column_position][column_position] + AC_g;*/
    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, column_position, column_position);
    temp_complex_dest = gsl_complex_add_real(temp_complex_src, AC_g);
    gsl_matrix_complex_set(AC_left_hand, column_position, column_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, row_position, row_position);
    temp_complex_dest = gsl_complex_add_real(temp_complex_src, AC_g);
    gsl_matrix_complex_set(AC_left_hand, row_position, row_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, row_position, column_position);
    temp_complex_dest = gsl_complex_sub_real(temp_complex_src, AC_g);
    gsl_matrix_complex_set(AC_left_hand, row_position, column_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, column_position, row_position);
    temp_complex_dest = gsl_complex_sub_real(temp_complex_src, AC_g);
    gsl_matrix_complex_set(AC_left_hand, column_position, row_position, temp_complex_dest);
}

void AC_voltagesource_LH(int posnode, int negnode){
    int pos_position, neg_position;

    int k_expression;

    pos_position = posnode - 1;
    neg_position = negnode - 1;

    gsl_complex temp_complex_dest;
    gsl_complex temp_complex_src;

    k_expression = curr_pos_oldname_table + k_LH;

    if(posnode == 0){
        /*array_A[k_expression][neg_position] = array_A[k_expression][neg_position] - 1;
        array_A[neg_position][k_expression] = array_A[neg_position][k_expression] - 1;
        k_variable++;*/

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, neg_position);
        temp_complex_dest = gsl_complex_sub_real(temp_complex_src, 1);
        gsl_matrix_complex_set(AC_left_hand, k_expression, neg_position, temp_complex_dest);

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, neg_position, k_expression);
        temp_complex_dest = gsl_complex_sub_real(temp_complex_src, 1);
        gsl_matrix_complex_set(AC_left_hand, neg_position, k_expression, temp_complex_dest);

        k_LH++;

        return;
    }
    
    if(negnode == 0){
        /*array_A[pos_position][k_expression] = array_A[pos_position][k_expression] + 1;
        array_A[k_expression][pos_position] = array_A[k_expression][pos_position] + 1;
        k_variable++;*/

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, pos_position, k_expression);
        temp_complex_dest = gsl_complex_add_real(temp_complex_src, 1);
        gsl_matrix_complex_set(AC_left_hand, pos_position, k_expression, temp_complex_dest);

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, pos_position);
        temp_complex_dest = gsl_complex_add_real(temp_complex_src, 1);
        gsl_matrix_complex_set(AC_left_hand, k_expression, pos_position, temp_complex_dest);

        k_LH++;

        return;
    }

    /*array_A[pos_position][k_expression] = array_A[pos_position][k_expression] + 1;
    array_A[k_expression][pos_position] = array_A[k_expression][pos_position] + 1;
    array_A[k_expression][neg_position] = array_A[k_expression][neg_position] - 1;
    array_A[neg_position][k_expression] = array_A[neg_position][k_expression] - 1;*/

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, pos_position, k_expression);
    temp_complex_dest = gsl_complex_add_real(temp_complex_src, 1);
    gsl_matrix_complex_set(AC_left_hand, pos_position, k_expression, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, pos_position);
    temp_complex_dest = gsl_complex_add_real(temp_complex_src, 1);
    gsl_matrix_complex_set(AC_left_hand, k_expression, pos_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, neg_position);
    temp_complex_dest = gsl_complex_sub_real(temp_complex_src, 1);
    gsl_matrix_complex_set(AC_left_hand, k_expression, neg_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, neg_position, k_expression);
    temp_complex_dest = gsl_complex_sub_real(temp_complex_src, 1);
    gsl_matrix_complex_set(AC_left_hand, neg_position, k_expression, temp_complex_dest);

    k_LH++;
}

void AC_capacitance_LH(int posnode, int negnode, double value, double f_value){
    int row_position, column_position;

    double AC_value = 0;

    gsl_complex temp_complex_dest;
    gsl_complex temp_complex_src;

    gsl_complex AC_result;

    row_position = posnode - 1;
    column_position = negnode - 1;

    AC_value = (2 * M_PI * f_value * value);

    AC_result = gsl_complex_rect(0, AC_value);

    if(posnode == 0){
        //array_C[column_position][column_position] = array_C[column_position][column_position] + AC_value;
        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, column_position, column_position);
        temp_complex_dest = gsl_complex_add(temp_complex_src, AC_result);
        gsl_matrix_complex_set(AC_left_hand, column_position, column_position, temp_complex_dest);
        
        return;
    }

    if(negnode == 0){
        //array_C[row_position][row_position] = array_C[row_position][row_position] + AC_value;
        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, row_position, row_position);
        temp_complex_dest = gsl_complex_add(temp_complex_src, AC_result);
        gsl_matrix_complex_set(AC_left_hand, row_position, row_position, temp_complex_dest);

        return;
    }

    // if((MNA_pos_node && MNA_neg_node ) != 0)

    /*array_C[row_position][row_position] = array_C[row_position][row_position] + AC_value;
    array_C[row_position][column_position] = array_C[row_position][column_position] - AC_value;
    array_C[column_position][row_position] = array_C[column_position][row_position] - AC_value;
    array_C[column_position][column_position] = array_C[column_position][column_position] + AC_value;*/

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, row_position, row_position);
    temp_complex_dest = gsl_complex_add(temp_complex_src, AC_result);
    gsl_matrix_complex_set(AC_left_hand, row_position, row_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, row_position, column_position);
    temp_complex_dest = gsl_complex_sub(temp_complex_src, AC_result);
    gsl_matrix_complex_set(AC_left_hand, row_position, column_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, column_position, row_position);
    temp_complex_dest = gsl_complex_sub(temp_complex_src, AC_result);
    gsl_matrix_complex_set(AC_left_hand, column_position, row_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, column_position, column_position);
    temp_complex_dest = gsl_complex_add(temp_complex_src, AC_result);
    gsl_matrix_complex_set(AC_left_hand, column_position, column_position, temp_complex_dest);
}

void AC_inductance_LH(int posnode, int negnode, double value, double f_value){
    int pos_position, neg_position;

    int k_expression;

    double AC_value = 0;
    gsl_complex AC_result;

    AC_value = (2 * M_PI * f_value * value);

    AC_result = gsl_complex_rect(0, AC_value);

    pos_position = posnode - 1;
    neg_position = negnode - 1;

    gsl_complex temp_complex_dest;
    gsl_complex temp_complex_src;

    k_expression = curr_pos_oldname_table + k_LH;

    if(posnode == 0){
        /*array_A[k_expression][neg_position] = array_A[k_expression][neg_position] - 1;
        array_A[neg_position][k_expression] = array_A[neg_position][k_expression] - 1;
        k_variable++;*/

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, neg_position);
        temp_complex_dest = gsl_complex_sub_real(temp_complex_src, 1);
        gsl_matrix_complex_set(AC_left_hand, k_expression, neg_position, temp_complex_dest);

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, neg_position, k_expression);
        temp_complex_dest = gsl_complex_sub_real(temp_complex_src, 1);
        gsl_matrix_complex_set(AC_left_hand, neg_position, k_expression, temp_complex_dest);

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, k_expression);
        temp_complex_dest = gsl_complex_sub(temp_complex_src, AC_result);
        gsl_matrix_complex_set(AC_left_hand, k_expression, k_expression, temp_complex_dest);

        k_LH++;

        return;
    }
    
    if(negnode == 0){
        /*array_A[pos_position][k_expression] = array_A[pos_position][k_expression] + 1;
        array_A[k_expression][pos_position] = array_A[k_expression][pos_position] + 1;
        k_variable++;*/

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, pos_position, k_expression);
        temp_complex_dest = gsl_complex_add_real(temp_complex_src, 1);
        gsl_matrix_complex_set(AC_left_hand, pos_position, k_expression, temp_complex_dest);

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, pos_position);
        temp_complex_dest = gsl_complex_add_real(temp_complex_src, 1);
        gsl_matrix_complex_set(AC_left_hand, k_expression, pos_position, temp_complex_dest);

        temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, k_expression);
        temp_complex_dest = gsl_complex_sub(temp_complex_src, AC_result);
        gsl_matrix_complex_set(AC_left_hand, k_expression, k_expression, temp_complex_dest);

        k_LH++;

        return;
    }

    /*array_A[pos_position][k_expression] = array_A[pos_position][k_expression] + 1;
    array_A[k_expression][pos_position] = array_A[k_expression][pos_position] + 1;
    array_A[k_expression][neg_position] = array_A[k_expression][neg_position] - 1;
    array_A[neg_position][k_expression] = array_A[neg_position][k_expression] - 1;*/

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, pos_position, k_expression);
    temp_complex_dest = gsl_complex_add_real(temp_complex_src, 1);
    gsl_matrix_complex_set(AC_left_hand, pos_position, k_expression, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, pos_position);
    temp_complex_dest = gsl_complex_add_real(temp_complex_src, 1);
    gsl_matrix_complex_set(AC_left_hand, k_expression, pos_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, neg_position);
    temp_complex_dest = gsl_complex_sub_real(temp_complex_src, 1);
    gsl_matrix_complex_set(AC_left_hand, k_expression, neg_position, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, neg_position, k_expression);
    temp_complex_dest = gsl_complex_sub_real(temp_complex_src, 1);
    gsl_matrix_complex_set(AC_left_hand, neg_position, k_expression, temp_complex_dest);

    temp_complex_src = gsl_matrix_complex_get(AC_left_hand, k_expression, k_expression);
    temp_complex_dest = gsl_complex_sub(temp_complex_src, AC_result);
    gsl_matrix_complex_set(AC_left_hand, k_expression, k_expression, temp_complex_dest);

    k_LH++;
}