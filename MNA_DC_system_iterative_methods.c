#include "MNA_DC_system_iterative_methods.h"
#include "MNA_DC_system_direct_methods.h"
#include "MNA_DC_system_creation.h"
#include "netlist_parser.h"

void choose_iterative_decomp(char *filename){
    if((spd_flag == 1) && (iter_flag == 1)){
        CG_decomp();

        CG_DC_op(filename);
    }

    if((spd_flag == 0) && (iter_flag == 1)){
        BiCG_decomp();

        BiCG_DC_op(filename);
    }
}

void CG_decomp(){
    int i, j;
    
    int iter;

    double alpha;
    double beta;

    double norm_r = 0;
    double norm_b = 0;

    double rho = 0;
    double rho1 = 0;

    double ptrans_q = 0;

    A_iterative = gsl_matrix_alloc(rowsA, columnsA);

    X_iterative = gsl_vector_alloc(rowsx);

    B_iterative = gsl_vector_alloc(rowsb);

    r_iterative = gsl_vector_alloc(rowsx);

    p_iterative = gsl_vector_alloc(rowsx);

    p_scale_beta = gsl_vector_alloc(rowsx);

    p_scale_alpha = gsl_vector_alloc(rowsx);

    q_iterative = gsl_vector_alloc(rowsx);

    q_scale_alpha = gsl_vector_alloc(rowsx);

    z_iterative = gsl_vector_alloc(rowsx);

    Ax_iterative = gsl_vector_alloc(rowsA);

    bminusAx_iterative = gsl_vector_alloc(rowsA);

    M_iterative = gsl_matrix_alloc(rowsA, columnsA);

    M_inverse_iterative = gsl_matrix_alloc(rowsA, columnsA);

    permut_iterative = gsl_permutation_alloc(rowsA);

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            gsl_matrix_set(A_iterative, i, j, array_A[i][j]);
        }
    }

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            gsl_matrix_set(M_iterative, i, j, 0);
        }
    }

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            gsl_matrix_set(M_inverse_iterative, i, j, 0);
        }
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(X_iterative, i, array_x[i]);
    }

    for(i = 0; i < rowsb; i++){
            gsl_vector_set(B_iterative, i, array_b[i]);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(r_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(p_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(p_scale_beta, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(p_scale_alpha, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(q_iterative, i, 0);
    }
    
    for(i = 0; i < rowsx; i++){
            gsl_vector_set(q_scale_alpha, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(z_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(Ax_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(bminusAx_iterative, i, 0);
    }

    //printf("\nB2 is: \n");
    //gsl_vector_fprintf(stdout, B_iterative, "%g");
///////////////////////////////////////////////////////////////////////////////////////////////////
    iter = 0;

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            if(i == j){
                if(array_A[i][j] != 0){
                    gsl_matrix_set(M_iterative, i, j, array_A[i][j]);
                }
                else{
                    gsl_matrix_set(M_iterative, i, j, 1);
                }
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1.0, A_iterative, X_iterative, 0.0, Ax_iterative);

    //printf("\nAx is: \n");
    //gsl_vector_fprintf(stdout, Ax_iterative, "%g");

    gsl_vector_memcpy(bminusAx_iterative, B_iterative);

    //printf("\nB-Ax is: \n");
    //gsl_vector_fprintf(stdout, bminusAx_iterative, "%g");

    gsl_vector_sub(bminusAx_iterative, Ax_iterative);

    gsl_vector_memcpy(r_iterative, bminusAx_iterative);

    //printf("\nr is: \n");
    //gsl_vector_fprintf(stdout, r_iterative, "%g");

    norm_r = gsl_blas_dnrm2(r_iterative);

    norm_b = gsl_blas_dnrm2(B_iterative);

    if(norm_b == 0){    // if B_iterative = 0
        norm_b = 1;
    }

    inverse_matrix();

    while(((norm_r / norm_b) > ITOL) && (iter < (curr_pos_oldname_table + G2_components))){
        iter = iter + 1;

        preconditioner_solve();

        gsl_blas_ddot(r_iterative, z_iterative, &rho);

        if(iter == 1){
            gsl_vector_memcpy(p_iterative, z_iterative);
        }
        else{
            beta = (rho / rho1);

            gsl_vector_memcpy(p_scale_beta, p_iterative);

            gsl_vector_scale(p_scale_beta, beta);

            gsl_vector_add(p_scale_beta, z_iterative);

            gsl_vector_memcpy(p_iterative, p_scale_beta);
        }

        rho1 = rho;

        gsl_blas_dgemv(CblasNoTrans, 1.0, A_iterative, p_iterative, 0.0, q_iterative);

        gsl_blas_ddot(p_iterative, q_iterative, &ptrans_q);

        alpha = (rho / ptrans_q);

        gsl_blas_daxpy(alpha, p_iterative, X_iterative);

        gsl_blas_daxpy(-alpha, q_iterative, r_iterative);

        norm_r = gsl_blas_dnrm2(r_iterative);
    }

    printf("\nCg solution vector is: \n");
    gsl_vector_fprintf(stdout, X_iterative, "%g");
}

void preconditioner_solve(){
    gsl_blas_dgemv(CblasNoTrans, 1.0, M_inverse_iterative, r_iterative, 0.0, z_iterative);
}

void inverse_matrix(){
    int signum;

    gsl_linalg_LU_decomp(M_iterative, permut_iterative, &signum);

    gsl_linalg_LU_invert(M_iterative, permut_iterative, M_inverse_iterative);
}

void BiCG_decomp(){
    int i, j;
    
    int iter;

    double alpha;
    double beta;

    double norm_r = 0;
    double norm_b = 0;

    double rho = 0;
    double rho1 = 0;
    
    double omega = 0;

    A_iterative = gsl_matrix_alloc(rowsA, columnsA);

    X_iterative = gsl_vector_alloc(rowsx);

    B_iterative = gsl_vector_alloc(rowsb);

    r_iterative = gsl_vector_alloc(rowsx);

    p_iterative = gsl_vector_alloc(rowsx);

    p_scale_beta = gsl_vector_alloc(rowsx);

    p_scale_alpha = gsl_vector_alloc(rowsx);

    q_iterative = gsl_vector_alloc(rowsx);

    q_scale_alpha = gsl_vector_alloc(rowsx);

    z_iterative = gsl_vector_alloc(rowsx);

    Ax_iterative = gsl_vector_alloc(rowsA);

    bminusAx_iterative = gsl_vector_alloc(rowsA);

    M_iterative = gsl_matrix_alloc(rowsA, columnsA);

    M_inverse_iterative = gsl_matrix_alloc(rowsA, columnsA);

    permut_iterative = gsl_permutation_alloc(rowsA);

    r_tilde_iterative = gsl_vector_alloc(rowsx);

    z_tilde_iterative = gsl_vector_alloc(rowsx);

    p_tilde_iterative = gsl_vector_alloc(rowsx);

    q_tilde_iterative = gsl_vector_alloc(rowsx);

    p_scale_beta_tilde = gsl_vector_alloc(rowsx);

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            gsl_matrix_set(A_iterative, i, j, array_A[i][j]);
        }
    }

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            gsl_matrix_set(M_iterative, i, j, 0);
        }
    }

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            gsl_matrix_set(M_inverse_iterative, i, j, 0);
        }
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(X_iterative, i, array_x[i]);
    }

    for(i = 0; i < rowsb; i++){
            gsl_vector_set(B_iterative, i, array_b[i]);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(r_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(p_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(p_scale_beta, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(p_scale_beta_tilde, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(p_scale_alpha, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(q_iterative, i, 0);
    }
    
    for(i = 0; i < rowsx; i++){
            gsl_vector_set(q_scale_alpha, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(z_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(Ax_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(bminusAx_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(r_tilde_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(z_tilde_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(p_tilde_iterative, i, 0);
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(q_tilde_iterative, i, 0);
    }

//////////////////////////////////////////////////////////////////
    iter = 0;

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            if(i == j){
                if(array_A[i][j] != 0){
                    gsl_matrix_set(M_iterative, i, j, array_A[i][j]);
                }
                else{
                    gsl_matrix_set(M_iterative, i, j, 1);
                }
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1.0, A_iterative, X_iterative, 0.0, Ax_iterative);

    gsl_vector_memcpy(bminusAx_iterative, B_iterative);

    gsl_vector_sub(bminusAx_iterative, Ax_iterative);

    gsl_vector_memcpy(r_iterative, bminusAx_iterative);

    gsl_vector_memcpy(r_tilde_iterative, bminusAx_iterative);

    norm_r = gsl_blas_dnrm2(r_iterative);

    norm_b = gsl_blas_dnrm2(B_iterative);

    if(norm_b == 0){    // if B_iterative = 0
        norm_b = 1;
    }

    inverse_matrix();

    while(((norm_r / norm_b) > ITOL) && (iter < (curr_pos_oldname_table + G2_components))){
        iter = iter + 1;

        preconditioner_solve();
        transpose_preconditioner_solve();

        gsl_blas_ddot(r_tilde_iterative, z_iterative, &rho);

        if((fabs(rho)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        if(iter == 1){
            gsl_vector_memcpy(p_iterative, z_iterative);

            gsl_vector_memcpy(p_tilde_iterative, z_tilde_iterative);
        }
        else{
            beta = (rho / rho1);

            gsl_vector_memcpy(p_scale_beta, p_iterative);

            gsl_vector_scale(p_scale_beta, beta);

            gsl_vector_add(p_scale_beta, z_iterative);

            gsl_vector_memcpy(p_iterative, p_scale_beta);

            gsl_vector_memcpy(p_scale_beta_tilde, p_tilde_iterative);

            gsl_vector_scale(p_scale_beta_tilde, beta);

            gsl_vector_add(p_scale_beta_tilde, z_tilde_iterative);

            gsl_vector_memcpy(p_tilde_iterative, p_scale_beta_tilde);
        }

        rho1 = rho;

        gsl_blas_dgemv(CblasNoTrans, 1.0, A_iterative, p_iterative, 0.0, q_iterative);

        gsl_blas_dgemv(CblasTrans, 1.0, A_iterative, p_tilde_iterative, 0.0, q_tilde_iterative);

        gsl_blas_ddot(p_tilde_iterative, q_iterative, &omega);

        if((fabs(omega)) < 1e-14){
            exit(10);
        }

        alpha = (rho / omega);

        gsl_blas_daxpy(alpha, p_iterative, X_iterative);

        gsl_blas_daxpy(-alpha, q_iterative, r_iterative);

        gsl_blas_daxpy(-alpha, q_tilde_iterative, r_tilde_iterative);

        norm_r = gsl_blas_dnrm2(r_iterative);
    }

    printf("\nBiCG solution vector is: \n");
    gsl_vector_fprintf(stdout, X_iterative, "%g");
}

void transpose_preconditioner_solve(){
    gsl_blas_dgemv(CblasTrans, 1.0, M_inverse_iterative, r_tilde_iterative, 0.0, z_tilde_iterative);
}


void CG_DC_op(char *filename){
    int i, j;

    FILE *fp;

    double vector_value;
    int vector_pos;

    char *node_name;
    char *node;

    fp = fopen(filename, "w");

    if(!fp){
        printf("\nOperating point file could not be opened.\n");
        exit(5);
    }

    for(i = 0; i < ht_size; i++){
        for(j = 0; j < ht_depth; j++){
            if(((ht + i) -> old_name[j]) != NULL){
                if(strcmp("0", ((ht + i) -> old_name[j])) == 0){
                    continue;
                }
                vector_pos = (ht + i) -> new_name_value[j];

                vector_value = gsl_vector_get(X_iterative, (vector_pos - 1));

                node_name = strdup(((ht + i) -> old_name[j]));

                //node = strdup("V");

                //fprintf(fp,"%s""(""%s"")"" ""%lf\n", node, node_name, vector_value);
                fprintf(fp,"%s"" ""%lf\n", node_name, vector_value);
                
                free(node_name);
                //free(node);
            }
        }
    }

    //free_CG();

    fclose(fp);
}

void BiCG_DC_op(char *filename){
    int i, j;

    FILE *fp;

    double vector_value;
    int vector_pos;

    char *node_name;
    char *node;

    fp = fopen(filename, "w");

    if(!fp){
        printf("\nOperating point file could not be opened.\n");
        exit(5);
    }

    for(i = 0; i < ht_size; i++){
        for(j = 0; j < ht_depth; j++){
            if(((ht+i) -> old_name[j]) != NULL){
                if(strcmp("0", ((ht + i) -> old_name[j])) == 0){
                    continue;
                }
                vector_pos = (ht + i) -> new_name_value[j];

                vector_value = gsl_vector_get(X_iterative, (vector_pos - 1));

                node_name = strdup(((ht + i) -> old_name[j]));

                //node = strdup("V");

                //fprintf(fp,"%s""(""%s"")"" ""%lf\n", node, node_name, vector_value);
                fprintf(fp,"%s"" ""%lf\n", node_name, vector_value);

                free(node_name);
                //free(node);
            }
        }
    }

    //free_BiCG();

    fclose(fp);
}

void free_CG(){
    gsl_permutation_free(permut_iterative);

    gsl_matrix_free(A_iterative);
    gsl_matrix_free(M_iterative);
    gsl_matrix_free(M_inverse_iterative);

    gsl_vector_free(X_iterative);
    gsl_vector_free(B_iterative);
    gsl_vector_free(r_iterative);
    gsl_vector_free(p_iterative);
    gsl_vector_free(p_scale_beta);
    gsl_vector_free(p_scale_alpha);
    gsl_vector_free(q_iterative);
    gsl_vector_free(q_scale_alpha);
    gsl_vector_free(z_iterative);
    gsl_vector_free(Ax_iterative);
    gsl_vector_free(bminusAx_iterative);
}

void free_BiCG(){
    gsl_permutation_free(permut_iterative);

    gsl_matrix_free(A_iterative);
    gsl_matrix_free(M_iterative);
    gsl_matrix_free(M_inverse_iterative);

    gsl_vector_free(X_iterative);
    gsl_vector_free(B_iterative);
    gsl_vector_free(r_iterative);
    gsl_vector_free(p_iterative);
    gsl_vector_free(p_scale_beta);
    gsl_vector_free(p_scale_alpha);
    gsl_vector_free(q_iterative);
    gsl_vector_free(q_scale_alpha);
    gsl_vector_free(z_iterative);
    gsl_vector_free(Ax_iterative);
    gsl_vector_free(bminusAx_iterative);

    gsl_vector_free(p_tilde_iterative);
    gsl_vector_free(z_tilde_iterative);
    gsl_vector_free(r_tilde_iterative);
    gsl_vector_free(q_tilde_iterative);
    gsl_vector_free(p_scale_beta_tilde);
}