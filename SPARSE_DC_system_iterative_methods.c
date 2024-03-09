#include "SPARSE_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_direct_methods.h"
#include "SPARSE_DC_system_creation.h"
#include "netlist_parser.h"


void choose_sparse_iterative_decomp(char *filename){
    if((spd_flag == 1) && (iter_flag == 1)){
        sparse_CG_decomp();

        sparse_CG_DC_op(filename);
    }

    if((spd_flag == 0) && (iter_flag == 1)){
        sparse_BiCG_decomp();

        sparse_BiCG_DC_op(filename);
    }
}

void sparse_BiCG_decomp(){
    int i, j, p;

    int sparse_iter;

    double sparse_rho = 0;
    double sparse_rho1 = 0;
    double sparse_beta = 0;
    double sparse_alpha = 0;

    double sparse_norm_b = 0;
    double sparse_norm_r = 0;

    double sparse_omega = 0;
    
    sparse_x_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_x_iterative == NULL){
        printf("\nMemory allocation for sparse x_iterative failed.\n");

        exit(8);
    }

    sparse_y_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_y_iterative == NULL){
        printf("\nMemory allocation for sparse y_iterative failed.\n");

        exit(8);
    }

    sparse_r_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_r_iterative == NULL){
        printf("\nMemory allocation for sparse r_iterative failed.\n");

        exit(8);
    }

    sparse_r_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_r_tilde_iterative == NULL){
        printf("\nMemory allocation for sparse r_tilde_iterative failed.\n");

        exit(8);
    }

    sparse_p_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_p_iterative == NULL){
        printf("\nMemory allocation for sparse p_iterative failed.\n");

        exit(8);
    }

    sparse_p_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_p_tilde_iterative == NULL){
        printf("\nMemory allocation for sparse p_tilde_iterative failed.\n");

        exit(8);
    }

    sparse_q_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_q_iterative == NULL){
        printf("\nMemory allocation for sparse q_iterative failed.\n");

        exit(8);
    }

    sparse_q_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_q_tilde_iterative == NULL){
        printf("\nMemory allocation for sparse q_tilde_iterative failed.\n");

        exit(8);
    }

    sparse_z_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_z_iterative == NULL){
        printf("\nMemory allocation for sparse z_iterative failed.\n");

        exit(8);
    }

    sparse_z_tilde_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_z_tilde_iterative == NULL){
        printf("\nMemory allocation for sparse z_tilde_iterative failed.\n");

        exit(8);
    }

    sparse_m_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_m_iterative == NULL){
        printf("\nMemory allocation for sparse m_iterative failed.\n");

        exit(8);
    }

    for(j = 0; j < sparse_rowsA; j++){
        for(p = sparse_C -> p[j] ; p < sparse_C -> p[j+1] ; p++){
            if(sparse_C -> i[p] == j){
                sparse_m_iterative[j] = sparse_C -> x[p];
            }
        }
    }

    sparse_M_iterative = gsl_matrix_alloc(sparse_rowsA, sparse_rowsA);
    sparse_M_inverse_iterative = gsl_matrix_alloc(sparse_rowsA, sparse_rowsA);

    for(i = 0; i < sparse_rowsA; i++){
        for(j = 0; j < sparse_rowsA; j++){
            gsl_matrix_set(sparse_M_iterative, i, j, 0);
            gsl_matrix_set(sparse_M_inverse_iterative, i, j, 0);
        }
    }

    sparse_permut_iterative = gsl_permutation_alloc(sparse_rowsA);

    // sparse_y_iterative = sparse_C * sparse_x_iterative (y = Ax)
    for (j = 0 ; j < sparse_rowsA ; j++){
        for (p = sparse_C -> p[j] ; p < sparse_C -> p[j+1] ; p++){
            sparse_y_iterative[sparse_C -> i[p]] = sparse_y_iterative[sparse_C -> i[p]] + sparse_C -> x[p] * sparse_x_iterative[j];
        }
    }

    // sparse_r_iterative = sparse_b - sparse_y_iterative (r1 = b-y)
    for(i = 0; i < sparse_rowsA; i++){
        sparse_r_iterative[i] = sparse_b[i] - sparse_y_iterative[i];
    }

    // sparse_r_tilde_iterative = sparse_r_iterative
    for(i = 0; i < sparse_rowsA; i++){
        sparse_r_tilde_iterative[i] = sparse_r_iterative[i];
    }

    sparse_norm_r = cblas_dnrm2(sparse_rowsA, sparse_r_iterative, 1);
    sparse_norm_b = cblas_dnrm2(sparse_rowsA, sparse_b, 1);

    if(sparse_norm_b == 0){
        sparse_norm_b = 1;
    }

    sparse_iter = 0;

    sparse_inverse_matrix();

    sparse_m_cmprsd_iterative = cs_spalloc(sparse_rowsA, sparse_columnsA, sparse_rowsA, 1, 1);
    sparse_m_cmprsd_iterative->nz = sparse_rowsA;

    for(i = 0; i < sparse_rowsA; i++){
        sparse_m_cmprsd_iterative -> i[i] = i;
        sparse_m_cmprsd_iterative -> p[i] = i;
        sparse_m_cmprsd_iterative -> x[i] = gsl_matrix_get(sparse_M_inverse_iterative, i, i);
    }

    sparse_m_cc_iterative = cs_compress(sparse_m_cmprsd_iterative);
    cs_dupl(sparse_m_cc_iterative);

    while(((sparse_norm_r / sparse_norm_b) > ITOL) && (sparse_iter < sparse_rowsA)){
        sparse_iter = sparse_iter + 1;

        sparse_preconditioner_solve();
        sparse_transpose_preconditioner_solve();

        sparse_rho = dot_product(sparse_r_tilde_iterative, sparse_z_iterative, sparse_rowsA);

        if((fabs(sparse_rho)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        if(sparse_iter == 1){
            for(i = 0; i < sparse_rowsA; i++){
                sparse_p_iterative[i] = sparse_z_iterative[i];

                sparse_p_tilde_iterative[i] = sparse_z_tilde_iterative[i];
            }
        }
        else{
            sparse_beta = (sparse_rho / sparse_rho1);
            for(i = 0; i < sparse_rowsA; i++){
                sparse_p_iterative[i] = (sparse_beta * sparse_p_iterative[i]) + sparse_z_iterative[i];

                sparse_p_tilde_iterative[i] = (sparse_beta * sparse_p_tilde_iterative[i]) + sparse_z_tilde_iterative[i];
            }
        }

        sparse_rho1 = sparse_rho;

        for(j = 0; j < sparse_rowsA; j++){
            for (p = sparse_C -> p[j]; p < sparse_C -> p[j + 1] ; p++){
                sparse_q_iterative[sparse_C -> i[p]] = sparse_q_iterative[sparse_C -> i[p]]
                                                        + (sparse_C -> x[p] * sparse_p_iterative[j]);
            }
        }

        for(j = 0; j < sparse_rowsA; j++){
            sparse_q_tilde_iterative[j] = 0;
            for (p = sparse_C -> p[j]; p < sparse_C -> p[j + 1] ; p++){
                sparse_q_tilde_iterative[j] = sparse_q_tilde_iterative[j]
                                                        + (sparse_C -> x[p] * sparse_p_tilde_iterative[sparse_C->i[p]]);
            }
        }

        sparse_omega = dot_product(sparse_p_tilde_iterative, sparse_q_iterative, sparse_rowsA);

        if((fabs(sparse_omega)) < 1e-14){    // EPS = 10 ^ (-14)
            exit(9);
        }

        sparse_alpha = (sparse_rho / sparse_omega);

        cblas_daxpy(sparse_rowsA, sparse_alpha, sparse_p_iterative, 1, sparse_x_iterative, 1);

        cblas_daxpy(sparse_rowsA, -sparse_alpha, sparse_q_iterative, 1, sparse_r_iterative, 1);

        cblas_daxpy(sparse_rowsA, -sparse_alpha, sparse_q_tilde_iterative, 1, sparse_r_tilde_iterative, 1);

        sparse_norm_r = cblas_dnrm2(sparse_rowsA, sparse_r_iterative, 1);

        // reset q_iterative vector to zero after each step
        for(i = 0; i < sparse_rowsA; i++){
            sparse_q_iterative[i] = 0;
        }
    }

    printf("\nSparse BiCG solution vector is: \n");
    for(i = 0; i < sparse_rowsA; i++){
        printf("%lf \n", sparse_x_iterative[i]);
    }
}

void sparse_BiCG_DC_op(char *filename){
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

                vector_value = sparse_x_iterative[vector_pos - 1];

                node_name = strdup(((ht + i) -> old_name[j]));

                //node = strdup("V");

                //fprintf(fp,"%s""(""%s"")"" ""%lf\n", node, node_name, vector_value);
                fprintf(fp,"%s"" ""%lf\n", node_name, vector_value);

                free(node_name);
                //free(node);
            }
        }
    }

    //sparse_free_BiCG();

    fclose(fp);

}


void sparse_transpose_preconditioner_solve(){
    int p, j;

    for(j = 0; j < sparse_rowsA; j++){
        sparse_z_tilde_iterative[j] = 0;
        for(p = sparse_m_cc_iterative -> p[j] ; p < sparse_m_cc_iterative -> p[j + 1] ; p++){
            sparse_z_tilde_iterative[j] = sparse_m_cc_iterative -> x[p] * sparse_r_iterative[sparse_m_cc_iterative->i[p]];
            //sparse_z_tilde_iterative[j] = sparse_m_cc_iterative -> x[p] * sparse_r_iterative[sparse_m_cc_iterative->i[p]];
        }
    }
}

void sparse_CG_decomp(){
    int i, j, p;

    int sparse_iter;

    double sparse_norm_r = 0;
    double sparse_norm_b = 0;

    double sparse_rho = 0;
    double sparse_rho1 = 0;
    double sparse_beta = 0;
    double sparse_alpha = 0;

    sparse_x_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_x_iterative == NULL){
        printf("\nMemory allocation for sparse x_iterative failed.\n");

        exit(8);
    }

    sparse_y_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_y_iterative == NULL){
        printf("\nMemory allocation for sparse y_iterative failed.\n");

        exit(8);
    }

    sparse_r_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_r_iterative == NULL){
        printf("\nMemory allocation for sparse r_iterative failed.\n");

        exit(8);
    }

    sparse_p_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_p_iterative == NULL){
        printf("\nMemory allocation for sparse p_iterative failed.\n");

        exit(8);
    }

    sparse_m_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_m_iterative == NULL){
        printf("\nMemory allocation for sparse m_iterative failed.\n");

        exit(8);
    }

    for(j = 0; j < sparse_rowsA; j++){
        for(p = sparse_C -> p[j] ; p < sparse_C -> p[j+1] ; p++){
            if(sparse_C -> i[p] == j){
                sparse_m_iterative[j] = sparse_C -> x[p];
            }
        }
    }

    sparse_z_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_z_iterative == NULL){
        printf("\nMemory allocation for sparse z_iterative failed.\n");

        exit(8);
    }

    sparse_q_iterative = (double *)calloc(sparse_rowsA, sizeof(double));
    if(sparse_q_iterative == NULL){
            printf("\nMemory allocation for sparse q_iterative failed.\n");

            exit(8);
    }


    sparse_M_iterative = gsl_matrix_alloc(sparse_rowsA, sparse_rowsA);
    sparse_M_inverse_iterative = gsl_matrix_alloc(sparse_rowsA, sparse_rowsA);

    for(i = 0; i < sparse_rowsA; i++){
        for(j = 0; j < sparse_rowsA; j++){
            gsl_matrix_set(sparse_M_iterative, i, j, 0);
            gsl_matrix_set(sparse_M_inverse_iterative, i, j, 0);
        }
    }

    sparse_permut_iterative = gsl_permutation_alloc(sparse_rowsA);

    // sparse_y_iterative = sparse_C * sparse_x_iterative (y = Ax)
    for (j = 0 ; j < sparse_rowsA ; j++){
        for (p = sparse_C -> p[j] ; p < sparse_C -> p[j+1] ; p++){
            sparse_y_iterative[sparse_C -> i[p]] = sparse_y_iterative[sparse_C -> i[p]] + sparse_C -> x[p] * sparse_x_iterative[j];
        }
    }

    // sparse_r_iterative = sparse_b - sparse_y_iterative (r1 = b-y)
    for(i = 0; i < sparse_rowsA; i++){
        sparse_r_iterative[i] = sparse_b[i] - sparse_y_iterative[i];
    }

    sparse_iter = 0;

    sparse_norm_r = cblas_dnrm2(sparse_rowsA, sparse_r_iterative, 1);
    sparse_norm_b = cblas_dnrm2(sparse_rowsA, sparse_b, 1);

    if(sparse_norm_b == 0){
        sparse_norm_b = 1;
    }

    sparse_inverse_matrix();

    sparse_m_cmprsd_iterative = cs_spalloc(sparse_rowsA, sparse_columnsA, sparse_rowsA, 1, 1);
    sparse_m_cmprsd_iterative->nz = sparse_rowsA;

    for(i = 0; i < sparse_rowsA; i++){
        sparse_m_cmprsd_iterative -> i[i] = i;
        sparse_m_cmprsd_iterative -> p[i] = i;
        sparse_m_cmprsd_iterative -> x[i] = gsl_matrix_get(sparse_M_inverse_iterative, i, i);
    }

    sparse_m_cc_iterative = cs_compress(sparse_m_cmprsd_iterative);
    cs_dupl(sparse_m_cc_iterative);

    while(((sparse_norm_r / sparse_norm_b) > ITOL) && (sparse_iter < sparse_rowsA)){
        sparse_iter = sparse_iter + 1;
        
        // sparse_z_iterative = sparse_r_iterative
        //for(i = 0; i < sparse_rowsA; i++){
            //sparse_z_iterative[i] = sparse_r_iterative[i];
        //}

        sparse_preconditioner_solve();

        sparse_rho = dot_product(sparse_r_iterative, sparse_z_iterative, sparse_rowsA);

        if(sparse_iter == 1){
            for(i = 0; i < sparse_rowsA; i++){
                sparse_p_iterative[i] = sparse_z_iterative[i];
            }
        }

        else{
            sparse_beta = (sparse_rho / sparse_rho1);

            for(i = 0; i < sparse_rowsA; i++){
                sparse_p_iterative[i] = (sparse_beta * sparse_p_iterative[i]) + sparse_z_iterative[i];
            }
        }

        sparse_rho1 = sparse_rho;

        /*sparse_q_iterative = (double *)calloc(sparse_rowsA, sizeof(double));

        if(sparse_q_iterative == NULL){
            printf("\nMemory allocation for sparse q_iterative failed.\n");

            exit(8);
        }*/

        //cs_gaxpy(sparse_C,sparse_p_iterative, sparse_q_iterative);
        for(j = 0; j < sparse_rowsA; j++){
            for (p = sparse_C -> p[j]; p < sparse_C -> p[j+1] ; p++){
                sparse_q_iterative[sparse_C -> i[p]] = sparse_q_iterative[sparse_C -> i[p]]
                                                        + (sparse_C -> x[p] * sparse_p_iterative[j]);
            }
        }

        //dot = cblas_ddot(sparse_rowsA, sparse_p_iterative, 1, sparse_q_iterative, 1);
        //sparse_alpha = (sparse_rho / dot);
        sparse_alpha = (sparse_rho / dot_product(sparse_p_iterative, sparse_q_iterative, sparse_rowsA));
        

        /*for(i = 0; i < sparse_rowsA; i++){
                sparse_x_iterative[i] = (sparse_alpha * sparse_p_iterative[i]) + sparse_x_iterative[i];
        }*/
        cblas_daxpy(sparse_rowsA, sparse_alpha, sparse_p_iterative, 1, sparse_x_iterative, 1);

        /*for(i = 0; i < sparse_rowsA; i++){
                sparse_r_iterative[i] = sparse_r_iterative[i] - (sparse_alpha * sparse_q_iterative[i]);
        }*/
        cblas_daxpy(sparse_rowsA, -sparse_alpha, sparse_q_iterative, 1, sparse_r_iterative, 1);

        sparse_norm_r = cblas_dnrm2(sparse_rowsA, sparse_r_iterative, 1);

        // reset q_iterative vector to zero after each step
        for(i = 0; i < sparse_rowsA; i++){
            sparse_q_iterative[i] = 0;
        }
    }

    printf("\nSparse CG solution vector is: \n");
    for(i = 0; i < sparse_rowsA; i++){
        printf("%lf \n", sparse_x_iterative[i]);
    }
}

void sparse_preconditioner_solve(){
    int p, j;

    //cs_gaxpy(sparse_m_cc_iterative, sparse_r_iterative, sparse_z_iterative);
    for(j = 0; j < sparse_rowsA; j++){
        for(p = sparse_m_cc_iterative -> p[j] ; p < sparse_m_cc_iterative -> p[j+1] ; p++){
            sparse_z_iterative[sparse_m_cc_iterative -> i[p]] = sparse_m_cc_iterative -> x[p] * sparse_r_iterative[j];
        }
    }
}

void sparse_inverse_matrix(){
    int i, j;
    int signum;

    for(i = 0; i < sparse_rowsA; i++){
        for(j = 0; j < sparse_rowsA; j++){
            if(i == j){
                if(sparse_C->x[i] != 0){
                    gsl_matrix_set(sparse_M_iterative, i, j, sparse_m_iterative[i]);
                }
                else{
                    gsl_matrix_set(sparse_M_iterative, i, j, 1);
                }
            }
        }
    }

    gsl_linalg_LU_decomp(sparse_M_iterative, sparse_permut_iterative, &signum);

    gsl_linalg_LU_invert(sparse_M_iterative, sparse_permut_iterative, sparse_M_inverse_iterative);
}

double dot_product(double *v, double *u, int n){
    double result = 0.0;

    for(int i = 0; i < n; i++){
        result += v[i] * u[i];
    }

    return(result);
}

void sparse_CG_DC_op(char *filename){
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

                vector_value = sparse_x_iterative[vector_pos - 1];

                node_name = strdup(((ht + i) -> old_name[j]));

                //node = strdup("V");

                //fprintf(fp,"%s""(""%s"")"" ""%lf\n", node, node_name, vector_value);
                fprintf(fp,"%s"" ""%lf\n", node_name, vector_value);

                free(node_name);
                //free(node);
            }
        }
    }

    //sparse_free_CG();

    fclose(fp);
}

void sparse_free_CG(){
    //destroy_sparse_matrices();

    gsl_permutation_free(sparse_permut_iterative);

    gsl_matrix_free(sparse_M_iterative);
    gsl_matrix_free(sparse_M_inverse_iterative);

    cs_spfree(sparse_m_cmprsd_iterative);
    cs_spfree(sparse_m_cc_iterative);
    
    free(sparse_x_iterative);
    free(sparse_y_iterative);
    free(sparse_r_iterative);
    free(sparse_q_iterative);
    free(sparse_p_iterative);
    free(sparse_m_iterative);
    free(sparse_z_iterative);
}

void sparse_free_BiCG(){
    //destroy_sparse_matrices();

    gsl_permutation_free(sparse_permut_iterative);

    gsl_matrix_free(sparse_M_iterative);
    gsl_matrix_free(sparse_M_inverse_iterative);

    cs_spfree(sparse_m_cmprsd_iterative);
    cs_spfree(sparse_m_cc_iterative);
    
    free(sparse_x_iterative);
    free(sparse_y_iterative);
    free(sparse_r_iterative);
    free(sparse_r_tilde_iterative);
    free(sparse_q_iterative);
    free(sparse_q_tilde_iterative);
    free(sparse_p_iterative);
    free(sparse_p_tilde_iterative);
    free(sparse_m_iterative);
    free(sparse_z_iterative);
    free(sparse_z_tilde_iterative);
}

void jacobiPreconditioner(const cs *A, double *diagonal){
    // Initialize the diagonal elements of A into 'diagonal'
    int n = A->n;
    for (int i = 0; i < n; ++i)
    {
        for (int p = A -> p[i]; p < A -> p[i + 1]; ++p)
        {
            if (A -> i[p] == i)
            {
                diagonal[i] = 1 / A -> x[p]; // Store the reciprocal of diagonal elements
                break;
            }
        }
    }
}