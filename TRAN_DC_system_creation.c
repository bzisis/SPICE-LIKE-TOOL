#include "TRAN_DC_system_creation.h"
#include "SPARSE_DC_system_creation.h"
#include "MNA_DC_system_creation.h"
#include "netlist_parser.h"
#include "MNA_DC_system_direct_methods.h"
#include "MNA_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_direct_methods.h"

void choose_tran_solution_method(){
    int i, j;

    if(sparse_flag == 0){
        init_tran_MNA_matrices();

        if((spd_flag == 1 ) && (iter_flag == 0)){
            // cholesky
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    gsl_matrix_set(tran_MNA_G, i, j, array_A[i][j]);
                }
            }

            for(i = 0; i < tran_MNA_rows; i++){
                gsl_vector_set(tran_MNA_X0, i,  gsl_vector_get(X, i));
            }

            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    gsl_matrix_set(tran_MNA_C, i, j, array_C[i][j]);
                }
            }
            
            /*printf("\nCholesky tran_MNA_G:\n");
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    printf("%lf ", gsl_matrix_get(tran_MNA_G, i, j));
                }
                printf("\n");
            }

            printf("\nCholesky tran_MNA_X0:\n");
            gsl_vector_fprintf(stdout, tran_MNA_X0, "%g");

            printf("\nCholesky tran_MNA_C:\n");
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    printf("%lf ", gsl_matrix_get(tran_MNA_C, i, j));
                }
                printf("\n");
            }*/

            gsl_vector_free(X);
        }

        if((spd_flag == 0) && (iter_flag == 0)){
            // LU
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    gsl_matrix_set(tran_MNA_G, i, j, array_A[i][j]);
                }
            }

            for(i = 0; i < tran_MNA_rows; i++){
                gsl_vector_set(tran_MNA_X0, i,  gsl_vector_get(X, i));
            }

            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    gsl_matrix_set(tran_MNA_C, i, j, array_C[i][j]);
                }
            }

            /*printf("\nLU tran_MNA_G:\n");
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    printf("%lf ", gsl_matrix_get(tran_MNA_G, i, j));
                }
                printf("\n");
            }

            printf("\nLU tran_MNA_X0:\n");
            gsl_vector_fprintf(stdout, tran_MNA_X0, "%g");

            printf("\nLU tran_MNA_C:\n");
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    printf("%lf ", gsl_matrix_get(tran_MNA_C, i, j));
                }
                printf("\n");
            }*/

            gsl_vector_free(X);
        }

        if((iter_flag == 1) && (spd_flag == 1)){
            // CG
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    gsl_matrix_set(tran_MNA_G, i, j, array_A[i][j]);
                }
            }

            for(i = 0; i < tran_MNA_rows; i++){
                gsl_vector_set(tran_MNA_X0, i,  gsl_vector_get(X_iterative, i));
            }

            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    gsl_matrix_set(tran_MNA_C, i, j, array_C[i][j]);
                }
            }

            /*printf("\nCG tran_MNA_G:\n");
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    printf("%lf ", gsl_matrix_get(tran_MNA_G, i, j));
                }
                printf("\n");
            }

            printf("\nCG tran_MNA_X0:\n");
            gsl_vector_fprintf(stdout, tran_MNA_X0, "%g");

            printf("\nCG tran_MNA_C:\n");
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    printf("%lf ", gsl_matrix_get(tran_MNA_C, i, j));
                }
                printf("\n");
            }*/

            free_CG();
        }

        if((iter_flag == 1) && (spd_flag == 0)){
            // BiCG
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    gsl_matrix_set(tran_MNA_G, i, j, array_A[i][j]);
                }
            }

            for(i = 0; i < tran_MNA_rows; i++){
                gsl_vector_set(tran_MNA_X0, i,  gsl_vector_get(X_iterative, i));
            }

            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    gsl_matrix_set(tran_MNA_C, i, j, array_C[i][j]);
                }
            }

            /*printf("\nBiCG tran_MNA_G:\n");
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    printf("%lf ", gsl_matrix_get(tran_MNA_G, i, j));
                }
                printf("\n");
            }

            printf("\nBiCG tran_MNA_X0:\n");
            gsl_vector_fprintf(stdout, tran_MNA_X0, "%g");

            printf("\nBiCG tran_MNA_C:\n");
            for(i = 0; i < tran_MNA_rows; i++){
                for(j = 0; j < tran_MNA_columns; j++){
                    printf("%lf ", gsl_matrix_get(tran_MNA_C, i, j));
                }
                printf("\n");
            }*/

            free_BiCG();
        }

        //free_tran_MNA_matrices();
    }

    if(sparse_flag == 1){
        init_tran_SPARSE_matrices();

        if((spd_flag == 1 ) && (iter_flag == 0)){
            // SPARSE cholesky
            // tran_sparse_tilde_C
            // tran_sparse_tilde_G

            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_SPARSE_X0[i] = sparse_b[i];
            }

            /*printf("\nSPARSE Cholesky tran_SPARSE_X0:\n");
            for(i = 0; i < tran_SPARSE_rows; i++){
                printf("%lf\n", tran_SPARSE_X0[i]);
            }*/
        }

        if((spd_flag == 0) && (iter_flag == 0)){
            // SPARSE LU
            // tran_sparse_tilde_C
            
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_SPARSE_X0[i] = sparse_b[i];
            }

            /*printf("\nSPARSE LU tran_SPARSE_X0:\n");
            for(i = 0; i < tran_SPARSE_rows; i++){
                printf("%lf\n", tran_SPARSE_X0[i]);
            }*/

        }

        if((iter_flag == 1) && (spd_flag == 1)){
            // SPARSE CG
            // tran_sparse_tilde_C
            // tran_sparse_tilde_C
            
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_SPARSE_X0[i] = sparse_x_iterative[i];
            }

            /*printf("\nSPARSE CG tran_SPARSE_X0:\n");
            for(i = 0; i < tran_SPARSE_rows; i++){
                printf("%lf\n", tran_SPARSE_X0[i]);
            }*/

            sparse_free_CG();
        }

        if((iter_flag == 1) && (spd_flag == 0)){
            // SPARSE BiCG
            // tran_sparse_tilde_C
            // tran_sparse_tilde_C
            
            for(i = 0; i < tran_SPARSE_rows; i++){
                tran_SPARSE_X0[i] = sparse_x_iterative[i];
            }

            /*printf("\nSPARSE BiCG tran_SPARSE_X0:\n");
            for(i = 0; i < tran_SPARSE_rows; i++){
                printf("%lf\n", tran_SPARSE_X0[i]);
            }*/

            sparse_free_BiCG();
        }
    }

    //free_tran_SPARSE_matrices();
}

void init_tran_MNA_matrices(){
    tran_MNA_C = NULL;
    tran_MNA_G = NULL;
    tran_MNA_X0 = NULL;
    tran_MNA_e = NULL;

    tran_MNA_rows = curr_pos_oldname_table + G2_components;
    tran_MNA_columns = curr_pos_oldname_table + G2_components;

    tran_MNA_G = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);

    tran_MNA_C = gsl_matrix_alloc(tran_MNA_rows, tran_MNA_columns);

    tran_MNA_X0 = gsl_vector_alloc(tran_MNA_rows);

    tran_MNA_e = gsl_vector_alloc(tran_MNA_rows);
}

void init_tran_SPARSE_matrices(){
    //tran_SPARSE_tilde_C = NULL;
    //tran_SPARSE_tilde_G = NULL;
    tran_SPARSE_X0 = NULL;
    tran_SPARSE_e = NULL;

    tran_SPARSE_rows = curr_pos_oldname_table + G2_components;
    tran_SPARSE_columns = curr_pos_oldname_table + G2_components;

    tran_SPARSE_X0 = (double *)calloc(tran_SPARSE_rows, sizeof(double));
    if(tran_SPARSE_X0 == NULL){
        printf("\nMemory allocation for TRANSIENT SPARSE MATRIX failed.\n");

        exit(10);
    }

    tran_SPARSE_e = (double *)calloc(tran_SPARSE_rows, sizeof(double));
    if(tran_SPARSE_e == NULL){
        printf("\nMemory allocation for TRANSIENT SPARSE MATRIX failed.\n");

        exit(10);
    }
}

void free_tran_SPARSE_matrices(){
    free(tran_SPARSE_X0);
    free(tran_SPARSE_e);
}

void free_tran_MNA_matrices(){
    gsl_matrix_free(tran_MNA_G);
    gsl_matrix_free(tran_MNA_C);
    gsl_vector_free(tran_MNA_X0);
    gsl_vector_free(tran_MNA_e);

    destroy_global_arrays(array_C, tran_MNA_rows);
}