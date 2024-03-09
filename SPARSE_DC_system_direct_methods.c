#include "SPARSE_DC_system_direct_methods.h"
#include "SPARSE_DC_system_creation.h"
#include "netlist_parser.h"
#include "MNA_DC_system_direct_methods.h"


void choose_sparse_decomp(char *filename){
    int i;

    sparse_b_plot = (double *)calloc(sparse_rowsA, sizeof(double));

    for(i = 0; i < sparse_rowsA; i++){
        sparse_b_plot[i] = sparse_b[i];
    }
    
    if((spd_flag == 1 && (iter_flag == 0))){
        sparse_cholesky_decomp();
        sparse_cholesky_DC_op(filename);
        
        return;
    }
    if(iter_flag == 0){
        sparse_LU_decomp();
        sparse_LU_DC_op(filename);
    }

    //free(sparse_b_plot);
}


void sparse_LU_decomp(){
    int i;

    sparse_S = cs_sqr(2, sparse_C, 0);
    sparse_N = cs_lu(sparse_C, sparse_S, 1);

    sparse_y = cs_malloc(sparse_rowsA, sizeof(double));

    cs_ipvec(sparse_N->pinv, sparse_b, sparse_y, sparse_rowsA);
    cs_lsolve(sparse_N->L, sparse_y);
    cs_usolve(sparse_N->U, sparse_y);
    cs_ipvec(sparse_S->q, sparse_y, sparse_b, sparse_rowsA);
    
    printf("\nSparse LU solution vector is: \n");
    for(i = 0; i < sparse_rowsA; i++){
        printf("%lf \n", sparse_b[i]);
    }

    cs_free(sparse_y);
    cs_sfree(sparse_S);
    cs_nfree(sparse_N);
}


void sparse_LU_DC_op(char *filename){
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
            if(((ht+i)->old_name[j]) != NULL){
                if(strcmp("0", ((ht+i)->old_name[j])) == 0){
                    continue;
                }
                vector_pos = (ht + i)->new_name_value[j];

                vector_value = sparse_b[vector_pos - 1];

                node_name = strdup(((ht+i)->old_name[j]));

                //node = strdup("V");

                //fprintf(fp,"%s""(""%s"")"" ""%lf\n", node, node_name, vector_value);
                fprintf(fp,"%s"" ""%lf\n", node_name, vector_value);

                free(node_name);
                //free(node);
            }
        }
    }

    fclose(fp);
    //destroy_sparse_matrices();
}

void sparse_cholesky_decomp(){
    int i;

    sparse_S = cs_schol(1, sparse_C);
    sparse_N = cs_chol(sparse_C, sparse_S);

    sparse_y = cs_malloc(sparse_rowsA, sizeof(double));

    cs_ipvec(sparse_S->pinv, sparse_b, sparse_y, sparse_rowsA);
    cs_lsolve(sparse_N->L, sparse_y);
    cs_ltsolve(sparse_N->L, sparse_y);
    cs_pvec(sparse_S->pinv, sparse_y, sparse_b, sparse_rowsA);

    printf("\nSparse Cholesky solution vector is: \n");
    for(i = 0; i < sparse_rowsA; i++){
        printf("%lf \n", sparse_b[i]);
    }

    cs_free(sparse_y);
    cs_sfree(sparse_S);
    cs_nfree(sparse_N);
}

void sparse_cholesky_DC_op(char *filename){
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
            if(((ht+i)->old_name[j]) != NULL){
                if(strcmp("0", ((ht+i)->old_name[j])) == 0){
                    continue;
                }
                vector_pos = (ht + i)->new_name_value[j];

                vector_value = sparse_b[vector_pos - 1];

                node_name = strdup(((ht+i)->old_name[j]));

                //node = strdup("V");

                //fprintf(fp,"%s""(""%s"")"" ""%lf\n", node, node_name, vector_value);
                fprintf(fp,"%s"" ""%lf\n", node_name, vector_value);

                free(node_name);
                //free(node);
            }
        }
    }

    fclose(fp);

    //destroy_sparse_matrices();
}

void sparse_dc_sweep(){
    int i;
    int s;
    double k;

    ptr_comp dcsweep_ptr;
    char dcsweeptype;
    char *dcsweepname;
    char *dcsweeptoken;

    int posofvsource;
    double resultofvsourse;

    double plotx = 0;
    double ploty = 0;

    int posofplotnode = 0;
    int posofprintnode = 0;

    FILE *fileptr = NULL;
    FILE *gnuplotpipe = NULL;
    char *gnuplotcommands = {"plot 'example.tmp'"};

    dcsweeptoken = strdup(DCsweep.dcsweep_component);

    dcsweeptype = DCsweep.dcsweep_component[0];

    dcsweeptoken++;

    dcsweepname = strdup(dcsweeptoken);
    
    dcsweep_ptr = list_head;

    dcsweeptoken--;
    free(dcsweeptoken);
    
    fileptr = fopen("example.tmp", "w");
    gnuplotpipe = popen("gnuplot -persistent", "w");

    if(plot_node != NULL){
        posofplotnode = return_newpos_node_name(plot_node);
    }
    
    if(print_node != NULL){
        posofprintnode = return_newpos_node_name(print_node);
    }

    if((DCsweep.dcsweep_flag == 1) && (plot_flag == 1)){
        if((DCsweep.dcsweep_component[0] == 'i') || (DCsweep.dcsweep_component[0] == 'I') || (DCsweep.dcsweep_component[0] == 'v') || 
           (DCsweep.dcsweep_component[0] == 'V')){
            
            while(dcsweep_ptr != NULL){
                if((dcsweeptype == dcsweep_ptr->comp_type) && (strcmp(dcsweepname, (dcsweep_ptr->comp_name)) == 0)){
                    for(k = DCsweep.dcsweep_startval; k < (DCsweep.dcsweep_endval + 0.000001); k = (k + DCsweep.dcsweep_step)){
                        sparse_b_plot_temp = (double *)calloc(sparse_rowsA, sizeof(double));
                        for(i = 0; i < sparse_rowsA; i++){
                            sparse_b_plot_temp[i] = sparse_b_plot[i];
                        }

                        posofvsource = (dcsweep_ptr -> vsource_pos) + curr_pos_oldname_table;
                        sparse_b_plot_temp[posofvsource] = k;

                        cs_lusol(2, sparse_C, sparse_b_plot_temp, 1);

                        if((posofplotnode != -1) && (posofplotnode != 0)){
                            resultofvsourse = sparse_b_plot_temp[posofplotnode - 1];
                        }

                        if((posofprintnode != -1) && (posofprintnode != 0)){
                            resultofvsourse = sparse_b_plot_temp[posofprintnode - 1];
                        }

                        fprintf(fileptr, "%f %f\n", plotx, ploty);

                        plotx = k + 0.000001;

                        ploty = resultofvsourse + 0.000001;

                        free(sparse_b_plot_temp);
                    }
                }

                dcsweep_ptr = dcsweep_ptr->nxt_comp;
            }
        }
    }

    //free(sparse_b_plot);

    fprintf(gnuplotpipe, "%s\n", gnuplotcommands);

    free(dcsweepname);

    fprintf(gnuplotpipe, "exit\n");

    fclose(fileptr);
    pclose(gnuplotpipe);
}