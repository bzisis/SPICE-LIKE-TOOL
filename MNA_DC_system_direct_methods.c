#include "MNA_DC_system_direct_methods.h"
#include "MNA_DC_system_creation.h"
#include "netlist_parser.h"
#include "MNA_DC_system_iterative_methods.h"

void LU_decomp(){
    int i, j;
    int s;

    A = gsl_matrix_alloc(rowsA, columnsA);

    X = gsl_vector_alloc(rowsx);

    B = gsl_vector_alloc(rowsb);

    P = gsl_permutation_alloc(rowsb);

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            gsl_matrix_set(A, i, j, array_A[i][j]);
        }
    }

    for(i = 0; i < rowsx; i++){
            gsl_vector_set(X, i, array_x[i]);
    }

    for(i = 0; i < rowsb; i++){
            gsl_vector_set(B, i, array_b[i]);
    }

    gsl_linalg_LU_decomp(A, P, &s);

    gsl_linalg_LU_solve(A, P, B, X);

    printf("\nLU solution vector is: \n");
    gsl_vector_fprintf(stdout, X, "%g");
}

void cholesky_decomp(){
    int i, j;

    A = gsl_matrix_alloc(rowsA, columnsA);

    X = gsl_vector_alloc(rowsx);

    B = gsl_vector_alloc(rowsb);

    for(i = 0; i < rowsA; i++){
        for(j = 0; j < columnsA; j++){
            gsl_matrix_set(A, i, j, array_A[i][j]);
        }
    }

    for(i = 0; i < rowsx; i++){
        gsl_vector_set(X, i, array_x[i]);
    }

    for(i = 0; i < rowsb; i++){
        gsl_vector_set(B, i, array_b[i]);
    }

    gsl_linalg_cholesky_decomp(A);

    gsl_linalg_cholesky_solve(A, B, X);

    printf("\nCholesky solution vector is: \n");
    gsl_vector_fprintf(stdout, X, "%g");
}

void LU_DC_op(char *filename){
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

                vector_value = gsl_vector_get(X, (vector_pos - 1));

                node_name = strdup(((ht+i)->old_name[j]));

                //node = strdup("V");

                //fprintf(fp,"%s""(""%s"")"" ""%lf\n", node, node_name, vector_value);
                fprintf(fp,"%s"" ""%lf\n", node_name, vector_value);

                free(node_name);
                //free(node);
            }
        }
    }
    
    gsl_permutation_free(P);
	gsl_matrix_free(A);
	//gsl_vector_free(X);
	gsl_vector_free(B);
    fclose(fp);
}

void cholesky_DC_op(char *filename){
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

                vector_value = gsl_vector_get(X, (vector_pos - 1));

                node_name = strdup(((ht+i)->old_name[j]));

                //node = strdup("V");

                //fprintf(fp,"%s""(""%s"")"" ""%lf\n", node, node_name, vector_value);
                fprintf(fp,"%s"" ""%lf\n", node_name, vector_value);

                free(node_name);
                //free(node);
            }
        }
    }
    
	gsl_matrix_free(A);
	//gsl_vector_free(X);
	gsl_vector_free(B);
    fclose(fp);
}

void choose_decomp(char *filename){
    if((spd_flag == 1 && (iter_flag == 0))){
        cholesky_decomp();
        cholesky_DC_op(filename);
        
        return;
    }
    if(iter_flag == 0){
        LU_decomp();
        LU_DC_op(filename);
    }
}

void dc_sweep(){
    int i, j;
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
                        //printf("%lf\n",k);
                        A = gsl_matrix_alloc(rowsA, columnsA);

                        X = gsl_vector_alloc(rowsx);

                        B = gsl_vector_alloc(rowsb);

                        P = gsl_permutation_alloc(rowsb);

                        for(i = 0; i < rowsA; i++){
                            for(j = 0; j < columnsA; j++){
                                gsl_matrix_set(A, i, j, array_A[i][j]);
                            }
                        }

                        for(i = 0; i < rowsx; i++){
                            gsl_vector_set(X, i, array_x[i]);
                        }

                        for(i = 0; i < rowsb; i++){
                            gsl_vector_set(B, i, array_b[i]);
                        }
/////////////////////////////// gsl_vector_set(B, posofvsource, k);//////////////////
                        posofvsource = (dcsweep_ptr -> vsource_pos) + curr_pos_oldname_table;
                        
                        gsl_vector_set(B, posofvsource, k);

                        gsl_linalg_LU_decomp(A, P, &s);

                        gsl_linalg_LU_solve(A, P, B, X);

                        //printf("\nSolution vector is: \n");
                        //gsl_vector_fprintf(stdout, X, "%g");

                        if((posofplotnode != -1) && (posofplotnode != 0)){
                            resultofvsourse = gsl_vector_get(X, (posofplotnode - 1));
                        }

                        if((posofprintnode != -1) && (posofprintnode != 0)){
                            resultofvsourse = gsl_vector_get(X, (posofprintnode - 1));
                        }

                        fprintf(fileptr, "%f %f\n", plotx, ploty);

                        plotx = k + 0.000001;

                        ploty = resultofvsourse + 0.000001;

                        gsl_matrix_free(A);
	                    gsl_vector_free(X);
	                    gsl_vector_free(B);
                        gsl_permutation_free(P);
                    }
                }

                dcsweep_ptr = dcsweep_ptr->nxt_comp;
            }
        }
    }

    fprintf(gnuplotpipe, "%s\n", gnuplotcommands);

    free(dcsweepname);

    fprintf(gnuplotpipe, "exit\n");

    fclose(fileptr);
    pclose(gnuplotpipe);
}

int return_newpos_node_name(char *oldnodename){
    int newposnodename;

    ptr_comp plotnodeptr;

    plotnodeptr = list_head;

    newposnodename = -1;

    while(plotnodeptr != NULL){
        if((strcmp((oldnodename), (plotnodeptr->pos_node.pos_node_name)) == 0)){
            
            newposnodename = plotnodeptr->pos_node.new_pos_node_name;

            break;
        }
        else if((strcmp((oldnodename), (plotnodeptr->neg_node.neg_node_name)) == 0)){

            newposnodename = plotnodeptr->neg_node.new_neg_node_name;

            break;
        }

        plotnodeptr = plotnodeptr->nxt_comp;
    }

    return(newposnodename);
}