/**************************************************/
/*         Circuit Simulation Algorithms          */
/**************************************************/ 
/* Authors: Georgia-Despoina Tzanetou             */
/*          Zisis Balatsos                        */
/*                                                */
/*       University of Thessaly Volos, 2023       */
/*                                                */
/**************************************************/



#include "netlist_parser.h"
#include "MNA_DC_system_creation.h"
#include "MNA_DC_system_direct_methods.h"
#include "MNA_DC_system_iterative_methods.h"
#include "SPARSE_DC_system_creation.h"
#include "SPARSE_DC_system_direct_methods.h"
#include "SPARSE_DC_system_iterative_methods.h"
#include "TRAN_DC_system_creation.h"
#include "TRAN_sweep.h"
#include "AC_TRAN_sweep.h"

void init_global_vars(){
    component_counter = 0;

    list_head = NULL;

    curr_pos_oldname_table = 0;

    G2_components = 0;

    spd_flag = 0;

    DCsweep.dcsweep_flag = 0;
    DCsweep.dcsweep_component = NULL;
    DCsweep.dcsweep_startval = 0;
    DCsweep.dcsweep_endval = 0;
    DCsweep.dcsweep_step = 0;

    plot_flag = 0;

    plot_node = NULL;

    iter_flag = 0;

    ITOL = 0;

    sparse_flag = 0;

    tran_flag = 0;

    tran_fin_time = 0;

    tran_time_step = 0;

    trapezoidal_flag = 0;
    backwardeuler_flag = 0;

    tran_m_parameter = 0;

    AC_flag = 0;
    AC_LIN_flag = 0;
    AC_LOG_flag = 0;

    AC_start_freq = 0;
    AC_end_freq = 0;

    AC_points = 0;
}

void init_components(ptr_comp curr_comp){
    curr_comp -> comp_type = 0;

    curr_comp -> comp_name = NULL;

    curr_comp -> pos_node.new_pos_node_name = 0;
    curr_comp -> pos_node.pos_node_name = NULL;
    curr_comp -> pos_node.pos_node_value = 0;

    curr_comp -> neg_node.new_neg_node_name = 0;
    curr_comp -> neg_node.neg_node_name = NULL;
    curr_comp -> neg_node.neg_node_value = 0;

    curr_comp -> value = 0;

    curr_comp -> nxt_comp = NULL;

    curr_comp -> vsource_pos = 0;

    curr_comp -> exp_flag = 0;
    curr_comp -> pulse_flag = 0;
    curr_comp -> sin_flag = 0;
    curr_comp -> pwl_flag = 0;
    curr_comp -> dc_flag = 0;

    curr_comp -> exp_field.exp_i1 = 0;
    curr_comp -> exp_field.exp_i2 = 0;
    curr_comp -> exp_field.exp_tc1 = 0;
    curr_comp -> exp_field.exp_tc2 = 0;
    curr_comp -> exp_field.exp_td1 = 0;
    curr_comp -> exp_field.exp_td2 = 0;

    curr_comp -> sin_field.sin_i1 = 0;
    curr_comp -> sin_field.sin_ia = 0;
    curr_comp -> sin_field.sin_fr = 0;
    curr_comp -> sin_field.sin_td = 0;
    curr_comp -> sin_field.sin_df = 0;
    curr_comp -> sin_field.sin_ph = 0;

    curr_comp -> pulse_field.pulse_i1 = 0;
    curr_comp -> pulse_field.pulse_i2 = 0;
    curr_comp -> pulse_field.pulse_td = 0;
    curr_comp -> pulse_field.pulse_tr = 0;
    curr_comp -> pulse_field.pulse_tf = 0;
    curr_comp -> pulse_field.pulse_pw = 0;
    curr_comp -> pulse_field.pulse_per = 0;

    curr_comp -> pwl_field = NULL;

    curr_comp -> comp_pwl_counter = 0;

    curr_comp -> AC_mag = 0;
    curr_comp -> AC_phase = 0;
    curr_comp -> AC_component = 0;
    
    /*curr_comp -> AC_complex_number.imag_part = 0;
    curr_comp -> AC_complex_number.real_part = 0;
    curr_comp -> AC_complex_number.complex_num = 0;*/
    
}

void add_components(FILE *fp, const char* word){
    init_global_vars();
    
    ptr_comp new_comp = NULL;

    ptr_comp prv_comp = NULL;

    char *ptr_line = NULL;

    char *ptr_token = NULL;

    ssize_t read;
    size_t line_size = 0;
    
    int token_flag = 0;

    int head_flag = 0;

    int prv_flag = 0;

    int pwl_counter = 0;

    char *spd_ptr = NULL;
    char *alt_spd_ptr = NULL;
    char *spd_token = NULL;
    char *spd_line = NULL;

    char *dcsweep_ptr = NULL;
    char *dcsweep_token = NULL;
    char *dcsweep_line = NULL;

    char *plot_ptr = NULL;
    char *plot_token = NULL;
    char *plot_line = NULL;

    char *iter_ptr = NULL;
    char *alt_iter_ptr = NULL;
    char *alternative_iter_ptr = NULL;
    char *iter_line = NULL;
    char *iter_token = NULL;

    char *itol_ptr = NULL;
    char *itol_line = NULL;
    char *itol_token = NULL;

    char *sparse_ptr = NULL;
    char *alt_sparse_ptr = NULL;
    char *alternative_sparse_ptr = NULL;
    char *sparse_line = NULL;
    char *sparse_token = NULL;

    char *plotnode_ptr = NULL;
    char *plotnode_line = NULL;
    char *plotnode_token = NULL;

    char *printnode_ptr = NULL;
    char *printnode_line = NULL;
    char *printnode_token = NULL;

    char *print_ptr = NULL;

    char *tran_ptr = NULL;
    char *tran_token = NULL;
    char *tran_line = NULL;

    char *tran_method_ptr = NULL;
    char *tran_method_line = NULL;
    char *tran_method_token = NULL;

    char *ACsweep_ptr = NULL;
    char *ACsweep_token = NULL;
    char *ACsweep_line = NULL;

    spd_ptr = "SPD";
    alt_spd_ptr = "SPD\n";

    iter_ptr = "ITER";
    alt_iter_ptr = "ITER\n";
    alternative_iter_ptr = "ITER\r\n";

    sparse_ptr = "SPARSE";
    alt_sparse_ptr = "SPARSE\n";
    alternative_sparse_ptr = "SPARSE\r\n";

    tran_method_ptr = ".OPTIONS METHOD=";

    itol_ptr = "ITOL=";

    dcsweep_ptr = ".DC";

    plot_ptr = ".PLOT";

    print_ptr = ".PRINT";

    plotnode_ptr = ".PLOT V(";

    printnode_ptr = ".PRINT V(";

    tran_ptr = ".TRAN";

    ACsweep_ptr = ".AC";

    fseek(fp, 0, SEEK_SET);

    while(!feof(fp)){
        while((read = getline(&ptr_line, &line_size, fp)) != -1){
            dcsweep_line = strdup(ptr_line);
            dcsweep_token = strtok(dcsweep_line, word);

            plot_line = strdup(ptr_line);
            plot_token = strtok(plot_line, word);

            itol_line = strdup(ptr_line);
            itol_token = strstr(itol_line, itol_ptr);

            plotnode_line = strdup(ptr_line);
            plotnode_token = strstr(plotnode_line, plotnode_ptr);

            printnode_line = strdup(ptr_line);
            printnode_token = strstr(printnode_line, printnode_ptr);

            tran_line = strdup(ptr_line);
            tran_token = strstr(tran_line, tran_ptr);

            ACsweep_line = strdup(ptr_line);
            ACsweep_token = strtok(ACsweep_line, word);

            if(tran_token != NULL){
                if(tran_line[0] != '*'){
                    tran_flag = 1;

                    tran_token = strtok(NULL, word);
                    
                    tran_time_step = strtod(tran_token, NULL);

                    tran_token = strtok(NULL, word);

                    tran_fin_time = strtod(tran_token, NULL);
                }
            }

            if(plotnode_token != NULL){
                if(plotnode_line[0] != '*'){
                    plotnode_token = strstr(plotnode_token, "(");

                    plotnode_token++;

                    plotnode_token = strtok(plotnode_token, ")");

                    plot_node = strdup(plotnode_token);

                    plotnode_token = NULL;
                }
            }

            if(printnode_token != NULL){
                if(printnode_line[0] != '*'){
                    printnode_token = strstr(printnode_token, "(");

                    printnode_token++;

                    printnode_token = strtok(printnode_token, ")");

                    print_node = strdup(printnode_token);

                    printnode_token = NULL;
                }
            }

            if(itol_token != NULL){
                if(itol_line[0] != '*'){
                    itol_token = strstr(itol_token, "=");

                    itol_token++;
                    itol_token = strtok(itol_token, "\n");

                    ITOL = atof(itol_token);
                }
            }

            if(strcmp(ACsweep_token, ACsweep_ptr) == 0){
                AC_flag = 1;

                ACsweep_token = strtok(NULL, word);

                if(strcmp(ACsweep_token, "LIN") == 0){
                    AC_LIN_flag = 1;
                }

                else if(strcmp(ACsweep_token, "LOG") == 0){
                    AC_LOG_flag = 1;
                }

                ACsweep_token = strtok(NULL, word);

                AC_points = strtod(ACsweep_token, NULL);

                ACsweep_token = strtok(NULL, word);

                AC_start_freq = strtod(ACsweep_token, NULL);

                ACsweep_token = strtok(NULL, word);

                ACsweep_token = strtok(ACsweep_token, "\n");

                AC_end_freq = strtod(ACsweep_token, NULL);
            }

            if(strcmp(dcsweep_token, dcsweep_ptr) == 0){
                DCsweep.dcsweep_flag = 1;

                dcsweep_token = strtok(NULL, word);

                DCsweep.dcsweep_component = strdup(dcsweep_token);

                dcsweep_token = strtok(NULL, word);

                DCsweep.dcsweep_startval = strtod(dcsweep_token, NULL);

                dcsweep_token = strtok(NULL, word);

                DCsweep.dcsweep_endval = strtod(dcsweep_token, NULL);

                dcsweep_token = strtok(NULL, word);

                DCsweep.dcsweep_step = strtod(dcsweep_token, NULL);

                dcsweep_token = NULL;
            }

            if((strcmp(plot_token, plot_ptr) == 0) || (strcmp(plot_token, print_ptr) == 0)){
                plot_flag = 1;
            }

            if(strstr(ptr_line, spd_ptr) != NULL){
                spd_line = strdup(ptr_line);

                spd_token = strtok(spd_line, word);

                if(spd_token[0] != '*'){
                    while(spd_token != NULL){
                        if((strcmp(spd_token, spd_ptr) == 0) || (strcmp(spd_token, alt_spd_ptr) == 0)){
                            spd_flag = 1;

                            break;
                        }

                        spd_token = strtok(NULL, word);
                    }
                    
                    free(spd_line);
                }
                else{
                    free(spd_line);
                }
            }

            if(strstr(ptr_line, tran_method_ptr) != NULL){
                tran_method_line = strdup(ptr_line);

                tran_method_token = strtok(tran_method_line, word);

                if(tran_method_token[0] != '*'){
                    while(tran_method_token != NULL){
                        tran_method_token = strtok(NULL, word);

                        if((strcmp(tran_method_token, "METHOD=TR") == 0) || strcmp(tran_method_token, "METHOD=TR\n") == 0){
                            trapezoidal_flag = 1;

                            break;
                        }
                        
                        if((strcmp(tran_method_token, "METHOD=BE") == 0) || (strcmp(tran_method_token, "METHOD=BE\n") == 0)){
                            backwardeuler_flag = 1;

                            break;
                        }
                    }
                }
                free(tran_method_line);
            }


            if(strstr(ptr_line, iter_ptr) != NULL){
                iter_line = strdup(ptr_line);

                iter_token = strtok(iter_line, word);

                if(iter_token[0] != '*'){
                    while(iter_token != NULL){
                        if((strcmp(iter_token, iter_ptr) == 0) || (strcmp(iter_token, alt_iter_ptr) == 0) ||
                            (strcmp(iter_token, alternative_iter_ptr) == 0)){
                            iter_flag = 1;

                            break;
                        }

                        iter_token = strtok(NULL, word);
                    }
                }
                free(iter_line);
            }

            if(strstr(ptr_line, sparse_ptr) != NULL){
                sparse_line = strdup(ptr_line);

                sparse_token = strtok(sparse_line, word);

                if(sparse_token[0] != '*'){
                    while(sparse_token != NULL){
                        if((strcmp(sparse_token, sparse_ptr) == 0) || (strcmp(sparse_token, alt_sparse_ptr) == 0) || 
                            (strcmp(sparse_token, alternative_sparse_ptr) == 0)){
                            sparse_flag = 1;

                            break;
                        }

                        sparse_token = strtok(NULL, word);
                    }
                }
                free(sparse_line);
            }


            if((ptr_line[0] == '*') || (ptr_line[0] == '\r') || (ptr_line[0] == '\n') || (ptr_line[0] == ' ') || (ptr_line[0] == '.')){
                free(dcsweep_line);
                free(plot_line);
                free(itol_line);
                free(plotnode_line);
                free(printnode_line);
                free(tran_line);
                free(ACsweep_line);
                //free(tran_method_line);

                continue;
            }

            new_comp = malloc(sizeof(component));

            if(!new_comp){
                printf("\nMemory allocation for parsing failed.\n");

                exit(1);
            }

            init_components(new_comp);

            if(head_flag == 0){
                list_head = new_comp;
            }

            ptr_token = strtok(ptr_line, word);

            while(ptr_token != NULL){
                if(token_flag == 0){
                    new_comp -> comp_type = ptr_token[0];

                    if((new_comp -> comp_type == 'v') || (new_comp -> comp_type == 'V') ||
                       (new_comp -> comp_type == 'l') || (new_comp -> comp_type == 'L')){
                        
                        new_comp -> vsource_pos = G2_components;

                        G2_components++;
                    }

                    ptr_token++;

                    new_comp -> comp_name = strdup(ptr_token);

                    ptr_token = strtok(NULL, word);
                    
                    token_flag++;
                }
                else if(token_flag == 1){
                    new_comp -> pos_node.pos_node_name = strdup(ptr_token);

                    ptr_token = strtok(NULL, word);
                    
                    token_flag++;
                }
                else if(token_flag == 2){
                    new_comp -> neg_node.neg_node_name = strdup(ptr_token);

                    ptr_token = strtok(NULL, word);
                    
                    token_flag++;
                }
                else if(token_flag == 3){
                    if((new_comp->comp_type == 'R') || (new_comp->comp_type == 'r') || (new_comp->comp_type == 'L') 
                        || (new_comp->comp_type == 'l') || (new_comp->comp_type == 'c') || (new_comp->comp_type == 'C')){
                        
                        new_comp -> value = strtod(ptr_token, NULL);

                        ptr_token = NULL;
                        token_flag = 0;

                        if(head_flag == 0){
                            prv_comp = new_comp;

                            head_flag = 1;
                        }

                        if(prv_flag == 1){
                            prv_comp -> nxt_comp = new_comp;

                            prv_comp = new_comp;
                        }
                    }
                    
                    else if((new_comp->comp_type == 'V') || (new_comp->comp_type == 'v') || (new_comp->comp_type == 'I')
                            || (new_comp->comp_type == 'i')){
                        
                        new_comp -> value = strtod(ptr_token, NULL);
                    
                        ptr_token = strtok(NULL, word);

                        if(ptr_token == NULL){
                            new_comp -> dc_flag = 1;

                            token_flag = 0;

                            if(head_flag == 0){
                                prv_comp = new_comp;

                                head_flag = 1;
                            }

                            if(prv_flag == 1){
                                prv_comp -> nxt_comp = new_comp;

                                prv_comp = new_comp;
                            }

                            break;
                        }

                        if(strcmp(ptr_token, "\n") == 0){
                            new_comp -> dc_flag = 1;

                            token_flag = 0;

                            if(head_flag == 0){
                                prv_comp = new_comp;

                                head_flag = 1;
                            }

                            if(prv_flag == 1){
                                prv_comp -> nxt_comp = new_comp;

                                prv_comp = new_comp;
                            }

                            break;
                        }

                        if(strcmp(ptr_token, "EXP") == 0){
                            new_comp -> exp_flag = 1;

                            token_flag++;
                        }

                        if(strcmp(ptr_token, "PULSE") == 0){
                            new_comp -> pulse_flag = 1;

                            token_flag++;
                        }

                        if(strcmp(ptr_token, "SIN") == 0){
                            new_comp -> sin_flag = 1;

                            token_flag++;
                        }

                        if(strcmp(ptr_token, "PWL") == 0){
                            new_comp -> pwl_flag = 1;

                            token_flag++;
                        }

                        //ptr_token = NULL;
                        //token_flag = 0;

                        if(strcmp(ptr_token, "AC") == 0){
                            new_comp -> AC_component = 1;

                            ptr_token = strtok(NULL, word);

                            new_comp -> AC_mag = strtod(ptr_token, NULL);

                            ptr_token = strtok(NULL, word);

                            new_comp -> AC_phase = strtod(ptr_token, NULL);

                            //new_comp -> AC_complex_number.real_part = (new_comp -> AC_mag) * cos(((new_comp -> AC_phase) * M_PI) / 180);

                            //new_comp -> AC_complex_number.imag_part = (new_comp -> AC_mag) * sin(((new_comp -> AC_phase) * M_PI) / 180);

                            //new_comp -> AC_complex_number.complex_num = CMPLX((new_comp -> AC_complex_number.real_part), 
                                                                              //(new_comp -> AC_complex_number.imag_part));
                            new_comp -> AC_complex_number = gsl_complex_rect((new_comp -> AC_mag) * cos(((new_comp -> AC_phase) * M_PI) / 180),
                                                                             (new_comp -> AC_mag) * sin(((new_comp -> AC_phase) * M_PI) / 180));
                            //token_flag++;
                        }
                    }

                    //ptr_token = NULL;
                    //token_flag = 0;

                    /*if(head_flag == 0){
                        prv_comp = new_comp;

                        head_flag = 1;
                    }

                    if(prv_flag == 1){
                        prv_comp -> nxt_comp = new_comp;

                        prv_comp = new_comp;
                    }*/
                }
                else if(token_flag == 4){
                    if(new_comp -> exp_flag == 1){
                        ptr_token = strtok(NULL, word);

                        ptr_token++;
                        new_comp -> exp_field.exp_i1 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> exp_field.exp_i2 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> exp_field.exp_td1 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> exp_field.exp_tc1 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> exp_field.exp_td2 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        strstr(ptr_token, ")");
                        new_comp -> exp_field.exp_tc2 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);

                        if((ptr_token != NULL) && (strcmp(ptr_token, "AC") == 0)){
                            new_comp -> AC_component = 1;

                            ptr_token = strtok(NULL, word);
                            new_comp -> AC_mag = strtod(ptr_token, NULL);
                            ptr_token = strtok(NULL, word);
                            new_comp -> AC_phase = strtod(ptr_token, NULL);
                            new_comp -> AC_complex_number = gsl_complex_rect((new_comp -> AC_mag) * cos(((new_comp -> AC_phase) * M_PI) / 180),
                                                                             (new_comp -> AC_mag) * sin(((new_comp -> AC_phase) * M_PI) / 180));
                        }

                        ptr_token = NULL;
                    }

                    if(new_comp -> pulse_flag == 1){
                        ptr_token = strtok(NULL, word);

                        ptr_token++;
                        new_comp -> pulse_field.pulse_i1 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> pulse_field.pulse_i2 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> pulse_field.pulse_td = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> pulse_field.pulse_tr = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> pulse_field.pulse_tf = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> pulse_field.pulse_pw = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        strstr(ptr_token, ")");
                        new_comp -> pulse_field.pulse_per = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);

                        if((ptr_token != NULL) && (strcmp(ptr_token, "AC") == 0)){
                            new_comp -> AC_component = 1;

                            ptr_token = strtok(NULL, word);
                            new_comp -> AC_mag = strtod(ptr_token, NULL);
                            ptr_token = strtok(NULL, word);
                            new_comp -> AC_phase = strtod(ptr_token, NULL);
                            new_comp -> AC_complex_number = gsl_complex_rect((new_comp -> AC_mag) * cos(((new_comp -> AC_phase) * M_PI) / 180),
                                                                             (new_comp -> AC_mag) * sin(((new_comp -> AC_phase) * M_PI) / 180));
                        }

                        ptr_token = NULL;
                    }

                    if(new_comp -> sin_flag == 1){
                        ptr_token = strtok(NULL, word);

                        ptr_token++;
                        new_comp -> sin_field.sin_i1 = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> sin_field.sin_ia = strtod(ptr_token, NULL);
                        
                        ptr_token = strtok(NULL, word);
                        new_comp -> sin_field.sin_fr = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> sin_field.sin_td = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        new_comp -> sin_field.sin_df = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        strstr(ptr_token, ")");
                        new_comp -> sin_field.sin_ph = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);

                        if((ptr_token != NULL) && (strcmp(ptr_token, "AC") == 0)){
                            new_comp -> AC_component = 1;

                            ptr_token = strtok(NULL, word);
                            new_comp -> AC_mag = strtod(ptr_token, NULL);
                            ptr_token = strtok(NULL, word);
                            new_comp -> AC_phase = strtod(ptr_token, NULL);
                            new_comp -> AC_complex_number = gsl_complex_rect((new_comp -> AC_mag) * cos(((new_comp -> AC_phase) * M_PI) / 180),
                                                                             (new_comp -> AC_mag) * sin(((new_comp -> AC_phase) * M_PI) / 180));
                        }

                        ptr_token = NULL;
                    }

                    if(new_comp -> pwl_flag == 1){
                        ptr_token = strtok(NULL, word);
                        ptr_token++;

                        new_comp -> pwl_field = (pwl_spec)malloc(sizeof(struct pwl_t));

                        new_comp -> pwl_field -> pwl_start = strtod(ptr_token, NULL);

                        ptr_token = strtok(NULL, word);
                        strstr(ptr_token, ")");

                        new_comp -> pwl_field -> pwl_end = strtod(ptr_token, NULL);

                        pwl_counter = 1;

                        while(ptr_token != NULL){
                            ptr_token = strtok(NULL, word);

                            if((ptr_token != NULL) && (strcmp(ptr_token, "AC") == 0)){
                                new_comp -> AC_component = 1;

                                ptr_token = strtok(NULL, word);

                                new_comp -> AC_mag = strtod(ptr_token, NULL);

                                ptr_token = strtok(NULL, word);

                                new_comp -> AC_phase = strtod(ptr_token, NULL);

                                new_comp -> AC_complex_number = gsl_complex_rect((new_comp -> AC_mag) * cos(((new_comp -> AC_phase) * M_PI) / 180),
                                                                             (new_comp -> AC_mag) * sin(((new_comp -> AC_phase) * M_PI) / 180));

                                ptr_token = NULL;
                            }
                            
                            if(ptr_token == NULL){
                                token_flag = 0;

                                /*if(head_flag == 0){
                                    prv_comp = new_comp;

                                    head_flag = 1;
                                }

                                if(prv_flag == 1){
                                    prv_comp -> nxt_comp = new_comp;

                                    prv_comp = new_comp;
                                }*/
                                
                                new_comp -> comp_pwl_counter = pwl_counter - 1;
                                //new_comp -> (pwl_field + pwl_counter) = NULL;
                                break;
                                //continue;
                            }

                            ptr_token++;
                            
                            new_comp -> pwl_field = (pwl_spec)realloc((new_comp -> pwl_field), (sizeof(struct pwl_t) * (pwl_counter + 1)));

                            new_comp -> pwl_field[pwl_counter].pwl_start = strtod(ptr_token, NULL);

                            ptr_token = strtok(NULL, word);
                            strstr(ptr_token, ")");

                            new_comp -> pwl_field[pwl_counter].pwl_end = strtod(ptr_token, NULL);

                            pwl_counter++;
                        }

                    }

                    token_flag = 0;

                    if(head_flag == 0){
                        prv_comp = new_comp;

                        head_flag = 1;
                    }

                    if(prv_flag == 1){
                        prv_comp -> nxt_comp = new_comp;

                        prv_comp = new_comp;
                    }
                }
            }
            prv_flag = 1;

            component_counter++;

            free(dcsweep_line);
            free(plot_line);
            free(itol_line);
            free(plotnode_line);
            free(printnode_line);
            free(tran_line);
            free(tran_method_line);
            free(ACsweep_line);
        }
    }
    
    free(ptr_line);
        
    G1_components = component_counter - G2_components;

    tran_m_parameter = (tran_fin_time / tran_time_step);

    // default value for ITOL
    if(ITOL == 0){
        ITOL = 1e-3;
    }
}

void print_components(){
    ptr_comp print_ptr;

    print_ptr = list_head;

    printf("Number of components: %lu",component_counter);

    while(print_ptr != NULL){
        printf("\nType: %c Name: %s pos: %s neg: %s value: %lf vsource: %d\n",print_ptr -> comp_type, print_ptr -> comp_name, 
                                                                    print_ptr -> pos_node.pos_node_name,
                                                                    print_ptr -> neg_node.neg_node_name, print_ptr -> value,
                                                                    print_ptr -> vsource_pos);
        print_ptr = print_ptr -> nxt_comp;
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ht_dimensions(){
    ht_size = (component_counter);
    ht_depth = HASHTABLE_DEPTH;
}

void ht_creation(){
    int i, j;
    ht_dimensions();

    ht = (hashtable *)malloc(ht_size * sizeof(hashtable));

    if(!ht){
        printf("\nMemory allocation for Hashtabled failed.\n");
        exit(2);
    }

//  hashtable initialization
    for(i = 0; i < ht_size; i++){
        for(j = 0; j < ht_depth; j++){
            ((ht + i) -> old_name[j]) = NULL;

            ((ht + i) -> new_name_value[j]) = 0;
        }
    }
}


unsigned long hash(char* key){
    long hashval = 0;

    long g;

	while (*key != '\0') {
		hashval = (hashval << 4) + *(key++);

		g = hashval & 0xF0000000L;

		if(g != 0){
            hashval ^= g >> 24;
        }
		
        hashval &= ~g;
	}

	return(hashval);
}
void add_ht_oldname(){
    ptr_comp oldname_ptr;

    unsigned long posnode_key;
    unsigned long negnode_key;

    unsigned long posnode_index;
    unsigned long negnode_index;

    int i;

    int posoldname_exists = 1;
    int negoldname_exists = 1;

    oldname_ptr = list_head;

    ht_creation();
    
    while(oldname_ptr != NULL){
        posnode_key = hash(oldname_ptr -> pos_node.pos_node_name);
        negnode_key = hash(oldname_ptr -> neg_node.neg_node_name);

        posnode_index = posnode_key % ht_size;
        negnode_index = negnode_key % ht_size;
    
        for(i = 0; i < ht_depth; i++){
            if(((ht + posnode_index) -> old_name[i]) != NULL){
                posoldname_exists = strcmp(((ht + posnode_index) -> old_name[i]), oldname_ptr -> pos_node.pos_node_name);

                if(posoldname_exists == 0){
                    break;
                }
            }

            if(((ht + posnode_index) -> old_name[i]) == NULL){
                ((ht + posnode_index) -> old_name[i]) = strdup(oldname_ptr -> pos_node.pos_node_name);

                if(strcmp(((ht + posnode_index) -> old_name[i]), "0") == 0){
                    ((ht + posnode_index) -> new_name_value[i]) = 0;

                    break;
                }
                
                ((ht + posnode_index) -> new_name_value[i]) = curr_pos_oldname_table + 1;

                curr_pos_oldname_table++;

                break;
            }
        }

        for(i = 0; i < ht_depth; i++){
            if(((ht + negnode_index) -> old_name[i]) != NULL){
                negoldname_exists = strcmp(((ht + negnode_index) -> old_name[i]), oldname_ptr -> neg_node.neg_node_name);

                if(negoldname_exists == 0){
                    break;
                }
            }

            if(((ht + negnode_index) -> old_name[i]) == NULL){
                ((ht + negnode_index) -> old_name[i]) = strdup(oldname_ptr -> neg_node.neg_node_name);

                if(strcmp(((ht + negnode_index) -> old_name[i]), "0") == 0){
                    ((ht + negnode_index) -> new_name_value[i]) = 0;

                    break;
                }

                ((ht + negnode_index) -> new_name_value[i]) = curr_pos_oldname_table + 1;

                curr_pos_oldname_table++;

                break;
            }
        }

        oldname_ptr = oldname_ptr -> nxt_comp;
    }
}

void add_oldname_pos(){
    ptr_comp add_oldname_ptr;

    unsigned long pos_index, neg_index;
    unsigned long pos_key, neg_key;

    int i;

    int position;

    add_oldname_ptr = list_head;

    while(add_oldname_ptr != NULL){


        pos_key = hash(add_oldname_ptr -> pos_node.pos_node_name);

        pos_index = pos_key % ht_size;
        
        for(i = 0; i < ht_depth; i++){
            if(strcmp(((ht + pos_index) -> old_name[i]), (add_oldname_ptr) -> pos_node.pos_node_name) == 0){
                position = (ht + pos_index) -> new_name_value[i];

                add_oldname_ptr -> pos_node.new_pos_node_name = position;

                break;
            }
        }
        
        neg_key = hash(add_oldname_ptr -> neg_node.neg_node_name);

        neg_index = neg_key % ht_size;

        for(i = 0; i < ht_depth; i++){
            if(strcmp(((ht+neg_index) -> old_name[i]), (add_oldname_ptr) -> neg_node.neg_node_name) == 0){
                position = (ht + neg_index) -> new_name_value[i];

                add_oldname_ptr -> neg_node.new_neg_node_name = position;

                break;
            }
        }

        add_oldname_ptr = add_oldname_ptr -> nxt_comp;
    }
}

/****************************************************************/
void print_ht(){
    int i, j;

    for(i = 0; i < ht_size; i++){
        printf("\n[%d] old_name: ", i);

        for(j = 0; j < ht_depth; j++){
            if((ht+i) -> old_name[j] != NULL){
                printf(" [%s->%d]  ", (ht + i) -> old_name[j], (ht + i) -> new_name_value[j]);
            }
        }
        printf("\n");
    }
}

void print_updated_components(){
    ptr_comp print_ptr;

    print_ptr = list_head;

    printf("\n*********************************\n");

    printf("Number of nodes except ground: %d\n", curr_pos_oldname_table);

    printf("\nG1 components: %d\nG2 components: %d\n", G1_components, G2_components);


    while(print_ptr != NULL){
        printf("\nType: %c Name: %s pos: %d neg: %d value: %lf\n", print_ptr -> comp_type, print_ptr -> comp_name, 
                                                                    print_ptr -> pos_node.new_pos_node_name,
                                                                    print_ptr -> neg_node.new_neg_node_name, print_ptr -> value);
        print_ptr = print_ptr -> nxt_comp;
    }
}

void free_ht(){
    int i, j;

    for(i = 0; i < ht_size; i++){
        for(j = 0; j < ht_depth; j++){
            if(((ht + i) -> old_name[j]) != NULL){
                free((ht + i) -> old_name[j]);
            }
        }
    }
    free(ht);
}

void destroy_components(){
    ptr_comp destroy_ptr;

    free_ht();

    while(list_head != NULL){
        destroy_ptr = list_head;

        list_head = list_head -> nxt_comp;

        free(destroy_ptr -> comp_name);
        free(destroy_ptr -> pos_node.pos_node_name);
        free(destroy_ptr -> neg_node.neg_node_name);

        free(destroy_ptr -> pwl_field);

        free(destroy_ptr);
    }

    destroy_global_arrays(array_A, rowsA);

    destroy_global_vectors(array_b);

    destroy_global_vectors(array_x);

    free(DCsweep.dcsweep_component);

    free(plot_node);

    free(print_node);

    free_tran_MNA_matrices();

    free_BE_hands();

    free_TR_hands();

    free_tran_SPARSE_matrices();
}

char *return_new_filename(char *filename){
    char *fileptr;

    fileptr = strtok(filename, ".");

    strcat(fileptr, ".op");

    return(fileptr);
}

char *return_solution_filename(char *filename){
    char *fileptr;

    fileptr = strtok(filename, ".");

    strcat(fileptr, ".sol");

    return(fileptr);
}

void choose_method(char *newfilename){
    if(AC_flag == 1){
        choose_AC_method();

        return;
    }

    if(sparse_flag == 0){
        global_arrays_allocation();

        MNA_algorithm();

        choose_decomp(newfilename);

        choose_iterative_decomp(newfilename);

        choose_tran_solution_method();

        choose_tran_method();

        if((plot_flag == 1) && (tran_flag == 0)){
            dc_sweep();
        }
    }
    else{
        triplet_format_method();

        choose_sparse_decomp(newfilename);

        choose_sparse_iterative_decomp(newfilename);

        choose_tran_solution_method();

        choose_tran_method();

        if((plot_flag == 1) && (tran_flag == 0)){
            sparse_dc_sweep();
        }
    }
}