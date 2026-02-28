#include "ctools.h"
int get_file_lines_LinuxSys(char * filename)
{
    char command[500];
    sprintf(command,"wc -l < %s",filename);
    FILE * fp=popen(command,"r");
    if(fp==NULL)
    {
        perror("popen failed");
    }
    int lines=0;
    if(fscanf(fp,"%d",&lines) != 1)
    {
        printf("Err! function get_file_lines_LinuxSys err!\n");
    }
    pclose(fp);
    
    return lines;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////functions 
void INPUT_free(struct PARAMETERS * input)
{
    struct MODE * tmp=input->modes;
    struct MODE * tmp2=NULL;
    while(tmp!=NULL)
    {
        tmp2=tmp->next;
        DATA_free(tmp->descripter);
        free(tmp);
        tmp=tmp2;
        
    }
    /////////////////////////////////////
    solidfication_free(input->so_check);
    /////////////////////////////////////
    
    free(input);
}
void PARA_toFILE(int method)
{
        FILE * fp=fopen("input.dat","w");

        if(method == 0) //single phase method
        {
            fprintf(fp,"method 1 \n");
            fprintf(fp,"solid_cutoff 0.8\n");
            fprintf(fp,"liquid_cutoff 0.2\n\n");

            fprintf(fp,"md_script_interpreter \n");
            fprintf(fp,"md_script_filename \n");
            fprintf(fp,"md_steps 2000\n");
            fprintf(fp,"md_steps_adaption 0\n\n");

            fprintf(fp,"run heatingup 500 1000 0.1\n");
            fprintf(fp,"max_step 20 \n");
            fprintf(fp,"max_sub_step 40 \n");
            fprintf(fp,"init_tem_check 1\n");
            fprintf(fp,"bisection_relax_struct 0\n\n");

            fprintf(fp,"supercell 5 5 5\n");
        }
        else if(method == 1) //void method
        {
            fprintf(fp,"method 1 \n");
            fprintf(fp,"solid_cutoff 0.8\n");
            fprintf(fp,"liquid_cutoff 0.2\n\n");

            fprintf(fp,"md_script_interpreter \n");
            fprintf(fp,"md_script_filename \n");
            fprintf(fp,"md_steps 2000\n");
            fprintf(fp,"md_steps_adaption 0\n\n");

            fprintf(fp,"run heatingup 500 1000 0.1\n");
            fprintf(fp,"max_step 20 \n");
            fprintf(fp,"max_sub_step 40 \n");
            fprintf(fp,"init_tem_check 1\n");
            fprintf(fp,"bisection_relax_struct 0\n\n");

            fprintf(fp,"supercell 5 5 5\n");
            fprintf(fp,"build defect 1 0 1.0\n");
            fprintf(fp,"size 0.5 0.5 0.5\n");
            fprintf(fp,"pos 0.5 0.5 0.5\n");
            fprintf(fp,"rot 0 0 0\n\n");

        }
        else if(method == 2) //modified void method
        {
            fprintf(fp,"method 1 \n");
            fprintf(fp,"solid_cutoff 0.8\n");
            fprintf(fp,"liquid_cutoff 0.2\n\n");

            fprintf(fp,"md_script_interpreter \n");
            fprintf(fp,"md_script_filename \n");
            fprintf(fp,"md_steps 20000\n");
            fprintf(fp,"md_steps_adaption 1\n\n");

            fprintf(fp,"run bisection 100 1000 2 0.1\n");
            fprintf(fp,"max_step 20 \n");
            fprintf(fp,"max_sub_step 40 \n");
            fprintf(fp,"bisection_relax_struct 0\n");
            fprintf(fp,"init_tem_check 1\n\n");


            fprintf(fp,"supercell 5 5 20\n");
            fprintf(fp,"build defect 1 0 0.8\n");
            fprintf(fp,"size 1.0 1.0 0.15\n");
            fprintf(fp,"pos 0.5 0.5 0.5\n");
            fprintf(fp,"rot 0 0 0\n\n");
                        
        }
        else if(method == 3) //two phase method
        {
            fprintf(fp,"method 1 \n");
            fprintf(fp,"solid_cutoff 0.8\n");
            fprintf(fp,"liquid_cutoff 0.2\n\n");

            fprintf(fp,"md_script_interpreter \n");
            fprintf(fp,"md_script_filename \n");
            fprintf(fp,"md_steps 2000\n");
            fprintf(fp,"md_steps_adaption 0\n\n");

            fprintf(fp,"run bisection 100 1000 1 100\n");
            fprintf(fp,"max_step 20 \n");
            fprintf(fp,"max_sub_step 40 \n");
            fprintf(fp,"bisection_relax_struct 1\n");
            fprintf(fp,"init_tem_check 1\n\n");

            fprintf(fp,"supercell 5 5 10\n");
            fprintf(fp,"build liquid 1 0\n");
            fprintf(fp,"size 1.0 1.0 0.5\n");
            fprintf(fp,"pos 0.5 0.5 0.25\n");
            fprintf(fp,"rot 0 0 0\n\n");

                      
        }
        else if(method == 4) //sanwitch method
        {
            fprintf(fp,"method 1 \n");
            fprintf(fp,"solid_cutoff 0.8\n");
            fprintf(fp,"liquid_cutoff 0.2\n\n");

            fprintf(fp,"md_script_interpreter \n");
            fprintf(fp,"md_script_filename \n");
            fprintf(fp,"md_steps 2000\n");
            fprintf(fp,"md_steps_adaption 0\n\n");

            fprintf(fp,"run bisection 100 1000 1 100\n");
            fprintf(fp,"max_step 20 \n");
            fprintf(fp,"max_sub_step 40 \n");
            fprintf(fp,"bisection_relax_struct 1\n");
            fprintf(fp,"init_tem_check 1\n\n");

            fprintf(fp,"supercell 4 4 12\n");
            fprintf(fp,"build liquid 1 0\n");
            fprintf(fp,"size 1.0 1.0 0.25\n");
            fprintf(fp,"pos 0.5 0.5 0.125\n");
            fprintf(fp,"rot 0 0 0\n");     
            fprintf(fp,"build liquid 1 0\n");
            fprintf(fp,"size 1.0 1.0 0.25\n");
            fprintf(fp,"pos 0.5 0.5 0.625\n");
            fprintf(fp,"rot 0 0 0\n\n");    
        }
        else if(method == 5) //Z method
        {
            fprintf(fp,"method 1 \n");
            fprintf(fp,"solid_cutoff 0.8\n");
            fprintf(fp,"liquid_cutoff 0.2\n\n");

            fprintf(fp,"md_script_interpreter \n");
            fprintf(fp,"md_script_filename \n");
            fprintf(fp,"md_steps 50000\n");
            fprintf(fp,"md_steps_adaption 1\n\n");

            fprintf(fp,"run bisection 100 1000 1 10\n");
            fprintf(fp,"max_step 40 \n");
            fprintf(fp,"max_sub_step 40 \n");
            fprintf(fp,"bisection_relax_struct 0\n");
            fprintf(fp,"init_tem_check 1\n\n");

            fprintf(fp,"supercell 5 5 5\n");         
        }
        else if(method == 6) //modified Z method
        {
            fprintf(fp,"method 1 \n");
            fprintf(fp,"solid_cutoff 0.7\n");
            fprintf(fp,"liquid_cutoff 0.3\n\n");

            fprintf(fp,"md_script_interpreter \n");
            fprintf(fp,"md_script_filename \n");
            fprintf(fp,"md_steps 50000\n");
            fprintf(fp,"md_steps_adaption 1\n\n");

            fprintf(fp,"run bisection 100 1000 2 0.1\n");
            fprintf(fp,"max_step 40 \n");
            fprintf(fp,"max_sub_step 40 \n");
            fprintf(fp,"bisection_relax_struct 0\n");
            fprintf(fp,"init_tem_check 1\n\n");

            fprintf(fp,"supercell 5 5 5\n"); 
        }

    fclose(fp);
}
struct PARAMETERS * read_parameters_from_file(char * filename)
{
    struct PARAMETERS *inputs=(struct PARAMETERS *)calloc(1,sizeof(struct PARAMETERS));

    inputs->vasp_NSW=1000;
    inputs->A=1;
    inputs->B=1;
    inputs->C=1;
    inputs->nsw_adaption_ratio=0;
    inputs->decision=0.8;
    inputs->decision_liquid=0.2;
    inputs->modes=NULL;
    inputs->max_step=20;
    inputs->max_sub_step=20;
    inputs->atom_method[0]=1;
    inputs->atom_method[1]=6;
    inputs->atom_method[1]=0.5;
    inputs->atom_method[1]=0.5;
    inputs->atom_method[1]=4.0;


    inputs->modes=NULL;
    inputs->so_check=NULL;
    //////////////////////////////////////////////////////
    FILE * fp=fopen(filename,"r");

    if(fp==NULL)
    {
        printf("Err! read file input.dat err!\n");
    }

    struct STRINGS * str;
    char buff[1000];
    int state=1;
    if(fgets(buff,sizeof(buff),fp)==NULL)
    {
        state=-100;
    }
    str=STRINGS_in(buff,100,100);
    int check_state=0;
    while(state)
    {

            if(!strcmp(STRINGS_get(str,0),"md_steps") && str->number_words ==2)
            {
                inputs->vasp_NSW=atoi(STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }
            else if(!strcmp(STRINGS_get(str,0),"md_steps_adaption") && str->number_words ==2)
            {
                inputs->nsw_adaption_ratio=atof(STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }
            else if(!strcmp(STRINGS_get(str,0),"solid_cutoff") && str->number_words ==2)
            {
                inputs->decision=atof(STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }
            else if(!strcmp(STRINGS_get(str,0),"liquid_cutoff") && str->number_words ==2)
            {
                inputs->decision_liquid=atof(STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }
            else if(!strcmp(STRINGS_get(str,0),"method") )
            {
                int methods=atoi(STRINGS_get(str,1));
                if(methods < 0)
                {
                    if(methods == -1 && str->number_words == 6)
                    {
                        inputs->atom_method[0]=1.0*methods;
                        inputs->atom_method[1]=atof(STRINGS_get(str,2));
                        inputs->atom_method[2]=atof(STRINGS_get(str,3));
                        inputs->atom_method[3]=atof(STRINGS_get(str,4));
                        inputs->atom_method[4]=atof(STRINGS_get(str,5));
                       
                    }
                    else if(methods == -1 && str->number_words != 6)
                    {
                        inputs->atom_method[0]=-1.0*methods;
                    }
                    else if(methods == -2 && str->number_words == 6)
                    {
                        inputs->atom_method[0]=1.0*methods;
                        inputs->atom_method[1]=atof(STRINGS_get(str,2));
                        inputs->atom_method[2]=atof(STRINGS_get(str,3));
                        inputs->atom_method[3]=atof(STRINGS_get(str,4));
                        inputs->atom_method[4]=atof(STRINGS_get(str,5));
                       
                    }
                    else if(methods == -2 && str->number_words != 6)
                    {
                        inputs->atom_method[0]=-1.0*methods;
                    }
                    else if(methods == -3 && str->number_words == 9)
                    {
                        inputs->atom_method[0]=1.0*methods;
                        inputs->atom_method[1]=atof(STRINGS_get(str,2));
                        inputs->atom_method[2]=atof(STRINGS_get(str,3));
                        inputs->atom_method[3]=atof(STRINGS_get(str,4));
                        inputs->atom_method[4]=atof(STRINGS_get(str,5));
                        inputs->atom_method[5]=atof(STRINGS_get(str,6));
                        inputs->atom_method[6]=atof(STRINGS_get(str,7));
                        inputs->atom_method[7]=atof(STRINGS_get(str,8));
                       
                    }
                    else if(methods == -3 && str->number_words != 9)
                    {
                        inputs->atom_method[0]=-1.0*methods;
                    }

                }
                else if(methods >0)
                {

                    inputs->atom_method[0]=1.0*methods;
                    
                }

                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
                
            }
            else if(!strcmp(STRINGS_get(str,0),"supercell") )
            {
                if(str->number_words == 4)
                {

                    inputs->A=atoi(STRINGS_get(str,1));
                    inputs->B=atoi(STRINGS_get(str,2));
                    inputs->C=atoi(STRINGS_get(str,3));
                }
                else
                {
                    inputs->A=1;
                    inputs->B=1;
                    inputs->C=1;
                }

                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
                
            }
            else if(!strcmp(STRINGS_get(str,0),"max_step") && str->number_words ==2)
            {
                inputs->max_step=atoi(STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }
            else if(!strcmp(STRINGS_get(str,0),"max_sub_step") && str->number_words ==2)
            {
                inputs->max_sub_step=atoi(STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }
            else if(!strcmp(STRINGS_get(str,0),"md_script_interpreter") && str->number_words ==2)
            {
                strcpy(inputs->commands,STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }   
            else if(!strcmp(STRINGS_get(str,0),"md_script_filename") && str->number_words ==2)
            {
                strcpy(inputs->script,STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }           
            else if(!strcmp(STRINGS_get(str,0),"init_tem_check") && str->number_words ==2)
            {
                inputs->init_tem_check=atoi(STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
                
            }            
            else if(!strcmp(STRINGS_get(str,0),"bisection_relax_struct") && str->number_words ==2)
            {
                inputs->dichotomy_structure_relax=atoi(STRINGS_get(str,1));
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
                
            }  
            else if(!strcmp(STRINGS_get(str,0),"build") && str->number_words >2)
            {
                struct MODE * tmp_mode=inputs->modes;
                if(  !strcmp(STRINGS_get(str,1),"defect")    )
                {
                    if(tmp_mode == NULL)
                    {
                        inputs->modes=(struct MODE *)calloc(1,sizeof(struct MODE));
                        inputs->modes->next=NULL;
                        tmp_mode=inputs->modes;

                    }
                    else
                    {
                        while(tmp_mode->next != NULL )
                        {
                            tmp_mode=inputs->modes->next;
                        }
                        tmp_mode->next=(struct MODE *)calloc(1,sizeof(struct MODE));
                        tmp_mode=tmp_mode->next;

                    }
                    tmp_mode->descripter=DATA2D_init(4,6);
                    //////////////////////////////////////////////////
                    //////////////////////////////////////////////////

                    DATA2D_set(tmp_mode->descripter,0,0,0); // type
                    DATA2D_set(tmp_mode->descripter,atoi(STRINGS_get(str,2)),0,1); // build scale
                    DATA2D_set(tmp_mode->descripter,atoi(STRINGS_get(str,3)),0,2); // build scale
                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,4)),0,3); // build scale
                    int sub_state=1;
                    while(sub_state)
                    {
                        STRINGS_free(str);
                        if(fgets(buff,sizeof(buff),fp)==NULL)
                        {
                            state=0;
                        }
                        else
                        {
                            str=STRINGS_in(buff,100,100);
                            if(!strcmp(STRINGS_get(str,0),"size"))
                            {
                                if(str->number_words == 4)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),1,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),1,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),1,2);
                                }
                                else if(str->number_words == 7)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),1,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),1,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),1,2);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,4)),1,3);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,5)),1,4);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,6)),1,5);
                                }
                                else
                                {
                                    DATA2D_set(tmp_mode->descripter,0,1,0);
                                    DATA2D_set(tmp_mode->descripter,0,1,1);
                                    DATA2D_set(tmp_mode->descripter,0,1,2);
                                    DATA2D_set(tmp_mode->descripter,0,1,3);
                                    DATA2D_set(tmp_mode->descripter,0,1,4);
                                    DATA2D_set(tmp_mode->descripter,0,1,5);
                                }

                            }
                            else if(!strcmp(STRINGS_get(str,0),"pos"))
                            {
                                if(str->number_words == 4)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),2,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),2,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),2,2);
                                }
                                else if(str->number_words == 7)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),2,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),2,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),2,2);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,4)),2,3);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,5)),2,4);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,6)),2,5);
                                }
                                else
                                {
                                    DATA2D_set(tmp_mode->descripter,0,2,0);
                                    DATA2D_set(tmp_mode->descripter,0,2,1);
                                    DATA2D_set(tmp_mode->descripter,0,2,2);
                                    DATA2D_set(tmp_mode->descripter,0,2,3);
                                    DATA2D_set(tmp_mode->descripter,0,2,4);
                                    DATA2D_set(tmp_mode->descripter,0,2,5);
                                }                                
                            }
                            else if(!strcmp(STRINGS_get(str,0),"rot"))
                            {
                                if(str->number_words == 4)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),3,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),3,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),3,2);
                                }
                                else if(str->number_words == 7)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),3,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),3,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),3,2);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,4)),3,3);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,5)),3,4);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,6)),3,5);
                                }
                                else
                                {
                                    DATA2D_set(tmp_mode->descripter,0,3,0);
                                    DATA2D_set(tmp_mode->descripter,0,3,1);
                                    DATA2D_set(tmp_mode->descripter,0,3,2);
                                    DATA2D_set(tmp_mode->descripter,0,3,3);
                                    DATA2D_set(tmp_mode->descripter,0,3,4);
                                    DATA2D_set(tmp_mode->descripter,0,3,5);
                                } 
                            }
                            else
                            {
                                sub_state=0;
                            }
                        }
                        
                    }                  


                }
                else if(  !strcmp(STRINGS_get(str,1),"liquid")    )
                {
                    inputs->isliquid=1;
                    if(tmp_mode == NULL)
                    {
                        inputs->modes=(struct MODE *)calloc(1,sizeof(struct MODE));
                        inputs->modes->next=NULL;
                        tmp_mode=inputs->modes;

                    }
                    else
                    {
                        while(tmp_mode->next != NULL )
                        {
                            tmp_mode=inputs->modes->next;
                        }
                        tmp_mode->next=(struct MODE *)calloc(1,sizeof(struct MODE));
                        tmp_mode=tmp_mode->next;

                    }
                    tmp_mode->descripter=DATA2D_init(4,6);
                    //////////////////////////////////////////////////
                    //////////////////////////////////////////////////

                    DATA2D_set(tmp_mode->descripter,1,0,0); // type
                    DATA2D_set(tmp_mode->descripter,atoi(STRINGS_get(str,2)),0,1); // build scale
                    DATA2D_set(tmp_mode->descripter,atoi(STRINGS_get(str,3)),0,2); // build scale
                    int sub_state=1;
                    while(sub_state)
                    {
                        STRINGS_free(str);
                        if(fgets(buff,sizeof(buff),fp)==NULL)
                        {
                            state=0;
                        }
                        else
                        {
                            str=STRINGS_in(buff,100,100);
                            if(!strcmp(STRINGS_get(str,0),"size"))
                            {
                                if(str->number_words == 4)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),1,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),1,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),1,2);
                                }
                                else if(str->number_words == 7)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),1,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),1,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),1,2);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,4)),1,3);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,5)),1,4);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,6)),1,5);
                                }
                                else
                                {
                                    DATA2D_set(tmp_mode->descripter,0,1,0);
                                    DATA2D_set(tmp_mode->descripter,0,1,1);
                                    DATA2D_set(tmp_mode->descripter,0,1,2);
                                    DATA2D_set(tmp_mode->descripter,0,1,3);
                                    DATA2D_set(tmp_mode->descripter,0,1,4);
                                    DATA2D_set(tmp_mode->descripter,0,1,5);
                                }

                            }
                            else if(!strcmp(STRINGS_get(str,0),"pos"))
                            {
                                if(str->number_words == 4)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),2,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),2,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),2,2);
                                }
                                else if(str->number_words == 7)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),2,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),2,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),2,2);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,4)),2,3);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,5)),2,4);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,6)),2,5);
                                }
                                else
                                {
                                    DATA2D_set(tmp_mode->descripter,0,2,0);
                                    DATA2D_set(tmp_mode->descripter,0,2,1);
                                    DATA2D_set(tmp_mode->descripter,0,2,2);
                                    DATA2D_set(tmp_mode->descripter,0,2,3);
                                    DATA2D_set(tmp_mode->descripter,0,2,4);
                                    DATA2D_set(tmp_mode->descripter,0,2,5);
                                }                                
                            }
                            else if(!strcmp(STRINGS_get(str,0),"rot"))
                            {
                                if(str->number_words == 4)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),3,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),3,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),3,2);
                                }
                                else if(str->number_words == 7)
                                {
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,1)),3,0);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,2)),3,1);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,3)),3,2);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,4)),3,3);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,5)),3,4);
                                    DATA2D_set(tmp_mode->descripter,atof(STRINGS_get(str,6)),3,5);
                                }
                                else
                                {
                                    DATA2D_set(tmp_mode->descripter,0,3,0);
                                    DATA2D_set(tmp_mode->descripter,0,3,1);
                                    DATA2D_set(tmp_mode->descripter,0,3,2);
                                    DATA2D_set(tmp_mode->descripter,0,3,3);
                                    DATA2D_set(tmp_mode->descripter,0,3,4);
                                    DATA2D_set(tmp_mode->descripter,0,3,5);
                                } 
                            }
                            else
                            {
                                sub_state=0;
                            }
                        }
                        
                    }                    

                }
                else
                {   
                    STRINGS_free(str);
                    if(fgets(buff,sizeof(buff),fp)==NULL)
                    {
                        state=0;
                    }
                    else str=STRINGS_in(buff,100,100);
                }

                
            }
            else if(!strcmp(STRINGS_get(str,0),"run") && str->number_words >=2)
            {
                if(!strcmp(STRINGS_get(str,1),"smallcell")  && str->number_words == 4)
                {
                    check_state=1;
                    inputs->solid_tem=atof(STRINGS_get(str,2));
                    inputs->melt_tem=atof(STRINGS_get(str,3));
                    STRINGS_free(str);
                    if(fgets(buff,sizeof(buff),fp)==NULL)
                    {
                        state=0;
                    }
                    else str=STRINGS_in(buff,100,100);
                    inputs->method=3;
                }
                else if(!strcmp(STRINGS_get(str,1),"heatingup") && str->number_words == 5 )
                {
                    check_state=1;
                    inputs->solid_tem=atof(STRINGS_get(str,2));
                    inputs->melt_tem=atof(STRINGS_get(str,3));
                    inputs->heat_ratio=atof(STRINGS_get(str,4));
                    STRINGS_free(str);
                    if(fgets(buff,sizeof(buff),fp)==NULL)
                    {
                        state=0;
                    }
                    else str=STRINGS_in(buff,100,100);
                    inputs->method=0;
                }
                else if(!strcmp(STRINGS_get(str,1),"bisection") && str->number_words == 6)
                {
                    check_state=1;
                    inputs->solid_tem=atof(STRINGS_get(str,2));
                    inputs->melt_tem=atof(STRINGS_get(str,3));
                    inputs->method=atoi(STRINGS_get(str,4));
                    inputs->prec_tem=atof(STRINGS_get(str,5));
                    if(inputs->method != 1 && inputs->method != 2)
                    {
                        printf("Err! bisection temperature control method err! the exit condition err!\n");
                        return NULL;
                    }
                    
                    STRINGS_free(str);
                    if(fgets(buff,sizeof(buff),fp)==NULL)
                    {
                        state=0;
                    }
                    else str=STRINGS_in(buff,100,100);
                }
                else if(!strcmp(STRINGS_get(str,1),"test") )
                {
                    check_state=1;
                    inputs->method=-1;
                    
                    STRINGS_free(str);
                    if(fgets(buff,sizeof(buff),fp)==NULL)
                    {
                        state=0;
                    }
                    else str=STRINGS_in(buff,100,100);
                }                
                else
                {
                    printf("Err! input.dat method err!\n");
                    INPUT_free(inputs);
                    state=0;
                    inputs=NULL;
                    break;
                }

            }
            else
            {
                STRINGS_free(str);
                if(fgets(buff,sizeof(buff),fp)==NULL)
                {
                    state=0;
                }
                else str=STRINGS_in(buff,100,100);
            }
    }

    if(check_state==0)
    {
        printf("err! function temperature control method cannt be found!\n");
        free(inputs);
        inputs=NULL;
    }
    fclose(fp);

    return inputs;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////functions for data structure
struct DATA ** DATAs_init(int numbers,int dimens, ...)
{
    int number_data=1;
    va_list arg;
    va_start(arg,dimens);
    int *args=NULL;
    if(dimens > 0)
    {
        args=(int * )calloc(dimens,sizeof(int));
        for(int i=0;i<dimens;i++)
        {
            args[i]=va_arg(arg,int);
            number_data=number_data*args[i];
        }
    }
    va_end(arg);
    //////////////////////////////////////////////////////////////////
    struct DATA ** tmp_out=NULL;
    if(numbers > 0 && number_data > 0)
    {
        tmp_out=(struct DATA **)calloc(numbers,sizeof(struct DATA *));
        
        for(int i=0;i<numbers;i++)
        {
            tmp_out[i]=(struct DATA *)calloc(1,sizeof(struct DATA));
            tmp_out[i]->dimen=dimens;
            tmp_out[i]->data_group=numbers;
            tmp_out[i]->all_data=number_data;
            
            tmp_out[i]->dimen_data=(int * )calloc(dimens,sizeof(int));
            memcpy(tmp_out[i]->dimen_data,args,dimens*sizeof(int));
            tmp_out[i]->data_data=(double *)calloc(number_data,sizeof(double));
            
        }
    }
    return tmp_out;
}
double DATA_get(struct DATA*data, ...)
{
    if((data->all_data) > 0)
    {
        int dim=data->dimen;
        va_list arg;
        va_start(arg,data);
        int *args=(int *)calloc(dim,sizeof(int));
        for(int i=0;i<dim;i++)
        {
            args[i]=va_arg(arg,int);
        }
        va_end(arg);
        //////////////////////////////////
        int pos=0;
        int pts=1;
        for(int i=0;i<dim;i++)
        {
            if(i==0)
            {
                pts=1;
            }
            else
            {
                pts=pts*(data->dimen_data[i-1]);
            }
            pos+=pts*args[i];
        }
        /////////////////////////////////
        free(args);
        return data->data_data[pos];
    }
    return 0;
}
void DATA_set(struct DATA*data,double set_double, ...)
{
    if((data->all_data) >0)
    {
        int dim=data->dimen;
        va_list arg;
        va_start(arg,set_double);
        int *args=(int *)calloc(dim,sizeof(int));
        for(int i=0;i<dim;i++)
        {
            args[i]=va_arg(arg,int);
        }
        va_end(arg);
        //////////////////////////////////
        int pos=0;
        int pts=1;
        for(int i=0;i<dim;i++)
        {
            if(i==0)
            {
                pts=1;
            }
            else
            {
                pts=pts*(data->dimen_data[i-1]);
            }
            pos+=pts*args[i];
        }
        /////////////////////////////////
        free(args);
        data->data_data[pos]=set_double;
    }
}
void DATA_free(struct DATA *data)
{
    if(data != NULL)
    {

        free(data->dimen_data);

        free(data->data_data);

        free(data);

    }
}
struct DATA * DATA_init(int number_int, ...)
{

    int num_datas=1;
    va_list arg;
    va_start(arg,number_int);
    int * args=NULL;
    if(number_int > 0 ) 
    {
        args=(int *)calloc(number_int,sizeof(int));
        for(int i=0;i<number_int;i++)
        {
            args[i]=va_arg(arg,int);
            num_datas=num_datas*args[i];
        }
    }
    va_end(arg);
    ////////////////////////////////////////////////////////////////////////////
    struct DATA* tmp_data=(struct DATA*)calloc(1,sizeof(struct DATA));

    if(num_datas > 0 && number_int > 0)
    {
        tmp_data->all_data=num_datas;
        tmp_data->data_data=(double *)calloc(num_datas,sizeof(double));
        tmp_data->dimen_data=args;
        tmp_data->dimen=number_int;
    }
    else
    {
        tmp_data->all_data=0;
        tmp_data->data_data=NULL;
        tmp_data->dimen_data=NULL;
        tmp_data->dimen=0;
    }
    return tmp_data;

}
void DATAs_free(struct DATA ** datas)
{
    int numbers=datas[0]->data_group;
    for(int i=0;i<numbers;i++)
    {
        free(datas[i]->data_data);
        free(datas[i]->dimen_data);
        free(datas[i]);
    }
    free(datas);
}
void DATA_getXYT(struct DATA * Tdata,int sel,double *x,double *y,double *z)
{
    if(Tdata->dimen==2 && Tdata->dimen_data[1]==3)
    {
        *x=DATA_get(Tdata,sel,0);
        *y=DATA_get(Tdata,sel,1);
        *z=DATA_get(Tdata,sel,2);

    }
    else
    {
        *x=-1;
        *y=-1;
        *z=-1;
    }
}
void DATA_setXYT(struct DATA * Tdata,int sel,double x,double y,double z)
{
    if(Tdata->dimen==2 && Tdata->dimen_data[1]==3)
    {
        DATA_set(Tdata,x,sel,0);
        DATA_set(Tdata,y,sel,1);
        DATA_set(Tdata,z,sel,2);
    }
}
struct DATA * DATA_cpy(struct DATA * org)
{
    struct DATA * tnp=NULL;
    if(org != NULL)
    {
        tnp=(struct DATA *)malloc(sizeof(struct DATA));
        tnp->dimen=org->dimen;
        tnp->data_group=org->data_group;
        tnp->all_data=org->all_data;
        if(org->dimen_data != NULL)
        {
             tnp->dimen_data=(int *)calloc(tnp->dimen,sizeof(int));
             memcpy(tnp->dimen_data,org->dimen_data,(tnp->dimen)*sizeof(int));
        }
        if(org->data_data !=NULL)
        {
            tnp->data_data=(double *)calloc(tnp->all_data,sizeof(double));
            memcpy(tnp->data_data,org->data_data,(tnp->all_data)*sizeof(double));
        }
    }
    return tnp;

}
struct DATA *DATA_shapeINIT_rect(int inside,int coor_type,double x0,double y0,double z0,double x1,double y1,double z1)
{
    struct DATA * out=DATA_init(2,3,3);
    DATA_set(out,1.0,0,0);
    DATA_set(out,1.0*inside,0,1);
    DATA_set(out,1.0*coor_type,0,2);

    DATA_set(out,x0,1,0);
    DATA_set(out,y0,1,1);
    DATA_set(out,z0,1,2);

    DATA_set(out,x1,2,0);
    DATA_set(out,y1,2,1);
    DATA_set(out,z1,2,2);

    return out;
}
double DATA4D_get(struct DATA * data,int dime1,int dime2,int dime3,int dime4)
{
    return data->data_data[dime4*(data->dimen_data[2])*(data->dimen_data[1])*(data->dimen_data[0]) + dime3*(data->dimen_data[1]) +dime2 * (data->dimen_data[0]) + dime1 ];
}
void DATA4D_set(struct DATA * data,double set_data, int dime1,int dime2,int dime3,int dime4)
{
    data->data_data[dime4*(data->dimen_data[2])*(data->dimen_data[1])*(data->dimen_data[0]) + dime3*(data->dimen_data[1]) +dime2 * (data->dimen_data[0]) + dime1 ]=set_data;
}
double DATA3D_get(struct DATA * data,int dime1,int dime2,int dime3)
{
    return data->data_data[ dime3*(data->dimen_data[1])*(data->dimen_data[0]) + dime2 * (data->dimen_data[0]) + dime1 ];
}
void DATA3D_set(struct DATA * data,double set_data, int dime1,int dime2,int dime3)
{
    data->data_data[ dime3*(data->dimen_data[1])*(data->dimen_data[0]) + dime2 * (data->dimen_data[0]) + dime1 ]=set_data;
}
double DATA2D_get(struct DATA * data,int dime1,int dime2)
{
    return data->data_data[  dime2 * (data->dimen_data[0]) + dime1     ];
}
void DATA2D_set(struct DATA * data,double set_data, int dime1,int dime2)
{
    data->data_data[dime2 * (data->dimen_data[0]) + dime1 ]=set_data;
}
double DATA1D_get(struct DATA * data,int dime1)
{
    return data->data_data[   dime1     ];
}
void DATA1D_set(struct DATA * data,double set_data, int dime1)
{
    data->data_data[ dime1 ]=set_data;
}

struct DATA * DATA1D_init(int dimen1)
{

    struct DATA* tmp_data=(struct DATA*)calloc(1,sizeof(struct DATA));

    if(dimen1 > 0)
    {
        tmp_data->all_data=dimen1;
        tmp_data->data_data=(double *)calloc(dimen1,sizeof(double));
        tmp_data->dimen_data=(int *)calloc(1,sizeof(int));
        tmp_data->dimen_data[0]=dimen1;
        tmp_data->dimen=1;
    }
    else
    {
        tmp_data->all_data=0;
        tmp_data->data_data=NULL;
        tmp_data->dimen_data=NULL;
        tmp_data->dimen=0;
    }
    return tmp_data;

}
struct DATA * DATA2D_init(int dimen1,int dimen2)
{

    struct DATA* tmp_data=(struct DATA*)calloc(1,sizeof(struct DATA));

    if(dimen1 > 0 && dimen2 >0)
    {
        tmp_data->all_data=dimen1*dimen2;
        tmp_data->data_data=(double *)calloc(dimen1*dimen2,sizeof(double));
        tmp_data->dimen_data=(int *)calloc(2,sizeof(int));
        tmp_data->dimen_data[0]=dimen1;
        tmp_data->dimen_data[1]=dimen2;
        tmp_data->dimen=2;
    }
    else
    {
        tmp_data->all_data=0;
        tmp_data->data_data=NULL;
        tmp_data->dimen_data=NULL;
        tmp_data->dimen=0;
    }
    return tmp_data;

}
struct DATA * DATA3D_init(int dimen1,int dimen2,int dimen3)
{

    struct DATA* tmp_data=(struct DATA*)calloc(1,sizeof(struct DATA));

    if(dimen1 > 0 && dimen2 >0 && dimen3 > 0)
    {
        tmp_data->all_data=dimen1*dimen2*dimen3;
        tmp_data->data_data=(double *)calloc(dimen1*dimen2*dimen3,sizeof(double));
        tmp_data->dimen_data=(int *)calloc(3,sizeof(int));
        tmp_data->dimen_data[0]=dimen1;
        tmp_data->dimen_data[1]=dimen2;
        tmp_data->dimen_data[2]=dimen3;
        tmp_data->dimen=3;
    }
    else
    {
        tmp_data->all_data=0;
        tmp_data->data_data=NULL;
        tmp_data->dimen_data=NULL;
        tmp_data->dimen=0;
    }
    return tmp_data;

}
struct DATA * DATA4D_init(int dimen1,int dimen2,int dimen3,int dimen4)
{

    struct DATA* tmp_data=(struct DATA*)calloc(1,sizeof(struct DATA));

    if(dimen1 > 0 && dimen2 >0 && dimen3 > 0 && dimen4>0)
    {
        tmp_data->all_data=dimen1*dimen2*dimen3*dimen4;
        tmp_data->data_data=(double *)calloc(dimen1*dimen2*dimen3*dimen4,sizeof(double));
        tmp_data->dimen_data=(int *)calloc(4,sizeof(int));
        tmp_data->dimen_data[0]=dimen1;
        tmp_data->dimen_data[1]=dimen2;
        tmp_data->dimen_data[2]=dimen3;
        tmp_data->dimen_data[2]=dimen4;
        tmp_data->dimen=4;
    }
    else
    {
        tmp_data->all_data=0;
        tmp_data->data_data=NULL;
        tmp_data->dimen_data=NULL;
        tmp_data->dimen=0;
    }
    return tmp_data;

}

struct DATA * XY_in(char *file_name)
{
    struct DATA * outputs=NULL;
    FILE *fp=fopen(file_name,"r");
    if(fp!=NULL)
    {
        char tmps[10000];
        fgets(tmps,10000,fp);
        struct STRINGS * lons=STRINGS_in(tmps,100,100);

        int row;
        row=lons->number_words;
        STRINGS_free(lons);
        fseek(fp,0,SEEK_SET);
        int col=get_file_lines_LinuxSys(file_name);
        ///////////////////////////
        outputs=DATA2D_init(col,row);
        for(int i=0;i<col;i++)
        {
            for(int j=0;j<row;j++)
            {
                double res=-1;
                fscanf(fp,"%lf",&res);
                DATA2D_set(outputs,res,i,j);
            }

        }
        fclose(fp);
    }
    else
    {
        printf("%s   XY_in function read err!\n",file_name);
        outputs=NULL;
    }
    return outputs;
}
void XY_toFILE(struct DATA * output,char *file_name)
{
    FILE *fp=fopen(file_name,"w");
    if(fp!=NULL)
    {
        int cow,row;
        cow=output->dimen_data[0];
        row=output->dimen_data[1];
        if(cow>0 & row >0)
        {
            for(int i=0;i<cow;i++)
            {
                double tmpdata;
               
                for(int j=0;j<row;j++)
                {
                     tmpdata=DATA_get(output,i,j);
                     fprintf(fp,"%lf\t",tmpdata);
                } 
                 fprintf(fp,"\n");
            }
        }
        else
        {
            printf("XY_toFILE function %d %d err!\n",cow,row);
        }
    }
    else
    {
        printf("XY_toFILE function out err!\n");
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////functions for crystal and POSCAR
int  MSD(struct POSCAR ** poscar,struct DATA **output,int All_fram,int init_fram,int cal_time,int skip_time,double time_step )
{
	
	int Time_sample=100; //this is the number of Time Average samples, which short the calculation time.

    ////////////////////////////////////////////////////////////////////////////////////
    int number_points=(cal_time-1)/skip_time + 1;

    struct DATA *out=DATA_init(2,number_points, 1 + (poscar[0]->num_type)*4);
    for(int i=0;i<All_fram;i++)
    {
        POSCAR_TranCoorType(poscar[i],1);
    }
    //////////////////////////
    int skip_T=0;
    int set_numbers=0;
	for(int step=0;step < cal_time; step+=skip_time)
	{
	    	/////////////////////////////////////////////////////
            DATA_set(out,1.0*step*time_step,set_numbers,0);// cal culate x coordinate

           skip_T=(All_fram-step-init_fram-1)/(Time_sample-1);
            if(skip_T < 1) skip_T=1;

	    	for(int num_type=0;num_type<(poscar[0]->num_type);num_type++)//每种原子的MSD
	    	{

                /////////////////////////////////////////////////////////////////////////////////
	    		int number_ion = (poscar[0]-> num_atoms[num_type]);
                int number_times=(All_fram-step-init_fram-1)/skip_T + 1;
	    		double para1=1.0/(1.0*number_ion*number_times  );
	    		double para2_sum=0;
	    		double para2_x=0;
	    		double para2_y=0;
	    		double para2_z=0;
	    		for(int j=init_fram;j<All_fram - step;j=j+skip_T)//Nstep - n循环
	    		{
	    			for(int i=0;i< number_ion;i++)
	    			{
	    				double dx1,dy1,dz1,dx2,dy2,dz2;
	    				POSCAR_getPOS(poscar[j+step],num_type,i,&dx1,&dy1,&dz1);
	    				POSCAR_getPOS(poscar[j],           num_type,i,&dx2,&dy2,&dz2);
    
	    				para2_sum+=(dx1-dx2)*(dx1-dx2)+(dy1-dy2)*(dy1-dy2)+(dz1-dz2)*(dz1-dz2);
	    				para2_x+=(dx1-dx2)*(dx1-dx2);
	    				para2_y+=(dy1-dy2)*(dy1-dy2);
	    				para2_z+=(dz1-dz2)*(dz1-dz2);
	    			}
	    		}
                DATA_set(out,para1*para2_sum ,set_numbers,num_type*4  +1);
                DATA_set(out,para1*para2_x        ,set_numbers,num_type*4  +2);
                DATA_set(out,para1*para2_y        ,set_numbers,num_type*4  +3);
                DATA_set(out,para1*para2_z        ,set_numbers,num_type*4  +4);

	    	}
            set_numbers++;
	}
	*output=out;
    return set_numbers;
}
int POSCAR_getAtomType(struct POSCAR * poscar,int atom_pos)
{
    for(int i=0;i<poscar->num_type;i++)
    {
        atom_pos -=  poscar->num_atoms[i];
        if(atom_pos < 0)
        {
            return i;
        }
    }
    return -1;
}
void POSCAR_getOrder(struct POSCAR * poscar,int atom_type,int *pre_order,int *end_order)
{
    int pre_num=0;
    int end_num=0;
    if(atom_type >=0)
    {
        
        for(int i=0;i<atom_type;i++)
        {
            pre_num+= poscar->num_atoms[i];
        }
        end_num=pre_num+(poscar->num_atoms[atom_type]);
    }
    else
    {
        pre_num=0;
        end_num=poscar->num_all;
    }
    *pre_order=pre_num;
    *end_order=end_num;
}
double POSCAR_getVolume(double Lattice[][3])
{
    double ax,ay,az,bx,by,bz,cx,cy,cz;
    ax=Lattice[0][0];
    ay=Lattice[0][1];
    az=Lattice[0][2];

    bx=Lattice[1][0];
    by=Lattice[1][1];
    bz=Lattice[1][2];

    cx=Lattice[2][0];
    cy=Lattice[2][1];
    cz=Lattice[2][2];

    return -az*by*cx + ay*bz*cx + az*bx*cy - ax*bz*cy - ay*bx*cz + ax*by*cz;

}
void POSCAR_getPOS_order(struct POSCAR *pos,int sel_atoms,double *x,double *y,double *z)
{
    *x=DATA2D_get(pos->atom_position,sel_atoms,0);
    *y=DATA2D_get(pos->atom_position,sel_atoms,1);
    *z=DATA2D_get(pos->atom_position,sel_atoms,2);
}
void POSCAR_setPOS_order(struct POSCAR *pos,int sel_atoms,double x,double y,double z)
{
    DATA_set(pos->atom_position, x,sel_atoms,0);
    DATA_set(pos->atom_position, y,sel_atoms,1);
    DATA_set(pos->atom_position, z,sel_atoms,2);
}

void POSCAR_free(struct POSCAR *poscar)
{
    if(poscar != NULL)
    {
        if(poscar->atom_position!=NULL)
        {
            DATA_free(poscar->atom_position);
        }
        if(poscar->external != NULL)
        {
            DATA_free(poscar->external);
        }
        if(poscar->external_2 != NULL)
        {
            DATA_free(poscar->external_2);
        }
        if(poscar->atom_position_ext !=NULL)
        {
            DATA_free(poscar->atom_position_ext);
        }
        free(poscar);
        poscar=NULL;
    }

}
void POSCARS_free(struct POSCAR **poscar,int numbers)
{
    for(int i=0;i<numbers;i++)
    {
        POSCAR_free(poscar[i]);
    }
    free(poscar);
}
struct POSCAR * XDATCAR_getAvgPOSCAR(struct POSCAR ** poscars,int numbers,int start_fram,int end_fram,int skip_fram)
{
    struct POSCAR * posn=POSCAR_cpy(poscars[0]);
    posn->coordinate_type=1;
    double ax=0,ay=0,az=0;
    double bx=0,by=0,bz=0;
    double cx=0,cy=0,cz=0;
    int  sf,ef;
    if(start_fram==end_fram&&start_fram==0)
    {
        sf=0;
        ef=numbers;
    }
    else
    {
        sf=start_fram;
        ef=end_fram;
    }
    int num_frams=((ef-sf)-1)/skip_fram + 1;

    for(int i=sf;i<ef;i=i+skip_fram)
    {
        
         perCrystal_full(poscars[i],0,0,0);
         POSCAR_TranCoorType(poscars[i],1);
         ax+=poscars[i]->Lattice[0][0];
         ay+=poscars[i]->Lattice[0][1];
         az+=poscars[i]->Lattice[0][2];
         bx+=poscars[i]->Lattice[1][0];
         by+=poscars[i]->Lattice[1][1];
         bz+=poscars[i]->Lattice[1][2];
         cx+=poscars[i]->Lattice[2][0];
         cy+=poscars[i]->Lattice[2][1];
         cz+=poscars[i]->Lattice[2][2];
    }
    ax=ax/num_frams;
    ay=ay/num_frams;
    az=az/num_frams;
    bx=bx/num_frams;
    by=by/num_frams;
    bz=bz/num_frams;
    cx=cx/num_frams;
    cy=cy/num_frams;
    cz=cz/num_frams;
    posn->Lattice[0][0]=ax;
    posn->Lattice[0][1]=ay;
    posn->Lattice[0][2]=az;
    posn->Lattice[1][0]=bx;
    posn->Lattice[1][1]=by;
    posn->Lattice[1][2]=bz;
    posn->Lattice[2][0]=cx;
    posn->Lattice[2][1]=cy;
    posn->Lattice[2][2]=cz;
    
    for(int at=0;at<poscars[0]->num_all;at++)
    {
            double sx=0,sy=0,sz=0;
            for(int i=sf;i<ef;i=i+skip_fram)
            {
                double x,y,z;
               
                x=DATA2D_get(poscars[i]->atom_position,at,0);
                y=DATA2D_get(poscars[i]->atom_position,at,1);
                z=DATA2D_get(poscars[i]->atom_position,at,2);
                sx+=x;
                sy+=y;
                sz+=z;
            }    
            sx=sx/num_frams;
            sy=sy/num_frams;
            sz=sz/num_frams;
            POSCAR_setPOS_order(posn,at,sx,sy,sz);
    }
    return posn;
}
void POSCAR_TranCoorType(struct POSCAR * pos,int type)//type == 0 crystal Else car
{
	if(type == 0)
	{
		if(pos->coordinate_type == 1)
		{
            if(pos->atom_position_ext == NULL)
            {
                pos->atom_position_ext=DATA_init(2,pos->num_all,3); 
                double x1,y1,z1,x2,y2,z2;
	    	    for(int i=0;i<(pos->num_all);i++)
	    	    {
	    	    	POSCAR_getPOS_order(pos,i,&x1,&y1,&z1);
	    	    	coorTurnC(pos->Lattice,1,&x2,&y2,&z2,x1,y1,z1);
	    	    	DATA_setXYT(pos->atom_position_ext,i,x2,y2,z2);
	    	    }
                struct DATA * tmp=pos->atom_position;
                pos->atom_position=pos->atom_position_ext;
                pos->atom_position_ext=tmp;
                pos->coordinate_type=0;
            }
            else
            {
                struct DATA * tmp=pos->atom_position;
                pos->atom_position=pos->atom_position_ext;
                pos->atom_position_ext=tmp;
                pos->coordinate_type=0;
            }

		}
	}
	else if(type == 1)
	{
		if(pos->coordinate_type == 0)
		{
            if(pos->atom_position_ext == NULL)
            {
                if(pos->atom_position_ext == NULL)
                {
                    pos->atom_position_ext=DATA_init(2,pos->num_all,3); 
                    double x1,y1,z1,x2,y2,z2;
	    	        for(int i=0;i<(pos->num_all);i++)
	    	        {
	    	        	POSCAR_getPOS_order(pos,i,&x1,&y1,&z1);
	    	        	coorTurnC(pos->Lattice,0,&x2,&y2,&z2,x1,y1,z1);
	    	        	DATA_setXYT(pos->atom_position_ext,i,x2,y2,z2);
	    	        }
                }
                struct DATA * tmp=pos->atom_position;
                pos->atom_position=pos->atom_position_ext;
                pos->atom_position_ext=tmp;
                pos->coordinate_type=1;
            }
            else
            {
                struct DATA * tmp=pos->atom_position;
                pos->atom_position=pos->atom_position_ext;
                pos->atom_position_ext=tmp;
                pos->coordinate_type=1;
            }
		}
	}

}

struct POSCAR * POSCAR_cpy(struct POSCAR * org)
{
    struct POSCAR * pos_new=(struct POSCAR *)malloc(sizeof(struct POSCAR));
    /////////////////////
    strcpy(pos_new->name,org->name);

    pos_new->scale=org->scale;


    memcpy(pos_new->Lattice,org->Lattice,3*3*sizeof(double));

    memcpy((pos_new->atoms),org->atoms,sizeof(int)*50);
    memcpy(pos_new->num_atoms,org->num_atoms,sizeof(int)*50);
    pos_new->num_type=org->num_type;
    pos_new->num_all=org->num_all;

    pos_new->coordinate_type=org->coordinate_type;
    if(org->atom_position != NULL)
    {
        pos_new->atom_position = DATA_init(2,org->atom_position->dimen_data[0],org->atom_position->dimen_data[1]);
        memcpy(pos_new->atom_position->data_data,org->atom_position->data_data,sizeof(double)*(org->atom_position->all_data));
    }

    pos_new->atom_fix=org->atom_fix;
    pos_new->atom_velocity=org->atom_velocity;
    pos_new->atom_external=org->atom_external;

    if(org->external != NULL)
    {
        pos_new->external=DATA_cpy(org->external);    
    }
    else
    {
        pos_new->external=NULL;
    }

    if(org->external_2 != NULL)
    {
        pos_new->external_2 = DATA_cpy(org->external_2);
    }
    else
    {
        pos_new->external_2=NULL;
    }

    if(org->atom_position_ext != NULL)
    {
        pos_new->atom_position_ext=DATA_cpy(org->atom_position_ext);
    }
    else
    {
        pos_new->atom_position_ext=NULL;
    }

    return pos_new;

}

void coorTurnC(double Lattice[][3],int type,double *d1,double *d2,double *d3,double x,double y,double z)//晶格，笛卡尔坐标转换
{
        double b1x=Lattice[0][0];
        double b1y=Lattice[0][1];
        double b1z=Lattice[0][2];

        double b2x=Lattice[1][0];
        double b2y=Lattice[1][1];
        double b2z=Lattice[1][2];

        double b3x=Lattice[2][0];
        double b3y=Lattice[2][1];
        double b3z=Lattice[2][2];

    if(type==0)//crystal to cartesion
    {
        *d1=x * b1x + y * b2x + z * b3x ;
        *d2=x * b1y + y * b2y + z * b3y ;
        *d3=x * b1z + y * b2z + z * b3z ;

    }
    else if(type==1)//cartesion to crystal
    {
        double LM=b1z * b2y * b3x - b1y * b2z * b3x - b1z * b2x * b3y + b1x * b2z * b3y + b1y * b2x * b3z - b1x * b2y * b3z;

        double cx=- y * b2z * b3x + x * b2z * b3y + y * b2x * b3z - x * b2y * b3z + b2y * b3x * z - b2x * b3y * z;
        double cy=  y * b1z * b3x - x * b1z * b3y - y * b1x * b3z + x * b1y * b3z - b1y * b3x * z + b1x * b3y * z;
        double cz=- y * b1z * b2x + x * b1z * b2y + y * b1x * b2z - x * b1y * b2z + b1y * b2x * z - b1x * b2y * z;

        *d1=cx / LM;
        *d2=cy / LM;
        *d3=cz / LM;
    }
}

void POSCAR_getPOS(struct POSCAR* pos,int atom,int number,double *x,double *y,double *z)
{
	//printf("-------------------POSCAR_getPOS function:\n");
	if(pos!=NULL)
	{
		if(atom >= (pos->num_type))
		{
			printf("POSCAR_getPOS: No such type of atoms");
		}
		else if(number >= (pos->num_atoms[atom]))
		{
			printf("POSCAR_getPOS: No such atom");
		}
		else
		{
			int numbers=0;
			for(int i=0;i<atom;i++)
			{
				numbers+=(pos->num_atoms[i]);
			}
			numbers+=number;
			*x=DATA_get(pos->atom_position,numbers,0);
			*y=DATA_get(pos->atom_position,numbers,1);
			*z=DATA_get(pos->atom_position,numbers,2);
		}
	}
	else
	{
		printf("POSCAR_getPOS: input data POSCAR is NULL!");
	}

	//printf("end-------------------\n \n");
}
void POSCAR_setPOS(struct POSCAR* pos,int atom,int number,double x,double y,double z)
{
	//printf("-------------------POSCAR_getPOS function:\n");
	if(pos!=NULL)
	{
		if(atom >= (pos->num_type))
		{
			printf("POSCAR_getPOS: No such type of atoms");
		}
		else if(number >= (pos->num_atoms[atom]))
		{
			printf("POSCAR_getPOS: No such atom");
		}
		else
		{
			int numbers=0;
			for(int i=0;i<atom;i++)
			{
				numbers+=(pos->num_atoms[i]);
			}
			numbers+=number;
			DATA_set(pos->atom_position,x,numbers,0);
			DATA_set(pos->atom_position,y,numbers,1);
			DATA_set(pos->atom_position,z,numbers,2);
		}
	}
	else
	{
		printf("POSCAR_getPOS: input data POSCAR is NULL!");
	}

	//printf("end-------------------\n \n");
}

void perCrystal(struct POSCAR * poscar,int A,int B,int C,double *rx,double *ry,double *rz,double x,double y,double z)
{
    double cx=0,cy=0,cz=0;

        if((poscar->coordinate_type)==0)//晶格坐标
        {
            cx=x;
            cy=y;
            cz=z;
            double resx=fmod(cx,1.0);
            double resy=fmod(cy,1.0);
            double resz=fmod(cz,1.0);
            resx+=1;
            resy+=1;
            resz+=1;
            resx=fmod(resx,1.0);
            resy=fmod(resy,1.0);
            resz=fmod(resz,1.0);

            *rx=resx + A;
            *ry=resy + B;
            *rz=resz + C;
        }
        else if((poscar->coordinate_type)==1)//笛卡尔坐标
        {
            coorTurnC(poscar->Lattice,1,&cx,&cy,&cz,x,y,z);
            double resx=fmod(cx,1.0);
            double resy=fmod(cy,1.0);
            double resz=fmod(cz,1.0);
            resx+=1;
            resy+=1;
            resz+=1;
            resx=fmod(resx,1.0);
            resy=fmod(resy,1.0);
            resz=fmod(resz,1.0);
            coorTurnC(poscar->Lattice,0,&cx,&cy,&cz,resx+A,resy+B,resz+C);
            *rx=cx;
            *ry=cy;
            *rz=cz;
        }



}
void perCrystal_full(struct POSCAR * poscar,int A,int B,int C)
{

        double cx,cy,cz,resx,resy,resz;

        if((poscar->coordinate_type)==0)//晶格坐标
        {
            for(int i=0;i<poscar->num_all;i++)
            {
                cx=DATA2D_get(poscar->atom_position, i ,0);
                cy=DATA2D_get(poscar->atom_position, i ,1);
                cz=DATA2D_get(poscar->atom_position, i ,2);
                resx=fmod(cx,1.0);
                resy=fmod(cy,1.0);
                resz=fmod(cz,1.0);
                resx+=1;
                resy+=1;
                resz+=1;
                resx=fmod(resx,1.0);
                resy=fmod(resy,1.0);
                resz=fmod(resz,1.0);
                DATA2D_set(poscar->atom_position,resx+A,i,0);
                DATA2D_set(poscar->atom_position,resy+B,i,1);
                DATA2D_set(poscar->atom_position,resz+C,i,2);

            }
        }
        else if((poscar->coordinate_type)==1)//笛卡尔坐标
        {
            double x,y,z;
            for(int i=0;i<poscar->num_all;i++)
            {
                x=DATA2D_get(poscar->atom_position, i ,0);
                y=DATA2D_get(poscar->atom_position, i ,1);
                z=DATA2D_get(poscar->atom_position, i ,2);
                coorTurnC(poscar->Lattice,1,&cx,&cy,&cz,x,y,z);

                resx=fmod(cx,1.0);
                resy=fmod(cy,1.0);
                resz=fmod(cz,1.0);
                resx+=1;
                resy+=1;
                resz+=1;
                resx=fmod(resx,1.0);
                resy=fmod(resy,1.0);
                resz=fmod(resz,1.0);
                coorTurnC(poscar->Lattice,0,&cx,&cy,&cz,resx+A,resy+B,resz+C);

                DATA2D_set(poscar->atom_position,cx,i,0);
                DATA2D_set(poscar->atom_position,cy,i,1);
                DATA2D_set(poscar->atom_position,cz,i,2);

            }
        }
  

}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////functions for strings 
int STRINGS_in_subfunction_check(char A)
{
    if(A == 32 || A == '\t')
    {
        return 1;
    }
    else if(A == '\0')
    {
        return -1;
    }
    else if( A >=48 && A<=57) //number
    {
        return 0;
    }
    else if(  (A>=97&&A<=122) || (A>=65&&A<=90) )//caps
    {
        return 0;
    }
    else if(A==')'||A=='_'||A=='-'||A=='+'||A=='='||A=='{'||A=='['||A=='}'||A==']'||A==92||A=='|' || A==';'||A==':'||A==39||A==34||A=='<'||A==','||A=='>'||A=='.'||A=='/'||A=='?' )
    {
        return 0;
    }
    else if(A == '~' || A=='`' || A=='!'|| A =='@'||A=='#'||A=='$'||A=='%'||A=='^'||A=='&'||A=='*' ||A=='('  )
    {
        return 0;
    }
    else
    {
        return 1;
    }
}
struct STRINGS * STRINGS_in(char * str,int sizea,int sizeb)
{

    struct STRINGS * tm_str=(struct STRINGS *) malloc(  sizeof( struct STRINGS)   );

    char ** res=(char **)calloc(sizea,sizeof(char *) ); //   res[sizea][sizeb]
    for(int i=0;i<sizea;i++)
    {
        res[i]=(char *)calloc(sizeb,sizeof(char));
    }

    int str_lens=(int)strlen(str);

    int state=0;
    int dstate=1;
    int number_words=0;
    int tmps_num=0;
    for(int i=0;i<str_lens-1;i++)
    {
        int st=STRINGS_in_subfunction_check(str[i]);
        if(st == 0)
        {
            state=0; // read words
            if(dstate == 1)
            {
                tmps_num=0;
                number_words++;
                dstate=0;
            }
        }
        else if(st == 1)
        {
            state=1;// read space
            if(dstate == 0 )
            {
                    dstate=1;
            }
            
        }
        else if(st == -1)
        {
            break;
        }
       ////////////////////////////////////
       if(state == 0 )
       {
            res[number_words-1][tmps_num]=str[i];
            res[number_words-1][tmps_num+1]='\0';
            tmps_num++;
       }

    }

    tm_str->strings=res;
    tm_str->number_words=number_words;
    tm_str->sizea=sizea;
    tm_str->sizeb=sizeb;
    return tm_str;
}
char * STRINGS_get(struct STRINGS * str,int n)
{
    if(n<str->number_words)
    {
        return str->strings[n];
    }
    else
    {
        return " ";
    }
    return " ";
    
}
void STRINGS_free(struct STRINGS * str)
{
    if(str->strings != NULL)
    {
        for(int i=0;i<str->sizea;i++)
        {
            if(str->strings[i] != NULL) free(str->strings[i]);
        }
        free(str->strings);
    }
    free(str);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////functions for statistics

int cTon(char *name)
{
    char element[119][5]={
    "H" ,
    "HE",
    "LI",
    "BE",
    "B" ,
    "C" ,
    "N" ,
    "O" ,
    "F" ,
    "NE",
    "NA",
    "MG",
    "AL",
    "SI",
    "P" ,
    "S" ,
    "CL",
    "AR",
    "K" ,
    "CA",
    "SC",
    "TI",
    "V" ,
    "CR",
    "MN",
    "FE",
    "CO",
    "NI",
    "CU",
    "ZN",
    "GA",
    "GE",
    "AS",
    "SE",
    "BR",
    "KR",
    "RB",
    "SR",
    "Y" ,
    "ZR",
    "NB",
    "MO",
    "TC",
    "RU",
    "RH",
    "PD",
    "AG",
    "CD",
    "IN",
    "SN",
    "SB",
    "TE",
    "I" ,
    "XE",
    "CS",
    "BA",
    "LA",
    "CE",
    "PR",
    "ND",
    "PM",
    "SM",
    "EU",
    "GD",
    "TB",
    "DY",
    "HO",
    "ER",
    "TM",
    "YB",
    "LU",
    "HF",
    "TA",
    "W" ,
    "RE",
    "OS",
    "IR",
    "PT",
    "AU",
    "HG",
    "TL",
    "PB",
    "BI",
    "PO",
    "AT",
    "RN",
    "FR",
    "RA",
    "AC",
    "TH",
    "PA",
    "U" ,
    "NP",
    "PU",
    "AM",
    "CM",
    "BK",
    "CF",
    "ES",
    "FM",
    "MD",
    "NO",
    "LR",
    "RF",
    "DB",
    "SG",
    "BH",
    "HS",
    "MT",
    "DS",
    "RG",
    "CN",
    "NH",
    "FL",
    "MC",
    "LV",
    "TS",
    "OG"
    "ALL"
    };

    char name_tmp[10];
    strcpy(name_tmp,name);
    for(int i=0;name_tmp[i]!='\0';i++)
    {
        name_tmp[i]=toupper(name_tmp[i]);
    }
    
    int ele=0;
    for(int i=0;i<119;i++)
    {
        int restt=strcmp(element[i],name_tmp);
        if(restt==0)
        {
            ele=i;
            break;
        }
    }
    if(ele == 118) ele=-1;
    else
    {
        ele++;
    }

    return ele;
}
void nToc(int n,char * re_name)
{
    char Ele[118][5]={
    "H ",
    "He",
    "Li",
    "Be",
    "B ",
    "C ",
    "N ",
    "O ",
    "F ",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P ",
    "S ",
    "Cl",
    "Ar",
    "K ",
    "Ca",
    "Sc",
    "Ti",
    "V ",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y ",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I ",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W ",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U ",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og"
    };


    if(n == -1)
    {
        strcpy(re_name,"ALL");
    }
    else if(n>0&& n<118)
    {
        
        strcpy(re_name,Ele[n-1]);
        
    }
    
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////functions for mathematics


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


struct POSCAR ** XDATCAR_in_advance(char *filename,int read_fram,int *number_poscar,int * MD_type)
{
    struct POSCAR **pos=NULL;
    FILE *fp=fopen(filename,"r");

    struct POSCAR *poscar=(struct POSCAR *)calloc(1,sizeof(struct POSCAR));
    struct POSCAR * tmp_pos=NULL;
    int all_atoms=0; 
    int mdtype=0;
    int numberPOS=0;
    int total_lines=0;
    int total_poscar=0;
    int skip=1;
    if(fp != NULL)
    {
        /////////////////////////////////////////////////////////////////get file lines
        total_lines=get_file_lines_LinuxSys(filename);
        /////////////////////////////////////////////////////////////////
            
        char tmps[200];
        char ts[20];

        double x,y,z;

        if(fgets(tmps,sizeof(tmps),fp)!=NULL)
        {
            sscanf(tmps,"%s",ts);
            strcpy(poscar->name,ts);
        }
        else
        {
            printf("Err! function XDATCAR_in err!\n");
        }

         if(fgets(tmps,sizeof(tmps),fp)!=NULL)
         {
            sscanf(tmps,"%lf",&x);
            poscar->scale=x;
         }
        else
        {
            printf("Err! function XDATCAR_in err!\n");
        }


        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[0][0]=x;
        poscar->Lattice[0][1]=y;
        poscar->Lattice[0][2]=z;
        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[1][0]=x;
        poscar->Lattice[1][1]=y;
        poscar->Lattice[1][2]=z;
        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[2][0]=x;
        poscar->Lattice[2][1]=y;
        poscar->Lattice[2][2]=z;

        struct STRINGS * str;
        if(fgets(tmps,sizeof(tmps),fp)!=NULL)// read atom type
        {
            str=STRINGS_in(tmps,10,5);
            poscar->num_type= str->number_words;
            for(int i=0;i<str->number_words;i++)
            {
                poscar->atoms[i] = cTon(STRINGS_get(str,i));
            }
            STRINGS_free(str);
        }
        else
        {
            printf("Err! function XDATCAR_in err!\n");
        }                            


         if(fgets(tmps,sizeof(tmps),fp)!=NULL)//  read atom numbers
         {
            str=STRINGS_in(tmps,10,8);
            for(int i=0;i<str->number_words;i++)
            {
                poscar->num_atoms[i] = atoi(STRINGS_get(str,i));
            }
            STRINGS_free(str);
         }
         else
        {
            printf("Err! function XDATCAR_in err!\n");
        }
            
                                                        // cal totle numbers
            for(int i=0;i<poscar->num_type;i++)
            {
                all_atoms+= poscar->num_atoms[i];
            }
            poscar->num_all= all_atoms;


            fgets(tmps,sizeof(tmps),fp); // read coordinate type
            if(tmps[0]=='D'||tmps[0]=='d')
            {
                poscar->coordinate_type=0;
            }
            else if(tmps[0]=='C'||tmps[0]=='c')
            {
                poscar->coordinate_type=1;
            }
            poscar->atom_position=DATA_init(2,poscar->num_all,3);
            tmp_pos=poscar;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////md type and
            
            if(total_lines%(all_atoms+8) == 0 && (total_lines-7)%(all_atoms+1)!=0)
            {
                total_poscar=total_lines/(all_atoms+8);
                mdtype = 2;
            }
            else if(total_lines%(all_atoms+8) != 0 && (total_lines-7)%(all_atoms+1)==0)
            {
                total_poscar=(total_lines-7)/(all_atoms+1);
                mdtype= 1;
            }
            else
            {
                fseek(fp,0,SEEK_SET);
                int n=6;
                while(n>0)
                {
                    fgets(tmps,sizeof(tmps),fp);
                    n--;
                }
                n=2+all_atoms+6;
                char tmp3[200];
                while(n>0)
                {
                    fgets(tmp3,sizeof(tmp3),fp);
                    n--;                    
                }
                if(strcmp(tmps,tmp3))
                {
                    total_poscar=(total_lines-7)/(all_atoms+1);
                    mdtype= 1;
                }
                else
                {
                    total_poscar=total_lines/(all_atoms+8);
                    mdtype = 2;                    
                }

            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////  calloc memery read
            fseek(fp,0,SEEK_SET);
            if(read_fram >= 0)
            {
                if(read_fram == 0)
                {
                    skip=1;
                }
                else
                {
                    skip=total_poscar/read_fram;
                    if(skip<1) skip=1;
                }
                int read_num=( total_poscar -1 )/skip +1;
                pos=(struct POSCAR**)calloc(read_num,sizeof(struct POSCAR *));
                for(int i=0;i<read_num;i++)
                {
                    pos[i]=(struct POSCAR *)calloc(1,sizeof(struct POSCAR));
                }

                numberPOS=0;
                for(int i=0;i<total_poscar;i+=skip)
                {
                    //////////////////////////////////
                    poscar=pos[numberPOS];
                    poscar->external=NULL;
                    poscar->external_2=NULL;
                    poscar->atom_position_ext=NULL;
                    if(i ==0 || mdtype== 2)
                    {
                        fgets(tmps,sizeof(tmps),fp);
                        strcpy(poscar->name,tmp_pos->name);

                        fgets(tmps,sizeof(tmps),fp);
                        sscanf(tmps,"%lf",&x);
                        poscar->scale=x;

                        fgets(tmps,sizeof(tmps),fp);
                        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
                        poscar->Lattice[0][0]=x;
                        poscar->Lattice[0][1]=y;
                        poscar->Lattice[0][2]=z;
                        fgets(tmps,sizeof(tmps),fp);
                        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
                        poscar->Lattice[1][0]=x;
                        poscar->Lattice[1][1]=y;
                        poscar->Lattice[1][2]=z;
                        fgets(tmps,sizeof(tmps),fp);
                        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
                        poscar->Lattice[2][0]=x;
                        poscar->Lattice[2][1]=y;
                        poscar->Lattice[2][2]=z;

                        struct STRINGS * str;
                        if(fgets(tmps,sizeof(tmps),fp)!=NULL)// read atom type
                        {
                            str=STRINGS_in(tmps,10,5);
                            poscar->num_type= str->number_words;
                            for(int i=0;i<str->number_words;i++)
                            {
                                poscar->atoms[i] = cTon(STRINGS_get(str,i));
                            }
                            STRINGS_free(str);
                        }
                        else
                        {
                            printf("Err! function XDATCAR_in err!\n");
                        }                            


                         if(fgets(tmps,sizeof(tmps),fp)!=NULL)//  read atom numbers
                         {
                            str=STRINGS_in(tmps,10,8);
                            for(int i=0;i<str->number_words;i++)
                            {
                                poscar->num_atoms[i] = atoi(STRINGS_get(str,i));
                            }
                            STRINGS_free(str);
                         }
                         else
                         {
                             printf("Err! function XDATCAR_in err!\n");
                         }

                            int all_atoms=0;                                             // cal totle numbers
                            for(int i=0;i<poscar->num_type;i++)
                            {
                                all_atoms+= poscar->num_atoms[i];
                            }
                            poscar->num_all= all_atoms;


                            fgets(tmps,sizeof(tmps),fp); // read coordinate type
                            if(tmps[0]=='D'||tmps[0]=='d')
                            {
                                poscar->coordinate_type=0;
                            }
                            else if(tmps[0]=='C'||tmps[0]=='c')
                            {
                                poscar->coordinate_type=1;
                            }
                    }
                    else
                    {
                        strcpy(poscar->name,tmp_pos->name);

                        poscar->scale=tmp_pos->scale;

                        poscar->Lattice[0][0]= tmp_pos->Lattice[0][0];
                        poscar->Lattice[0][1]= tmp_pos->Lattice[0][1];
                        poscar->Lattice[0][2]= tmp_pos->Lattice[0][2];

                        poscar->Lattice[1][0]= tmp_pos->Lattice[1][0];
                        poscar->Lattice[1][1]= tmp_pos->Lattice[1][1];
                        poscar->Lattice[1][2]= tmp_pos->Lattice[1][2];

                        poscar->Lattice[2][0]= tmp_pos->Lattice[2][0];
                        poscar->Lattice[2][1]= tmp_pos->Lattice[2][1];
                        poscar->Lattice[2][2]= tmp_pos->Lattice[2][2];

                        memcpy(poscar->atoms,tmp_pos->atoms,sizeof(tmp_pos->atoms));

                        memcpy(poscar->num_atoms,tmp_pos->num_atoms,sizeof(tmp_pos->num_atoms));

                        poscar->num_type= tmp_pos-> num_type;
                        poscar->num_all = tmp_pos-> num_all;
                        poscar->coordinate_type = tmp_pos-> coordinate_type;
                        fgets(tmps,sizeof(tmps),fp);

                    }
                    /////////////////////////////////////////////////////////////////////////
                    poscar->atom_position=DATA2D_init(poscar->num_all,3);
                    for(int j=0;j<poscar->num_all;j++)
                    {
                            if(fgets(tmps,sizeof(tmps),fp)!=NULL)
                            {
                                sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
                                DATA2D_set(poscar->atom_position,x,j,0);
                                DATA2D_set(poscar->atom_position,y,j,1);
                                DATA2D_set(poscar->atom_position,z,j,2);
                            }
                            else
                            {
                                POSCAR_free(poscar);
                                printf("Err! function XDATCAR_read err!    some poscar read err!\n");
                                goto ENDS;
                            }

                    }
                    /////////////////////////////////////////////////////////////////////////skip
                    int sk_lines=0;
                    if(mdtype==1) // nvt
                    {

                        sk_lines=(skip-1)*(all_atoms + 1);
                    }
                    else if(mdtype == 2) //npt
                    {
                        sk_lines=(skip-1)*(all_atoms + 8);
                    }
                    while(sk_lines-- )
                    {
                        if(fgets(tmps,sizeof(tmps),fp)==NULL)
                        {
                            goto ENDS;                                
                        }
                    }
                    /////////////////////////////////////////////////////////////////////////
                    numberPOS++;
                }
            }
            else
            {
                numberPOS=1;
                pos=(struct POSCAR**)calloc(1,sizeof(struct POSCAR *));
                pos[0]=(struct POSCAR *)calloc(1,sizeof(struct POSCAR));
                poscar=pos[0];
                int sk_lines=0;
                if(mdtype==1) // nvt
                {
                    sk_lines=7+(abs(read_fram)-1)*(all_atoms + 1);
                }
                else if(mdtype == 2) //npt
                {
                    sk_lines=(abs(read_fram)-1)*(all_atoms + 7);
                }
                while(sk_lines-- )
                {
                    if(fgets(tmps,sizeof(tmps),fp)==NULL)
                    {
                        printf("Err! function XDATCAR_in_acvance err! input read_fram does not exist!\n");
                        goto ENDS;                                
                    }
                }
                if(mdtype == 1)//NVT
                {
                        strcpy(poscar->name,tmp_pos->name);

                        poscar->scale=tmp_pos->scale;

                        poscar->Lattice[0][0]= tmp_pos->Lattice[0][0];
                        poscar->Lattice[0][1]= tmp_pos->Lattice[0][1];
                        poscar->Lattice[0][2]= tmp_pos->Lattice[0][2];

                        poscar->Lattice[1][0]= tmp_pos->Lattice[1][0];
                        poscar->Lattice[1][1]= tmp_pos->Lattice[1][1];
                        poscar->Lattice[1][2]= tmp_pos->Lattice[1][2];

                        poscar->Lattice[2][0]= tmp_pos->Lattice[2][0];
                        poscar->Lattice[2][1]= tmp_pos->Lattice[2][1];
                        poscar->Lattice[2][2]= tmp_pos->Lattice[2][2];

                        memcpy(poscar->atoms,tmp_pos->atoms,sizeof(tmp_pos->atoms));

                        memcpy(poscar->num_atoms,tmp_pos->num_atoms,sizeof(tmp_pos->num_atoms));

                        poscar->num_type= tmp_pos-> num_type;
                        poscar->num_all = tmp_pos-> num_all;
                        poscar->coordinate_type = tmp_pos-> coordinate_type;
                        fgets(tmps,sizeof(tmps),fp);
                }
                else //NPT
                {
                        
                        fgets(tmps,sizeof(tmps),fp);
                        strcpy(poscar->name,tmps);
                        if(fgets(tmps,sizeof(tmps),fp)!=NULL)
                        {
                           sscanf(tmps,"%lf",&x);
                           poscar->scale=x;
                        }
                        else
                        {
                            printf("Err! function XDATCAR_ib err!\n");
                        }
                        fgets(tmps,sizeof(tmps),fp);
                        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
                        poscar->Lattice[0][0]=x;
                        poscar->Lattice[0][1]=y;
                        poscar->Lattice[0][2]=z;
                        fgets(tmps,sizeof(tmps),fp);
                        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
                        poscar->Lattice[1][0]=x;
                        poscar->Lattice[1][1]=y;
                        poscar->Lattice[1][2]=z;
                        fgets(tmps,sizeof(tmps),fp);
                        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
                        poscar->Lattice[2][0]=x;
                        poscar->Lattice[2][1]=y;
                        poscar->Lattice[2][2]=z;

                        struct STRINGS * str;
                        if(fgets(tmps,sizeof(tmps),fp)!=NULL)// read atom type
                        {
                            str=STRINGS_in(tmps,10,5);
                            poscar->num_type= str->number_words;
                            for(int i=0;i<str->number_words;i++)
                            {
                                poscar->atoms[i] = cTon(STRINGS_get(str,i));
                            }
                            STRINGS_free(str);
                        }
                        else
                        {
                            printf("Err! function XDATCAR_in err!\n");
                        }                            


                         if(fgets(tmps,sizeof(tmps),fp)!=NULL)//  read atom numbers
                         {
                            str=STRINGS_in(tmps,10,8);
                            for(int i=0;i<str->number_words;i++)
                            {
                                poscar->num_atoms[i] = atoi(STRINGS_get(str,i));
                            }
                            STRINGS_free(str);
                         }
                         else
                        {
                            printf("Err! function XDATCAR_in err!\n");
                        }

                        int all_atoms=0;                                             // cal totle numbers
                        for(int i=0;i<poscar->num_type;i++)
                        {
                            all_atoms+= poscar->num_atoms[i];
                        }
                        poscar->num_all= all_atoms;

                        fgets(tmps,sizeof(tmps),fp); // read coordinate type
                        if(tmps[0]=='D'||tmps[0]=='d')
                        {
                            poscar->coordinate_type=0;
                        }
                        else if(tmps[0]=='C'||tmps[0]=='c')
                        {
                            poscar->coordinate_type=1;
                        }
                }

                poscar->atom_position=DATA_init(2,poscar->num_all,3);
                for(int i=0;i<poscar->num_all;i++)
                {
                        if(fgets(tmps,sizeof(tmps),fp)!=NULL)
                        {
                            sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
                            DATA2D_set(poscar->atom_position,x,i,0);
                            DATA2D_set(poscar->atom_position,y,i,1);
                            DATA2D_set(poscar->atom_position,z,i,2);
                        }
                        else
                        {
                            POSCAR_free(poscar);
                            printf("Err! function XDATCAR_read err!    some poscar read err!\n");
                            goto ENDS;
                        }
                }                
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////  
    }
    else
    {
        printf("XDATCAR_read function: err file open err!");
    }
    POSCAR_free(tmp_pos);
    ENDS:
    *number_poscar=numberPOS;
    *MD_type=mdtype;
    return pos;
}


struct POSCAR *  POSCAR_in(char * file_name)
{
    double x,y,z;
    //printf("-------------------POSCAR_in function:\n");
    FILE *fp=fopen(file_name,"r");
    struct POSCAR *poscar=NULL;
    if(fp==NULL)
    {
        printf("POSCAR_in:err read err!");
    }
    else
    {
        char tmps[200];
        char ts[20];
        poscar=(struct POSCAR *)malloc(sizeof(struct POSCAR));
    //////////////////////////////////////////////////
        fgets(tmps,sizeof(tmps),fp); 
        sscanf(tmps,"%s",ts);
        strcpy(poscar->name,ts);
    ////////////////////////////////////////////////
        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf",&x);
        poscar->scale=x;
    ////////////////////////////////////////////////

        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[0][0]=x;
        poscar->Lattice[0][1]=y;
        poscar->Lattice[0][2]=z;
        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[1][0]=x;
        poscar->Lattice[1][1]=y;
        poscar->Lattice[1][2]=z;
        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[2][0]=x;
        poscar->Lattice[2][1]=y;
        poscar->Lattice[2][2]=z;
    ////////////////////////////////////////////////

                fgets(tmps,sizeof(tmps),fp);                             // read atom type
                struct STRINGS * str=STRINGS_in(tmps,10,5);
                poscar->num_type= str->number_words;
                for(int i=0;i<str->number_words;i++)
                {
                    poscar->atoms[i] = cTon(STRINGS_get(str,i));
                }
                STRINGS_free(str);

                fgets(tmps,sizeof(tmps),fp); //  read atom numbers
                str=STRINGS_in(tmps,10,8);
                for(int i=0;i<str->number_words;i++)
                {
                    poscar->num_atoms[i] = atoi(STRINGS_get(str,i));
                }
                STRINGS_free(str);

                fgets(tmps,sizeof(tmps),fp); // read coordinate type
                str=STRINGS_in(tmps,10,30);
                char * tl=STRINGS_get(str,0);
                if(tl[0] == 'D' || tl[0]=='d') poscar->coordinate_type = 0;
                else if( tl[0] == 'C' || tl[0]=='c')  poscar->coordinate_type = 1;
                else
                {
                    printf("Err! function XDATCAR_in err ! for NPT ensemble the coordinate type cant be found!\n");
                    STRINGS_free(str);
                    return NULL;
                }
                STRINGS_free(str);

                int all_atoms=0;                                             // cal totle numbers
                for(int i=0;i<poscar->num_type;i++)
                {
                    all_atoms+= poscar->num_atoms[i];
                }
                poscar->num_all= all_atoms;

                poscar->atom_external=0; // set the datatype
                poscar->atom_fix=0;
                poscar->atom_velocity=0;
                poscar->external=NULL;
                poscar->external_2=NULL;
                poscar->atom_position_ext=NULL;
                
    //////////////////////////////////////////////////
        poscar->atom_position=DATA_init(2,poscar->num_all,3);
        for(int i=0;i<poscar->num_all;i++)
        {
            fgets(tmps,sizeof(tmps),fp);
            sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
            DATA_set(poscar->atom_position,x,i,0);
            DATA_set(poscar->atom_position,y,i,1);
            DATA_set(poscar->atom_position,z,i,2);
        }

        fclose(fp);
        poscar->external=NULL;

        
    }
    return poscar;
    //printf("end-------------------\n \n");
}
void POSCAR_toFILE(struct POSCAR *poscar,char *filename)
{
    char names[10];
	//printf("-------------------POSCAR_in function:\n");
	if(poscar!=NULL)
	{
		FILE *fp=fopen(filename,"w");
		fprintf(fp,"%s\n",poscar->name);
		fprintf(fp,"  %lf\n",poscar->scale);
		fprintf(fp,"%.8f\t%.8f\t%.8f\n",poscar->Lattice[0][0],poscar->Lattice[0][1],poscar->Lattice[0][2]);
		fprintf(fp,"%.8f\t%.8f\t%.8f\n",poscar->Lattice[1][0],poscar->Lattice[1][1],poscar->Lattice[1][2]);
		fprintf(fp,"%.8f\t%.8f\t%.8f\n",poscar->Lattice[2][0],poscar->Lattice[2][1],poscar->Lattice[2][2]);
		for(int i=0;i< (poscar->num_type);i++)
		{
            nToc(poscar->atoms[i],names);
			fprintf(fp,"%s ",names);
		}
		fprintf(fp,"\n");
		for(int i=0;i< (poscar->num_type);i++)
		{
			fprintf(fp,"%d ",poscar->num_atoms[i]);
		}
		fprintf(fp,"\n");
		if(poscar->coordinate_type == 0) fprintf(fp,"Dir\n");
		else fprintf(fp,"car\n");
		for(int i=0;i<poscar->num_all;i++)
		{
			fprintf(fp," %.8f\t%.8f\t%.8f\n",DATA_get(poscar->atom_position,i,0),DATA_get(poscar->atom_position,i,1),DATA_get(poscar->atom_position,i,2));
		}
		fclose(fp);
	}
	else
	{
		printf("POSCAR_toFILE err! POSCAR is NULL!\n");
	}
	 //printf("end-------------------\n \n");
}
struct POSCAR *  FPOSCAR_in(char * file_name)
{
    double x,y,z;
    //printf("-------------------POSCAR_in function:\n");
    FILE *fp=fopen(file_name,"r");
    struct POSCAR *poscar=NULL;
    if(fp==NULL)
    {
        printf("POSCAR_in:err read err!");
    }
    else
    {
        char tmps[200];
        char ts[20];
        poscar=(struct POSCAR *)malloc(sizeof(struct POSCAR));

        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%s",ts);
        strcpy(poscar->name,ts);

        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf",&x);
        poscar->scale=x;

        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[0][0]=x;
        poscar->Lattice[0][1]=y;
        poscar->Lattice[0][2]=z;
        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[1][0]=x;
        poscar->Lattice[1][1]=y;
        poscar->Lattice[1][2]=z;
        fgets(tmps,sizeof(tmps),fp);
        sscanf(tmps,"%lf %lf %lf",&x,&y,&z);
        poscar->Lattice[2][0]=x;
        poscar->Lattice[2][1]=y;
        poscar->Lattice[2][2]=z;

            fgets(tmps,sizeof(tmps),fp);                          // read atom type
            struct STRINGS * str=STRINGS_in(tmps,10,5);
            poscar->num_type= str->number_words;
            for(int i=0;i<str->number_words;i++)
            {
                poscar->atoms[i] = cTon(STRINGS_get(str,i)); 
            }
            STRINGS_free(str);

            fgets(tmps,sizeof(tmps),fp); //  read atom numbers
            str=STRINGS_in(tmps,10,8);
            for(int i=0;i<str->number_words;i++)
            {
                poscar->num_atoms[i] = atoi(STRINGS_get(str,i));
            }
            STRINGS_free(str);

            fgets(tmps,sizeof(tmps),fp); //  read Dy
            str=STRINGS_in(tmps,10,50);
            char *tmp_read=STRINGS_get(str,0);
            if(tmp_read[0]=='S'||tmp_read[0]=='s')
            {
                STRINGS_free(str);
            }
            else
            {
                printf("Err !!   function FPOSCAR_in err! Not dynamic POSCAR!\n");
                goto ends;
            }
            
            fgets(tmps,sizeof(tmps),fp); // read coordinate type
            str=STRINGS_in(tmps,10,30);
            char * tl=STRINGS_get(str,0);
            if(tl[0] == 'D' || tl[0]=='d') poscar->coordinate_type = 0;
            else if( tl[0] == 'C' || tl[0]=='c')  poscar->coordinate_type = 1;
            else
            {
                printf("Err! function FPOSCAR_in err !  coordinate type cant be found!\n");
                STRINGS_free(str);
                goto ends;
            }
            STRINGS_free(str);
            int all_atoms=0;                                             // cal totle numbers
            for(int i=0;i<poscar->num_type;i++)
            {
                all_atoms+= poscar->num_atoms[i];
            }
            poscar->num_all= all_atoms;

    //////////////////////////////////////////////////
        
        poscar->atom_fix = 1;
        poscar->atom_velocity =0;
        poscar->atom_external =0;


        char A,B,C;
        poscar->atom_position=DATA_init(2,poscar->num_all,3);
        poscar->external=DATA_init(2,poscar->num_all,3);
        for(int i=0;i<poscar->num_all;i++)
        {
            fgets(tmps,sizeof(tmps),fp);
            sscanf(tmps,"%lf %lf %lf %c %c %c",&x,&y,&z,&A,&B,&C);
            DATA_set(poscar->atom_position,x,i,0);
            DATA_set(poscar->atom_position,y,i,1);
            DATA_set(poscar->atom_position,z,i,2);

            if(A == 'T' || A == 't')             DATA_set(poscar->external,1,i,0);
            else if (A == 'F' || A == 'f')   DATA_set(poscar->external,0,i,0);

            if(B == 'T' || B == 't')             DATA_set(poscar->external,1,i,1);
            else if (B == 'F' || B == 'f')   DATA_set(poscar->external,0,i,1);

            if(C == 'T' || C == 't')             DATA_set(poscar->external,1,i,2);
            else if (C == 'F' || C == 'f')   DATA_set(poscar->external,0,i,2);

        }
        poscar->external_2=NULL;
        poscar->atom_position_ext=NULL;
        fclose(fp);
    }
    return poscar;
    ends:
    free(poscar);
    return NULL;
    //printf("end-------------------\n \n");
}
void FPOSCAR_toFILE(struct POSCAR *poscar,char *filename)
{
    char names[10];
	//printf("-------------------POSCAR_in function:\n");
	if(poscar!=NULL && poscar->atom_fix == 1)
	{
		FILE *fp=fopen(filename,"w");
		fprintf(fp,"%s\n",poscar->name);
		fprintf(fp,"  %lf\n",poscar->scale);
		fprintf(fp,"%lf  %lf  %lf\n",poscar->Lattice[0][0],poscar->Lattice[0][1],poscar->Lattice[0][2]);
		fprintf(fp,"%lf  %lf  %lf\n",poscar->Lattice[1][0],poscar->Lattice[1][1],poscar->Lattice[1][2]);
		fprintf(fp,"%lf  %lf  %lf\n",poscar->Lattice[2][0],poscar->Lattice[2][1],poscar->Lattice[2][2]);
		for(int i=0;i< (poscar->num_type);i++)
		{
            nToc(poscar->atoms[i],names);
			fprintf(fp," %s ",names);
		}
		fprintf(fp,"\n");
		for(int i=0;i< (poscar->num_type);i++)
		{
			fprintf(fp," %d ",poscar->num_atoms[i]);
		}
		fprintf(fp,"\n");
        fprintf(fp,"Selective dynamics\n");
		if(poscar->coordinate_type == 0) fprintf(fp,"direct\n");
		else fprintf(fp,"cartesian\n");
 
        char A,B,C;
        int xm,ym,zm;
		for(int i=0;i<poscar->num_all;i++)
		{
   
            xm=DATA_get(poscar->external,i,0);
            ym=DATA_get(poscar->external,i,1);
            zm=DATA_get(poscar->external,i,2);
            if(xm == 1) A='T';
            else A='F';
            if(ym == 1) B='T';
            else B='F';
            if(zm == 1) C='T';
            else C='F';
			fprintf(fp," %.8f  %.8f  %.8f  %c  %c  %c\n",DATA_get(poscar->atom_position,i,0),DATA_get(poscar->atom_position,i,1),DATA_get(poscar->atom_position,i,2),A,B,C);
		}
		fclose(fp);
	}
	else if(poscar!=NULL && poscar->atom_fix == 0)
    {
        POSCAR_toFILE(poscar,filename);
    }
    else
	{
		printf("POSCAR_toFILE err! POSCAR is NULL!\n");
	}
	 //printf("end-------------------\n \n");
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



















