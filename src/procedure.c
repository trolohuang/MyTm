#include "ctools.h"

double func(double x,double * para,int n_para)
{
    return 1.0/( 1.0 + exp(-para[0]*(x-para[1])) );
}
void G_func(double x,double *para,int n_para,double * res)
{
    double y=func(x,para,n_para);
    res[0]=(x-para[1])*y*(1-y);
    res[1]=-para[0]*y*(1-y);
}
double functionsFit_MSE(struct DATA * input,double * para, int n_para ,double (* fit_function)(double x, double * para, int n_para), void (* gradient_fit_function)(double x, double * para, int n_para, double * cal_res), int iter_steps, double stop_prec)
{
    double * G_para=(double *)calloc(n_para,sizeof(double));
    double * LGM=(double *)calloc(n_para,sizeof(double));
    double * LGV=(double *)calloc(n_para,sizeof(double));
    double * re_G_para=(double *)calloc(n_para,sizeof(double));

    double * tmp_result=(double *)calloc(input->dimen_data[0],sizeof(double));
    double final_loss=-1;

    double init_lr=0.001;
    double end_lr=0.0000001;
    double alpha=2*log(1.0/end_lr)/(iter_steps);
    double lr=init_lr;
    int steps=1;
    while(steps>0)
    {
        memset(G_para,0,n_para*sizeof(double));
        ///////////////////////////////////////
        double loss=0;
        for(int i=0;i<input->dimen_data[0];i++)
        {
            tmp_result[i]=fit_function(DATA2D_get(input,i,0),para,n_para);
            loss+=0.5*pow(tmp_result[i] - DATA2D_get(input,i,1),2);
        }
        loss/=(input->dimen_data[0]);
        if(loss - stop_prec <0)
        {
            final_loss=loss;
            break;
        }
        else if(steps>iter_steps)
        {
            final_loss=loss;
            printf("Warning function functionsFit_MSE has not converged!\n");
            break;
        }
        else
        {
            /////////////////////////////////////////////////////////////for gradient
            for(int i=0;i<input->dimen_data[0];i++)
            {
                gradient_fit_function(DATA2D_get(input,i,0),para,n_para,re_G_para);
                for(int n=0;n<n_para;n++)
                {
                    G_para[n] += -1.0*(DATA2D_get(input,i,1) - tmp_result[i])*re_G_para[n];
                }
            }
            for(int n=0;n<n_para;n++)
            {
                G_para[n] /=2.0/(input->dimen_data[0]);
            }
            /////////////////////////////////////////////////////////////for adam update parameters
            for(int n=0;n<n_para;n++)
            {
                LGM[n] = 0.9*LGM[n] + 0.1*G_para[n];
                LGV[n] = 0.999*LGV[n] + 0.001*pow(G_para[n],2);
                double tutam=LGM[n]/(1 - pow(0.9,steps));
                double tutav=LGV[n]/(1 - pow(0.999,steps) );
                para[n]=para[n] - lr*(tutam/( 0.0000001 +  sqrt(tutav)  ));
            }
            /////////////////////////////////////////////////////////////
        }
        lr = init_lr*exp(-alpha * steps);
        steps++;
    }


    free(G_para);
    free(LGM);
    free(LGV);
    free(re_G_para);
    free(tmp_result);



    return final_loss;
}
void bubble_sort(struct DATA * input, int n)
{
    for(int i = 0; i < n - 1; i++)
    {
        for(int j = 0; j < n - 1 - i; j++)
        {
            double a1=DATA_get(input,j,0);
            double a2=DATA_get(input,j+1,0);
            double b1=DATA_get(input,j,1);
            double c1=DATA_get(input,j,2);
            double b2=DATA_get(input,j+1,1);
            double c2=DATA_get(input,j+1,2);
            if ((a1 - a2) > 0.0001)
            {
                DATA_set(input,a2,j,0);
                DATA_set(input,b2,j,1);
                DATA_set(input,c2,j,2);

                DATA_set(input,a1,j+1,0);
                DATA_set(input,b1,j+1,1);
                DATA_set(input,c1,j+1,2);                
            }
        }
    }
}
double heatingUP(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions )
{
    ///////////////////////////////////////////////////////////
    int mtypes=(int)inputs->atom_method[0];
    if(mtypes == 0) mtypes=1;
    if(mtypes == 3 || mtypes == -3)
    {
        inputs->so_check=solidfication_init(1);
    }
    else
    {
        inputs->so_check=solidfication_init(mtypes);
    }
        
    ///////////////////////////////////////////////////////////
    FILE * fp=NULL;
    int MD_steps_back=inputs->vasp_NSW;
    ////////////////////////////////////////////////////////////build mode
    fp=fopen("MyTm_result.dat","a");
    printf("2.build mode...\n"); 
    fprintf(fp,"2.build mode...\n");
    fprintf(fp,"----------------------------------------------------------------------------\n"); 
    struct POSCAR * pos_solid=NULL;
    struct  POSCAR * pos_liquid=NULL;
    struct POSCAR * pos=POSCAR_cpy(init_poscar);
    double melting_tem=0;
    //crystal type transitions
    //set supercell
    if(inputs->A >0 && inputs->B>0 && inputs->C>0)
    {
        POSCAR_setSuperCell(pos,inputs->A,inputs->B,inputs->C);
    }
    /////////////////////////////////////////////////////////////learn

    if(mtypes == 3 || mtypes == -3)
    {
        inputs->so_check->solid=pos;
        inputs->so_check->run_type=0;
        melting_tem=solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1; 

        inputs->init_tem_check=1;
    }
    else
    {
        inputs->so_check->solid=pos;
        inputs->so_check->run_type=0;
        melting_tem=solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1;         
    }
    
    /////////////////////////////////////////////////////////////

    printf(".\tRelaxing...\n");
    fprintf(fp,".\tRelaxing...\n");

    struct DATA * tem;
    struct DATA * pre;
    struct DATA * vol;
    struct POSCAR * pos_mode;
    double solidfication;
    struct POSCAR * pos_pre=NULL;
    /////////////////////////////////////////////////////////////init temperature check and get solid and liquid struct
    
    
    if(inputs->init_tem_check != 0)
    {
        int nsw_back=inputs->vasp_NSW;
 
        printf(".\tRecheck the solid temperature: \n");
        fprintf(fp,".\tRecheck the solid temperature: \n");
        fflush(stdout);
    
        int temper_count=0;
        
        int quit_state=0;
        double pre_tmp=-1000;
        int vas=-1;
        pos_solid=NULL;
        pos_pre=NULL;
        while(quit_state == 0 && inputs->init_tem_check != 0)//get solid 
        {
            printf(".\t\t%.2f ",inputs->solid_tem);
            fprintf(fp,".\t\t%.2f ",inputs->solid_tem);
            fflush(stdout);
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->solid_tem ,inputs->solid_tem,inputs,-1,1);
            POSCAR_free(pos_solid);
            pos_solid=functions->read_MDstructure();

            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            //issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_solid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            /////////////////////////////////////////////////////////////////
            // stable
            printf("solid fraction:%.2f ",solidfication);
            if(solidfication - inputs->decision < 0) // liquid
            {
                printf("liquid\n");
                fprintf(fp,"liquid\n");
                if(inputs->init_tem_check != 0 &&vas == -1)
                {
                    inputs->solid_tem -= 500;
                    if(inputs->solid_tem <= 0)
                    {
                        inputs->solid_tem = (inputs->solid_tem + 500) * 0.5;
                    }
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 &&vas == 0)
                {
                    inputs->solid_tem -= 500;
                    if(inputs->solid_tem <= 0)
                    {
                        inputs->solid_tem = (inputs->solid_tem + 500) * 0.5;
                    }
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 1)
                {
                    inputs->solid_tem=pre_tmp;
                    quit_state=-1;
                }
                else
                {
                    printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                    return -1;
                }
                
            }
            else if(solidfication - inputs->decision > 0) //solid
            {

                    if(pos_pre != NULL) POSCAR_free(pos_pre);
                    pos_pre=POSCAR_cpy(pos_solid);
                    printf("solid\n");
                    fprintf(fp,"solid\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=1; 
                        pre_tmp=inputs->solid_tem;

                        inputs->solid_tem += 500;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==0)
                    {
                        quit_state=-1;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==1)
                    {
                        vas=1; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem += 500;
    
                    }
                    else
                    {
                        quit_state=-1;
                    }
                    
            }
            else //uncertain
            {
                    printf("\n");
                    fprintf(fp,"\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=0; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem -= 500;
                        if(inputs->solid_tem <=0)
                        {
                            inputs->solid_tem=(inputs->solid_tem+500) * 0.5;
                        } 
     
                    }
                    else if(inputs->init_tem_check != 0 && vas==0)
                    {
                        vas=0; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem -= 500;
                        if(inputs->solid_tem <=0)
                        {
                            inputs->solid_tem=(inputs->solid_tem+500) * 0.5;
                        } 
          
                    }
                    else if(inputs->init_tem_check != 0 && vas==1)
                    {
                        inputs->solid_tem=pre_tmp;
                        quit_state=-1;
                    }
                    else
                    {
                        printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                        return -1;
                    }
                    

            }                

            DATA_free(tem);
            DATA_free(vol);
            DATA_free(pre);
            fflush(stdout);

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                return -1;
            }
        }
        POSCAR_free(pos_solid);
        pos_solid=pos_pre;
        pos_pre=NULL;

    
        printf("\n");
        fprintf(fp,"\n");
    
        temper_count=0;
        printf(".\tRecheck the liquid temperature \n");

        fprintf(fp,".\tRecheck the liquid Temperature \n");

        fflush(stdout);

        inputs->vasp_NSW =MD_steps_back;
        quit_state=0;
        pre_tmp=-1000;
        vas=-1;
        pos_liquid=NULL;
        pos_pre=NULL;
        if(inputs->melt_tem < inputs->solid_tem) inputs->melt_tem=inputs->solid_tem;
        
        while(quit_state == 0 && inputs->init_tem_check != 0)//get liquid
        {
            printf(".\t\t%.2f ",inputs->melt_tem);
            fprintf(fp,".\t\t %.2f ",inputs->melt_tem);
            fflush(stdout);
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->melt_tem ,inputs->melt_tem,inputs,-1,1);
            POSCAR_free(pos_liquid);
            pos_liquid=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            //issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_liquid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            /////////////////////////////////////////////////////////////////
            // stable
            printf("solid fraction:%.2f ",solidfication);
            if(solidfication -inputs->decision_liquid > 0) // solid
            {
                printf("solid\n");
                fprintf(fp,"solid\n");
                if(inputs->init_tem_check != 0 && vas == -1)
                {
                    inputs->melt_tem += 500;
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 0)
                {
                    inputs->melt_tem += 500;
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 1)
                {
                    inputs->melt_tem=pre_tmp;
                    quit_state=-1;
                }
                else
                {
                    printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                    return -1;
                }
                
            }
            else if(solidfication - inputs->decision_liquid < 0) //liquid
            {
                    if(pos_pre != NULL) POSCAR_free(pos_pre);
                    pos_pre=POSCAR_cpy(pos_liquid);
                    printf("liquid\n");
                    fprintf(fp,"liquid\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=1; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem -=500;
                        if(inputs->melt_tem <=0)
                        {
                            inputs->melt_tem=(inputs->melt_tem + 500) * 0.5;
                        }
                
                    }
                    else  if(inputs->init_tem_check != 0 && vas==0)
                    {

                        quit_state=-1;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==1)
                    {
                        vas=1; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem -=500;
                        if(inputs->melt_tem <=0)
                        {
                            inputs->melt_tem=(inputs->melt_tem + 500) * 0.5;
                        }
                   
                    }
                    else if(inputs->melt_tem <= 500)
                    {
                        quit_state=-1;
                    }
                    else
                    {
                        quit_state=-1;
                    }
                    
                
            }
            else //uncertain
            {
                printf("\n");
                fprintf(fp,"\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=0; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem += 500;
                    
                    }
                    else if(inputs->init_tem_check != 0 && vas==0)
                    {
                        vas=0; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem += 500;
                   
                    }
                    else if(inputs->init_tem_check != 0 && vas==1)
                    {
                        inputs->melt_tem=pre_tmp;
                        quit_state=-1;
                    }
                    else
                    {
                        printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                        return -1;
                    }
                    

            }                

            DATA_free(tem);
            DATA_free(vol);
            DATA_free(pre);
            fflush(stdout);

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                return -1;
            }
        }
    
        POSCAR_free(pos_liquid);
        pos_liquid=pos_pre;
        pos_pre=NULL;


        inputs->vasp_NSW=nsw_back;

    }
    else
    {
        int issolid;
        printf(".\tRelax solid structure: ");
        fprintf(fp,".\tRelax solid structure: ");
        fflush(stdout);
        int quit_state=0;
        int temper_count=0;
        while(quit_state == 0 ) // for solid
        {
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->solid_tem ,inputs->solid_tem,inputs,-1,1);
            pos_solid=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_solid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            DATA_free(tem);
            DATA_free(pre);
            DATA_free(vol);

            if(solidfication - inputs->decision <0)
            {
                printf("Warning! get solid and liquid struct err!\ninit liquid temperature is too high!\nget solid struct failed.\ndecrease startTem or open init_tem_check.\nif solid phase is required, the result might err!\n");
                quit_state=1;                  
            }
            else
            {
                if(issolid == 1)
                {
                    quit_state=1;
                }
                else
                {
                    POSCAR_free(pos_solid);
                    if(inputs->vasp_NSW <10000) inputs->vasp_NSW *=2.0;
                    else inputs->vasp_NSW +=MD_steps_back;
                    temper_count++;                     
                }
            }

           

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err!  adjust init temperature or increase max_step in run command.\n");
                return -1;
            }
        }



        if(inputs->isliquid == 1) // for liquid
        {

            printf(".\tRelax liquid structure: ");
            fprintf(fp,".\tRelax liquid structure: ");
            fflush(stdout);
            quit_state=0;
            temper_count=0;
            while(quit_state == 0 ) 
            {
                functions->write_runstructure(pos);
                functions->MD_relax(inputs->melt_tem ,inputs->melt_tem,inputs,-1,1);
                pos_liquid=functions->read_MDstructure();
                tem=XY_in("Temperature_MyTm");
                pre=XY_in("Pressure_MyTm");
                vol=XY_in("Volume_MyTm");
                issolid=solidfication_check_stable(tem,vol,pre,0.1);
                inputs->so_check->solid=pos_liquid;
                solidfication=solidfication_check(inputs->so_check,inputs);
                DATA_free(tem);
                DATA_free(pre);
                DATA_free(vol);

                if(issolid == 1)
                {
                    if(solidfication - inputs->decision_liquid > 0)
                    {
                        printf("Warning! get solid and liquid struct err!\ninit liquid temperature is too low!\nget liquid struct failed.\nincrease endTem or open init_tem_check.\nif liuqid phase is required, the result might err!\n");
                        quit_state=-1;
                    }
                    else
                    {
                        quit_state=-1;
                    }
                 
                }
                else
                {
                        POSCAR_free(pos_liquid);
                        if(inputs->vasp_NSW <10000) inputs->vasp_NSW *=2.0;
                        else inputs->vasp_NSW +=MD_steps_back;
                        temper_count++; 
                }

                if(temper_count >= inputs->max_step )
                {
                    printf("Err! get solid and liquid struct err!  adjust init temperature or increase max_step in run command.\n");
                    return -1;
                }
            }

        }
        else
        {
            pos_liquid=POSCAR_cpy(pos_solid);
        }      

    }
    


   
    pos_mode=mode_build(pos_solid,pos_liquid,inputs->modes);
    POSCAR_toFILE(pos_solid,"POSCAR_solid");
    POSCAR_toFILE(pos_liquid,"POSCAR_liquid");
    
    fclose(fp);



    if(mtypes == 3 || mtypes == -3)
    {
        printf(".\tML classification selected! collecting data...\n");
        solidfication_free(inputs->so_check);
        inputs->so_check=solidfication_init(mtypes);

        int run_num=inputs->vasp_NSW;
        inputs->vasp_NSW=1000;
        int numbers;
        functions->write_runstructure(pos_solid);
        functions->MD_relax(inputs->solid_tem,inputs->solid_tem,inputs,-2,-2);
        inputs->so_check->poss_solid=functions->read_MDtrajectory(&numbers);
        inputs->so_check->num_solid=numbers;

        functions->write_runstructure(pos_liquid);
        functions->MD_relax(inputs->melt_tem,inputs->melt_tem,inputs,-2,-2);
        inputs->so_check->poss_liquid=functions->read_MDtrajectory(&numbers);
        inputs->so_check->num_liquid=numbers;


        
        inputs->so_check->solid=pos;
        printf(".\t\tcollection completed. train start.\n");
        

        inputs->so_check->run_type=0;
        solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1;

        POSCARS_free(inputs->so_check->poss_liquid,numbers);
        POSCARS_free(inputs->so_check->poss_solid,numbers);
        inputs->so_check->poss_liquid=NULL;
        inputs->so_check->poss_solid=NULL;

        inputs->vasp_NSW=run_num;
    }


    POSCAR_free(pos);

    functions->write_modestructure(pos_mode);
    POSCAR_free(pos_mode);
    /////////////////////////////////////////////////////////////heating up
    int state=1;
    double stTem=inputs->solid_tem;
    double etTem=inputs->melt_tem;
    struct POSCAR ** poscars=NULL;
    struct DATA * res_tmps;

    
    fp=fopen("MyTm_result.dat","a");
    fprintf(fp,"3. Melting Temperature Search Start.\n");
    fprintf(fp,"----------------------------------------------------------------------------\n");
    printf("3. Melting Temperature Search Start.\n");  

    pos=functions->read_modestructure();

    ////////////////////////////////////
    functions->write_runstructure(pos);
    POSCAR_free(pos);

    int run_steps=(int)((etTem-stTem)/(inputs->heat_ratio));
    if(run_steps <= 0)
    {
        printf("Err! the heating up steps err! steps:%d\n",run_steps);
        return -1;
    }
    inputs->vasp_NSW=run_steps;
    printf(".\theating up start. MD Steps:%d From:%.2fK to %.2fK\n",run_steps,stTem,etTem);
    fprintf(fp,".\theating up start. MD Steps:%d From:%.2fK to %.2fK\n",run_steps,stTem,etTem);
    functions->MD_heatingUP(stTem,etTem,inputs,state,0);
        ////////////////////////////////////

        res_tmps=functions->getMDthermodynamic();

        fprintf(fp,"\nTP_%d_%d\t\tsteps:\ttemperature:\tpressure:\tvolume:\n",0,0);

        for(int i=0;i<res_tmps->dimen_data[0];i++)
        {

            fprintf(fp,"TP_%d_%d\t\t%d\t%lf\t%lf\t%lf\n",1,0,i,DATA2D_get(res_tmps,i,0),DATA2D_get(res_tmps,i,1),DATA2D_get(res_tmps,i,2));
        }
        DATA_free(res_tmps);
        ////////////////////////////////////
    if(mtypes == -1 ||mtypes == 1 ||mtypes == 0)
    {
        printf(".\t\tqmli OP classification start.\n");
        fprintf(fp,".\t\tqlmi OP classification start.\n");
    }
    else if(mtypes == -2 ||mtypes == 2 )
    {
        printf(".\t\tLD classification start.\n");
        fprintf(fp,".\t\tLD classification start.\n");
    }
    else if(mtypes == -3 ||mtypes == 3 )
    {
        printf(".\t\tML classification start.\n");
        fprintf(fp,".\t\tML classification start.\n");
    }




    int nums;

    poscars=functions->read_MDtrajectory(&nums);
    int skip=nums/100;
    if(skip<1) skip=1;
    double melt_scale=0;
    melting_tem=-1;  
  
    fprintf(fp,"\nSF_%d_%d\t\tsteps:\tsolidfication\n",1,0);
    double heat_ratios=(etTem - stTem)/nums;

  



    int sti=0;
    int eti=nums;
    while( (eti - sti)*heat_ratios - 100 > 0   )
    {
        int i=0.5*(eti + sti);
        inputs->so_check->solid=poscars[i];
        melt_scale=solidfication_check(inputs->so_check,inputs);


        printf(".\t\ttime step:%d [%d - %d]\t\tsolid fraction:%.2f\tmelting tem:%.2f\n",i,sti,eti,melt_scale,stTem + 0.5*(eti +sti)*(heat_ratios));
        fprintf(fp,"SF_%d_%d.\t\ttime step:%d [%d - %d]\t\tsolid fraction:%.2f\n",i,0,i,sti,eti,melt_scale);

        if(melt_scale - 0.5*(inputs->decision_liquid + inputs->decision) < 0 )
        {
            eti=i;
        }
        else
        {
            sti=i;
        }


    }


    POSCARS_free(poscars,nums);POSCAR_free(pos_solid);pos_solid=NULL;POSCAR_free(pos_liquid);pos_liquid=NULL;
    poscars=NULL;

    if(eti == nums)
    {

        melting_tem=stTem + 0.5*(eti +sti)*(heat_ratios);
        printf("Unable to find melting point, the value of parameter liquid_cutoff may be set too small or there may be a bug in MyTm. The result may not be accurate!!!\n");      
    }
    else
    {
        melting_tem=stTem + 0.5*(eti + sti)*(heat_ratios);
    }
    
    printf(".done\n\n");
    fprintf(fp,".done\n\n");

    printf("4. Melting point calculation: %lf (%lf+%lf)\n",melting_tem,stTem,melting_tem-stTem);
    fprintf(fp,"\n4. Melting point calculation: %lf (%lf+%lf)\n",melting_tem,stTem,melting_tem-stTem);
    fprintf(fp,"----------------------------------------------------------------------------\n");
    fprintf(fp,"\tTemperature: %lf\n",melting_tem);
    fclose(fp);

    return 0;
}

double bisection_tem(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions )    
{
    ///////////////////////////////////////////////////////////
    int mtypes=(int)inputs->atom_method[0];
    if(mtypes == 0) mtypes=1;
    if(mtypes == 3 || mtypes == -3)
    {
        inputs->so_check=solidfication_init(1);
    }
    else
    {
        inputs->so_check=solidfication_init(mtypes);
    }
    ///////////////////////////////////////////////////////////
    FILE * fp=NULL;
    int MD_steps_back=inputs->vasp_NSW;
    ////////////////////////////////////////////////////////////build mode
    fp=fopen("MyTm_result.dat","a");
    printf("2.build mode...\n"); 
    fprintf(fp,"2.build mode...\n");
    fprintf(fp,"----------------------------------------------------------------------------\n"); 
    struct POSCAR * pos_solid;
    struct  POSCAR * pos_liquid;
    struct POSCAR * pos=POSCAR_cpy(init_poscar);
    double melting_tem=0;
    //crystal type transitions

    //set supercell
    if(inputs->A >0 && inputs->B>0 && inputs->C>0)
    {
        POSCAR_setSuperCell(pos,inputs->A,inputs->B,inputs->C);
    }
    /////////////////////////////////////////////////////////////learn

    if(mtypes == 3 || mtypes == -3)
    {
        inputs->so_check->type=1;

        inputs->so_check->solid=pos;
        inputs->so_check->run_type=0;
        melting_tem=solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1; 

        inputs->init_tem_check=1;
    }
    else
    {
        inputs->so_check->solid=pos;
        inputs->so_check->run_type=0;
        melting_tem=solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1;         
    }
    
    /////////////////////////////////////////////////////////////

    printf(".\tRelaxing...\n");
    fprintf(fp,".\tRelaxing...\n");

    struct DATA * tem;
    struct DATA * pre;
    struct DATA * vol;
    struct POSCAR * pos_mode;
    double solidfication;
    struct POSCAR * pos_pre=NULL;
    /////////////////////////////////////////////////////////////init temperature check and get solid and liquid struct
    
    if(inputs->init_tem_check != 0)
    {
        int nsw_back=inputs->vasp_NSW;
 
        printf(".\tRecheck the solid temperature: \n");
        fprintf(fp,".\tRecheck the solid temperature: \n");
        fflush(stdout);
    
        int temper_count=0;
        
        int quit_state=0;
        double pre_tmp=-1000;
        int vas=-1;
        pos_solid=NULL;
        pos_pre=NULL;
        while(quit_state == 0 && inputs->init_tem_check != 0)//get solid 
        {
            printf(".\t\t%.2f ",inputs->solid_tem);
            fprintf(fp,".\t\t%.2f ",inputs->solid_tem);
            fflush(stdout);
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->solid_tem ,inputs->solid_tem,inputs,-1,1);
            POSCAR_free(pos_solid);
            pos_solid=functions->read_MDstructure();

            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            //issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_solid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            /////////////////////////////////////////////////////////////////
            // stable
            printf("solid fraction:%.2f ",solidfication);
            if(solidfication - inputs->decision < 0) // liquid
            {
                printf("liquid\n");
                fprintf(fp,"liquid\n");
                if(inputs->init_tem_check != 0 &&vas == -1)
                {
                    inputs->solid_tem -= 500;
                    if(inputs->solid_tem <= 0)
                    {
                        inputs->solid_tem = (inputs->solid_tem + 500) * 0.5;
                    }
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 &&vas == 0)
                {
                    inputs->solid_tem -= 500;
                    if(inputs->solid_tem <= 0)
                    {
                        inputs->solid_tem = (inputs->solid_tem + 500) * 0.5;
                    }
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 1)
                {
                    inputs->solid_tem=pre_tmp;
                    quit_state=-1;
                }
                else
                {
                    printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                    return -1;
                }
                
            }
            else if(solidfication - inputs->decision > 0) //solid
            {

                    if(pos_pre != NULL) POSCAR_free(pos_pre);
                    pos_pre=POSCAR_cpy(pos_solid);
                    printf("solid\n");
                    fprintf(fp,"solid\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=1; 
                        pre_tmp=inputs->solid_tem;

                        inputs->solid_tem += 500;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==0)
                    {
                        quit_state=-1;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==1)
                    {
                        vas=1; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem += 500;
    
                    }
                    else
                    {
                        quit_state=-1;
                    }
                    
            }
            else //uncertain
            {
                    printf("\n");
                    fprintf(fp,"\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=0; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem -= 500;
                        if(inputs->solid_tem <=0)
                        {
                            inputs->solid_tem=(inputs->solid_tem+500) * 0.5;
                        } 
     
                    }
                    else if(inputs->init_tem_check != 0 && vas==0)
                    {
                        vas=0; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem -= 500;
                        if(inputs->solid_tem <=0)
                        {
                            inputs->solid_tem=(inputs->solid_tem+500) * 0.5;
                        } 
          
                    }
                    else if(inputs->init_tem_check != 0 && vas==1)
                    {
                        inputs->solid_tem=pre_tmp;
                        quit_state=-1;
                    }
                    else
                    {
                        printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                        return -1;
                    }
                    

            }                

            DATA_free(tem);
            DATA_free(vol);
            DATA_free(pre);
            fflush(stdout);

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                return -1;
            }
        }
        POSCAR_free(pos_solid);
        pos_solid=pos_pre;
        pos_pre=NULL;

    
        printf("\n");
        fprintf(fp,"\n");
    
        temper_count=0;
        printf(".\tRecheck the liquid temperature \n");

        fprintf(fp,".\tRecheck the liquid Temperature \n");

        fflush(stdout);

        inputs->vasp_NSW =MD_steps_back;
        quit_state=0;
        pre_tmp=-1000;
        vas=-1;
        pos_liquid=NULL;
        pos_pre=NULL;
        if(inputs->melt_tem < inputs->solid_tem) inputs->melt_tem=inputs->solid_tem;
        
        while(quit_state == 0 && inputs->init_tem_check != 0)//get liquid
        {
            printf(".\t\t%.2f ",inputs->melt_tem);
            fprintf(fp,".\t\t %.2f ",inputs->melt_tem);
            fflush(stdout);
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->melt_tem ,inputs->melt_tem,inputs,-1,1);
            POSCAR_free(pos_liquid);
            pos_liquid=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            //issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_liquid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            /////////////////////////////////////////////////////////////////
            // stable
            printf("solid fraction:%.2f ",solidfication);
            if(solidfication -inputs->decision_liquid > 0) // solid
            {
                printf("solid\n");
                fprintf(fp,"solid\n");
                if(inputs->init_tem_check != 0 && vas == -1)
                {
                    inputs->melt_tem += 500;
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 0)
                {
                    inputs->melt_tem += 500;
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 1)
                {
                    inputs->melt_tem=pre_tmp;
                    quit_state=-1;
                }
                else
                {
                    printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                    return -1;
                }
                
            }
            else if(solidfication - inputs->decision_liquid < 0) //liquid
            {
                    if(pos_pre != NULL) POSCAR_free(pos_pre);
                    pos_pre=POSCAR_cpy(pos_liquid);
                    printf("liquid\n");
                    fprintf(fp,"liquid\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=1; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem -=500;
                        if(inputs->melt_tem <=0)
                        {
                            inputs->melt_tem=(inputs->melt_tem + 500) * 0.5;
                        }
                
                    }
                    else  if(inputs->init_tem_check != 0 && vas==0)
                    {

                        quit_state=-1;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==1)
                    {
                        vas=1; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem -=500;
                        if(inputs->melt_tem <=0)
                        {
                            inputs->melt_tem=(inputs->melt_tem + 500) * 0.5;
                        }
                   
                    }
                    else if(inputs->melt_tem <= 500)
                    {
                        quit_state=-1;
                    }
                    else
                    {
                        quit_state=-1;
                    }
                    
                
            }
            else //uncertain
            {
                printf("\n");
                fprintf(fp,"\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=0; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem += 500;
                    
                    }
                    else if(inputs->init_tem_check != 0 && vas==0)
                    {
                        vas=0; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem += 500;
                   
                    }
                    else if(inputs->init_tem_check != 0 && vas==1)
                    {
                        inputs->melt_tem=pre_tmp;
                        quit_state=-1;
                    }
                    else
                    {
                        printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                        return -1;
                    }
                    

            }                

            DATA_free(tem);
            DATA_free(vol);
            DATA_free(pre);
            fflush(stdout);

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                return -1;
            }
        }
    
        POSCAR_free(pos_liquid);
        pos_liquid=pos_pre;
        pos_pre=NULL;


        inputs->vasp_NSW=nsw_back;

    }
    else
    {
        int issolid;
        printf(".\tRelax solid structure: ");
        fprintf(fp,".\tRelax solid structure: ");
        fflush(stdout);
        int quit_state=0;
        int temper_count=0;
        while(quit_state == 0 ) // for solid
        {
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->solid_tem ,inputs->solid_tem,inputs,-1,1);
            pos_solid=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_solid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            DATA_free(tem);
            DATA_free(pre);
            DATA_free(vol);

            if(solidfication - inputs->decision <0)
            {
                printf("Warning! get solid and liquid struct err!\ninit liquid temperature is too high!\nget solid struct failed.\ndecrease startTem or open init_tem_check.\nif solid phase is required, the result might err!\n");
                quit_state=1;                  
            }
            else
            {
                if(issolid == 1)
                {
                    quit_state=1;
                }
                else
                {
                    POSCAR_free(pos_solid);
                    if(inputs->vasp_NSW <10000) inputs->vasp_NSW *=2.0;
                    else inputs->vasp_NSW +=MD_steps_back;
                    temper_count++;                     
                }
            }

           

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err!  adjust init temperature or increase max_step in run command.\n");
                return -1;
            }
        }



        if(inputs->isliquid == 1) // for liquid
        {

            printf(".\tRelax liquid structure: ");
            fprintf(fp,".\tRelax liquid structure: ");
            fflush(stdout);
            quit_state=0;
            temper_count=0;
            while(quit_state == 0 ) 
            {
                functions->write_runstructure(pos);
                functions->MD_relax(inputs->melt_tem ,inputs->melt_tem,inputs,-1,1);
                pos_liquid=functions->read_MDstructure();
                tem=XY_in("Temperature_MyTm");
                pre=XY_in("Pressure_MyTm");
                vol=XY_in("Volume_MyTm");
                issolid=solidfication_check_stable(tem,vol,pre,0.1);
                inputs->so_check->solid=pos_liquid;
                solidfication=solidfication_check(inputs->so_check,inputs);
                DATA_free(tem);
                DATA_free(pre);
                DATA_free(vol);

                if(issolid == 1)
                {
                    if(solidfication - inputs->decision_liquid > 0)
                    {
                        printf("Warning! get solid and liquid struct err!\ninit liquid temperature is too low!\nget liquid struct failed.\nincrease endTem or open init_tem_check.\nif liuqid phase is required, the result might err!\n");
                        quit_state=-1;
                    }
                    else
                    {
                        quit_state=-1;
                    }
                 
                }
                else
                {
                        POSCAR_free(pos_liquid);
                        if(inputs->vasp_NSW <10000) inputs->vasp_NSW *=2.0;
                        else inputs->vasp_NSW +=MD_steps_back;
                        temper_count++; 
                }

                if(temper_count >= inputs->max_step )
                {
                    printf("Err! get solid and liquid struct err!  adjust init temperature or increase max_step in run command.\n");
                    return -1;
                }
            }

        }
        else
        {
            pos_liquid=POSCAR_cpy(pos_solid);
        }      

    }
    



    fprintf(fp,"\n.done\n\n");
    printf("\n.done\n\n");
    POSCAR_toFILE(pos_solid,"POSCAR_solid");
    POSCAR_toFILE(pos_liquid,"POSCAR_liquid");
    pos_mode=mode_build(pos_solid,pos_liquid,inputs->modes);
    POSCAR_toFILE(pos_mode,"POSCAR_initmodel");
    fclose(fp);
    
    functions->write_modestructure(pos_mode);
    POSCAR_free(pos_mode);

    if(mtypes == 3 || mtypes == -3)
    {
        printf(".\tML classification selected! collecting data...\n");
        solidfication_free(inputs->so_check);
        inputs->so_check=solidfication_init(mtypes);

        int run_num=inputs->vasp_NSW;
        inputs->vasp_NSW=1000;
        int numbers;
        functions->write_runstructure(pos_solid);
        functions->MD_relax(inputs->solid_tem,inputs->solid_tem,inputs,-2,-2);
        inputs->so_check->poss_solid=functions->read_MDtrajectory(&numbers);
        inputs->so_check->num_solid=numbers;

        functions->write_runstructure(pos_liquid);
        functions->MD_relax(inputs->melt_tem,inputs->melt_tem,inputs,-2,-2);
        inputs->so_check->poss_liquid=functions->read_MDtrajectory(&numbers);
        inputs->so_check->num_liquid=numbers;


        
        inputs->so_check->solid=pos;
        printf(".\t\tcollection completed. train start.\n");
        

        inputs->so_check->run_type=0;
        solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1;

        POSCARS_free(inputs->so_check->poss_liquid,numbers);
        POSCARS_free(inputs->so_check->poss_solid,numbers);
        inputs->so_check->poss_liquid=NULL;
        inputs->so_check->poss_solid=NULL;

        inputs->vasp_NSW=run_num;
    }

    POSCAR_free(pos);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////bisection
   
    inputs->vasp_NSW=MD_steps_back;
    int state=1;
    double stTem=inputs->solid_tem;
    double etTem=inputs->melt_tem;
    double tmp_ratio=0.5;


    struct POSCAR * pos_tmp;
    inputs->vasp_NSW=MD_steps_back;

    fp=fopen("MyTm_result.dat","a");
    printf("3. Melting Temperature Search Start.\n");  
    fprintf(fp,"3. Melting Temperature Search Start.\n");
    fprintf(fp,"----------------------------------------------------------------------------\n");
    printf(".\tbisection method start_tem:%lf\tend_tem:%lf\n",stTem,etTem);
    fprintf(fp,".\tbisection method start temperature:%lf\tend temperature:%lf\n",stTem,etTem);
    if(mtypes == -1 ||mtypes == 1 ||mtypes == 0)
    {
        printf(".\t\tqmli OP classification start.\n");
        fprintf(fp,".\t\tqlmi OP classification start.\n");
    }
    else if(mtypes == -2 ||mtypes == 2 )
    {
        printf(".\t\tLD classification start.\n");
        fprintf(fp,".\t\tLD classification start.\n");
    }
    else if(mtypes == -3 ||mtypes == 3 )
    {
        printf(".\t\tML classification start.\n");
        fprintf(fp,".\t\tML classification start.\n");
    }
    fclose(fp);


    
    struct POSCAR * pos1;
    struct POSCAR * pos2;

    if(inputs->nsw_adaption_ratio >0) inputs->dichotomy_stable=1;

    stTem=0;
    while(state >= 0)
    {
        fp=fopen("MyTm_result.dat","a");
        /////////////////////////////////////////init MD model set
        double set_tmp=stTem+(etTem - stTem)*tmp_ratio;
        printf(".\tStep:%d\tTem:%lf [%lf-%lf,%lf]\n",state,set_tmp,stTem,etTem,tmp_ratio);
        fprintf(fp,"IF.\tStep:%d\tTemperature:%lf [from: %lf to %lf, ratio: %lf]\n",state,set_tmp,stTem,etTem,tmp_ratio);
        printf(".\t\tMode Rebuilding...\n");
        fprintf(fp,".\t\tMode Rebuilding...\n");

        if(inputs->dichotomy_structure_relax != 0)
        {
            functions->write_runstructure(pos_solid);
            functions->MD_relax(set_tmp,set_tmp,inputs,state,-1);
            pos1=functions->read_MDstructure();
            functions->write_runstructure(pos_liquid);
            functions->MD_relax(set_tmp,set_tmp,inputs,state,-2);
            pos2=functions->read_MDstructure();
            pos_mode=mode_build(pos1,pos2,inputs->modes);
            POSCAR_free(pos1);
            POSCAR_free(pos2);
        }
        else
        {
            pos_mode=mode_build(pos_solid,pos_liquid,inputs->modes);
        }
        functions->write_modestructure(pos_mode);
        functions->write_runstructure(pos_mode);
        POSCAR_free(pos_mode);
        //////////////////////////////////////////
        int sub_state=1;
        printf(".\t\tRebuild done\n");
        fprintf(fp,".\t\tRebuild done\n");
        int issolid;
        while(sub_state > 0)
        {
            functions->MD_relax(set_tmp,set_tmp,inputs,state,sub_state);
            //////////////////////////////////////////////
            pos_tmp=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_tmp;
            solidfication=solidfication_check(inputs->so_check,inputs);
            double now_tem=functions->get_MDavgTemperature();
            double now_pre=functions->get_MDavgPressure();
            if(inputs->dichotomy_stable == 1) // stable needed
            {

                if(issolid == 1 && solidfication - inputs->decision >0)
                {


                    printf(".\t\tsub_step:%d\tsolidfication:%.2f\tsolid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tsolid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    stTem=set_tmp;
                    inputs->vasp_NSW=MD_steps_back;
                    tmp_ratio=0.5;
                    POSCAR_free(pos_tmp);
                    sub_state=-1000;
                    DATA_free(tem);
                    DATA_free(pre);
                    DATA_free(vol);
                    break;
                }
                else if(solidfication - inputs->decision_liquid <0)
                {

                    printf(".\t\tsub_step:%d\tsolidfication:%.2f\tliquid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tliquid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    etTem=set_tmp;
                    inputs->vasp_NSW=MD_steps_back;
                    tmp_ratio=0.5;
                    POSCAR_free(pos_tmp);
                    sub_state=-1000;
                    
                    DATA_free(tem);
                    DATA_free(pre);
                    DATA_free(vol);
                    break;
                }
                else
                {
                    printf(".\t\tsub_step:%d\tsolidfication:%.2f\tunstable\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tunstable\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    
                    if(   (inputs->nsw_adaption_ratio > 0 ))
                    {
                        
                        inputs->vasp_NSW *=1.0+inputs->nsw_adaption_ratio;
                    }
                    else
                    {
                        functions->write_runstructure(pos_tmp);
                    }
                    POSCAR_free(pos_tmp);
                    
                }
                
            }
            else //
            {
                if( solidfication - inputs->decision >0)
                {
                    printf(".\t\tsub_step:%d\tsolidfication:%.2f\tsolid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tsolid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    stTem=set_tmp;
                    inputs->vasp_NSW=MD_steps_back;
                    tmp_ratio=0.5;
                    POSCAR_free(pos_tmp);
                    sub_state=-1000;
                    DATA_free(tem);
                    DATA_free(pre);
                    DATA_free(vol);
                    break;
                }
                else if( solidfication - inputs->decision_liquid <0)
                {


                    printf(".\t\tsub_step:%d\tsolidfication:%.2f\tliquid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tliquid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    etTem=set_tmp;
                    inputs->vasp_NSW=MD_steps_back;
                    tmp_ratio=0.5;
                    POSCAR_free(pos_tmp);
                    sub_state=-1000;
                    DATA_free(tem);
                    DATA_free(pre);
                    DATA_free(vol);
                    break;
                }
                else
                {
                    printf(".\t\tsub_step:%d\tsolidfication:%.2f\tunstable\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tunstable\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                    if(   (inputs->nsw_adaption_ratio > 0 ))
                    {
                        
                        inputs->vasp_NSW *=1.0+inputs->nsw_adaption_ratio;
                    }
                    else
                    {
                        functions->write_runstructure(pos_tmp);
                    }

                    POSCAR_free(pos_tmp);
                }


            }
            ////////////////////////////////////////////// output thermodynamics parameters
            fprintf(fp,"\nTP_%d_%d\t\tsteps:\ttemperature:\tpressure:\tvolume:\n",state,sub_state);
            if(tem->dimen_data[0] == pre->dimen_data[0] && tem->dimen_data[0] == vol->dimen_data[0])
            {
                for(int i=0;i<tem->dimen_data[0];i++)
                {

                    fprintf(fp,"TP_%d_%d\t\t%d\t%lf\t%lf\t%lf\n",state,sub_state,i,DATA2D_get(tem,i,1),DATA2D_get(pre,i,1),DATA2D_get(vol,i,1));
                }
            }
            else
            {
                printf("Err! read the Temperature_MyTm or Pressure_MyTm or Volume_MyTm file err!\n");
            }

            DATA_free(tem);
            DATA_free(pre);
            DATA_free(vol);
            //////////////////////////////////////////////
            sub_state++;
            if(sub_state>=inputs->max_sub_step)
            {
                if(inputs->decision_liquid>0 && inputs->decision-1.0<0)
                {
                    if(solidfication - 0.5<=0)
                    {
                        tmp_ratio=tmp_ratio + (1.0-tmp_ratio)*0.5;
                    }
                    else
                    {
                        tmp_ratio=tmp_ratio + (0-tmp_ratio)*0.5;
                    }
                    printf(".\t\tsub_step reach the max_sub_step! reset temperature ratio as %lf\n",tmp_ratio);
                    fprintf(fp,"IF.\t\tsub_step reach the max_sub_step! reset temperature ratio as %lf\n",tmp_ratio);
                    sub_state=-1000;
                    break;
                }
                else
                {
                    stTem=set_tmp;
                    sub_state=-1000;
                    break;
                }

                
            }
            
        }
        //////////////////////////////////////////

        if(etTem - stTem < inputs->prec_tem)
        {
            state=-1000;
        } 
        if(state > inputs->max_step)
        {
            printf("Err! reached the maximum number of max_sub_steps!\n");
            fclose(fp);
            return -1;            
        }
        //////////////////////////////////////////
        state++; 
        fclose(fp); 
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    fp=fopen("MyTm_result.dat","a");
    fprintf(fp,"IF.\tbisection method done. Temperature: %lf-%lf\n",stTem,etTem);
    melting_tem=0.5*(stTem + etTem);
    ////////////////////////////////////
    POSCAR_free(pos_solid);
    POSCAR_free(pos_liquid);
    
    printf(".done\n\n");
    fprintf(fp,".done\n\n");

    printf("4. Melting point calculation: %lf (%lf %lf)\n",melting_tem,stTem,etTem);
    fprintf(fp,"----------------------------------------------------------------------------\n");
    fprintf(fp,"IF.\tTemperature: %lf\n",melting_tem);
    fclose(fp);

    return 0;





    return 0;
}

double bisection_sta(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions )    
{
    ///////////////////////////////////////////////////////////
    int mtypes=(int)inputs->atom_method[0];
    if(mtypes == 0) mtypes=1;
    if(mtypes == 3 || mtypes == -3)
    {
        inputs->so_check=solidfication_init(1);
    }
    else
    {
        inputs->so_check=solidfication_init(mtypes);
    }
    ///////////////////////////////////////////////////////////
    FILE * fp=NULL;
    int MD_steps_back=inputs->vasp_NSW;
    ////////////////////////////////////////////////////////////build mode
    fp=fopen("MyTm_result.dat","a");
    printf("2.build mode...\n"); 
    fprintf(fp,"2.build mode...\n");
    fprintf(fp,"----------------------------------------------------------------------------\n"); 
    struct POSCAR * pos_solid;
    struct  POSCAR * pos_liquid;
    struct POSCAR * pos=POSCAR_cpy(init_poscar);
    double melting_tem=0;
    //crystal type transitions

    //set supercell
    if(inputs->A >0 && inputs->B>0 && inputs->C>0)
    {
        POSCAR_setSuperCell(pos,inputs->A,inputs->B,inputs->C);
    }
    /////////////////////////////////////////////////////////////learn

    if(mtypes == 3 || mtypes==-3)
    {
        inputs->so_check->type=1;

        inputs->so_check->solid=pos;
        inputs->so_check->run_type=0;
        melting_tem=solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1; 

        inputs->init_tem_check=1;
    }
    else
    {
        inputs->so_check->solid=pos;
        inputs->so_check->run_type=0;
        melting_tem=solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1;         
    }
 
    
    /////////////////////////////////////////////////////////////

    printf(".\tRelaxing...\n");
    fprintf(fp,".\tRelaxing...\n");

    struct DATA * tem;
    struct DATA * pre;
    struct DATA * vol;
    struct POSCAR * pos_mode;
    double solidfication;
    struct POSCAR * pos_pre=NULL;
    /////////////////////////////////////////////////////////////init temperature check and get solid and liquid struct
    

    
    if(inputs->init_tem_check != 0)
    {
        int nsw_back=inputs->vasp_NSW;
 
        printf(".\tRecheck the solid temperature: \n");
        fprintf(fp,".\tRecheck the solid temperature: \n");
        fflush(stdout);
    
        int temper_count=0;
        
        int quit_state=0;
        double pre_tmp=-1000;
        int vas=-1;
        pos_solid=NULL;
        pos_pre=NULL;
        while(quit_state == 0 && inputs->init_tem_check != 0)//get solid 
        {
            printf(".\t\t%.2f ",inputs->solid_tem);
            fprintf(fp,".\t\t%.2f ",inputs->solid_tem);
            fflush(stdout);
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->solid_tem ,inputs->solid_tem,inputs,-1,1);
            POSCAR_free(pos_solid);
            pos_solid=functions->read_MDstructure();

            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            //issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_solid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            /////////////////////////////////////////////////////////////////
            // stable
            printf("solid fraction:%.2f ",solidfication);
            if(solidfication - inputs->decision < 0) // liquid
            {
                printf("liquid\n");
                fprintf(fp,"liquid\n");
                if(inputs->init_tem_check != 0 &&vas == -1)
                {
                    inputs->solid_tem -= 500;
                    if(inputs->solid_tem <= 0)
                    {
                        inputs->solid_tem = (inputs->solid_tem + 500) * 0.5;
                    }
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 &&vas == 0)
                {
                    inputs->solid_tem -= 500;
                    if(inputs->solid_tem <= 0)
                    {
                        inputs->solid_tem = (inputs->solid_tem + 500) * 0.5;
                    }
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 1)
                {
                    inputs->solid_tem=pre_tmp;
                    quit_state=-1;
                }
                else
                {
                    printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                    return -1;
                }
                
            }
            else if(solidfication - inputs->decision > 0) //solid
            {

                    if(pos_pre != NULL) POSCAR_free(pos_pre);
                    pos_pre=POSCAR_cpy(pos_solid);
                    printf("solid\n");
                    fprintf(fp,"solid\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=1; 
                        pre_tmp=inputs->solid_tem;

                        inputs->solid_tem += 500;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==0)
                    {
                        quit_state=-1;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==1)
                    {
                        vas=1; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem += 500;
    
                    }
                    else
                    {
                        quit_state=-1;
                    }
                    
            }
            else //uncertain
            {
                    printf("\n");
                    fprintf(fp,"\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=0; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem -= 500;
                        if(inputs->solid_tem <=0)
                        {
                            inputs->solid_tem=(inputs->solid_tem+500) * 0.5;
                        } 
     
                    }
                    else if(inputs->init_tem_check != 0 && vas==0)
                    {
                        vas=0; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem -= 500;
                        if(inputs->solid_tem <=0)
                        {
                            inputs->solid_tem=(inputs->solid_tem+500) * 0.5;
                        } 
          
                    }
                    else if(inputs->init_tem_check != 0 && vas==1)
                    {
                        inputs->solid_tem=pre_tmp;
                        quit_state=-1;
                    }
                    else
                    {
                        printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                        return -1;
                    }
                    

            }                

            DATA_free(tem);
            DATA_free(vol);
            DATA_free(pre);
            fflush(stdout);

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                return -1;
            }
        }
        POSCAR_free(pos_solid);
        pos_solid=pos_pre;
        pos_pre=NULL;

    
        printf("\n");
        fprintf(fp,"\n");
    
        temper_count=0;
        printf(".\tRecheck the liquid temperature \n");

        fprintf(fp,".\tRecheck the liquid Temperature \n");

        fflush(stdout);

        inputs->vasp_NSW =MD_steps_back;
        quit_state=0;
        pre_tmp=-1000;
        vas=-1;
        pos_liquid=NULL;
        pos_pre=NULL;
        if(inputs->melt_tem < inputs->solid_tem) inputs->melt_tem=inputs->solid_tem;
        
        while(quit_state == 0 && inputs->init_tem_check != 0)//get liquid
        {
            printf(".\t\t%.2f ",inputs->melt_tem);
            fprintf(fp,".\t\t %.2f ",inputs->melt_tem);
            fflush(stdout);
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->melt_tem ,inputs->melt_tem,inputs,-1,1);
            POSCAR_free(pos_liquid);
            pos_liquid=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            //issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_liquid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            /////////////////////////////////////////////////////////////////
            // stable
            printf("solid fraction:%.2f ",solidfication);
            if(solidfication -inputs->decision_liquid > 0) // solid
            {
                printf("solid\n");
                fprintf(fp,"solid\n");
                if(inputs->init_tem_check != 0 && vas == -1)
                {
                    inputs->melt_tem += 500;
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 0)
                {
                    inputs->melt_tem += 500;
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 1)
                {
                    inputs->melt_tem=pre_tmp;
                    quit_state=-1;
                }
                else
                {
                    printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                    return -1;
                }
                
            }
            else if(solidfication - inputs->decision_liquid < 0) //liquid
            {
                    if(pos_pre != NULL) POSCAR_free(pos_pre);
                    pos_pre=POSCAR_cpy(pos_liquid);
                    printf("liquid\n");
                    fprintf(fp,"liquid\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=1; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem -=500;
                        if(inputs->melt_tem <=0)
                        {
                            inputs->melt_tem=(inputs->melt_tem + 500) * 0.5;
                        }
                
                    }
                    else  if(inputs->init_tem_check != 0 && vas==0)
                    {

                        quit_state=-1;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==1)
                    {
                        vas=1; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem -=500;
                        if(inputs->melt_tem <=0)
                        {
                            inputs->melt_tem=(inputs->melt_tem + 500) * 0.5;
                        }
                   
                    }
                    else if(inputs->melt_tem <= 500)
                    {
                        quit_state=-1;
                    }
                    else
                    {
                        quit_state=-1;
                    }
                    
                
            }
            else //uncertain
            {
                printf("\n");
                fprintf(fp,"\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=0; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem += 500;
                    
                    }
                    else if(inputs->init_tem_check != 0 && vas==0)
                    {
                        vas=0; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem += 500;
                   
                    }
                    else if(inputs->init_tem_check != 0 && vas==1)
                    {
                        inputs->melt_tem=pre_tmp;
                        quit_state=-1;
                    }
                    else
                    {
                        printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                        return -1;
                    }
                    

            }                

            DATA_free(tem);
            DATA_free(vol);
            DATA_free(pre);
            fflush(stdout);

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                return -1;
            }
        }
    
        POSCAR_free(pos_liquid);
        pos_liquid=pos_pre;
        pos_pre=NULL;


        inputs->vasp_NSW=nsw_back;

    }
    else
    {
        int issolid;
        printf(".\tRelax solid structure: ");
        fprintf(fp,".\tRelax solid structure: ");
        fflush(stdout);
        int quit_state=0;
        int temper_count=0;
        while(quit_state == 0 ) // for solid
        {
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->solid_tem ,inputs->solid_tem,inputs,-1,1);
            pos_solid=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_solid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            DATA_free(tem);
            DATA_free(pre);
            DATA_free(vol);

            if(solidfication - inputs->decision <0)
            {
                printf("Warning! get solid and liquid struct err!\ninit liquid temperature is too high!\nget solid struct failed.\ndecrease startTem or open init_tem_check.\nif solid phase is required, the result might err!\n");
                quit_state=1;                  
            }
            else
            {
                if(issolid == 1)
                {
                    quit_state=1;
                }
                else
                {
                    POSCAR_free(pos_solid);
                    if(inputs->vasp_NSW <10000) inputs->vasp_NSW *=2.0;
                    else inputs->vasp_NSW +=MD_steps_back;
                    temper_count++;                     
                }
            }

           

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err!  adjust init temperature or increase max_step in run command.\n");
                return -1;
            }
        }



        if(inputs->isliquid == 1) // for liquid
        {

            printf(".\tRelax liquid structure: ");
            fprintf(fp,".\tRelax liquid structure: ");
            fflush(stdout);
            quit_state=0;
            temper_count=0;
            while(quit_state == 0 ) 
            {
                functions->write_runstructure(pos);
                functions->MD_relax(inputs->melt_tem ,inputs->melt_tem,inputs,-1,1);
                pos_liquid=functions->read_MDstructure();
                tem=XY_in("Temperature_MyTm");
                pre=XY_in("Pressure_MyTm");
                vol=XY_in("Volume_MyTm");
                issolid=solidfication_check_stable(tem,vol,pre,0.1);
                inputs->so_check->solid=pos_liquid;
                solidfication=solidfication_check(inputs->so_check,inputs);
                DATA_free(tem);
                DATA_free(pre);
                DATA_free(vol);

                if(issolid == 1)
                {
                    if(solidfication - inputs->decision_liquid > 0)
                    {
                        printf("Warning! get solid and liquid struct err!\ninit liquid temperature is too low!\nget liquid struct failed.\nincrease endTem or open init_tem_check.\nif liuqid phase is required, the result might err!\n");
                        quit_state=-1;
                    }
                    else
                    {
                        quit_state=-1;
                    }
                 
                }
                else
                {
                        POSCAR_free(pos_liquid);
                        if(inputs->vasp_NSW <10000) inputs->vasp_NSW *=2.0;
                        else inputs->vasp_NSW +=MD_steps_back;
                        temper_count++; 
                }

                if(temper_count >= inputs->max_step )
                {
                    printf("Err! get solid and liquid struct err!  adjust init temperature or increase max_step in run command.\n");
                    return -1;
                }
            }

        }
        else
        {
            pos_liquid=POSCAR_cpy(pos_solid);
        }      

    }
    




    fprintf(fp,"\n.done\n\n");
    printf("\n.done\n\n");
    POSCAR_toFILE(pos_solid,"POSCAR_solid");
    POSCAR_toFILE(pos_liquid,"POSCAR_liquid");
    pos_mode=mode_build(pos_solid,pos_liquid,inputs->modes);
    POSCAR_toFILE(pos_mode,"POSCAR_initmodel");
  
    fclose(fp);
    
    functions->write_modestructure(pos_mode);
    POSCAR_free(pos_mode);
    ///////////////////////////////////////////////////////////////
    if(mtypes == 3 || mtypes == -3)
    {
        printf(".\tML classification selected! collecting data...\n");
        solidfication_free(inputs->so_check);
        inputs->so_check=solidfication_init(mtypes);

        int run_num=inputs->vasp_NSW;
        inputs->vasp_NSW=1000;
        int numbers;
        functions->write_runstructure(pos_solid);
        functions->MD_relax(inputs->solid_tem,inputs->solid_tem,inputs,-2,-2);
        inputs->so_check->poss_solid=functions->read_MDtrajectory(&numbers);
        inputs->so_check->num_solid=numbers;

        functions->write_runstructure(pos_liquid);
        functions->MD_relax(inputs->melt_tem,inputs->melt_tem,inputs,-2,-2);
        inputs->so_check->poss_liquid=functions->read_MDtrajectory(&numbers);
        inputs->so_check->num_liquid=numbers;


        
        inputs->so_check->solid=pos;
        printf(".\t\tcollection completed. train start.\n");
        

        inputs->so_check->run_type=0;
        solidfication_check(inputs->so_check,inputs);
        inputs->so_check->run_type=1;

        POSCARS_free(inputs->so_check->poss_liquid,numbers);
        POSCARS_free(inputs->so_check->poss_solid,numbers);
        inputs->so_check->poss_liquid=NULL;
        inputs->so_check->poss_solid=NULL;

        inputs->vasp_NSW=run_num;
    }
    POSCAR_free(pos);




    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////bisection
   
    inputs->vasp_NSW=MD_steps_back;
    int state=1;
    double stTem=inputs->solid_tem;
    double etTem=inputs->melt_tem;
    double tmp_ratio=0.5;


    struct POSCAR * pos_tmp;
    inputs->vasp_NSW=MD_steps_back;

    fp=fopen("MyTm_result.dat","a");
    printf("3. Melting Temperature Search Start.\n");  
    fprintf(fp,"3. Melting Temperature Search Start.\n");
    fprintf(fp,"----------------------------------------------------------------------------\n");
    printf(".\tbisection method start_tem:%lf\tend_tem:%lf\n",stTem,etTem);
    fprintf(fp,".\tbisection method start temperature:%lf\tend temperature:%lf\n",stTem,etTem);
    if(mtypes == -1 ||mtypes == 1 ||mtypes == 0)
    {
        printf(".\t\tqmli OP classification start.\n");
        fprintf(fp,".\t\tqlmi OP classification start.\n");
    }
    else if(mtypes == -2 ||mtypes == 2 )
    {
        printf(".\t\tLD classification start.\n");
        fprintf(fp,".\t\tLD classification start.\n");
    }
    else if(mtypes == -3 ||mtypes == 3 )
    {
        printf(".\t\tML classification start.\n");
        fprintf(fp,".\t\tML classification start.\n");
    }
    fclose(fp);

    


    struct POSCAR * pos1;
    struct POSCAR * pos2;

    struct POSCAR ** poss;
    int num_poss;

    stTem=0;
    while(state >= 0)
    {
        fp=fopen("MyTm_result.dat","a");
        /////////////////////////////////////////init MD model set
        double set_tmp=stTem+(etTem - stTem)*tmp_ratio;
        printf(".\tStep:%d\tTem:%lf [%lf-%lf,%lf]\n",state,set_tmp,stTem,etTem,tmp_ratio);
        fprintf(fp,"IF.\tStep:%d\tTemperature:%lf [from: %lf to %lf, ratio: %lf]\n",state,set_tmp,stTem,etTem,tmp_ratio);
        printf(".\t\tMode Rebuilding...\n");
        fprintf(fp,".\t\tMode Rebuilding...\n");

        if(inputs->dichotomy_structure_relax != 0)
        {
            functions->write_runstructure(pos_solid);
            functions->MD_relax(set_tmp,set_tmp,inputs,state,-1);
            pos1=functions->read_MDstructure();


            functions->write_runstructure(pos_liquid);
            functions->MD_relax(set_tmp,set_tmp,inputs,state,-2);
            pos2=functions->read_MDstructure();
            pos_mode=mode_build(pos1,pos2,inputs->modes);
            POSCAR_free(pos1);
            POSCAR_free(pos2);
        }
        else
        {
            pos_mode=mode_build(pos_solid,pos_liquid,inputs->modes);
        }
        functions->write_modestructure(pos_mode);
        functions->write_runstructure(pos_mode);
        POSCAR_free(pos_mode);
        //////////////////////////////////////////
        int sub_state=1;
        printf(".\t\tRebuild done\n");
        fprintf(fp,".\t\tRebuild done\n");
        int issolid;

        while(sub_state > 0)
        {
            functions->MD_relax(set_tmp,set_tmp,inputs,state,sub_state);
            
            //////////////////////////////////////////////
            pos_tmp=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_tmp;
            solidfication=solidfication_check(inputs->so_check,inputs);
            double now_tem=functions->get_MDavgTemperature();
            double now_pre=functions->get_MDavgPressure();


            if(solidfication - inputs->decision_liquid < 0)
            {
                printf(".\t\tsub_step:%d\tsolidfication:%.2f\tliquid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tliquid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                etTem=stTem+(etTem - stTem)*tmp_ratio;
                inputs->vasp_NSW=MD_steps_back;
                tmp_ratio=0.5;
                POSCAR_free(pos_tmp);
                sub_state=-1000;
                break;
            }
            else
            {
                if(issolid == 1 || sub_state > 2)
                {
                    if(solidfication - inputs->decision >0)
                    {
                        printf(".\t\tsub_step:%d\tsolidfication:%.2f\tsolid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                        fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tsolid\tTemperature:%.2f\tPressure:%.2f\n",sub_state,solidfication,now_tem,now_pre);
                        stTem=stTem+(etTem - stTem)*tmp_ratio;
                        inputs->vasp_NSW=MD_steps_back;
                        tmp_ratio=0.5;
                        POSCAR_free(pos_tmp);
                        sub_state=-1000;
                        break;                        
                    }
                    else
                    {
                        printf(".\t\tsub_step:%d\tsolidfication:%.2f\tTemperature:%.2f\tPressure:%.2f stable\n",sub_state,solidfication,now_tem,now_pre);
                        fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tliquid\tTemperature:%.2f\tPressure:%.2f stable\n",sub_state,solidfication,now_tem,now_pre);
                        melting_tem=now_tem;
                        inputs->vasp_NSW=MD_steps_back;
                        tmp_ratio=0.5;
                        POSCAR_free(pos_tmp);
                        sub_state=-1000;
                        state=-1000;
                        break;                        
                    }
                }
                else
                {
                    printf(".\t\tsub_step:%d\tsolidfication:%.2f\tTemperature:%.2f\tPressure:%.2f unstable\n",sub_state,solidfication,now_tem,now_pre);
                    fprintf(fp,"IF.\t\tsub_step:%d\tsolidfication:%.2f\tliquid\tTemperature:%.2f\tPressure:%.2f unstable\n",sub_state,solidfication,now_tem,now_pre);
                    
                    if( (inputs->nsw_adaption_ratio > 0))
                    {
                        inputs->vasp_NSW *=1.0+inputs->nsw_adaption_ratio;
                    }
                    else
                    {
                        functions->write_runstructure(pos_tmp);
                    }
                    
                    POSCAR_free(pos_tmp);
                }
            }
            
            
            ////////////////////////////////////////////// output thermodynamics parameters
            fprintf(fp,"\nTP_%d_%d\t\tsteps:\ttemperature:\tpressure:\tvolume:\n",state,sub_state);
            if(tem->dimen_data[0] == pre->dimen_data[0] && tem->dimen_data[0] == vol->dimen_data[0])
            {
                for(int i=0;i<tem->dimen_data[0];i++)
                {

                    fprintf(fp,"TP_%d_%d\t\t%d\t%lf\t%lf\t%lf\n",state,sub_state,i,DATA2D_get(tem,i,1),DATA2D_get(pre,i,1),DATA2D_get(vol,i,1));
                }
            }
            else
            {
                printf("Err! read the Temperature_MyTm or Pressure_MyTm or Volume_MyTm file err!\n");
            }

            DATA_free(tem);
            DATA_free(pre);
            DATA_free(vol);
            //////////////////////////////////////////////
            sub_state++;
            if(sub_state>=inputs->max_sub_step)
            {
                stTem=set_tmp;
                sub_state=-1000;
                printf("Err! reached the maximum number of max_sub_steps!\n");
                return -1;
                break;   
            }
            
        }
        //////////////////////////////////////////

        if(etTem - stTem < inputs->prec_tem)
        {
            state=-1000;
        } 
        if(state > inputs->max_step)
        {
            printf("Err! reached the maximum number of max_steps!\n");
            fclose(fp);
            return -1;            
        }
        //////////////////////////////////////////
        state++; 
        fclose(fp); 
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    poss=functions->read_MDtrajectory(&num_poss);
    int skip=1;
    if(num_poss > 200) skip=num_poss/100;
    fp=fopen("MyTm_result.dat","a");
    fprintf(fp,"\n");
    fprintf(fp,"SF.\tsteps\tsolidfraction\n");
    for(int i=0;i<num_poss;i+=skip)
    {
        inputs->so_check->solid=poss[i];
        solidfication=solidfication_check(inputs->so_check,inputs);
        fprintf(fp,"SF.\t%d\t%.4f\n",i,solidfication);
    }
    fprintf(fp,"\n");
    fclose(fp);
    POSCARS_free(poss,num_poss);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    fp=fopen("MyTm_result.dat","a");
    fprintf(fp,"IF.\tbisection method done. Temperature Range: %lf-%lf\n",stTem,etTem);
    ////////////////////////////////////
    POSCAR_free(pos_solid);
    POSCAR_free(pos_liquid);
    
    printf(".done\n\n");
    fprintf(fp,".done\n\n");

    printf("4. Melting point calculation: %lf (%lf %lf)\n",melting_tem,stTem,etTem);
    fprintf(fp,"----------------------------------------------------------------------------\n");
    fprintf(fp,"IF.\tTemperature: %lf\n",melting_tem);
    fclose(fp);

    return 0;





    return 0;
}

double smallcell(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions )    
{
    ///////////////////////////////////////////////////////////
    int mtypes=(int)inputs->atom_method[0];
    if(mtypes == 0) mtypes=1;
    if(mtypes == 3 || mtypes == -3)
    {
        inputs->so_check=solidfication_init(1);
    }
    else
    {
        inputs->so_check=solidfication_init(mtypes);
    }
    ///////////////////////////////////////////////////////////
    inputs->init_tem_check=1;
    //fprintf(fp,"----------------------------------------------------------------------------\n");
  
    ////////////////////////////////////////////////////////////build mode
    FILE * fp=fopen("MyTm_result.dat","a");
    fprintf(fp,"2.build mode...\n"); 
    fprintf(fp,"----------------------------------------------------------------------------\n");
    printf("2.build mode...\n");  
    struct POSCAR * pos_solid;
    struct  POSCAR * pos_liquid;
    struct POSCAR * pos=POSCAR_cpy(init_poscar);
    double melting_tem=0;
    //crystal type transitions
    //set supercell
    if(inputs->A >0 && inputs->B>0 && inputs->C>0)
    {
        POSCAR_setSuperCell(pos,inputs->A,inputs->B,inputs->C);
    }
    /////////////////////////////////////////////////////////////learn
    inputs->so_check->solid=pos;
    inputs->so_check->run_type=0;
    melting_tem=solidfication_check(inputs->so_check,inputs);
    inputs->so_check->run_type=1;


    double pre_st_tem=inputs->solid_tem;
    double pre_et_tem=inputs->melt_tem;
    /////////////////////////////////////////////////////////////
    printf(".\tRelaxing...\n");
    fprintf(fp,".\tRelaxing...\n");

    int MD_steps_back=inputs->vasp_NSW;
    inputs->vasp_NSW=1000;

    /////////////////////////////////////////////////////////////init temperature check and get solid and liquid struct
    struct DATA * tem=NULL;
    struct DATA * pre=NULL;
    struct DATA * vol=NULL;
    struct POSCAR * pos_mode=NULL;
    struct POSCAR * pos_pre;
    double solidfication;
    /////////////////////////////////////////////////////////////init temperature check and get solid and liquid struct
    
    if(1)
    {
        int nsw_back=inputs->vasp_NSW;
        inputs->vasp_NSW=1000;
       
        printf(".\tRecheck the solid temperature: \n");
        fprintf(fp,".\tRecheck the solid temperature: \n");
        fflush(stdout);
    
        int temper_count=0;
        
        int quit_state=0;
        double pre_tmp=-1000;
        int vas=-1;
        pos_solid=NULL;
        pos_pre=NULL;
        while(quit_state == 0 && inputs->init_tem_check != 0)//get solid 
        {
            printf(".\t\t%.2f ",inputs->solid_tem);
            fprintf(fp,".\t\t%.2f ",inputs->solid_tem);
            fflush(stdout);
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->solid_tem ,inputs->solid_tem,inputs,-1,1);
            POSCAR_free(pos_solid);
            pos_solid=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");

            inputs->so_check->solid=pos_solid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            printf("solid fraction:%.2f ",solidfication);
            /////////////////////////////////////////////////////////////////
            // stable

            if(solidfication - inputs->decision_liquid < 0) // liquid
            {
                printf("liquid\n");
                fprintf(fp,"liquid\n");
                if(inputs->init_tem_check != 0 &&vas == -1)
                {
                    inputs->solid_tem -= 500;
                    if(inputs->solid_tem <= 0)
                    {
                        inputs->solid_tem = (inputs->solid_tem + 500) * 0.5;
                    }
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 &&vas == 0)
                {
                    inputs->solid_tem -= 500;
                    if(inputs->solid_tem <= 0)
                    {
                        inputs->solid_tem = (inputs->solid_tem + 500) * 0.5;
                    }
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 1)
                {
                    inputs->solid_tem=pre_tmp;
                    quit_state=-1;
                }
                else
                {
                    printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                    return -1;
                }
                
            }
            else if(solidfication -inputs->decision > 0) //solid
            {

                    if(pos_pre != NULL) POSCAR_free(pos_pre);
                    pos_pre=POSCAR_cpy(pos_solid);
                    printf("solid\n");
                    fprintf(fp,"solid\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=1; 
                        pre_tmp=inputs->solid_tem;

                        inputs->solid_tem += 500;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==0)
                    {
                        quit_state=-1;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==1)
                    {
                        vas=1; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem += 500;
    
                    }
                    else
                    {
                        quit_state=-1;
                    }
                    
            }
            else //uncertain
            {
                    printf("\n");
                    fprintf(fp,"\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=0; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem -= 500;
                        if(inputs->solid_tem <=0)
                        {
                            inputs->solid_tem=(inputs->solid_tem+500) * 0.5;
                        } 
     
                    }
                    else if(inputs->init_tem_check != 0 && vas==0)
                    {
                        vas=0; 
                        pre_tmp=inputs->solid_tem;
                        inputs->solid_tem -= 500;
                        if(inputs->solid_tem <=0)
                        {
                            inputs->solid_tem=(inputs->solid_tem+500) * 0.5;
                        } 
          
                    }
                    else if(inputs->init_tem_check != 0 && vas==1)
                    {
                        inputs->solid_tem=pre_tmp;
                        quit_state=-1;
                    }
                    else
                    {
                        printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                        return -1;
                    }
                    

            }                

            DATA_free(tem);
            DATA_free(vol);
            DATA_free(pre);
            fflush(stdout);

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                return -1;
            }
        }
        POSCAR_free(pos_solid);
        pos_solid=pos_pre;
        pos_pre=NULL;

    
        printf("\n");
        fprintf(fp,"\n");
    
        temper_count=0;
        printf(".\tRecheck the liquid temperature \n");

        fprintf(fp,".\tRecheck the liquid Temperature \n");

        fflush(stdout);

        inputs->vasp_NSW =MD_steps_back;
        quit_state=0;
        pre_tmp=-1000;
        vas=-1;
        pos_liquid=NULL;
        pos_pre=NULL;
        if(inputs->melt_tem < inputs->solid_tem) inputs->melt_tem=inputs->solid_tem;
        
        while(quit_state == 0 && inputs->init_tem_check != 0)//get liquid
        {
            printf(" %lf_",inputs->melt_tem);
            fprintf(fp," %lf_",inputs->melt_tem);
            fflush(stdout);
            functions->write_runstructure(pos);
            functions->MD_relax(inputs->melt_tem ,inputs->melt_tem,inputs,-1,1);
            POSCAR_free(pos_liquid);
            pos_liquid=functions->read_MDstructure();
            tem=XY_in("Temperature_MyTm");
            pre=XY_in("Pressure_MyTm");
            vol=XY_in("Volume_MyTm");
            //issolid=solidfication_check_stable(tem,vol,pre,0.1);
            inputs->so_check->solid=pos_liquid;
            solidfication=solidfication_check(inputs->so_check,inputs);
            printf("solid fraction:%.2f ",solidfication);
            /////////////////////////////////////////////////////////////////
            // stable

            if(solidfication - inputs->decision > 0) // solid
            {
                printf("solid\n");
                fprintf(fp,"solid\n");
                if(inputs->init_tem_check != 0 && vas == -1)
                {
                    inputs->melt_tem += 500;
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 0)
                {
                    inputs->melt_tem += 500;
                    vas=0;
                }
                else if(inputs->init_tem_check != 0 && vas == 1)
                {
                    inputs->melt_tem=pre_tmp;
                    quit_state=-1;
                }
                else
                {
                    printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                    return -1;
                }
                
            }
            else if(solidfication -inputs->decision_liquid < 0) //liquid
            {
                    if(pos_pre != NULL) POSCAR_free(pos_pre);
                    pos_pre=POSCAR_cpy(pos_liquid);
                    printf("liquid\n");
                    fprintf(fp,"liquid\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=1; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem -=500;
                        if(inputs->melt_tem <=0)
                        {
                            inputs->melt_tem=(inputs->melt_tem + 500) * 0.5;
                        }
                
                    }
                    else  if(inputs->init_tem_check != 0 && vas==0)
                    {

                        quit_state=-1;
                    }
                    else  if(inputs->init_tem_check != 0 && vas==1)
                    {
                        vas=1; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem -=500;
                        if(inputs->melt_tem <=0)
                        {
                            inputs->melt_tem=(inputs->melt_tem + 500) * 0.5;
                        }
                   
                    }
                    else if(inputs->melt_tem <= 500)
                    {
                        quit_state=-1;
                    }
                    else
                    {
                        quit_state=-1;
                    }
                    
                
            }
            else //uncertain
            {
                    printf("\n");
                    fprintf(fp,"\n");
                    if(inputs->init_tem_check != 0 && vas==-1)
                    {
                        vas=0; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem += 500;
                    
                    }
                    else if(inputs->init_tem_check != 0 && vas==0)
                    {
                        vas=0; 
                        pre_tmp=inputs->melt_tem;
                        inputs->melt_tem += 500;
                   
                    }
                    else if(inputs->init_tem_check != 0 && vas==1)
                    {
                        inputs->melt_tem=pre_tmp;
                        quit_state=-1;
                    }
                    else
                    {
                        printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                        return -1;
                    }
                    

            }                

            DATA_free(tem);
            DATA_free(vol);
            DATA_free(pre);
            fflush(stdout);

            if(temper_count >= inputs->max_step )
            {
                printf("Err! get solid and liquid struct err! please open: init_tem_check or adjust init temperature in run command.\n");
                return -1;
            }
        }
    
        POSCAR_free(pos_liquid);
        pos_liquid=pos_pre;
        pos_pre=NULL;


        inputs->vasp_NSW=nsw_back;

    }

    fprintf(fp,"\n.done\n\n");
    printf("\n.done\n\n");
    POSCAR_toFILE(pos_solid,"POSCAR_solid");
    POSCAR_toFILE(pos_liquid,"POSCAR_liquid");
    pos_mode=mode_build(pos_solid,pos_liquid,inputs->modes);
    POSCAR_toFILE(pos_mode,"POSCAR_initmodel");
    POSCAR_free(pos);
    fclose(fp);

    inputs->vasp_NSW=MD_steps_back;
    /////////////////////////////////////////////////////////////dichotomy
    inputs->solid_tem=pre_st_tem;
    inputs->melt_tem=pre_et_tem;
    double stTem=inputs->solid_tem;
    double etTem=inputs->melt_tem;


    struct DATA *liquid_ratio=DATA2D_init(inputs->max_step,2);
    double tem_prec=(etTem-stTem)/(inputs->max_step - 1);

    struct POSCAR ** poscars=NULL;

    struct POSCAR * pos_tmp;

    fp=fopen("MyTm_result.dat","a");
    printf("3. Melting Temperature Search Start.\n");  
    fprintf(fp,"3. Melting Temperature Search Start.\n");
    fprintf(fp,"----------------------------------------------------------------------------\n");
    printf(".\tsmallcell method start_tem:%lf\tend_tem:%lf\n",stTem,etTem);
    fprintf(fp,".\tsmallcell method start temperature:%lf\tend temperature:%lf\n",stTem,etTem);
    if(mtypes == -1 ||mtypes == 1 ||mtypes == 0)
    {
        printf(".\t\tqmli OP classification start.\n");
        fprintf(fp,".\t\tqlmi OP classification start.\n");
    }
    else if(mtypes == -2 ||mtypes == 2 )
    {
        printf(".\t\tLD classification start.\n");
        fprintf(fp,".\t\tLD classification start.\n");
    }
    else if(mtypes == -3 ||mtypes == 3 )
    {
        printf(".\t\tML classification start.\n");
        fprintf(fp,".\t\tML classification start.\n");
    }
    fclose(fp);

  

    struct POSCAR * pos1;
    struct POSCAR * pos2;





    fp=fopen("MyTm_result.dat","a");


    double min_temperature=-1;
    double max_temperature=-1;
    for(int state=0;state<inputs->max_step;state++)
    {
        
        /////////////////////////////////////////MD run
        double set_tmp=stTem + state*tem_prec;
        printf(".\tStep:%d\tTem:%lf [%lf-%lf]\n",state,set_tmp,stTem,etTem);
        fprintf(fp,"IF.\tStep:%d\tTemperature:%lf [from: %lf to %lf]\n",state,set_tmp,stTem,etTem);
        printf(".\t\tMode Rebuilding...\n");
        fprintf(fp,".\t\tMode Rebuilding...\n");

        functions->write_runstructure(pos_solid);
        functions->MD_relax(set_tmp,set_tmp,inputs,state,-1);
        pos1=functions->read_MDstructure();


        functions->write_runstructure(pos_liquid);
        functions->MD_relax(set_tmp,set_tmp,inputs,state,-2);
        pos2=functions->read_MDstructure();

        pos_mode=mode_build(pos1,pos2,inputs->modes);
        functions->write_modestructure(pos_mode);

        POSCAR_free(pos1);
        POSCAR_free(pos2);
     
        printf(".\t\tRebuild done\n");
        fprintf(fp,".\t\tRebuild done\n");
        int num_solids=0;
        int num_liquids=0;

  
        double now_tem;
        double now_pre;

        pos_tmp=pos_mode;



        for(int sub_state = 0;sub_state<inputs->max_sub_step;sub_state++)
        {
            while(1)
            {
                functions->write_runstructure(pos_tmp);
                POSCAR_free(pos_tmp);
                
                functions->MD_relax(set_tmp,set_tmp,inputs,state,sub_state);
                /////////////////////////////////////////////////////////////////
                int nums;
                double tmp_sij;
                poscars=functions->read_MDtrajectory(&nums);

                pos_tmp=XDATCAR_getAvgPOSCAR(poscars,nums,0,nums-1,1);
     
                inputs->so_check->solid=pos_tmp;
                tmp_sij=solidfication_check(inputs->so_check,inputs);  
    
                POSCARS_free(poscars,nums);
   
                POSCAR_free(pos_tmp);
        
                if(tmp_sij - inputs->decision > 0  )
                {
                    now_tem=functions->get_MDavgTemperature();
                    now_pre=functions->get_MDavgPressure();

                    printf(".\t\tsub_step:%d\tTemperature:%lf\tPressure:%lf solidfication:%lf\tfin\n",sub_state,set_tmp,now_pre,tmp_sij);
                    fprintf(fp,".\t\tsub_step:%d\tTemperature:%lf\tPressure:%lf solidfication:%lf\tfin\n",sub_state,set_tmp,now_pre,tmp_sij);

                    num_solids++;

                    pos_tmp=functions->read_modestructure();
                    break;
                }
                else if(tmp_sij-inputs->decision_liquid<0)
                {
                    now_tem=functions->get_MDavgTemperature();
                    now_pre=functions->get_MDavgPressure();
                 
                    num_liquids++;

                    printf(".\t\tsub_step:%d\tTemperature:%lf\tPressure:%lf solidfication:%lf\tfin\n",sub_state,set_tmp,now_pre,tmp_sij);
                    fprintf(fp,".\t\tsub_step:%d\tTemperature:%lf\tPressure:%lf solidfication:%lf\tfin\n",sub_state,set_tmp,now_pre,tmp_sij);

                    pos_tmp=functions->read_modestructure();
                    break;
                }
                else
                {
                    pos_tmp=functions->read_MDstructure();
                    printf(".\t\tsub_step:%d\tTemperature:%lf\tPressure:%lf solidfication:%lf\tunstable\n",sub_state,set_tmp,now_pre,tmp_sij);
                    fprintf(fp,".\t\tsub_step:%d\tTemperature:%lf\tPressure:%lf solidfication:%lf\tunstable\n",sub_state,set_tmp,now_pre,tmp_sij);
                }
            }        
            ////////////////////////////////////////////////////////////
            
        }
        if(num_liquids <= 0)
        {
            min_temperature=now_tem;
        }
        if(num_solids <=0 && max_temperature <0)
        {
            max_temperature=now_tem;
        }
        DATA2D_set(liquid_ratio,now_tem,state,0);
        DATA2D_set(liquid_ratio,1.0*num_liquids/(num_solids+num_liquids),state,1);
        /////////////////////////////////////////    
    }
    fclose(fp);

    ///////////////////////////////////////////
    double * para=(double *)calloc(2,sizeof(double));
    para[0]=1.0;
    para[1]=0.5*(min_temperature + max_temperature);
    functionsFit_MSE(liquid_ratio,para,2,func,G_func,500000,0.00001);
    melting_tem=para[1];
    fp=fopen("MyTm_result.dat","a");
    fprintf(fp,"small cell liquid ratio data\n");
    fprintf(fp,"----------------------------------------------------------------------------\n");
    fprintf(fp,"Temperature(K)\tLiquid_ratio\n");
    for(int i=0;i<liquid_ratio->dimen_data[0];i++)
    {
        fprintf(fp,"%lf\t%lf\n",DATA2D_get(liquid_ratio,i,0),DATA2D_get(liquid_ratio,i,1));
    }

    fprintf(fp,"\n\nsmallcell liquid ratio Fit data\n");
    fprintf(fp,"Temperature(K)\tLiquid_ratio\tcentral value:%lf (melting temperature)\n",para[1]);
    int set_num=2000;
    double min_tem=DATA2D_get(liquid_ratio,0,0);
    double max_tem=DATA2D_get(liquid_ratio,liquid_ratio->dimen_data[0] - 1,0);
    double precs_tems=(max_tem - min_tem)/set_num;
    for(int i=0;i<set_num;i++)
    {
        double x=min_tem + i*precs_tems;
        double y=1.0/(1.0 + exp(-para[0]*(x-para[1])));
        fprintf(fp,"%lf\t%lf\n",x,y);
    }

    //////////////////////////////////////////
    free(para);
    DATA_free(liquid_ratio);


    ////////////////////////////////////
    POSCAR_free(pos_solid);
    POSCAR_free(pos_liquid);
    printf(".done\n\n");
    fprintf(fp,".done\n\n");

    printf("4. Melting point calculation: %lf\n",melting_tem);
    fprintf(fp,"----------------------------------------------------------------------------\n");
    fprintf(fp,"IF\tTemperature: %lf\n",melting_tem);
    fclose(fp);

    return 0;
}



double test_run(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions )    
{
   
printf("1\n");
    inputs->so_check=solidfication_init(-1);

    inputs->so_check->solid=init_poscar;
    inputs->so_check->run_type=0;

    double melt=solidfication_check(inputs->so_check,inputs);
printf("2\n");
    inputs->so_check->run_type=1;

    struct POSCAR * pos=POSCAR_in("POSCAR_check");
printf("3\n");
    inputs->so_check->solid=pos;

    melt=solidfication_check(inputs->so_check,inputs);
    printf("%lf\n",melt);

    return 0;
}
