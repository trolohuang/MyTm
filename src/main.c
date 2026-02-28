#include "ctools.h"


int Calculate_melting_point(struct PARAMETERS *inputs);
struct FUNCTIONS * MD_interface_init(void);
void help_message(void);
void version(void);
int main(int argc,char * argv[])
{

    char * optstring="hm:rv";
    int opt;
    opt=getopt(argc,argv,optstring);
    if( opt == 'h' )
    {
        help_message();
    }
    else if( opt == 'm' )
    {
        int method;
        method=atoi(optarg);
        printf("%d\n",method);
        PARA_toFILE(method);
        printf(".\t\tThe input.dat file has been written, please modify the relevant parameters!\n");
        printf("done\n");
    }
    else if(opt == 'r' )
    {
        struct PARAMETERS *inputs=read_parameters_from_file("input.dat");

        Calculate_melting_point(inputs);

        INPUT_free(inputs);

    }
    else if(opt == 'v')
    {
        version();
    }
    else
    {
        help_message();
    }

    return 0;    
}
int Calculate_melting_point(struct PARAMETERS *inputs)
{
    
    int state=0;
    struct FUNCTIONS * functions=MD_interface_init();
    struct POSCAR * pos=POSCAR_in("POSCAR_init");

   
    ///////////////////////////////////////////////// 
    FILE * fp=fopen("MyTm_result.dat","w");
    fprintf(fp,"1.Parameters:\n");
    fprintf(fp,"----------------------------------------------------------------------------\n");
        fprintf(fp,".\tmd_steps %d\n",inputs->vasp_NSW);
        fprintf(fp,".\tA %d\n",inputs->A);
        fprintf(fp,".\tB %d\n",inputs->B);
        fprintf(fp,".\tC %d\n",inputs->C);
        fprintf(fp,".\tmax_step %d\n",inputs->max_step);
        fprintf(fp,".\tmax_sub_step %d\n",inputs->max_sub_step);
        if(fabs(inputs->atom_method[0]) - 0.01 < 0  || fabs(inputs->atom_method[0]   - 1) - 0.01 <0)
        {
            fprintf(fp,".\tmethod 1 qlmi \n"); 
        }
        else if(fabs(inputs->atom_method[0]   - 2) - 0.01 <0)
        {
            fprintf(fp,".\tmethod 2 local atomic density\n"); 
        }
        else if(fabs(inputs->atom_method[0]   - 3) - 0.01 <0)
        {
            fprintf(fp,".\tmethod 3 ML\n"); 
        }
        else if(fabs(inputs->atom_method[0]   - 4) - 0.01 <0)
        {
            fprintf(fp,".\tmethod 4 SOAP\n"); 
        }
        fprintf(fp,".\tsolid_cutoff %lf\n",inputs->decision);
        fprintf(fp,".\tliquid_cutoff %lf\n",inputs->decision_liquid); 
     

    printf("1.Parameters:\n");

        printf(".\tMD infomation:\n");
        printf(".\tmd_steps %d\n",inputs->vasp_NSW);
        printf(".\tinterface file name: %s\n",inputs->script);
        printf(".\tinterface interpreter: %s\n\n",inputs->commands);

        printf(".\tTemperature control:\n");
        printf(".\tmax_step %d\n",inputs->max_step);
        printf(".\tmax_sub_step %d\n",inputs->max_sub_step);

        if(inputs->method == 0)
        {
            printf(".\ttype: heatingUp\n");
        }
        else if(inputs->method == 1 || inputs->method == 2)
        {
            printf(".\ttype: bisection\n");
        }        
        else if(inputs->method == 3)
        {
            printf(".\ttype: smallcell\n");
        }
        printf("\n");



        printf(".\tAtomic classification:\n");
        if(fabs(inputs->atom_method[0]) - 0.01 < 0  || fabs(inputs->atom_method[0]   - 1) - 0.01 <0)
        {
            printf(".\tmethod 1 qlmi OP\n"); 
        }
        else if(fabs(inputs->atom_method[0]   - 2) - 0.01 <0)
        {
            printf(".\tmethod 2 local atomic density\n"); 
        }
        else if(fabs(inputs->atom_method[0]   - 3) - 0.01 <0)
        {
            printf(".\tmethod 3 ML\n"); 
        }
        else if(fabs(inputs->atom_method[0]   - 4) - 0.01 <0)
        {
            printf(".\tmethod 4 SOAP\n"); 
        }
        printf(".\tsolid_cutoff %lf\n",inputs->decision);
        printf(".\tliquid_cutoff %lf\n\n",inputs->decision_liquid); 

    ///////////////////////////////////////////////////////////////////////
    if(1)
    {
        fprintf(fp,".\tMode\n");
        struct MODE * tmp_mode=inputs->modes;
        while(tmp_mode != NULL)
        {
            int type=DATA2D_get(tmp_mode->descripter,0,0);
            if(  type == 0    )
            {
                fprintf(fp,".\tbuild defect %d %d %d\n",(int)DATA2D_get(tmp_mode->descripter,0,1),(int)DATA2D_get(tmp_mode->descripter,0,2),(int)DATA2D_get(tmp_mode->descripter,0,3));
            }
            else if(type == 1)
            {
                fprintf(fp,".\tbuild liquid %d %d\n",(int)DATA2D_get(tmp_mode->descripter,0,1),(int)DATA2D_get(tmp_mode->descripter,0,2));
            }
            fprintf(fp,".\tsize\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",DATA2D_get(tmp_mode->descripter,1,0),DATA2D_get(tmp_mode->descripter,1,1),DATA2D_get(tmp_mode->descripter,1,2),DATA2D_get(tmp_mode->descripter,1,3),DATA2D_get(tmp_mode->descripter,1,4),DATA2D_get(tmp_mode->descripter,1,5));
            fprintf(fp,".\tpos\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",DATA2D_get(tmp_mode->descripter,2,0),DATA2D_get(tmp_mode->descripter,2,1),DATA2D_get(tmp_mode->descripter,2,2),DATA2D_get(tmp_mode->descripter,2,3),DATA2D_get(tmp_mode->descripter,2,4),DATA2D_get(tmp_mode->descripter,2,5));
            fprintf(fp,".\trot\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",DATA2D_get(tmp_mode->descripter,3,0),DATA2D_get(tmp_mode->descripter,3,1),DATA2D_get(tmp_mode->descripter,3,2),DATA2D_get(tmp_mode->descripter,3,3),DATA2D_get(tmp_mode->descripter,3,4),DATA2D_get(tmp_mode->descripter,3,5));
            tmp_mode=tmp_mode->next;
        }

        printf(".\tModel:\n");
        printf(".\tA %d\n",inputs->A);
        printf(".\tB %d\n",inputs->B);
        printf(".\tC %d\n",inputs->C);
        tmp_mode=inputs->modes;
        while(tmp_mode != NULL)
        {
            int type=DATA2D_get(tmp_mode->descripter,0,0);
            if(  type == 0    )
            {
                printf(".\tbuild defect %d %d %d\n",(int)DATA2D_get(tmp_mode->descripter,0,1),(int)DATA2D_get(tmp_mode->descripter,0,2),(int)DATA2D_get(tmp_mode->descripter,0,3));
            }
            else if(type == 1)
            {
                printf(".\tbuild liquid %d %d\n",(int)DATA2D_get(tmp_mode->descripter,0,1),(int)DATA2D_get(tmp_mode->descripter,0,2));
            }
            printf(".\tsize\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",DATA2D_get(tmp_mode->descripter,1,0),DATA2D_get(tmp_mode->descripter,1,1),DATA2D_get(tmp_mode->descripter,1,2),DATA2D_get(tmp_mode->descripter,1,3),DATA2D_get(tmp_mode->descripter,1,4),DATA2D_get(tmp_mode->descripter,1,5));
            printf(".\tpos\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",DATA2D_get(tmp_mode->descripter,2,0),DATA2D_get(tmp_mode->descripter,2,1),DATA2D_get(tmp_mode->descripter,2,2),DATA2D_get(tmp_mode->descripter,2,3),DATA2D_get(tmp_mode->descripter,2,4),DATA2D_get(tmp_mode->descripter,2,5));
            printf(".\trot\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",DATA2D_get(tmp_mode->descripter,3,0),DATA2D_get(tmp_mode->descripter,3,1),DATA2D_get(tmp_mode->descripter,3,2),DATA2D_get(tmp_mode->descripter,3,3),DATA2D_get(tmp_mode->descripter,3,4),DATA2D_get(tmp_mode->descripter,3,5));
            tmp_mode=tmp_mode->next;
        }
    }
    printf(".done\n\n");
    fprintf(fp,".done\n\n");
    //////////////////////////////////////////////////////////////////
    fclose(fp);
    if(inputs->method == 0)
    {
        heatingUP(pos,inputs,functions);
    }
    else if(inputs->method == 1)
    {
        bisection_tem(pos,inputs,functions);
    }
    else if(inputs->method == 2)
    {
        bisection_sta(pos,inputs,functions);
    }
    else if(inputs->method == 3)
    {
        smallcell(pos,inputs,functions);
    } 
    else if(inputs->method == -1)
    {
        test_run(pos,inputs,functions);
    }   
    else
    {
        printf("No Temperature control!\n");
    }
   
    free(functions);

    POSCAR_free(pos);

    return state;
    
}
struct FUNCTIONS * MD_interface_init(void)
{
    struct FUNCTIONS * func=(struct FUNCTIONS *)calloc(1,sizeof(struct FUNCTIONS));
        func->read_initstructure=VASP_read_initstructure;//read init structure
        func->read_modestructure=VASP_read_modestructure;//read mode structure
        func->read_MDstructure=VASP_read_MDstructure;//read afterMD structure
        func->read_MDtrajectory=VASP_read_MDtrajectory;//read MD trajectory
        func->write_modestructure=VASP_write_modestructure;//write mode to file
        func->write_runstructure=VASP_write_runstructure;//write runing structure
        func->MD_relax=VASP_MD_relax;// run MD relax

        func->MD_heatingUP=VASP_MD_heatingUP;// run MD temperure up or down
        func->get_MDavgPressure=VASP_get_MDavgPressure;
        func->get_MDavgTemperature=VASP_get_MDavgTemperature;
        func->get_MDavgEnthalpy=VASP_get_MDavgEnthalpy;
        func->getMDthermodynamic=VASP_getMDthermodynamic;

    return func;
}
void help_message(void)
{
    printf(".Welcome to MyTM\n\n");
    printf(".\t1.What files must be prepared:\n");
    printf(".\t\t(1). Model file containing the smallest unit. (POSCAR_init)\n");
    printf(".\t\t(2). MyTm's control file. (input.dat)\n");
    printf(".\t\t(3). The script for MD program.\n\n");
    printf("\n");

    printf(".\t2.How to generate input.dat file:\n");
    printf(".\t\t$ MyTm -m para1\n");
    printf(".\t\t\tpara1 is the integer number representing the melting point calculation method.\n");
    printf(".\t\t\t0:Direct-heating method;\n");
    printf(".\t\t\t1:void method;\n");
    printf(".\t\t\t2:modified void method;\n");
    printf(".\t\t\t3:two-phase method;\n");
    printf(".\t\t\t4:sandwich method;\n");
    printf(".\t\t\t5:Z method;\n");
    printf(".\t\t\t6:modified Z method;\n");
    printf("\n");

    printf(".\t3.Which variables in the input.dat file need to be manually set:\n");
    printf(".\t\tsupercell a b c;  a, b and c are the expansion multiples of supercells.\n");
    printf(".\t\tmd_script_interpreter and md_script_filename; If you use a bash script named MD.sh to call MD code, please set:\n");
    printf("\t\tmd_script_interpreter bash\n");
    printf("\t\tmd_script_filename ./MD.sh\n");
    printf("\n");

    printf(".\t4.Please refer to usage.pdf for the meanings of all variables and the format of calling MD script\n");

}
void version(void)
{
    printf("MyTm 1.0.0\n");
}