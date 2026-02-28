#include "ctools.h"
struct POSCAR * VASP_read_initstructure(void)//read init structure
{
    return POSCAR_in("POSCAR_init");
}
struct POSCAR * VASP_read_modestructure(void)//read mode structure
{
    return POSCAR_in("POSCAR_mode");
}
struct POSCAR * VASP_read_MDstructure(void)//read afterMD structure
{
    int numbers,type;
    struct POSCAR ** poss=XDATCAR_in_advance("XDATCAR_MyTm",10000,&numbers,&type);
    struct POSCAR * pos=POSCAR_cpy(poss[numbers-1]);
    POSCARS_free(poss,numbers);
    return pos;
}
struct POSCAR ** VASP_read_MDtrajectory(int * number_frames)//read MD trajectory
{
    int numbers,type;
    struct POSCAR ** poss=XDATCAR_in_advance("XDATCAR_MyTm",1000,&numbers,&type);
    *number_frames=numbers;
    return poss;
}

void VASP_write_modestructure(struct POSCAR * poscar)//write mode to file
{
    POSCAR_toFILE(poscar,"POSCAR_mode");
}
void VASP_write_runstructure(struct POSCAR * poscar)//write runing structure
{
    POSCAR_toFILE(poscar,"POSCAR_MyTm");
}


void VASP_MD_relax(double TemBeg,double TemEnd,struct PARAMETERS *input,int state,int sub_state)// run MD relax
{
    char commands[500];
    sprintf(commands,"%s %s %lf %lf %d %d %d",input->commands,input->script,TemBeg,TemEnd,input->vasp_NSW,state,sub_state);
    system(commands);
}
void VASP_MD_heatingUP(double TemBeg,double TemEnd,struct PARAMETERS *input,int state,int sub_state)// run MD temperure up or down
{
    char commands[500];
    sprintf(commands,"%s %s %lf %lf %d %d %d",input->commands,input->script,TemBeg,TemEnd,input->vasp_NSW,state,sub_state);
    system(commands);
}

double VASP_get_MDavgPressure(void)
{

    struct DATA * pre=XY_in("Pressure_MyTm");
    double res=0;
    int numbers=0;
    for(int i= (pre->dimen_data[0])*3/4;i<pre->dimen_data[0];i++)
    {
        res+=DATA2D_get(pre,i,1);
        numbers++;
    }
    return res/numbers;
}
double VASP_get_MDavgTemperature(void)
{
    struct DATA * pre=XY_in("Temperature_MyTm");
    double res=0;
    int numbers=0;
    for(int i= (pre->dimen_data[0])*3/4;i<pre->dimen_data[0];i++)
    {
        res+=DATA2D_get(pre,i,1);
        numbers++;
    }
    return res/numbers;
}
double VASP_get_MDavgEnthalpy()
{
    FILE * fp=fopen("Enthalpy_MyTm","r");
    double res=0;
    int steps;
    fscanf(fp,"%d %lf",&steps,&res);
    fclose(fp);
    return res;
}
struct DATA * VASP_getMDthermodynamic()
{
    struct DATA * tem=XY_in("Temperature_MyTm");
    struct DATA * pre=XY_in("Pressure_MyTm");
    struct DATA * vol=XY_in("Volume_MyTm");
    struct DATA * res=DATA_init(2,tem->dimen_data[0],3);
    for(int i=0;i<tem->dimen_data[0];i++)
    {

        DATA2D_set(res,DATA2D_get(tem,i,1),i,0);
        DATA2D_set(res,DATA2D_get(pre,i,1),i,1);
        DATA2D_set(res,DATA2D_get(vol,i,1),i,2);

    }
    DATA_free(tem);
    DATA_free(pre);
    DATA_free(vol);
    return res;
}

