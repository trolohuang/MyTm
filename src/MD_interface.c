#include "ctools.h"
struct POSCAR * VASP_read_initstructure(void)//read init structure
{
    return POSCAR_in("POSCAR_init");
}
struct POSCAR * VASP_read_modestructure(void)//read mode structure
{
    return POSCAR_in("POSCAR_mode");
}
double sum_line(const char *filename, int target_line)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
        return 0.0;

    char buf[4096];
    int line = 0;

    while (fgets(buf, sizeof(buf), fp))
    {
        line++;
        if (line == target_line)
            break;
    }

    fclose(fp);

    if (line != target_line)
        return 0.0;      // 文件没有这么多行

    double sum = 0.0;
    char *p = buf;
    char *end;

    while (1)
    {
        double x = strtod(p, &end);

        if (p == end)
            break;       // 没有更多数字

        sum += x;
        p = end;
    }

    return sum;
}
int compare_lines(const char *filename, int line1, int line2)
{
    char cmd[1024];

    snprintf(cmd, sizeof(cmd),
        "awk 'NR==%d{a=$0} NR==%d{exit !(a==$0)}' \"%s\"",
        line1, line2, filename);

    int status = system(cmd);

    if (status == -1)
        return -1;      // system()失败

    if (WIFEXITED(status))
        return WEXITSTATUS(status);

    return -1;
}
struct POSCAR * VASP_read_MDstructure(void)//read afterMD structure
{

    struct POSCAR * pos=NULL;

    int atoms=sum_line("XDATCAR_MyTm",7);

    int type=compare_lines("XDATCAR_MyTm",2,atoms+9);
    if(type == 0)
    {
        char names[256];
        int tat=atoms+8;
        sprintf(names,"tail -n %d XDATCAR_MyTm > POSCAR_tmps",tat);
        system(names);
    }
    else if(type == 1)
    {
        char names[256];

        sprintf(names,"head -n 8 XDATCAR_MyTm > POSCAR_tmps");
        system(names);  

        sprintf(names,"tail -n %d XDATCAR_MyTm >> POSCAR_tmps",atoms);
        system(names);
    }
    pos=POSCAR_in("POSCAR_tmps");
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
    DATA_free(pre);
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
    DATA_free(pre);
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
    struct DATA * res=DATA2D_init(tem->dimen_data[0],3);
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

