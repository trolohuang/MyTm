#include "ctools.h"

void POSCAR_delAtoms_with_distance(struct POSCAR * pos,struct DATA * atom_protect,double min)
{
    perCrystal_full(pos,0,0,0);
    int state=0;
    struct DATA ** nets=build_coordinations_table_byRadius(pos,-1,min);
    struct DATA * del_table=DATA1D_init(pos->num_all);
    
    int num_d=1;
    while(state==0)
    {
  
        ///////////////////////////////////////////
        int max_at=-99999;
        int max_val=-99999;
        for(int i=0;i<pos->num_all;i++)
        {
            if(del_table->data_data[i] - 0.5 < 0  && atom_protect->data_data[i] - 0.5 <0)
            {
                int n_at=0;
                for(int nt=1;nt<nets[i]->dimen_data[0];nt++)
                {
                    int xui=(int)DATA2D_get(nets[i],nt,7);
                    if(del_table->data_data[xui] - 0.5<0)
                    {
                        n_at++;
                    }
                }

                if(n_at > max_val )
                {
                    max_val=n_at;
                    max_at=i;
                }
            }
        }
        if(max_val > 0)
        {
            del_table->data_data[max_at]=1.0;
        }
        
        /////////////////////////////////////////// 
        num_d=0;
        for(int i=0;i<pos->num_all;i++)
        {
            if( del_table->data_data[i] - 0.5 <0 && atom_protect->data_data[i] - 0.5 <0)
            {
                for(int at=1;at<nets[i]->dimen_data[0];at++)
                {
                    int nei_n=(int)DATA2D_get(nets[i],at,7);
                    if(del_table->data_data[nei_n] - 0.5 <0)
                    {
                        num_d++;
                    }
                }
            }
        }
        ///////////////////////////////////////////
        if(num_d<=0)
        {
            state=1000;
        }
    }

    ///////////////////////////////////////////////
    struct DATA ** atoms=(struct DATA **)calloc(pos->num_type,sizeof(struct DATA *));
    for(int i=0;i<pos->num_type;i++)
    {
        atoms[i]=DATA_init(2,8*(pos->num_all),3);
    }
    int * nums=(int *)calloc(pos->num_type,sizeof(int));
    //////
    int pre_num,end_num;
    int all_atoms_at=0;
    for(int i=0;i<pos->num_type;i++)
    {
        POSCAR_getOrder(pos,i,&pre_num,&end_num);
        for(int at=pre_num;at<end_num;at++)
        {
            double x,y,z;
            POSCAR_getPOS_order(pos,at,&x,&y,&z);
            if(del_table->data_data[at] - 0.5 <0)
            {
                DATA_set(atoms[i],x,nums[i],0);
                DATA_set(atoms[i],y,nums[i],1);
                DATA_set(atoms[i],z,nums[i],2);
                nums[i]++;
                all_atoms_at++;
            }
        }
    } 

    ///////////////
    DATA_free(pos->atom_position);
    pos->atom_position=DATA2D_init(all_atoms_at,3);
    if(pos->atom_position !=NULL)
    {
        DATA_free(pos->atom_position_ext);
        pos->atom_position_ext=NULL;
    }
    ///////////////
    int all_atoms=0;
    for(int i=0;i<pos->num_type;i++)
    {
        pos->num_atoms[i]=nums[i];
        for(int at=0;at<nums[i];at++)
        {
            double x,y,z;
            x=DATA_get(atoms[i],at,0);
            y=DATA_get(atoms[i],at,1);
            z=DATA_get(atoms[i],at,2);
            DATA_set(pos->atom_position,x,all_atoms,0);
            DATA_set(pos->atom_position,y,all_atoms,1);
            DATA_set(pos->atom_position,z,all_atoms,2);
            all_atoms++;
        }
    }
    pos->num_all=all_atoms;
    /////

    for(int i=0;i<pos->num_type;i++)
    {
        DATA_free(atoms[i]);
    }
    free(atoms);
    free(nums);
    DATA_free(del_table);
    DATAs_free(nets);
}

void rote_positions(double x0, double y0,double z0,double rotx,double roty,double rotz,double *rx,double *ry,double *rz)
{
    
    double x1=x0*cos(rotz) - y0*sin(rotz);
    double y1=x0*sin(rotz) + y0*cos(rotz);
    double z1=z0;


    double x2=x1*cos(roty)+z1*sin(roty);
    double y2=y1;
    double z2=-x1*sin(roty)+z1*cos(roty);

    double x3=x2;
    double y3=y2*cos(rotx)-z2*sin(rotx);
    double z3=y2*sin(rotx)+z2*cos(rotx);

    *rx=x3;
    *ry=y3;
    *rz=z3;
}

struct POSCAR * build_defect(struct POSCAR * poscar, struct DATA * descripter)
{
    struct POSCAR * pos=POSCAR_cpy(poscar);
    perCrystal_full(pos,0,0,0);
    POSCAR_TranCoorType(pos,0);
    DATA_free(pos->atom_position_ext);
    pos->atom_position_ext=NULL;
    /////////////////////////////////////////////
    struct DATA * del_table=DATA1D_init(poscar->num_all);

    int shape1=(int )DATA_get(descripter,0,1);
    int shape2=(int )DATA_get(descripter,0,2);
    double defect_scale=DATA_get(descripter,0,3);

    double size1x=DATA_get(descripter,1,0);
    double size1y=DATA_get(descripter,1,1);
    double size1z=DATA_get(descripter,1,2);
    double size2x=DATA_get(descripter,1,3);
    double size2y=DATA_get(descripter,1,4);
    double size2z=DATA_get(descripter,1,5);

    double pos1x=DATA_get(descripter,2,0);
    double pos1y=DATA_get(descripter,2,1);
    double pos1z=DATA_get(descripter,2,2);
    double pos2x=DATA_get(descripter,2,3);
    double pos2y=DATA_get(descripter,2,4);
    double pos2z=DATA_get(descripter,2,5);


    
    double rote1x=DATA_get(descripter,3,0);
    double rote1y=DATA_get(descripter,3,1);
    double rote1z=DATA_get(descripter,3,2);

    double rote2x=DATA_get(descripter,3,3);
    double rote2y=DATA_get(descripter,3,4);
    double rote2z=DATA_get(descripter,3,5);
    

    ///////////////////////////////////////////
    int num_del=0;
    srand(time(NULL));
    for(int at=0;at<pos->num_all;at++)
    {
        double x,y,z;
        POSCAR_getPOS_order(pos,at,&x,&y,&z);
        ///////////////////////////////////////
        int check=0;
        if(shape1 > 0)
        {
            double x1=x-pos1x;
            double y1=y-pos1y;
            double z1=z-pos1z;

            double finx,finy,finz;
            rote_positions(x1,y1,z1,rote1x,rote1y,rote1z,&finx,&finy,&finz);
            if(shape1== 1)//rect
            {
                if(finx + size1x*0.5 > 0 && finx - size1x*0.5 < 0 &&  finy + size1y*0.5 > 0 && finy - size1y*0.5 < 0 && finz + size1z*0.5 > 0 && finz - size1z*0.5 < 0   )
                {
                    check=1;
                }
            }
            else if(shape1 == 2)
            {
                double radius = pow(finx,2)/pow(size1x*0.5,2) + pow(finy,2)/pow(size1y*0.5,2) + pow(finz,2)/pow(size1z*0.5,2) -1;
                if(radius <0)
                {
                    check=1;
                }
            }
            else if(shape1 == 3)
            {
                double radius = pow(finx,2)/pow(size1x*0.5,2) + pow(finy,2)/pow(size1y*0.5,2)-1;
                if(radius < 0 && finz + size1z*0.5 > 0 && finz - size1z*0.5 < 0 )
                {
                    check=1;
                }
            }
            else
            {
                check=0;
            }
        }
        else
        {
            check=0;
        }
        //////////
        if(check == 1)
        {
            if(shape2 > 0)
            {
                double x2=x-pos2x;
                double y2=y-pos2y;
                double z2=z-pos2z;
    
                double finx,finy,finz;
                rote_positions(x2,y2,z2,rote2x,rote2y,rote2z,&finx,&finy,&finz);
                if(shape2== 1)//rect
                {
                    if(!(finx + size2x*0.5 > 0 && finx - size2x*0.5 < 0 &&  finy + size2y*0.5 > 0 && finy - size2y*0.5 < 0 && finz + size2z*0.5 > 0 && finz - size2z*0.5 < 0 )  )
                    {
                        check=2;
                    }
                }
                else if(shape2 == 2)
                {
                    double radius = pow(finx,2)/pow(size2x*0.5,2) + pow(finy,2)/pow(size2y*0.5,2) + pow(finz,2)/pow(size2z*0.5,2) -1;
                    if(radius >0)
                    {
                        check=2;
                    }
                }
                else if(shape2 == 3)
                {
                    double radius = pow(finx,2)/pow(size2x*0.5,2) + pow(finy,2)/pow(size2y*0.5,2) -1;
                    if( !(radius < 0 && finz + size2z*0.5 > 0 && finz - size2z*0.5 < 0 ))
                    {
                        check=2;
                    }
                }
                else
                {
                    check=0;
                }
            }
            else
            {
                check=2;
            }
        }
        ///////////////////////////////////////
        if(check == 2)
        {
            int rgf=rand()%100 + 1;
         
            if( rgf <= (int)(100*defect_scale)  )
            {
                DATA_set(del_table,1.0,at);
                num_del++;
            }

        }
    }   
    
 
    ////////////////////////////////////////////
    struct DATA ** atoms=DATAs_init(pos->num_type,2,pos->num_all,3);
    
    int num_atoms[1000];
    for(int ty=0;ty<pos->num_type;ty++)
    {
        num_atoms[ty]=0;
    }

    for(int ty=0;ty<pos->num_type;ty++)
    {
        int st,et;
        POSCAR_getOrder(pos,ty,&st,&et);
        for(int i=st;i<et;i++)
        {
            if(del_table->data_data[i] - 0.5 <0)
            {
                double x,y,z;
                POSCAR_getPOS_order(pos,i,&x,&y,&z);
                DATA2D_set(atoms[ty],x,num_atoms[ty],0);
                DATA2D_set(atoms[ty],y,num_atoms[ty],1);
                DATA2D_set(atoms[ty],z,num_atoms[ty],2);
                num_atoms[ty]++;
            }
        }
    }

    int num_all=0;
 
    struct DATA * pos_new=DATA2D_init(pos->num_all - num_del,3);
    
    for(int ty=0;ty<pos->num_type;ty++)
    {
        pos->num_atoms[ty]=num_atoms[ty];
        for(int i=0;i<num_atoms[ty];i++)
        {

                double x,y,z;
                
                x=DATA2D_get(atoms[ty],i,0);
                y=DATA2D_get(atoms[ty],i,1);
                z=DATA2D_get(atoms[ty],i,2);
                DATA2D_set(pos_new,x,num_all,0);
                DATA2D_set(pos_new,y,num_all,1);
                DATA2D_set(pos_new,z,num_all,2);
                num_all++;
        }
    }

    /////////////////////////////////////////////
    DATA_free(pos->atom_position);
    pos->atom_position=pos_new;

    pos->num_all=num_all;


    /////////////////////////////////////////////
    DATA_free(del_table);
    DATAs_free(atoms);

    /////////////////////////////////////////////
    return pos;
}

struct POSCAR * build_liquid(struct POSCAR * pos_solid, struct POSCAR * pos_liquid, struct DATA * descripter)
{

    struct POSCAR * pos=POSCAR_cpy(pos_solid);
    struct POSCAR * liquid= POSCAR_cpy(pos_liquid);
    ////////////////////////////////////////////////find the minus distances
    double rdfss=99999;
    perCrystal_full(liquid,0,0,0);
    POSCAR_TranCoorType(liquid,1);
    for(int i=0;i<liquid->num_all;i++)
    {
        double x0,y0,z0;
        POSCAR_getPOS_order(liquid,i,&x0,&y0,&z0);
        for(int j=i+1;j<liquid->num_all;j++)
        {
            double x1,y1,z1;
            POSCAR_getPOS_order(liquid,j,&x1,&y1,&z1);
            double radius=sqrt(  (x1-x0)*(x1-x0) +  (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)    );
            if(radius - rdfss < 0 && radius - 0.1 >0)
            {
                rdfss=radius;
            }
        }
    }
    perCrystal_full(pos,0,0,0);
    POSCAR_TranCoorType(pos,1);
    for(int i=0;i<pos->num_all;i++)
    {
        double x0,y0,z0;
        POSCAR_getPOS_order(pos,i,&x0,&y0,&z0);
        for(int j=i+1;j<pos->num_all;j++)
        {
            double x1,y1,z1;
            POSCAR_getPOS_order(pos,j,&x1,&y1,&z1);
            double radius=sqrt(  (x1-x0)*(x1-x0) +  (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)    );
            if(radius - rdfss < 0 && radius - 0.1 >0)
            {
                rdfss=radius;
            }
        }
    }
    rdfss=rdfss*0.8;
    ////////////////////////////////////////////////

    POSCAR_TranCoorType(pos,0);
    DATA_free(pos->atom_position_ext);
    pos->atom_position_ext=NULL;

    POSCAR_TranCoorType(liquid,0);
    DATA_free(liquid->atom_position_ext);
    liquid->atom_position_ext=NULL;
    /////////////////////////////////////////////// get the shape information 
    struct DATA * del_table=DATA1D_init(pos_solid->num_all);
    struct DATA * in_table=DATA1D_init(pos_liquid->num_all);

    int shape1=(int )DATA2D_get(descripter,0,1);
    int shape2=(int )DATA2D_get(descripter,0,2);

    double size1x=DATA2D_get(descripter,1,0);
    double size1y=DATA2D_get(descripter,1,1);
    double size1z=DATA2D_get(descripter,1,2);
    double size2x=DATA2D_get(descripter,1,3);
    double size2y=DATA2D_get(descripter,1,4);
    double size2z=DATA2D_get(descripter,1,5);

    double pos1x=DATA2D_get(descripter,2,0);
    double pos1y=DATA2D_get(descripter,2,1);
    double pos1z=DATA2D_get(descripter,2,2);
    double pos2x=DATA2D_get(descripter,2,3);
    double pos2y=DATA2D_get(descripter,2,4);
    double pos2z=DATA2D_get(descripter,2,5);

    double rote1x=DATA2D_get(descripter,3,0);
    double rote1y=DATA2D_get(descripter,3,1);
    double rote1z=DATA2D_get(descripter,3,2);
    double rote2x=DATA2D_get(descripter,3,3);
    double rote2y=DATA2D_get(descripter,3,4);
    double rote2z=DATA2D_get(descripter,3,5);


    //////////////////////////////////////////////get the atoms that we need to delet in poscar_solid
    for(int at=0;at<pos->num_all;at++)
    {
        double x,y,z;
        POSCAR_getPOS_order(pos,at,&x,&y,&z);
        ///////////////////////////////////////
        int check=0;
        if(shape1 > 0)
        {
            double x1=x-pos1x;
            double y1=y-pos1y;
            double z1=z-pos1z;

            double finx,finy,finz;
            rote_positions(x1,y1,z1,rote1x,rote1y,rote1z,&finx,&finy,&finz);
            if(shape1== 1)//rect
            {
                if(finx + size1x*0.5 > 0 && finx - size1x*0.5 < 0 &&  finy + size1y*0.5 > 0 && finy - size1y*0.5 < 0 && finz + size1z*0.5 > 0 && finz - size1z*0.5 < 0   )
                {
                    check=1;
                }
            }
            else if(shape1 == 2)
            {
                double radius = pow(finx,2)/pow(size1x*0.5,2) + pow(finy,2)/pow(size1y*0.5,2) + pow(finz,2)/pow(size1z*0.5,2) -1;
                if(radius <0)
                {
                    check=1;
                }
            }
            else if(shape1 == 3)
            {
                double radius = pow(finx,2)/pow(size1x*0.5,2) + pow(finy,2)/pow(size1y*0.5,2)-1;
                if(radius < 0 && finz + size1z*0.5 > 0 && finz - size1z*0.5 < 0 )
                {
                    check=1;
                }
            }
            else
            {
                check=0;
            }
        }
        else
        {
            check=0;
        }
        //////////
        if(check == 1)
        {
            if(shape2 > 0)
            {
                double x2=x-pos2x;
                double y2=y-pos2y;
                double z2=z-pos2z;
    
                double finx,finy,finz;
                rote_positions(x2,y2,z2,rote2x,rote2y,rote2z,&finx,&finy,&finz);
                if(shape2== 1)//rect
                {
                    if(!(finx + size2x*0.5 > 0 && finx - size2x*0.5 < 0 &&  finy + size2y*0.5 > 0 && finy - size2y*0.5 < 0 && finz + size2z*0.5 > 0 && finz - size2z*0.5 < 0 )  )
                    {
                        check=2;
                    }
                }
                else if(shape2 == 2)
                {
                    double radius = pow(finx,2)/pow(size2x*0.5,2) + pow(finy,2)/pow(size2y*0.5,2) + pow(finz,2)/pow(size2z*0.5,2) -1;
                    if(radius >0)
                    {
                        check=2;
                    }
                }
                else if(shape2 == 3)
                {
                    double radius = pow(finx,2)/pow(size2x*0.5,2) + pow(finy,2)/pow(size2y*0.5,2) -1;
                    if( !(radius < 0 && finz + size2z*0.5 > 0 && finz - size2z*0.5 < 0 ))
                    {
                        check=2;
                    }
                }
                else
                {
                    check=0;
                }
            }
            else
            {
                check=2;
            }
        }
        ///////////////////////////////////////
        if(check == 2)
        {

            DATA_set(del_table,1.0,at);

        }
    }  

    //pos - liquid
    //del_table - in_table

    /////////////////////////////////////////////get the atoms that we need to keep in poscar_liquid


    for(int at=0;at<liquid->num_all;at++)
    {
        double x,y,z;
        POSCAR_getPOS_order(liquid,at,&x,&y,&z);
        ///////////////////////////////////////
        int check=0;
        if(shape1 > 0)
        {
            double x1=x-pos1x;
            double y1=y-pos1y;
            double z1=z-pos1z;

            double finx,finy,finz;
            rote_positions(x1,y1,z1,rote1x,rote1y,rote1z,&finx,&finy,&finz);
            if(shape1== 1)//rect
            {
                if(finx + size1x*0.5 > 0 && finx - size1x*0.5 < 0 &&  finy + size1y*0.5 > 0 && finy - size1y*0.5 < 0 && finz + size1z*0.5 > 0 && finz - size1z*0.5 < 0   )
                {
                    check=1;
                }
            }
            else if(shape1 == 2)
            {
                double radius = pow(finx,2)/pow(size1x*0.5,2) + pow(finy,2)/pow(size1y*0.5,2) + pow(finz,2)/pow(size1z*0.5,2) -1;
                if(radius <0)
                {
                    check=1;
                }
            }
            else if(shape1 == 3)
            {
                double radius = pow(finx,2)/pow(size1x*0.5,2) + pow(finy,2)/pow(size1y*0.5,2)-1;
                if(radius < 0 && finz + size1z*0.5 > 0 && finz - size1z*0.5 < 0 )
                {
                    check=1;
                }
            }
            else
            {
                check=0;
            }
        }
        else
        {
            check=0;
        }
        //////////
        if(check == 1)
        {
            if(shape2 > 0)
            {
                double x2=x-pos2x;
                double y2=y-pos2y;
                double z2=z-pos2z;
    
                double finx,finy,finz;
                rote_positions(x2,y2,z2,rote2x,rote2y,rote2z,&finx,&finy,&finz);
                if(shape2== 1)//rect
                {
                    if(!(finx + size2x*0.5 > 0 && finx - size2x*0.5 < 0 &&  finy + size2y*0.5 > 0 && finy - size2y*0.5 < 0 && finz + size2z*0.5 > 0 && finz - size2z*0.5 < 0 )  )
                    {
                        check=2;
                    }
                }
                else if(shape2 == 2)
                {
                    double radius = pow(finx,2)/pow(size2x*0.5,2) + pow(finy,2)/pow(size2y*0.5,2) + pow(finz,2)/pow(size2z*0.5,2) -1;
                    if(radius >0)
                    {
                        check=2;
                    }
                }
                else if(shape2 == 3)
                {
                    double radius = pow(finx,2)/pow(size2x*0.5,2) + pow(finy,2)/pow(size2y*0.5,2) -1;
                    if( !(radius < 0 && finz + size2z*0.5 > 0 && finz - size2z*0.5 < 0 ))
                    {
                        check=2;
                    }
                }
                else
                {
                    check=0;
                }
            }
            else
            {
                check=2;
            }
        }
        ///////////////////////////////////////
        if(check == 2)
        {

            DATA_set(in_table,1.0,at);

        }
    }  
    ///////////////////////////////////////////// combin the solid and liquid atoms
    
    struct DATA ** atoms=DATAs_init(pos->num_type,2,4*(pos->num_all),3);             
    struct DATA ** atoms2=DATAs_init(liquid->num_type,2,4*(liquid->num_all),3);
    int num_atoms[1000];
    int in_atoms[1000];
    for(int ty=0;ty<pos->num_type;ty++)
    {
        num_atoms[ty]=0;
        in_atoms[ty]=0;
    }
    
    int pos_total=0;
    for(int ty=0;ty<pos->num_type;ty++)
    {
        int st,et;
        POSCAR_getOrder(pos,ty,&st,&et);
        for(int i=st;i<et;i++)
        {
            if(del_table->data_data[i] - 0.5 <0)
            {
                double x,y,z;
                POSCAR_getPOS_order(pos,i,&x,&y,&z);
                DATA2D_set(atoms[ty],x,num_atoms[ty],0);
                DATA2D_set(atoms[ty],y,num_atoms[ty],1);
                DATA2D_set(atoms[ty],z,num_atoms[ty],2);
                num_atoms[ty]++;
                pos_total++;
            }
        }
    }

    for(int ty=0;ty<liquid->num_type;ty++)
    {
        int st,et;
        POSCAR_getOrder(liquid,ty,&st,&et);
        for(int i=st;i<et;i++)
        {
            if(in_table->data_data[i] - 0.5 >0)
            {
                double x,y,z;
                POSCAR_getPOS_order(liquid,i,&x,&y,&z);
                DATA2D_set(atoms2[ty],x,in_atoms[ty],0);
                DATA2D_set(atoms2[ty],y,in_atoms[ty],1);
                DATA2D_set(atoms2[ty],z,in_atoms[ty],2);
                in_atoms[ty]++;
                pos_total++;
       
              
            }
        }
    }
    /////////////////////////////////////////////
    int num_all=0;
 
    struct DATA * pos_new=DATA2D_init(pos_total,3);
    struct DATA * pos_protect=DATA1D_init(pos_total);
    
    for(int ty=0;ty<pos->num_type;ty++)
    {
        ////////////////////////////////////////solid
        for(int i=0;i<num_atoms[ty];i++)
        {
                double x,y,z;
                x=DATA2D_get(atoms[ty],i,0);
                y=DATA2D_get(atoms[ty],i,1);
                z=DATA2D_get(atoms[ty],i,2);

                DATA2D_set(pos_new,x,num_all,0);
                DATA2D_set(pos_new,y,num_all,1);
                DATA2D_set(pos_new,z,num_all,2);

                pos_protect->data_data[num_all]=1.0;
                num_all++;
        }
        ///////////////////////////////////////liquid
        for(int i=0;i<in_atoms[ty];i++)
        {
                double x,y,z;
                x=DATA2D_get(atoms2[ty],i,0);
                y=DATA2D_get(atoms2[ty],i,1);
                z=DATA2D_get(atoms2[ty],i,2);

                DATA2D_set(pos_new,x,num_all,0);
                DATA2D_set(pos_new,y,num_all,1);
                DATA2D_set(pos_new,z,num_all,2);

                pos_protect->data_data[num_all]=0;

                num_all++;
        }
        pos->num_atoms[ty]=num_atoms[ty] + in_atoms[ty];
    }
    


    pos->atom_position=pos_new;
    pos->num_all=num_all;
    /////////////////////////////////////////////

    POSCAR_delAtoms_with_distance(pos,pos_protect,rdfss);


    /////////////////////////////////////////////
    POSCAR_free(liquid);

    DATA_free(del_table);
    DATA_free(in_table);

    DATA_free(pos_protect);

    DATAs_free(atoms);
    DATAs_free(atoms2);
    /////////////////////////////////////////////

    return pos;    
}


/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
long long int gcd(long long int a, long long int b)
{
    while (b != 0) {
        long long int tmp = b;
        b = a % b;
        a = tmp;
    }
    return a;
}
// 计算最小公倍数
long long int lcm(long long int a, long long int b)
{
    return (a * b) / gcd(a, b);
}
// 主函数：近似 LCM（忽略 0，并返回倍数因子）
double approx_lcm(double a1, double a2, double a3, int precision, int *m1, int *m2, int *m3)
{
    int scale = (int)pow(10, precision);

    long long int i1 = (a1 != 0.0) ? (long long int)round(a1 * scale) : 0;
    long long int i2 = (a2 != 0.0) ? (long long int)round(a2 * scale) : 0;
    long long int i3 = (a3 != 0.0) ? (long long int)round(a3 * scale) : 0;

    // 统计有效数
    long long int lcm_result = 0;
    if (i1 != 0) lcm_result = i1;
    if (i2 != 0) lcm_result = (lcm_result == 0) ? i2 : lcm(lcm_result, i2);
    if (i3 != 0) lcm_result = (lcm_result == 0) ? i3 : lcm(lcm_result, i3);

    if (lcm_result == 0) {
        // 全部为 0
        *m1 = *m2 = *m3 = 0;
        return 0.0;
    }

    // 计算倍数因子（基于整数）
    *m1 = (i1 != 0) ? (lcm_result / i1) : 0;
    *m2 = (i2 != 0) ? (lcm_result / i2) : 0;
    *m3 = (i3 != 0) ? (lcm_result / i3) : 0;

    return (double)lcm_result / scale;
}
void mat_mult(double * vec,double phi,double theta)
{
    double x=vec[0];
    double y=vec[1];
    double z=vec[2];

    double x1=x*cos(phi)-y*sin(phi);
    double y1=x*sin(phi)+y*cos(phi);
    double z1=z;

    vec[0]=x1*cos(theta)  + z1*sin(theta);
    vec[1]=y1;
    vec[2]=-x1*sin(theta) + z1 * cos(theta);



}

////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
void POSCAR_CrystalType2Orthogonal(struct POSCAR * poscar)
{
    POSCAR_TranCoorType(poscar,0);
    DATA_free(poscar->atom_position_ext);
    poscar->atom_position_ext=NULL;

    double min_volume=9999999999;

    double *vec=(double *)calloc(9,sizeof(double));

    double X1,X2,X3;

    int n1x,n1y,n1z,n2x,n2y,n2z,n3x,n3y,n3z;

    double vec_tmp[9];
    vec_tmp[0]=poscar->Lattice[0][0];
    vec_tmp[1]=poscar->Lattice[0][1];
    vec_tmp[2]=poscar->Lattice[0][2];
    vec_tmp[3]=poscar->Lattice[1][0];
    vec_tmp[4]=poscar->Lattice[1][1];
    vec_tmp[5]=poscar->Lattice[1][2];
    vec_tmp[6]=poscar->Lattice[2][0];
    vec_tmp[7]=poscar->Lattice[2][1];
    vec_tmp[8]=poscar->Lattice[2][2];

    double vec_tmp2[9];

    for(int t=0;t<200;t++)
    {
        for(int p=0;p<200;p++)
        {
            vec[0]=vec_tmp[0];
            vec[1]=vec_tmp[1];
            vec[2]=vec_tmp[2];
            vec[3]=vec_tmp[3];
            vec[4]=vec_tmp[4];
            vec[5]=vec_tmp[5];
            vec[6]=vec_tmp[6];
            vec[7]=vec_tmp[7];
            vec[8]=vec_tmp[8];
            double theta=(2*3.1415926535/200.0) * t;
            double phi=(2*3.1415926535/200.0) * p;
            //////////////////////////////////////
            mat_mult(vec  ,phi,theta);
            mat_mult(vec+3,phi,theta);
            mat_mult(vec+6,phi,theta);

            vec_tmp2[0]=vec[0];
            vec_tmp2[1]=vec[1];
            vec_tmp2[2]=vec[2];
            vec_tmp2[3]=vec[3];
            vec_tmp2[4]=vec[4];
            vec_tmp2[5]=vec[5];
            vec_tmp2[6]=vec[6];
            vec_tmp2[7]=vec[7];
            vec_tmp2[8]=vec[8];           
            ///////////////////////////////////////
            matrix_inverse(vec,3);
            double a1,a2,a3;
            int n1,n2,n3,n4,n5,n6,n7,n8,n9;
            a1=vec[0];
            a2=vec[1];
            a3=vec[2];
            X1=approx_lcm(a1,a2,a3,4,&n1,&n2,&n3);
            a1=(double)n1;
            a2=(double)n2;
            a3=(double)n3;
            double Xn=approx_lcm(a1,a2,a3,0,&n1,&n2,&n3);
            X1=fabs(Xn/X1);

            n1=round(X1*vec[0]);
            n2=round(X1*vec[1]);
            n3=round(X1*vec[2]);
            ///////////////////////////
            a1=vec[3];
            a2=vec[4];
            a3=vec[5];
            X2=approx_lcm(a1,a2,a3,4,&n4,&n5,&n6);
            a1=(double)n4;
            a2=(double)n5;
            a3=(double)n6;
            Xn=approx_lcm(a1,a2,a3,0,&n4,&n5,&n6);
            X2=fabs(Xn/X2);
            n4=round(X2*vec[3]);
            n5=round(X2*vec[4]);
            n6=round(X2*vec[5]);
            ///////////////////////////////
            a1=vec[6];
            a2=vec[7];
            a3=vec[8];
            X3=approx_lcm(a1,a2,a3,4,&n7,&n8,&n9);
            a1=(double)n7;
            a2=(double)n8;
            a3=(double)n9;
            Xn=approx_lcm(a1,a2,a3,0,&n7,&n8,&n9);
            X3=fabs(Xn/X3);

            n7=round(X3*vec[6]);
            n8=round(X3*vec[7]);
            n9=round(X3*vec[8]);

            if(min_volume - X1*X2*X3 >0)
            {

                min_volume=X1*X2*X3;
                poscar->Lattice[0][0]=vec_tmp2[0];
                poscar->Lattice[0][1]=vec_tmp2[1];
                poscar->Lattice[0][2]=vec_tmp2[2];
                poscar->Lattice[1][0]=vec_tmp2[3];
                poscar->Lattice[1][1]=vec_tmp2[4];
                poscar->Lattice[1][2]=vec_tmp2[5];
                poscar->Lattice[2][0]=vec_tmp2[6];
                poscar->Lattice[2][1]=vec_tmp2[7];
                poscar->Lattice[2][2]=vec_tmp2[8];

                n1x=n1;
                n1y=n2;
                n1z=n3;

                n2x=n4;
                n2y=n5;
                n2z=n6;

                n3x=n7;
                n3y=n8;
                n3z=n9;

            }
        }
    }

    double Lattices[3][3];

    Lattices[0][0]=n1x*(poscar->Lattice[0][0]) + n1y*(poscar->Lattice[1][0]) + n1z*(poscar->Lattice[2][0]); 
    Lattices[1][0]=n2x*(poscar->Lattice[0][0]) + n2y*(poscar->Lattice[1][0]) + n2z*(poscar->Lattice[2][0]); 
    Lattices[2][0]=n3x*(poscar->Lattice[0][0]) + n3y*(poscar->Lattice[1][0]) + n3z*(poscar->Lattice[2][0]);

    Lattices[0][1]=n1x*(poscar->Lattice[0][1]) + n1y*(poscar->Lattice[1][1]) + n1z*(poscar->Lattice[2][1]); 
    Lattices[1][1]=n2x*(poscar->Lattice[0][1]) + n2y*(poscar->Lattice[1][1]) + n2z*(poscar->Lattice[2][1]); 
    Lattices[2][1]=n3x*(poscar->Lattice[0][1]) + n3y*(poscar->Lattice[1][1]) + n3z*(poscar->Lattice[2][1]);

    Lattices[0][2]=n1x*(poscar->Lattice[0][2]) + n1y*(poscar->Lattice[1][2]) + n1z*(poscar->Lattice[2][2]); 
    Lattices[1][2]=n2x*(poscar->Lattice[0][2]) + n2y*(poscar->Lattice[1][2]) + n2z*(poscar->Lattice[2][2]); 
    Lattices[2][2]=n3x*(poscar->Lattice[0][2]) + n3y*(poscar->Lattice[1][2]) + n3z*(poscar->Lattice[2][2]);
   
    //////////////////////////////////////////////
    POSCAR_TranCoorType(poscar,1);
    perCrystal_full(poscar,0,0,0);
    struct DATA ** positions=DATAs_init(poscar->num_type,2,10000*(poscar->num_all),3);
    int numtypes[50];
    int num_all=0;
    for(int i=0;i<poscar->num_type;i++)
    {
        numtypes[i]=0;
    }
    ///////
    for(int i=0;i<poscar->num_type;i++)
    {
        int st,et;
        POSCAR_getOrder(poscar,i,&st,&et);
        for(int at=st;at<et;at++)
        {
            double x,y,z;
            POSCAR_getPOS_order(poscar,at,&x,&y,&z);
            double crx,cry,crz;
            coorTurnC(Lattices,1,&crx,&cry,&crz,x,y,z);
            if( crx>=0 && crx-1.0<0 && cry>=0 && cry-1.0<0 && crz>=0&&crz-1.0 <0)
            {
                DATA_set(positions[i],x,numtypes[i],0);
                DATA_set(positions[i],y,numtypes[i],1);
                DATA_set(positions[i],z,numtypes[i],2);
                numtypes[i]++;
                num_all++;
            }
        }
    }
    for(int n=1;n>0;n++)
    {
        int checks=0;
        int sA,eA,sB,eB,sC,eC;
        ////////////////////////////up and down
        sA=-n;
        eA=n;
        sB=-n;
        eB=n;
        sC=n;
        eC=n;
        for(int a=sA;a<=eA;a++)
        {
            for(int b=sB;b<=eB;b++)
            {
                for(int c=sC;c<=eC;c++)
                {
                    perCrystal_full(poscar,a,b,c);
                    for(int i=0;i<poscar->num_type;i++)
                    {
                        int st,et;
                        POSCAR_getOrder(poscar,i,&st,&et);
                        for(int at=st;at<et;at++)
                        {
                            double x,y,z;
                            POSCAR_getPOS_order(poscar,at,&x,&y,&z);
                            double diss=9999;
                            for(int La=-1;La<=1;La++)
                            {
                                for(int Lb=-1;Lb<=1;Lb++)
                                {
                                    for(int Lc=-1;Lc<=1;Lc++)
                                    {
                                        for(int att=0;att<numtypes[i];att++)
                                        {
                                            double tx=DATA_get(positions[i],att,0) + La*X1;
                                            double ty=DATA_get(positions[i],att,1) + Lb*X2;
                                            double tz=DATA_get(positions[i],att,2) + Lc*X3;
                                        
                                            double tmpdiss=sqrt( (x-tx)*(x-tx) +  (y-ty)*(y-ty) + (z-tz)*(z-tz)  );
                                            if(diss-tmpdiss > 0 ) diss=tmpdiss;

                                        }                        
                                    }
                                }
                            }
                            double crx,cry,crz;
                            coorTurnC(Lattices,1,&crx,&cry,&crz,x,y,z);
                            if( crx>=0 && crx-1.0<0 && cry>=0 && cry-1.0<0 && crz>=0&&crz-1.0 <0 && diss-0.5>0)
                            {
                                DATA_set(positions[i],x,numtypes[i],0);
                                DATA_set(positions[i],y,numtypes[i],1);
                                DATA_set(positions[i],z,numtypes[i],2);
                                numtypes[i]++;
                                num_all++;
                                checks++;
                            }
                        }
                    }
                }
            }
        }
        sA=-n;
        eA=n;
        sB=-n;
        eB=n;
        sC=-n;
        eC=-n;
        for(int a=sA;a<=eA;a++)
        {
            for(int b=sB;b<=eB;b++)
            {
                for(int c=sC;c<=eC;c++)
                {
                    perCrystal_full(poscar,a,b,c);
                    for(int i=0;i<poscar->num_type;i++)
                    {
                        int st,et;
                        POSCAR_getOrder(poscar,i,&st,&et);
                        for(int at=st;at<et;at++)
                        {
                            double x,y,z;
                            POSCAR_getPOS_order(poscar,at,&x,&y,&z);
                            double diss=9999;
                            for(int La=-1;La<=1;La++)
                            {
                                for(int Lb=-1;Lb<=1;Lb++)
                                {
                                    for(int Lc=-1;Lc<=1;Lc++)
                                    {
                                        for(int att=0;att<numtypes[i];att++)
                                        {
                                            double tx=DATA_get(positions[i],att,0) + La*X1;
                                            double ty=DATA_get(positions[i],att,1) + Lb*X2;
                                            double tz=DATA_get(positions[i],att,2) + Lc*X3;
                                        
                                            double tmpdiss=sqrt( (x-tx)*(x-tx) +  (y-ty)*(y-ty) + (z-tz)*(z-tz)  );
                                            if(diss-tmpdiss > 0) diss=tmpdiss;

                                        }                        
                                    }
                                }
                            }
                            double crx,cry,crz;
                            coorTurnC(Lattices,1,&crx,&cry,&crz,x,y,z);
                            if( crx>=0 && crx-1.0<0 && cry>=0 && cry-1.0<0 && crz>=0&&crz-1.0 <0 && diss-0.5>0)
                            {
                                DATA_set(positions[i],x,numtypes[i],0);
                                DATA_set(positions[i],y,numtypes[i],1);
                                DATA_set(positions[i],z,numtypes[i],2);
                                numtypes[i]++;
                                num_all++;
                                checks++;
                            }
                        }
                    }
                }
            }
        }
        /////////////////////////////////left and right
        sA=n;
        eA=n;
        sB=-n;
        eB=n;
        sC=-(n-1);
        eC=n-1;
        for(int a=sA;a<=eA;a++)
        {
            for(int b=sB;b<=eB;b++)
            {
                for(int c=sC;c<=eC;c++)
                {
                    perCrystal_full(poscar,a,b,c);
                    for(int i=0;i<poscar->num_type;i++)
                    {
                        int st,et;
                        POSCAR_getOrder(poscar,i,&st,&et);
                        for(int at=st;at<et;at++)
                        {
                            double x,y,z;
                            POSCAR_getPOS_order(poscar,at,&x,&y,&z);
                            double diss=9999;
                            for(int La=-1;La<=1;La++)
                            {
                                for(int Lb=-1;Lb<=1;Lb++)
                                {
                                    for(int Lc=-1;Lc<=1;Lc++)
                                    {
                                        for(int att=0;att<numtypes[i];att++)
                                        {
                                            double tx=DATA_get(positions[i],att,0) + La*X1;
                                            double ty=DATA_get(positions[i],att,1) + Lb*X2;
                                            double tz=DATA_get(positions[i],att,2) + Lc*X3;
                                        
                                            double tmpdiss=sqrt( (x-tx)*(x-tx) +  (y-ty)*(y-ty) + (z-tz)*(z-tz)  );
                                            if(diss-tmpdiss > 0) diss=tmpdiss;

                                        }                        
                                    }
                                }
                            }
                            double crx,cry,crz;
                            coorTurnC(Lattices,1,&crx,&cry,&crz,x,y,z);
                            if( crx>=0 && crx-1.0<0 && cry>=0 && cry-1.0<0 && crz>=0&&crz-1.0 <0 && diss-0.5>0)
                            {
                                DATA_set(positions[i],x,numtypes[i],0);
                                DATA_set(positions[i],y,numtypes[i],1);
                                DATA_set(positions[i],z,numtypes[i],2);
                                numtypes[i]++;
                                num_all++;
                                checks++;
                            }
                        }
                    }
                }
            }
        }
        sA=-n;
        eA=-n;
        sB=-n;
        eB=n;
        sC=-(n-1);
        eC=n-1;
        for(int a=sA;a<=eA;a++)
        {
            for(int b=sB;b<=eB;b++)
            {
                for(int c=sC;c<=eC;c++)
                {
                    perCrystal_full(poscar,a,b,c);
                    for(int i=0;i<poscar->num_type;i++)
                    {
                        int st,et;
                        POSCAR_getOrder(poscar,i,&st,&et);
                        for(int at=st;at<et;at++)
                        {
                            double x,y,z;
                            POSCAR_getPOS_order(poscar,at,&x,&y,&z);
                            double diss=9999;
                            for(int La=-1;La<=1;La++)
                            {
                                for(int Lb=-1;Lb<=1;Lb++)
                                {
                                    for(int Lc=-1;Lc<=1;Lc++)
                                    {
                                        for(int att=0;att<numtypes[i];att++)
                                        {
                                            double tx=DATA_get(positions[i],att,0) + La*X1;
                                            double ty=DATA_get(positions[i],att,1) + Lb*X2;
                                            double tz=DATA_get(positions[i],att,2) + Lc*X3;
                                        
                                            double tmpdiss=sqrt( (x-tx)*(x-tx) +  (y-ty)*(y-ty) + (z-tz)*(z-tz)  );
                                            if(diss-tmpdiss > 0) diss=tmpdiss;

                                        }                        
                                    }
                                }
                            }
                            double crx,cry,crz;
                            coorTurnC(Lattices,1,&crx,&cry,&crz,x,y,z);
                            if( crx>=0 && crx-1.0<0 && cry>=0 && cry-1.0<0 && crz>=0&&crz-1.0 <0 && diss-0.5>0)
                            {
                                DATA_set(positions[i],x,numtypes[i],0);
                                DATA_set(positions[i],y,numtypes[i],1);
                                DATA_set(positions[i],z,numtypes[i],2);
                                numtypes[i]++;
                                num_all++;
                                checks++;
                            }
                        }
                    }
                }
            }
        }     
        ///////////////////////front and behand   
        sA=-(n-1);
        eA=n-1;
        sB=n;
        eB=n;
        sC=-(n-1);
        eC=n-1;
        for(int a=sA;a<=eA;a++)
        {
            for(int b=sB;b<=eB;b++)
            {
                for(int c=sC;c<=eC;c++)
                {
                    perCrystal_full(poscar,a,b,c);
                    for(int i=0;i<poscar->num_type;i++)
                    {
                        int st,et;
                        POSCAR_getOrder(poscar,i,&st,&et);
                        for(int at=st;at<et;at++)
                        {
                            double x,y,z;
                            POSCAR_getPOS_order(poscar,at,&x,&y,&z);
                            double diss=9999;
                            for(int La=-1;La<=1;La++)
                            {
                                for(int Lb=-1;Lb<=1;Lb++)
                                {
                                    for(int Lc=-1;Lc<=1;Lc++)
                                    {
                                        for(int att=0;att<numtypes[i];att++)
                                        {
                                            double tx=DATA_get(positions[i],att,0) + La*X1;
                                            double ty=DATA_get(positions[i],att,1) + Lb*X2;
                                            double tz=DATA_get(positions[i],att,2) + Lc*X3;
                                        
                                            double tmpdiss=sqrt( (x-tx)*(x-tx) +  (y-ty)*(y-ty) + (z-tz)*(z-tz)  );
                                            if(diss-tmpdiss > 0) diss=tmpdiss;

                                        }                        
                                    }
                                }
                            }
                            double crx,cry,crz;
                            coorTurnC(Lattices,1,&crx,&cry,&crz,x,y,z);
                            if( crx>=0 && crx-1.0<0 && cry>=0 && cry-1.0<0 && crz>=0&&crz-1.0 <0 && diss-0.5>0)
                            {
                                DATA_set(positions[i],x,numtypes[i],0);
                                DATA_set(positions[i],y,numtypes[i],1);
                                DATA_set(positions[i],z,numtypes[i],2);
                                numtypes[i]++;
                                num_all++;
                                checks++;
                            }
                        }
                    }
                }
            }
        }
        sA=-(n-1);
        eA=n-1;
        sB=-n;
        eB=-n;
        sC=-(n-1);
        eC=n-1;
        for(int a=sA;a<=eA;a++)
        {
            for(int b=sB;b<=eB;b++)
            {
                for(int c=sC;c<=eC;c++)
                {
                    perCrystal_full(poscar,a,b,c);
                    for(int i=0;i<poscar->num_type;i++)
                    {
                        int st,et;
                        POSCAR_getOrder(poscar,i,&st,&et);
                        for(int at=st;at<et;at++)
                        {
                            double x,y,z;
                            POSCAR_getPOS_order(poscar,at,&x,&y,&z);
                            double diss=9999;
                            for(int La=-1;La<=1;La++)
                            {
                                for(int Lb=-1;Lb<=1;Lb++)
                                {
                                    for(int Lc=-1;Lc<=1;Lc++)
                                    {
                                        for(int att=0;att<numtypes[i];att++)
                                        {
                                            double tx=DATA_get(positions[i],att,0) + La*X1;
                                            double ty=DATA_get(positions[i],att,1) + Lb*X2;
                                            double tz=DATA_get(positions[i],att,2) + Lc*X3;
                                        
                                            double tmpdiss=sqrt( (x-tx)*(x-tx) +  (y-ty)*(y-ty) + (z-tz)*(z-tz)  );
                                            if(diss-tmpdiss > 0) diss=tmpdiss;

                                        }                        
                                    }
                                }
                            }
                            double crx,cry,crz;
                            coorTurnC(Lattices,1,&crx,&cry,&crz,x,y,z);
                            if( crx>=0 && crx-1.0<0 && cry>=0 && cry-1.0<0 && crz>=0&&crz-1.0 <0 && diss-0.5>0)
                            {
                                DATA_set(positions[i],x,numtypes[i],0);
                                DATA_set(positions[i],y,numtypes[i],1);
                                DATA_set(positions[i],z,numtypes[i],2);
                                numtypes[i]++;
                                num_all++;
                                checks++;
                            }
                        }
                    }
                }
            }
        }  


        ////////////////////////////////
        if(checks==0)
        {
            n=-1000;
        }
    }
    ///////////////


    //////////////
    struct DATA * positions_new=DATA_init(2,num_all,3);
    int n=0;
    for(int i=0;i<poscar->num_type;i++)
    {
        poscar->num_atoms[i]=numtypes[i];
        for(int at=0;at<numtypes[i];at++)
        {
            double x,y,z;
            x=DATA_get(positions[i],at,0);
            y=DATA_get(positions[i],at,1);
            z=DATA_get(positions[i],at,2);
            DATA_set(positions_new,x,n,0);
            DATA_set(positions_new,y,n,1);
            DATA_set(positions_new,z,n,2);
            n++;
        }
    }
    DATAs_free(positions);
    DATA_free(poscar->atom_position);
    DATA_free(poscar->atom_position_ext);
    poscar->atom_position=positions_new;
    poscar->atom_position_ext=NULL;
    poscar->num_all=num_all;

    poscar->Lattice[0][0]=Lattices[0][0];
    poscar->Lattice[0][1]=Lattices[0][1];
    poscar->Lattice[0][2]=Lattices[0][2];

    poscar->Lattice[1][0]=Lattices[1][0];
    poscar->Lattice[1][1]=Lattices[1][1];
    poscar->Lattice[1][2]=Lattices[1][2];

    poscar->Lattice[2][0]=Lattices[2][0];
    poscar->Lattice[2][1]=Lattices[2][1];
    poscar->Lattice[2][2]=Lattices[2][2];

}
struct POSCAR * mode_build(struct POSCAR * poscar_solid,struct POSCAR * poscar_liquid, struct MODE * mode_string)
{
    struct MODE * mode_tmp=mode_string;
    struct POSCAR * pos_tmp1=NULL;
    struct POSCAR * pos_build=POSCAR_cpy(poscar_solid);

    while(mode_tmp != NULL)
    {

        int type=(int)DATA_get(mode_tmp->descripter,0,0);
        if(type == 0)
        {
            pos_tmp1=build_defect(pos_build,mode_tmp->descripter);
            POSCAR_free(pos_build);
            pos_build=pos_tmp1;
        }
        else if(type == 1 )
        {
            pos_tmp1=build_liquid(pos_build,poscar_liquid,mode_tmp->descripter);
            POSCAR_free(pos_build);
            pos_build=pos_tmp1;            
        }
        else
        {
            printf("Fun Err! mode_build mode type err!\n");
        }
        mode_tmp=mode_tmp->next;
    }
    return pos_build;
}

void POSCAR_setSuperCell(struct POSCAR * pos,int A,int B,int C)
{

	int all_nums=(pos->num_all)*A*B*C;
  
	struct DATA * positions=DATA2D_init( all_nums ,3   );
	double x1,y1,z1,x2,y2,z2;
	int atoms=0;
	for(int num=0;num<(pos->num_all);num++)
	{

		POSCAR_getPOS_order(pos,num,&x1,&y1,&z1);
		for(int i=0;i<A;i++)
		{
			for(int j=0;j<B;j++)
			{
				for(int k=0;k<C;k++)
				{
                    perCrystal(pos,i,j,k,&x2,&y2,&z2,x1,y1,z1);
					DATA2D_set(positions,x2,atoms,0);
					DATA2D_set(positions,y2,atoms,1);
					DATA2D_set(positions,z2,atoms,2);
					atoms++;
				}
			}
		}

	}
	free(pos->atom_position);
	pos->atom_position=positions;
	pos->num_all=all_nums;

	for(int type=0;type<(pos->num_type);type++)
	{
		pos->num_atoms[type] = A*B*C*(pos->num_atoms[type]);
	}
	/////////
    POSCAR_TranCoorType(pos,1);
    if(pos->atom_position_ext != NULL)
    {
        DATA_free(pos->atom_position_ext);

        pos->atom_position_ext=NULL;
    }
	pos->Lattice[0][0]= (pos->Lattice[0][0])*(1.0*A);
	pos->Lattice[0][1]= (pos->Lattice[0][1])*(1.0*A);
	pos->Lattice[0][2]= (pos->Lattice[0][2])*(1.0*A);

	pos->Lattice[1][0]= (pos->Lattice[1][0])*(1.0*B);
	pos->Lattice[1][1]= (pos->Lattice[1][1])*(1.0*B);
	pos->Lattice[1][2]= (pos->Lattice[1][2])*(1.0*B);

	pos->Lattice[2][0]= (pos->Lattice[2][0])*(1.0*C);
	pos->Lattice[2][1]= (pos->Lattice[2][1])*(1.0*C);
	pos->Lattice[2][2]= (pos->Lattice[2][2])*(1.0*C);

	
}
