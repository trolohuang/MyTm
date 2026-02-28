#include "ctools.h"
double sum_arithmetic_1(double a1, double d, int n)
{
    return 1.0*n/2*(2*a1 + (n-1)*d);
}
double sum_arithmetic_2(double a1, double an, int n)
{
    return 1.0*n/2*( a1 + an );
}
double integrade_trapezoid(struct DATA * input)
{
    int numbers=input->dimen_data[0];
    double res=0;
    for(int i=0;i<numbers-1;i++)
    {
        double x1=DATA2D_get(input,i,1);
        double x2=DATA2D_get(input,i+1,1);
        double delx=DATA2D_get(input,i+1,0) - DATA2D_get(input,i,0);
        res+= (x1+x2)*delx *0.5; 
    }
    return res;
}
double points2surface(double desx,double desy,double desz,double Vec3_1[],double  Vec3_2[])
{
    double Nvec[3];
    Vec_cross(Nvec,Vec3_1,Vec3_2);
    double nor[3];
    Vec_Nor(nor,Nvec);

    double des[3];
    des[0]=desx;
    des[1]=desy;
    des[2]=desz;
    
    double res=Vec_dot(des,nor);
    if(res<0) res=-res;
    return res;

}
int Qsort_DATA_sub(struct DATA ** data,int left,int right,int dest)
{

    double checks=DATA1D_get(data[right],dest);

    int pos=left;
    
    for(int i=left;i<right;i++)
    {
        double tmps=DATA1D_get(data[i],dest);

        if(tmps-checks < 0)
        {
            struct DATA *tmdata=data[i];
            data[i]=data[pos];
            data[pos]=tmdata;
            pos++;
        }
    }
    struct DATA *tmdata=data[pos];
    data[pos]=data[right];
    data[right]=tmdata;
    return  pos;

}
void Qsort_DATA(struct DATA ** data,int left,int right,int dest)
{
    if(left < right)
    {
        int pos=Qsort_DATA_sub(data,left,right,dest);
        Qsort_DATA(data,left,pos-1,dest);
        Qsort_DATA(data,pos+1,right,dest);
    }

}
struct DATA *  smear_Gaussian(struct DATA * input,double sigma) // sigma < 0.1
{
    int nums=0;
    double delx=DATA_get(input,1,0)-DATA_get(input,0,0);
     nums=(int)(2*sqrt(sigma/2)/delx   );
    /////////////////////////////
    struct DATA * outs=DATA_init(2,input->dimen_data[0],2);
    for(int i=0;i<input->dimen_data[0];i++)
    {
        double num_all=0;
        double average=0;
        for(int ms=-nums;ms<=nums;ms++)
        {
            int num_get= i + ms;
            if(num_get >=0 && num_get < input->dimen_data[0])
            {
                average=average + DATA_get(input,num_get,1)*Gaussian_functions(1.0,ms*delx,0,sigma);
                num_all=num_all + Gaussian_functions(1.0,ms*delx,0,sigma);
            }
        }
        average=average/num_all;
        DATA_set(outs,average,i,1);
        DATA_set(outs,DATA_get(input,i,0),i,0);
    }
    return outs;
}
struct DATA *  smear_linear(struct DATA * input,double sigma) // sigma < 0.1
{
    int nums=(int)((input->dimen_data[0])*sigma);
    struct DATA * outs=DATA_init(2,input->dimen_data[0],2);
    for(int i=0;i<input->dimen_data[0];i++)
    {
        int num_all=0;
        double average=0;
        for(int ms=-nums;ms<=nums;ms++)
        {
            int num_get= i + ms;
            if(num_get >=0 && num_get < input->dimen_data[0])
            {
                num_all++;
                average=average + DATA_get(input,num_get,1);
            }
        }
        average=average/(1.0*num_all);
        DATA_set(outs,average,i,1);
        DATA_set(outs,DATA_get(input,i,0),i,0);
    }
    return outs;
}
struct DATA  * Interpolation_liner(struct DATA * input)
{
    struct DATA * tmp=DATA_init(2,input->dimen_data[0],2);
    int nums=0;
    for(int i=1;i<input->dimen_data[0];i++)
    {
        double dx=0.5*(DATA_get(input,i,0) + DATA_get(input,i-1,0) );
        double dy=0.5*(DATA_get(input,i,1) + DATA_get(input,i-1,1) );

        DATA_set(tmp,dx,nums,0);
        DATA_set(tmp,dy,nums,1);
        nums++;
    }
    struct DATA ** datas=DATAs_init(input->dimen_data[0] + nums,1,2);
    for(int i=0;i<input->dimen_data[0];i++)
    {
        DATA_set(datas[i],DATA_get(input,i,0) ,0);
        DATA_set(datas[i],DATA_get(input,i,1) ,1);
    }
    for(int i=0;i<nums;i++)
    {
        DATA_set(datas[i + input->dimen_data[0]],DATA_get(tmp,i,0) ,0);
        DATA_set(datas[i + input->dimen_data[0]],DATA_get(tmp,i,1) ,1);       
    }
    DATA_free(tmp);
    Qsort_DATA(datas,0,nums+input->dimen_data[0] -1,0  );
    tmp=DATA_init(2,nums+input->dimen_data[0],2);
    for(int i=0;i<nums+input->dimen_data[0];i++)
    {
        DATA_set(tmp,DATA_get(datas[i],0),i,0);
        DATA_set(tmp,DATA_get(datas[i],1),i,1);
    }
    DATAs_free(datas);
    return tmp;

}
void cartesian2spherical(double x,double y,double z,double * r,double *th,double *phi)
{
    *r= sqrt(x*x+y*y+z*z);
    *phi=acos(z / (*r) );
    *th=atan2(y,x);
}
void spherical2cartesian(double *x,double *y,double *z,double  r,double th,double phi)
{
    *x=r*sin(phi)*cos(th);
    *y=r*sin(phi)*sin(th);
    *z=r*cos(phi);
}
double factorial(int n)
{
    if(n == 0 ) return 1;
    else return n*factorial(n-1);
}
double Lengendre(int L,double x)
{
    double re=0;
    if(L==0) re = 1.0;

    else if(L==1) re = x;

    else if(L>1)  re= ((2.0*L -1) *x *Lengendre(L-1,x) -  1.0*(L-1) *Lengendre(L-2,x))/(1.0*L);

    else
    {
        printf("Err: function Legendre err! err L: %d\n",L);
        re=-1;
    } 
    return re;
}
double P_Legendre(int L,int M,double x)
{
    double res=0;
    int MS=abs(M);
    if(MS==0&&L>=0)
    {
        res= Lengendre(L,x);
    }
    else if(L>MS && MS>0)
    {
        double Pmm= pow(-1,MS) * factorial(2*MS)*pow(1-x*x,1.0*MS/2.0)/ (pow(2,MS) *factorial(MS)   )   ; 
        double Pmm1=x*(2*MS + 1)*Pmm;

        double Pmm2=Pmm1;

        for(int K=MS+2;K<=L;K++)
        {
            Pmm2=(    (2*K-1)*x*Pmm1    - 1.0*(K + MS -1)*Pmm                )/  (   1.0*(K-MS)  );
            Pmm=Pmm1;
            Pmm1=Pmm2;
        }
        res=Pmm2;
    }
    else if(L==MS && MS>0)
    {
        res=pow(-1,MS) * factorial(2*MS)*pow(1-x*x,1.0*MS/2.0)/ (pow(2,MS) *factorial(MS)   )   ; 
    }
    else
    {
        printf("Err: function P_Legendre err! err L and M: %d %d\n",L,M);
    }
    ///////////////////////
    if(M<0)
    {
        res=res*pow(-1,MS)*factorial(L-MS)  /   factorial(L+MS);
    }
    return res;
}
void Sph_Har(int L,int M,double theta,double phi,double *re_real,double *re_imag )
{
    *re_real= sqrt(     (2.0*L+1)/(4.0*3.1415926535)  *   factorial(L-M)/factorial(L+M)               ) * P_Legendre(L,M,cos(theta)) * cos(M*phi);
    *re_imag=sqrt(     (2.0*L+1)/(4.0*3.1415926535)  * factorial(L-M)/factorial(L+M)               ) * P_Legendre(L,M,cos(theta)) * sin(M*phi);
}
double Vec_distance(double dx,double dy,double dz)
{
    return sqrt(dx*dx+dy*dy+dz*dz);
}
void Vec_cross(double Vec3_C[],double Vec3_A[],double Vec3_B[])
{
        double ax,ay,az;
        double bx,by,bz;
        ax=Vec3_A[0];
        ay=Vec3_A[1];
        az=Vec3_A[2];
        bx=Vec3_B[0];
        by=Vec3_B[1];
        bz=Vec3_B[2];

        Vec3_C[0]=-az*by+ay*bz;
        Vec3_C[1]=az*bx-ax*bz;
        Vec3_C[2]=-ay*bx+ax*by;

}
double Vec_dot(double Vec3_A[],double Vec3_B[])
{
        double ax,ay,az;
        double bx,by,bz;
        ax=Vec3_A[0];
        ay=Vec3_A[1];
        az=Vec3_A[2];
        bx=Vec3_B[0];
        by=Vec3_B[1];
        bz=Vec3_B[2];
        return ax*bx+ay*by+az*bz;
}
void Vec_Nor(double Vec3_C[],double Vec3_A[])
{

        double ax,ay,az;
        ax=Vec3_A[0];
        ay=Vec3_A[1];
        az=Vec3_A[2];
        double nor=sqrt(ax*ax+ay*ay+az*az);

        Vec3_C[0]=ax/nor;
        Vec3_C[1]=ay/nor;
        Vec3_C[2]=az/nor;
}
double Gaussian_functions(double a,double x,double x0,double t)
{
    return a*pow(M_E,-(x-x0)*(x-x0)/t);
}
void matrix_inverse(double* A, int n)
{
    int *ipiv = (int*)malloc(n * sizeof(int));
    int lwork = n * n;
    double *work = (double*)malloc(lwork * sizeof(double));
    int info;
    
    // LU 分解
    dgetrf(&n, &n, A, &n, ipiv, &info);
    
    if (info != 0) {
        printf("LU factorization failed with code %d\n", info);
        free(ipiv);
        free(work);
        return;
    }
    
    // 计算逆矩阵
    dgetri(&n, A, &n, ipiv, work, &lwork, &info);
    
    if (info != 0) {
        printf("Matrix inversion failed with code %d\n", info);
    }
    
    free(ipiv);
    free(work);
}
struct DATA *  AUC_ROC_cure(struct DATA * fx,struct DATA * fx_state,double min_la,double max_la,int bin_la,double * auc)
{

    double prec = (max_la - min_la)/(1.0*bin_la);
   
    struct DATA * ret=DATA2D_init(bin_la,3);
   
    struct DATA ** tmps=DATAs_init(bin_la,1,3);

    for(int i=0;i<bin_la;i++)
    {
        double TP=0,FP=0,TN=0,FN=0;
        double la_val=min_la + prec*i;
        for(int m=0;m<fx->dimen_data[0];m++)
        {
            double f_val=DATA2D_get(fx,m,1);
            double state=DATA1D_get(fx_state,m);
            if(f_val - la_val >0 && state - 0.5 >0)
            {
                TP+=1.0;
            }
            if(f_val - la_val >0 && state - 0.5 <0)
            {
                FP+=1.0;
            }
            if(f_val - la_val <= 0 && state - 0.5 < 0)
            {
                TN+=1.0;
            }
            if(f_val - la_val <= 0 && state - 0.5 >0)
            {
                FN+=1.0;
            }
        }
        //////////////////////////////
        double roc_x=FP/(FP+TN);
        double roc_y=TP/(TP+FN);

        DATA_set(tmps[i],roc_x,0);
        DATA_set(tmps[i],roc_y,1);
        DATA_set(tmps[i],la_val,2);
    }

    Qsort_DATA(tmps,0,bin_la-1,0);

    for(int i=0;i<bin_la;i++)
    {
        DATA2D_set(ret,tmps[i]->data_data[0],i,0);
        DATA2D_set(ret,tmps[i]->data_data[1],i,1);
        DATA2D_set(ret,tmps[i]->data_data[2],i,2);
    }
    DATAs_free(tmps);

    if(auc != NULL)
    {
        double val=integrade_trapezoid(ret);
        *auc=val;
    }

    return ret;


}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// sjk check
struct DATA ** build_coordinations_table_bySANN(struct POSCAR * input,int atom_type)
{
    int MAX_atoms=10000;
    perCrystal_full(input,0,0,0);
    POSCAR_TranCoorType(input,1);
    double radius=3.0;
    /////////////////////////////////////////////////////////////////////////////// target atoms
    int pre_num=0;
    int end_num=0;
    if(atom_type >=0)
    {
        
        for(int i=0;i<atom_type;i++)
        {
            pre_num+= input->num_atoms[i];
        }
        end_num=pre_num+(input->num_atoms[atom_type]);
    }
    else
    {
        pre_num=0;
        end_num=input->num_all;
    }
    ////////////////////////////////////////////////////////////////////////////cal max radius
    double desx,desy,desz;
    //////////////////////////////////////////////////////////////////////////////init table
   
    struct DATA ** nets=(struct DATA **)calloc(end_num-pre_num,sizeof(struct DATA *));
    struct DATA **  data=(struct DATA ** )calloc(MAX_atoms,sizeof(struct DATA *));

    for(int i=0;i<MAX_atoms;i++)
    {
        data[i]=DATA1D_init(8); //neighbor    x y z radius A B C num
    }

    double vec3_a[3];
    double vec3_b[3];
    double vec3_c[3];
    vec3_a[0]=input->Lattice[0][0];
    vec3_a[1]=input->Lattice[0][1];
    vec3_a[2]=input->Lattice[0][2];
    vec3_b[0]=input->Lattice[1][0];
    vec3_b[1]=input->Lattice[1][1];
    vec3_b[2]=input->Lattice[1][2];
    vec3_c[0]=input->Lattice[2][0];
    vec3_c[1]=input->Lattice[2][1];
    vec3_c[2]=input->Lattice[2][2];
    double maxA=points2surface(input->Lattice[0][0],input->Lattice[0][1],input->Lattice[0][2],vec3_b,vec3_c);
    double maxB=points2surface(input->Lattice[1][0],input->Lattice[1][1],input->Lattice[1][2],vec3_a,vec3_c);
    double maxC=points2surface(input->Lattice[2][0],input->Lattice[2][1],input->Lattice[2][2],vec3_a,vec3_b);


    for(int atom=pre_num;atom<end_num;atom++)
    {
        
        int number_neighbor=0;
        int num=atom-pre_num;
        int quites=1;
        int AS,AE,BS,BE,CS,CE;
        POSCAR_getPOS_order(input,atom,&desx,&desy,&desz);

       //////////////////////////////////////////////////////////////////////////////  cal the minimal A B C

         //for AB
        double min1=points2surface(desx,desy,desz,vec3_a,vec3_b);
        double min2=maxC-min1;
        //forBC
         double min3=points2surface(desx,desy,desz,vec3_b,vec3_c);
         double min4=maxA-min3;
        //forAC
        double min5=points2surface(desx,desy,desz,vec3_a,vec3_c);
        double min6=maxB-min5;
       
        
       while(quites == 1)
        {
           if(min3-radius<=0) AS=-1;
           else AS=0;
           if(min4-radius<=0) AE=1;
           else AE=0;
           if(min5 -radius<=0) BS=-1;
           else BS=0;
           if(min6-radius<=0) BE=1;
           else BE=0;
           if(min1-radius<=0) CS=-1;
           else CS=0;
           if(min2-radius<=0) CE=1;
           else CE=0; 

           number_neighbor=0;
           ////////////////////////////////////////////////////////////////////////////// cal neighbor atoms
           
           for(int A=AS;A<=AE;A++) //A
           {
               for(int B=BS;B<=BE;B++) // B
               {
                   for(int C=CS;C<=CE;C++) // C
                   {
                           for(int i=pre_num;i<end_num;i++)
                           {
                               double tx,ty,tz,mtx,mty,mtz,tr;
                               POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                               perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                               tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );
                               if(tr-radius <= 0 && tr-0.05 >0)
                               {
                                   DATA1D_set(data[number_neighbor],tx,0);
                                   DATA1D_set(data[number_neighbor],ty,1);
                                   DATA1D_set(data[number_neighbor],tz,2);
                                   DATA1D_set(data[number_neighbor],tr,3);

                                   DATA1D_set(data[number_neighbor],(double)A,4);
                                   DATA1D_set(data[number_neighbor],(double)B,5);
                                   DATA1D_set(data[number_neighbor],(double)C,6);

                                   DATA1D_set(data[number_neighbor],(double)(i-pre_num),7);
                                   number_neighbor++;
                               }
                           }
                   }
               }
           }
           //////////////////////////////////////////////////////////////////////////////// cal coordinates atoms
           
           if(number_neighbor >3)
           {
               Qsort_DATA(data,0, number_neighbor - 1, 3  );
               /////////////////////////////////////////////////
               int numbers=0;
               double Rmi;
               for(int i=3;i<number_neighbor;i++)
               {
                   numbers=i;
                   Rmi=0;
                   for(int m=0;m<i;m++)
                   {
                       Rmi=Rmi+DATA1D_get(data[m],3);
                   }
                   Rmi=Rmi/(numbers-2);

                   double checks=DATA1D_get(data[i],3);
                   if(Rmi - checks < 0)
                   {    
                       break;
                   }
               }

               quites=0;
               if(numbers    <   number_neighbor - 1)
               {
                        nets[num]=DATA2D_init(numbers + 1 ,8);
                        nets[num]->data_group=end_num-pre_num;
                        DATA2D_set(nets[num],desx,0,0);
                        DATA2D_set(nets[num],desy,0,1);
                        DATA2D_set(nets[num],desz,0,2);
                        DATA2D_set(nets[num],numbers,0,3);
                       for(int i=1;i< numbers+1;i++)
                       {
                           double x,y,z,r,A,B,C,n;
                           x = DATA1D_get(data[i-1],0    );
                           y = DATA1D_get(data[i-1],1    );
                           z = DATA1D_get(data[i-1],2    );
                           r = DATA1D_get(data[i-1],3    );
                           A = DATA1D_get(data[i-1],4    );
                           B = DATA1D_get(data[i-1],5    );
                           C = DATA1D_get(data[i-1],6    );
                           n = DATA1D_get(data[i-1],7    );
                           DATA2D_set(nets[num],x,i,0);
                           DATA2D_set(nets[num],y,i,1);
                           DATA2D_set(nets[num],z,i,2);
                           DATA2D_set(nets[num],r,i,3);
                           DATA2D_set(nets[num],A,i,4);
                           DATA2D_set(nets[num],B,i,5);
                           DATA2D_set(nets[num],C,i,6);
                           DATA2D_set(nets[num],n,i,7);
                       }
                        quites=0;
               }
               else
               {
                    quites=1;
                    radius=radius+2;
                    //printf("%d  Warning function build_coordinations_table Err!  max_radius is too small, we have set max_radius=2.0*max_radius!\n",atom);
               }
           }
           else
           {
                quites=1;
                radius=radius+2;
                //printf("%d  Warning function build_coordinations_table Err!  max_radius is too small, we have set max_radius=2.0*max_radius!\n",atom);
           }

           
        }
        

    }
     ///////////////////////////////////////////////////free memery;
    
    
    
    for(int i=0;i<MAX_atoms;i++)
    {
        DATA_free(data[i]);
    }
    free(data);

    return nets;
}

struct DATA ** build_coordinations_table_byRadius_typedecompose(struct POSCAR * input,double Rcutoff)
{

    if(input->coordinate_type ==0)
    {
        perCrystal_full(input,0,0,0);
        POSCAR_TranCoorType(input,1);
        perCrystal_full(input,0,0,0);
    }
    else
    {
        perCrystal_full(input,0,0,0);
        POSCAR_TranCoorType(input,0);
        perCrystal_full(input,0,0,0);
        POSCAR_TranCoorType(input,1);
    }
    ///////////////////////////////////////////////////////////////////////////////
    double mindis=9999999;
    double disa=points2surface(input->Lattice[0][0],input->Lattice[0][1],input->Lattice[0][2],input->Lattice[1],input->Lattice[2]);
    if(disa - mindis <0) mindis=disa;
    double disb=points2surface(input->Lattice[1][0],input->Lattice[1][1],input->Lattice[1][2],input->Lattice[0],input->Lattice[2]);
    if(disb - mindis <0) mindis=disb;
    double disc=points2surface(input->Lattice[2][0],input->Lattice[2][1],input->Lattice[2][2],input->Lattice[0],input->Lattice[1]);
    if(disc - mindis <0) mindis=disc;

    ///////////////////////////////////////////////////////////////////////////////
    int num_type=input->num_type;
    int num_all=input->num_all;
    double radius=Rcutoff;
    struct DATA ** nets=(struct DATA **)calloc(num_type * num_all,sizeof(struct DATA *));
    ///////////////////////////////////////////////////////////////////////////////
    double volume=POSCAR_getVolume(input->Lattice);
    double density=volume/input->num_all;
    int MAX_atoms=(int)(2*pow(2.5*Rcutoff,3)/density);
    //struct DATA *  data=DATA2D_init(MAX_atoms,8);
    
    int sizea=(int)disa/Rcutoff;
    int sizeb=(int)disa/Rcutoff;
    int sizec=(int)disa/Rcutoff;    
    ///////////////////////////////////////////////////////////////////////////////
    // disa - 2*Rcutoff < 0 && disb - 2*Rcutoff < 0 && disc - 2*Rcutoff < 0
    if(sizea < 2  || sizeb< 2 ||sizec<  2   )
    {

        
        double * data=(double *)calloc(9*MAX_atoms,sizeof(double));

        ////////////////////////////////////////////////////////////////////////////cal max radius
        //////////////////////////////////////////////////////////////////////////////init table
        /////////////////////////////////////////////////////////////////////////////cal MAX A, B and C
        int A,B,C;
        for(int atom=0;atom<num_all;atom++)
        {
            
            double desx,desy,desz;
            POSCAR_getPOS_order(input,atom,&desx,&desy,&desz);
            for(int types=0;types<num_type;types++)
            {
                
                int number_neighbor=1;
                data[0]=desx;
                data[1]=desy;
                data[2]=desz;
                data[8]=POSCAR_getAtomType(input,atom);
                int pre_num=0;
                int end_num=0;
                
                POSCAR_getOrder(input,types,&pre_num,&end_num);
                //////////////////////////////////////////////

                for(int i=pre_num;i<end_num;i++)
                {
                    double mtx,mty,mtz,tr;
                    POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                    tr=sqrt(      (mtx-desx)*(mtx-desx)   + (mty-desy)*(mty-desy)   + (mtz-desz)*(mtz-desz)     );
                    if(tr-radius < 0 && tr-0.005 >0)
                    {
                        
                        data[number_neighbor*9 + 0]=mtx;
                        data[number_neighbor*9 + 1]=mty;
                        data[number_neighbor*9 + 2]=mtz;
                        data[number_neighbor*9 + 3]=tr;
                        data[number_neighbor*9 + 4]=(double)A;
                        data[number_neighbor*9 + 5]=(double)B;
                        data[number_neighbor*9 + 6]=(double)C;
                        data[number_neighbor*9 + 7]=(double)(i*num_type);
                        data[number_neighbor*9 + 8]=(double)types;
                        number_neighbor++;
                    }
                }
                int Nn=1;
                while(Nn>0)
                {
                    int num_single=0;
                    ////////////////////////////////////////////////////////////////////////////// cal neighbor atoms
                    /////for up C=n
                    C=Nn;
                    for(A=-Nn;A<=Nn;A++) 
                    {
                        for(B=-Nn;B<=Nn;B++) 
                        {
                                for(int i=pre_num;i<end_num;i++)
                                {
                                    double tx,ty,tz,mtx,mty,mtz,tr;
                                    POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                                    perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                                    tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );
                                    if(tr-radius < 0 && tr-0.05 >0)
                                    {
                                        data[number_neighbor*9 + 0]=tx;
                                        data[number_neighbor*9 + 1]=ty;
                                        data[number_neighbor*9 + 2]=tz;
                                        data[number_neighbor*9 + 3]=tr;
                                        data[number_neighbor*9 + 4]=(double)A;
                                        data[number_neighbor*9 + 5]=(double)B;
                                        data[number_neighbor*9 + 6]=(double)C;
                                        data[number_neighbor*9 + 7]=(double)(i*num_type);
                                        data[number_neighbor*9 + 8]=(double)types;
                                        number_neighbor++;
                                        num_single++;
                                    }
                                }
                        }
                    }
                    /////for down C=-n
                    C=-Nn;
                    for(A=-Nn;A<=Nn;A++) 
                    {
                        for(B=-Nn;B<=Nn;B++) 
                        {
                                for(int i=pre_num;i<end_num;i++)
                                {
                                    double tx,ty,tz,mtx,mty,mtz,tr;
                                    POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                                    perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                                    tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                                    if(tr-radius < 0 && tr-0.05 >0)
                                    {
                                        data[number_neighbor*9 + 0]=tx;
                                        data[number_neighbor*9 + 1]=ty;
                                        data[number_neighbor*9 + 2]=tz;
                                        data[number_neighbor*9 + 3]=tr;
                                        data[number_neighbor*9 + 4]=(double)A;
                                        data[number_neighbor*9 + 5]=(double)B;
                                        data[number_neighbor*9 + 6]=(double)C;
                                        data[number_neighbor*9 + 7]=(double)(i*num_type);
                                        data[number_neighbor*9 + 8]=(double)types;
                                        number_neighbor++;
                                        num_single++;
                                    }
                                }
                        }
                    }
                    /////for left B=-n
                    B=-Nn;
                    for(A=-Nn;A<=Nn;A++) 
                    {
                        for(C=-(Nn-1);C<=(Nn-1);C++) 
                        {
                                for(int i=pre_num;i<end_num;i++)
                                {
                                    double tx,ty,tz,mtx,mty,mtz,tr;
                                    POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                                    perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                                    tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                                    if(tr-radius < 0 && tr-0.05 >0)
                                    {
                                        data[number_neighbor*9 + 0]=tx;
                                        data[number_neighbor*9 + 1]=ty;
                                        data[number_neighbor*9 + 2]=tz;
                                        data[number_neighbor*9 + 3]=tr;
                                        data[number_neighbor*9 + 4]=(double)A;
                                        data[number_neighbor*9 + 5]=(double)B;
                                        data[number_neighbor*9 + 6]=(double)C;
                                        data[number_neighbor*9 + 7]=(double)(i*num_type);
                                        data[number_neighbor*9 + 8]=(double)types;
                                        number_neighbor++;
                                        num_single++;
                                    }
                                }
                        }
                    }
                    /////for right B=-n
                    B=Nn;
                    for(A=-Nn;A<=Nn;A++) 
                    {
                        for(C=-(Nn-1);C<=(Nn-1);C++) 
                        {
                                for(int i=pre_num;i<end_num;i++)
                                {
                                    double tx,ty,tz,mtx,mty,mtz,tr;
                                    POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                                    perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                                    tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                                    if(tr-radius < 0 && tr-0.05 >0)
                                    {
                                        data[number_neighbor*9 + 0]=tx;
                                        data[number_neighbor*9 + 1]=ty;
                                        data[number_neighbor*9 + 2]=tz;
                                        data[number_neighbor*9 + 3]=tr;
                                        data[number_neighbor*9 + 4]=(double)A;
                                        data[number_neighbor*9 + 5]=(double)B;
                                        data[number_neighbor*9 + 6]=(double)C;
                                        data[number_neighbor*9 + 7]=(double)(i*num_type);
                                        data[number_neighbor*9 + 8]=(double)types;
                                        number_neighbor++;
                                        num_single++;
                                    }
                                }
                        }
                    }      
                    /////for front A=n
                    A=Nn;
                    for(B=-(Nn-1);B<=(Nn-1);B++) 
                    {
                        for(C=-(Nn-1);C<=(Nn-1);C++) 
                        {
                                for(int i=pre_num;i<end_num;i++)
                                {
                                    double tx,ty,tz,mtx,mty,mtz,tr;
                                    POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                                    perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                                    tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                                    if(tr-radius < 0 && tr-0.05 >0)
                                    {
                                        data[number_neighbor*9 + 0]=tx;
                                        data[number_neighbor*9 + 1]=ty;
                                        data[number_neighbor*9 + 2]=tz;
                                        data[number_neighbor*9 + 3]=tr;
                                        data[number_neighbor*9 + 4]=(double)A;
                                        data[number_neighbor*9 + 5]=(double)B;
                                        data[number_neighbor*9 + 6]=(double)C;
                                        data[number_neighbor*9 + 7]=(double)(i*num_type);
                                        data[number_neighbor*9 + 8]=(double)types;
                                        number_neighbor++;
                                        num_single++;
                                    }
                                }
                        }
                    }     
                    /////for behand A=n
                    A=-Nn;
                    for(B=-(Nn-1);B<=(Nn-1);B++) 
                    {
                        for(C=-(Nn-1);C<=(Nn-1);C++) 
                        {
                                for(int i=pre_num;i<end_num;i++)
                                {
                                    double tx,ty,tz,mtx,mty,mtz,tr;
                                    POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                                    perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                                    tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                                    if(tr-radius < 0 && tr-0.05 >0)
                                    {
                                        data[number_neighbor*9 + 0]=tx;
                                        data[number_neighbor*9 + 1]=ty;
                                        data[number_neighbor*9 + 2]=tz;
                                        data[number_neighbor*9 + 3]=tr;
                                        data[number_neighbor*9 + 4]=(double)A;
                                        data[number_neighbor*9 + 5]=(double)B;
                                        data[number_neighbor*9 + 6]=(double)C;
                                        data[number_neighbor*9 + 7]=(double)(i*num_type);
                                        data[number_neighbor*9 + 8]=(double)types;
                                        number_neighbor++;
                                        num_single++;
                                    }
                                }
                        }
                    }   

                    ////////////////////////////////////////////////////////////////////////////// reset atom pos
                    if(num_single==0) Nn=-1000;  
                    Nn++;
                }
                
                nets[atom*num_type + types]=DATA2D_init(number_neighbor ,9);
                nets[atom*num_type + types]->data_group=num_all*num_type;
                data[3]=number_neighbor-1;
                for(int i=0;i< number_neighbor;i++)
                {
                    DATA2D_set(nets[atom*num_type + types],data[i*9+0],i,0);
                    DATA2D_set(nets[atom*num_type + types],data[i*9+1],i,1);
                    DATA2D_set(nets[atom*num_type + types],data[i*9+2],i,2);
                    DATA2D_set(nets[atom*num_type + types],data[i*9+3],i,3);
                    DATA2D_set(nets[atom*num_type + types],data[i*9+4],i,4);
                    DATA2D_set(nets[atom*num_type + types],data[i*9+5],i,5);
                    DATA2D_set(nets[atom*num_type + types],data[i*9+6],i,6);
                    DATA2D_set(nets[atom*num_type + types],data[i*9+7],i,7);
                    DATA2D_set(nets[atom*num_type + types],data[i*9+8],i,8);
                }
                
            }
        
        }
        
        free(data);
    }
    else
    {
        POSCAR_TranCoorType(input,0);
        POSCAR_TranCoorType(input,1);
        ////////////////////////////////////////
        struct DATA ** table=DATAs_init(sizea*sizeb*sizec,2,num_all,2); // cell 
        struct DATA *  n_table=DATA1D_init(sizea*sizeb*sizec);        // number of cell
        struct DATA *  atom_pos_cell=DATA2D_init(num_all,3);          //atoms in cell pos
        for(int i=0;i<num_all;i++)
        {
            double desx,desy,desz;
            POSCAR_getPOS_order(input,i,&desx,&desy,&desz);
            int posa=(int)(   DATA2D_get(input->atom_position_ext,i,0)  *sizea      );   
            int posb=(int)(   DATA2D_get(input->atom_position_ext,i,1)  *sizeb      );
            int posc=(int)(   DATA2D_get(input->atom_position_ext,i,2)  *sizec      );   

            DATA2D_set(atom_pos_cell,posa,i,0); 
            DATA2D_set(atom_pos_cell,posb,i,1); 
            DATA2D_set(atom_pos_cell,posc,i,2);  
            if(posa >=sizea) posa=sizea-1;
            if(posb >=sizeb) posb=sizeb-1;
            if(posc >=sizec) posc=sizec-1;
            DATA2D_set(table[ posc*sizea*sizeb + posb*sizea + posa ], i, (int)(n_table->data_data[posc*sizea*sizeb + posb*sizea + posa] ) , 0      )   ;
            int type=POSCAR_getAtomType(input,i);
            DATA2D_set(table[ posc*sizea*sizeb + posb*sizea + posa ], type, (int)(n_table->data_data[posc*sizea*sizeb + posb*sizea + posa] ) , 1      )   ;
            n_table->data_data[posc*sizea*sizeb + posb*sizea + posa]+=1;
        }
        ////////////////////////////////////////
        struct DATA ** data=DATAs_init( num_type,2,num_all,9);
        int *          n_data=(int *)calloc(num_type,sizeof(int));
        int trans_a=0;
        int trans_b=0;
        int trans_c=0;
        int posa;
        int posb;
        int posc;
        for(int at=0;at<num_all;at++)
        {
            memset(n_data,0,num_type*sizeof(int));
            double desx,desy,desz;
            POSCAR_getPOS_order(input,at,&desx,&desy,&desz);
            for(int dt=0;dt<num_type;dt++)
            {
                DATA2D_set(data[dt],desx,n_data[dt],0);
                DATA2D_set(data[dt],desy,n_data[dt],1);
                DATA2D_set(data[dt],desz,n_data[dt],2);  
                DATA2D_set(data[dt],POSCAR_getAtomType(input,at),n_data[dt],8); 
                n_data[dt]++;             
            }

            posa=DATA2D_get(atom_pos_cell,at,0);
            posb=DATA2D_get(atom_pos_cell,at,1);
            posc=DATA2D_get(atom_pos_cell,at,2);        

            for(int ca=-1;ca<=1;ca++)
            {
                for(int cb=-1;cb<=1;cb++)
                {
                    for(int cc=-1;cc<=1;cc++)
                    {
                        trans_a=0;
                        trans_b=0;
                        trans_c=0;
                        int nowa=posa+ca;
                        int nowb=posb+cb;
                        int nowc=posc+cc;
                        if(nowa  >=sizea)
                        {
                            nowa=0;
                            trans_a=1;
                        }
                        else if(nowa<0)
                        {
                            nowa=sizea-1;
                            trans_a=-1;
                        }
                        if(nowb  >=sizeb)
                        {
                            nowb=0;
                            trans_b=1;
                        }
                        else if(nowb<0)
                        {
                            nowb=sizeb-1;
                            trans_b=-1;
                        }
                        if(nowc  >=sizec)
                        {
                            nowc=0;
                            trans_c=1;
                        }
                        else if(nowc<0)
                        {
                            nowc=sizec-1;
                            trans_c=-1;
                        }

                        int real_posn=nowc*sizea*sizeb + nowb*sizea + nowa;
                        for(int nei=0;nei<(int)(n_table->data_data[real_posn] ) ;nei++ )
                        {
                            int posn=DATA2D_get(table[ real_posn ], nei, 0      );
                            int type=DATA2D_get(table[ real_posn ], nei, 1      );

                            double tx,ty,tz,mtx,mty,mtz,tr;
                            POSCAR_getPOS_order(input,posn,&mtx,&mty,&mtz);
                            perCrystal(input,trans_a,trans_b,trans_c,&tx,&ty,&tz,mtx,mty,mtz);
                            tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                            if(tr-radius < 0 && tr-0.05 >0)
                            {
                                DATA2D_set(data[type],tx,                     n_data[type],0);
                                DATA2D_set(data[type],ty,                     n_data[type],1);
                                DATA2D_set(data[type],tz,                     n_data[type],2);
                                DATA2D_set(data[type],tr,                     n_data[type],3);
                                DATA2D_set(data[type],(double)trans_a,        n_data[type],4);
                                DATA2D_set(data[type],(double)trans_b,        n_data[type],5);
                                DATA2D_set(data[type],(double)trans_c,        n_data[type],6);
                                DATA2D_set(data[type],(double)(posn*num_type),n_data[type],7);
                                DATA2D_set(data[type],(double)type,           n_data[type],8);
                                n_data[type]++;
                            }
                        }



                    }
                }
            }
            
            for(int dt=0;dt<num_type;dt++)
            {
                nets[at*num_type + dt]=DATA2D_init(n_data[dt] ,9);
                nets[at*num_type + dt]->data_group=num_all*num_type;
                DATA2D_set(data[dt],n_data[dt]-1,0,3);

                for(int i=0;i< n_data[dt];i++)
                {
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,0),i,0);
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,1),i,1);
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,2),i,2);
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,3),i,3);
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,4),i,4);
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,5),i,5);
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,6),i,6);
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,7),i,7);
                    DATA2D_set(nets[at*num_type+ dt],DATA2D_get(data[dt],i,8),i,8);
                }  
            }      
            
        }

        DATAs_free(table);
        DATA_free(n_table);
        DATA_free(atom_pos_cell);
        DATAs_free(data);
        free(n_data);
    }
    
    return nets;
}

struct DATA ** build_coordinations_table_byRadius(struct POSCAR * input,int atom_type,double Rcutoff)
{
    int MAX_atoms=100000;
    perCrystal_full(input,0,0,0);
    POSCAR_TranCoorType(input,1);
    
    double radius=0;
    /////////////////////////////////////////////////////////////////////////////// target atoms
    int pre_num=0;
    int end_num=0;
    if(atom_type >=0)
    {
        
        for(int i=0;i<atom_type;i++)
        {
            pre_num+= input->num_atoms[i];
        }
        end_num=pre_num+(input->num_atoms[atom_type]);
    }
    else
    {
        pre_num=0;
        end_num=input->num_all;
    }
    ////////////////////////////////////////////////////////////////////////////cal max radius
    radius=Rcutoff;
    //////////////////////////////////////////////////////////////////////////////init table
    struct DATA ** nets=(struct DATA **)calloc(end_num-pre_num,sizeof(struct DATA *));
  

    struct DATA *  data=DATA2D_init(MAX_atoms,8);
    /////////////////////////////////////////////////////////////////////////////cal MAX A, B and C
    int A,B,C;
    for(int atom=pre_num;atom<end_num;atom++)
    {

        double desx,desy,desz;
        POSCAR_getPOS_order(input,atom,&desx,&desy,&desz);
        int number_neighbor=1;
        DATA2D_set(data,desx,0,0);
        DATA2D_set(data,desy,0,1);
        DATA2D_set(data,desz,0,2);
        for(int i=pre_num;i<end_num;i++)
        {
            double mtx,mty,mtz,tr;
            POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
            tr=sqrt(      (mtx-desx)*(mtx-desx)   + (mty-desy)*(mty-desy)   + (mtz-desz)*(mtz-desz)     );
            if(tr-radius < 0 && tr-0.05 >0)
            {
                DATA2D_set(data,mtx,number_neighbor,0);
                DATA2D_set(data,mty,number_neighbor,1);
                DATA2D_set(data,mtz,number_neighbor,2);
                DATA2D_set(data,tr,number_neighbor,3);
                DATA2D_set(data,(double)A,number_neighbor,4);
                DATA2D_set(data,(double)B,number_neighbor,5);
                DATA2D_set(data,(double)C,number_neighbor,6);
                DATA2D_set(data,(double)(i-pre_num),number_neighbor,7);
                number_neighbor++;
            }
        }
        int Nn=1;
        while(Nn>0)
        {
            int num_single=0;
            ////////////////////////////////////////////////////////////////////////////// cal neighbor atoms
            /////for up C=n
            C=Nn;
            for(A=-Nn;A<=Nn;A++) 
            {
                for(B=-Nn;B<=Nn;B++) 
                {
                        for(int i=pre_num;i<end_num;i++)
                        {
                            double tx,ty,tz,mtx,mty,mtz,tr;
                            POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                            perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                            tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                            if(tr-radius < 0 && tr-0.05 >0)
                            {
                                DATA2D_set(data,tx,number_neighbor,0);
                                DATA2D_set(data,ty,number_neighbor,1);
                                DATA2D_set(data,tz,number_neighbor,2);
                                DATA2D_set(data,tr,number_neighbor,3);
                                DATA2D_set(data,(double)A,number_neighbor,4);
                                DATA2D_set(data,(double)B,number_neighbor,5);
                                DATA2D_set(data,(double)C,number_neighbor,6);
                                DATA2D_set(data,(double)(i-pre_num),number_neighbor,7);
                                number_neighbor++;
                                num_single++;
                            }
                        }
                }
            }
            /////for down C=-n
            C=-Nn;
            for(A=-Nn;A<=Nn;A++) 
            {
                for(B=-Nn;B<=Nn;B++) 
                {
                        for(int i=pre_num;i<end_num;i++)
                        {
                            double tx,ty,tz,mtx,mty,mtz,tr;
                            POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                            perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                            tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                            if(tr-radius < 0 && tr-0.05 >0)
                            {
                                DATA2D_set(data,tx,number_neighbor,0);
                                DATA2D_set(data,ty,number_neighbor,1);
                                DATA2D_set(data,tz,number_neighbor,2);
                                DATA2D_set(data,tr,number_neighbor,3);
                                DATA2D_set(data,(double)A,number_neighbor,4);
                                DATA2D_set(data,(double)B,number_neighbor,5);
                                DATA2D_set(data,(double)C,number_neighbor,6);
                                DATA2D_set(data,(double)(i-pre_num),number_neighbor,7);
                                number_neighbor++;
                                num_single++;
                            }
                        }
                }
            }
            /////for left B=-n
            B=-Nn;
            for(A=-Nn;A<=Nn;A++) 
            {
                for(C=-(Nn-1);C<=(Nn-1);C++) 
                {
                        for(int i=pre_num;i<end_num;i++)
                        {
                            double tx,ty,tz,mtx,mty,mtz,tr;
                            POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                            perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                            tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                            if(tr-radius < 0 && tr-0.05 >0)
                            {
                                DATA2D_set(data,tx,number_neighbor,0);
                                DATA2D_set(data,ty,number_neighbor,1);
                                DATA2D_set(data,tz,number_neighbor,2);
                                DATA2D_set(data,tr,number_neighbor,3);
                                DATA2D_set(data,(double)A,number_neighbor,4);
                                DATA2D_set(data,(double)B,number_neighbor,5);
                                DATA2D_set(data,(double)C,number_neighbor,6);
                                DATA2D_set(data,(double)(i-pre_num),number_neighbor,7);
                                number_neighbor++;
                                num_single++;
                            }
                        }
                }
            }
            /////for right B=-n
            B=Nn;
            for(A=-Nn;A<=Nn;A++) 
            {
                for(C=-(Nn-1);C<=(Nn-1);C++) 
                {
                        for(int i=pre_num;i<end_num;i++)
                        {
                            double tx,ty,tz,mtx,mty,mtz,tr;
                            POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                            perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                            tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                            if(tr-radius < 0 && tr-0.05 >0)
                            {
                                DATA2D_set(data,tx,number_neighbor,0);
                                DATA2D_set(data,ty,number_neighbor,1);
                                DATA2D_set(data,tz,number_neighbor,2);
                                DATA2D_set(data,tr,number_neighbor,3);
                                DATA2D_set(data,(double)A,number_neighbor,4);
                                DATA2D_set(data,(double)B,number_neighbor,5);
                                DATA2D_set(data,(double)C,number_neighbor,6);
                                DATA2D_set(data,(double)(i-pre_num),number_neighbor,7);
                                number_neighbor++;
                                num_single++;
                            }
                        }
                }
            }      
            /////for front A=n
            A=Nn;
            for(B=-(Nn-1);B<=(Nn-1);B++) 
            {
                for(C=-(Nn-1);C<=(Nn-1);C++) 
                {
                        for(int i=pre_num;i<end_num;i++)
                        {
                            double tx,ty,tz,mtx,mty,mtz,tr;
                            POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                            perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                            tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                            if(tr-radius < 0 && tr-0.05 >0)
                            {
                                DATA2D_set(data,tx,number_neighbor,0);
                                DATA2D_set(data,ty,number_neighbor,1);
                                DATA2D_set(data,tz,number_neighbor,2);
                                DATA2D_set(data,tr,number_neighbor,3);
                                DATA2D_set(data,(double)A,number_neighbor,4);
                                DATA2D_set(data,(double)B,number_neighbor,5);
                                DATA2D_set(data,(double)C,number_neighbor,6);
                                DATA2D_set(data,(double)(i-pre_num),number_neighbor,7);
                                number_neighbor++;
                                num_single++;
                            }
                        }
                }
            }     
            /////for behand A=n
            A=-Nn;
            for(B=-(Nn-1);B<=(Nn-1);B++) 
            {
                for(C=-(Nn-1);C<=(Nn-1);C++) 
                {
                        for(int i=pre_num;i<end_num;i++)
                        {
                            double tx,ty,tz,mtx,mty,mtz,tr;
                            POSCAR_getPOS_order(input,i,&mtx,&mty,&mtz);
                            perCrystal(input,A,B,C,&tx,&ty,&tz,mtx,mty,mtz);
                            tr=sqrt(      (tx-desx)*(tx-desx)   + (ty-desy)*(ty-desy)   + (tz-desz)*(tz-desz)     );

                            if(tr-radius < 0 && tr-0.05 >0)
                            {
                                DATA2D_set(data,tx,number_neighbor,0);
                                DATA2D_set(data,ty,number_neighbor,1);
                                DATA2D_set(data,tz,number_neighbor,2);
                                DATA2D_set(data,tr,number_neighbor,3);
                                DATA2D_set(data,(double)A,number_neighbor,4);
                                DATA2D_set(data,(double)B,number_neighbor,5);
                                DATA2D_set(data,(double)C,number_neighbor,6);
                                DATA2D_set(data,(double)(i-pre_num),number_neighbor,7);
                                number_neighbor++;
                                num_single++;
                            }
                        }
                }
            }   
            
            ////////////////////////////////////////////////////////////////////////////// reset atom pos
            if(num_single==0) Nn=-1000;  
            Nn++;
        }
                    
        if(number_neighbor > 1)
        {
            nets[atom-pre_num]=DATA2D_init(number_neighbor ,8);
            nets[atom-pre_num]->data_group=end_num-pre_num;
            DATA2D_set(data,number_neighbor-1,0,3);
            for(int i=0;i< number_neighbor;i++)
            {
                DATA2D_set(nets[atom-pre_num],DATA2D_get(data,i,0),i,0);
                DATA2D_set(nets[atom-pre_num],DATA2D_get(data,i,1),i,1);
                DATA2D_set(nets[atom-pre_num],DATA2D_get(data,i,2),i,2);
                DATA2D_set(nets[atom-pre_num],DATA2D_get(data,i,3),i,3);
                DATA2D_set(nets[atom-pre_num],DATA2D_get(data,i,4),i,4);
                DATA2D_set(nets[atom-pre_num],DATA2D_get(data,i,5),i,5);
                DATA2D_set(nets[atom-pre_num],DATA2D_get(data,i,6),i,6);
                DATA2D_set(nets[atom-pre_num],DATA2D_get(data,i,7),i,7);
            }
        }
        else
        {
            nets[atom-pre_num]=DATA2D_init(1 ,8);
            DATA2D_set(nets[atom-pre_num],DATA2D_get(data,0,0),0,0);
            DATA2D_set(nets[atom-pre_num],DATA2D_get(data,0,1),0,1);
            DATA2D_set(nets[atom-pre_num],DATA2D_get(data,0,2),0,2);
            DATA2D_set(nets[atom-pre_num],DATA2D_get(data,0,3),0,3);
            DATA2D_set(nets[atom-pre_num],DATA2D_get(data,0,4),0,4);
            DATA2D_set(nets[atom-pre_num],DATA2D_get(data,0,5),0,5);
            DATA2D_set(nets[atom-pre_num],DATA2D_get(data,0,6),0,6);
            DATA2D_set(nets[atom-pre_num],DATA2D_get(data,0,7),0,7);
        }
    }
     ///////////////////////////////////////////////////free memery;

    DATA_free(data);
    return nets;
}

struct  DATA * RDF(struct POSCAR * poscar,double max_radius,double prec,int sample_atoms)
{

    struct DATA ** nets=build_coordinations_table_byRadius(poscar,-1,max_radius);

    ////////////////////////cal atom number
    int skip=poscar->num_all / sample_atoms;
    if(skip <= 1)
    {
        skip=1;
        sample_atoms=poscar->num_all;
    }
    else
    {
        sample_atoms=(poscar->num_all - 1)/skip + 1;
    }
    /////////////////////// number bins
    int number_bins=0;
    double dis_bins=0;
    if(prec > 0)
    {
        dis_bins=prec;
        number_bins = max_radius/dis_bins;
    }
    else if(prec <0)
    {
        number_bins=-prec;
        dis_bins=max_radius/number_bins;
    }
    struct DATA * out=DATA_init(2, number_bins ,2);
    ////////////////////////
    for(int at=0;at<poscar->num_all;at+=skip)
    {
        for(int i=1;i<nets[at]->dimen_data[0];i++)
        {
            double r=DATA2D_get(nets[at],i,3);
            int set_bins=(int)(r/dis_bins);
            if(set_bins>=number_bins) set_bins=number_bins-1;
            DATA_set(out,  DATA_get(out,set_bins,1) + 1  , set_bins,1 );            
        }
    }
    double density=1.0*(poscar->num_all)/POSCAR_getVolume(poscar->Lattice);
    //printf("%lf     %lf\n",density,  1.0*(nets[0]->dimen_data[0])/(  4.0/3.0*M_PI*pow(max_radius,3)        )          );
    for(int br=0;br<out->dimen_data[0];br++)
    {
        double para=1.0/(  sample_atoms * ( 4.0/3.0*M_PI*pow(br*dis_bins + dis_bins,3) - 4.0/3.0*M_PI*pow(br*dis_bins ,3)    )   );

        para=para/density;

        DATA2D_set(out,DATA2D_get(out,br,1)*para,br,1);
        DATA2D_set(out,br*dis_bins,br,0);
    }
    DATAs_free(nets);
    return out;

}

struct DATA * find_peak(struct DATA * input,int type,double sigma,int find_type)
{
    struct DATA * resu=NULL;
    if(input != NULL)
    {
        struct DATA *afterdata=NULL;
        if(type ==1)
        {
            afterdata=input;
        }
        else if(type ==2)
        {
            afterdata=smear_Gaussian(input,sigma);
        } 
        //////////////////////////////////////////
        int checks,peak_type;
        if(find_type >=0) //peak
        {
             peak_type=1;
             checks=-1;
        }     
        else if(find_type < 0)
        {
            peak_type=-1;
            checks=1;
        }
        struct DATA * outputs=DATA_init(2, afterdata->dimen_data[0],2);


        int number_peaks=0;
        for(int i=0;i<( afterdata->dimen_data[0])-1;i++)
        {
            double x,xp;
            x=DATA_get( afterdata,i,1);
            xp=DATA_get( afterdata,i+1,1);
            if(  (xp-x)*checks  < 0    ) //diff
            {
                if(  (xp-x)*peak_type  < 0  ) //find
                {
                        DATA_set(outputs,DATA_get(afterdata,i,0),number_peaks,0);
                        DATA_set(outputs,x,number_peaks,1);
                        number_peaks++;
                        checks=-peak_type;
                }
                else
                {
                        checks=peak_type;
                }
            }
        }
        if(number_peaks > 1)
        {
            resu=DATA_init(2,number_peaks,2);
            for(int i=0;i<number_peaks;i++)
            {
                DATA_set(resu,DATA_get(outputs,i,0),i,0);
                DATA_set(resu,DATA_get(outputs,i,1),i,1);
            }
            
        }
        else
        {
            printf("Err! functions find_peak err! No peaks find!\n");
        }

        DATA_free(outputs);
        if(type == 2)
        {
            DATA_free(afterdata);
        }
        //////////////////////////////////////////

    }
    else
    {
        printf("Err! function find_peak err!  DATA is empty!\n");
    }
    return resu;
}
void statistic(struct DATA * input, struct DATA ** output ,double prec)
{
    if(input->dimen == 1 )
    {
        double max,min;
        max=-99999999999999;
        min=999999999999999;
        for(int i=0;i<(input->dimen_data[0]);i++)
        {
            double tmps=DATA_get(input,i);
            if(   tmps- min < 0) min=  tmps;
            if(  tmps - max > 0) max=  tmps;
        }
        //////////////////////////////////
        if (prec<0)
        {
            prec = (max - min)/(- (int)prec);
        }
        int numbers = (int)( (max-min)/prec  );
        numbers += 1;
        struct DATA * sta_res=DATA_init(2,numbers,2);
        
        //////////////////////////////////
        for(int i=0;i<(input->dimen_data[0]);i++)
        {
            double tmps=DATA_get(input,i);
            int pos= (int)floor((tmps   - min )/prec );
            DATA_set(sta_res, DATA_get(sta_res, pos ,1) + 1 , pos ,1);
        }
        for(int i=0;i<numbers;i++)
        {
            DATA_set(sta_res, min + 1.0*i*prec , i , 0   );
        }
        ///////////////////////////////////

        *output= sta_res;
    }
    else
    {
        printf("function Err: statistic err: the dimension is not one!\n");
    }
}
double PolynomialFit_LSM(struct DATA * input,struct DATA ** para,int n)
{
    struct DATA * R_par=DATA1D_init(n);
    struct DATA * L_par=DATA1D_init( 2*n - 1);

    //////////////////////////////////////////////// calculate the vector and matrix's parmeters
    for(int i=0;i<(input->dimen_data[0]);i++)
    {
        /////////R_parameters
        for(int m=0;m<n;m++)
        {
            DATA_set(  R_par, DATA1D_get(R_par,m) +  DATA2D_get(input,i,1)*  pow(DATA2D_get(input,i,0), m ), m       );
        }
        /////////L_parmeters
        for(int m=0;m< (2*n-1);m++)
        {
            DATA1D_set(  L_par, DATA1D_get(L_par,m) +  pow(DATA2D_get(input,i,0), m ),    m       );
        }
    }
    //////////////////////////////////////////////// set the matrix and vector: A*X = B
    struct DATA * B=R_par;
    struct DATA * A=DATA_init(2,n,n);
    for(int row =0 ;row<n;row++)
    {
        for(int col=0;col<n;col++)
        {
            DATA2D_set(A,           DATA_get(L_par,row + col)          ,row,col);
        }
    }
    ////////////////////////////////////////////////mkl solve the liner functions
    int * ipivs=(int *)calloc(n,sizeof(int));
    int num=n,nrhs=1,lda=n,ldb=n,info;
    dgesv(&num,&nrhs,A->data_data,&lda,ipivs,B->data_data,&ldb,&info);

    double chanca=0;
    if(info != 0)
    {
        chanca=-1;
        printf("Err!  function PolynomialFit_LSM err! dgesv solve err!\n");
    } 
    else
    {
        for(int i=0;i< (input->dimen_data[0]) ;i++   )
        {
            double sumsss=0;
            for(int j=0;j<n;j++)
            {
                sumsss+= DATA1D_get(B,j)  *  pow(DATA2D_get(input,i,0),j);
            }
            chanca+=pow(  DATA2D_get(input,i,1) - sumsss  ,2);
        }
        *para = B;
    }
    DATA_free(A);
    DATA_free(L_par);
    free(ipivs);
    return chanca;
    ////////////////////////////////////////////////


}


int solidfication_check_stable(struct DATA * temperature, struct DATA * volume, struct DATA * pressure,double  precs)
{

    int unstable=1;
    struct DATA * para=NULL;
    if(temperature != NULL)
    {
        struct DATA * ntemperature=DATA2D_init(temperature->dimen_data[0]/2,2);
        int start=temperature->dimen_data[0]/2 - 1;
        for(int i=0;i<ntemperature->dimen_data[0];i++)
        {
            DATA2D_set(ntemperature,DATA2D_get(temperature,i+start,0),i,0);
            DATA2D_set(ntemperature,DATA2D_get(temperature,i+start,1),i,1);
        }

        PolynomialFit_LSM(ntemperature,&para,2);
        double sts=fabs(para->data_data[1]) ;

        if( sts - precs > 0) unstable=0;
        DATA_free(para);
        DATA_free(ntemperature);

    }

    if(volume != NULL)
    {
        struct DATA * nvolume=DATA2D_init(volume->dimen_data[0]/2,2);
        int start=volume->dimen_data[0]/2 - 1;
        for(int i=0;i<nvolume->dimen_data[0];i++)
        {
            DATA2D_set(nvolume,DATA2D_get(volume,i+start,0),i,0);
            DATA2D_set(nvolume,DATA2D_get(volume,i+start,1),i,1);
        }


        PolynomialFit_LSM(nvolume,&para,2);
        double sts=fabs(para->data_data[1]) ;

        if( sts - precs > 0) unstable=0;
        DATA_free(para);
        DATA_free(nvolume);

    }
    if(pressure != NULL)
    {
        struct DATA * npressure=DATA2D_init(pressure->dimen_data[0]/2,2);
        int start=pressure->dimen_data[0]/2 - 1;
        for(int i=0;i<npressure->dimen_data[0];i++)
        {
            DATA2D_set(npressure,DATA2D_get(pressure,i+start,0),i,0);
            DATA2D_set(npressure,DATA2D_get(pressure,i+start,1),i,1);
        }

        PolynomialFit_LSM(npressure,&para,2);
        double sts=fabs(para->data_data[1]) ;

        if( sts - precs > 0) unstable=0;
        DATA_free(para);
        DATA_free(npressure);

    }

   
    return unstable;

}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// sjk check
void Qlmi_Nor(double L,double M,struct DATA * neighbor, double *real,double *imag)
{
    double sum_qlmr=0;
    double sum_qlmi=0;
    struct DATA * final_atoms=neighbor;

    double desx,desy,desz;
    desx=DATA2D_get(final_atoms,0,0);
    desy=DATA2D_get(final_atoms,0,1);
    desz=DATA2D_get(final_atoms,0,2);
    

    for(int j=1;j<final_atoms->dimen_data[0];j++)
    {
        double rad,theta,phi;
        double ymlr,ymli;
        double rx,ry,rz;
        rx=DATA2D_get(final_atoms,j,0);
        ry=DATA2D_get(final_atoms,j,1);
        rz=DATA2D_get(final_atoms,j,2);

        cartesian2spherical(desx-rx, desy-ry,desz-rz,&rad,&theta,&phi);
        Sph_Har(L,M,phi,theta,&ymlr,&ymli);

        sum_qlmr+=ymlr;
        sum_qlmi+=ymli;
    }
    double numbers=1.0*(final_atoms->dimen_data[0])-1;
    sum_qlmr=sum_qlmr/numbers;
    sum_qlmi=sum_qlmi/numbers;
    ///////////////////////////////
    double sumss=0;
    for(int m=-L;m<=L;m++)
    {
        double suL=0;
        double suR=0;
        for(int j=1;j<final_atoms->dimen_data[0];j++)
        {
            double rad,theta,phi;
            double ymlr,ymli;
            double rx,ry,rz;
            rx=DATA2D_get(final_atoms,j,0);
            ry=DATA2D_get(final_atoms,j,1);
            rz=DATA2D_get(final_atoms,j,2);

            cartesian2spherical(desx-rx, desy-ry,desz-rz,&rad,&theta,&phi);
            Sph_Har(L,m,phi,theta,&ymlr,&ymli);

            suR+=ymlr;
            suL+=ymli;
        }  
        suL=suL/numbers;
        suR=suR/numbers;  

        sumss=sumss+ ( suL*suL + suR*suR   );
    }

    sumss=sqrt(sumss);

    *real=sum_qlmr/sumss;
    *imag=sum_qlmi/sumss;
}
struct DATA * Qlm_full_Nor(struct DATA ** nets,int L)
{
    int numbers=nets[0]->data_group;
    struct DATA * qlms=DATA3D_init(numbers,2*L+1,2);
    
    for(int i=0;i<numbers;i++)
    {
        for(int m=-L;m<=L;m++)
        {
            double re,img;
            Qlmi_Nor(L,m,nets[i],&re,&img);
            DATA3D_set(qlms,re,i,m+L,0);
            DATA3D_set(qlms,img,i,m+L,1);

        }
    }
    return qlms;
}
double Sjk_Nor(struct DATA * qlm,int sel_atom1,int sel_atom2)
{
    ///////////////////////////
    double real=0;
    double jr,ji;
    double kr,ki;
    int L=qlm->dimen_data[1];
    for(int m=0;m<L;m++)
    {
        jr=DATA3D_get(qlm,sel_atom1,m,0);
        ji=DATA3D_get(qlm,sel_atom1,m,1);

        kr=DATA3D_get(qlm,sel_atom2,m,0);
        ki=DATA3D_get(qlm,sel_atom2,m,1);

        real+=   jr*kr+ji*ki;

    }

    return real;
}
struct DATA * Sjk_Avg_Nor_full(struct DATA * qml,struct DATA ** nets)
{
    struct DATA * out=DATA1D_init(qml->dimen_data[0]);
    for(int at=0;at<(qml->dimen_data[0]);at++)
    {
        double sum=0;
        for(int i=1;i<   (nets[at]->dimen_data[0]);i++   )
        {
            int sel_atom2=(int)DATA2D_get(nets[at],i,7);
            sum+=Sjk_Nor(qml,at,sel_atom2);
        }
        sum=sum/(  (nets[at]->dimen_data[0]) -1    );
        out->data_data[at]=sum;
    }
    return out;   
}
struct DATA * solid_liquid_bySjk_atom_from_best_table(struct POSCAR * poscar,struct DATA * learn_res,double neigh_ratio,double *solidfication)
{
    struct DATA ** nets_par=NULL;
    struct DATA * qml_par=NULL;
    struct DATA * sjk_par=NULL;


    int ground_L=(int)DATA2D_get(learn_res,0,1);
    double ground_prec=DATA2D_get(learn_res,0,2);
    struct DATA ** nets_all=build_coordinations_table_bySANN(poscar,-1);
    struct DATA * qml_all=Qlm_full_Nor(nets_all,ground_L);
    struct DATA * sjk_all=Sjk_Avg_Nor_full(qml_all,nets_all);


    struct DATA * ret_table=DATA1D_init(poscar->num_all);
    int solid_number=0;
    int liquid_number=0;

    double ground_auc=DATA2D_get(learn_res,0,3);
    
    for(int ty=0;ty<poscar->num_type;ty++)
    {
        double ty_auc=DATA2D_get(learn_res,ty+1,3);
        int pre_num,end_num;
        POSCAR_getOrder(poscar,ty,&pre_num,&end_num);
        if(ty_auc - ground_auc >0)
        {
            double par_L=(int)DATA2D_get(learn_res,ty+1,1);
            double par_prec=DATA2D_get(learn_res,ty+1,2);
            nets_par=build_coordinations_table_bySANN(poscar,ty);
            qml_par=Qlm_full_Nor(nets_par,par_L);
            
            sjk_par=Sjk_Avg_Nor_full(qml_par,nets_par);
            for(int i=0;i<qml_par->dimen_data[0];i++)
            {
                if(sjk_par->data_data[i] - par_prec >=0)
                {
                    solid_number++;
                    ret_table->data_data[i+pre_num] = 1.0;
                }
                else
                {
                    liquid_number++;
                    ret_table->data_data[i+pre_num] = 0;
                }
            }
            DATA_free(sjk_par);
            DATA_free(qml_par);
            DATAs_free(nets_par);
        }
        else 
        {
            for(int i=pre_num;i<end_num;i++)
            {
                if(sjk_all->data_data[i] - ground_prec >=0)
                {
                    solid_number++;
                    ret_table->data_data[i]=1.0;
                }
                else
                {
                    liquid_number++;
                    ret_table->data_data[i]=0;
                }
            }
        }
    }



    for(int at=0;at<poscar->num_all;at++)
    {
        if(ret_table->data_data[at] - 0.5 <0)//liquid
        {
            int solid_num=0;
            for(int i=1;i<nets_all[at]->dimen_data[0];i++)
            {
                int nei=DATA2D_get(nets_all[at],i,7);
                if(ret_table->data_data[nei] - 0.5 >0)
                {
                    solid_num++;
                }
            }
            double ratio=1.0*solid_num/(  nets_all[at]->dimen_data[0] -1     );
            if(ratio - neigh_ratio > 0)
            {
                ret_table->data_data[at]=1.0;
            }
        }
        else                                  //solid
        {
            int liquid_num=0;
            for(int i=1;i<nets_all[at]->dimen_data[0];i++)
            {
                int nei=DATA2D_get(nets_all[at],i,7);
                if(ret_table->data_data[nei] - 0.5 <0)
                {
                    liquid_num++;
                }
            }
            double ratio=1.0*liquid_num/(  nets_all[at]->dimen_data[0] -1     );
            if(ratio - neigh_ratio > 0)
            {
                ret_table->data_data[at]=0;
            }            
        }
    }


    DATA_free(qml_all);
    DATA_free(sjk_all);
    DATAs_free(nets_all);


    *solidfication=1.0*solid_number/(solid_number+liquid_number);
    return ret_table;
}
struct DATA * max_AUC_SJK_learn_from_POSCAR(struct POSCAR * poscar_solid,struct POSCAR * poscar_liquid,int solid_sample_number,int liquid_sample_number,int *per_L,int * per_atom_type,double * sjk_perc,struct DATA ** out_result)
{
    struct DATA ** nets_solid=NULL;
    struct DATA ** nets_liquid=NULL;
    struct DATA * qml_solid=NULL;
    struct DATA * qml_liquid=NULL;
    struct DATA * data=NULL;
    struct DATA * state=NULL;

    struct DATA * outs=DATA2D_init(poscar_solid->num_type + 1,4);


    double max_auc=-99999;
    struct DATA * max_auc_cure=NULL;
    double m_ty=-1;
    double m_L=0;
    double m_prec=-1;

    struct DATA * tmp_auc=NULL;

    for(int ty=-1;ty<poscar_solid->num_type;ty++)
    {

        nets_solid=build_coordinations_table_bySANN(poscar_solid,ty);
        nets_liquid=build_coordinations_table_bySANN(poscar_liquid,ty);

        int solid_numbers=0;
        int liquid_numbers=0;
        if(ty==-1)
        {
            solid_numbers=poscar_solid->num_all;
            liquid_numbers=poscar_liquid->num_all;
        }
        else
        {   
            solid_numbers=poscar_solid->num_atoms[ty];
            liquid_numbers=poscar_liquid->num_atoms[ty];
        }

        int s_s_n=solid_sample_number;
        int l_s_n=liquid_sample_number;
        int solid_skip=solid_numbers/s_s_n;
        int liquid_skip=liquid_numbers/l_s_n;

        if(solid_skip <1) solid_skip=1;
        if(liquid_skip <1) liquid_skip=1;
        int solid_cal_numbers=(solid_numbers-1)/solid_skip +1;

        int liquid_cal_numbner=(liquid_numbers -1)/liquid_skip +1;

        int total_data_number=solid_cal_numbers+liquid_cal_numbner;

        data=DATA2D_init(total_data_number,2);
        state=DATA1D_init(total_data_number);


        double tmp_max_auc=-99999;
        int tmp_max_L;
        double tmp_max_prec;
        for(int L=4;L<=8;L+=2)
        {
            qml_solid=Qlm_full_Nor(nets_solid,L);
            qml_liquid=Qlm_full_Nor(nets_liquid,L);
            ///////////////////////////////////////////////////////////data_preper
            double maxs_sij=-9999;
            double mins_sij=9999;
            double maxl_sij=-9999;
            double minl_sij=9999;
            int set_numbers=0;
            for(int i=0;i<solid_numbers;i+=solid_skip)
            {
                double sum=0;
                for(int ni=1;ni<   (nets_solid[i]->dimen_data[0]);ni++   )
                {
                    int sel_atom2=(int)DATA2D_get(nets_solid[i],ni,7);
                    sum+=Sjk_Nor(qml_solid,i,sel_atom2);
                }
                sum=sum/(  (nets_solid[i]->dimen_data[0]) -1    );
                
                
                DATA2D_set(data,sum,set_numbers,1);
                state->data_data[set_numbers]=1.0;
                if(sum - maxs_sij >0) maxs_sij=sum;
                if(sum - mins_sij <0) mins_sij=sum;
                set_numbers++;
            }
     
            for(int i=0;i<liquid_numbers;i+=liquid_skip)
            {


                double sum=0;
                for(int ni=1;ni<   (nets_liquid[i]->dimen_data[0]);ni++   )
                {
                    int sel_atom2=(int)DATA2D_get(nets_liquid[i],ni,7);
                    sum+=Sjk_Nor(qml_liquid,i,sel_atom2);
                }
                sum=sum/(  (nets_liquid[i]->dimen_data[0]) -1    );


                DATA2D_set(data,sum,set_numbers,1);
                state->data_data[set_numbers]=0;
                if(sum - maxl_sij >0) maxl_sij=sum;
                if(sum - minl_sij <0) minl_sij=sum;
                set_numbers++;
            }
     
            ///////////////////////////////////////////////////////////
            double aucs=0;
            tmp_auc=AUC_ROC_cure(data,state,0,1,20,&aucs);


            double per_prec=-1;

            if( aucs - max_auc >0   )
            {
                max_auc=aucs;
                if(max_auc_cure != NULL)
                {
                    DATA_free(max_auc_cure);
                }
                max_auc_cure=tmp_auc;
                m_L=L;
                m_ty=ty;
                m_prec=per_prec;
            }
            else
            {
                DATA_free(tmp_auc);
            }

            if(aucs -tmp_max_auc >0)
            {
                tmp_max_auc=aucs;
                tmp_max_L=L;
                tmp_max_prec=  0.5*(0.5*(maxs_sij+mins_sij)  +  0.5*(maxl_sij+minl_sij) )  ;
            }

            DATA_free(qml_solid);
            DATA_free(qml_liquid);

        }
        
        DATA2D_set(outs,1.0*ty,       ty+1,0);
        DATA2D_set(outs,1.0*tmp_max_L,ty+1,1);
        DATA2D_set(outs,tmp_max_prec, ty+1,2);
        DATA2D_set(outs,tmp_max_auc, ty+1,3);


        DATAs_free(nets_solid);
        DATAs_free(nets_liquid);
        DATA_free(data);
        DATA_free(state);

    }


    *per_L=m_L;
    *per_atom_type=m_ty;
    *sjk_perc=m_prec;
    *out_result=outs;
    return max_auc_cure;

}
struct DATA * solid_liquid_bySjk(struct DATA * qml,struct DATA ** nets,double compare_sij,double solid_ratio,double neigh_ratio,double * solidfication)
{
    int solid=0;
    int liquid=0;
    struct DATA * out=DATA1D_init(qml->dimen_data[0]);
    for(int at=0;at<(qml->dimen_data[0]);at++)
    {
        int sim=0;
        for(int i=1;i<   (nets[at]->dimen_data[0]);i++   )
        {
            int sel_atom2=(int)DATA_get(nets[at],i,7);
            double sijs=Sjk_Nor(qml,at,sel_atom2);
            if(sijs - compare_sij>0)
            {
                sim++;
            }
        }
        double ratio=1.0*sim/((nets[at]->dimen_data[0])-1);
        if(ratio - solid_ratio>0)
        {
            DATA1D_set(out,1,at);
            solid++;
        }
        else
        {
            DATA1D_set(out,0,at);
            liquid++;
        }
    }


    for(int at=0;at<(qml->dimen_data[0]);at++)
    {
        if(out->data_data[at] - 0.5 <0)//liquid
        {
            int solid_num=0;
            for(int i=1;i<nets[at]->dimen_data[0];i++)
            {
                int nei=DATA2D_get(nets[at],i,7);
                if(out->data_data[nei] - 0.5 >0)
                {
                    solid_num++;
                }
            }
            double ratio=1.0*solid_num/(  nets[at]->dimen_data[0] -1     );
            if(ratio - neigh_ratio > 0)
            {
                out->data_data[at]=1.0;
            }
        }
        else                                  //solid
        {
            int liquid_num=0;
            for(int i=1;i<nets[at]->dimen_data[0];i++)
            {
                int nei=DATA2D_get(nets[at],i,7);
                if(out->data_data[nei] - 0.5 <0)
                {
                    liquid_num++;
                }
            }
            double ratio=1.0*liquid_num/(  nets[at]->dimen_data[0] -1     );
            if(ratio - neigh_ratio > 0)
            {
                out->data_data[at]=0;
            }            
        }
    }


    *solidfication=1.0*solid/(solid+liquid);
    return out;
}
struct DATA ** solid_liquid_bySjkStatistic_get_reference(struct DATA ** nets_decompose,int number_atoms,int number_atom_types, int * number_references, int L)
{
    struct DATA * tmp_qml=DATA3D_init(number_atom_types,2*L+1,2);//qml for different atom types
    
    struct DATA ** tmp_qml_total=(struct DATA **)calloc(number_atoms,sizeof(struct DATA *));
    int n_tmp_qml_total=0;

    for(int at=0;at<number_atoms;at++)
    { 
        //////////////////////////////////////for qml
        for(int ty=0;ty<number_atom_types;ty++)
        {
            for(int m=-L;m<=L;m++)
            {
                double real,img;
                Qlmi_Nor(L,m,nets_decompose[at*number_atom_types + ty],&real,&img);
                DATA3D_set(tmp_qml,real,ty,m+L,0);
                DATA3D_set(tmp_qml,img,ty,m+L,1);
            }
        }
        //////////////////////////////////////for similar
        
        if(n_tmp_qml_total==0)
        {
            tmp_qml_total[0]=tmp_qml;
            tmp_qml=DATA3D_init(number_atom_types,2*L+1,2);
            n_tmp_qml_total++;
        }
        else
        {
            int check=0;
            for(int i=0;i<n_tmp_qml_total;i++)
            {
                double sums=0;
                ////////////////
                for(int ty=0;ty<number_atom_types;ty++)
                {
                    double real=0;
                    double real_a=0;
                    double real_b=0;
                    double jr,ji;
                    double kr,ki;

                    for(int m=0;m<2*L + 1;m++)
                    {
                        jr=DATA3D_get(tmp_qml,ty,m,0);
                        ji=DATA3D_get(tmp_qml,ty,m,1);
                    
                        kr=DATA3D_get(tmp_qml_total[i],ty,m,0);
                        ki=DATA3D_get(tmp_qml_total[i],ty,m,1);
                    
                        real+=   jr*kr+ji*ki;
                        real_a += jr*jr+ji*ji;
                        real_b +=kr*kr + ki*ki;
                    
                    }
                    real/=(  pow(real_a,0.5) * pow(real_b,0.5)   )   ;
                    
                    sums+=real;
                }
                sums/=number_atom_types;
                ////////////////
                if(sums>0.95)
                {
                    check=1;
                    break;
                }
            }
            if(check == 0)
            {
                tmp_qml_total[n_tmp_qml_total]=tmp_qml;
                tmp_qml=DATA3D_init(number_atom_types,2*L+1,2);
                n_tmp_qml_total++;
            }
        }

    }


    struct DATA ** rets=(struct DATA **)calloc(n_tmp_qml_total,sizeof(struct DATA *));
    for(int i=0;i<n_tmp_qml_total;i++)
    {
        rets[i]=tmp_qml_total[i];
        rets[i]->data_group=n_tmp_qml_total;
    }

    free(tmp_qml_total);
    DATA_free(tmp_qml);
    *number_references=n_tmp_qml_total;

    return rets;
}
struct DATA *  solid_liquid_bySjkStatistic_from_reference(struct DATA ** nets_decompose,int number_atoms, int number_atoms_type ,struct DATA ** reference, int number_refer, double sjk_cutoff, double neigh_cutoff, double * solidfication)
{
    struct DATA * ret=DATA1D_init(number_atoms);

    int L=reference[0]->dimen_data[1];
    L=(L-1)/2;

    struct DATA * tmp_qml=DATA3D_init(number_atoms_type,2*L+1,2);
    int solid=0;
    for(int at=0;at<number_atoms;at++)
    { 
        //////////////////////////////////////for qml
        for(int ty=0;ty<number_atoms_type;ty++)
        {
            for(int m=-L;m<=L;m++)
            {
                double real,img;
                Qlmi_Nor(L,m,nets_decompose[at*number_atoms_type + ty],&real,&img);
                DATA3D_set(tmp_qml,real,ty,m+L,0);
                DATA3D_set(tmp_qml,img,ty,m+L,1);
            }
        }
        //////////////////////////////////////for similar
        double max_sjk=-9999;
            for(int i=0;i<number_refer;i++)
            {
                double sums=0;
                ////////////////
                for(int ty=0;ty<number_atoms_type;ty++)
                {
                    double real=0;
                    double real_a=0;
                    double real_b=0;
                    double jr,ji;
                    double kr,ki;

                    for(int m=0;m<2*L + 1;m++)
                    {
                        jr=DATA3D_get(tmp_qml,ty,m,0);
                        ji=DATA3D_get(tmp_qml,ty,m,1);
                    
                        kr=DATA3D_get(reference[i],ty,m,0);
                        ki=DATA3D_get(reference[i],ty,m,1);
                    
                        real+=   jr*kr+ji*ki;
                        real_a += jr*jr+ji*ji;
                        real_b +=kr*kr + ki*ki;
                    
                    }
                    real/=(  pow(real_a,0.5) * pow(real_b,0.5)   )   ;
                    
                    sums+=real;
                }
                sums/=number_atoms_type;
                ////////////////
                if(max_sjk -sums <0)
                {
                    max_sjk=sums;
                }
            }
        
        if(max_sjk - sjk_cutoff > 0)
        {
            ret->data_data[at] = 1.0;
            solid++;
        }

    }
    double melt=1.0*solid/number_atoms;
    ////////////////////////////////////////////////////////////
    for(int at=0;at<number_atoms;at++)
    {
        int total=0;
        int r_total=0;
        for(int ty=0;ty<number_atoms_type;ty++)
        {
            for(int ne=1;ne<nets_decompose[at*number_atoms_type + ty]->dimen_data[0];ne++)
            {
                int n_nei=(int)DATA2D_get(nets_decompose[at*number_atoms_type + ty],ne,7);
                n_nei=n_nei/number_atoms_type;
                int n_st=ret->data_data[n_nei];

                if(   (int)(ret->data_data[at]) != n_st    )
                {
                    r_total++;
                }
                total++;
            }
        }
        if(1.0*r_total/total - neigh_cutoff > 0)
        {
            if(ret->data_data[at] < 0.5 ) ret->data_data[at]=1.0;
            else ret->data_data[at]=0;
        }
    }
    solid=0;
    for(int at=0;at<number_atoms;at++)
    {
        if(ret->data_data[at] - 0.5 >0) solid++;
    }
    melt=1.0*solid/number_atoms;
    ////////////////////////////////////////////////////////////
    DATA_free(tmp_qml);
    tmp_qml=NULL;
    *solidfication=melt;
    return ret;

}
struct DATA *  solid_liquid_bySjkStatistic_getSjk(struct DATA ** nets_decompose,int number_atoms, int number_atoms_type ,struct DATA ** reference, int number_refer)
{
    struct DATA * ret=DATA1D_init(number_atoms);

    int L=reference[0]->dimen_data[1];
    L=(L-1)/2;

    struct DATA * tmp_qml=DATA3D_init(number_atoms_type,2*L+1,2);
    for(int at=0;at<number_atoms;at++)
    { 
        //////////////////////////////////////for qml
        for(int ty=0;ty<number_atoms_type;ty++)
        {
            for(int m=-L;m<=L;m++)
            {
                double real,img;
                Qlmi_Nor(L,m,nets_decompose[at*number_atoms_type + ty],&real,&img);
                DATA3D_set(tmp_qml,real,ty,m+L,0);
                DATA3D_set(tmp_qml,img,ty,m+L,1);
            }
        }
        //////////////////////////////////////for similar
        double max_sjk=-9999;
            for(int i=0;i<number_refer;i++)
            {
                double sums=0;
                ////////////////
                for(int ty=0;ty<number_atoms_type;ty++)
                {
                    double real=0;
                    double real_a=0;
                    double real_b=0;
                    double jr,ji;
                    double kr,ki;

                    for(int m=0;m<2*L + 1;m++)
                    {
                        jr=DATA3D_get(tmp_qml,ty,m,0);
                        ji=DATA3D_get(tmp_qml,ty,m,1);
                    
                        kr=DATA3D_get(reference[i],ty,m,0);
                        ki=DATA3D_get(reference[i],ty,m,1);
                    
                        real+=   jr*kr+ji*ki;
                        real_a += jr*jr+ji*ji;
                        real_b +=kr*kr + ki*ki;
                    
                    }
                    real/=(  pow(real_a,0.5) * pow(real_b,0.5)   )   ;
                    
                    sums+=real;
                }
                sums/=number_atoms_type;
                ////////////////
                if(max_sjk -sums <0)
                {
                    max_sjk=sums;
                }
            }
        

        ret->data_data[at] = max_sjk ;
        

    }
    DATA_free(tmp_qml);
    tmp_qml=NULL;
    return ret;

}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// local density check
double local_density(struct DATA * net_reference,struct DATA * net_cal,double sigma)
{
    double sums=0;

    double posx,posy,posz;
    posx=DATA2D_get(net_cal,0,0);
    posy=DATA2D_get(net_cal,0,1);
    posz=DATA2D_get(net_cal,0,2);

    double posrx,posry,posrz;
    posrx=DATA2D_get(net_reference,0,0);
    posry=DATA2D_get(net_reference,0,1);
    posrz=DATA2D_get(net_reference,0,2);

        for(int n1=1;n1<net_cal->dimen_data[0];n1++)
        {
            double x1,y1,z1;
            x1=DATA2D_get(net_cal,n1,0);
            y1=DATA2D_get(net_cal,n1,1);
            z1=DATA2D_get(net_cal,n1,2);
            x1-=posx;
            y1-=posy;
            z1-=posz;
            for(int n2=1;n2<net_reference->dimen_data[0];n2++)
            {
                double x2,y2,z2;
                x2=DATA2D_get(net_reference,n2,0);
                y2=DATA2D_get(net_reference,n2,1);
                z2=DATA2D_get(net_reference,n2,2);
                x2-=posrx;
                y2-=posry;
                z2-=posrz;
                double radius=sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)  );
                sums+=Gaussian_functions(1.0,radius,0,4*sigma*sigma);
            }
 
        }

    return sums;
}
struct DATA ** solid_liquid_byLD_get_reference(struct DATA ** nets_decompose,int number_atoms,int number_atom_types, int * number_references, struct DATA ** normal , double sigma)
{
    struct DATA * tmp_qml=DATA1D_init(number_atom_types*number_atoms);//qml for different atom types

    int * tmp_qml_total=(int*)calloc(number_atoms,sizeof(int));
    int n_tmp_qml_total=0;

    for(int at=0;at<number_atoms*number_atom_types;at++)
    {
        double val=local_density(nets_decompose[at],nets_decompose[at],sigma);
        tmp_qml->data_data[at]=val;
    }


    for(int at=0;at<number_atoms;at++)
    { 

        //////////////////////////////////////for similar
        
        if(n_tmp_qml_total==0)
        {
            tmp_qml_total[n_tmp_qml_total] = at;
            n_tmp_qml_total++;
        }
        else
        {
            int check=0;
            for(int i=0;i<n_tmp_qml_total;i++)
            {
                double sums=0;
                ////////////////
                for(int ty=0;ty<number_atom_types;ty++)
                {

                    double real=local_density(nets_decompose[tmp_qml_total[i] + ty],nets_decompose[at*number_atom_types + ty],sigma);
                    real/=tmp_qml->data_data[tmp_qml_total[i]*number_atom_types + ty];
                    sums+=real;
                }
                sums/=number_atom_types;
                ////////////////
                if(sums>0.95)
                {
                    check=1;
                    break;
                }
            }
            if(check == 0)
            {
                tmp_qml_total[n_tmp_qml_total] = at;
                n_tmp_qml_total++;
            }
        }

    }


    struct DATA ** rets=(struct DATA **)calloc(n_tmp_qml_total*number_atom_types,sizeof(struct DATA *));
    struct DATA *  rets_normal=DATA1D_init(n_tmp_qml_total*number_atom_types);
    for(int i=0;i<n_tmp_qml_total;i++)
    {
        for(int ty=0;ty<number_atom_types;ty++)
        {
            rets[i*number_atom_types + ty] = DATA_cpy(nets_decompose[tmp_qml_total[i] + ty]);
            rets_normal->data_data[ i*number_atom_types + ty] = tmp_qml->data_data[tmp_qml_total[i]*number_atom_types + ty];
        }
        rets[i]->data_group=n_tmp_qml_total;
    }

    free(tmp_qml_total);
    DATA_free(tmp_qml);
    *number_references=n_tmp_qml_total;
    *normal=rets_normal;

    return rets;
}
struct DATA * solid_liquid_by_LD_from_reference(struct DATA **nets_decompose,int number_atoms,int number_atoms_types,struct DATA ** reference,struct DATA * reference_normal,int number_refer,double ld_cutoff,double neigh_cutoff, double sigma,double * solidfication)
{
    struct DATA * ret=DATA1D_init(number_atoms);


    int solid=0;
    for(int at=0;at<number_atoms;at++)
    { 
        //////////////////////////////////////for qml

        //////////////////////////////////////for similar
        double max_sjk=-9999;
            for(int i=0;i<number_refer;i++)
            {
                double sums=0;
                ////////////////
                for(int ty=0;ty<number_atoms_types;ty++)
                {

                    double real=local_density(reference[i*number_atoms_types + ty],nets_decompose[at*number_atoms_types + ty],sigma);
                    real/=reference_normal->data_data[i*number_atoms_types + ty]; 

                    sums+=real;
                }
                
                sums/=number_atoms_types;
                
                ////////////////
                if(max_sjk -sums <0)
                {
                    max_sjk=sums;
                }
            }
        //printf("%d %lf\n",at,max_sjk);
        
        if(max_sjk - ld_cutoff > 0)
        {
            ret->data_data[at] = 1.0;
            solid++;
        }

    }
    double melt=1.0*solid/number_atoms;
    ////////////////////////////////////////////////////////////
    for(int at=0;at<number_atoms;at++)
    {
        int total=0;
        int r_total=0;
        for(int ty=0;ty<number_atoms_types;ty++)
        {
            for(int ne=1;ne<nets_decompose[at*number_atoms_types + ty]->dimen_data[0];ne++)
            {
                int n_nei=(int)DATA2D_get(nets_decompose[at*number_atoms_types + ty],ne,7);
                n_nei=n_nei/number_atoms_types;
                int n_st=ret->data_data[n_nei];

                if(   (int)(ret->data_data[at]) != n_st    )
                {
                    r_total++;
                }
                total++;
            }
        }
        if(1.0*r_total/total - neigh_cutoff > 0)
        {
            if(ret->data_data[at] < 0.5 ) ret->data_data[at]=1.0;
            else ret->data_data[at]=0;
        }
    }
    solid=0;
    for(int at=0;at<number_atoms;at++)
    {
        if(ret->data_data[at] - 0.5 >0) solid++;
    }
    melt=1.0*solid/number_atoms;
    ////////////////////////////////////////////////////////////

    *solidfication=melt;
    return ret;

}
struct DATA * solid_liquid_by_LD_getLD(struct DATA **nets_decompose,int number_atoms,int number_atoms_types,struct DATA ** reference,struct DATA * reference_normal,int number_refer, double sigma)
{
    struct DATA * ret=DATA1D_init(number_atoms);
    for(int at=0;at<number_atoms;at++)
    { 
        //////////////////////////////////////for qml

        //////////////////////////////////////for similar
        double max_sjk=-9999;
            for(int i=0;i<number_refer;i++)
            {
                double sums=0;
                ////////////////
                for(int ty=0;ty<number_atoms_types;ty++)
                {

                    double real=local_density(reference[i*number_atoms_types + ty],nets_decompose[at*number_atoms_types + ty],sigma);
                    real/=reference_normal->data_data[i*number_atoms_types + ty]; 

                    sums+=real;
                }
                
                sums/=number_atoms_types;
                
                ////////////////
                if(max_sjk -sums <0)
                {
                    max_sjk=sums;
                }
            }
        
        ret->data_data[at] = max_sjk;

    }
    ////////////////////////////////////////////////////////////
   
    return ret;

}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// atom_check func

int solid_liquid_ml_getNEP_and_Qlmi_LD_Neigh_all(struct POSCAR * pos_perfect,struct POSCAR ** poss_solid,struct POSCAR ** poss_liquid,int max_frame_solid,int max_frame_liquid,int nep_N,int nep_L,int qlmi_L,double ld_sigma,double rcut,int max_config_solid,int max_config_liquid,double ** inputs,double ** outputs)
{

    struct DATA **refer_qlmi;
    int max_refer_qlmi;

    struct DATA ** refer_ld;
    struct DATA * refer_ld_nor;
    int max_refer_ld;


    struct DATA ** nets_per=build_coordinations_table_byRadius_typedecompose(pos_perfect,rcut);

    refer_qlmi=solid_liquid_bySjkStatistic_get_reference(nets_per,pos_perfect->num_all,pos_perfect->num_type,&max_refer_qlmi,qlmi_L);
    refer_ld=solid_liquid_byLD_get_reference(nets_per,pos_perfect->num_all,pos_perfect->num_type,&max_refer_ld,&refer_ld_nor,ld_sigma);

    /////////////////////////////////////////////////////
    int num_all_s=poss_solid[0]->num_all;
    int num_type_s=poss_solid[0]->num_type;
    int num_all_l=poss_liquid[0]->num_all;
    int num_type_l=poss_liquid[0]->num_type;


    int cal_solid_fram=max_config_solid/num_all_s + 1;
    int cal_liquid_fram=max_config_liquid/num_all_l +1;

    int skip_s=1;
    if(max_frame_solid > 2*cal_solid_fram ) skip_s=max_frame_solid/cal_solid_fram;

    int skip_l=1;
    if(max_frame_liquid > 2*cal_liquid_fram) skip_l=max_frame_liquid/cal_liquid_fram;

    int real_solid=0;
    int real_liquid=0;

    struct NEP * nep;
    struct DATA ** nets;
    int size_input=nep_N*nep_N*nep_L;
    double *input=(double *)calloc(size_input*max_config_solid + size_input*max_config_liquid,sizeof(double));
    double *output=(double *)calloc(2*max_config_solid + 2*max_config_liquid,sizeof(double));

    struct DATA * tmp_qlmi;
    struct DATA * tmp_ld;

    struct DATA * tmp_nep;
    ////////////////////////////////////////////////for solid describter

    ////////////////////////////////////
    nep=nep_init(nep_N,nep_L,rcut);
    
    for(int fra=0;fra<max_frame_solid;fra+=skip_s)
    {
        nets=build_coordinations_table_byRadius_typedecompose(poss_solid[fra],rcut);
        ////cal qlmi ld neigh and nep;
        tmp_nep=nep_generate_nepdescribtor(nets,num_type_s,num_all_s,nep,rcut,false,true,false,true);
        tmp_qlmi=solid_liquid_bySjkStatistic_getSjk(nets,num_all_s,num_type_s,refer_qlmi,max_refer_qlmi);
        tmp_ld=solid_liquid_by_LD_getLD(nets,num_all_s,num_type_s,refer_ld,refer_ld_nor,max_refer_ld,ld_sigma);

        ////keep the minimal nep , qlmi , ld and neigh
        for(int at=0;at<num_all_s;at++)
        {
                //copy nep
                for(int i=0;i<size_input;i++)
                {
                    input[real_solid * size_input + i]=DATA2D_get(tmp_nep,at,i);
                }
                //copy qlmi
                output[real_solid * 2 + 0] =tmp_qlmi->data_data[at];
                output[real_solid * 2 + 1] =tmp_ld->data_data[at];

                real_solid++;
           
            if(real_solid >=max_config_solid)
            {
                break;
            }

        }   
        //////////
        DATAs_free(nets); DATAs_free(nets_per);DATA_free(tmp_qlmi);
        DATA_free(tmp_nep);
        DATA_free(tmp_ld);
        if(real_solid >=max_config_solid)
        {
            break;
        }
    }

    ////////////////////////////////////////////////////for liquid describter

    for(int fra=0;fra<max_frame_liquid;fra+=skip_l)
    {
        nets=build_coordinations_table_byRadius_typedecompose(poss_liquid[fra],rcut);
        ////cal qlmi ld neigh and nep;
        tmp_nep=nep_generate_nepdescribtor(nets,num_type_l,num_all_l,nep,rcut,false,true,false,true);
        tmp_qlmi=solid_liquid_bySjkStatistic_getSjk(nets,num_all_l,num_type_l,refer_qlmi,max_refer_qlmi);
        tmp_ld=solid_liquid_by_LD_getLD(nets,num_all_l,num_type_l,refer_ld,refer_ld_nor,max_refer_ld,ld_sigma);

        ////keep the minimal nep , qlmi , ld and neigh
        for(int at=0;at<num_all_l;at++)
        {

                //copy nep
                for(int i=0;i<size_input;i++)
                {
                    input[(real_solid + real_liquid) * size_input + i]=DATA2D_get(tmp_nep,at,i);
                }
                //copy qlmi
                output[(real_solid+real_liquid) * 2 + 0] =tmp_qlmi->data_data[at];
                output[(real_solid+real_liquid) * 2 + 1] =tmp_ld->data_data[at];

                real_liquid++;
        
            if(real_liquid >=max_config_liquid)
            {
                break;
            }

        }   
        //////////
        DATAs_free(nets);
        DATA_free(tmp_nep);
        DATA_free(tmp_ld);
        DATA_free(tmp_qlmi);
        if(real_liquid >=max_config_liquid)
        {
            break;
        }
    }

    nep_free(nep);
    ///////////////////////////////////////////////////for min check
    if(real_solid + real_liquid < max_config_liquid + max_config_solid)
    {
        double * tmp_input=(double *)calloc(size_input*real_solid+ size_input*real_liquid,sizeof(double));
        double * tmp_output=(double *)calloc(2*real_solid + 2*real_liquid,sizeof(double));   
        for(int i=0;i<real_solid + real_liquid;i++)
        {
            for(int j=0;j<size_input;j++)
            {
                tmp_input[i*size_input + j] = input[i*size_input + j];
            }
            for(int j=0;j<3;j++)
            {
                tmp_output[i*2 + j] = output[i*2 + j];
            }
        } 
        free(input);
        free(output);
        input=tmp_input;
        output=tmp_output;    
    }
    DATAs_free(refer_qlmi);DATA_free(refer_ld_nor);DATAs_free(refer_ld);
    *inputs=input;
    *outputs=output;

    return real_solid + real_liquid;
}

struct NEURAL * solid_liquid_ml_train(double * input,double * output,int number_data,int size_input,int size_output)
{
    struct NEURAL * neu=neural_init_MLP(size_input,size_output,3,15,15,15);

    neural_set_data(neu,input,output,number_data);

    neural_set_activate_function(neu,activate_function_leaky_ReLU,G_activate_function_leaky_ReLU,0);
    neural_set_activate_function(neu,activate_function_default,G_activate_function_default,-1);

    neural_set_loss_func(neu,loss_func_MSE,G_loss_func_MSE);

    neural_set_init_w_b_MLP(neu,1);


    data_selector_randall_init(neu,200.0/number_data);

    forward_calculator_MLP_init(neu);

    gradient_calculator_MLP_init(neu,10,10000);

    weight_matrix_updator_Adam_init(neu,0.9,0.999);

    learn_ratio_exp_init(neu,0.005,0.0000001,500);

    exitor_avg_init(neu,10,0.001,1000);

    print_function_init(neu,1,10);

    neural_train(neu);



    return neu;

}

struct DATA * solid_liquid_ml(struct POSCAR * pos,struct NEURAL * neu,int N,int L,double ml_cutoff,double alpha,double beta,double neigh_cutoff,double rcut,double * melt)
{

    /////////////////////////////////////////

    struct DATA ** nets=build_coordinations_table_byRadius_typedecompose(pos,rcut);

    struct NEP * nep=nep_init(N,L,rcut);

    int num_all=pos->num_all;
    int num_type=pos->num_type;
    struct DATA * tmp_nep=nep_generate_nepdescribtor(nets,num_type,num_all,nep,rcut,false,true,false,true);


    double * input=(double *)calloc(num_all*N*N*L,sizeof(double));

    for(int i=0;i<num_all;i++)
    {
        for(int j=0;j<N*N*L;j++)
        {
            input[i*N*N*L + j] = DATA2D_get(tmp_nep,i,j);
        }
    }
    DATA_free(tmp_nep);
    

    neural_set_data(neu,input,NULL,num_all);

    neural_run(neu);

    ////////////////////////////////



    int solid=0;
    struct DATA * ret=DATA1D_init(num_all);
    for(int i=0;i<num_all;i++)
    {
        double pa=neu->output_data[ i*2 + 0    ];
        double pb=neu->output_data[ i*2 + 1    ];
        double avg=alpha*pa/(  alpha*pa + beta*(1-pb)      );
        if( avg- ml_cutoff > 0)
        {
            ret->data_data[i]=1.0;
            solid++;
        }

    }
    *melt=1.0*solid/num_all;

    ////////////////////////////////neigh
    for(int at=0;at<num_all;at++)
    {
        int total=0;
        int r_total=0;
        for(int ty=0;ty<num_type;ty++)
        {
            for(int ne=1;ne<nets[at*num_type + ty]->dimen_data[0];ne++)
            {
                int n_nei=(int)DATA2D_get(nets[at*num_type + ty],ne,7);
                n_nei=n_nei/num_type;
                int n_st=ret->data_data[n_nei];

                if(   (int)(ret->data_data[at]) != n_st    )
                {
                    r_total++;
                }
                total++;
            }
        }
        if(1.0*r_total/total - neigh_cutoff > 0)
        {
            if(ret->data_data[at] < 0.5 ) ret->data_data[at]=1.0;
            else ret->data_data[at]=0;
        }
    }
    solid=0;
    for(int at=0;at<num_all;at++)
    {
        if(ret->data_data[at] - 0.5 >0) solid++;
    }
    *melt=1.0*solid/num_all;


    nep_free(nep);
    DATAs_free(nets);
    free(input);
    input=NULL;
    return ret;


}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// 

double solidfication_check(struct SOLID_CHECK * so_check, struct PARAMETERS * inputs)
{
    double melt=-1;
    if(so_check->type == 1 || so_check->type == -1)//for qlmi similarity
    {    
        if(so_check->run_type == 0) // learn
        {
            if(so_check->best_table != NULL)
            {
                DATAs_free(so_check->best_table);
            }

            int L;
            double radius;
            if(so_check->type == 1)
            {
                L=6;
                radius=6.0;
            }
            else
            {
                L=(int)(inputs->atom_method[1]);
                radius=inputs->atom_method[4];
            }
            ///////////////////////////////////////////

            struct DATA ** nets=build_coordinations_table_byRadius_typedecompose(so_check->solid,radius);

            so_check->best_table=solid_liquid_bySjkStatistic_get_reference(nets,so_check->solid->num_all,so_check->solid->num_type,&(so_check->num_ref),L);

            DATAs_free(nets);

        }
        else if(so_check->run_type == 1)
        {
    
            double cutoff;
            double neigh_cutoff;
            double radius;
            if(so_check->type == 1)
            {
            
                cutoff=0.35;
                neigh_cutoff=0.5;
                radius=6.0;
            }
            else
            {
                cutoff=inputs->atom_method[2];
                neigh_cutoff=inputs->atom_method[3];      
                radius=inputs->atom_method[4];          
            }

            struct DATA ** nets=build_coordinations_table_byRadius_typedecompose(so_check->solid,radius);
            struct DATA * ret=solid_liquid_bySjkStatistic_from_reference(nets,so_check->solid->num_all,so_check->solid->num_type,so_check->best_table,so_check->num_ref,cutoff,neigh_cutoff,&melt);
            DATA_free(ret);
            DATAs_free(nets);
        }
    }
    else if(so_check->type == 2 || so_check->type == -2) //local density
    {
        if(so_check->run_type == 0) // learn
        {
            double sigma;
            double radius;
            if(so_check->type == 2)
            {
                sigma=0.4;
                radius=6.0;
            }
            else
            {
                sigma=inputs->atom_method[1];
                radius=inputs->atom_method[4];
            }
            struct DATA ** nets=build_coordinations_table_byRadius_typedecompose(so_check->solid,radius);
            so_check->refer=solid_liquid_byLD_get_reference(nets,so_check->solid->num_all,so_check->solid->num_type,&(so_check->num_ld_ref),&(so_check->refer_normal),sigma);
            DATAs_free(nets);
        }
        else if(so_check->run_type == 1)
        {
            double sigma;
            double ld_cutoff;
            double neigh_cutoff;
            double radius;


            if(so_check->type == 2)
            {
                sigma=0.4;
                ld_cutoff=0.3;
                neigh_cutoff=0.5;
                radius=6.0;
            }
            else
            {
                sigma=inputs->atom_method[1];
                ld_cutoff=inputs->atom_method[2];
                neigh_cutoff=inputs->atom_method[3];
                radius=inputs->atom_method[4];                
            }


            struct DATA ** nets=build_coordinations_table_byRadius_typedecompose(so_check->solid,radius);

            struct DATA * ret=solid_liquid_by_LD_from_reference(nets,so_check->solid->num_all,so_check->solid->num_type,so_check->refer,so_check->refer_normal,so_check->num_ld_ref,ld_cutoff,neigh_cutoff,sigma,&melt);
            DATA_free(ret);
            DATAs_free(nets);
        }        
    }
    else if(so_check->type == 3 || so_check->type == -3) //ml
    {
        if(so_check->run_type == 0) // learn
        {
            double L;
            double sigma;
            double rcut;

            if(so_check->type == 3)
            {
                L=6;
                sigma=0.4;
                rcut=6.0;
            }
            else
            {
                L=(int)(inputs->atom_method[5]);
                sigma=inputs->atom_method[3];  
                rcut=inputs->atom_method[7];              
            }

            int numbers=solid_liquid_ml_getNEP_and_Qlmi_LD_Neigh_all(so_check->solid,so_check->poss_solid,so_check->poss_liquid,so_check->num_solid,so_check->num_liquid,6,16,L,sigma,rcut,1000,1000,&(so_check->input_data),&(so_check->output_data));
            printf("........................................................................................start.\n");
            so_check->neu=solid_liquid_ml_train(so_check->input_data,so_check->output_data,numbers,6*6*16,2);
            free(so_check->input_data);
            free(so_check->output_data);
            so_check->input_data=NULL;
            so_check->output_data=NULL;
        }
        else if(so_check->run_type == 1)
        {
            double ml_cutoff;
            double neigh_cutoff;
            double alpha;
            double beta;
            double rcut;
            if(so_check->type == 3)
            {
                ml_cutoff=0.35;
                neigh_cutoff=0.5;
                alpha=0.5;
                beta=0.5;      
                rcut=6.0;          
            }
            else
            {
                ml_cutoff=inputs->atom_method[5];
                neigh_cutoff=inputs->atom_method[6];
                alpha=inputs->atom_method[1];
                beta=inputs->atom_method[2];  
                rcut=inputs->atom_method[7];               
            }

            struct DATA * ret=solid_liquid_ml(so_check->solid,so_check->neu,6,16,ml_cutoff,alpha,beta,neigh_cutoff,rcut,&melt);
            DATA_free(ret);
        }        
    }

    return melt;
}
struct SOLID_CHECK * solidfication_init(int type)
{
    struct SOLID_CHECK * out=(struct SOLID_CHECK *)calloc(1,sizeof(struct SOLID_CHECK));
    out->type=type;
    out->run_type=-1;
    out->solid=NULL;
    out->liquid=NULL;
    out->poss_solid=NULL;
    out->poss_liquid=NULL;


    out->best_table=NULL;

    out->input_data=NULL;
    out->output_data=NULL;


    out->refer=NULL;
    out->refer_normal=NULL;
    
    out->neu=NULL;
    return out;
}
void solidfication_free(struct SOLID_CHECK * so_check)
{

    if(so_check->best_table != NULL)
    {
        DATAs_free(so_check->best_table);
    }
        
    

    if(so_check->refer_normal != NULL) DATA_free(so_check->refer_normal);
    if(so_check->refer != NULL) DATAs_free(so_check->refer);
    

    if(so_check->neu != NULL)
    {
        neural_free(so_check->neu);
    }
    
    free(so_check);
    so_check=NULL;
}



