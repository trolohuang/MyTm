#include "nep.h"
double Heaviside(double x)
{
    if(x>0) return 1;
    else if(x==0) return 0.5;
    else if(x<0) return 0;
    
    return -1;
}







double nep_fc_function(double r,struct NEP * nep)
{
    return pow(1-r/(nep->Rcut),nep->NP) * Heaviside(nep->Rcut - r);
}

double nep_gr_function(int n,double r,struct NEP * nep)
{
    return pow(r,n)*nep_fc_function(r,nep);
}

////////////////////////////////////////////////////////////////////////
// basic calculation for NEP
struct NEP * nep_init(int N, int L,double rcut)
{
    struct NEP * out=(struct NEP *)calloc(1,sizeof(struct NEP));

    out->N=N;
    out->L=L;
    out->NP=5;
    out->Rcut=rcut;
    /////////////////
    out->cnlm=DATA4D_init(N,L,2*(L - 1)+1,2) ;
    out->Gi2=DATA1D_init(N);
    out->Gi3=DATA3D_init(N,N,L);
    out->neigh=NULL;
    return out;
}
void nep_free(struct NEP * nep)
{
    if(nep->Gi2 != NULL) DATA_free(nep->Gi2);
    if(nep->Gi3 != NULL) DATA_free(nep->Gi3);
    if(nep->cnlm != NULL) DATA_free(nep->cnlm);
    free(nep);
    nep=NULL;
}
void nep_generate_Gi2(struct NEP * nep,double scaled)
{
    if(nep->neigh != NULL)
    {
        for(int i=0;i<nep->N;i++)
        {
            double sum=0;
            for(int at=1;at<nep->neigh->dimen_data[0];at++)
            {
                sum+=nep_gr_function(i,  DATA2D_get(nep->neigh,at,3)/scaled   ,nep);
            }
            nep->Gi2->data_data[i]=sum;
        }
    }
}
void nep_generate_Gi3(struct NEP * nep,double scaled)
{
    if(nep->neigh != NULL)
    {

        //////////////////////////////////////////////////////////////
        double x0,y0,z0;
        x0=DATA2D_get(nep->neigh,0,0);
        y0=DATA2D_get(nep->neigh,0,1);
        z0=DATA2D_get(nep->neigh,0,2);
        

        //////////////////////////////////////////////////////////////
        struct DATA * cnlm=nep->cnlm;
        //for cnlm
        for(int n=0;n<nep->N;n++)
        {
            for(int l=0;l<nep->L;l++)
            {
                for(int m=-l;m<=l;m++)
                {
                    double real=0;
                    double img=0;
                    for(int at=1;at<nep->neigh->dimen_data[0];at++)
                    {
                        double x1=DATA2D_get(nep->neigh,at,0);
                        double y1=DATA2D_get(nep->neigh,at,1);
                        double z1=DATA2D_get(nep->neigh,at,2);
                        double r,theta,phi;
                        cartesian2spherical(x1-x0,y1-y0,z1-z0,&r,&theta,&phi);
                        r=r/scaled;
                        /////////////////////////
                        double rasf=nep_gr_function(n,r,nep);
                        double yhr,yhi;
                        Sph_Har(l,m,theta,phi,&yhr,&yhi);
                        real+=yhr*rasf;
                        img+=yhi*rasf;
                    }
                    DATA4D_set(cnlm,real,n,l,m+l,0);
                    DATA4D_set(cnlm,img,n,l,m+l,1);
                }
            }
        }

        // for Gi3
        for(int n1=0;n1<nep->N;n1++)
        {
            for(int n2=n1;n2<nep->N;n2++)
            {
                for(int l=0;l<nep->L;l++)
                {
                    double sum_real=0;
                    /////////////////////////////
                    for(int m= -l ;m<=l ;m++)
                    {
                        double r1=DATA4D_get(cnlm,n1,l,m+l,0);
                        double i1=DATA4D_get(cnlm,n1,l,m+l,1);
                        double r2=DATA4D_get(cnlm,n2,l,m+l,0);
                        double i2=DATA4D_get(cnlm,n2,l,m+l,1);
                        sum_real+=r1*r2-i1*i2;
                    }
                    ////////////////////////////
                    if(n1 == n2)
                    {

                        DATA3D_set(nep->Gi3,sum_real,n1,n2,l);                    
                    }
                    else
                    {
                        DATA3D_set(nep->Gi3,sum_real,n1,n2,l);
                        DATA3D_set(nep->Gi3,sum_real,n2,n1,l);
                    }
                }
            }
        }
    }
}


////////////////////////////////////////////////////////////////////////
//generator describtor
struct DATA * nep_generate_nepdescribtor(struct DATA ** nets_decompose,int atom_num_type,int atom_num_all,struct NEP * nep,double scaled,bool gen_Gi2,bool gen_Gi3,bool normal_Gi2,bool normal_Gi3)
{
    int N=nep->N;
    int L=nep->L;
    struct  DATA * output=NULL;
    if(gen_Gi2 && gen_Gi3)
    {
        output=DATA2D_init(atom_num_all,  atom_num_type*(N+N*N*L)   );
        if(nep->Gi2 == NULL) nep->Gi2=DATA1D_init(N);
        if(nep->Gi3 == NULL) nep->Gi3=DATA3D_init(N,N,L);
        ///////////////////////////////////////////////////////////////
        for(int at=0;at<atom_num_all;at++)
        {
            for(int dt=0;dt<atom_num_type;dt++)
            {

                nep->neigh=nets_decompose[at*atom_num_type + dt];
                if(nets_decompose[at*atom_num_type + dt]->dimen_data[0] <=1)
                {       
                    memset(nep->Gi2->data_data,0,N*sizeof(double));
                    memset(nep->Gi3->data_data,0,N*N*L*sizeof(double));
                }
                else
                {
                    nep_generate_Gi2(nep,scaled);
                    nep_generate_Gi3(nep,scaled);
                }

                if(normal_Gi2)
                {
                    nep_Gi2_normal(nep->Gi2);
                }
                if(normal_Gi3)
                {
                    nep_Gi3_normal(nep->Gi3);
                }
                
                for(int i=0;i<N;i++)
                {
                    DATA2D_set(output,nep->Gi2->data_data[i],at,dt*(N+N*N*L) + i);
                }

                for(int i=0;i<N*N*L;i++)
                {
                    DATA2D_set(output,nep->Gi3->data_data[i],at,dt*(N+N*N*L) + N + i);
                }
            }
        }
        ////////////////////////////////////////////////////////////////
    }
    else if(gen_Gi2 && !gen_Gi3 )
    {
        output=DATA2D_init(atom_num_all,  atom_num_type*(N)   );
        if(nep->Gi2 == NULL) nep->Gi2=DATA1D_init(N);
        ///////////////////////////////////////////////////////////////
        for(int at=0;at<atom_num_all;at++)
        {
            for(int dt=0;dt<atom_num_type;dt++)
            {

                nep->neigh=nets_decompose[at*atom_num_type + dt];
                if(nets_decompose[at*atom_num_type + dt]->dimen_data[0] <=1)
                {       
                    memset(nep->Gi2->data_data,0,N*sizeof(double));

                }
                else
                {
                    nep_generate_Gi2(nep,scaled);
                }

                if(normal_Gi2)
                {
                    nep_Gi2_normal(nep->Gi2);
                }

                for(int i=0;i<N;i++)
                {
                    DATA2D_set(output,nep->Gi2->data_data[i],at,dt*(N) + i);
                }

            }
        }
        ////////////////////////////////////////////////////////////////
    }
    else if(!gen_Gi2 && gen_Gi3)
    {
        output=DATA2D_init(atom_num_all,  atom_num_type*(N*N*L)   );
        if(nep->Gi3 == NULL) nep->Gi3=DATA3D_init(N,N,L);
        ///////////////////////////////////////////////////////////////
        for(int at=0;at<atom_num_all;at++)
        {
            for(int dt=0;dt<atom_num_type;dt++)
            {

                nep->neigh=nets_decompose[at*atom_num_type + dt];
                if(nets_decompose[at*atom_num_type + dt]->dimen_data[0] <=1)
                {       
                    memset(nep->Gi3->data_data,0,N*N*L*sizeof(double));
                }
                else
                {

                    nep_generate_Gi3(nep,scaled);
                }


                if(normal_Gi3)
                {
                    nep_Gi3_normal(nep->Gi3);
                }
                
                for(int i=0;i<N*N*L;i++)
                {
                    DATA2D_set(output,nep->Gi3->data_data[i],at,dt*(N*N*L) + i);
                }

            }
        }
        ////////////////////////////////////////////////////////////////
    }
    return output;
}

////////////////////////////////////////////////////////////////////////
double nep_kernal_Gi2(struct DATA * Gi2_1,struct DATA * Gi2_2)
{
    double sum=0;
    for(int i=0;i<Gi2_1->all_data;i++)
    {
        sum+=Gi2_1->data_data[i] * Gi2_2->data_data[i];
    }
    return sum;
}   
double nep_kernal_Gi3(struct DATA * Gi3_1,struct DATA * Gi3_2)
{
    double sum=0;
    for(int i=0;i<Gi3_1->all_data;i++)
    {
        sum+=Gi3_1->data_data[i] * Gi3_2->data_data[i];
    }
    return sum;
}  

void nep_Gi2_normal(struct DATA * Gi2)
{
    double sum=0;
    for(int i=0;i<Gi2->all_data;i++)
    {
        sum+=Gi2->data_data[i] * Gi2->data_data[i];
    }    
    if(sum>0)
    {
        sum=sqrt(sum);
        for(int i=0;i<Gi2->all_data;i++)
        {
            Gi2->data_data[i] = Gi2->data_data[i]/sum;
        }   
    } 
}
void nep_Gi3_normal(struct DATA * Gi3)
{
    double sum=0;
    for(int i=0;i<Gi3->all_data;i++)
    {
        sum+=Gi3->data_data[i] * Gi3->data_data[i];
    }    
    if(sum>0)
    {
        sum=sqrt(sum);
        for(int i=0;i<Gi3->all_data;i++)
        {
            Gi3->data_data[i] = Gi3->data_data[i]/sum;
        }   
    }
}






