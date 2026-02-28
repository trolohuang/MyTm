#include "soap.h"
//径向函数
//高斯轨道基
double soap_gr_function_gauss(struct SOAP * soap)
{
    double sigma=soap->gr_para->data_data[0];
    int l=(int)(soap->gr_para->data_data[1]);
    int n=(int)(soap->gr_para->data_data[2]);
    double r=(soap->gr_para->data_data[3]);
    double Nn=DATA2D_get(soap->gr_data,n,l);

    return Nn*pow(r,l)*exp( -(  pow(r-  1.0*(n+1)/(soap->N_max+1)*soap->Rcut ,2 )/(2*sigma*sigma)   )    )*soap_fc_function_cos(r,soap->Rcut);
}

//截断函数
double soap_fc_function_cos(double r,double rcut)
{
    if(r - rcut <=0)  return 0.5*(  cos(M_PI*r/rcut) + 1   );
    else return 0; 
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

struct SOAP * soap_init(int N,int L,double Rcut,double sigma_density)
{
    struct SOAP * out=(struct SOAP *)calloc(1,sizeof(struct SOAP));
    out->N_max=N;
    out->L_max=L;
    out->Rcut=Rcut;
    out->sigma_density=sigma_density;

    ///////////////////////////////////////for gr
    out->gr_data=NULL;
    out->gr_para=NULL;
    ///////////////////////////////////////for result


    out->pnnl=DATA3D_init(N,N,L);
    out->cnlm=DATA4D_init(N,L,2*L-1 ,2);
    out->integrade_table=NULL;
    return out;
}
void soap_gr_init(struct SOAP * soap,int gr_type,double para1,double para2)
{
    if(gr_type == 1) // for gauss type
    {
        if(soap->gr_para == NULL)
        {
            soap->gr_para=DATA1D_init(4); // sigma, l , n , r
            soap->gr_para->data_data[0]=para1; // N x L
        } 
        else
        {
            DATA_free(soap->gr_para);
            soap->gr_para=DATA1D_init(4);
            soap->gr_para->data_data[0]=para1; // N x L
        }
        if(soap->gr_data == NULL)  
        {
            soap->gr_data=DATA2D_init(soap->N_max,soap->L_max); // for Nn
        }
        else
        {
            DATA_free(soap->gr_data);
            soap->gr_data=DATA2D_init(soap->N_max,soap->L_max); // for Nn
        }
        
        soap->soap_gr=soap_gr_function_gauss;

        int n_inte=(soap->Rcut )/0.01;
        for(int n=0;n<soap->N_max;n++)
        {
            for(int l=0;l<soap->L_max;l++)
            {
                double n_val=0;
                /////////////////////
                double now_val;
                double pre_val;
                for(int i=0;i<n_inte;i++)
                {
                    double r=i*0.01;
                    now_val=pow(r,2*l+2)*exp(   -2*( r- (n+1)/(soap->N_max + 1)*soap->Rcut    )/( 2*para1*para1  )  ) * pow( soap_fc_function_cos(r,soap->Rcut),2);
                    if(i == 0) pre_val=now_val;
                    n_val+=(pre_val+now_val)*0.01*0.5;
                }
                n_val=pow(n_val,-0.5);
                /////////////////////
                DATA2D_set(soap->gr_data,n_val,n,l);
            }
        }

    }
}
void soap_free(struct SOAP * soap)
{
    if(soap->integrade_table != NULL) DATA_free(soap->integrade_table);
    if(soap->gr_data != NULL) DATA_free(soap->gr_data);
    if(soap->gr_para != NULL) DATA_free(soap->gr_para);
    if(soap->cnlm != NULL) DATA_free(soap->cnlm);
    if(soap->pnnl != NULL) DATA_free(soap->pnnl);
    free(soap);
    soap=NULL;
}




void soap_integration_table(struct SOAP * soap,int bin_r,int bin_integration)
{
    if(soap->integrade_table != NULL) DATA_free(soap->integrade_table);
    struct DATA * output=DATA3D_init(soap->N_max ,soap->L_max ,bin_r);
    soap->integrade_table=output;

    double prec_r=soap->Rcut/(bin_r);

    double inte_prec_r=soap->Rcut/(bin_integration);
    for(int nl=0;nl<soap->L_max;nl++)
    {
        for(int nn=0;nn<soap->N_max;nn++)
        {
            for(int nr=0;nr<bin_r;nr++)
            {
                ////////////////////////////
                double rij=nr*prec_r;

                double inte_sum=0;
                double pre_val=0;
                double now_val=0;
                for(int i=0;i<bin_integration;i++)
                {
                    double r=i*inte_prec_r;
                    /////////////////////
                    soap->gr_para->data_data[1]=nl;
                    soap->gr_para->data_data[2]=nn;
                    soap->gr_para->data_data[3]=r;
                    double gr= soap->soap_gr(soap);
                    /////////////////////
                    now_val= pow(r,2) * gr * gsl_sf_bessel_il_scaled( nl , r*rij/pow( soap->sigma_density,2)    )*exp(    (2*r*rij - (r*r + rij*rij) )/(2*pow(soap->sigma_density,2)) );
                    if(i==0) pre_val=now_val;
                    inte_sum+=(now_val + pre_val)*inte_prec_r*0.5;
                    pre_val=now_val;
                }
                ////////////////////////
                DATA3D_set(output,inte_sum,nn,nl,nr);

            }

            ///////////////////////////////

        }
    }


}

void soap_cnlm_generate_withScaled(struct SOAP * soap,double scale_dis)
{
    double x0,y0,z0;
    double x,y,z;
    double r,theta,phi;
    double val_real=0;
    double val_img=0;
    x0=DATA2D_get(soap->neig_atoms,0,0);
    y0=DATA2D_get(soap->neig_atoms,0,1);
    z0=DATA2D_get(soap->neig_atoms,0,2);    
    for(int n=0;n<soap->N_max;n++)
    {
        for(int l=0;l<soap->L_max;l++)
        {
            for(int m= -l ;m<=l ;m++)
            {
                
                val_real=0;
                val_img=0;
                ////////////////////////
                for(int at=1;at<soap->neig_atoms->dimen_data[0];at++)
                {
                    x=DATA2D_get(soap->neig_atoms,at,0);
                    y=DATA2D_get(soap->neig_atoms,at,1);
                    z=DATA2D_get(soap->neig_atoms,at,2);
                    
                    ///////////////
                    cartesian2spherical(x-x0,y-y0,z-z0,&r,&theta,&phi);
                    r=r/scale_dis;
                    ///
                    double tmp_val=4*M_PI*soap_fc_function_cos(r,soap->Rcut);
                    //
                    int N_R=(int)(    (r/soap->Rcut )*soap->integrade_table->dimen_data[2]        );
                    tmp_val*=DATA3D_get(soap->integrade_table,n,l,N_R);
                    //
                    double real,img;

                    Sph_Har(l,m,theta,phi,&real,&img);

                    real=real*tmp_val;
                    img=img*tmp_val;
                    //
                    val_real+=real;
                    val_img+=img;
                }
                ////////////////////////
                DATA4D_set(soap->cnlm,val_real,n,l,m+(soap->L_max-1),0);
                DATA4D_set(soap->cnlm,val_img,n,l,m+(soap->L_max-1),1);
            }
        }
    }
}

void soap_pnnl_generate(struct SOAP * soap)
{

    double total=0;
    for(int n1=0;n1<soap->N_max;n1++)
    {
        for(int n2=n1;n2<soap->N_max;n2++)
        {
            for(int l=0;l<soap->L_max;l++)
            {
                double sum_real=0;
                /////////////////////////////
                for(int m= -(soap->L_max-1) ;m<=(soap->L_max-1) ;m++)
                {
                    double r1=DATA4D_get(soap->cnlm,n1,l,m+(soap->L_max-1),0);
                    double i1=DATA4D_get(soap->cnlm,n1,l,m+(soap->L_max-1),1);
                    double r2=DATA4D_get(soap->cnlm,n2,l,m+(soap->L_max-1),0);
                    double i2=DATA4D_get(soap->cnlm,n2,l,m+(soap->L_max-1),1);
                    sum_real+=r1*r2-i1*i2;
                }
                ////////////////////////////
                if(n1 == n2)
                {
                    total+=sum_real*sum_real;
                    DATA3D_set(soap->pnnl,sum_real,n1,n2,l);                    
                }
                else
                {
                    total+=2*sum_real*sum_real;
                    DATA3D_set(soap->pnnl,sum_real,n1,n2,l);
                    DATA3D_set(soap->pnnl,sum_real,n2,n1,l);
                }
            }
        }
    }
    total=pow(total,0.5);

    for(int i=0;i<soap->pnnl->all_data;i++)
    {
        soap->pnnl->data_data[i]= (soap->pnnl->data_data[i])/total;
    }   

}

double soap_kernal(struct DATA * pnnl_1,struct DATA * pnnl_2)
{
    double sum=0;
    for(int i=0;i<pnnl_1->all_data;i++)
    {
        sum+=(pnnl_1->data_data[i])*(pnnl_2->data_data[i]);
    }
    return sum;
}







