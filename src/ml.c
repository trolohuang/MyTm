#include "ml.h"
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//神经网络使用的矩阵函数
// A x B + C
void matrix_A_X_B_add_C(double * A,int ma,int na,double *B,int mb,int nb,double *C,double * res)//WXa+B
{
    if(na == mb)
    {
        for(int m1=0;m1<ma;m1++)
        {
            for(int n2=0;n2<nb;n2++)
            {
                res[m1*nb + n2]=0;
                for(int i=0;i<na;i++)
                {
                    res[m1*nb + n2]+=A[ m1*na + i  ] * B[ i*nb + n2 ];
                }
                res[m1*nb + n2]+=C[m1*nb + n2];
            }
        }
    }

}

// A x B 
void matrix_A_X_B(double * A,int ma,int na,double *B,int mb,int nb,double * res)//WXa+B
{
    if(na == mb)
    {
        for(int m1=0;m1<ma;m1++)
        {
            for(int n2=0;n2<nb;n2++)
            {
                res[m1*nb + n2]=0;
                for(int i=0;i<na;i++)
                {
                    res[m1*nb + n2]+=A[ m1*na + i  ] * B[ i*nb + n2 ];
                }
            }
        }
    }

}

// W转至 x B 点乘 C
void matrix_WT_X_B_adX_C(double *w,int mw,int nw,double *B,int mb,int nb,double *C,double *res)//Wt X B .O C
{
    if(mw == mb)
    {
        for(int n1=0;n1<nw;n1++)
        {
            for(int n2=0;n2<nb;n2++)
            {
                res[n1*nb + n2]=0;
                for(int i=0;i<mw;i++)
                {
                    res[n1*nb + n2]+=w[   i*nw + n1      ] * B[  i*nb + n2 ];
                }
                res[n1*nb + n2]=res[n1*nb + n2]*C[n1*nb + n2];
            }
        }
    }
}

// A x B转至
void matrix_A_X_BT(double * A,int ma,int na,double *B,int mb,int nb,double * res) // A X BT
{
    if(na == nb)
    {
        for(int m1=0;m1<ma;m1++)
        {
            for(int m2=0;m2<mb;m2++)
            {
                res[m1*mb + m2]=0;
                for(int i=0;i<na;i++)
                {
                    res[m1*mb + m2]+=A[ m1*na + i  ] * B[ m2*nb + i ];
                }
            }
        }
    }   
}

//矩阵A全部元素调用函数activate_func
void matrix_A_ac(double * A,int ma,int na,double *res,double (*activate_func)(double))
{
    for(int m=0;m<ma;m++)
    {
        for(int n=0;n<na;n++)
        {
            res[m*na+n]=activate_func(A[m*na+n]);
        }
    }
}

//矩阵A全部元素调用函数G_activate_func
void matrix_A_Gac(double * A,int ma,int na,double *res,double (*G_activate_func)(double))
{
    for(int m=0;m<ma;m++)
    {
        for(int n=0;n<na;n++)
        {
            res[m*na+n]=G_activate_func(A[m*na+n]);
        }
    }
}

// vector A . B
double matrix_A_dot_B(double * A, double *B, int numbers) // vector A.B
{
    double sums=0;
    for(int i=0;i<numbers;i++)
    {
        sums+=A[i]*B[i];
    }
    return sums;
}

// vector C = A .. B  ；C=A，B按元素相乘
void matrix_A_dd_B(double *A,double *B,int numbers,double *C)
{
    for(int i=0;i<numbers;i++)
    {
        C[i]=A[i]*B[i];
    }
}

//常数 a. Vector B
void matrix_a_dot_B(double a,double *B,int numbers,double *C)
{
    for(int i=0;i<numbers;i++)
    {
        C[i]=a*B[i];
    }
}

double box_muller(double x0,double sigma)
{
    return x0 + sigma*(     sqrt(    -2*(log(   0.000000001 + 1.0*(rand()%10000)/10000        )  )        ) * cos(2*3.1415926535*(  0.000000001 + 1.0*(rand()%10000)/10000  ))       );
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void neu_printf(struct NEURAL * neu)
{
    printf("\n******************************************************************************\n");
    struct NEURAL_NETWORK * tmp=neu->neu;
    int num=0;
while(tmp != NULL)
{
    printf("-----------------------------%d\n",num);
    num++;
    printf("\nw\n");
    for(int i=0;i<tmp->m*tmp->n;i++)
    {
        printf("%lf ",tmp->w[i]);
    }
    printf("\n");
    printf("\nlw\n");
    for(int i=0;i<tmp->m*tmp->n;i++)
    {
        printf("%lf ",tmp->lw[i]);
    }
    printf("\n");
    printf("\nb\n");
    for(int i=0;i<tmp->m;i++)
    {
        printf("%lf ",tmp->b[i]);
    }
    printf("\n");
    printf("\nlb\n");
    for(int i=0;i<tmp->m;i++)
    {
        printf("%lf ",tmp->lb_sigma[i]);
    }
    printf("\n");
    printf("\nz\n");
    for(int i=0;i<tmp->m;i++)
    {
        printf("%lf ",tmp->z[i]);
    }
    printf("\n");
    printf("\na\n");
    for(int i=0;i<tmp->m;i++)
    {
        printf("%lf ",tmp->a[i]);
    }
    printf("\n");
    printf("\naa\n");
    for(int i=0;i<tmp->m;i++)
    {
        printf("%lf ",tmp->aa[i]);
    }
    printf("\n");
    tmp=tmp->next;
}
}
void neural_network_free(struct NEURAL_NETWORK * neus)
{
    if(neus->a != NULL) free(neus->a);
    if(neus->aa != NULL) free(neus->aa);
    if(neus->z != NULL) free(neus->z);
    if(neus->w != NULL) free(neus->w);

    if(neus->lw != NULL) free(neus->lw);
    if(neus->lb_sigma != NULL) free(neus->lb_sigma);

    if(neus->mbt != NULL) free(neus->mbt);
    if(neus->vbt != NULL) free(neus->vbt);
    if(neus->mt != NULL) free(neus->mt);
    if(neus->vt != NULL) free(neus->vt);

    free(neus);
    neus=NULL;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

struct NEURAL * neural_init_MLP(int input_dimen,int output_dimen,int layers,...)
{
    va_list arg;
    va_start(arg,layers);
    int * args=NULL;
    args=(int *)calloc(layers,sizeof(int));
    for(int i=0;i<layers;i++)
    {
        args[i]=va_arg(arg,int);
        
    }
    va_end(arg);
    /////////////////////////////////

    struct NEURAL * out=(struct NEURAL *)calloc(1,sizeof(struct NEURAL));

    ///////////////////////////////
    struct NEURAL_NETWORK * tmp;
    struct NEURAL_NETWORK * tmp2;
    if(layers>=1)
    {
        for(int i=0;i<layers;i++)
        {
            int units=args[i];
            tmp=(struct NEURAL_NETWORK *)calloc(1,sizeof(struct NEURAL_NETWORK ));
            tmp->next=NULL;
            tmp->previous=NULL;

            //////////////////////

            if(i==0)
            {
                out->neu=tmp;
                tmp2=out->neu;

                tmp2->m=units;
                tmp2->n=input_dimen;

                tmp2->w=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->b=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->lw=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->lb_sigma=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->z=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->a=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->aa=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->mt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->vt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->mbt=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->vbt=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->activate_function=activate_function_default;
                tmp2->G_activate_function=G_activate_function_default;


            }
            else
            {
                tmp2->next=tmp;
                tmp->previous=tmp2;
                tmp2=tmp2->next;

                tmp2->m=units;
                tmp2->n=args[i-1];

                tmp2->w=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->b=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->lw=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->lb_sigma=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->z=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->a=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->aa=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->mt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->vt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->mbt=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->vbt=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->activate_function=activate_function_default;
                tmp2->G_activate_function=G_activate_function_default;
            }
        }

        //output layers
        tmp=(struct NEURAL_NETWORK *)calloc(1,sizeof(struct NEURAL_NETWORK ));
        tmp->next=NULL;
        tmp->previous=NULL;

        tmp2->next=tmp;
        tmp->previous=tmp2;
        tmp2=tmp2->next;

        tmp2->m=output_dimen;
        tmp2->n=args[layers-1];

        tmp2->w=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
        tmp2->b=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->lw=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
        tmp2->lb_sigma=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->z=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->a=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->aa=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->mt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
        tmp2->vt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
        tmp2->mbt=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->vbt=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->activate_function=activate_function_default;
        tmp2->G_activate_function=G_activate_function_default;

    }
    else if(layers == 0)
    {
        tmp2=(struct NEURAL_NETWORK *)calloc(1,sizeof(struct NEURAL_NETWORK ));
        tmp2->next=NULL;
        tmp2->previous=NULL;

        tmp2->m=output_dimen;
        tmp2->n=input_dimen;

        tmp2->w=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
        tmp2->b=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->lw=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
        tmp2->lb_sigma=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->z=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->a=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->aa=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->mt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
        tmp2->vt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
        tmp2->mbt=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->vbt=(double *)calloc(tmp2->m,sizeof(double));
        tmp2->activate_function=activate_function_default;
        tmp2->G_activate_function=G_activate_function_default;      
        out->neu=tmp2;  
    }
    ///////////////////////////////



    out->input_units=input_dimen;
    out->output_units=output_dimen;
    out->layers=layers+1;


    out->in_m=-1;
    out->in_n=-1;
    out->input_data=NULL;
    out->out_m=-1;
    out->out_n=-1;
    out->output_data=NULL;
    out->tr_m=-1;
    out->tr_n=-1;
    out->train_data=NULL;

    out->para_data_sel_func=NULL;
    out->para_exitor=NULL;
    out->para_forward_calculator=NULL;
    out->para_gradient_calculator=NULL;
    out->para_learn_function=NULL;
    out->para_weight_matrix_updator=NULL;
    out->loss_data=NULL;
    out->para_print_function=NULL;
    out->para_wither_birth_creator=NULL;

    out->data_sel_func=NULL;
    out->exitor=NULL;
    out->forward_calculator=NULL;
    out->gradient_calculator=NULL;
    out->learn_function=NULL;
    out->weight_matrix_updator=NULL;
    out->print_function=NULL;
    out->wither_birth_creator=NULL;
    


    free(args);
    return out;

}

struct NEURAL * neural_init_ASLP(int input_dimen,int layers)
{
    /////////////////////////////////

    struct NEURAL * out=(struct NEURAL *)calloc(1,sizeof(struct NEURAL));

    ///////////////////////////////
    struct NEURAL_NETWORK * tmp;
    struct NEURAL_NETWORK * tmp2;

        for(int i=0;i<layers;i++)
        {

            tmp=(struct NEURAL_NETWORK *)calloc(1,sizeof(struct NEURAL_NETWORK ));
            tmp->next=NULL;
            tmp->previous=NULL;

            //////////////////////

            if(i==0)
            {
                out->neu=tmp;
                tmp2=out->neu;

                tmp2->m=1;
                tmp2->n=input_dimen;

                tmp2->w=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->b=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->lw=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->lb_sigma=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->z=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->a=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->aa=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->mt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->vt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->mbt=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->vbt=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->activate_function=activate_function_default;
                tmp2->G_activate_function=G_activate_function_default;


            }
            else
            {
                tmp2->next=tmp;
                tmp->previous=tmp2;
                tmp2=tmp2->next;

                tmp2->m=1;
                tmp2->n=input_dimen;

                tmp2->w=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->b=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->lw=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->lb_sigma=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->z=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->a=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->aa=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->mt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->vt=(double *)calloc((tmp2->m)*(tmp2->n),sizeof(double));
                tmp2->mbt=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->vbt=(double *)calloc(tmp2->m,sizeof(double));
                tmp2->activate_function=activate_function_default;
                tmp2->G_activate_function=G_activate_function_default;
            }
        }


    ///////////////////////////////



    out->input_units=input_dimen;
    out->output_units=1;
    out->layers=layers;


    out->in_m=-1;
    out->in_n=-1;
    out->input_data=NULL;
    out->out_m=-1;
    out->out_n=-1;
    out->output_data=NULL;
    out->tr_m=-1;
    out->tr_n=-1;
    out->train_data=NULL;

    out->para_data_sel_func=NULL;
    out->para_exitor=NULL;
    out->para_forward_calculator=NULL;
    out->para_gradient_calculator=NULL;
    out->para_learn_function=NULL;
    out->para_weight_matrix_updator=NULL;
    out->loss_data=NULL;
    out->para_print_function=NULL;
    out->para_wither_birth_creator=NULL;

    out->data_sel_func=NULL;
    out->exitor=NULL;
    out->forward_calculator=NULL;
    out->gradient_calculator=NULL;
    out->learn_function=NULL;
    out->weight_matrix_updator=NULL;
    out->print_function=NULL;
    out->wither_birth_creator=NULL;
    
    


    return out;

}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void neural_set_init_w_b_MLP(struct NEURAL * neu,int type)
{
    struct NEURAL_NETWORK * tmp_neu=neu->neu;
    srand(time(NULL));
    while(tmp_neu!=NULL)
    {
        double set=0;
        if(type == 0)// He /kaiming 
        {
            set=1.0/sqrt(tmp_neu->m);
            for(int i=0;i<(tmp_neu->m)*(tmp_neu->n);i++)
            {
                tmp_neu->w[i]=-set + 2*rand()%10000/( 10000/set  );
            }
        }
        else if(type == 1) // Xavier/ Glorot
        {
            set=1.0/sqrt(tmp_neu->m + tmp_neu->n);
            for(int i=0;i<(tmp_neu->m)*(tmp_neu->n);i++)
            {
                tmp_neu->w[i]=-set + 2*rand()%10000/( 10000/set  );
            }
        }
        tmp_neu=tmp_neu->next;
    }
}
void neural_set_init_w_b_const(struct NEURAL * neu,double val)
{
    struct NEURAL_NETWORK * tmp_neu=neu->neu;
    while(tmp_neu!=NULL)
    {
        for(int i=0;i<(tmp_neu->m)*(tmp_neu->n);i++)
        {
            tmp_neu->w[i]=val;
        }
        tmp_neu=tmp_neu->next;
    }
}
void neural_set_activate_function(struct NEURAL * neu,double (*activate_function)(double),double (*G_activate_function)(double),int set_layers)
{
    if(set_layers==0)
    {
        struct NEURAL_NETWORK * tmp=neu->neu;
        while(tmp!=NULL)
        {
            tmp->activate_function=activate_function;
            tmp->G_activate_function=G_activate_function;
            tmp=tmp->next;
        }
    }
    else if(set_layers>0)
    {
        struct NEURAL_NETWORK * tmp=neu->neu;
        set_layers=set_layers-1;
        while(set_layers--)
        {
            tmp=tmp->next;
        }
        tmp->activate_function=activate_function;
        tmp->G_activate_function=G_activate_function;
    }
    else if(set_layers <0)
    {
        struct NEURAL_NETWORK * tmp=neu->neu;
        while(tmp->next != NULL)
        {
            tmp=tmp->next;
        }        
        ///////////////////
        set_layers=abs(set_layers);
        set_layers--;
        while(set_layers--)
        {
            tmp=tmp->previous;
        }        
        tmp->activate_function=activate_function;
        tmp->G_activate_function=G_activate_function;        


    }
}
void neural_set_loss_func(struct NEURAL * neu,  double (*loss_func)(double * , double *, int , double ,double) , void (* G_loss_func)(double * , double * , int ,double * ,double ,double))
{
    neu->loss_func=loss_func;
    neu->G_loss_func=G_loss_func;
}
void neural_set_data(struct NEURAL * neu,double * input_data, double * reference_data,int data_numbers)
{
    neu->input_data=input_data;
    neu->in_m=data_numbers;
    neu->in_n=neu->input_units;

    neu->train_data=reference_data;
    neu->tr_m=data_numbers;
    neu->tr_n=neu->output_units;

    if(neu->output_data != NULL)
    {
        free(neu->output_data);
    }
    neu->output_data=(double *)calloc((neu->output_units)*data_numbers,sizeof(double));
    neu->out_m=data_numbers;
    neu->out_n=neu->output_units;
    
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void neural_train(struct NEURAL * neu)
{
    int state=0;
    neu->now_steps=0;
    while(state >= 0)
    {
        neu->now_steps=state;
        ////////////////////////////////////////数据分批和采样
        neu->data_sel_func(neu);  //printf("%d\n",1);
        ////////////////////////////////////////向前计算 并 梯度计算
        neu->gradient_calculator(neu);//printf("%d\n",2);
        ////////////////////////////////////////权重矩阵更新
        neu->weight_matrix_updator(neu);///printf("%d\n",3);
        ////////////////////////////////////////学习率更新
        neu->learn_function(neu);//printf("%d\n",4);
        ////////////////////////////////////////收敛检查
        neu->now_steps=neu->exitor(neu);//printf("%d\n",5);
        state=neu->now_steps;
        ////////////////////////////////////////
        if(neu->print_function != NULL)
        {
            neu->print_function(neu);
        }//printf("%d\n",6);
        if(neu->wither_birth_creator!=NULL)
        {
            neu->wither_birth_creator(neu);
        }
        ////////////////////////////////////////
        if(state%1000 == 0)
        {
            srand(time(NULL));
        }
        state++;

    }
}

void neural_run(struct NEURAL * neu)
{
    data_selector_total_init(neu);
    neu->data_sel_func(neu);
    neu->forward_calculator(neu);
}


void neural_free(struct NEURAL * neu)
{
    struct NEURAL_NETWORK * neu_tmp1=neu->neu;
    struct NEURAL_NETWORK * neu_tmp2=neu->neu;
    while(neu_tmp1 != NULL)
    {
        neu_tmp2=neu_tmp1->next;
        neural_network_free(neu_tmp1);
        neu_tmp1=neu_tmp2;
    }
    /////////
    if(neu->input_data != NULL) free(neu->input_data);
    if(neu->output_data != NULL) free(neu->output_data);
    if(neu->train_data != NULL) free(neu->train_data);

    if(neu->data_sel_table != NULL) free(neu->data_sel_table);

    if(neu->para_data_sel_func != NULL) free(neu->para_data_sel_func);
    if(neu->para_exitor != NULL) free(neu->para_exitor);
    if(neu->para_forward_calculator != NULL) free(neu->para_forward_calculator);
    if(neu->para_gradient_calculator != NULL) free(neu->para_gradient_calculator);
    if(neu->para_learn_function != NULL) free(neu->para_learn_function);
    if(neu->para_print_function != NULL) free(neu->para_print_function);
    if(neu->para_weight_matrix_updator != NULL) free(neu->para_weight_matrix_updator);
    
    if(neu->loss_data != NULL) free(neu->loss_data);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//////////////////////////////////////////// data selector
void data_selector_total(struct NEURAL * neu)
{
    for(int i=0;i<neu->in_m;i++)
    {
        neu->data_sel_table[i]=i;
    }
    
}
void data_selector_total_init(struct NEURAL * neu)
{
    /////////////////////
    if(neu->data_sel_table != NULL)
    {
        free(neu->data_sel_table);
    }
    neu->data_sel_table=(int *)calloc(neu->in_m,sizeof(int));
    if(neu->para_data_sel_func != NULL) free(neu->para_data_sel_func);
    neu->para_data_sel_func=NULL;
    neu->data_sel_func=data_selector_total;
    neu->data_cal_number=neu->in_m;
}

void data_selector_rand(struct NEURAL * neu)
{
    int start=rand()%(neu->in_m);
    for(int i=0;i<neu->data_cal_number;i++)
    {
        
        neu->data_sel_table[i]=i + start;
        if(neu->data_sel_table[i] >=neu->in_m) 
        {
            neu->data_sel_table[i]=neu->data_sel_table[i] % (neu->in_m);
        }
    }
}
void data_selector_rand_init(struct NEURAL * neu,double sel_scale)
{
    int sel_num;
    if(sel_scale > 0) sel_num=sel_scale*(neu->in_m);
    else sel_num=(int)fabs(sel_scale);

    
    /////////////////////
    if(neu->data_sel_table != NULL)
    {
        free(neu->data_sel_table);
    }
    neu->data_sel_table=(int *)calloc(sel_num,sizeof(int));
    neu->para_data_sel_func=NULL;
    neu->data_sel_func=data_selector_rand;
    neu->data_cal_number=sel_num;
}

void data_selector_randall(struct NEURAL * neu)
{
    for(int i=0;i<neu->data_cal_number;i++)
    {
        neu->data_sel_table[i]=rand()%(neu->in_m);
    }
}
void data_selector_randall_init(struct NEURAL * neu,double sel_scale)
{
    int sel_num;
    if(sel_scale > 0) sel_num=sel_scale*(neu->in_m);
    else sel_num=(int)fabs(sel_scale);
    /////////////////////
    if(neu->data_sel_table != NULL)
    {
        free(neu->data_sel_table);
    }
    neu->data_sel_table=(int *)calloc(sel_num,sizeof(int));
    neu->para_data_sel_func=NULL;
    neu->data_sel_func=data_selector_randall;
    neu->data_cal_number=sel_num;
}
//////////////////////////////////////////// forward calculator
void forward_calculator_MLP(struct NEURAL * neu)
{
    double * input_data=NULL;
    struct NEURAL_NETWORK * tmp_neu;
    for(int i=0;i<neu->data_cal_number;i++)
    {
        int sel_n=neu->data_sel_table[i];
        input_data=(neu->input_data + sel_n*(neu->in_n) );
        tmp_neu=neu->neu;
        while(tmp_neu!=NULL)
        {
            matrix_A_X_B_add_C(tmp_neu->w,tmp_neu->m,tmp_neu->n,input_data,tmp_neu->n,1,tmp_neu->b,tmp_neu->z);
            matrix_A_ac(tmp_neu->z,tmp_neu->m,1,tmp_neu->a,tmp_neu->activate_function);
            matrix_A_Gac(tmp_neu->z,tmp_neu->m,1,tmp_neu->aa,tmp_neu->G_activate_function);
            input_data=tmp_neu->a;
            tmp_neu=tmp_neu->next; 
        }
        memcpy(neu->output_data + sel_n*(neu->out_n),input_data,neu->output_units * sizeof(double));
    } 

}
void forward_calculator_MLP_init(struct NEURAL * neu)
{
    if(neu->para_forward_calculator != NULL) free(neu->para_forward_calculator);
    neu->para_forward_calculator=NULL;
    neu->forward_calculator=forward_calculator_MLP;
}

void forward_calculator_ASLP(struct NEURAL * neu)
{
    double * input_data=NULL;
    struct NEURAL_NETWORK * tmp_neu;
    struct NEURAL_NETWORK * max_neu;

    tmp_neu=neu->neu;
    while(tmp_neu!=NULL)
    {
        tmp_neu->a[0]=0;
        tmp_neu->aa[0]=0;
        tmp_neu=tmp_neu->next;
    }

    for(int i=0;i<neu->data_cal_number;i++)
    {
        double max_val=-999;
        int sel_n=neu->data_sel_table[i];
        input_data=(neu->input_data + sel_n*(neu->in_n) );
        tmp_neu=neu->neu;
        while(tmp_neu!=NULL)
        {
            matrix_A_X_B_add_C(tmp_neu->w,tmp_neu->m,tmp_neu->n,input_data,tmp_neu->n,1,tmp_neu->b,tmp_neu->z);
            //matrix_A_ac(tmp_neu->z,tmp_neu->m,1,tmp_neu->a,tmp_neu->activate_function);
            //matrix_A_Gac(tmp_neu->z,tmp_neu->m,1,tmp_neu->aa,tmp_neu->G_activate_function);
            if(max_val - tmp_neu->z[0] < 0)
            {
                max_val = tmp_neu->z[0];
                max_neu=tmp_neu;
            }
            tmp_neu=tmp_neu->next; 
        }
        if(max_val - neu->para_gradient_calculator[1] >0)
        {
            max_neu->a[0] += 1.0;
        }
        else
        {
            max_neu->aa[0] += 1.0;
        }
        neu->output_data[sel_n]=max_val;
    } 

}
void forward_calculator_ASLP_init(struct NEURAL * neu)
{
    if(neu->para_forward_calculator != NULL) free(neu->para_forward_calculator);
    neu->para_forward_calculator=NULL;
    neu->forward_calculator=forward_calculator_ASLP;
}

//////////////////////////////////////////// gradient calculator
void gradient_calculator_MLP(struct NEURAL * neu)
{

    double * ref_data=NULL;
    double * input_data=NULL;
    struct NEURAL_NETWORK * tmp_neu=neu->neu;
    struct NEURAL_NETWORK * tmp_neu2;
    double * tmp_sigma=NULL;
    double * tmp_lw=NULL;
    double * tmp_sigma1=NULL;
    while(tmp_neu!=NULL)
    {
        memset(tmp_neu->lb_sigma,0,tmp_neu->m*sizeof(double));
        memset(tmp_neu->lw,0,(tmp_neu->m)*(tmp_neu->n)*sizeof(double));
        tmp_neu=tmp_neu->next;
    }
    for(int i=0;i<neu->data_cal_number;i++)
    {
        int sel_n=neu->data_sel_table[i];
        //forward_calculator
        
        ///////////////////////////////////////////////
        input_data=(neu->input_data + sel_n*(neu->in_n) );
        tmp_neu=neu->neu;
        while(tmp_neu!=NULL)
        {
            tmp_neu2=tmp_neu;
            matrix_A_X_B_add_C(tmp_neu->w,tmp_neu->m,tmp_neu->n,input_data,tmp_neu->n,1,tmp_neu->b,tmp_neu->z);
            matrix_A_ac(tmp_neu->z,tmp_neu->m,1,tmp_neu->a,tmp_neu->activate_function);
            matrix_A_Gac(tmp_neu->z,tmp_neu->m,1,tmp_neu->aa,tmp_neu->G_activate_function);
            input_data=tmp_neu->a;
            tmp_neu=tmp_neu->next; 
        }
        //graden back cal lw and lb
        ///////////////////////////////////////////////
        ref_data=neu->train_data + sel_n*(neu->tr_n);
        
        tmp_sigma=NULL;
        tmp_lw=NULL;
        tmp_sigma1=NULL;

        

        while(tmp_neu2!=NULL)
        {
            
            tmp_lw=(double *)calloc((tmp_neu2->m)*(tmp_neu2->n),sizeof(double));
            tmp_sigma=(double *)calloc((tmp_neu2->m),sizeof(double));

            if(tmp_neu2->next == NULL)// for l_sigma
            {
                neu->G_loss_func(tmp_neu2->a,ref_data, tmp_neu2->m , tmp_sigma, neu->loss_para[0], neu->loss_para[1]);
                for(int lm=0;lm<tmp_neu2->m;lm++) 
                {
                   tmp_sigma[lm]*=tmp_neu2->aa[lm];
                }
            }
            else
            {
                matrix_WT_X_B_adX_C(tmp_neu2->next->w,tmp_neu2->next->m,tmp_neu2->next->n,tmp_sigma1,tmp_neu2->next->m,1,tmp_neu2->aa,tmp_sigma);
            }

            if(tmp_neu2->previous == NULL) // for l_w
            {
                matrix_A_X_BT(tmp_sigma,tmp_neu2->m,1,neu->input_data + sel_n*(neu->in_n),neu->input_units,1,tmp_lw);
            }
            else
            {
                matrix_A_X_BT(tmp_sigma,tmp_neu2->m,1,tmp_neu2->previous->a,tmp_neu2->previous->m,1,tmp_lw);
            }


            for(int lm=0;lm<tmp_neu2->m;lm++)
            {
                tmp_neu2->lb_sigma[lm]+=tmp_sigma[lm];
            }

            for(int i=0;i<(tmp_neu2->m)*(tmp_neu2->n);i++)
            {
                tmp_neu2->lw[i]+=tmp_lw[i];
            }

            free(tmp_sigma1);
            tmp_sigma1=tmp_sigma;

           
            free(tmp_lw);
            tmp_neu2=tmp_neu2->previous;
            
        }
       
        free(tmp_sigma1);
    }
    tmp_neu=neu->neu;
      
    while(tmp_neu!=NULL)
    {
        double sum_w=0;
        double sum_b=0;
        for(int i=0;i<tmp_neu->m*tmp_neu->n;i++)
        {
            sum_w+=pow(tmp_neu->lw[i],2);
        }
        sum_w=sqrt(sum_w);
        if(sum_w - neu->para_gradient_calculator[0] >0)
        {
            for(int i=0;i<tmp_neu->m*tmp_neu->n;i++)
            {
                tmp_neu->lw[i] *=(neu->para_gradient_calculator[0])/sum_w;
            }            
        }
        
        for(int i=0;i<tmp_neu->m;i++)
        {
            sum_b+=pow(tmp_neu->lb_sigma[i],2);
        }
        if(sum_b - neu->para_gradient_calculator[1] >0)
        {
            for(int i=0;i<tmp_neu->m;i++)
            {
                tmp_neu->lb_sigma[i] *=(neu->para_gradient_calculator[1])/sum_b;
            }
        }
        tmp_neu=tmp_neu->next;
    }
        

}
void gradient_calculator_MLP_init(struct NEURAL *  neu, double gradient_w_clip,double gradient_b_clip)
{
    neu->para_gradient_calculator=(double *)calloc(2,sizeof(double));
    neu->para_gradient_calculator[0]=gradient_w_clip;
    neu->para_gradient_calculator[1]=gradient_b_clip;
    neu->gradient_calculator=gradient_calculator_MLP;
}

void gradient_calculator_MLP_GaussSmear(struct NEURAL * neu)
{

    double * ref_data=NULL;
    double * input_data=NULL;
    struct NEURAL_NETWORK * tmp_neu=neu->neu;
    struct NEURAL_NETWORK * tmp_neu2;
    double * tmp_sigma=NULL;
    double * tmp_lw=NULL;
    double * tmp_sigma1=NULL;
    //////////////////////////////
    double nor=0;
    double * tmp_vec=(double *)calloc(neu->input_units,sizeof(double));

   
    /////////////////////////////
    while(tmp_neu!=NULL)
    {
        memset(tmp_neu->lb_sigma,0,tmp_neu->m*sizeof(double));
        memset(tmp_neu->lw,0,(tmp_neu->m)*(tmp_neu->n)*sizeof(double));
        tmp_neu=tmp_neu->next;
    }
    for(int i=0;i<neu->data_cal_number;i++)
    {
        int sel_n=neu->data_sel_table[i];
        //forward_calculator
        
        ///////////////////////////////////////////////Gauss_smaer_data
        input_data=(neu->input_data + sel_n*(neu->in_n) );

        nor=box_muller(0,neu->para_gradient_calculator[2]);
        double sums=0;
        for(int sv=0;sv<neu->input_units;sv++)
        {
            tmp_vec[sv]=box_muller(0,1);
            sums+=pow(tmp_vec[sv],2);
        }
        sums=pow(sums,0.5);
        for(int sv=0;sv<neu->input_units;sv++)
        {
            tmp_vec[sv]=tmp_vec[sv]*nor/sums + input_data[sv];
        } 
        ///////////////////////////////////////////////
        input_data=tmp_vec;

        tmp_neu=neu->neu;
        while(tmp_neu!=NULL)
        {
            tmp_neu2=tmp_neu;
            matrix_A_X_B_add_C(tmp_neu->w,tmp_neu->m,tmp_neu->n,input_data,tmp_neu->n,1,tmp_neu->b,tmp_neu->z);
            matrix_A_ac(tmp_neu->z,tmp_neu->m,1,tmp_neu->a,tmp_neu->activate_function);
            matrix_A_Gac(tmp_neu->z,tmp_neu->m,1,tmp_neu->aa,tmp_neu->G_activate_function);
            input_data=tmp_neu->a;
            tmp_neu=tmp_neu->next; 
        }
        //graden back cal lw and lb
        ///////////////////////////////////////////////
        ref_data=neu->train_data + sel_n*(neu->tr_n);
        
        tmp_sigma=NULL;
        tmp_lw=NULL;
        tmp_sigma1=NULL;

        

        while(tmp_neu2!=NULL)
        {
            
            tmp_lw=(double *)calloc((tmp_neu2->m)*(tmp_neu2->n),sizeof(double));
            tmp_sigma=(double *)calloc((tmp_neu2->m),sizeof(double));

            if(tmp_neu2->next == NULL)// for l_sigma
            {
                neu->G_loss_func(tmp_neu2->a,ref_data, tmp_neu2->m , tmp_sigma, neu->loss_para[0], neu->loss_para[1]);
                for(int lm=0;lm<tmp_neu2->m;lm++) 
                {
                   tmp_sigma[lm]*=tmp_neu2->aa[lm];
                }
            }
            else
            {
                matrix_WT_X_B_adX_C(tmp_neu2->next->w,tmp_neu2->next->m,tmp_neu2->next->n,tmp_sigma1,tmp_neu2->next->m,1,tmp_neu2->aa,tmp_sigma);
            }

            if(tmp_neu2->previous == NULL) // for l_w
            {
                matrix_A_X_BT(tmp_sigma,tmp_neu2->m,1,neu->input_data + sel_n*(neu->in_n),neu->input_units,1,tmp_lw);
            }
            else
            {
                matrix_A_X_BT(tmp_sigma,tmp_neu2->m,1,tmp_neu2->previous->a,tmp_neu2->previous->m,1,tmp_lw);
            }


            for(int lm=0;lm<tmp_neu2->m;lm++)
            {
                tmp_neu2->lb_sigma[lm]+=tmp_sigma[lm];
            }

            for(int i=0;i<(tmp_neu2->m)*(tmp_neu2->n);i++)
            {
                tmp_neu2->lw[i]+=tmp_lw[i];
            }

            free(tmp_sigma1);
            tmp_sigma1=tmp_sigma;

           
            free(tmp_lw);
            tmp_neu2=tmp_neu2->previous;
            
        }
       
        free(tmp_sigma1);
    }
    tmp_neu=neu->neu;
    
    while(tmp_neu!=NULL)
    {
        double sum_w=0;
        double sum_b=0;
        for(int i=0;i<tmp_neu->m*tmp_neu->n;i++)
        {
            sum_w+=pow(tmp_neu->lw[i],2);
        }
        sum_w=sqrt(sum_w);
        if(sum_w - neu->para_gradient_calculator[0] >0)
        {
            for(int i=0;i<tmp_neu->m*tmp_neu->n;i++)
            {
                tmp_neu->lw[i] *=(neu->para_gradient_calculator[0])/sum_w;
            }            
        }
        
        for(int i=0;i<tmp_neu->m;i++)
        {
            sum_b+=pow(tmp_neu->lb_sigma[i],2);
        }
        if(sum_b - neu->para_gradient_calculator[1] >0)
        {
            for(int i=0;i<tmp_neu->m;i++)
            {
                tmp_neu->lb_sigma[i] *=(neu->para_gradient_calculator[1])/sum_b;
            }
        }
        tmp_neu=tmp_neu->next;
    }
        
    free(tmp_vec);
}
void gradient_calculator_MLP_GaussSmear_init(struct NEURAL *  neu, double gradient_w_clip,double gradient_b_clip,double sigma)
{
    neu->para_gradient_calculator=(double *)calloc(3,sizeof(double));
    neu->para_gradient_calculator[0]=gradient_w_clip;
    neu->para_gradient_calculator[1]=gradient_b_clip;
    neu->para_gradient_calculator[2]=sigma;
    neu->gradient_calculator=gradient_calculator_MLP;
}


void gradient_calculator_ASLP(struct NEURAL * neu)
{
    //0 w clip
    //1 similar cutoff
    //2 cluster parameters
    double * ref_data=NULL;
    double * input_data=NULL;
    struct NEURAL_NETWORK * tmp_neu=neu->neu;
    double * tmp_sigma=NULL;
    double * tmp_lw=NULL;

    while(tmp_neu!=NULL)
    {
        memset(tmp_neu->lb_sigma,0,tmp_neu->m*sizeof(double));
        memset(tmp_neu->lw,0,(tmp_neu->m)*(tmp_neu->n)*sizeof(double));
        memset(tmp_neu->a,0,tmp_neu->m*sizeof(double));
        
        //////////////////////////////
        double sum_w=0;
        for(int i=0;i<tmp_neu->n;i++)
        {
            sum_w+=pow(tmp_neu->w[i],2);
        }
        sum_w=sqrt(sum_w);
        if(sum_w >0.00001)
        {
            for(int i=0;i<tmp_neu->m*tmp_neu->n;i++)
            {
                tmp_neu->w[i] /=sum_w;
            }  
        }      
        //////////////////////////////
        tmp_neu=tmp_neu->next;
    }
    tmp_neu=neu->neu;
    struct NEURAL_NETWORK * max_neural;
    double max_val;   
    tmp_lw=(double *)calloc((tmp_neu->m)*(tmp_neu->n),sizeof(double));
    tmp_sigma=(double *)calloc(1,sizeof(double));
    for(int i=0;i<neu->data_cal_number;i++)
    {
        int sel_n=neu->data_sel_table[i];
        //forward_calculator
        max_val=-99999;
        ///////////////////////////////////////////////
        input_data=(neu->input_data + sel_n*(neu->in_n) );
        tmp_neu=neu->neu;
        while(tmp_neu!=NULL)
        {
            matrix_A_X_B_add_C(tmp_neu->w,tmp_neu->m,tmp_neu->n,input_data,tmp_neu->n,1,tmp_neu->b,tmp_neu->z);
            //matrix_A_ac(tmp_neu->z,tmp_neu->m,1,tmp_neu->a,tmp_neu->activate_function);
            //matrix_A_Gac(tmp_neu->z,tmp_neu->m,1,tmp_neu->aa,tmp_neu->G_activate_function);
            if(max_val - tmp_neu->z[0] < 0)
            {
                max_val = tmp_neu->z[0];
                max_neural=tmp_neu;
            }
            tmp_neu=tmp_neu->next; 
        }
        if(max_val - neu->para_gradient_calculator[1] >0)
        {
            max_neural->a[0]+=1.0;
        }

        //graden back cal lw and lb
        ///////////////////////////////////////////////
        ref_data=neu->train_data + sel_n*(neu->tr_n);
        
        tmp_neu=neu->neu;
        
        while(tmp_neu!=NULL)
        {
            neu->G_loss_func(tmp_neu->z,ref_data,1,tmp_sigma,neu->loss_para[0],neu->loss_data[1]);
            matrix_A_X_BT(tmp_sigma,1,1,neu->input_data + sel_n*(neu->in_n),neu->input_units,1,tmp_lw);
            for(int nw=0;nw<tmp_neu->n;nw++)
            {
                tmp_neu->lw[nw]+=tmp_lw[nw];
            }
            tmp_neu=tmp_neu->next;
        }
        

    
    }
    free(tmp_lw);
    free(tmp_sigma);
    



    //gradient clip
    tmp_neu=neu->neu;
    while(tmp_neu!=NULL)
    {
        double sum_w=0;
        double clu_coef=exp(-1.0*pow(tmp_neu->a[0],2)/pow(neu->para_gradient_calculator[2],2)) ; 
        for(int i=0;i<tmp_neu->m*tmp_neu->n;i++)
        {
            sum_w+=pow(tmp_neu->lw[i],2);
        }
        sum_w=sqrt(sum_w);
        if(sum_w - neu->para_gradient_calculator[0] >0)
        {
            for(int i=0;i<tmp_neu->m*tmp_neu->n;i++)
            {
                tmp_neu->lw[i] *=(neu->para_gradient_calculator[0]  * clu_coef  )/sum_w;
            }            
        }
        else
        {
            for(int i=0;i<tmp_neu->m*tmp_neu->n;i++)
            {
                tmp_neu->lw[i] *= clu_coef ;
            }            
        }        
        
        tmp_neu=tmp_neu->next;
    }
}
void gradient_calculator_ASLP_init(struct NEURAL *  neu, double gradient_w_clip,double similarity_cutoff,double clustering_coeff)
{
    neu->para_gradient_calculator=(double *)calloc(3,sizeof(double));
    neu->para_gradient_calculator[0]=gradient_w_clip;
    neu->para_gradient_calculator[1]=similarity_cutoff;
    neu->para_gradient_calculator[2]=clustering_coeff;
    neu->gradient_calculator=gradient_calculator_ASLP;
}

//////////////////////////////////////////// weight matrix updator
void weight_matrix_updator_Adam(struct NEURAL * neu)
{



    struct NEURAL_NETWORK * tmp_neu=neu->neu;

    ////////////////////////////////////////////////////////////////////////////
    tmp_neu=neu->neu;
    double study_len=neu->learn_ratio;
    while(tmp_neu!= NULL)
    {
        //////////////////////////////////////////
        for(int i=0;i<(tmp_neu->m )*(tmp_neu->n);i++)
        {
            tmp_neu->mt[i] = (neu->para_weight_matrix_updator[0]) * (tmp_neu->mt[i]) + (1 - neu->para_weight_matrix_updator[0])*(tmp_neu->lw[i]);
            tmp_neu->vt[i] = (neu->para_weight_matrix_updator[1]) * (tmp_neu->vt[i]) + (1 - neu->para_weight_matrix_updator[1])*pow((tmp_neu->lw[i]),2);
            tmp_neu->w[i]=  tmp_neu->w[i] -   study_len*(  (tmp_neu->mt[i]/(1+  pow(neu->para_weight_matrix_updator[0],neu->now_steps + 1)   ))/(0.00001 +  sqrt(  tmp_neu->vt[i]/(1-  pow(neu->para_weight_matrix_updator[1],neu->now_steps + 1)   )    )     )    );
        
        }                      
        /////////////////////////////////////////
        for(int i=0;i<(tmp_neu->m );i++)
        {
            tmp_neu->mbt[i] = (neu->para_weight_matrix_updator[0]) * (tmp_neu->mbt[i]) + (1 - neu->para_weight_matrix_updator[0])*(tmp_neu->lb_sigma[i]);
            tmp_neu->vbt[i] = (neu->para_weight_matrix_updator[1]) * (tmp_neu->vbt[i]) + (1 - neu->para_weight_matrix_updator[1])*pow((tmp_neu->lb_sigma[i]),2);
            tmp_neu->b[i]=  tmp_neu->b[i] -   study_len*(  (tmp_neu->mbt[i]/(1+  pow(neu->para_weight_matrix_updator[0],neu->now_steps + 1)   ))/(0.00001 +  sqrt(  tmp_neu->vbt[i]/(1-  pow(neu->para_weight_matrix_updator[1],neu->now_steps + 1)   )    )     )    );
        }
        /////////////////////////////////////////
        tmp_neu=tmp_neu->next;
    }
}
void weight_matrix_updator_Adam_init(struct NEURAL * neu,double beta1,double beta2)
{
    neu->para_weight_matrix_updator=(double *)calloc(2,sizeof(double));
    neu->para_weight_matrix_updator[0]=beta1;
    neu->para_weight_matrix_updator[1]=beta2;
    neu->weight_matrix_updator=weight_matrix_updator_Adam;
}


void weight_matrix_updator_GradientDown(struct NEURAL * neu)
{



    struct NEURAL_NETWORK * tmp_neu=neu->neu;

    ////////////////////////////////////////////////////////////////////////////

    tmp_neu=neu->neu;

    double study_len=neu->learn_ratio;
    while(tmp_neu!= NULL)
    {
        //////////////////////////////////////////
        for(int i=0;i<(tmp_neu->m )*(tmp_neu->n);i++)
        {
            tmp_neu->w[i]-=study_len*tmp_neu->lw[i];
        
        }                      
        /////////////////////////////////////////
        for(int i=0;i<(tmp_neu->m );i++)
        {
            tmp_neu->b[i]-=study_len*tmp_neu->lb_sigma[i];
        }
        /////////////////////////////////////////
        tmp_neu=tmp_neu->next;
    }
}
void weight_matrix_updator_GradientDown_init(struct NEURAL * neu)
{
    if(neu->para_weight_matrix_updator!=NULL)
    {
        free(neu->para_weight_matrix_updator);
        neu->para_weight_matrix_updator=NULL;
    }
    neu->weight_matrix_updator=weight_matrix_updator_GradientDown;
}

//////////////////////////////////////////// learn functions
void learn_ratio_stage(struct NEURAL * neu)
{    
    if(  neu->now_steps %(  (int)(neu->para_learn_function[2])   )   == 0  )
    {
        
        if(neu->learn_ratio - neu->para_learn_function[1] <=0)
        {
             neu->learn_ratio = neu->para_learn_function[1];
        }
        else
        {
            neu->learn_ratio= (neu->para_learn_function[0]) * pow(neu->para_learn_function[3] , (int)(neu->now_steps /((neu->para_learn_function[2])   ))     )   ;
        }
    }

}
void learn_ratio_stage_init(struct NEURAL * neu,double  init_lr,double end_lr,int flush_steps,double ratio)
{
    neu->para_learn_function=(double *)calloc(4,sizeof(double));
    neu->para_learn_function[0]=init_lr;
    neu->para_learn_function[1]=end_lr;
    neu->para_learn_function[2]=1.0*flush_steps;
    neu->para_learn_function[3]=ratio;
    neu->learn_ratio=init_lr;
    neu->learn_function=learn_ratio_stage;
}

void learn_ratio_exp(struct NEURAL * neu)
{    

    neu->learn_ratio=(neu->para_learn_function[1]) + (neu->para_learn_function[0] - neu->para_learn_function[1])*exp(- (neu->para_learn_function[2]) * pow(neu->now_steps,2)   );

}
void learn_ratio_exp_init(struct NEURAL * neu,double  init_lr,double end_lr,int half_steps)
{
    neu->para_learn_function=(double *)calloc(3,sizeof(double));
    neu->para_learn_function[0]=init_lr;
    neu->para_learn_function[1]=end_lr;
    neu->para_learn_function[2]=-log(0.5)/pow(half_steps,2);
    neu->learn_ratio=init_lr;
    neu->learn_function=learn_ratio_exp;
}

//////////////////////////////////////////// exitor
int exitor_max(struct NEURAL * neu)
{
    
    int state=neu->now_steps;
    if(neu->now_steps %  ( (int)(neu->para_exitor[0])     ) == 0)
    {
        
        double * input;
        double * refer;
        //////////////////////////
        int * tmp_sel_table=neu->data_sel_table;
        int cal_num=neu->data_cal_number;
        int * mew_data_sel_table=(int *)calloc(neu->in_m,sizeof(int));
        for(int i=0;i<neu->in_m;i++)
        {
            mew_data_sel_table[i]=i;
        }
       
        neu->data_sel_table=mew_data_sel_table;
        neu->data_cal_number=neu->in_m;
        neu->forward_calculator(neu); 
        ///////////////////////////
        double avg=-99999999;
        for(int i=0;i<neu->in_m;i++)
        {
            input=neu->output_data + i*(neu->out_n);
            refer=neu->train_data + i*(neu->tr_n);
            
            neu->loss_data[i]=neu->loss_func(input,refer,neu->tr_n,neu->loss_para[0],neu->loss_para[1]);

            if(neu->loss_data[i] - avg >0) avg=neu->loss_data[i];

            
        }
        if(avg- neu->para_exitor[1] < 0)
        {
            state=-1000;
        }
        free(mew_data_sel_table);
        neu->data_sel_table=tmp_sel_table;
        neu->data_cal_number=cal_num;
        
    }

    if(neu->now_steps - neu->para_exitor[2]>0)
    {
        state=-1000;
    }

    return state;
    
}
void exitor_max_init(struct NEURAL * neu,int flush_steps,double cutoff_prec,int max_train_steps)
{
    neu->para_exitor=(double *)calloc(3,sizeof(double));
    neu->para_exitor[0]=1.0*flush_steps;
    neu->para_exitor[1]=cutoff_prec;   
    neu->para_exitor[2]=1.0*max_train_steps; 
    if(neu->loss_data != NULL) free(neu->loss_data);
    neu->loss_data=(double *)calloc(neu->in_m,sizeof(double));
    for(int i=0;i<neu->in_m;i++)
    {
        neu->loss_data[i]=9999999;
    }
    neu->exitor=exitor_max;
}


int exitor_avg(struct NEURAL * neu)
{
    
    int state=neu->now_steps;
    if(neu->now_steps %  ( (int)(neu->para_exitor[0])     ) == 0)
    {
        
        double * input;
        double * refer;
        //////////////////////////
        int * tmp_sel_table=neu->data_sel_table;
        int cal_num=neu->data_cal_number;
        int * mew_data_sel_table=(int *)calloc(neu->in_m,sizeof(int));
        for(int i=0;i<neu->in_m;i++)
        {
            mew_data_sel_table[i]=i;
        }
       
        neu->data_sel_table=mew_data_sel_table;
        neu->data_cal_number=neu->in_m;
        neu->forward_calculator(neu); 
        ///////////////////////////
        double avg=0;
        for(int i=0;i<neu->in_m;i++)
        {
            input=neu->output_data + i*(neu->out_n);
            refer=neu->train_data + i*(neu->tr_n);
            
            neu->loss_data[i]=neu->loss_func(input,refer,neu->tr_n,neu->loss_para[0],neu->loss_para[1]);

            avg+=neu->loss_data[i];
        }
        avg=avg/(neu->in_m);
        if(avg- neu->para_exitor[1] < 0)
        {
            state=-1000;
        }
        free(mew_data_sel_table);
        neu->data_sel_table=tmp_sel_table;
        neu->data_cal_number=cal_num;
        
    }

    if(neu->now_steps - neu->para_exitor[2]>0)
    {
        state=-1000;
    }

    return state;
    
}
void exitor_avg_init(struct NEURAL * neu,int flush_steps,double cutoff_prec,int max_train_steps)
{
    neu->para_exitor=(double *)calloc(3,sizeof(double));
    neu->para_exitor[0]=1.0*flush_steps;
    neu->para_exitor[1]=cutoff_prec;   
    neu->para_exitor[2]=1.0*max_train_steps; 
    if(neu->loss_data != NULL) free(neu->loss_data);
    neu->loss_data=(double *)calloc(neu->in_m,sizeof(double));
    for(int i=0;i<neu->in_m;i++)
    {
        neu->loss_data[i]=9999999;
    }
    neu->exitor=exitor_avg;
}


//////////////////////////////////////////// printor
void print_function(struct NEURAL * neu)
{
    if( (int)(neu->para_print_function[0] )  == 1 && neu->now_steps % ( (int)(neu->para_print_function[1])  ) == 0) // for loss function
    {
        double sum=0;
        for(int i=0;i<neu->in_m;i++)
        {
            sum+=neu->loss_data[i];
        }
        sum/=1.0*neu->in_m;
        printf("steps: %d\tloss: %.6f\tlr: %.8f\n",neu->now_steps,sum,neu->learn_ratio);
    }
    else if( (int)(neu->para_print_function[0] )  == 2 && neu->now_steps % ( (int)(neu->para_print_function[1])  ) == 0) // for matrix w and b
    {
        struct NEURAL_NETWORK * tmp_neu=neu->neu;
        printf("%d\t",neu->now_steps);
        while(tmp_neu != NULL)
        {
            double sumw=0;
            double sumb=0;
            for(int i=0;i<tmp_neu->m * tmp_neu->n;i++)
            {
                sumw+=pow(tmp_neu->w[i],2);
            }
            sumw=sqrt(sumw);
            for(int i=0;i<tmp_neu->m ;i++)
            {
                sumb+=pow(tmp_neu->b[i],2);
            }
            sumb=sqrt(sumb);
            printf("%.4f:%.4f\t",sumw,sumb);
            tmp_neu=tmp_neu->next;
        }
        printf("\n");
    }
    else if( (int)(neu->para_print_function[0] )  == 3 && neu->now_steps % ( (int)(neu->para_print_function[1])  ) == 0) // for matrix lw and lb
    {
        struct NEURAL_NETWORK * tmp_neu=neu->neu;
        printf("%d\t",neu->now_steps);
        while(tmp_neu != NULL)
        {
            double sumw=0;
            double sumb=0;
            for(int i=0;i<tmp_neu->m * tmp_neu->n;i++)
            {
                sumw+=pow(tmp_neu->lw[i],2);
            }
            sumw=sqrt(sumw);
            for(int i=0;i<tmp_neu->m ;i++)
            {
                sumb+=pow(tmp_neu->lb_sigma[i],2);
            }
            sumb=sqrt(sumb);
            printf("%.4f:%.4f\t",sumw,sumb);
            tmp_neu=tmp_neu->next;
        }
        printf("\n");
    }
    if( (int)(neu->para_print_function[0] )  == 4 && neu->now_steps % ( (int)(neu->para_print_function[1])  ) == 0) // for loss function and number layers for ASLP
    {
        double sum=0;
        for(int i=0;i<neu->in_m;i++)
        {
            sum+=neu->loss_data[i];
        }
        sum/=1.0*neu->in_m;
        printf("steps: %d\tloss: %.6f\tlr: %.8f layers:%d\n",neu->now_steps,sum,neu->learn_ratio,neu->layers);
    }
}
void print_function_init(struct NEURAL * neu,int print_type,int flush_steps)
{
    neu->para_print_function=(double *)calloc(2,sizeof(double));
    neu->para_print_function[0]=1.0*print_type;
    neu->para_print_function[1]=1.0*flush_steps;

    neu->print_function=print_function;
}
////////////////////////////////////////////wither and birth creator

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
////////////////////////////

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//mathematica
//the activate functions
double activate_function_default(double x)//无
{
    return x;
}
double activate_function_leaky_ReLU(double x)//拟合分段线性函数
{
    if(x<0) return alpha_LRU*x;
    else return x;
}
double activate_function_sigmoid(double x)// 拟合，S型，渐近，压缩
{
    return 1.0/(1.0+exp(x));
}
double activate_function_Tanh(double x) //拟合光滑有界，连续函数
{
    return tanh(x);
}
double activate_function_softplus(double x) // 和ReLu相似，但二阶连续可导
{
    return log(1+exp(x));
}
double activate_function_SINE(double x) //周期振荡函数
{
    return sin(x);
}
double activate_function_swish(double x)// 无限可微，光滑，复杂函数
{
    return x/(1 + exp(-x));
}


double G_activate_function_default(double x)
{
    return 1;
}
double G_activate_function_leaky_ReLU(double x)
{
    if(x >0) return 1;
    else return alpha_LRU;
}
double G_activate_function_sigmoid(double x)
{
    return 1.0/(1.0 + exp(x)) * ( 1- 1.0/(1.0+exp(x))    );
}
double G_activate_function_Tanh(double x) 
{
    return 1-pow(tanh(x),2);
}
double G_activate_function_softplus(double x) 
{
    return 1.0/(1.0+exp(-x));
}
double G_activate_function_SINE(double x) 
{
    return cos(x);
}
double G_activate_function_swish(double x)
{
    return 1/(1 + exp(-x))  + x*(1/(1 + exp(-x)))*( 1- 1/(1 + exp(-x))   ) ;
}


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//loss functions
double loss_func_MSE(double * yi,double * yit, int numbers, double para1,double para2)
{
    double loss=0;
    for(int i=0;i<numbers;i++)
    {
        loss+=pow(yi[i] - yit[i], 2);
    }
    return 1.0/(2.0*numbers)    *loss;
}
void G_loss_func_MSE( double * yi , double * yit , int numbers, double * output, double para1,double para2)
{
    for(int i=0;i<numbers;i++)
    {
        output[i]= 1.0/numbers*(  yi[i] - yit[i]    );
    }

}


double loss_func_Huber(double * yi,double * yit, int numbers, double para1,double para2)
{
    double loss=0;
    for(int i=0;i<numbers;i++)
    {
        if(  fabs(yi[i] - yit[i])  - para1 <=0  ) loss+=0.5*pow(yi[i] - yit[i], 2);
        else loss+=para1 * (  fabs(yi[i] - yit[i]) - 0.5*para1     );
        
    }
    return loss;    
}
void G_loss_func_Huber( double * yi , double * yit , int numbers, double * output, double para1,double para2)
{
    for(int i=0;i<numbers;i++)
    {
        if(  fabs(yi[i] - yit[i])  - para1 <=0  ) output[i]= (  yi[i] - yit[i]    );
        else
        {
            if(yi[i] - yit[i] > 0)
            {
                output[i]=para1;
            }
            else if( yi[i] - yit[i] < 0)
            {
                output[i]=-para1;
            }
            else
            {   
                output[i]=0;
            }
        }
        
    }    
}


double loss_func_BCEwithLogitsLoss(double * yi,double * yit, int numbers, double para1,double para2)
{
    double loss=0;
    for(int i=0;i<numbers;i++)
    {
        if(yi[i] <= 0)
        {
            loss+=0-yi[i]*yit[i] + log10(1+exp(-fabs(yi[i])));
        }
        else
        {
            loss+=yi[i]-yi[i]*yit[i] + log10(1+exp(-fabs(yi[i])));
        }
    }    
    return loss;
}
void G_loss_func_BCEwithLogitsLoss( double * yi , double * yit , int numbers, double * output, double para1,double para2)
{
    for(int i=0;i<numbers;i++)
    {
        output[i]= yi[i]-yit[i];
    }   
}














