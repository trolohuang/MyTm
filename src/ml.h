#ifndef ML_H
#define ML_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <stdarg.h>
#define alpha_LRU 0.001   //激活函数参数 leaky ReLU 参数alpha

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//MLP build process
//
//struct NEURAL * neu=neural_init_MLP()
//neural_set_data() //set train data befor data_selector
//neural_set_activate_function
//neural_set_loss_func
//neural_init_w_b
//data_selector
//forward_calculator_MLP_init
//gradient_calculator_MLP_init
//weight_matrix_updator_Adam_init
//learn_ratio_stage_init
//exitor_avg_init
//print_function_init //optional
//neural_train






//神经网络结构
struct NEURAL_NETWORK // a=activate_func(z); z=w*a_1 + b ; aa=G_activate_func(z)
{
    int m;         //权重矩阵 w： mxn 
    int n;

    double alpha; //神经元幅值
    double lalpha; //神经元幅值梯度
    double mat;   //动量一阶矩
    double vat;   //动量二阶矩

    double * w;   //机器学习 权重矩阵
    double * lw;  //机器学习 权重矩阵梯度
    double * mt; // for adam：权重矩阵w的动量一阶矩
    double * vt; // for adam：权重矩阵w的动量二阶矩


    double * b;   //机器学习      权重增幅
    double * lb_sigma; //机器学习 权重增幅梯度
    double * mbt;// for adam：权重矩阵b的动量一阶矩
    double * vbt;// for adam：权重矩阵b的动量二阶矩    

    double * z;   //  z=w*a_1 + b       机器学习 临时输出值
    double * a; //a=f(z)    临时输出值的激活函数值
    double * aa; //aa=f'(z) 临时输出值的激活函数梯度值 


    double (*activate_function)(double); // 激活函数
    double (*G_activate_function)(double);//激活函数梯度

    struct NEURAL_NETWORK * next; 
    struct NEURAL_NETWORK * previous;
};
struct NEURAL
{

    //neural struct 
    int layers;              
    int input_units;
    int output_units;
    struct NEURAL_NETWORK * neu;

    //data 
    int in_m;              
    int in_n;
    double * input_data;   // 输入数据， in_m x in_n ， 数据格式，1,2,3,4\n 5,6,7,8\n 9 10 ...
    int out_m;
    int out_n;
    double * output_data;  // 输出数据 out_m x out_n   
    int tr_m;
    int tr_n;
    double * train_data;   // 


    ////////////////////////////////////
    int now_steps;


    // data selector
    int data_cal_number;
    int * data_sel_table;
    void (* data_sel_func)(struct NEURAL * neu);
    double * para_data_sel_func;
    // calculator
    void (* forward_calculator)(struct NEURAL * neu);
    double * para_forward_calculator;
    // gradient calculator
    void (* gradient_calculator)(struct NEURAL * neu);
    double * para_gradient_calculator;
    //weight updating 
    void (* weight_matrix_updator)(struct NEURAL * neu);
    double * para_weight_matrix_updator;
    //learn function updator
    double learn_ratio;
    void (* learn_function)(struct NEURAL * neu);
    double *para_learn_function;
    //exit checker
    int (* exitor)(struct NEURAL * neu);
    double * para_exitor;
    // printor
    void (* print_function)(struct NEURAL * neu);
    double * para_print_function;
    //neural dead and birth
    void (* wither_birth_creator)(struct NEURAL * neu);
    double * para_wither_birth_creator;
    /////////////////////////////////////

    double loss_para[10];    //损失函数参数
    double * loss_data;
    double (*loss_func)(double * , double *, int , double ,double); // 损失函数
    void (* G_loss_func)(double * , double * , int ,double * ,double ,double); //损失函数对a导数


};

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//调用函数



//设置激活函数
//参数：神经网络，激活函数，激活函数导数，设置层
//set_layers 设置第几层，当set_layers = 0 时设置全部层；set_layers > 0 设置对应层的激活函数 1是指第0层；set_layers < 0 从后往前第几层，-1是指最后一层
void neural_set_activate_function(struct NEURAL * neu,double (*activate_function)(double),double (*G_activate_function)(double),int set_layers);
//设置损失函数以及损失函数对a求导数的函数
void neural_set_loss_func(struct NEURAL * neu,  double (*loss_func)(double * , double *, int , double ,double) , void (* G_loss_func)(double * , double * , int ,double * ,double ,double));

void neural_set_data(struct NEURAL * neu,double * input_data, double * reference_data,int data_numbers);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//初始化一个神经网络
//参数：输出神经网络结构，输入层维度，输出层维度，隐藏层层数，第一隐藏层神经元数，第二第一隐藏层神经元数，...
struct NEURAL * neural_init_MLP(int input_dimen,int output_dimen,int layers,...);
struct NEURAL * neural_init_ASLP(int input_dimen,int layers);
//////////////////////////////////////////////////////////////////////////////////
//初始化权重矩阵，训练开始前必须调用
//参数：神经网络结构，初始化类型
//初始化类型：type=0 He/Kaiming 初始化 ； type=1 Xavier初始化
void neural_set_init_w_b_MLP(struct NEURAL * neu,int type);
void neural_set_init_w_b_const(struct NEURAL * neu,double val);
//////////////////////////////////////////////////////////////////////////////////
//数据选择器，控制每次迭代的数据分批
//数据选择器，迭代时，分配data_sel_table内存，存储数据编号，在向前计算中依次计算其中对应的数据；设置data_cal_number，代表当前计算的全部数据个数
//
//
//全部选择，选择所有数据，通常用于训练后直接计算
void data_selector_total_init(struct NEURAL * neu);
//随机初始选择，随机选择数据起始点，然后计算（sel_scale*全部数据）个数，当sel_scale<0，每次计算abs(sel_scale)个数据
//参数：神经网络,选择比例
//sel_scale<0，每次计算abs(sel_scale)个数据；当sel_scale>0每次计算（sel_scale*全部数据）个数据
void data_selector_rand_init(struct NEURAL * neu,double sel_scale);
//
//
void data_selector_randall_init(struct NEURAL * neu,double sel_scale);
//////////////////////////////////////////////////////////////////////////////////
//网络向前计算器，将数据选择器中的所有数据，根据网络类型向前计算的出的结果放到output_data中
//
//
//MLP网络向前计算其
void forward_calculator_MLP_init(struct NEURAL * neu);

void forward_calculator_ASLP_init(struct NEURAL * neu);
//////////////////////////////////////////////////////////////////////////////////
//梯度计算器，根据向前计算器的计算结果和网络类型计算权重矩阵w和偏执矩阵b的梯度 对应 NEURAl_STRUCT 中的lw和lb_sigma
//
//
//MLP网络的梯度计算
//参数：神经网络：w梯度裁减值，b梯度裁减值
void gradient_calculator_MLP_init(struct NEURAL *  neu, double gradient_w_clip,double gradient_b_clip);

//
//MLP网络的梯度计算, 带输入数据作gauss smear
//参数：神经网络：w梯度裁减值，b梯度裁减值, Gauss smear sigma
//对输入描述符作一个高斯smear以增强对噪声的处理
void gradient_calculator_MLP_GaussSmear_init(struct NEURAL *  neu, double gradient_w_clip,double gradient_b_clip,double sigma);

void gradient_calculator_ASLP_init(struct NEURAL *  neu, double gradient_w_clip,double similarity_cutoff,double clustering_coeff);
//////////////////////////////////////////////////////////////////////////////////
//权重矩阵更新器，控制如何通过梯度矩阵 更新 权重矩阵
//
//
//Adam权重更新器
///参数：神经网络，beta1,beta2
void weight_matrix_updator_Adam_init(struct NEURAL * neu,double beta1,double beta2);


//梯度下降法更新梯度
void weight_matrix_updator_GradientDown_init(struct NEURAL * neu);
//////////////////////////////////////////////////////////////////////////////////
//学习率更新器，控制学习率 neu->learn_ratio 如何跟随学习步数neu->steps变化
//
//
//台阶更新，每隔flush步数，学习率从初始init_lr变为 init_lr * pow(ratio,n),n是更新了多少次
void learn_ratio_stage_init(struct NEURAL * neu,double  init_lr,double end_lr,int flush_steps,double ratio);

//exp（-x^2）衰减
//参数：neu，初始学习率，最小学习率，多少步数衰减到一半
void learn_ratio_exp_init(struct NEURAL * neu,double  init_lr,double end_lr,int half_steps);
//////////////////////////////////////////////////////////////////////////////////
//退出器，控制机器学习如何截断，控制机器学习停止
//
//
//平均损失函数截断，当所有数据平均损失函数小于 cutoff_prec时，停止
//参数：神经网络，多少步数判断一次，截断值，最大迭代步数
void exitor_avg_init(struct NEURAL * neu,int flush_steps,double cutoff_prec,int max_train_steps);

//
//
//
void exitor_max_init(struct NEURAL * neu,int flush_steps,double cutoff_prec,int max_train_steps);
////////////////////////////////////////////////////////////////////////////////
void print_function_init(struct NEURAL * neu,int print_type,int flush_steps);
///////////////////////////////////////////////////////////////////////////////
void wither_birth_creator_ASLP_init(struct NEURAL * neu,int flush_steps);


//////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//执行一次任务
//参数：神经网络neural，打印内容
void neural_run(struct NEURAL * neu);
void neural_train(struct NEURAL * neu);
void neural_free(struct NEURAL * neu);

void neu_printf(struct NEURAL * neu);
void neural_network_free(struct NEURAL_NETWORK * neus);
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//激活函数以及导数
double activate_function_default(double x);    //无
double activate_function_leaky_ReLU(double x);//拟合分段线性函数
double activate_function_sigmoid(double x);// 拟合，S型，渐近，压缩
double activate_function_Tanh(double x); //拟合光滑有界，连续函数
double activate_function_softplus(double x); // 和ReLu相似，但二阶连续可导
double activate_function_SINE(double x); //周期振荡函数
double activate_function_swish(double x);       // 无限可微，光滑，复杂函数


double G_activate_function_default(double x);
double G_activate_function_leaky_ReLU(double x);
double G_activate_function_sigmoid(double x);
double G_activate_function_Tanh(double x);
double G_activate_function_softplus(double x);
double G_activate_function_SINE(double x);
double G_activate_function_swish(double x);

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//率函损失函数

// MSE 均方误差函数
// 计算的yi数组，期望yi数组，yi数组多少个数，损失函数参数1,损失函数参数2
double loss_func_MSE(double * yi,double * yit, int numbers, double para1,double para2);
//MSE对a的导数，结果输出到数组output中
void G_loss_func_MSE( double * yi , double * yit , int numbers, double * output, double para1,double para2);


double loss_func_Huber(double * yi,double * yit, int numbers, double para1,double para2);
void G_loss_func_Huber( double * yi , double * yit , int numbers, double * output, double para1,double para2);

double loss_func_BCEwithLogitsLoss(double * yi,double * yit, int numbers, double para1,double para2);
void G_loss_func_BCEwithLogitsLoss( double * yi , double * yit , int numbers, double * output, double para1,double para2);
#endif