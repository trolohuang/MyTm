#ifndef SOAP_H
#define SOAP_H

#include "ctools.h"




struct SOAP
{

    int N_max;  //径向个数
    int L_max;  //角向个数
    double Rcut; //截断距离
    double sigma_density; //原子密度展开系数
    struct DATA * neig_atoms; // 邻近原子配位表，只读取前三列坐标
    
    //////radius functions
    double (*soap_gr)(struct SOAP * soap); // radial function
    struct DATA * gr_para;  // 径向函数参数
    struct DATA * gr_data; //预计算数据
    //////
    struct DATA * cnlm;    //cnlm  4D data; cnlm[N][L][M][real/img] 最后一个维度是复数的实部和虚部
    struct DATA * pnnl;    // pnnl 3D data； pnnl[N][N][L]
    struct DATA * integrade_table; // cnlm计算积分表  3D data; integrade_table[N][L][R_bin]
};
//对于 gr_data 和gr_para
//1. 对于高斯轨道基函数：gr_para 储存参数和要计算的数据；1D data[4]; sigma l n r; 
//                  ：gr_data 2D data[Nmax][Lmax]; 储存归一化系数Nn 
//
//
//
//
double soap_gr_function_gauss(struct SOAP * soap);
double soap_fc_function_cos(double r,double rcut);
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

//初始化SOAP算子
//参数：SOAP算子，径向个数，角向个数，截断距离，原子密度展开系数
struct SOAP * soap_init(int N,int L,double Rcut,double sigma_density);
void soap_free(struct SOAP * soap);

//初始化径向基函数
//参数：soap算子，径向函数类型，径向函数可选参数，径向函数可可可选参数
//径向函数类型：1：gauss轨道基函数；para1=高斯轨道基的展开系数
//径向函数类型：2：
void soap_gr_init(struct SOAP * soap,int gr_type,double para1,double para2);


//构造soap描述符p快速计算的积分表
//参数：soap算子，radius方向计算点的个数，radius积分点个数，径向函数参数1,径向函数参数2
//计算SOAP->integrade_table; 3D data; integrade_table[N][L][R_bin]
void soap_integration_table(struct SOAP * soap,int bin_r,int bin_integration);

//计算原子密度展开系数Cnlm
//参数：SOAP算子
//计算SOAP->cnlm; 4D data; cnlm[N][L][M][real/img] 最后一个维度是复数的实部和虚部
void soap_cnlm_generate(struct SOAP * soap);

void soap_cnlm_generate_withScaled(struct SOAP * soap,double scale_dis);

//计算展开系数功率谱 
//参数：SOAP算子
//计算SOAP->pnnl; 3D data； pnnl[N][N][L]
void soap_pnnl_generate(struct SOAP * soap);


// pnnl_1.pnnl_2
double soap_kernal(struct DATA * pnnl_1,struct DATA * pnnl_2);

#endif