#ifndef NEP_H
#define NEP_H
#include "ctools.h"
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//调用流程
//nep_init
//nep_generate_Gi2
//nep_generate_Gi3
struct NEP
{

    int N;     //径向个数
    int L;     //角向个数
    double Rcut; //截断距离
    struct DATA * neigh; //邻近原子表
    //////////////////
    int NP; // default = 5; // 径向函数 P
    //////////////////
    struct DATA * cnlm;   // Gi3 的cnlm[N,L,M]
    struct DATA * Gi2;    //二阶
    struct DATA * Gi3;    //三阶
};


//
//初始化一个NEP算子
//参数：径向个数，角向个数，邻近原子截断半径
struct NEP * nep_init(int N, int L,double rcut);

//
//释放NEP内存
void nep_free(struct NEP * nep);

//构建二阶NEP算子，结果储存到NEP->Gi2; 一维DATA Gi2[N]
//参数：nep算子，缩放系数
//缩放系数是对邻近原子与中心原子距离除以scaled
void nep_generate_Gi2(struct NEP * nep,double scaled);

//构建三阶NEP算子，结果储存到NEP->Gi3; 三维DATA Gi2[N，N，L] ;
//参数：nep算子，缩放系数
//缩放系数是对邻近原子与中心原子距离除以scaled
void nep_generate_Gi3(struct NEP * nep,double scaled);



//二阶算子Gi2的kernal，向量内积
double nep_kernal_Gi2(struct DATA * Gi2_1,struct DATA * Gi2_2);  
//三阶算子Gi3的kernal，向量内积
double nep_kernal_Gi3(struct DATA * Gi3_1,struct DATA * Gi3_2);
//二阶算子归一化
void nep_Gi2_normal(struct DATA * Gi2);
//三阶算子归一化
void nep_Gi3_normal(struct DATA * Gi3);


//
//生成NEP机器学习描述符
//参数：输出描述符，元素分解配位表，元素类型数量，原子个数数量，nep算子，缩放系数，是否计算二阶NEP，是否计算三阶NEP，是否归一化二阶NEP，是否归一化三阶NEP
//输出描述符：是一个2D DATA[num_atom x describe],每一行记录一个原子的nep描述符，每一行由各个元素通道的描述符组合
//缩放系数是对邻近原子与中心原子距离除以scaled
struct DATA * nep_generate_nepdescribtor(struct DATA ** nets_decompose,int atom_num_type,int atom_num_all,struct NEP * nep,double scaled,bool gen_Gi2,bool gen_Gi3,bool normal_Gi2,bool normal_Gi3);


#endif 