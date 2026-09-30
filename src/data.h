#ifndef DATA_H
#define DATA_H

struct MODE 
{
    struct DATA * descripter;
    struct MODE * next;
};
struct SOLID_CHECK
{
    int type;            // 1=qlm, 2=local density, 3=ML , 4=soap,
    int run_type;        // 0=learn, 1=check
    struct POSCAR * solid;
    struct POSCAR * liquid;
    struct POSCAR ** poss_solid;
    struct POSCAR ** poss_liquid;
    int num_solid,num_liquid;

    //////////////for qlmi
    struct DATA ** best_table; //REFERENCES
    int num_ref;
    /////////////for local density

    struct DATA ** refer;
    struct DATA * refer_normal;
    int num_ld_ref;

    /////////////for ML
    double * input_data;
    double * output_data;
    struct NEURAL * neu;

};
struct PARAMETERS
{
    ///////////////////////////MD runing steps
    int vasp_NSW;  // the steps of MD running 
    double nsw_adaption_ratio;
    //int vasp_NSW_adaptive;// automatic adjust the vasp_NSW?


    //////////////////////////temperature adjustment method 
    double solid_tem,melt_tem,prec_tem,heat_ratio; // solid temperature , liquid temperature, the difference of solid and liquid, heating up ratio
    int method; //runing type //heating up or bisection // small cell
    int dichotomy_stable; // in the bisection method, we need thermodynamic stable for every steps?
    int dichotomy_structure_relax; // in the bisection method, relax the mode for every temperature?
    
    int max_step,max_sub_step; // runing steps of melt calculation
    int init_tem_check; // chek the init solid_tem and melt_tem?
    
    ////////////////////////// model building 
    int A,B,C; // supercell for A, B and C
    //int orthogonal; // cell to cube?
    struct MODE * modes;  //the mode describter
    int isliquid; //do we need to build liquid phase?


    ////////////////////////// MD interface
    char commands[500]; // MD interface type: python, bash and etc.
    char script[500];   // MD interface script path and file name

    ///////////////////////// solid liquid check
    double decision; // the precies of solid checking cutoff
    double decision_liquid; //the precies of liquid checking cutoff
    double atom_method[10];// solid and liquid check method type and parameters
                           // for qlm op: 1 L qlm_cutoff neigh_cutoff atom_neigh_cutoff
                           // for qlm op auto adjustment para: -1
                           // for LD op: 2 sigma ld_cutoff neigh_cutoff atom_neigh_cutoff
                           // for qlm op auto adjustment para: -2
                           // for ml   : 3 alpha beta sigma L ml_cutoff ml_neigh rcut


    struct SOLID_CHECK * so_check; // the atom check data
};
typedef struct POSCAR * (*read_structure)(void);
typedef struct POSCAR **(*read_structures)(int * );
typedef double (*get_parameters)(void);
typedef struct DATA * (*get_values)(void);
typedef void (*write_structure)(struct POSCAR * );
typedef void (*MD_run)(double , double , struct PARAMETERS *, int ,int);
struct FUNCTIONS
{
    read_structure read_initstructure;//read init structure
    read_structure read_modestructure;//read mode structure
    read_structure read_MDstructure;//read afterMD structure
    read_structures read_MDtrajectory;//read MD trajectory

    get_parameters get_MDavgTemperature;
    get_parameters get_MDavgPressure;
    get_parameters get_MDavgEnthalpy;
    get_values getMDthermodynamic;

    write_structure write_modestructure;//write mode to file
    write_structure write_runstructure;//write runing structure


    MD_run MD_relax;// run MD relax
    MD_run MD_heatingUP;// run MD temperure up or down
    MD_run MD_relax_ensemble;// run MD relax with specifix ensemble
    MD_run MD_FILEbackup;

};
struct DATA
{
    int dimen;      //// data的唯度
    int all_data;   //// data的数据个数
    int * dimen_data; // data各个唯度大小
    double * data_data; // data

    int data_group; //data数组启用时，记录data数组大小；通常在程序中使用第0个DATA中的data_group判断个数
    
};
/////////////////////////////////////////////////////////////////////////////
struct STRINGS
{
    int sizea; //最大单词个数
    int sizeb; // 单个单词最大长度
    int number_words;//词的总个数
    char ** strings;//字符串数组，每一个组元代表一个识别到的单词，
};
////////////// POSCAR/CHGCAR/ELFCAR 适用 
struct POSCAR
{
    char name[50];//POSCAR注释

    double scale;//坐标缩放
    double Lattice[3][3];//晶格坐标 [5]

    int atoms[50];             //  原子序数
    int num_atoms[50]; //  原子个数
    int num_type;            // 多少种原子
    int num_all;                // 原子总个数

    int coordinate_type;// 输出坐标类型： 0 分数坐标 1 笛卡尔
    struct DATA *atom_position; //原子位置 晶格坐标
    struct DATA * atom_position_ext;//原子位置笛卡尔坐标

    int  atom_fix ;              // 是否有原子固定的信息
    int  atom_velocity;   //是否有原子速度信息
    int  atom_external;  //是否有ELF或者charge等网格信息

    struct DATA *external;       // 额外数据； 储存电荷密度或者ELF等；分子动力学储存原子固定信息
    struct DATA *external_2;  //额外数据，储存分子动力学原子速度;
};

#endif