#ifndef SOLID_DETECTION_H
#define SOLID_DETECTION_H

struct DATA ** build_coordinations_table_bySANN(struct POSCAR * input,int atom_type);
struct DATA ** build_coordinations_table_byRadius(struct POSCAR * input,int atom_type,double Rcutoff);
struct  DATA * RDF(struct POSCAR * poscar,double max_radius,double prec,int sample_atoms);
struct DATA * find_peak(struct DATA * input,int type,double sigma,int find_type);
void statistic(struct DATA * input, struct DATA ** output ,double prec);
int cTon(char *name);
void nToc(int n,char * re_name);


double solidfication_check(struct SOLID_CHECK * so_check, struct PARAMETERS * inputs);
struct SOLID_CHECK * solidfication_init(int type);
void solidfication_free(struct SOLID_CHECK * so_check);
int solidfication_check_stable(struct DATA * temperature, struct DATA * volume, struct DATA * pressure,double precs);
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double points2surface(double desx,double desy,double desz,double Vec3_1[],double  Vec3_2[]);
int Qsort_DATA_sub(struct DATA ** data,int left,int right,int dest);
void Qsort_DATA(struct DATA ** data,int left,int right,int dest);
struct DATA *  smear_Gaussian(struct DATA * input,double sigma);
struct DATA *  smear_linear(struct DATA * input,double sigma);
struct DATA  * Interpolation_liner(struct DATA * input);
void cartesian2spherical(double x,double y,double z,double * r,double *th,double *phi);
void spherical2cartesian(double *x,double *y,double *z,double  r,double th,double phi);
double factorial(int n);
double Lengendre(int L,double x);
double P_Legendre(int L,int M,double x);
void Sph_Har(int L,int M,double theta,double phi,double *re_real,double *re_imag );
double Vec_distance(double dx,double dy,double dz);
void Vec_cross(double Vec3_C[],double Vec3_A[],double Vec3_B[]);
double Vec_dot(double Vec3_A[],double Vec3_B[]);
void Vec_Nor(double Vec3_C[],double Vec3_A[]);
double Gaussian_functions(double a,double x,double x0,double t);
void matrix_inverse(double* A, int n);
#endif