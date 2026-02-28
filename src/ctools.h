#ifndef ctools_H
#define ctools_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdbool.h>
#include <time.h>
#include <mkl.h>
#include <getopt.h>
#include <gsl/gsl_sf_bessel.h>
#include <complex.h>
#include <stdbool.h>
#include "ml.h"
#include "nep.h"
#include "data.h"


#include "solid_detection.h"
#include "mode_build.h"
#include "procedure.h"
#include "MD_interface.h"
#include "soap.h"

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void INPUT_free(struct PARAMETERS * input);
void PARA_toFILE(int method);
struct PARAMETERS * read_parameters_from_file(char * filename);
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void DATA_getXYT(struct DATA * Tdata,int sel,double *x,double *y,double *z);
void DATA_setXYT(struct DATA * Tdata,int sel,double x,double y,double z);
double DATA_get(struct DATA*data, ...);
void DATA_set(struct DATA*data,double set_double, ...);
void DATA_free(struct DATA *data);
struct DATA * DATA_init(int number_int, ...);
void DATAs_free(struct DATA ** datas);
struct DATA * DATA_cpy(struct DATA * org);
struct DATA *DATA_shapeINIT_rect(int inside,int coor_type,double x0,double y0,double z0,double x1,double y1,double z1);
struct DATA ** DATAs_init(int numbers,int dimens, ...);
double read_parameters(char * filename,char * comds,char * split);
double DATA4D_get(struct DATA * data,int dime1,int dime2,int dime3,int dime4);
void DATA4D_set(struct DATA * data,double set_data, int dime1,int dime2,int dime3,int dime4);
double DATA3D_get(struct DATA * data,int dime1,int dime2,int dime3);
void DATA3D_set(struct DATA * data,double set_data, int dime1,int dime2,int dime3);
double DATA2D_get(struct DATA * data,int dime1,int dime2);
void DATA2D_set(struct DATA * data,double set_data, int dime1,int dime2);
double DATA1D_get(struct DATA * data,int dime1);
void DATA1D_set(struct DATA * data,double set_data, int dime1);
struct DATA * DATA1D_init(int dimen1);
struct DATA * DATA2D_init(int dimen1,int dimen2);
struct DATA * DATA3D_init(int dimen1,int dimen2,int dimen3);
struct DATA * DATA4D_init(int dimen1,int dimen2,int dimen3,int dimen4);
struct DATA * XY_in(char *file_name);
void XY_toFILE(struct DATA * output,char *file_name);
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int  MSD(struct POSCAR ** poscar,struct DATA **output,int All_fram,int init_fram,int cal_time,int skip_time,double time_step );
void POSCAR_getOrder(struct POSCAR * poscar,int atom_type,int *pre_order,int *end_order);
double POSCAR_getVolume(double Lattice[][3]);
void POSCAR_getPOS_order(struct POSCAR *pos,int sel_atoms,double *x,double *y,double *z);
void POSCAR_setPOS_order(struct POSCAR *pos,int sel_atoms,double x,double y,double z);

void POSCAR_free(struct POSCAR *poscar);
void POSCARS_free(struct POSCAR **poscar,int numbers);

struct POSCAR * XDATCAR_getAvgPOSCAR(struct POSCAR ** poscars,int numbers,int start_fram,int end_fram,int skip_fram);
void POSCAR_TranCoorType(struct POSCAR * pos,int type);//type == 0 crystal Else car
void perCrystal(struct POSCAR * poscar,int A,int B,int C,double *rx,double *ry,double *rz,double x,double y,double z);
void perCrystal_full(struct POSCAR * poscar,int A,int B,int C);
struct POSCAR * POSCAR_cpy(struct POSCAR * org);
void coorTurnC(double Lattice[][3],int type,double *d1,double *d2,double *d3,double x,double y,double z);
int POSCAR_setFrenkelDefect_byRand(struct POSCAR * poscar,int atom_type, struct DATA * shape_rect, double defect_scale);

void POSCAR_getPOS(struct POSCAR* pos,int atom,int number,double *x,double *y,double *z);
void POSCAR_setPOS(struct POSCAR* pos,int atom,int number,double x,double y,double z);
double PolynomialFit_LSM(struct DATA * input,struct DATA ** para,int n);
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int STRINGS_in_subfunction_check(char A);
struct STRINGS * STRINGS_in(char * str,int sizea,int sizeb);
char * STRINGS_get(struct STRINGS * str,int n);
void STRINGS_free(struct STRINGS * str);
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
struct POSCAR *  POSCAR_in(char * file_name);
void POSCAR_toFILE(struct POSCAR *poscar,char *filename);
struct POSCAR *  FPOSCAR_in(char * file_name);
void FPOSCAR_toFILE(struct POSCAR *poscar,char *filename);
struct POSCAR ** XDATCAR_in_advance(char *filename,int read_fram,int *number_poscar,int * MD_type);
int POSCAR_getAtomType(struct POSCAR * poscar,int atom_pos);


void POSCAR_toFILE(struct POSCAR *poscar,char *filename);
#endif