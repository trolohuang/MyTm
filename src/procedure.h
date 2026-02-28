#ifndef PROCEDURE_H
#define PROCEDURE_H
double heatingUP(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions );
double bisection_tem(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions );
double bisection_sta(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions );
double smallcell(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions ) ;


double test_run(struct POSCAR * init_poscar,struct PARAMETERS *inputs,struct FUNCTIONS *functions ) ;
#endif