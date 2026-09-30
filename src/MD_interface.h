#ifndef MD_INTERFACE_H
#define MD_INTERFACE_H

struct POSCAR * VASP_read_initstructure(void);
struct POSCAR * VASP_read_modestructure(void);
struct POSCAR * VASP_read_MDstructure(void);
struct POSCAR ** VASP_read_MDtrajectory(int * number_frames);
void VASP_write_modestructure(struct POSCAR * poscar);
void VASP_write_runstructure(struct POSCAR * poscar);
void VASP_MD_relax(double TemBeg,double TemEnd,struct PARAMETERS *input,int state,int sub_state);
void VASP_MD_heatingUP(double TemBeg,double TemEnd,struct PARAMETERS *input,int state,int sub_state);
double VASP_get_MDavgPressure();
double VASP_get_MDavgTemperature();
double VASP_get_MDavgEnthalpy();
struct DATA * VASP_getMDthermodynamic();

#endif