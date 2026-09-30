#ifndef MODE_BUILD_H
#define MODE_BUILD_H

void POSCAR_CrystalType2Orthogonal(struct POSCAR * poscar);
void POSCAR_setSuperCell(struct POSCAR * pos,int A,int B,int C);
struct POSCAR * mode_build(struct POSCAR * poscar_solid,struct POSCAR * poscar_liquid, struct MODE * mode_string);
#endif