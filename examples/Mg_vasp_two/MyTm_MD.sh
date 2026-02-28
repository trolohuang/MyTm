#!/bin/bash

# get parameters for MD  
tem_start=$1
tem_end=$2
nsw=$3
steps=$4
sub_steps=$5

# run MD by VASP code
cp INCAR_init INCAR
cp POSCAR_MyTm POSCAR
sed -i "1i TEBEG=${tem_start}" ./INCAR
sed -i "1i TEEND=${tem_end}" ./INCAR
sed -i "1i NSW=${nsw}" ./INCAR
mpirun -n 32 vasp6.4_gam > vasp.log

#get MD trajectory to file XDATCAR_MyTm and temperature, pressure and volume to file Temperature_MyTm, Volume_MyTm and Pressure_MyTm, 
cp XDATCAR XDATCAR_MyTm
grep 'T=' OSZICAR |awk '{print NR,$3}' > Temperature_MyTm
grep 'volume of cell :' OUTCAR |awk '{print NR,$5}' > Volume_MyTm
sed -i '1,2d' Volume_MyTm
grep 'total pressure  =' OUTCAR |awk '{print NR,$4}' > Pressure_MyTm

#MD result backup
mv XDATCAR XDATCAR_${steps}_${sub_steps}
mv OUTCAR OUTCAR_${steps}_${sub_steps}
mv OSZICAR OSZICAR_${steps}_${sub_steps}





