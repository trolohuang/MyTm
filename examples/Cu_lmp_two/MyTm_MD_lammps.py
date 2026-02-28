#!/home/trolo/anaconda3/envs/deepmd/bin/python

import sys
import os
import subprocess
from ase.io import read, write
from ase.data import atomic_masses

# ----------------------------
# 1. read the MD parameters
# ----------------------------

T_start = float(sys.argv[1])
T_end   = float(sys.argv[2])
steps   = int(sys.argv[3])

poscar = "POSCAR_MyTm"
lammps_data = "data_mytm.lammps"
dump_file = "dump_mytm.lammpstrj"
xdatcar_file = "XDATCAR_MyTm"

# ============================
# 2. Convert POSCAR to LAMMPS data and the Mass section
# ============================
atoms = read(poscar, format="vasp")
# set charge
		

# 写入原始 LAMMPS data 文件（无 Masses）
write(lammps_data, atoms, format="lammps-data", atom_style="atomic")

# 提取元素信息
symbols = atoms.get_chemical_symbols()
unique_symbols = []
for s in symbols:
    if s not in unique_symbols:
        unique_symbols.append(s)

# 生成 Masses 字段
masses_text = "Masses\n\n"
for i, sym in enumerate(unique_symbols, start=1):
    mass = atomic_masses[atoms[symbols.index(sym)].number]
    masses_text += f"{i} {mass:.6f}\n"
masses_text += "\n"

# 正确插入 Masses
with open(lammps_data, "r") as f:
    lines = f.readlines()

# 找到 Atoms 行的索引
atoms_index = None
for i, line in enumerate(lines):
    if line.strip().lower().startswith("atoms"):
        atoms_index = i
        break

if atoms_index is None:
    raise ValueError("未找到 'Atoms' 行，无法插入 Masses")

# 插入 Masses 前保持 header
new_lines = lines[:atoms_index] + [masses_text] + lines[atoms_index:]
with open(lammps_data, "w") as f:
    f.writelines(new_lines)


# ============================
# 3.set LAMMPS input files and run lammps
# ============================
lmp_in = "in.mytm"

with open(lmp_in, "w") as f:
    f.write(f"""
units           metal
atom_style      atomic
boundary        p p p

read_data       {lammps_data}


pair_style	eam/fs

pair_coeff	* *  Cu.eam.fs Cu




neighbor	2.0 bin
neigh_modify	delay 10 every 1 check yes

velocity	all create {T_start} 134563

# Thermodynamic variables
variable        t equal temp
variable        p equal press
variable        v equal vol
variable        s equal step

# 输出热力学量：两列（step value）
fix tprint all print 1 "${{s}} ${{t}}" file Temperature_MyTm screen no
fix pprint all print 1 "${{s}} ${{p}}" file Pressure_MyTm   screen no
fix vprint all print 1 "${{s}} ${{v}}" file Volume_MyTm     screen no

# 动力学
timestep        0.001
thermo          100
thermo_style    custom step temp press vol pe ke etotal

# NVT 升温过程
fix myNVT all npt temp {T_start} {T_end}  0.1 iso 0.0001 0.0001 1
#fix myNVT all nph iso 0.0001 0.0001 1
dump mydump all custom 10 {dump_file} id type x y z
dump_modify mydump sort id

run {steps}

unfix myNVT
""")


with open("lammps.output","w") as log:
	subprocess.run(["/opt/intel/oneapi/mpi/2021.16/bin/mpirun","-n","32","lmp_oneapi", "-in", lmp_in], stdout=log, stderr=log)


# ============================
# 5. Convert LAMMPS dump to XDATCAR and get thermodynamic parameters
# ============================
symbles_oputput = {1:"Cu"}
atoms_list = read(dump_file, format="lammps-dump-text", index=":")

for at in atoms_list:
	symbols = [symbles_oputput[t] for t in at.get_atomic_numbers()  ]
	at.set_chemical_symbols(symbols)

write(xdatcar_file, atoms_list, format="vasp-xdatcar")



with open("Temperature_MyTm","r") as f:
	lines = f.readlines()
lines = [l for l in lines if not l.lstrip().startswith("#")]
with open("Temperature_MyTm","w") as f:
	f.writelines(lines)

with open("Pressure_MyTm","r") as f:
	lines = f.readlines()
lines = [l for l in lines if not l.lstrip().startswith("#")]
with open("Pressure_MyTm","w") as f:
	f.writelines(lines)
	
with open("Volume_MyTm","r") as f:
	lines = f.readlines()
lines = [l for l in lines if not l.lstrip().startswith("#")]
with open("Volume_MyTm","w") as f:
	f.writelines(lines)	

