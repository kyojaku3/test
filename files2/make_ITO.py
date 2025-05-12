#!/usr/bin/env python3

import sys
import numpy as np
from ase.io import read, write
from ase.visualize import view

poscar_filepath="POSCAR_In2O3/POSCAR_In2O3_conventional_2218907"
atoms=read(poscar_filepath)

n_Sn=3
#n_Sn=4

outfile=f"POSCAR_In2O3_Sn{n_Sn}"
print(f"outfile={outfile}")

conc=n_Sn/atoms.get_volume()*1e24
print(f"n_Sn={n_Sn}, conc={conc} 1/cm3")

In_indices = [i for i, atom in enumerate(atoms) if atom.symbol == "In"]
O_indices = [i for i, atom in enumerate(atoms) if atom.symbol == "O"]
n_atoms = len(atoms)
n_In = len(In_indices)
n_O = len(O_indices)
In_ratio = (n_In - n_Sn)/n_atoms*100
Sn_ratio = n_Sn/n_atoms*100
O_ratio = n_O/n_atoms*100
print(f"In={In_ratio}%, Sn={Sn_ratio}%, O={O_ratio}%")
#sys.exit()

positions = atoms.get_positions()[In_indices]
selected_indices = []

for i in range(n_Sn):
    if i == 0:
        selected_indices.append(In_indices[0])
    else:
        remaining_indices = [idx for idx in In_indices 
            if idx not in selected_indices]
        msd_values = []
        for idx in remaining_indices:
            distances = atoms.get_distances(idx, selected_indices, mic=True)
            msd = np.mean(distances**2)
            msd_values.append((idx, msd))
        
        best_idx = max(msd_values, key=lambda x: x[1])[0]
        selected_indices.append(best_idx)
print(selected_indices)

for idx in selected_indices:
    atoms[idx].symbol = "Sn"

view(atoms)
write(outfile, atoms)

