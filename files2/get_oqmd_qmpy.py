#!/usr/bin/env python3

import argparse
import csv
import sys
import qmpy
import qmpy.io.ase_mapper
from ase.io import write

parser = argparse.ArgumentParser(description="get property of oqmd entry")
parser.add_argument("oqmd_ids", type=int, nargs='*', help="Entry ID of OQMD")
parser.add_argument('-m', '--mode', type=str, choices=["new","add"], default="new")
parser.add_argument('-f', '--outcsvfile', type=str, default="oqmd_qmpy.csv")
parser.add_argument('--save_poscar', action="store_true")

args = parser.parse_args()


def main():
    data = []
    for oqmd_id in args.oqmd_ids:
        entry = OQMD_Interface(oqmd_id)
#        sys.exit()
        data.append({
                     "OQMD ID" : entry.get_oqmd_id(),
                     "ICSD ID" : entry.get_icsd_id(),
                     "Name" : entry.get_name(),
                     "Generic" : entry.get_generic(),
                     "Spacegroup" : entry.get_spacegroup(),
                     "Prototype" : entry.get_prototype(),
                     "# of Element Types" : entry.get_ntype(),
                     "# of Atoms" : entry.get_natom(),
                     "Volume [A^3]" : entry.get_volume(),
                     "Formation Energy [eV/atom]" : entry.get_formation_energy(),
                     "Band Gap [eV]" : entry.get_bandgap(),
                     "Stability [eV/atom]" : entry.get_stability(),
                     "Energy [eV/atom]" : entry.get_energy_pa(),
                    })
        

        print(f"OQMD ID = {entry.get_oqmd_id()}")
        print(f"ICSD ID = {entry.get_icsd_id()}")
        print(f"Name = {entry.get_name()}")
        print(f"Generic = {entry.get_generic()}")
        print(f"Spacegroup = {entry.get_spacegroup()}")
        print(f"Prototype = {entry.get_prototype()}")
        print(f"# of Element Types = {entry.get_ntype()}")
        print(f"# of Atoms = {entry.get_natom()}")
        print(f"Volume [A^3] = {entry.get_volume()}")
        print(f"Formation Energy [eV/atom] = {entry.get_formation_energy()}")
        print(f"Band Gap [eV] = {entry.get_bandgap()}")
        print(f"Stability [eV/atom] = {entry.get_stability()}")
        print(f"Energy [eV/atom] = {entry.get_energy_pa()}")
#        print("atoms = {}".format(oqmd_entry.get_ase_atoms()))
        print("")

        if args.save_poscar:
            poscar_prefix = f"{entry.entry.composition.name}_OQMD{oqmd_id}"            
            print(f"For oqmd_id={oqmd_id}, save structure as {poscar_prefix}")
            prim_cell = entry.structure.make_primitive(in_place=False)
            conv_cell = entry.structure.make_conventional(in_place=False)
            atoms = qmpy.io.ase_mapper.structure_to_atoms(entry.structure)
            atoms.symbols = [str(atom).split()[0] for atom in entry.structure.atoms]
#            atoms_prim = qmpy.io.ase_mapper.structure_to_atoms(prim_cell)
#            atoms_conv = qmpy.io.ase_mapper.structure_to_atoms(conv_cell)
            atoms_prim = entry.get_ase_atoms(is_primitive=True)
            atoms_conv = entry.get_ase_atoms(is_primitive=False)
            write(f"{poscar_prefix}.poscar",      atoms)
            write(f"{poscar_prefix}_prim.poscar", atoms_prim)
            write(f"{poscar_prefix}_conv.poscar", atoms_conv)

    field_name = list(data[0].keys())
    print("field_name = ",field_name)

    if args.mode == "new":
        with open(args.outcsvfile, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames = field_name)
            writer.writeheader()
            writer.writerows(data)

    elif args.mode == "add":
        with open(args.outcsvfile, "a", newline='') as f:
            writer = csv.DictWriter(f, fieldnames = field_name)
            writer.writerows(data)
        

class OQMD_Interface:

    def __init__(self, entry_id):
        self.entry = qmpy.Entry.objects.get(id=entry_id)
        if "static" in self.entry.calculations:
            self.calculation: qmpy.Calculation = self.entry.calculations["static"]
            self.structure:   qmpy.Structure   = self.entry.calculations["static"].output
            self.compositoin: qmpy.Composition = self.entry.calculations["static"].composition
        elif "standard" in self.entry.calculations:
            self.calculation: qmpy.Calculation = self.entry.calculations["standard"]
            self.structure:   qmpy.Structure   = self.entry.calculations["standard"].output
            self.composition: qmpy.Composition = self.entry.calculations["standard"].composition
        else:
            oself.calculation = None
            self.structure = None
            self.composition = None

        self.prim_cell = self.structure
        self.conv_cell = self.prim_cell.make_conventional(in_place=False)


    def get_ase_atoms(self, is_primitive=False):
        if is_primitive:
            atoms = qmpy.io.ase_mapper.structure_to_atoms(self.prim_cell)
            atoms.symbols = [str(atom).split()[0] for atom in self.prim_cell.atoms]
        else:
            atoms = qmpy.io.ase_mapper.structure_to_atoms(self.conv_cell)
            atoms.symbols = [str(atom).split()[0] for atom in self.conv_cell.atoms]
        atoms.set_initial_magnetic_moments(None)
        return atoms


    def get_oqmd_id(self):
        return self.entry.id

    def get_icsd_id(self):
        if self.entry.label == None:
            return self.entry.label
        else:
            return self.entry.label.split('-')[1]

    def get_name(self):
        return self.entry.name

    def get_generic(self):
        return None

    def get_spacegroup(self):
        return self.structure.spacegroup.hm  

    def get_prototype(self):
        if self.entry.prototype:
            return str(self.entry.prototype).split(' - ')[0]
        else:
            return None

    def get_ntype(self):
        return self.entry.ntypes

    def get_natom(self):
        return self.entry.natoms
  
    def get_volume(self):
        return self.calculation.volume

    def get_formation_energy(self):
        return self.calculation.formation_energy()

    def get_bandgap(self):
        return self.calculation.band_gap

    def get_stability(self):
        forms = self.calculation.formationenergy_set.filter(fit="standard")
        forms = forms.exclude(stability=None)
        if not forms.exists():
            return None
        return forms[0].stability

    def get_energy_pa(self):
        return self.calculation.energy_pa

if __name__ == "__main__":
    main()
