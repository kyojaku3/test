# ini:20210228, update:20250108
# /usr/bin/env python3
# usage:

import sys
import os
import shutil
from monty.os import cd
from time import perf_counter

from ase.io import read, write
from ase.visualize import view
from ase import units
from ase import Atoms
from ase.calculators.calculator import Calculator
# For Opt
from ase.optimize.optimize import Optimizer
from ase.constraints import FixSymmetry, ExpCellFilter, StrainFilter, UnitCellFilter
from ase.optimize import QuasiNewton, LBFGS, BFGS, FIRE
# For MD
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.nptberendsen import NPTBerendsen
from ase.md import MDLogger

#mattersim_pth = os.path.join(os.environ['MATTERSIM_MODEL_PATH'], 'mattersim-v1.0.0-1M.pth')
mattersim_pth = os.path.join(os.environ['MATTERSIM_MODEL_PATH'], 'mattersim-v1.0.0-5M.pth')


class ASECalculator:

    def __init__(
        self,
        atoms_in: Atoms,
        ):

        self.atoms = atoms_in.copy()
        self.cf = None
        self.opt = None
        self.dyn = None

    def set_CellFilter(
        self,
        optcell: int = 3,  # 1: StrainFilter, 2: ExpCellFilter, 3: UnitCellFilter
        hydrostatic_strain: bool = False,
        pressure=0.0,
        ):

        if optcell == 1:
            self.cf = StrainFilter(self.atoms, hydrostatic_strain=hydrostatic_strain)

        elif optcell == 2:
            self.cf = ExpCellFilter(self.atoms, hydrostatic_strain=hydrostatic_strain)

        elif optcell == 3:
            self.cf = UnitCellFilter(self.atoms, scalar_pressure=pressure * units.GPa)

        else:
            self.cf = self.atoms

    def set_Optimizer(
        self,
        optimizer: str = "LBFGS", # (BFGS, LBFGS, FIRE, QuasiNewton)
    ):
       
        if optimizer == "BFGS":
            self.opt = BFGS(self.cf,
                trajectory = self.traj_filename,
                logfile = self.log_filename)

        elif optimizer == "LBFGS":
            self.opt = LBFGS(self.cf, 
                trajectory = self.traj_filename,
                logfile = self.log_filename)

        elif optimizer == "FIRE":
            self.opt = FIRE(self.cf,
                trajectory = self.traj_filename,
                logfile = self.log_filename)

        elif optimizer == "LBFGS":
            self.opt = QuasiNewton(self.cf,
                trajectory = self.traj_filename,
                logfile = self.log_filename)

    def set_MD(
        self,
        ion_dyn,
        time_step = 1.0,   # fsec
        num_interval = 1,
        temperature_K = 0,  # Kelvin
        # NVT by Berendsen thermostat
        taut = 1.0,  # fs
        # NVT by Langevin thermostat
        friction_coeff = 0.002,
        # NVT by Nose-Hoover thermostat
        ttime  = 20.0,     # tau_T thermostat time constant in fsec
        # NPT by Nose-Hoover thermostat with Parrinello-Rahman barostat
        sigma   = 1.0,     # External pressure in bar
        pfactor = 2e6,     # Barostat parameter in GPa,     
    ):

        # Set the momenta corresponding to the given "temperature"
        MaxwellBoltzmannDistribution(
            self.atoms,
            temperature_K = temperature_K,
            force_temp = True)
#        Stationary(self.atoms)  # Set zero total momentum to avoid drifting
#        sys.exit()

# NVE
# https://docs.matlantis.com/atomistic-simulation-tutorial/ja/6_1_md-nve.html
        if ion_dyn == "nve-verlet" or ion_dyn == "nve":
            print(f"nve: {time_step}, {num_interval}, {self.log_filename}, {self.traj_filename}")
#            sys.exit()
            self.dyn = VelocityVerlet(
                self.atoms,
                time_step * units.fs,
                loginterval = num_interval,
#                logfile = self.log_filename,
                trajectory = self.traj_filename,
            )

# NVT
# https://docs.matlantis.com/atomistic-simulation-tutorial/ja/6_2_md-nvt.html
        elif ion_dyn == "nvt-berendsen":
            self.dyn = NVTBerendsen(
                self.atoms,
                time_step * units.fs,
                loginterval = num_interval,
                logfile = self.log_filename,
                trajectory = self.traj_filename,
                temperature_K = temperature_K,
                taut = taut * units.fs,
            )

        elif ion_dyn == "nvt-langevin":
            self.dyn = Langevin(
                self.atoms,
                time_step * units.fs,
                loginterval = self.num_interval,
                logfile = self.log_filename,
                trajectory = traj_filename,
                temperature_K = temperature_K,
                friction = friction_coeff,
            )


        elif ion_dyn == "nvt-nose" or ion_dyn == "nvt":
            self.dyn = NPT(
                self.atoms,
                time_step * units.fs,
                loginterval = num_interval,
                logfile = self.log_filename,
                trajectory = self.traj_filename,
                temperature_K = temperature_K,
                externalstress = 0.1e-6 * units.GPa,  # Ignored in NVT
                ttime = ttime * units.fs,
                pfactor = None,   # None for NVT
            )

# NPT
# https://docs.matlantis.com/atomistic-simulation-tutorial/ja/6_3_md-npt.html
        elif ion_dyn == "npt-nose" or ion_dyn == "npt":
            self.dyn = NPT(
                self.atoms,
                time_step * units.fs,
                loginterval = num_interval,
                logfile = self.log_filename,
                trajectory = self.traj_filename,
                temperature_K = temperature_K,
                externalstress = sigma * units.bar,
                ttime = ttime * units.fs,
                pfactor = pfactor * units.GPa * (units.fs**2),
            )

        elif ion == "npt-berendsen":
            self.dyn = NPTBerendsen(
                self.atoms,
                time_step * units.fs,
                loginterval = num_interval,
                logfile = self.log_filename,
                trajectory = self.traj_filename,
                temperature_K = temperature_K,
                pressure_au = 1.0 * units.bar,
                taut = 5.0 * units.fs,
                taup = 500.0 * units.fs,
                compressibility_au = 5e-7 / units.bar,
            )



    def set_Calculator(
        self,
        calculator: Calculator,
        fix_symmetry: bool = True,
        # For CellFilter
        optcell: int = 0,
        hydrostatic_strain: bool = False,
        pressure = 0.0,
        # For Optimizer
        optimizer = None,
        # For MD
        ion_dyn = None,
        time_step = 1.0,  # fsec
        temperature_K = 0,  # Kelvin
        # printout
        num_interval = 1,
        log_filename = "ase.log",
        traj_filename = "ase.traj",
        ):

        self.log_filename = log_filename
        self.traj_filename = traj_filename

        self.atoms.calc = calculator
        if fix_symmetry:
            self.atoms.set_constraint([FixSymmetry(self.atoms)])

        self.set_CellFilter(optcell,
             hydrostatic_strain, pressure)

        if calctype == 'scf' or calctype == 'relax' or calctype == 'vc-relax':
            self.set_Optimizer(optimizer)

        if calctype == 'md' and ion_dyn:
            self.set_MD(ion_dyn,
                time_step = time_step,
                num_interval = num_interval,
                temperature_K = temperature_K,
            )


    def run(self,
            calctype,
            fmax: float = 0.005,
            num_md_steps = 1,
        ):

        print(f"Current directory is {os.getcwd()}")
        print(f"calctype = {calctype}, fmax = {fmax}, num_md_steps = {num_md_steps}")
        if calctype == "scf":
            self.atoms.get_potential_energy()
#            self.atoms.calc.calculate(self.atoms)

        elif calctype == "relax" or calctype == "vc-relax":
            self.opt.run(fmax=0.005)

        elif calctype == "md":
            start_time = perf_counter()
            print(f"    imd     Etot(eV)    T(K)    stress(mean,xx,yy,zz,yz,xz,xy)(units.GPa)  elapsed_time(sec)")
            self.dyn.run(num_md_steps)



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('structfile', type=str)
    parser.add_argument('calctype', type=str,
        choices=["scf", "relax", "vc-relax", "md"], default="scf")
    parser.add_argument('-o', '--outdir', type=str, help='if not given, outdir=calcdir')
    parser.add_argument('-r', '--MATSIM_RUN', type=str, default='true')
    parser.add_argument('-nosym', '--is_nosym', action="store_true")
    parser.add_argument('-oc', '--optcell', type=int, choices=[0, 1, 2, 3], default=0, help="0: no cell relax, 1: StrainFilter, 2: ExpCellFilter, 3: UnitCellFilter")
    parser.add_argument('-hydro', '--is_hydrostatic', action="store_true")
    parser.add_argument('-opt', '--optimizer', choices=[None, "BFGS", "LBFGS", "FIRE", "QuasiNewton"], default='LBFGS')
    parser.add_argument('-fm', '--fmax', type=float, default=0.0)
    parser.add_argument('-pr', '--pressure', type=float, default=0.0)
    parser.add_argument('-ion', '--ion_dyn', choices=[None, "nve", "nve-verlet",
                        "nvt", "nvt-nose", "nvt-berendsen", "nvt-langevin",
                        "npt", "npt-nose", "npt-berendsen"], default=None)
    parser.add_argument('-dt', '--time_step', type=float, default=1.0, help="time step in fsec")
    parser.add_argument('-t', '--temp_K', type=float, default=0, help="Temperature in K")
    parser.add_argument('-ni', '--n_interval', type=int, default=1)
    parser.add_argument('-ns', '--n_steps', type=int, default=1)

    args = parser.parse_args()
    structfile = os.path.abspath(args.structfile)
    calctype = args.calctype
    MATSIM_RUN = args.MATSIM_RUN.lower() == 'true'
    print(f"structfile={structfile}")
    print(f"calctype={calctype}")
    print(f"MATSIM_RUN={MATSIM_RUN}")
#    sys.exit()

    outdir = calctype
    if args.outdir:
        outdir = args.outdir
    print(f"outdir={outdir}")

    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    atoms = read(structfile)
#    view(atoms)
#    sys.exit()

    # prepare potential
    import torch
    from mattersim.forcefield import MatterSimCalculator
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Running MatterSim on {device}")
    calculator = MatterSimCalculator(load_path=mattersim_pth, device=device)

    with cd(outdir):
        asecalc = ASECalculator(atoms)
        asecalc.set_Calculator(
            calculator,
            fix_symmetry = not args.is_nosym,
            optcell = args.optcell,
            hydrostatic_strain = args.is_hydrostatic,
            pressure = args.pressure,
            optimizer = args.optimizer,
            ion_dyn = args.ion_dyn,
            time_step = args.time_step,
            temperature_K = args.temp_K,
            num_interval = args.n_interval,
            )

        if MATSIM_RUN:
            print("Calculation is to be done") 
            asecalc.run(calctype,
                fmax = args.fmax,
                num_md_steps = args.n_steps,
            )


            if calctype == 'scf' or calctype == 'relax':

#                etot = asecalc.atoms.calc.get_total_energy()
                epot = asecalc.atoms.calc.get_potential_energy()
#                ekin = asecalc.atoms.calc.get_kinetic_energy()

#                print(f"Total Energy     : {etot:f} eV")
                print(f"Potential Energy : {epot:f} eV")
#                print(f"Kinetic Energy   : {ekin:f} eV")

            if calctype == 'relax' or calctype == 'vc-relax':
                write('final.poscar', asecalc.atoms)

