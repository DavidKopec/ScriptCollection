import os
import sys

import pandas as pd

GROMACS_V = 'gromacs:2022.3-plumed'


PATH_SYSTEM='' # path to system for which we want to calculate density
PATH_DENSITIES='' # path where to run this script and store calculated data
PATH_DEP='' # path to additional dependencies (plumed template etc.)

RADII = [1.2]       #[1.0,1.2,1.5] # radiuses to use

# calculate the average bead density from plumed driver densities for whole trajectory
def calculate_mean(radius:float):
    data = pd.read_csv(f'density_{radius}.dat', header=0, sep=' ')
    average = data['density'].mean()

    with open('output.dat','w') as f:
        f.write('%.1f\t%.0f\n' % (radius,average)) 
    pass


def main():
    
    # name of the system
    folder_name = PATH_SYSTEM.split('/')[-1]
    
    os.chdir(PATH_DENSITIES)
    
    # make a folder of this system in directory for densitites (all in one place :))
    os.system(f"mkdir {folder_name}")
    os.chdir(folder_name)
    
    # concatenate trajectory from md of the system
    os.system(f"""module add {GROMACS_V}
    gmx trjcat -f {PATH_SYSTEM}/SIM/storage/*.xtc -o traj.xtc""") # ! maybe start frm md0002 (md0001 may influence results in unwanted way)
    
    os.system(f"cp {PATH_SYSTEM}/SIM/storage/md_0004."+"{gro,tpr} .")
    os.system(f"cp {PATH_DEP}/make_ndx.py .")
    os.system("python3 make_ndx.py") # make ndx groups (important for plumed calculation)
    
    
    # creating a plumed for the calculation of the constant 
    for radius in RADII:
        # open a plumed-teplate
        fill=''
        
        # read the template
        with open("plumed-template",'r') as f:
            fill = f.read()
        
        # fill edited template to actual plumed file
        fill.replace("<RADIUS>",radius).replace("<FILE>",f"density_{radius}.dat")
        with open(f'plumed{radius}.dat','w') as f:
            f.write(fill)
        
        # run plumed computation
        os.system(f"""module add {GROMACS_V}
        export PLUMED_NUM_THREADS=12
        plumed driver --mf_xtc traj.xtc --plumed plumed.dat""")
        
        pass
    
    pass



if __name__=="__main__":
    main()