import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import os

GROMACS_V ="gromacs:2022.3-plumed"

def main():
    os.system("gmx trjcat -f *.xtc -o completeTraj.xtc -dt 100")

    u = mda.Universe("md0001.gro","completeTraj.xtc")
    prot1 = u.atoms[:38]
    prot2 = u.atoms[39:77]

    for ts in u.trajectory[::100]:
        prot1 = u.atoms[:38].center_of_mass()
        prot2 = u.atoms[39:77].center_of_mass()
        time = u.trajectory.time
        print(distances.dist(prot1.position,prot2.position,box=u.dimensions))


    

if __name__=="__main__":
    main()