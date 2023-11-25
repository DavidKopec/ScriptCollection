# This script edits, minimizes, equilibrates and finally runs MD of a CHARMMGUI generated M2.2P systems

import os
import sys
import MDAnalysis as mda
from MDAnalysis.selections.gromacs import SelectionWriter

# Path to dependencies files (may differ depending on your setup)
dependencies = ["../../002*", "../../001*"]
gromacs_version = "gromacs:2022.3-plumed"

#Edits topology
def TopologyEdit():
    headers = """#include "FF/martini_v2.3P.itp"
    #include "FF/martini_v2.0_lipids_all_201506.itp"
    #include "FF/martini_v2.0_ions.itp"
    """
    fill=""
    f = open("system.top","r")
    for line in f:
        if ("#" not in line):
            fill+=line
    f.close()
    
    fill= headers+fill
    
    f = open("system.top","w")
    f.write(fill)
    f.close()

def Minimization():
    os.system("mkdir minim")
    os.system(f"cp -r system.* {dependencies[0]}/FF {dependencies[1]}/relax.mdp {dependencies[1]}/minim.mdp ./minim")
    os.chdir(os.getcwd()+"/minim")
    
    #os.system(f"module add {gromacs_version} \n gmx grompp -f minim.mdp -c system.gro -p system.top -o minim.tpr \n gmx convert -tpr -f minim.tpr -o minim.trr")
    
    # universe creation
    u = mda.Universe("system.gro")
    
    # selections
    membrane = u.select_atoms("resname DOPC POPS DOPE")
    plumed_selec = u.select_atoms("(resname DOPC POPS DOPE) and (name C1A D2A C3A C4A C1B C2B C3B C4B)")
    water = u.select_atoms("resname PW")
    ions = u.select_atoms("name NA CL")
    SOLU = u.select_atoms("resname POPS DOPC DOPE")
    SOLV = u.select_atoms("name CL NA or resname PW")
    SOLV2 = u.select_atoms("name CL NA or resname PW")
    
    with mda.selections.gromacs.SelectionWriter("index.ndx", mode="w") as ndx:
        ndx.write(membrane, name="Membrane")
        ndx.write(plumed_selec, name="Selec")
        ndx.write(water, name="Water")
        ndx.write(ions, name="Ions")
        ndx.write(SOLU, name="SOLU")
        ndx.write(SOLV, name="SOLV")
        ndx.write(SOLV2, name="W_ION")        
    
    os.system(f"module add {gromacs_version} \n gmx grompp -f relax.mdp -c system.gro -p system.top -n index.ndx -o relax.tpr -maxwarn 1 \n gmx mdrun -deffnm relax -v")
    #input("Make groups [Membrane] --> lipids and [W_ION] --> solvent \n press enter to continue...")
    #os.system(f"module add {gromacs_version} \n gmx make_ndx -f system.gro -o index.ndx \n gmx grompp -f relax.mdp -c system.gro -p system.top -n index.ndx -o minim.tpr \n gmx mdrun -deffnm minim -v")
    
    os.chdir(os.getcwd()+"/../")
    
    


def EQ_MD():
    os.system("mkdir SIM")
    os.system(f"cp -r ./minim/relax.gro ./minim/index.ndx system.top {dependencies[0]}/FF {dependencies[0]}/martini_eq {dependencies[1]}/*.mdp ./SIM")
    os.chdir(os.getcwd()+"/SIM")
    
    os.system("psubmit cpu martini_eq ncpus=8 walltime=1d -y")
    
def change_dependencies(new_dep:str):
    new_dep = new_dep.replace("[","")
    new_dep = new_dep.replace("]","")
    
    paths = new_dep.split(",")
    
    dependencies = [paths[1],paths[0]]
    

def main(*args):
    if("-h" in args):
        print("Usage: place your script in the folder with .gro and .top files")
        print("All dependencies like .mdp files or .itp files have to be stored in specific folders")
        print("Default folder for .mdp files is ../001, and ../002 for the rest ")
        print("To change location of dependencies folders, use: python3 CG_MinimToMD.py -dep[path_to_mdps,path_to_others]")
        sys.exit()
    elif(len(args)>=2):
        if("-dep" in args[1]):
            change_dependencies(args[1])
    
    TopologyEdit()
    Minimization()
    EQ_MD()
  
  
if __name__=="__main__":
    main()
