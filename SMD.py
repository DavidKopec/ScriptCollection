import os
import MDAnalysis as mda
from MDAnalysis.selections.gromacs import SelectionWriter

dependencies = ["../../002Dependencies", "../../001*"]
gromacs_version = "gromacs:2022.3-plumed"
ndx = 0

def SMD():
    os.system("mkdir SMD")
    os.system(f"cp -r {dependencies[0]}/FF {dependencies[0]}/plumed.dat {dependencies[0]}/precycle_md ./SMD")
    os.system(f"cp {dependencies[1]}/smd.mdp ./SMD")
    os.system("cp ./SIM/md.gro ./SIM/md.cpt ./SIM/md.tpr ./SIM/system.top ./SIM/index.ndx ./SMD")
    
    if(ndx!=0):
        #print(f"####### \n\n SMD --> ndx !=0 ndx={ndx} \n\n #######")
        #input()
        make_ndx()
    
    os.chdir(os.getcwd()+"/SMD")
    os.system("mv md.gro system.gro")
    os.system("mv md.cpt system.cpt")
    
    os.system("psubmit cpu precycle_md ncpus=16 walltime=1d -y")
    
def make_ndx():
    #os.chdir(os.getcwd()+"/SMD")
    # universe creation
    u = mda.Universe("md.tpr", "md.gro")
    
    # selections
    membrane = u.select_atoms("resname DOPC POPS DOPE")
    plumed_selec = u.select_atoms("(resname DOPC POPS DOPE) and (name C1A D2A C3A C4A C1B C2B C3B C4B)")
    water = u.select_atoms("resname PW")
    ions = u.select_atoms("resname ION")
    SOLU = u.select_atoms("resname POPS DOPC DOPE")
    SOLV = u.select_atoms("resname PW ION")
    
    with mda.selections.gromacs.SelectionWriter("index.ndx", mode="w") as ndx:
        ndx.write(membrane, name="Membrane")
        ndx.write(plumed_selec, name="Selec")
        ndx.write(water, name="Water")
        ndx.write(ions, name="Ions")
        ndx.write(SOLU, name="SOLU")
        ndx.write(SOLV, name="SOLV")



def main(*args):
    if(len(args)>2):
        #print("####### \n\n MAIN --> args > 2 \n\n #######")
        #input()
        ndx=1        

    SMD()
    
if __name__=="__main__":
    main()
