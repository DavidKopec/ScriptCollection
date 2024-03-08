import os
import sys
import MDAnalysis as mda
from MDAnalysis.selections.gromacs import SelectionWriter

dependencies = ["/storage/brno14-ceitec/shared/softmatter/kopec/Masters/Diploma/02MP1_CG/01MembranesOnly/002Dependencies", 
                "/storage/brno14-ceitec/shared/softmatter/kopec/Masters/Diploma/02MP1_CG/01MembranesOnly/001MDPs"]
gromacs_version = "gromacs:2022.3-plumed"

RADIUS = '1.2'
#D_MAX = '1.5'
DENS = '164'
FILENAME = f'output{os.getcwd().split('/')[-1]}'

CV1_ONLY  =False

def prepare_plumed_file(filename:str):
    fill=''
    with open(filename,'r') as template:
        fill = template.read()
    
    fill = fill.replace("<RADIUS>",RADIUS).replace("<DENS>",DENS).replace("<FILENAME>",FILENAME)
    
    with open('plumed.dat','w') as plumed:
        plumed.write(fill)
    # replace:
        # radius --> 1.2 currently 
        # CUTOFF -->  (1.5)
        # !!! DENS !!! --> calculated density average
        # BETA -->  15
        # DELTA --> 20
        # GAMMA --> 19
        # FILENAME on the end --> file where to print the data calculated 
    pass

def SMD():
    plumed_template_file = ''
    os.system("mkdir SMD")
    os.system(f"cp -r {dependencies[0]}/FF {dependencies[0]}/precycle_md ./SMD")

    if(CV1_ONLY==False):
        plumed_template_file = 'plumed-fixed.dat.template'
        os.system(f"cp {dependencies[0]}/SMD/plumed-fixed.dat.template .")
    else:
        plumed_template_file = 'CV1_only_template.dat'
        os.system(f"cp {dependencies[0]}/SMD/CV1_only_template.dat .")

    os.system(f"cp {dependencies[1]}/smd.mdp ./SMD")
    os.system("cp ./SIM/storage/md0004.gro ./SIM/storage/md0004.cpt ./SIM/storage/md0004.tpr ./SIM/system.top ./SIM/index.ndx ./SMD")
    
    os.chdir(os.getcwd()+"/SMD")
    os.system("mv md0004.gro system.gro")
    os.system("mv md0004.cpt system.cpt")


    prepare_plumed_file(plumed_template_file)
    
    os.system("psubmit cpu precycle_md ncpus=16 walltime=1d -y")

def main():
    if(len(sys.argv)>1 and sys.argv[2]!='1'):
        print("USAGE: python3 SMD.py [1]\n default will use plumed with CV1 and CV2\n with argument 1 --> only use CV1")

    proceed = input("Have you checked the values of denity and radius ???, press 'y' to continue")
    if(proceed!='y'):
        sys.exit()      
    SMD()
    
if __name__=="__main__":
    main()
