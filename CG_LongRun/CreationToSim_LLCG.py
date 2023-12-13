import os
import sys
import MDAnalysis as mda

GROMACS_V:str = "gromacs:2022.3-plumed"
DEPENDENCIES:str = "../A_Dependencies"
PROT = "GLY ALA VAL CYS PRO LEU ILE MET TRP PHE LYS ARG HIS SER THR TYR ASN GLN ASP GLU"


def prepare_peptides(count:int) ->None:
    coords = [1,1,0] # coordinates on start
    shifts = [[0,0],[1,0],[0,1],[1,1],[1,-1]] # shifts we want to use
    
    peptide_num = 1 # how many peptides were generated already
    while(peptide_num<count):
        for i,shift in enumerate(shifts):
            X,Y = shift
            for _ in range(2): # to do the shift in both (positive and negative direction)
                X=-X
                Y=-Y
                os.system(f"""module add {GROMACS_V}
                            gmx editconf -f MP1_rdy.pdb -O MP1_{i}.pdb -translate {X} {Y} 0""")#
                peptide_num+=1
                if(i==0 or peptide_num==count):
                    break
            if(peptide_num==count):break

def assemble_system(count:int):

    fill=""
    prot_fill=""
    first_atom = False

    for i in range(1,count+1):
        with open(f"MP1_{i}.pdb") as prot_pdb:
            for line in prot_pdb:
                if("ATOM" in line):
                    prot_fill+=line
    
    with open("system.pdb","r") as memb_pdb:
        for line in memb_pdb:
            if(first_atom==False):
                if("ATOM" in line):
                    fill+=prot_fill
                    first_atom = True
            fill+=line
    with open("sys_gen.pdb","w") as f:
        f.write(fill)
    os.system(f"""module add {GROMACS_V}
              gmx editconf -f sys_gen.pdb -o system.gro""")
    pass

def ionize():
    os.system(f"cp {DEPENDENCIES}/ions.mdp")
    os.system(f"""module add {GROMACS_V}
              gmx grompp -f ions.mdp -c system.gro -p system.top -o ions.tpr
              gmx genion -s ions.tpr -pname NA -nname CL -neutral -o system_ion.gro -p system.top""")
    
def index_file():
    u = mda.Universe("system_ion.gro")
    
    # selections
    membrane = u.select_atoms("resname DOPC POPS DOPE")
    water = u.select_atoms("resname SOL or resname W")
    ions = u.select_atoms("name NA CL CLA SOD or resname ION")
    SOLU = u.select_atoms(f"resname POPS DOPC DOPE or {PROT}")
    SOLV = u.select_atoms("name CL NA SOL SOD CLA W")
    
    with mda.selections.gromacs.SelectionWriter("index.ndx", mode="w") as ndx:
        ndx.write(membrane, name="Membrane")
        ndx.write(water, name="Water")
        ndx.write(ions, name="Ions")
        ndx.write(SOLU, name="MEMB_PROT")
        ndx.write(SOLV, name="W_ION")        
    pass

def std_sim():
    def minim():
        os.system("mkdir MIN")
        os.system("cp system_ion.gro index.ndx system.top ./MIN")
        os.system(f"cp {DEPENDENCIES}/minim*.mdp ./MIN")
        os.chdir(os.getcwd()+"/MIN")
        os.system(f"""module add {GROMACS_V}
                   gmx grompp -f minim1.mdp -c system_ion.gro -n index.ndx -p system.top -o minim1.gro -maxwarn 1
                   gmx mdrun -deffnm minim1 -v
                   gmx grompp -f minim2.mdp -c minim1.gro -n index.ndx -p system.top -o minim2.gro -maxwarn 1
                   gmx mdrun -deffnm minim2 -v""")
        os.chdir(os.getcwd()+"/..")
    def sim():
        os.system("mkdir SIM")
        os.system("cp minim2.gro index.ndx system.top ./SIM")
        os.system(f"cp {DEPENDENCIES}/*.mdp ./SIM")
        os.system(f"cp {DEPENDENCIES}/SIM ./SIM")
        os.chdir(os.getcwd()+"/SIM")
        os.system("psubmit cpu SIM ncpus=8 walltime=1d -y")
        os.chdir(os.getcwd()+"/..")

    minim()
    sim()

def main():
    if(len(sys.argv)!=2):
        print("Usage: python3 script.py peptide_count")
        sys.exit()
    
    prepare_peptides()
    assemble_system()
    ionize()
    index_file()
    std_sim()
    
    pass

if __name__ =="__main__":
    main()