import os
import MDAnalysis as mda

GROMACS_V:str = "gromacs:2022.3-plumed"
PATH_DEP:str ="/storage/brno14-ceitec/shared/softmatter/kopec/Masters/Diploma/01MP1_AA/01_Affinities/A_Dependencies" #TO DO

SYSTEM:str="system.pdb"
PEPTIDE:str ="MP1_done.pdb"
TOPOLOGY:str = "system.top"
INDEX:str = "index.ndx"
MAIN_DIR:str = ''

def make_ndx():
    u = mda.Universe("system.gro")

    membrane = u.select_atoms("resname DOPC POPS DOPE")
    water = u.select_atoms("resname SOL")
    ions = u.select_atoms("resname NA CL SOD CLA")
    SOLU = u.select_atoms("resname POPS DOPC DOPE GLY ALA VAL CYS PRO LEU ILE MET TRP PHE LYS ARG HIS SER THR TYR ASN GLN ASP GLU") # Myabe select also the protein ?
    SOLV = u.select_atoms("resname CL NA SOL SOD CLA")
    
    with mda.selections.gromacs.SelectionWriter("index.ndx", mode="w") as ndx:
        ndx.write(membrane, name="Membrane")
        #ndx.write(plumed_selec, name="Selec")
        ndx.write(water, name="Water")
        ndx.write(ions, name="Ions")
        ndx.write(SOLU, name="MEMB_PROT")
        ndx.write(SOLV, name="W_ION")

def create_system():
    peptide_pdb:str=""
    system_pdb:str=""

    with open(f"{PEPTIDE}","r") as f:
        lines = f.readlines()
        for line in lines:
            if "ATOM" in line:
                peptide_pdb+=line
    
    with open(f"{SYSTEM}","r") as f:
        first_atom_line:bool = True

        lines = f.readlines()
        for line in lines:
            if "ATOM" in line and first_atom_line==True:
                system_pdb+="\n"+peptide_pdb+"\n"
                first_atom_line = False
            system_pdb+=line
            
    with open(SYSTEM,'w') as f:
        f.write(system_pdb)
    
    os.system(f"""module add {GROMACS_V}
                gmx editconf -f {SYSTEM} -o system.gro""")
    

def ionization():
    os.chdir(MAIN_DIR)

    os.system("mkdir 01_IONIZATION")
    os.system("cp -r {system.gro,system.top,FF} ./01_IONIZATION")
    os.system(f"cp -r {PATH_DEP}/ions.mdp ./01_IONIZATION")

    os.chdir("./01_IONIZATION")

    os.system(f"""module add {GROMACS_V}
                gmx grompp -f ions.mdp -c system.gro -p system.top -o ions.tpr -maxwarn 1
                gmx genion -s ions.tpr -pname NA -nname CL -neutral -p system.top -o system.gro""")
    make_ndx()

def simulation():
    os.chdir(MAIN_DIR)
    
    os.system("mkdir 02_MINIMIZATION")
    os.system("cp -r ./01_IONIZATION/{system.gro,system.top,FF,index.ndx} ./02_MINIMIZATION")
    os.system(f"cp -r {PATH_DEP}/minim.mdp ./02_MINIMIZATION")

    os.chdir(os.getcwd()+"/02_MINIMIZATION")

    os.system(f"""module add {GROMACS_V}
            gmx grompp -f minim.mdp -c system.gro -p system.top -r system.gro -n index.ndx -o minim.tpr
            gmx mdrun -deffnm minim -v""")
    
    os.chdir(MAIN_DIR)
    os.system('mkdir 03_SIM')
    os.system("cp -r ./02_MINIMIZATION/{minim.gro,system.top,FF,index.ndx} ./03_SIM")
    os.system(f"cp -r {PATH_DEP}/"+"{*.mdp,SIM}" +" ./03_SIM")

    os.chdir(os.getcwd()+'/03_SIM')
    os.system('psubmit gpu SIM ncpus=16 ngpus=2 walltime=1d -y')

    pass

def main():
    global MAIN_DIR
    MAIN_DIR = os.getcwd()
    
    create_system()
    ionization()
    simulation()

if __name__=="__main__":
    main()
