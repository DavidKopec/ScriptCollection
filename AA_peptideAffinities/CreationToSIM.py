import os
import MDAnalysis as mda

GROMACS_V:str = "gromacs.2022.3-plumed"
PATH_DEP:str ="../../" #TO DO

SYSTEM:str="system.pdb"
PEPTIDE:str ="MP1_rdy.pdb"
TOPOLOGY:str = "system.top"
INDEX:str = "index.ndx"
MAIN_DIR:str = ''

def create_system():
    peptide_pdb:str=""
    system_pdb:str=""

    with open(f"{PEPTIDE}","r") as f:
        lines = f.readlines()
        for line in lines:
            if "ATOM" in line:
                peptide_pdb+=line
    
    with open(f"{SYSTEM}","x") as f:
        first_atom_line:bool = True

        lines = f.readlines()
        for line in lines:
            if "ATOM" in line and first_atom_line==True:
                system_pdb+="\n"+peptide_pdb+"\n"
                first_atom_line = False
            system_pdb+=line
        
        f.write(system_pdb)
    
    os.system(f"""module add {GROMACS_V}
                gmx editconf -f {SYSTEM} -o system.gro""")
    

def ionization():
    os.chdir(MAIN_DIR)

    os.system("mkdir 01_IONIZATION")
    os.system("cp -r {system.gro,system.top,index.ndx,FF} ./01_IONIZATION")
    os.system(f"cp -r {PATH_DEP}/ion.mdp ./01_IONIZATION")

    os.chdir("./01_IONIZATION")

    os.system(f"""module add {GROMACS_V}
                gmx grompp -f ion.mdp -c system.gro -p system.top -o ion.tpr -maxwarn 1
                gmx genion -s ion.tpr -pname NA -nname CL -neutral -p system.top -o system.gro""")

def simulation():
    os.chdir(MAIN_DIR)
    
    os.system("mdkir 02_MINIMIZATION")
    os.system("cp -r ./01_IONIZATION/{system.gro,system.top,FF} ./02_MINIMIZATION")
    os.system(f"cp r{PATH_DEP}/minim.mdp ./02_MINIMIZATION")

    os.chdir(os.getcwd()+"/02_MINIMIZATION")

    ### create index.ndx groups
    u = mda.Universe("system.gro")

    membrane = u.select_atoms("resname DOPC POPS DOPE")
    water = u.select_atoms("resname PW")
    ions = u.select_atoms("name NA CL")
    SOLU = u.select_atoms("resname POPS DOPC DOPE GLY ALA VAL CYS PRO LEU ILE MET TRP PHE LYS ARG HIS SER THR TYR ASN GLN ASP GLU") # Myabe select also the protein ?
    SOLV = u.select_atoms("name CL NA or resname PW")
    
    with mda.selections.gromacs.SelectionWriter("index.ndx", mode="w") as ndx:
        ndx.write(membrane, name="Membrane")
        #ndx.write(plumed_selec, name="Selec")
        ndx.write(water, name="Water")
        ndx.write(ions, name="Ions")
        ndx.write(SOLU, name="MEMB_PROT")
        ndx.write(SOLV, name="W_ION")

    os.system(f"""module add {GROMACS_V}
            gmx grompp -f minim.mdp -c system.gro -p system.top -r system.gro -n index.ndx -o minim.tpr
            gmx mdrun -deffnm minim -v""")
    
    os.chdir(MAIN_DIR)
    os.system('mdkir 03_SIM')
    os.system("cp -r ./02_MINIMIZATION/{minim.gro,system.top,FF} ./03_SIM")
    os.system(f"cp -r {PATH_DEP}/"+"{*.mdp,SIM}" +" ./03_SIM")

    os.system('psubmit gpu SIM ncpus=16 ngpus=2 walltime=1d -y')


    pass


def main():
    global MAIN_DIR
    MAIN_DIR = os.getcwd()
    pass

if __name__=="__main__":
    main
