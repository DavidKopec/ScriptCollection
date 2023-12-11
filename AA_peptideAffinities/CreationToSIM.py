import os
import MDAnalysis as mda

GROMACS_V:str = "gromacs.2022.3-plumed"
PATH_DEP:str ="../../" #TO DO

SYSTEM:str="system.pdb"
PEPTIDE:str ="MP1_rdy.pdb"
TOPOLOGY:str = "system.top"
INDEX:str = "index.ndx"

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
                gmx editconf -f {SYSTEM} -system.gro""")
    

def ionization():
    os.system("mkdir 01_IONIZATION")
    os.system("cp -r {system.gro,system.top,index.ndx,FF} ./01_IONIZATION")
    os.system(f"cp -r {PATH_DEP}/ion.mdp ./01_IONIZATION")

    os.chdir("./01_IONIZATION")

    os.system(f"""module add {GROMACS_V}
                gmx grompp -f ion.mdp -c system.gro -p system.top -o ion.tpr
                gmx genion -s ion.tpr -pname NA -nname CL -neutral -p system.top -o system.gro""")

def simulation():
    if(os.getcwd().split("/")[-1] =="01_IONIZATION"):
        os.chdir(os.getcwd()+"../")
    
    os.system("mdkir 02_MINIMIZATION")
    os.system("cp -r ./01_IONIZATION/{system.gro,system.top,FF} ./02_MINIMIZATION")
    os.system(f"cp r{PATH_DEP}/minim.mdp ./02_MINIMIZATION")

    os.chdir(os.getcwd()+"/02_MINIMIZATION")

    ### create index.ndx groups
    u = mda.Universe("system.gro")

    membrane = u.select_atoms("resname DOPC POPS DOPE")
    #pumed_selec = u.select_atoms("(resname DOPC POPS DOPE) and (name C1A D2A C3A C4A C1B C2B C3B C4B)")
    water = u.select_atoms("resname PW")
    ions = u.select_atoms("name NA CL")
    SOLU = u.select_atoms("resname POPS DOPC DOPE")
    SOLV = u.select_atoms("name CL NA or resname PW")
    SOLV2 = u.select_atoms("name CL NA or resname PW")
    
    with mda.selections.gromacs.SelectionWriter("index.ndx", mode="w") as ndx:
        ndx.write(membrane, name="Membrane")
        #ndx.write(plumed_selec, name="Selec")
        ndx.write(water, name="Water")
        ndx.write(ions, name="Ions")
        ndx.write(SOLU, name="MEMB_PROT")
        ndx.write(SOLV, name="W_ION")
        ndx.write(SOLV2, name="W_ION")  

    os.system(f"""module add {GROMACS_V}
            gmx grompp -f minim.mdp -c system.gro -p system.top -r system.gro -n index.ndx -o minim.tpr
            gmx mdrun -deffnm minim -v""")

    pass


def main():
    pass

if __name__=="__main__":
    main