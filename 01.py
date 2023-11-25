import os

import numpy as np
import pandas as pd
from math import *
import MDAnalysis as mda

GROMACS_V:str = "gromacs:2022.3-plumed"
PATH_SMD:str = f"{os.getcwd()}/SMD" 
PATH_DEP:str = "/mnt/storage-brno14-ceitec/shared/softmatter/kopec/Masters/Diploma/02MP1_CG/01MembranesOnly/001MDPs" 

# directory prep
# trajectory prep
def phase01() -> None:
    # assuming we are in the concrete membrane directory
    os.system("mkdir FE")
    os.system(f"cp -r {PATH_SMD}/"+"{system.top,index.ndx,FF} ./FE")
    os.system(f"cp -r {PATH_SMD}/storage/md0004.tpr ./FE")
    os.system(f"cp -r {PATH_DEP}/md.mdp ./FE")

    os.system(f"""module add {GROMACS_V}
            gmx trjcat -f {PATH_SMD}/storage/*.xtc -o FE/trjcat.xtc""")
    
    os.system(f"""module add {GROMACS_V}
            gmx trjconv -f FE/trjcat.xtc -s FE/md0004.tpr -ur rect -pbc mol -b 100000 -o FE/trajout.xtc""")

#preparation of smd_out.dat file    
def phase02() ->None:
    #prepare smd_out.dat
    os.system(f"mv {PATH_SMD}/output.dat {PATH_SMD}/storage/md0004.output.dat")
    
    fill=""
    for i in range(1,5):
        with open(f"{PATH_SMD}/storage/md000{i}.output.dat","r") as f:
            for line in f:
                if "#" in line:
                    continue
                fill+=line

    os.chdir(os.getcwd()+"/FE")          
    with open("tmp.dat","x") as f:
        f.write(fill)

    os.system("""sed '/^#/d; s/^ //g' tmp.dat | cut -f 1,2 -d ' ' > smd_out.dat 
              sed -i ' 1,5000d' smd_out.dat 
              sed -i '1 i\time cv' smd_out.dat""")
    

# creates windows equally distributed regarding frames 
def phase03_custom() ->None:

    if(os.getcwd().split("/")[-1] == "FE"):
        os.chdir(os.getcwd()+"/../")

    os.mkdir("FE/inputs")
    os.chdir(os.getcwd()+"/FE/inputs")
    
    tprFile = '../md0001.tpr'
    xtcFile = '../trajout.xtc'-0.1
    
    input_file:str = "../smd_out.dat"
    
    data = pd.read_csv(input_file, header=0, skip_blank_lines=True, sep=' ', dtype=np.float64)
    
    u = mda.Universe(tprFile,xtcFile)
    sys = u.select_atoms('all')
    
    traj_len = len(u.trajectory)
    div = 120 #space between each frame
    
    frames_arr = np.arange(0,traj_len,div)
    
    Out:str = ""
    for win_nr,frame in enumerate(frames_arr):
        if(frame>traj_len):
            print("reached the end of the trajectory !!!")
            break
        temp:str=""
        win_name = "win-"+str(win_nr) + ".gro"
        time = int(data['time'][frame])
        
        temp = '%s\t%f\t%i\t%i' % (win_name,time,int(round(frame)),win_nr)
        Out += temp+"\n"
        print(temp)
        sys.write(win_name,frames= u.trajectory[[int(round(frame))]])
        
    with open("windows","x") as f:
        f.write(Out)
    
    os.chdir(os.getcwd()+"/../")    
    
#windows generation
def phase03() ->None:
    os.mkdir("FE/inputs")
    os.chdir(os.getcwd()+"/FE/inputs")
    
    tprFile = '../md0001.tpr'
    xtcFile = '../trajout.xtc'

    input_file:str = "../smd_out.dat"

    dmax = 2.20
    dmin = -0.10
    step = 0.04

    windows = np.arange(dmin,dmax,step)

    data = pd.read_csv(input_file, header=0, skip_blank_lines=True, sep=' ', dtype=np.float64)

    u = mda.Universe(tprFile,xtcFile)

    sys = u.select_atoms('all')

    Out:str = ""
    for win_nr,window in enumerate(windows):
        temp:str=""
        frame = round(np.abs(data['cv']-window).argmin()/5.5)
        time = int(data['time'][frame])
        win_name = "win-"+str(win_nr) + ".gro"
        
        temp = '%s\t%f\t%i\t%f\t%i' % (win_name,time,int(round(frame)), window, win_nr)
        Out += temp+"\n"
        print(temp)
        sys.write(win_name,frames= u.trajectory[[int(round(frame))]])
    
    with open("windows","x") as f:
        f.write(Out)
    
    os.chdir(os.getcwd()+"/../")


#prepare the simulation folders
#manipulation with plumed parameters for each window
def phase04() -> None:

    if(os.getcwd().split("/")[-1]!="FE"):
        os.chdir(os.getcwd()+"/FE")

    os.system(f"cp {PATH_DEP}/"+"{md.mdp,plumed-template.dat} .")

    pos_init:float=-0.1 #-0.1 is default
    step:float = 0.042 #0.035 is default
    pos:float=0.0
    plumed:str = ""
    with open("plumed-template.dat","r") as f:
        plumed = f.read()

    for iter_n, _ in enumerate(open("inputs/windows",mode="r")):
        os.system(f"mkdir win_{iter_n}")
        pos = pos_init+step*iter_n
        with open(f"./win_{iter_n}/plumed.dat",mode="x") as f:
            fill = plumed.replace("<AT>",f"{pos}")
            f.write(fill)
        
        #copy all necesary files into simulation folder
        os.system(f"cp inputs/win-{iter_n}.gro ./win_{iter_n}")
        os.system(f"mv win_{iter_n}/win-{iter_n}.gro win_{iter_n}/system.gro")
        os.system("cp -r {system.top,index.ndx,md.mdp,FF}"+f" ./win_{iter_n}")
        os.system(f"cp {PATH_DEP}/precycle_md_FE ./win_{iter_n}")
    
    #plumed = plumed.replace()

def main(*args):
    phase01()
    phase02()
    phase03_custom()
    phase04()

if __name__=="__main__":
    main()

#plumed params
# R_0=<RADIUS> D_MAX=<CUTOFF> 
# cv_tails: FUNC=-d/<DENS>+1
# dist_tails: ALT_MIN={BETA=<BETA>}
# cv definition FUNC ... <DELTA> <GAMMA>
#umbrella restrain AT=<AT>

#diff between win 0-30
# <RADIUS> 1.2 - 1.2 --> no diff
# <DMAX> 1.5 - 1.5 --> no diff
# <DENS> 114 - 114 --> no diff
# <BETA> 15 - 15 --> no diff
# <DELTA> 20 - 20 --> no diff
# <GAMMA> 19 - 19 --> no diff

# <AT> -100 - 950 --> DIFF !!!!!
#   mathematical formula: posinit+step*win_num
