import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import os
from pathlib import Path

import matplotlib.pyplot as plt

GROMACS_V ="gromacs:2022.3-plumed"

def main():
    traj_file = Path(os.getcwd()+"/completeTraj.xtc")
    if traj_file.is_file():
        os.system(f"rm -fr {traj_file}")
        
    os.system(f"""module add {GROMACS_V}
              gmx trjcat -f *.xtc -o completeTraj.xtc -dt 500""")

    u = mda.Universe("md0001.gro","completeTraj.xtc")
    sys = u.select_atoms('all')
    
    traj_len = len(u.trajectory)
    distance_threshold = 4
    os.mkdir("CloseContacts")
    os.system("cp {center,index.ndx} ./CloseContacts")
    
    prot1 = u.atoms[:38]
    prot2 = u.atoms[39:77]
    fill=""
    info_fill=""

    for i,ts in enumerate(u.trajectory[::100]):
        prot1 = u.atoms[:38]
        prot2 = u.atoms[39:77]
        time = u.trajectory.time
        frame = u.trajectory.frame
        d = min(distances.dist(prot1,prot2,box=u.dimensions)[2])
        
        if(d<=distance_threshold):
            sys.write(f"cc{i}.gro",frames=u.trajectory[[int(round(frame))]])
            os.system(f"mv cc{i}.gro CloseContacts/")
            os.system(f"./center -c CloseContacts/cc{i}.gro -n index.ndx -o CloseContacts/cc{i}.gro")
            
            ucc = mda.Universe(f"CloseContacts/cc{i}.gro")
            prot1 = ucc.atoms[:38].center_of_mass(compound="residues")
            prot2 = ucc.atoms[39:77].center_of_mass(compound="residues")
            dist_matrix = distances.distance_array(prot1,prot2,box=u.dimensions)
            np.savetxt(f"CloseContacts/distance_matrix{i}.csv",dist_matrix,delimiter=",")
                      
            
        info = "Frame: %f\t Distance: %f\t" % (frame,d)
        info_fill+=info+"\n"
        print(info)
        fill+=str(d)+"\n"
    
    with open("distances_info","w") as f:
        f.write(info_fill)
    
    with open("2_p_distances.xvg","w") as f:
        f.write(fill)
        
        
    
if __name__=="__main__":
    main()