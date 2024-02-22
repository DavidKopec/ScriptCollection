import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd

gro_file = ''
xtc_file = ''

u = mda.Universe(gro_file,xtc_file)

distance_threshold = 4
peptide_lenght=14

# aggregate identification

# how long is each peptide part of the aggreagate ?
# create a list / dictionary for each peptide, with values 'yes' and 'no' 
# ccorresponding to frame or time of the trajectory

aggregates_over_time = []

AAs = u.select_atoms('resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL')
peptides = [AAs.residues[i:i+14] for i in range(0,len(AAs.residues),14)]

data = dict()

for ts in u.trajectory[::100]:
    centroids = np.array([peptide.center_of_mass() for peptide in peptides])

    for i in range(len(centroids)):
        fill = ''
        for x in range(len(centroids)):
            data[f'p{i}'] = []
            if(i==x):
                continue
            distance = min(distances.dist(centroids[i],centroids[x],box=u.dimensions)[2])
            info = "Frame: %f\t Distance: %f\t" % (u.trajectory.frame,distance)
            print(info)

            if(distance<distance_threshold):
                data[f'p{i}'].append('YES')
            else:
                data[f'p{i}'].append('NO')
                
                pass

df = pd.DataFrame(data)
df.to_csv('ResidenceTime.csv')

