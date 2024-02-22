import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd

gro_file = 'storage/md0031.gro'
xtc_file = 'storage/md0031.xtc'

u = mda.Universe(gro_file,xtc_file)

distance_threshold = 4
peptide_lenght=14


AAs = u.select_atoms('resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL')
peptides = [AAs.residues[i:i+14] for i in range(0,len(AAs.residues),14)]

data = dict()

for ts in u.trajectory[::100]:
    COMs = [peptide.center_of_mass() for peptide in peptides]

    for i in range(len(COMs)):
        data[f'p{i}'] = []
        
        for x in range(len(COMs)):
            if(x<=i):
                continue
            distance = distances.distance_array(COMs[i],COMs[x],box=u.dimensions)
            info = "Frame: %f\t Distance: %f\t" % (u.trajectory.frame,distance)
            print(f'p{i}-p{x}')
            print(info)

            if(distance<distance_threshold):
                data[f'p{i}'].append('YES')
            else:
                data[f'p{i}'].append('NO')

for key,val in data:
    print(key,len(val))

#df = pd.DataFrame(data)
#df.hist()
#df.to_csv('ResidenceTime.csv')

