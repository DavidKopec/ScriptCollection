import MDAnalysis as mda
from MDAnalysis.analysis import leaflet

import pandas as pd
import matplotlib.pyplot as plt

data = dict()

data['Frame'] = []
data['upper'] = []
data['lower'] = []

gro_file = 'md0005.gro'
xtc_file = 'md0005.xtc'
peptide_selection = ('resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL')
membrane_selection = 'resname DOPC DOPE POPS'
phoshate_selection = 'name PO4'

traj_step = 10

u = mda.Universe(gro_file,xtc_file)

membraneCOG = u.select_atoms(f'{membrane_selection}',updating=True).center_of_geometry()
z_cutoff = membraneCOG[2]
upper_leaflet = u.select_atoms(f'{phoshate_selection} and ( prop z > {z_cutoff})',updating=True)
proximity_selection = u.select_atoms(f'(cyzone 8 5 -5 {peptide_selection}) and {phoshate_selection}',updating=True)

for ts in u.trajectory[::traj_step]:

    lower_leaflet = u.select_atoms(f'{phoshate_selection}').residues - upper_leaflet.residues - proximity_selection.residues

    num_lipids_upper = len(upper_leaflet.residues)
    num_lipids_lower = len(lower_leaflet.residues)

    box = u.dimensions[:2]
    box_area = box[0]*box[1]

    area_per_lipid_upper = box_area / num_lipids_upper
    area_per_lipid_lower = box_area / num_lipids_lower

    data['Frame'].append(u.trajectory.frame)
    data['upper'].append(area_per_lipid_upper)
    data['lower'].append(area_per_lipid_lower)

    print(f"Area per lipid for upper leaflet: {area_per_lipid_upper:.2f} Å^2")
    print(f"Area per lipid for lower leaflet: {area_per_lipid_lower:.2f} Å^2")

df=pd.DataFrame(data)
df.to_csv('lipid_area.csv')