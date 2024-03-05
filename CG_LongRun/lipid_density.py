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

z_cutoff = 0
peptide_cutoff = 3.5

traj_step = 10

u = mda.Universe(gro_file,xtc_file)

membraneCOG = u.select_atoms(f'{membrane_selection}',updating=True).center_of_geometry()
z_cutoff = membraneCOG[2]
upper_leaflet = u.select_atoms(f'{phoshate_selection} and ( prop z > {z_cutoff})',updating=True)
proximity_selection = u.select_atoms(f'(around 10 {peptide_selection}) and {phoshate_selection}',updating=True)

for ts in u.trajectory[::traj_step]:
    lower_leaflet = u.select_atoms(f'{phoshate_selection}') - upper_leaflet

    data['Frame'].append(u.trajectory.frame)
    data['upper'].append(len(upper_leaflet.residues))
    data['lower'].append(len(lower_leaflet.residues-proximity_selection.residues))

    print(str(len(lower_leaflet.residues-proximity_selection.residues))+f"\t Frame: {u.trajectory.frame}")

df = pd.DataFrame(data)
df.to_csv('leaflet_density.csv')