import os
import sys

import matplotlib.pyplot as plt
import numpy as np

GROMACS_V = "gromacs:2022.3-plumed"
analyses = ('Potential','Pressure','Temperature','Density')

if(len(sys.argv)<=1):
    print("Usage: python/ std_analysis.py file.edr")
    sys.exit()
edr_file = sys.argv[1]

for analysis in analyses:
    os.system(f"""module add {GROMACS_V}
              echo '{analysis}' | gmx energy -f {edr_file} -o {analysis}.xvg""")
    a = np.loadtxt(f'{analysis}.xvg', comments=['@','#'])
    # avg = np.average(a)
    # plt.plot(a[:,0],a[:,1],'r.')
    # plt.plot(a[:,0],avg,"b-")
