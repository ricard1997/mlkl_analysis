import MDAnalysis as mda
from mlkl_analysis import Protein
import matplotlib as mpl
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import MDAnalysis.analysis.encore as encore
import seaborn as sns
import re
import subprocess
import matplotlib.patches as patches
import networkx as nx
from scipy.sparse import coo_matrix




gro = "../2pmlkl/consistent_plumed_r0/com_plumed_relaxation/test_umbrella/out_concatenated.gro"
xtc = "../2pmlkl/consistent_plumed_r0/com_plumed_relaxation/test_umbrella/out_concatenated.xtc"




protein = Protein(gro,xtc, selection_string = "protein")
angles = []
frames = []

temporal_dependent_angles = []
temp = []
count = 0
for ts in protein.u.trajectory:
    
    count = True
    value = protein.angle_4hb_psk(ref=True)
    temp.append(value)
    if int(ts.time) == 0:
        temporal_dependent_angles.append(temp)
        temp = []

    angles.append(value)
    frames.append(ts.time)
plt.plot(frames, angles)
plt.savefig("test_angles.png")

plt.close()

i = 0
for rep in temporal_dependent_angles[:6]:
    len_val = len(rep)
    sns.kdeplot(rep[-int(0.5*len_val):], label = f"number {i}")
    sns.kdeplot(rep[:-int(0.5*len_val)], label = f"number_final {i}")
    i += 1

plt.legend()
plt.savefig("kdes_angles_lats.png")
plt.close()
i = 0
for rep in temporal_dependent_angles[6:]:
    len_val = len(rep)
    sns.kdeplot(rep[-int(0.5*len_val):], label = f"number {i}")
    sns.kdeplot(rep[:-int(0.5*len_val)], label = f"number_final {i}")
    i += 1

plt.legend()
plt.savefig("kdes_angles_last.png")

#print(angles)
