import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlb.pyplot as plt



############################### Code developed to analysis result from simulationf of protein in water, specifically MLKL #######################


# Defines the class protein by specifying the residues belonging by the protein and setups of the simulation such as timestep, among others
class Protein:
    def __init(self, gro, traj, selection_string, timestep, tpr = False):
        if tpr:
            self.u = mda.Universe(tpr, traj)
        else:
            self.u = mda.Universe(gro, traj)
        self.protein = self.u.select_atoms(selection_string)
        selt.timestep = timestep







    


