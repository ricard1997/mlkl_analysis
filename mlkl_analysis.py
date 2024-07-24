import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlb.pyplot as plt



############################### Code developed to analysis result from simulationf of protein in water, specifically MLKL #######################


# Defines the class protein by specifying the residues belonging by the protein and setups of the simulation such as timestep, among others
class Protein:
    def __init(self, gro, traj, selection_string, timestep = False, tpr = False):
        if tpr:
            self.u = mda.Universe(tpr, traj)
        else:
            self.u = mda.Universe(gro, traj)
        self.protein = self.u.select_atoms(selection_string)
        if timestep:
            selt.timestep = timestep
        self.selection_string = selection_string





    def extract_protein(self, start = 0, stop = -1, step = 1,
                            print_info = False,
                            custom_protein = False,
                            output_xtc = "only_protein.xtc",
                            output_gro = "only_protein.gro",
                            ):
        protein = self.protein
        protein.write(output_gro)
        
        seelction = self.selection_string
        if not custom_protein:
            selection = custom_protein
            protein = self.u.select_atoms(custom_protein)
            
        if not print_info:
            print(f"Writting atoms corresponding to protein in residues {selection}, start frame: {start}, final frame: {stop}, step: {step}")
            n_frames = len(self.u.trajectory[start:stop:step])
            print(f"Number of frames to print: {len(self.u.trajectory[start:stop:step])}")
            print(f"This correspond to {n_frames*timestep} timeunits")


        with mda.Writer(output_xtc, n_atoms = protein.n_atoms) as W:
	        for ts in self.u.trajectory[start:stop:step]:
		        W.write(protein)


        




    


