import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms


############################### Code developed to analysis result from simulationf of protein in water, specifically MLKL #######################


# Defines the class protein by specifying the residues belonging by the protein and setups of the simulation such as timestep, among others
class Protein:
    def __init__(self, gro, traj, selection_string, timestep = False, tpr = False):
        if tpr:
            self.u = mda.Universe(tpr, traj)
        else:
            self.u = mda.Universe(gro, traj)
        self.protein = self.u.select_atoms(selection_string)
        if timestep:
            self.timestep = timestep
        self.selection_string = selection_string




    ################### Code to extract protein ##########################
    def extract_protein(self, start = 0, stop = -1, step = 1,
                            print_info = False,
                            custom_protein = False,
                            output_xtc = "only_protein.xtc",
                            output_gro = "only_protein.gro",
                            ):
        protein = self.protein
        protein.write(output_gro)
        selection = self.selection_string
        if custom_protein:
            selection = custom_protein
            protein = self.u.select_atoms(custom_protein)
            
        if print_info:
            print(f"Writting atoms corresponding to protein in residues {selection}, start frame: {start}, final frame: {stop}, step: {step}")
            n_frames = len(self.u.trajectory[start:stop:step])
            print(f"Number of frames to print: {len(self.u.trajectory[start:stop:step])}")
            if self.timestep:
                print(f"This correspond to {n_frames*self.timestep} timeunits")

        with mda.Writer(output_xtc, n_atoms=protein.n_atoms) as W:
            for ts in self.u.trajectory[start:stop:step]:
                #print(ts.frame, protein.n_atoms)
                W.write(protein)

    # Compute RMSD of the full protein and subgroups (In case of MLKL: 4HB, brace and PsK) and save it in a txt file 
    def get_rmsd(self, selection,start = 0, stop = -1, step=1, subgroups = None):
        if subgroups:
            subgroups += [selection]
        rmsd = rms.RMSD(self.u, select = selection, groupselections = subgroups)


        rmsd.run(start = start, stop= stop, step= step)
        rmsd_o = rmsd.results.rmsd

        columns = ["frame", "time", "full_rmsd"]
        n_columns = rmsd_o.shape[1]
        for i in range(n_columns-3):
            columns += [f"group{i}"]

        rmsd_df = pd.DataFrame(rmsd_o, columns = columns)
        rmsd_df.to_csv("rmsd.dat",index = False)
        


    

        
    def get_rmsf(self, selection = None, start = 0, stop = -1, step = 1):
        protein = self.protein.select_atoms("name CA")
        if selection:
            protein = self.u.select_atoms(selection)
        
        R = rms.RMSF(protein)
        rmsf = R.run(start = start, stop = stop, step = step)

        rmsf_values = R.results.rmsf
        rmsf_df = pd.DataFrame()
        rmsf_df["resnum"] = protein.resnums
        rmsf_df["rmsf"] = rmsf_values
        plt.plot(rmsf_df["resnum"], rmsf_df["rmsf"])
        plt.savefig("test.png")
        rmsf_df.to_csv("rmsf.dat", index = False)










