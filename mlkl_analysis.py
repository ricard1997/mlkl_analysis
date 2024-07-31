import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align

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
        self.gro = gro



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
    def get_rmsd(self, selection,start = 0, stop = -1, step=1, subgroups = None, sufix=""):
        if subgroups:
            subgroups += [selection]
        print(selection, subgroups)
        rmsd = rms.RMSD(self.u, select = selection, groupselections = subgroups)


        rmsd.run(start = start, stop= stop, step= step)
        rmsd_o = rmsd.results.rmsd

        columns = ["frame", "time", "full_rmsd"]
        n_columns = rmsd_o.shape[1]
        for i in range(n_columns-3):
            columns += [f"group{i}"]

        rmsd_df = pd.DataFrame(rmsd_o, columns = columns)
        rmsd_df.to_csv(f"rmsd{sufix}.dat",index = False)
        

    # Function to align protein with respect to the first frame or a reference CA structure (CA number of atoms must match)
    def align_prot(self, selection, ref_file = None, sufix = ""):
        mobile = self.u
        ref = mda.Universe(self.gro)
        if ref_file:
            ref = mda.Universe(ref_file)
        ref_at = ref.select_atoms(selection)
        aligner = align.AlignTraj(mobile ,ref_at ,select= selection ,filename = f'aligned_prot{sufix}.xtc').run()
        self.u.atoms.write(f"aligned_prot{sufix}.gro")

        temp_u = mda.Universe(f"aligned_prot{sufix}.gro", f"aligned_prot{sufix}.xtc")
        for ts in temp_u.trajectory[1:2]:
            temp_u.atoms.write(f"aligned_prot{sufix}.gro")
    


        

        


    # Compute RMSF of all the CA atoms and write it to a file
    def get_rmsf(self, selection = None, start = 0, stop = -1, step = 1, sufix = ""):
        protein = self.protein.select_atoms("name CA")
        if selection:
            protein = self.u.select_atoms(selection)
        
        R = rms.RMSF(protein)
        rmsf = R.run(start = start, stop = stop, step = step)

        rmsf_values = R.results.rmsf
        rmsf_df = pd.DataFrame()
        rmsf_df["resnum"] = protein.resnums
        rmsf_df["rmsf"] = rmsf_values
        #plt.plot(rmsf_df["resnum"], rmsf_df["rmsf"])
        rmsf_df.to_csv(f"rmsf{sufix}.dat", index = False)


    # Compute temporal RMSF of all the cA atom and write it to a file
    def get_time_rmsf(self, selection = None, interval = 20, start = 0, step = 1, stop =-1, sufix = ""):
        N = len(self.u.trajectory[start:stop:step])
        n_blocks = N // interval
        
        protein = self.protein.select_atoms("name CA")
        if selection:
            protein = self.u.select_atoms(selection)

        R = rms.RMSF(protein)

        rmsf_df = pd.DataFrame()
        rmsf_df["resnum"] = protein.resnums
        rmsf_dict = {}
        for i in range(n_blocks):
            middle = int(0.5*(i*interval+(i+1)*interval))
            R.run(start = i * interval, stop = ( i+1 ) * interval, step = 1)
            rmsf_dict[f"{middle}"] = R.results.rmsf
        rmsf_dict = pd.DataFrame(rmsf_dict)
        rmsf_df = pd.concat([rmsf_df, rmsf_dict], axis = 1)
        rmsf_df.to_csv(f"time_rmsf{sufix}.dat", index = False)

    
    # Time dependence distance between two groups of atoms
    def distance_two(self, selection1, selection2, write = False, sufix = ""):
        group1 = self.protein.select_atoms(selection1)
        group2 = self.protein.select_atoms(selection2)
        
        data = []
        for ts in self.u.trajectory:
            vector = group2.center_of_mass()-group1.center_of_mass()
            data.append([ts.frame, np.linalg.norm(vector)])

        data = np.array(data)
        if write:
            data = pd.DataFrame(data, columns =["frame", "distance"])
            data.to_csv(f"dist_{sufix}.dat")
        return data
        





        






