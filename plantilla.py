import MDAnalysis as mda
from mlkl_analysis import Protein
import os
import sys




gro = "vscratch/grp-vmonje/ricardox/c-phos-project/2ubpmlkl/rep1/production/centered_prot.gro"
xtc = "vscratch/grp-vmonje/ricardox/c-phos-project/2ubpmlkl/rep1/production/centered_prot.xtc"


d_dir = "/vscratch/grp-vmonje/ricardox/d-phos-project/"
c_dir = "/vscratch/grp-vmonje/ricardox/c-phos-project/"

directories = {
                f"normalmlkl": [d_dir, "rep0", "rep1", "rep2"],
                f"345mlkl": [d_dir, "rep0", "rep1", "rep1"],
                f"347mlkl": [d_dir, "rep0", "rep1", "rep2"],
                f"2pmlkl":[d_dir, "rep0", "rep1", "rep2"],
                f"s345d": [c_dir,"rep0", "rep1"],
                f"q343a": [c_dir,"rep0", "rep1"],
                f"q343a_s345d": [c_dir,"rep0", "rep1"],
                f"2ubpmlkl":[c_dir, "rep1"],
}


# Set up working directory
home = "/vscratch/grp-vmonje/ricardox/d-phos-project/mlkl_analysis"
os.chdir(home)


def extract_xtc(directories):
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            try:
                os.chdir(f"{directories[key][0]}{key}/{rep}/production/")
                os.system(f"sbatch {home}/get_centered_protein.csh")
            except:
                print(f"Error while working on {key}/{rep}")


#extract_xtc(directories)


def extractions():
    os.chdir(home)
    os.makedirs("data", exist_ok = True)
    for key in list(directories.keys()):
        os.chdir(f"{home}/data")
        os.makedirs(f"{key}", exist_ok =True)
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}")
            os.makedirs(f"{rep}", exist_ok=True)
            os.chdir(f"{home}/data/{key}/{rep}")
            print(f"working on {os.getcwd()}")

            gro = f"{directories[key][0]}{key}/{rep}/production/centered_prot.gro"
            xtc = f"{directories[key][0]}{key}/{rep}/production/centered_prot.xtc"
            protein = Protein(gro, xtc, "protein", timestep = 0.1)
            # Extract only proteins
            protein.extract_protein(step = 10)

            ref = '/vscratch/grp-vmonje/ricardox/d-phos-project/ref_structure.gro'
            # Align the extractions
            selections = {"psk":"((resid 182-351 or resid 365-460) and name CA)",
                            "found" : "(resid 7-83 or resid 100-122 or resid 134-175 or resid 182-460) and name CA",
                            "4hbbrace": "((resid 7-83 or resid 100-122 or resid 134-175) and name CA)"
                            }
            new_gro = "only_protein.gro"
            new_xtc = "only_protein.xtc"
            protein_o = Protein(new_gro, new_xtc, "protein", timestep = 0.1)
            for part in list(selections.keys()):
                protein.align_prot(selections[part], ref_file = ref, sufix = part)
            




    
extractions()



#ref = '../ref_structure.gro'

#protein = Protein(gro , xtc, "protein", timestep = 0.1)


#selection = "(resid 7-83 or resid 100-122 or resid 134-175 or resid 182-460) and name CA"
#selection = "((resid 182-351 or resid 365-460) and name CA)"
#protein.align_prot(selection, ref_file = ref)



