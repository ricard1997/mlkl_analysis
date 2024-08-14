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

gro = "vscratch/grp-vmonje/ricardox/c-phos-project/2ubpmlkl/rep1/production/centered_prot.gro"
xtc = "vscratch/grp-vmonje/ricardox/c-phos-project/2ubpmlkl/rep1/production/centered_prot.xtc"


d_dir = "/vscratch/grp-vmonje/ricardox/d-phos-project/"
c_dir = "/vscratch/grp-vmonje/ricardox/c-phos-project/"
e_dir = "/vscratch/grp-vmonje/ricardox/e-phos-project/"

directories = {
    #            f"normalmlkl": [d_dir, "rep0", "rep1", "rep2"],
    #            f"345mlkl": [d_dir, "rep0", "rep1", "rep2"],
    #            f"347mlkl": [d_dir, "rep0", "rep1", "rep2"],
    #            f"2pmlkl":[d_dir, "rep0", "rep1", "rep2"],
    #            f"s345d": [e_dir,"rep0", "rep1"],
    #            f"s345ds347d":[e_dir, "rep0"],
                f"4btfalpha": [e_dir,"rep0"],
#                f"q343a": [e_dir,"rep0", "rep1"],
#                f"q343a_s345d": [e_dir,"rep0", "rep1"],
#                f"2ubpmlkl":[e_dir, "rep1"],
}


# Set up working directory
home = "/vscratch/grp-vmonje/ricardox/d-phos-project/mlkl_analysis"
os.chdir(home)


def extract_xtc(directories, batch = True):
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            try:
                os.chdir(f"{directories[key][0]}{key}/{rep}/production/")
                if batch:
                    os.system(f"sbatch {home}/get_centered_protein.csh")
                else:
                    os.system(f"{home}/get_centered_protein.csh")
                    
            except:
                print(f"Error while working on {key}/{rep}")





def check_files(directories):
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            test = os.path.isfile(f"{directories[key][0]}{key}/{rep}/production/centered_prot.xtc")
            if test:
                print(f"File centered_prot.xtc exist in {directories[key][0]}{key}/{rep}/production")
            else:
                print(f"File centered_prot.xtc DOES NOT exist in {directories[key][0]}{key}/{rep}/production")
    
    for key in list(directories.keys()):
        for rep in directories[key][1:]:

            file = f"{home}/data/{key}/{rep}/only_protein.xtc"
            test = os.path.isfile(f"{home}/data/{key}/{rep}/only_protein.xtc")
            if test:
                print(f"File only_protein.xtc exist in {home}/data/{key}/{rep}/production")
                u = mda.Universe(file.replace(".xtc", ".gro"), file)
                print(u.atoms.n_atoms)
                print(f"The lenght of the trajectory is {len(u.trajectory)}")
            else:
                print(f"File only_protein.xtc DOES NOT exist in {home}/data/{key}/{rep}/production")

            file = f"{home}/data/{key}/{rep}/aligned_protpsk.xtc"
            test = os.path.isfile(f"{home}/data/{key}/{rep}/aligned_protpsk.xtc")
            if test:
                print(f"File aligned_protpsk.xtc exist in {home}/data/{key}/{rep}/production")
                u = mda.Universe(file.replace(".xtc", ".gro"), file)
                print(u.atoms.n_atoms)
                print(f"The lenght of the trajectory is {len(u.trajectory)}")
            else:
                print(f"File aligned_protpsk.xtc DOES NOT exist in {home}/data/{key}/{rep}/production")

                
    



#extract_xtc(directories)


def extractions(directories, step = 10):
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
            protein.extract_protein(step = step)

            ref = '/vscratch/grp-vmonje/ricardox/d-phos-project/ref_structure.gro'
            # Align the extractions
            selections = {"psk":"((resid 182-351 or resid 365-460) and name CA)",
                            "found" : "(resid 7-83 or resid 100-122 or resid 134-175 or resid 182-460) and name CA",
                            "4hbbrace": "((resid 7-83 or resid 100-122 or resid 134-175) and name CA)"
                            }
            selections2 = {"psk":"((resid 177-347 or resid 360-455) and name CA)",
                            "found" : "(resid 2-78 or resid 95-117 or resid 129-170 or resid 177-455) and name CA",
                            "4hbbrace": "((resid 2-78 or resid 95-117 or resid 129-170) and name CA)"
                            }
            new_gro = "only_protein.gro"
            new_xtc = "only_protein.xtc"
            protein_o = Protein(new_gro, new_xtc, "protein", timestep = 0.1)
            for part in list(selections.keys()):
                
                if "alpha" in key:
                    print("here")
                    protein_o.align_prot(selections[part], ref_file = ref, selection2 = selections2[part], sufix = part)
                else:
                    protein_o.align_prot(selections[part], ref_file = ref, sufix = part)
            




    
def rmsds(directories):
    files = ["aligned_prot4hbbrace", "aligned_protfound", "aligned_protpsk"]
    for key in list(directories.keys()):
        os.chdir(f"{home}/data")
        os.chdir(f"{home}/data/{key}")
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}")
            print(f"working on {os.getcwd()}")
            for file in files:
                gro = f"{file}.gro"
                xtc = f"{file}.xtc"
                protein = Protein(gro, xtc,selection_string="resid 1-469", timestep = 1)

                # rmsd groups
                subgroups = ["resid 6-124 and name CA", "resid 125-180 and name CA", "resid 181-469 and name CA"]
                protein.get_rmsd("resid 1-469 and name CA",subgroups =  subgroups, sufix = f"_{file.replace('aligned_prot', '')}")
            

def rmsfs(directories):
    files = ["aligned_prot4hbbrace", "aligned_protfound", "aligned_protpsk"]
    for key in list(directories.keys()):
        os.chdir(f"{home}/data")
        os.chdir(f"{home}/data/{key}")
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}")
            print(f"working on {os.getcwd()}")
            for file in files:
                gro = f"{file}.gro"
                xtc = f"{file}.xtc"
                protein = Protein(gro, xtc,selection_string = "resid 1-469", timestep = 1)
                # rmsf 
                protein.get_rmsf("resid 1-469 and name CA", sufix = f"_{file.replace('aligned_prot', '')}")
            

def time_rmsfs(directories):
    files = ["aligned_prot4hbbrace", "aligned_protfound", "aligned_protpsk"]
    for key in list(directories.keys()):
        os.chdir(f"{home}/data")
        os.chdir(f"{home}/data/{key}")
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}")
            print(f"working on {os.getcwd()}")
            for file in files:
                gro = f"{file}.gro"
                xtc = f"{file}.xtc"
                protein = Protein(gro, xtc,selection_string = "resid 1-469", timestep = 1)
                # rmsf 
                protein.get_time_rmsf("resid 1-469 and name CA", sufix = f"_{file.replace('aligned_prot', '')}")
            
 
def plot_rmsds(directories, filename):

    fig = plt.figure(figsize =(20,30), dpi = 500)
    layout = []
    for key in list(directories.keys()):
        layout.append([f"{key}_rep0", f"{key}_rep1", f"{key}_rep2"])
    
    ax_dict = fig.subplot_mosaic(layout, sharey=True,
                                gridspec_kw={"hspace":0.5,"wspace":0.15})
    
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            data = pd.read_csv(f"{home}/data/{key}/{rep}/{filename}")
            columns = list(data.columns)
            for column in columns[2:-1]:
                ax_dict[f"{key}_{rep}"].plot(data["frame"], data[column], linewidth = 0.2, alpha = 0.5, label = column) 
                ax_dict[f"{key}_{rep}"].set_title(f"{key}_{rep}") 
                ax_dict[f"{key}_{rep}"].legend() 
    #plt.tight_layout()
    plt.savefig(f"plot_{filename.replace('.dat', '.png')}")
    plt.close()






def plot_main_rmsds(directories, filename):
    fig = plt.figure(figsize =(20,30), dpi = 500)
    layout = []
    for key in list(directories.keys()):
        layout.append([f"{key}", f"{key}"])
    
    ax_dict = fig.subplot_mosaic(layout, sharey=True,
                                gridspec_kw={"hspace":0.5,"wspace":0.15})
    
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            data = pd.read_csv(f"{home}/data/{key}/{rep}/{filename}")
            columns = list(data.columns)
            ax_dict[f"{key}"].plot(data["frame"], data[columns[2]], linewidth = 0.2, alpha = 0.5, label = rep) 
            ax_dict[f"{key}"].set_title(f"{key}") 
        ax_dict[f"{key}"].legend() 
    #plt.tight_layout()
    plt.savefig(f"plot_prindipal_{filename.replace('.dat', '.png')}")
    plt.close()



def plot_rmsfs(directories, filename):

    fig = plt.figure(figsize =(10,15), dpi = 500)
    layout = []
    for key in list(directories.keys()):
        layout.append([f"{key}", f"{key}"])
    
    ax_dict = fig.subplot_mosaic(layout, sharey=True,sharex = True,
                                gridspec_kw={"hspace":0.5,"wspace":0.15})
    
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            data = pd.read_csv(f"{home}/data/{key}/{rep}/{filename}")
            columns = list(data.columns)
            ax_dict[f"{key}"].plot(data[columns[0]], data[columns[1]], linewidth = 0.5, label = rep) 
            ax_dict[f"{key}"].set_title(f"{key}")
            ax_dict[f"{key}"].set_xlabel(columns[0])
            ax_dict[f"{key}"].set_ylabel(columns[1])
            ax_dict[f"{key}"].legend() 
        ax_dict[f"{key}"].axvspan(6, 125, color = "blue", alpha = 0.2)
        ax_dict[f"{key}"].axvspan(126, 128, color = "red", alpha = 0.2)
        ax_dict[f"{key}"].axvspan(343, 350, color = "black", alpha = 0.2)
    #plt.tight_layout()
    plt.savefig(f"plot_{filename.replace('.dat', '.png')}")
    plt.close()
           


# Plot all the rmsfs in the same plot
def plot_all_rmsfs(directories, filename):
    plt.close()
    axvspan = [(7,26,"blue"),
                (29,55, "blue"),
                (59,83, "blue"),
                (101,121, "blue"),
                (136,161, "red"),
                (163,175, "red"),
                (192,196, "black"),
                (200,228, "yellow"),
                (233,248, "black"),
                (259,278, "yellow"),
                (284,292, "black"),
                (296,317, "black"),
                (326,340, "yellow"),
                (344,352, "green"),
                (365,369, "green"),
                (370,376, "black"),
                (382,400, "black"),
                (408,420, "black"),
                (430,441, "black"),
                (444,448, "black"),
                (450,459, "black"),


    ]
    #fig = plt.figure(figsize =(10,15), dpi = 500)
    layout = []
    for key in list(directories.keys()):
        layout.append([f"{key}", f"{key}"])
    
    #ax_dict = fig.subplot_mosaic(layout, sharey=True,sharex = True,
     #                           gridspec_kw={"hspace":0.5,"wspace":0.15})
    plt.figure(figsize = (15,5))
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            data = pd.read_csv(f"{home}/data/{key}/{rep}/{filename}")
            columns = list(data.columns)
            plt.plot(data[columns[0]], data[columns[1]], linewidth = 0.5, label = f"{key}_{rep}") 
            plt.xlabel(columns[0])
            plt.ylabel(columns[1])
    plt.legend(bbox_to_anchor = (1.01,0.9))
    for axv in axvspan:
        plt.axvspan(axv[0], axv[1], color = axv[2], alpha = 0.2)
    plt.tight_layout()
    plt.savefig(f"plot_all_{filename.replace('.dat', '.png')}", dpi=1000)
    plt.close()
            

def plot_time_rmsfs(directories, filename):

    fig = plt.figure(figsize =(20,30), dpi = 500)
    layout = []
    for key in list(directories.keys()):
        layout.append([f"{key}_rep0", f"{key}_rep1", f"{key}_rep2"])

    #layout = [[row[i] for row in layout] for i in range(len(layout[0]))]
    
    ax_dict = fig.subplot_mosaic(layout, sharey=True,sharex= True,
                                gridspec_kw={"hspace":0.5,"wspace":0.15})
    plot = []    
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            data = pd.read_csv(f"{home}/data/{key}/{rep}/{filename}")
            columns = list(data.columns)

            resnums = data[columns[0]]
            data = data[columns[1:]]
            data_v = data.values
            resnums = list(resnums)
            #print(resnums)
            plot.append(ax_dict[f"{key}_{rep}"].imshow(data_v,cmap = "Spectral", extent = [int(columns[1]), int(columns[-1]), int(resnums[0]), int(resnums[-1])]))
            
            ax_dict[f"{key}_{rep}"].set_title(f"{key}_{rep}")
    #plt.tight_layout()
    fig.colorbar(plot[0])
    plt.savefig(f"plot_{filename.replace('.dat', '.png')}")
    plt.close()

def plot_time_individual_rmsfs(directories, filename):

    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            data = pd.read_csv(f"{home}/data/{key}/{rep}/{filename}")
            columns = list(data.columns)

            resnums = data[columns[0]].iloc[190:445]

            data = data[columns[1:]].iloc[190:445]
            dt = int(columns[1])-int(columns[2])
            dt2 = dt*0.5
            print(dt2)
            data_v = data.values
            data_v = data_v.T
            data_v = np.rot90(data_v)
            resnums = list(resnums)
            #print(resnums)
            plt.imshow(data_v,cmap = "Spectral", aspect = "auto",extent = [int(columns[1]), int(columns[-1]), int(resnums[0]), int(resnums[-1])])
            
            plt.title(f"{key}_{rep}")
    #plt.tight_layout()
            plt.colorbar()
            plt.savefig(f"plot_{key}_{rep}{filename.replace('.dat', '.png')}")
            plt.close()



def plot_dist_two(directories, selection1, selection2, sufix = ""):

    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}/")
            prot = Protein("./aligned_protpsk.gro", "./aligned_protpsk.xtc", "resid 1-469")
            dist = prot.distance_two(selection1, selection2)
            dist = pd.DataFrame(dist, columns = ["frame", "dist"])
            dist["rolling"] = dist["dist"].rolling(window = 50, center = True).mean()
            plt.plot(dist["frame"], dist["rolling"], label = f"{key}_{rep}", alpha = 0.5) 
    
    os.chdir(home)
    plt.legend(bbox_to_anchor=(1.01,0.9))
    plt.tight_layout()
    plt.savefig(f"plot_dist_two_{sufix}.png")
    plt.close()


def pca_for_all(directories, selection):

    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}/")
            prot = Protein("./aligned_protfound.gro", "./aligned_protfound.xtc", "resid 1-469")
            dist = prot.pca(selection = selection)
    return dist 
    os.chdir(home)


def tica_for_all(directories, selection):

    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}/")
            prot = Protein("./aligned_protpsk.gro", "./aligned_protpsk.xtc", "resid 1-469")
            dist = prot.tica(selection = selection)
    
    os.chdir(home)

def ttica_for_all(directories, selection):

    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}/")
            prot = Protein("./aligned_protpsk.gro", "./aligned_protpsk.xtc", "resid 1-469")
            dist = prot.ttica(selection = selection)
    
    os.chdir(home)


def tttica_for_all(directories, selection, lag):

    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}/")
            prot = Protein("./aligned_protpsk.gro", "./aligned_protpsk.xtc", "resid 1-469")
            dist = prot.tttica(selection = selection, lag = lag)
    return dist
    os.chdir(home)


def cluster_encore(directories):
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}/")
            prot = Protein("./aligned_protpsk.gro", "./aligned_protpsk.xtc", "resid 1-469")
            clusters = encore.cluster(prot.u)
            print(clusters)



def recorrer_string(string):
    def minus_five(match):
        return str(int(match.group())-5)
    modified_string = re.sub(r"\d+", minus_five, string)
    return modified_string
    

    
def applycluster(directories, selection, projector, classifier, selection2= None, original_data = None):

    dict_results = {}
    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            os.chdir(f"{home}/data/{key}/{rep}/")
            print(f"working on {home}/data/{key}/{rep}/")
            prot = Protein("./aligned_protpsk.gro", "./aligned_protpsk.xtc", "resid 1-469")
            
            data = prot.get_features(selection)
            if "alpha" in key:
                print(key)
                selection_s = recorrer_string(selection)
            else:
                selection_s = selection
            
            data = prot.get_features(selection_s)
            print(selection)

                
            projection = projector.transform(data)
            clusterization = classifier.transform(projection)

            
            if len(original_data) >0 :
                sns.kdeplot(x = original_data[:,0],y = original_data[:,1], levels = 1, thresh = 0.05)


            plt.scatter(*projection.T, c= clusterization, alpha = 0.1)
            plt.scatter(*classifier.cluster_centers.T, marker = "o", c = "black")


            for i in range(classifier.n_clusters):
                plt.annotate(f"{i}", classifier.cluster_centers[i], xytext=classifier.cluster_centers[i]+.1)


            plt.savefig(f"{key}{rep}classification.png")
            plt.close()
            unique, counts = np.unique(clusterization, return_counts=True)
            percentages = counts*100/len(clusterization)
            dict_clusters = {}
            for value, percentage in zip(unique, percentages):
                dict_clusters[f"{value}"] = percentage

            dict_results[f"{key}_{rep}"] = dict_clusters

    return dict_results


def check_file(file):
    value = os.path.isfile(file)
    if value:
        print(f"File {file} exists")
    else:
        print(f"Files {file} does not exist, make sure it does before rerun this code")
        return value

def generate_topol(file):
    file_w = open(file, "r")
    words = ["CLA", "POT", "TIP3"]
    new_topol = open(file.replace("topol.top", "new_topol.top"))
    lines = file_w.readlines()
    for line in lines:
        if not any(line.startswith(word) for word in words):
            new_topol.write(line)



    # Function to extract the tpr for only the protein if all the files are present
def generate_tprs(directories):
    files_needed = ["topol.top", "new_prod.mdp", "aligned_protpsk.gro"]

    for key in list(directories.keys()):
        for rep in directories[key][1:]:
            test = True
            os.chdir(f"{directories[key][0]}/{key}/{rep}/")
            print(f"Working on the tpr file for {directories[key][0]}/{key}/{rep}/")
            test = check_file(f"{directories[key][0]}/{key}/{rep}/topol.top")
            if test:
                generate_topol(f"{directories[key][0]}/{key}/{rep}/topol.top")
            test = check_file(f"{home}/new_prod.mdp")
            test = check_file(f"{home}/data/{key}/{rep}/aligned_protpsk.gro")

            if test:
                subprocess.run(f"gmx grompp -f {home}/new_prod.mdp -o new_prod.tpr -c {home}/data/{key}/{rep}/aligned_protpsk.gro -r {home}/data/{key}/{rep}/aligned_protpsk.gro -p {directories[key][0]}/{key}/{rep}/new_topol.top")




#####################################################################################################################################

# Run section

#####################################################################################################################################



#------ Check files and run rmsd, rmsf, and time rmsf ----------

#check_files(directories)
#extractions(directories)
#rmsds(directories)
#rmsfs(directories)
#time_rmsfs(directories)



# Provide a filenames for rmsd and rmsf with (i) 4hb and brace alignment
# (ii) all the atoms found in the original pdb
# (iii) psk domain

files = ["rmsd_4hbbrace.dat", "rmsd_found.dat", "rmsd_psk.dat"]
files = ["rmsd_4hbbrace.dat", "rmsd_found.dat", "rmsd_psk.dat"]

#for file in files:
#    plot_rmsds(directories, file)
#    plot_main_rmsds(directories, file)
#    plot_rmsfs(directories, file.replace("rmsd", "rmsf"))
#    plot_time_rmsfs(directories, file.replace("rmsd", "time_rmsf"))
#    plot_time_individual_rmsfs(directories, file.replace("rmsd", "time_rmsf"))






# ---- Dictionary to plot distances betwoeen to grups of atoms
dict_dist = {
    "alpha_5_6": ["resid 344-352 and name CA", "resid 365-369 and name CA"],
    "alpha_5_7": ["resid 344-352 and name CA", "resid 370-376 and name CA"],
    "alpha_5_8": ["resid 344-352 and name CA", "resid 383-399 and name CA"],
    "alpha_5_9": ["resid 344-352 and name CA", "resid 408-420 and name CA"],
    "alpha_5_6_7_9": ["resid 344-352 and name CA", "(resid 408-420 or resid 370-376 or resid 365-369) and name CA"],
}
#for key in list(dict_dist.keys()):
#    plot_dist_two(directories, dict_dist[key][0], dict_dist[key][1], sufix = key)

#for file in files:
    #plot_rmsds(directories, file)
    #plot_main_rmsds(directories, file)
#    plot_all_rmsfs(directories, file.replace("rmsd", "rmsf"))











# ----- Set up the directories we are working with ----------

directories = {
#                f"normalmlkl": [d_dir, "rep0", "rep1", "rep2"],
#                f"345mlkl": [d_dir, "rep0", "rep1", "rep2"],
#                f"347mlkl": [d_dir, "rep0", "rep1", "rep2"],
#                f"2pmlkl":[d_dir, "rep0", "rep1", "rep2"],
#                f"s345d": [e_dir,"rep1"],
#                f"s345ds347d":[e_dir, "rep0"],
                f"s345ds347dalpha":[e_dir, "rep0"],
                f"4btfalpha": [e_dir,"rep0"],
                f"4btfalpha_2pmlkl": [e_dir,"rep0"],
#                f"q343a": [e_dir,"rep0", "rep1"],
#                f"q343a_s345d": [e_dir,"rep0", "rep1"],
#                f"2ubpmlkl":[e_dir, "rep1"],
}



#extract_xtc(directories, batch = True)
extractions(directories, step = 1)
check_files(directories)



# ---------------









########### The following code is th eworkflow for a markov state model with 50 clusters and then its 
########### reduction to three metastable states

# ------- Reset the directories we are wotking with ---------
directories = {
#                f"normalmlkl": [d_dir, "rep0"],
#                f"345mlkl": [d_dir, "rep0"],
#                #f"347mlkl": [d_dir, "rep0"],
                f"2pmlkl":[e_dir, "rep0"]#, "rep1"],
#                f"s345d": [e_dir,"rep1"],
#                #f"q343a": [c_dir,"rep0", "rep1"],
#                f"q343a_s345d": [c_dir,"rep0"],
#                f"2ubpmlkl":[c_dir, "rep1"],
#                f"s345danton":[e_dir, "rep0", "rep2", "repnew"],
#                f"s345ds347d":[e_dir, "rep0"],
}



# ------ Residues corresponding to either a alpha-helix or to a B-sheet --------
axvspan = [(7,26,"blue"),
                (29,55, "blue"),
                (59,83, "blue"),
                (101,121, "blue"),
                (136,161, "red"),
                (163,175, "red"),
                (192,196, "black"),
                (200,228, "yellow"),
                (233,248, "black"),
                (259,278, "yellow"),
                (284,292, "black"),
                (296,317, "black"),
                (326,340, "yellow"),
                (344,352, "green"),
                (365,369, "green"),
                (370,376, "black"),
                (382,400, "black"),
                (408,420, "black"),
                (430,441, "black"),
                (444,448, "black"),
                (450,459, "black"),
                ]



# ------ Set up a selection string for the clustering ---------
selection_string = f"(resid {axvspan[0][0]}-{axvspan[0][1]}"
for item in axvspan[1:]:
    selection_string += f" or resid {item[0]}-{item[1]}"
selection_string += ") and name CA"

#selection_string = "(resid 7-83 or resid 100-122 or resid 129-170 or resid 177-455) and name CA"
#selection_string = "(resid 7-83 or resid 100-122 or resid 129-170 or resid 177-455) and name CA"

# -------- Reset the selection string dependeing on needs ------ 
selection_string = "(resid 6-469 and name CA)"

print(selection_string)









# ------ Do PCA and return the data and the pca model -----
# ------ Only one element should be in the directory ------
# In this case we are using only 2pmlkl rep0
data, pca = pca_for_all(directories, selection_string)
print(data, "pca", pca)
data_fitted = pca.transform(data)
print(data_fitted.shape)



# ----- Move to directory
os.chdir(f"{home}/data/2pmlkl/rep0/")



# ------- Generates the classifier model for 2pmlkl rep0 -----------
protein = Protein("aligned_protpsk.gro", "aligned_protpsk.xtc", selection_string)
classifier = protein.cluster(data_fitted, n_clusters = 50,sufix = "tica")

assignments = classifier.transform(data_fitted)


# ------- plot the cluster and its classification --------------
plt.close()
fig, ax = plt.subplots(1,1,figsize = (18,5))
ax.scatter(*data_fitted.T, c = assignments)
ax.scatter(*classifier.cluster_centers.T, marker = "o", c = "black")
for i in range(classifier.n_clusters):
    ax.annotate(f"{i}", classifier.cluster_centers[i], xytext=classifier.cluster_centers[i]+.1)
fig.savefig("clusterannotate.png")



# ----- Generates Markov State Model based on the clusters passed -----------
from deeptime.markov.msm import MaximumLikelihoodMSM
msm = MaximumLikelihoodMSM().fit(assignments, lagtime = 1).fetch_model()
print(f"Number of states: {msm.n_states}")
























# ------ Plot the trnaistion matrix neglecting the connectivity of low transition rates -----------
import networkx as nx
plt.close()
fig, ax = plt.subplots(1, 1, figsize=(10, 10))

threshold = 1e-2
title = f"Transition matrix with connectivity threshold {threshold:.0e}"
G = nx.DiGraph()
ax.set_title(title)
for i in range(msm.n_states):
    G.add_node(i, title=f"{i+1}")
for i in range(msm.n_states):
    for j in range(msm.n_states):
        if msm.transition_matrix[i, j] > threshold:
            G.add_edge(i, j, title=f"{msm.transition_matrix[i, j]:.3e}")

edge_labels = nx.get_edge_attributes(G, 'title')
pos = nx.fruchterman_reingold_layout(G)
nx.draw_networkx_nodes(G, pos, ax=ax)
nx.draw_networkx_labels(G, pos, ax=ax, labels=nx.get_node_attributes(G, 'title'));
nx.draw_networkx_edges(G, pos, ax=ax, arrowstyle='-|>',
                       connectionstyle='arc3, rad=0.3');
print("plotgi")
fig.savefig("graph.png")




#--------- Use pccs to coarse grain the MSM, here i am using 3 states ---------
pcca = msm.pcca(n_metastable_sets=3)
print(pcca.coarse_grained_transition_matrix)




# ---- Use membership probabilities to compute the membership for each metastable state
# ---- In this case the criteria is the max probability is enough to belong to the cluster
memberships = pcca.memberships

mapping_colors = ["b", "g", "r"]
max_indices = [np.argmax(row) for row in memberships]
max_colors = [mapping_colors[np.argmax(row)] for row in memberships] # Dim number of clusters
map_assignments = [max_colors[value] for value in assignments]
    
    

# ----- Plot the metastable states in the principal components
plt.close()
plt.scatter(*data_fitted.T, c = map_assignments)
plt.scatter(*classifier.cluster_centers.T, c = max_colors)
for i in range(classifier.n_clusters):
    plt.annotate(f"{i}", classifier.cluster_centers[i], xytext=classifier.cluster_centers[i]+.1)
plt.savefig("clusterannotate_1.png")


# ------ Print the metastable stated over time
plt.close()
plt.plot([max_indices[value] for value in assignments], label = "cl1")
plt.legend()
plt.savefig("temporal.png")
plt.close()



# ----- Plot the probability of belonging to a particular metastable state
fig, axes = plt.subplots(1, 2, figsize=(15, 10))
for i in range(len(axes)):
    ax = axes[i]
    ax.set_title(f"Metastable set {i+1} assignment probabilities")
    print(data_fitted.shape, pcca.memberships[assignments,i])
    ax.scatter(*data_fitted.T, c=pcca.memberships[assignments, i], cmap=plt.cm.Blues)
norm = mpl.colors.Normalize(vmin=0, vmax=1)
fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.Blues), ax=axes, shrink=.8);
fig.savefig("newplot.png")
plt.close()




# ---------- Directories to be included in the barplot of th emetastable states ------------
directories = {
                f"normalmlkl": [e_dir, "rep0", "rep1", "rep2"],
                f"345mlkl": [e_dir, "rep0", "rep1", "rep2"],
                f"347mlkl": [e_dir, "rep0", "rep1", "rep2"],
                f"2pmlkl":[e_dir, "rep0", "rep1", "rep2"],
                f"s345d": [e_dir,"rep0", "rep1"],
                f"s345ds347d":[e_dir, "rep0"],
                f"s345ds347dalpha":[e_dir, "rep0"],
                f"4btfalpha": [e_dir,"rep0"],
                f"4btfalpha_2pmlkl": [e_dir,"rep0"],
                f"q343a": [e_dir,"rep0", "rep1"],
                f"q343a_s345d": [e_dir,"rep0", "rep1"],
                f"2ubpmlkl":[e_dir, "rep1"],
}





# If not defined here, the code will continue using the selection string from the beggining which is the ideal
#selection_string = "(resid 7-83 or resid 100-122 or resid 134-175 or resid 182-460) and name CA"

plt.close()
perce = applycluster(directories, selection_string, pca, classifier, original_data = data_fitted)
perce = pd.DataFrame(perce)
perce = perce.transpose()


perce_melted = perce.reset_index().melt(id_vars=["index"], var_name = "Cluster", value_name = "Value")
perce_melted["pcca"] = perce_melted["Cluster"].apply(lambda x: max_indices[int(x)])


print(perce_melted)
perce_melted.rename(columns={"index":"Condition"}, inplace=True)
perce_melted = perce_melted.groupby(["Condition", "pcca"], as_index=False)["Value"].sum()
print(perce_melted)
plt.close()
plt.figure(figsize=(15, 8))
sns.barplot(data = perce_melted, x="Condition", y = "Value", hue = "pcca")
print(perce_melted)
os.chdir(home)
plt.xlabel("Replica")
plt.ylabel("Percentage")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("barplot.png")



print(perce)



#ref = '../ref_structure.gro'
#protein = Protein(gro , xtc, "protein", timestep = 0.1)



#selection = "(resid 7-83 or resid 100-122 or resid 134-175 or resid 182-460) and name CA"
#selection = "((resid 182-351 or resid 365-460) and name CA)"
#protein.align_prot(selection, ref_file = ref)



