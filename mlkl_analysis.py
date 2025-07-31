import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
from MDAnalysis.analysis.hydrogenbonds import (HydrogenBondAnalysis as HBA)
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
from scipy.sparse.linalg import eigs
from scipy import stats
from deeptime.decomposition import TICA
from deeptime.clustering import KMeans
############################### Code developed to analysis result from simulationf of protein in water, specifically MLKL #######################


# Defines the class protein by specifying the residues belonging by the protein and setups of the simulation such as timestep, among others
class Protein:
    def __init__(self, gro, traj, selection_string, timestep = False, tpr = False):
        if tpr:
            print(f"Reading {tpr}, {traj}")
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
    def align_prot(self, selection, ref_file = None, selection2 = None, sufix = ""):
        mobile = self.u

        ref = mda.Universe(self.gro)
        if ref_file:
            ref = mda.Universe(ref_file)
        if selection2:
            ref_at = ref.select_atoms("resid 6-469")
            ref_at.write(f"ref_structure.gro") 
            ref_at = ref_at.select_atoms("name CA")
            ref_at.residues.resids = list(range(1,ref_at.n_atoms+1))
            print(ref_at.select_atoms("resid 1-50").resnames)
        else:
            ref_at = ref.select_atoms(selection)
        if selection2:
            print(self.u.select_atoms("resid 1-50 and name CA").resnames, "u", selection2)
            aligner = align.AlignTraj(mobile ,ref_at ,select= selection2 ,filename = f'aligned_prot{sufix}.xtc').run()
        else:
            aligner = align.AlignTraj(mobile ,ref_at ,select= selection ,filename = f'aligned_prot{sufix}.xtc').run()
        self.u.atoms.write(f"aligned_prot{sufix}.gro")

        temp_u = mda.Universe(f"aligned_prot{sufix}.gro", f"aligned_prot{sufix}.xtc")
        for ts in temp_u.trajectory[1:2]:
            temp_u.atoms.write(f"aligned_prot{sufix}.gro")
    


        

    def write_cluster_trajs(self, cluster_list, step):
        if len(cluster_list) != len(self.u.trajectory):
            print("Warning: The cluster list and the trajectory have different length")
            print(len(cluster_list), len(self.u.trajectory))
            print("Well try to make it match:")
            #print(len(list(range(0, len(self.u.trajectory)-1, step))), len(cluster_list))

        else:
            print("######################## The cluster list and the trajectory have the same lenght, we can work with thta######")
            

        cluster = set(cluster_list)
        cluster_dict = {}
        for value in cluster:
            cluster_dict[value] = [i for i, n in zip(list(range(0, len(self.u.trajectory)-1, step)), cluster_list) if n == value]

            with mda.Writer(f"clustered_traj_{value}.xtc") as W:
                count = 0
                for ts in self.u.trajectory[cluster_dict[value]]:
                    W.write(self.protein)
                    #print(ts.frame)
                    if count == 0:
                        self.protein.write(f"clustered_traj_{value}.gro")
                    count += 1
            

        
        


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
        
    # Time dependence distance between two groups of atoms
    def angle_4hb_psk(self, write = False, sufix = "", alpha = None, ref = None):
        

        if alpha is not None:
            groups = {"h1": ["resid 16-19", "resid 1-3"],
                    "h2":["resid 25-28", "resid 45-48"],
                    "h3": ["resid 74-77", "resid 55-58"],
                    "h4": ["resid 97-100", "resid 112-115"]}
            psk = ["resid 431-434", "resid 188-191"]
        else:
            groups = {"h1": ["resid 21-24", "resid 5-8"],
                    "h2":["resid 30-33", "resid 50-53"],
                    "h3": ["resid 79-82", "resid 60-63"],
                    "h4": ["resid 102-105", "resid 117-120"]}
            psk = ["resid 436-439", "resid 193-196"]

        fourhb = []
        for group in groups:
            vect = self.protein.select_atoms(groups[group][1]).center_of_mass() - self.protein.select_atoms(groups[group][0]).center_of_mass()
            fourhb.append(vect/np.linalg.norm(vect))
            
            if alpha:
                print(vect)



        fourhb = np.array(fourhb)
        fourhb = np.mean(fourhb, axis=0)
    
            
        print("at least start #########################")
        sel = [groups[group][1] for group in groups]
        print("at least start #########################")
        selection = sel[0]
        for atoms in sel[1:]:
            selection += f" or {atoms}"
        top_4hb = self.protein.select_atoms(selection).center_of_mass() - self.protein.select_atoms(psk[0]).center_of_mass()
            
        print(selection)

        vect = self.protein.select_atoms(psk[1]).center_of_mass() - self.protein.select_atoms(psk[0]).center_of_mass()
        print(vect) 

        if ref:
            print("test")
            fourhb = top_4hb
        angle = np.arccos(np.dot(vect,fourhb)/(np.linalg.norm(vect)*np.linalg.norm(fourhb)))
        angle = np.rad2deg(angle)
        return angle

    def get_angles(self, ref_selection = None, selection1 = None, selection2 = None):
        ref_selection = "resid 436-439" if ref_selection is None else ref_selection
        selection1 = "(resid 7-26 or resid 29-55 or resid 59-83 or resid 101-121) and name CA" # 4hb residues
        selection2 = "(resid 192-196 or resid 200-228 or resid 259-278 or resid 284-292 or resid 296-317 or resid 326-340) and name CA" # psk residues (only on one side and without the ones that move)
        
        sel1 = self.protein.select_atoms(selection1)
        sel2 = self.protein.select_atoms(selection2)
        ref = self.protein.select_atoms(ref_selection)
        vector1 = sel1.positions - ref.center_of_mass()
        vector2 = sel2.positions - ref.center_of_mass()

        vector1 = vector1 / np.linalg.norm(vector1, axis = 1, keepdims = True)
        vector2 = vector2 / np.linalg.norm(vector2, axis = 1, keepdims = True)
        angle = np.arccos(np.dot(vector1, vector2.T)).flatten()
        #print(angle.shape)
        return angle.reshape(-1)



    def get_features(self, selection = None,dist = "full",raw_pos = True,angles = False, start=0, stop=-1, step =1):
        
        sel = self.protein
        sel = sel.select_atoms("name CA")

        
        if isinstance(selection, list):
            #print(selection)
            dist_sel = []
            selection_string = f"({selection[0]} and name CA) "
            for i,sele in enumerate(selection):
                dist_sel.append(self.u.select_atoms(f"({sele}) and name CA"))
                if i > 0:
                    selection_string += f" or ({sele} and name CA)"
            if len(selection) == 1:
                selection_string = f"({selection[0]} and name CA)"
            sel = self.u.select_atoms(selection_string)
        else:
            sel = self.u.select_atoms(f"({selection}) and name CA")
            
        data = []
        from scipy.spatial.distance import cdist

        
        for ts in self.u.trajectory[start:stop:step]:

            feats = []
            if raw_pos:
                pos = sel.positions
                zero_centered = pos - sel.center_of_mass()
                raw_pos_feat = zero_centered.flatten()
                feats.append(raw_pos_feat)
            if dist == "full":
                
                pos = sel.positions
                dists = cdist(pos, pos)
                i, j = np.triu_indices_from(dists, k=1)
                full_interdist_feat = dists[i,j]
                feats.append(full_interdist_feat)
                #print(full_interdist_feat.shape)
            elif dist == "inter":
                dist_feats = []
                for i in range(len(dist_sel)):
                    for j in range(len(dist_sel)):

                        #if i == j:
                        #    posi = dist_sel[i].positions
                        #    posj = dist_sel[j].positions
                        #    disti = cdist(posi,posj)
                        #    k, l = np.triu_indices_from(dists, k=1)
                        #    dists = dists[k,l]
                        #    dist_feats.append(dist)

                        if i > j:
                            posi = dist_sel[i].positions
                            posj = dist_sel[j].positions
                            
                            disti = cdist(posi,posj)
                            flat_dist = disti.flatten()
                            dist_feats.append(flat_dist)
                partial_inter_feat = np.concatenate(dist_feats)
                feats.append(partial_inter_feat)
            if angles:
                angle = self.get_angles()
                
                feats.append(angle)

            #print("hetettete",feats, feats[1].shape, feats[0].shape)
            feats = np.concatenate(feats)
            #print(feats.shape)




            data.append(feats) #generates a vector (feature vector containing NX3 features, or NX3 + NXN features if dist=True)
        #print(len(data))
        data = np.array(data)
        #print(data.shape, "##########333")
        return data
        



    # Function computes the PCA for CA atoms for the protein or for the atoms given
    def pca(self, selection = None, start = 0, stop = -1, step = 1, data = False):
        print("Make sure that your trajectory is aligned before running this function other calculation")


        #print(data.shape)
        
        if data.shape or data:
            data = data
        else:
            sel = self.protein
            sel = sel.select_atoms("name CA")
            if selection:
                sel = self.u.select_atoms(f"({selection}) and name CA")
        #data = []
        #for ts in self.u.trajectory[start:stop:step]:
        #    zero_centered = sel.positions - sel.center_of_mass()
        #    data.append(zero_centered.flatten()) #generates a vector (feature vector containing NX3 features)

        #data = np.array(data)


            test_number = len(sel.atoms)
            test_number = test_number * 3 + (test_number*(test_number-1))/2 
            data = self.get_features(selection = selection,dist = True, start = start, stop = stop, step  = step)
            #print(data.shape, test_number)
        #std = np.std(data, axis = 0)
        #scaler = StandardScaler()
        #scaler.fit(data)
        #data = scaler.transform(data)


        pca = PCA(n_components = 2)


        pca.fit(data) # Find the  eigenvectors and then use them to compute the projections using the transfrom
        proj = pca.transform(data) # Projected data into the eigenvectors, in this case [nframes, 10] we can plot the first two columns as the principal component projections
        plt.plot(proj[:,0], label = "PC1")
        plt.plot(proj[:,1], label = "PC2")
        plt.savefig("testi.png")
        plt.close()
        #print(pca.explained_variance_ratio_)
        #print(pca.singular_values_)
        return data, pca        



    def tica(self, selection = None, start = 0, stop = -1, step = 1):
        print("Make sure that your trajectory is aligned before running this function other calculation")
        sel = self.protein
        sel = sel.select_atoms("name CA")
        if selection:
            sel = self.u.select_atoms(f"({selection}) and name CA")
        data = []
        for ts in self.u.trajectory[start:stop:step]:
            zero_centered = sel.positions - sel.center_of_mass()
            data.append(zero_centered.flatten()) #generates a vector (feature vector containing NX3 features)

        data = np.array(data)
        tica = FastICA(n_components = 2)
        tica.fit(data) # Find the  eigenvectors and then use them to compute the projections using the transfrom
        proj = tica.transform(data) # Projected data into the eigenvectors, in this case [nframes, 10] we can plot the first two columns as the principal component projections
        plt.plot(proj[:,0], label = "TIC1")
        plt.plot(proj[:,1], label = "TIC2")
        plt.legend()
        plt.savefig("test.png")
        #print(tica.explained_variance_ratio_)
        #print(tica.singular_values_)
        
        

    def compute_lag_cov(self,data, lag):
        n_samples, n_features = data.shape
        if lag!=0:
            data_lagged = data[lag:]
            data_original = data[:-lag]
        else:
            data_lagged = data
            data_original = data
            

        cov = np.dot(data_original.T, data_lagged)/(n_samples-lag)

        return cov


        

    def ttica(self,n_components = 2, selection = None,lag = 10, start = 0, stop = -1, step = 1):
        print("Make sure that your trajectory is aligned before running this function other calculation")
        sel = self.protein
        sel = sel.select_atoms("name CA")
        if selection:
            sel = self.u.select_atoms(f"({selection}) and name CA")
        data = []
        for ts in self.u.trajectory[start:stop:step]:
            zero_centered = sel.positions - sel.center_of_mass()
            data.append(zero_centered.flatten()) #generates a vector (feature vector containing NX3 features)

        data = np.array(data)

        cov = self.compute_lag_cov(data, lag)
        cov_o = self.compute_lag_cov(data, 0)
        mat = np.dot(np.linalg.inv(cov_o), cov)
        print(mat.shape)
        eigenvalues, eigenvectors = eigs(mat, k=3, which="LM")


        print(eigenvalues)

        print(eigenvectors.shape, data.shape)
        transformed = np.dot(data,eigenvectors)

        

        plt.plot(transformed[:,0], label = "PC1")
        plt.plot(transformed[:,1], label = "PC2")
        plt.plot(transformed[:,2], label = "PC2")
        plt.savefig("test_ttica.png")
        #print(pca.explained_variance_ratio_)
        #print(pca.singular_values_)



    def tttica(self,n_components = 3, selection = None,lag = 10, start = 0, stop = -1, step = 1):
        print("Make sure that your trajectory is aligned before running this function other calculation")
        sel = self.protein
        sel = sel.select_atoms("name CA")
        if selection:
            sel = self.u.select_atoms(f"({selection}) and name CA")
        data = []
        for ts in self.u.trajectory[start:stop:step]:
            zero_centered = sel.positions - sel.center_of_mass()
            data.append(zero_centered.flatten()) #generates a vector (feature vector containing NX3 features)

        data = np.array(data)


        tica = TICA(lagtime=lag, dim=n_components)
        tica.fit(data, lagtime = lag)
        model = tica.fetch_model()
        transformed = model.transform(data)

        x = transformed[:,0]
        y = transformed[:,1]
#print(values.shape, x, y)

        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        
        print(transformed)
        plt.plot(transformed[:,0], label = "PC1")
        plt.plot(transformed[:,1], label = "PC2")
        plt.plot(transformed[:,2], label = "PC3")
        plt.legend()
        plt.savefig(f"test_ttttica_{lag}.png")
        plt.close()



        X,Y = np.mgrid[xmin:xmax:500j,ymin:ymax:500j]
        positions = np.vstack([X.ravel(),Y.ravel()])
        values = np.vstack([x,y])
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T,X.shape)



        fig, ax = plt.subplots()

        ax.imshow(np.rot90(Z), extent = [xmin,xmax,ymin,ymax])




        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.savefig(f"final_tica_{lag}.png")
        plt.close()
        print("tica", transformed)
        return  data, model



    # Function to cluster data using deeptime
    def cluster(self, data,n_clusters = 3, sufix = ""):
        centers = KMeans(n_clusters = n_clusters,
                                init_strategy = "kmeans++",
                                max_iter = 0,
                                fixed_seed = 13,
                                n_jobs = 8,
                                )
        print(data)    
        clustering = centers.fit(data).fetch_model()
        clustering.transform(data)
        centers.initial_centers = clustering.cluster_centers
        centers.max_iter = 5000

        cluster_optimization = centers.fit(data).fetch_model()

        

        assignments = cluster_optimization.transform(data)

        plt.scatter(data[:,0], data[:,1],alpha = 0.1, c=assignments)
        plt.savefig(f"clust_{sufix}.png")
        plt.close()
        plt.plot(assignments)
        plt.xlabel("time $ns$")
        plt.ylabel("cluster")
        #print("values clusering", assignments)
        plt.savefig(f"tempclus{sufix}.png")

        return cluster_optimization
            









    def hb_protein_protein(self, sufix = ""):
        
        hbonds = HBA(universe=self.u,
		    between = ['protein', 'protein'],
		    d_a_cutoff = 3.2, d_h_a_angle_cutoff=150, update_selections=False)

        hbonds.run()

        np.save(f'numpyhb_{sufix}.dat', hbonds.results.hbonds)
        df = pd.DataFrame(hbonds.results.hbonds[:,:4].astype(int), columns = ['Frame', 'Donor_ix', 'Hydrogen_ix', 'Acceptor_ix',])
        df["Distances"] = hbonds.results.hbonds[:,4]
        df["Angles"] = hbonds.results.hbonds[:,5]
        df["Donor_resname"] = self.u.atoms[df.Donor_ix].resnames
        df["Acceptor_resname"] = self.u.atoms[df.Acceptor_ix].resnames
        df["Donor_resid"] = self.u.atoms[df.Donor_ix].resids
        df["Acceptor_resid"] = self.u.atoms[df.Acceptor_ix].resids
        df["Donor_name"] = self.u.atoms[df.Donor_ix].names
        df["Aceptor_name"] = self.u.atoms[df.Acceptor_ix].names
        df.to_csv(f'hbonds_data{sufix}.dat')










        #print(pca.explained_variance_ratio_)
        #print(pca.singular_values_)


