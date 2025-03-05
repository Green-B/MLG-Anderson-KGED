import numpy as np
import pickle
import sys
import os

w = float(sys.argv[1]) # Anderson disorder strength - identifies files to be processed
ne = int(sys.argv[2]) # Steps for integrating over e
ns = int(sys.argv[3]) # Resolution for plotting s
n_to_combine = 0
if len(sys.argv)>4:
	n_to_combine = int(sys.argv[4]) # Number of files to combine
# Get a list of everything in the path and pick out only the ones that are files
jes_indiv_path = "jes_individual"
jes_comb_path = "jes_combined"
if not os.path.exists(jes_comb_path):
	os.makedirs(jes_comb_path)
onlyfiles = [f for f in os.listdir(jes_indiv_path) if os.path.isfile(os.path.join(jes_indiv_path, f))]
search = ("jes-w%s-ne%s-ns%s" % (w, ne, ns))
w_ne_ns_files = [f for f in onlyfiles if search in f]
# Prepare to average the data
jes = np.zeros((ne,ns))
num_deltas = np.zeros((ne,ns))
dos = np.zeros((ne))
# Make note of the number of unit cells and save the least number of unit cells used among all jes files
search = "ne%s-ns%s-n" % (ne, ns)
lowest_Ntot = int(w_ne_ns_files[0][w_ne_ns_files[0].index(search)+len(search):w_ne_ns_files[0].index("-pbc")])
# Cut down the list to take only as many as we want
if n_to_combine>0 and len(w_ne_ns_files)>n_to_combine:
	w_ne_ns_files = w_ne_ns_files[:n_to_combine]
n_files = len(w_ne_ns_files)
# Combine the files
for file in w_ne_ns_files:
	with open(jes_indiv_path+"/"+file, "rb") as infile:
		# Compare Ntot
		this_Ntot = int(file[file.index(search)+len(search):file.index("-pbc")])
		if this_Ntot!=lowest_Ntot:
			print("\nWarning: combining files with different Ntot - saving the lowest Ntot\n")
			if this_Ntot<lowest_Ntot:
				lowest_Ntot = this_Ntot
		# Now we actually average over jes
		# Old jes files do not have num_deltas, so be prepared to load a file without it
		this_jes, emin, emax, smin, smax, this_num_deltas, this_dos = pickle.load(infile)
		jes += this_jes
		num_deltas += this_num_deltas
		dos += this_dos
jes = jes/n_files
num_deltas = num_deltas/n_files
dos = dos/n_files

data = [jes, emin, emax, smin, smax, num_deltas, dos]
savefile = "jes-combined-w%s-ne%s-ns%s-n%s-nc%s.pkl" % (w, ne, ns, lowest_Ntot, n_files)
with open(jes_comb_path+"/"+savefile, "wb") as outfile:
	pickle.dump(data, outfile, pickle.HIGHEST_PROTOCOL)

print("Done combining jes")
