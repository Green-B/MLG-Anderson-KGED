import numpy as np
import pickle
import sys

w = float(sys.argv[1]) # Anderson disorder strength - identifies files to be processed
ne_old = int(sys.argv[2]) # Steps for integrating over e
ns_old = int(sys.argv[3]) # Resolution for plotting s
Ntot = int(sys.argv[4]) # Total number of unit cells in the supercell
nconfigs = int(sys.argv[5]) # Number of disorder configurations
ne_combine = int(sys.argv[6]) # Factor by which to reduce ne by joining bins
ns_combine = int(sys.argv[7]) # Factor by which to reduce ns by joining bins

jes_prefix = "combined"#"samesite"

ne_new = ne_old/ne_combine
if (ne_new-int(ne_new)) != 0:
	raise Exception("ne_old must be divisible by ne_combine.")
ne_new = int(ne_new)

ns_new = ns_old/ns_combine
if (ns_new-int(ns_new)) != 0:
	raise Exception("ns_old must be divisible by ns_combine.")
ns_new = int(ns_new)

jes_comb_path = "jes_"+jes_prefix
loadfile = "jes-%s-w%s-ne%s-ns%s-n%s-nc%s.pkl" % (jes_prefix, w, ne_old, ns_old, Ntot, nconfigs)
with open(jes_comb_path+"/"+loadfile, "rb") as infile:
	jes_old, emin, emax, smin, smax, num_deltas_old, dos_old = pickle.load(infile)

jes_intermediate = np.zeros((ne_new,ns_old))
num_deltas_intermediate = np.zeros((ne_new,ns_old))
jes_new = np.zeros((ne_new,ns_new))
num_deltas_new = np.zeros((ne_new,ns_new))
dos_new = np.zeros((ne_new))
for e_index in np.arange(ne_new):
	for e_combine_index in np.arange(ne_combine):
		jes_intermediate[e_index,:] += jes_old[ne_combine*e_index+e_combine_index,:]
		num_deltas_intermediate[e_index,:] += num_deltas_old[ne_combine*e_index+e_combine_index,:]
		dos_new[e_index] += dos_old[ne_combine*e_index+e_combine_index]
for s_index in np.arange(ns_new):
	for s_combine_index in np.arange(ns_combine):
		jes_new[:,s_index] += jes_intermediate[:,ns_combine*s_index+s_combine_index]
		num_deltas_new[:,s_index] += num_deltas_intermediate[:,ns_combine*s_index+s_combine_index]
# There is no need to adjust normalization because the normalization by ne,ns is done in the Jupyter notebook

data = [jes_new, emin, emax, smin, smax, num_deltas_new, dos_new]
savefile = "jes-%s-w%s-ne%s-ns%s-n%s-nc%s.pkl" % (jes_prefix, w, ne_new, ns_new, Ntot, nconfigs)
with open(jes_comb_path+"/"+savefile, "wb") as outfile:
	pickle.dump(data, outfile, pickle.HIGHEST_PROTOCOL)

print("Done rebinning jes")
