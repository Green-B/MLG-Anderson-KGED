import numpy as np
import pickle
import sys
import os
import time

w = float(sys.argv[1]) # Anderson disorder strength - identifies files to be processed
ne = int(sys.argv[2]) # Steps for integrating over e
ns = int(sys.argv[3]) # Resolution for plotting s
# Get a list of everything in the path and pick out only the ones that are files
ev_path = "en_v2"
jes_path = "jes_individual"
if not os.path.exists(jes_path):
	os.makedirs(jes_path)
onlyfiles = [f for f in os.listdir(ev_path) if os.path.isfile(os.path.join(ev_path, f))]
# Load minimum and maximum energies, and be sure there is only one such file
search = "e-min-max-w%s" % w
thisw_minmax_file_list = [k for k in onlyfiles if search in k]
if (len(thisw_minmax_file_list) !=1):
	raise Exception("There should only be one e-min-max file for each w.")
thisw_minmax_file = thisw_minmax_file_list[0]
with open(ev_path+"/"+ thisw_minmax_file, "rb") as infile:
	emin, emax = pickle.load(infile) 
# Set minimum and maximum epsilon and s
ewidth = emax-emin
smin = 0.0+np.finfo(float).eps
smax = 1.0
swidth = smax-smin

# Load each file with energies and velocity matrix elements and calculate J(e,s)
search = "e-v-w%s" % w
thisw_ev_files = [k for k in onlyfiles if search in k]
nfiles = len(thisw_ev_files)
whichfile = 1
print("Number of files to process: %d" % nfiles)
for file in thisw_ev_files:
	print("Now processing "+file)
	#n_jes_included = 0 # Count the total number of jes included in the range [smin, smax] for this file
	#n_jes_tot = 0 # Count the total number of jes possible for this file
	t0 = time.time()
	with open(ev_path+"/"+file, "rb") as infile:
		# Calculate jes for each run and save it separately
		en, v2, Ntot = pickle.load(infile)
		[enen, emem] = np.meshgrid(en, en)
		ss = enen - emem
		ee = ( enen + emem )/2
		v2 = v2.flatten()
		ee = ee.flatten()
		ss = ss.flatten()
		# Remove deltas with s<smin; removal from ss must come last
		v2 = v2[ ss>=smin ]
		ee = ee[ ss>=smin ]
		ss = ss[ ss>=smin ]
		# Remove deltas with s>smax; removal from ss must come last
		v2 = v2[ ss<=smax ]
		ee = ee[ ss<=smax ]
		ss = ss[ ss<=smax ]
		[jes, xedges, yedges] = np.histogram2d(ee, ss, bins=[ne,ns], range=[[emin,emax],[smin,smax]], weights=v2)
		[num_deltas, xedges, yedges] = np.histogram2d(ee, ss, bins=[ne,ns], range=[[emin,emax],[smin,smax]])
		# Also find the density of states, using ne bins for consistency
		[dos, edges] = np.histogram(en, bins=ne, range=(emin,emax))
		# Save the data in files with the same name except with the prefix "jes", and ne, ns noted
		savefile = ("jes-w%s-ne%s-ns%s" % (w, ne, ns))+file[file.index("-n"):]
		data = [jes, emin, emax, smin, smax, num_deltas, dos]
		with open(jes_path+"/"+savefile, "wb") as outfile:
			pickle.dump(data, outfile, pickle.HIGHEST_PROTOCOL)
	t = time.time()-t0
	ts = t % 60
	tm = int((t-ts)/60 % 60)
	th = int((t-tm*60-ts)/3600 % 60)
	print("Finished processing %s of %s in time %d:%d:%.1f" % (whichfile, nfiles, th, tm, ts))
	whichfile += 1

print("Done processing all files")
