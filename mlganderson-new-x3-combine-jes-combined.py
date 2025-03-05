import numpy as np
import pickle
import sys
import os
import re

# If two jes-combined files match in everything except the tags, then 
nc1 = int(sys.argv[1])
tag1 = sys.argv[2]
nc2 = int(sys.argv[3])
tag2 = sys.argv[4]
outputtag = sys.argv[5]
nctot = nc1+nc2

jes_comb_path = "jes_combined"
onlyfiles = [f for f in os.listdir(jes_comb_path) if os.path.isfile(os.path.join(jes_comb_path, f))]

for file1 in onlyfiles:
	result = re.search("jes-combined-(.*)-nc(.*)\.(.*).pkl", file1)
	if result!=None: # Only rename matching entries
		etc = str(result.group(1))
		this_nconfig = int(result.group(2))
		this_tag = str(result.group(3))
		# If two files with different tags and nconfigs but identical w, ne, ns exist, combine them
		if this_nconfig==nc1 and this_tag==tag1:
			try:
				file2 = "jes-combined-%s-nc%s.%s.pkl" % (etc, nc2, tag2)
				print(file1)
				print(file2)
				with open(jes_comb_path+"/"+file1, "rb") as infile:
					jes1, emin1, emax1, smin1, smax1, num_deltas1, dos1 = pickle.load(infile)
				with open(jes_comb_path+"/"+file2, "rb") as infile:
					jes2, emin2, emax2, smin2, smax2, num_deltas2, dos2 = pickle.load(infile)
				# Combine them and save the result
				if emin1==emin2 and emax1==emax2 and smin1==smin2 and smax1==smax2:
					jes_recb = jes1*nc1/nctot + jes2*nc2/nctot
					num_deltas_recb = num_deltas1*nc1/nctot + num_deltas2*nc2/nctot
					dos_recb = dos1*nc1/nctot + dos2*nc2/nctot
					savefile = "jes-combined-%s-nc%s.%s.pkl" % (etc, nc1+nc2, outputtag)
					data = [jes_recb, emin1, emax1, smin1, smax1, num_deltas_recb, dos_recb]
					with open(jes_comb_path+"/"+savefile, "wb") as outfile:
						pickle.dump(data, outfile, pickle.HIGHEST_PROTOCOL)
				else:
					raise Exception("Some jes-combined files do not have the same energy ranges! Check %s" % file1)
			except:
				pass

print("Done recombining jes")
