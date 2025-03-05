import sys
import os
import re

w = float(sys.argv[1]) # Anderson disorder strength - identifies files to be processed
ne = int(sys.argv[2]) # Steps for integrating over e
ns = int(sys.argv[3]) # Resolution for plotting s
repshift = int(sys.argv[4]) # Shift to apply to rep numbers

jes_path = "jes_individual"
onlyfiles = [f for f in os.listdir(jes_path) if os.path.isfile(os.path.join(jes_path, f))]

n_renamed = 0
for file in onlyfiles:
	result = re.search("jes-w%s-ne%s-ns%s-(.*)-rep(.*).pkl" % (w, ne, ns), file)
	if result!=None: # Only rename matching entries
		etc = str(result.group(1))
		whichrep = int(result.group(2))
		newname = ("jes-w%s-ne%s-ns%s-%s-rep%s.pkl" % (w, ne, ns, etc, (whichrep+repshift)))
		os.rename(jes_path+"/"+file, jes_path+"/"+newname)
		n_renamed +=1

if n_renamed==1:
	print("Renamed %s file" % n_renamed)
else:
	print("Renamed %s files" % n_renamed)