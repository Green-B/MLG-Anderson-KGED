
echo "Starting time:"
date
echo "---"
echo "w = "$wdisorder
Ntot=$((9*$Nred*$Nred))
echo "Ntot = "$Ntot
echo "These are the reps starting with "$(($repspernode*$whichrepset))

##########################################
# ACTIVATE PARALLELIZATION PACKAGES HERE #
##########################################

# We will keep files in the working directory until we're done
home= # ENTER HOME DIRECTORY HERE
work= # ENTER WORKING DIRECTORY HERE
# Copy the necessary code to the working directory
cp $home/mlganderson-new-1-diag.py $work/mlganderson-new-1-diag.py
cp $home/mlganderson-new-2-calc-j.py $work/mlganderson-new-2-calc-j.py
cp $home/mlganderson-new-3-combine-j.py $work/mlganderson-new-3-combine-j.py

cd $work

# Clean out data from old runs if any is present
if [ -d $work/en_v2 ]; then
	rm -rf $work/en_v2
fi
if [ -d $work/jes_individual ]; then
	rm -rf $work/jes_individual
fi
if [ -d $work/jes_combined ]; then
	rm -rf $work/jes_combined
fi

# Do the diagonalization
python $work/mlganderson-new-1-diag.py $wdisorder $Nred 1 $repspernode
wait
echo "Done diagonalizing"

# Do the binning
python $work/mlganderson-new-2-calc-j.py $wdisorder 1000 1000
wait
echo "Done binning"
rm -rf $work/en_v2
wait

# Combine the jes over runs
python $work/mlganderson-new-3-combine-j.py $wdisorder 1000 1000
wait
echo "Done combining"
rm -rf $work/jes_individual
wait

# Make the number string appear to be a float if it has no decimal
if [[ $wdisorder != *.* ]]; then
	wdisorder=$wdisorder.0
fi
# Move the combined file to the network drive
if [ ! -d $home/jes_combined ]; then
	mkdir $home/jes_combined
fi
mv $work/jes_combined/jes-combined-w$wdisorder-ne1000-ns1000-n$Ntot-nc$repspernode.pkl $home/jes_combined/jes-combined-w$wdisorder-ne1000-ns1000-n$Ntot-nc$repspernode.repset$whichrepset.pkl
wait
rm -rf $work/jes_combined
wait

# Remove the working code
rm $work/mlganderson-new-1-diag.py
rm $work/mlganderson-new-2-calc-j.py
rm $work/mlganderson-new-3-combine-j.py

echo "Ending time:"
date

