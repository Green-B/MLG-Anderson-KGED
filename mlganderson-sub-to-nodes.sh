
# Do reps in sets of repspernode
wdisorder=$1
Nred=$2
repspernode=$3
howmanyrepsets=$4

for i in $( seq 1 $howmanyrepsets )
do
	qsub -pe impi 8 -v wdisorder=$wdisorder -v Nred=$Nred -v repspernode=$repspernode -v whichrepset=$i mlganderson-node-script.sh
done
