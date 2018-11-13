#!/bin/bash
cp ../script/Makefile ./
make

for ((i = 1; i<=20; i=$i+1)) 
do
	for ((j = 5; j<=5; j=$j+1))
	do
		for k in 2000 
		do
			#	sed "s/^functionToRun [0-9]*$/functionToRun $i/" configure.ini > temp
			#	cat temp > configure.ini
			echo funcID $i, learnStrat $j, dimension $k
			#	cat configure.ini | grep '^functionToRun [0-9]*$'
			# ./submit.sh $i
			qsub ccvil.pbs -v fun=$i,strat=$j,dim=$k
			sleep 1
		done
	done
done

echo 'All job have been submitted successfully'
