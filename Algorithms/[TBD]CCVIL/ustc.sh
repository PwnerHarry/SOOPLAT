for ((fun = 1; fun<=20; fun=$fun+1))
do
	for ((strat = 5; strat<=5; strat=$strat+1))
	do
		for dim in 2000 
		do
			sed "s/^functionToRun [0-9]*$/functionToRun $fun/" configure.ini > F$fun-S$strat-D$dim
			sed -i "s/^learnStrategy [0-9]*$/learnStrategy $strat/" F$fun-S$strat-D$dim
			sed -i "s/^dimension [0-9]*$/dimension $dim/" F$fun-S$strat-D$dim

			# check whether the modified result is correct or not
			cat F$fun-S$strat-D$dim | grep '^functionToRun [0-9]*$'
			cat F$fun-S$strat-D$dim | grep '^learnStrategy [0-9]*$'
			cat F$fun-S$strat-D$dim | grep '^dimension [0-9]*$'

			bsub -q serial -o out.F$fun-S$strat-D$dim -e err.F$fun-S$strat-D$dim ./main.out F$fun-S$strat-D$dim
			sleep 1
		done
	done
done

echo 'All job have been submitted successfully'
