for P in 16 36 49 64
do
	touch $P
	for N in 256 1024 4096 16384 65536 262144 1048576
	do
		# echo $N
		# echo $P
		for m in 1 2 3 4 5
		do
			mpirun -np $P -hosts csews20,csews23,csews24,csews25,csews26,csews27,csews28,csews29,csews30 ./src $N 50 >> $P 
		done
	done
done
