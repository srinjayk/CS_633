

rm ~/Assignment1/hostsimproved
rm ~/Assignment1/plot*.png
rm ~/Assignment1/out.csv
rm ~/Assignment1/out

python3 ~/UGP/eagle/monitor/monitord.py start

# sleep 1m
g++ -o ~/UGP/allocator/src/allocator_improved.out ~/UGP/allocator/src/allocator_improved.cpp 
echo "Entering home directory"
#cd ~
./UGP/allocator/src/allocator_improved.out 64 8
#echo "Waiting starts .."
# sleep 1m
# echo "Waiting ends .. "
# python3 UGP/eagle/monitor/monitord.py stop
cp ~/UGP/allocator/src/hostsimproved ~/Assignment1/

echo "Entering Assignment directory"
cd ~/Assignment1/
echo "Compilation begins . . "
mpicc -o src src.c -lm

pwd

touch out.csv
touch out

echo "p,dp,time" >> out.csv
echo "Execution starting .. "
for P in 16 25 36 49 64
do
	# touch $P
	echo "running for P = $P"
	for N in 256 1024 4096 16384 65536 262144 1048576
	do
		# echo $N
		# echo $P
		for m in 1 2 3 4 5
		do
			BOUT=$(mpirun -np $P -f hostsimproved ./src $N 50)
			# echo "${BOUT}"
			for A in ${BOUT}
			do
				echo $P,$N,$A >> out.csv
				echo $A >> out
			done
		done
	done
done

cd ~

python3 UGP/eagle/monitor/monitord.py stop

cd ~/Assignment1/

echo "Plotting begins .. "
sleep 1m
python3 plot.py

echo "Plotting ends .. "
