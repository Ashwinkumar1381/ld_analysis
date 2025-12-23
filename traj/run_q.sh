#! /bin/bash

dirpath="/home/ethaya_lab_c1/Desktop/ashwin_md/LD/lmp/psps/tau_1e-1"
qpath="$dirpath/collisions_queue.dat"

count=0
while read -r -a line; do

	count=$((count+1))

	if [ $count -ne 1 ]; then

		Fd=${line[0]}
		frameW=${line[1]}

		fpathI="$dirpath/Fd$Fd/traj2.cfg"
		fpathO="$dirpath/Fd$Fd/collisions.dat"

		if [ -f $fpathI ]; then

			printf "\n\nProcessing file %d, Fd = %d\n\n" $((count-1)) $Fd

			g++ -w -g -std=c++17 collisions.cpp analysis.cpp -o ./bin/coll.o
			./bin/coll.o $fpathI $fpathO $frameW

		else
			echo "\nInput error! File $fpathI does not exist.\n"
			exit 1
		fi
	fi

done <$qpath

exit 0