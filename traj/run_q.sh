#! /bin/bash

dirpath="//media/ashwin/ASH_DRIVE_3/ashwin_md/lane/Aug_Nov2025/Fd100"
qpath="$dirpath/msd_queue.dat"

count=0
while read -r -a line; do

	count=$((count+1))

	if [ $count -ne 1 ]; then

		Tau=${line[0]}
		frameStart=${line[1]}
		frameEnd=${line[2]}

		fpathI="$dirpath/tau_$Tau/traj2.cfg"
		fpathO="$dirpath/tau_$Tau/msd_A.dat"

		if [ -f $fpathI ]; then

			printf "\n\nProcessing file %d, Tau = %s\n\n" $((count-1)) $Tau

			g++ -w -g -std=c++17 msd.cpp analysis.cpp -o ./bin/msd.o
			./bin/msd.o $fpathI $fpathO $frameStart $frameEnd

		else
			echo "\nInput error! File $fpathI does not exist.\n"
			exit 1
		fi
	fi

done <$qpath

exit 0