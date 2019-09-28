#!/bin/bash
function cpu_usage()
{
if [ -z "$1" ]
then
	NR_TOP_SMP="2"
else
	NR_TOP_SMP="$1"
fi
top -bn$NR_TOP_SMP | grep "Cpu(s)" | tail -1 | \
           sed "s/.*, *\([0-9]*\)[0-9.]*%* id.*/\1/" | \
           awk '{print 100 - $1}'
}

CPU_BOUND=40
MAX=10
counter=$MAX
while true
do
	cpu=$(cpu_usage 2)
	echo $cpu"%"
	if [ $cpu -le $CPU_BOUND ]
	then
		counter=$((counter-1))
		cpu=$(cpu_usage 3)
		echo $cpu"%"
		if [ $counter -le 0 ] && [ $cpu -le $CPU_BOUND ]
		then
			echo cpu usage too low
		#	cat code | ./send_mail.py s1243221@gmail.com mail.pub "CPU ($cpu%) usage is lower than $CPU_BOUND at $(date)"
			exit 0
		fi
	else
		counter=$MAX
	fi
	sleep 120
done
