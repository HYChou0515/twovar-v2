#!/bin/bash
function cpu_usage()
{
top -bn2 | grep "Cpu(s)" | tail -1 | \
           sed "s/.*, *\([0-9]*\)[0-9.]*%* id.*/\1/" | \
           awk '{print 100 - $1}'
}

CPU_BOUND=50
MAX=10
counter=$MAX
while true
do
	sleep 300
	cpu=$(cpu_usage)
	echo $cpu"%"
	if [ $cpu -le $CPU_BOUND ]
	then
		counter=$((counter-1))
		if [ $counter -le 0 ]
		then
			echo cpu usage too low
			cat code | ./send_mail.py s1243221@gmail.com mail.pub "CPU ($cpu%) usage is lower than $CPU_BOUND at $(date)"
			exit 1
		fi
	else
		counter=$MAX
	fi
done
