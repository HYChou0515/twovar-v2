#!/bin/bash

read -r -a fields <<< "iter t obj decr_rate actsize sucsize nr_n_ops ops_per_sucs updsize sucs_rate cdsteps nSV nBSV nFree nNonSV PGmax PGmin PGdiff Gmax Gmin Gdiff n_exchange alpha_diff nr_pos_y nr_neg_y"

function exit_with_help
{
	echo "Usage: viewlog [options] logfile"
	echo "Options:"
	echo "-f field: shows only field"
	echo "-v: vim mode"
	echo "-o: horizontal split"
	for i in "${!fields[@]}"; do 
		printf "%s\t%s\n" "$(($i+1))" "${fields[$i]}"
	done
	exit 1
}


if [[ $# < 1 ]]; then
	exit_with_help
fi

files=()

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
	-f)
		cutfields="$2"
		shift # past argument
		shift # past value
		;;
	-v)
		vimmode="1"
		shift
		;;
	-o)
		spmode="1"
		shift
		;;
	*)# unknown option
		files+=("$1")
		shift
		;;
esac
done
TMPDIR="$HOME/.viewlog"
TMPSUFFIX=".tmp"
mkdir -p $TMPDIR

for file in "${files[@]}"; do
	if [ ! -f $file ]; then
		echo "no such file: $file"
		exit_with_help
	fi
done
for i in "${!files[@]}"; do
	file=${files[$i]}
	loglen=$(awk '$0==""{print NR-1}' $file)
	summarylen=$(($(wc -l < $file)-$loglen))
	if [ -z $vimmode ]; then
		out="/dev/stdout"
	else
		out="$TMPDIR/$(basename $file)$TMPSUFFIX$i"
		tmps+=("$out")
	fi
	if [ -z $cutfields ]; then
		head -${loglen} ${file} | column -t > $out
		tail -${summarylen} ${file} >> $out 
	else
		echo $file
		head -${loglen} ${file} | sed 's/\ \([^\ ]*\)\ /\ \1\|/g' | cut -d'|' -f ${cutfields} | sed 's/|/\ /g' | column -t > $out 
		tail -${summarylen} ${file} >> $out
	fi
done

if [ $vimmode ]; then
	if [ $spmode ]
	then
		vim -o ${tmps[*]}
	else
		vim -O ${tmps[*]}
	fi
fi
