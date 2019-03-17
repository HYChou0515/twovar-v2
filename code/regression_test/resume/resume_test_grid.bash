#!/bin/bash
if [ -z $1 ]
then 
	exit -1
else
	STYPE=$1
fi
SMRY=12

if [ "$(wc -l < $STYPE.2log)" -lt $(($SMRY+1)) ]; then
	>&2 echo "$STYPE $STYPE.2log empty"
	exit 1
fi
LOG2_ITER_START="$(head -1 $STYPE.2log | cut -d' ' -f 2)"
head -n$(($LOG2_ITER_START-1)) $STYPE.1log > $STYPE.1log.tmp
cat $STYPE.2log >> $STYPE.1log.tmp
viewlog.bash -f 1,3- $STYPE.1log.tmp > $STYPE.1logcmp
viewlog.bash -f 1,3- $STYPE.3log > $STYPE.3logcmp
if [ -n "$(diff $STYPE.1logcmp $STYPE.3logcmp)" ]
then
	>&2 echo $STYPE
fi
