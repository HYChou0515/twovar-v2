#!/bin/bash
DATA="../../../data/heart_scale"
TRAIN="../../train"
GRID="../../grid"
TEST_OPTS="grids"
CORES=8

LOG=".log"
GRID_LOG=".grid$LOG"
TRAIN_LOG=".train$LOG"

GRIDJOB="grid.job"
TRAINJOB="train.job"
rm -rf $GRIDJOB
rm -rf $TRAINJOB
while read opt; do
	NAME=""
	for token in $opt; do
		if [ ${token:0:1} == "-" ]; then
			token=${token:1}
		else
			token="${token}_"
		fi
		NAME="${NAME}${token}"
	done
	NAME=${NAME::-1}
	echo "$TRAIN $opt $DATA $NAME$TRAIN_LOG" >> $TRAINJOB
	echo "$opt $NAME$GRID_LOG" >> $GRIDJOB
done < $TEST_OPTS
parallel --jobs $CORES < $TRAINJOB 
$GRID $DATA < $GRIDJOB

for log in *"${LOG}"; do
	viewlog.bash -f 1,3- $log > $log.tmp
	mv $log.tmp $log
done

DIFF=""
ONLY_TRAIN=""
ONLY_GRID=""
PASS=""
for trainlog in *"${TRAIN_LOG}"; do
	if [ ! -f $trainlog ]; then
		continue
	fi
	base=${trainlog::-$(printf "$TRAIN_LOG" | wc -c)}
	gridlog="$base$GRID_LOG"
	if [ ! -f $gridlog ]; then
		ONLY_TRAIN="$ONLY_TRAIN $trainlog"
		continue
	fi
	if [ -n "$(diff $trainlog $gridlog)" ]; then
		DIFF="$DIFF $base"
		continue
	else
		PASS="$PASS $base"
	fi
	rm $gridlog $trainlog
done
for gridlog in *"${GRID_LOG}"; do
	if [ -f $gridlog ]; then
		ONLY_GRID="$ONLY_GRID $gridlog"
	fi
done

if [ -n "$PASS" ]; then
	for base in $PASS; do
		echo "[pass] $base"
	done
fi

if [ -n "$ONLY_TRAIN" ]; then
	for log in $ONLY_TRAIN; do
		>&2 echo "[only train] $log"
	done
fi

if [ -n "$ONLY_GRID" ]; then
	for log in $ONLY_GRID; do
		>&2 echo "[only grid] $log"
	done
fi

if [ -n "$DIFF" ]; then
	for base in $DIFF; do
		>&2 echo "[diff] $base"
	done
fi
