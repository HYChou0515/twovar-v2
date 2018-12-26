#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd $DIR

PLOT=./plot.py
OLD=log_old
NEW=log_new
TIME=.time
CD=.cd

cd ..
$PLOT time > $DIR/$NEW$TIME
if [ $? -gt 0 ]
then
	echo time error
	exit
fi
$PLOT cd > $DIR/$NEW$CD
if [ $? -gt 0 ]
then
	echo cd error
	exit
fi

cd $DIR

diff ../config.py config.py > config.diff
if [ $(wc -l < config.diff) -gt 0 ]
then
	cp ../config.py ./
	$PLOT time > $OLD$TIME
	$PLOT cd > $OLD$CD
fi

diff $OLD$TIME $NEW$TIME > diff.time
diff $OLD$CD $NEW$CD > diff.cd

echo test finished

if [ $(wc -l < diff.time) -gt 0 ]
then
	echo time is different
	echo see diff.time
fi

if [ $(wc -l < diff.cd) -gt 0 ]
then
	echo cd is different
	echo see diff.cd
fi
