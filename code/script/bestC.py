from config import *
import subprocess
import itertools
rangeC = [2**(k) for k in range(-8,9)]

ROOT_PATH="../"
DATA_PATH="../../../data/"
LOG_PATH= ROOT_PATH + "L2bestc_log/"
for data in dataset:
	print("data:"+data)
	maxc = -1000
	maxr = -1
	for c in rangeC:
		filename = LOG_PATH + data + "_" + str(c)
		f = open(filename, 'r')
		line = f.readline()
		num = float(line.split(' ')[4].strip('\n%'))
		if num > maxr:
			maxr = num
			maxc = c
	print(maxc,maxr)

		
