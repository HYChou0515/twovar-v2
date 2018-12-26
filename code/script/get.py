import math
import sys
import os
import numpy as np
from config import *

LOGPATH = '../logpatty/'
stype = [runs[k] for k in runtype]
if stype[0] in L1:
	loss = "L1"
else:
	loss = "L2"

def getlog(data,s,c,e):
	activesize = 0
	notup = 0
	logname = "%s_s%d_c%g_e%g"%(data, s, c, e)
	print(logname)
	fl = open(LOGPATH+logname,'r')
	if s in [32,35]:
		upd = 'udp'
	else:
		upd = 'upd'
	for line in fl:
		line = line.split(' ')
		if upd in line:
			iteration = float(line[line.index('iter')+1])
			if iteration > 200:
				break
			elif iteration < 0:
				continue
			activesize = activesize + float(line[line.index('actsize')+1])
			notup = notup + float(line[line.index(upd)+1])
	if s in [32,35]:
		print(data, s, notup/activesize)
	else:
		print(data, s, (activesize-notup)/activesize)
for data in dataset:
	for s in stype:
		for c in clist:
			for e in elist:
				getlog(data,s,c,e)
