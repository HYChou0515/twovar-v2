#!/usr/bin/env python
from liblrconf import *
from plotconst import *
from expconf import *
import os
import math
import subprocess
import itertools
import time
import sys
ROOT_PATH="../"
DATA_PATH="../../data/"
LOG_PATH=ROOT_PATH+logfolder+"/"
train="train"

def exit_with_help():
	print("USAGE: " + sys.argv[0] + " [option]")
	print("option: [run,print]")
	exit(1)

if len(sys.argv) != 2:
	exit_with_help()
mode = -1
if sys.argv[1] == "run":
	mode = 0
elif sys.argv[1] == "print":
	process_count = 1
	ind = 0
	bashfile_name="out_run.bash"
	bashfile = open(bashfile_name, 'w')
	mode = 1
else:
	exit_with_help()

for data in dataset:
	print("Running "+ data)
	for tp in runtype:
		if is_semigd(runs[tp]):
			need_r = True
			rlist_real = rlist
		else:
			rlist_real = [0]
			need_r = False
		if is_shrink(runs[tp]):
			need_e = True
			elist_real = elist
		else:
			elist_real = [0.1]
			need_e = False
		if is_oneclass(runs[tp]):
			need_n = True
			nlist_real = nlist
			clist_real = [1]
		else:
			need_n = False
			nlist_real = [0.1]
			clist_real = clist
		kk = itertools.product(clist_real, elist_real, nlist_real, rlist_real)
		for tc, e, n, r in kk:
			datapath = "%s%s" % (DATA_PATH, data)
			filename = "%s%s_s%d_c%g_e%g" %(LOG_PATH, data, runs[tp], tc, e)
			if need_n:
				filename = filename + "_n%g" %(n)
			if need_r:
				filename = filename + "_r%g" %(r)
			print(filename)
			loginfo = LogInfo(filename)
			#try:
			#	opt_val = loginfo.get_obj_minimal()
			#	opt_val += math.fabs(tolerance * opt_val)
			#except Exception as e:
			#	sys.stderr.write('error: ' + filename+os.linesep)
			#	sys.stderr.write(str(e)+os.linesep)
			#	continue
			#cmd = "%s%s -s %d -c %g -e %g -m %d -t %d -o %.16g" % (ROOT_PATH, train, runs[tp], tc, e, m, timeout, opt_val)
			cmd = "%s%s -s %d -c %g -e %g -m %d -t %d" % (ROOT_PATH, train, runs[tp], tc, e, m, timeout)
			if need_n:
				cmd = cmd + " -n %g" % (n)
			if need_r:
				cmd = cmd + " -r %g" % (r)
			cmd = cmd+ " %s" % datapath
			print(cmd)
			if mode == 0:
				with open(filename, "w") as logfile:
					p = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE)
					out, err = p.communicate()
					logfile.write(out)
			elif mode == 1:
				cmd = cmd+ " %d.model > %s\n" % (ind, filename)
				ind += 1
				bashfile.write(cmd)

if mode == 1:
	bashfile.close()
