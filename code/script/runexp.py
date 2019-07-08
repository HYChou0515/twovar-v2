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
RESUME_PATH=ROOT_PATH+resumefolder+"/"
train="train"
grid="grid"

def exit_with_help():
	print("USAGE: " + sys.argv[0] + " [option]")
	print("option: [run,print,grid]")
	exit(1)

if len(sys.argv) != 2:
	exit_with_help()
mode = -1
ind = 0
if sys.argv[1] == "run":
	mode = 0
elif sys.argv[1] == "print":
	bashfile_name="out_run.bash"
	bashfile = open(bashfile_name, 'w')
	mode = 1
elif sys.argv[1] == "grid":
	bashfile_name="out_run.bash"
	bashfile = open(bashfile_name, 'w')
	gridjob_name_format="grid_%s.job.%d"
	mode = 2
else:
	exit_with_help()

for data_count, data in enumerate(dataset):
	print("Running "+ data)
	if mode == 2:
		process_count = PROCESS_MAX
		gridjob_count = 0
		gridjob_name = gridjob_name_format % (data, gridjob_count)
		gridjob_file = None
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
		if is_oneclass(runs[tp])==1:
			need_n = True
			nlist_real = nlist
			clist_real = [1]
		else:
			need_n = False
			nlist_real = [0.1]
			clist_real = clist
		kk = itertools.product(clist_real, elist_real, nlist_real, rlist_real)
		datapath = "%s%s" % (DATA_PATH, data)
		for tc, e, n, r in kk:
			basename = "%s_s%d_c%g_e%g" %(data, runs[tp], tc, e)
			if need_n:
				basename = basename + "_n%g" %(n)
			if need_r:
				basename = basename + "_r%g" %(r)
			filename = "%s%s" % (LOG_PATH, basename)
			resumename = "%s%s" % (RESUME_PATH, basename)
			print(filename)
			loginfo = LogInfo(filename)
			#	try:
			#		opt_val = loginfo.get_obj_minimal()
			#		opt_val += math.fabs(tolerance * opt_val)
			#	except Exception as e:
			#		sys.stderr.write('error: ' + filename+os.linesep)
			#		sys.stderr.write(str(e)+os.linesep)
			#		continue
			#cmd = "%s%s -s %d -c %g -e %g -m %d -t %d -o %.16g" % (ROOT_PATH, train, runs[tp], tc, e, m, timeout, opt_val)
			cmd = "%s%s -s %d -c %g -e %g -m %d -t %d" % (ROOT_PATH, train, runs[tp], tc, e, m, timeout)
			if need_n:
				cmd = cmd + " -n %g" % (n)
			if need_r:
				cmd = cmd + " -r %g" % (r)
			cmd = cmd+ " %s" % datapath
			print(str(ind) + ': ' + cmd)
			if mode == 0:
				cmd = cmd + " %s %s" % (filename, resumename)
				subprocess.Popen(cmd.split(' '))
				ind += 1
			elif mode == 1:
				cmd = cmd+ " > %s\n" % (filename)
				ind += 1
				bashfile.write(cmd)
			elif mode == 2:
				grid_opt = ' '.join(cmd.split(' ')[1:-1]) + " %s %s\n" % (filename, resumename)
				ind += 1
				if not gridjob_file or gridjob_file.closed:
					gridjob_file = open(gridjob_name, 'w')
				gridjob_file.write(grid_opt)
				process_count -= 1
				if process_count <= 0:
					gridjob_count += 1
					process_count = PROCESS_MAX
					gridjob_file.close()
					gridjob_name = gridjob_name_format % (data, gridjob_count)

	if mode == 2:
		if gridjob_file and not gridjob_file.closed:
			gridjob_file.close()
			gridjob_count += 1
		for i in range(gridjob_count):
			bashfile.write("%s%s %s < %s && echo \"job %d finished\"\n" % (ROOT_PATH, grid, datapath, gridjob_name_format % (data, i), data_count*gridjob_count+i))

if mode == 1 or mode == 2:
	bashfile.close()
