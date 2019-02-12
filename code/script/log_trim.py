#!/usr/bin/env python
from subprocess import Popen, PIPE
import math, bisect
import sys, os
from liblrconf import *
from plotconst import *
from expconf import *

LOGDIR = '../log_notime/tmp/'
TOLERANCE = 1e-9
LOG_SUMMARY_LINES = 10
LINE_BUFFER = 10

def column_n_of_space(filename, f, fout):
	if fout is None:
		fout = PIPE
	p = Popen(['cut', '-d ', '-f'+str(f), filename], shell=False, stderr=PIPE, stdout=fout)
	res, err = p.communicate()
	if err:
		print(err)
		raise Exception('cut command error')
	if fout is None:
		if sys.version_info[0] >= 3: #for python 3
			res = res.decode()
		return res
	return None

def mv(fromfile, tofile):
	p = Popen(['mv', fromfile, tofile], shell=False, stderr=PIPE, stdout=PIPE)
	res, err = p.communicate()
	if err:
		raise Exception('mv command error')

def head(filename, n, fout):
	if fout is None:
		fout = PIPE
	p = Popen(['head', '-n'+str(n), filename], shell=False, stderr=PIPE, stdout=fout)
	res, err = p.communicate()
	if err:
		raise Exception('head command error')
	if fout is None:
		if sys.version_info[0] >= 3: #for python 3
			res = res.decode()
		return res
	return None

def tail(filename, n, fout):
	if fout is None:
		fout = PIPE
	p = Popen(['tail', '-n'+str(n), filename], shell=False, stderr=PIPE, stdout=fout)
	res, err = p.communicate()
	if err:
		raise Exception('head command error')
	if fout is None:
		if sys.version_info[0] >= 3: #for python 3
			res = res.decode()
		return res
	return None

def last_n_line_of_file(filename, n):
	p = Popen(['tail', '-'+str(n), filename], shell=False, stderr=PIPE, stdout=PIPE)
	res, err = p.communicate()
	if err:
		raise Exception('tail command error')
	if sys.version_info[0] >= 3: #for python 3
		res = res.decode()
	return res.split(os.linesep)[0]

def line_number(filename):
	p = Popen(['wc', '-l', filename], shell=False, stderr=PIPE, stdout=PIPE)
	res, err = p.communicate()
	if err:
		raise Exception('wc command error')
	if sys.version_info[0] >= 3: #for python 3
		res = res.decode()
	return int(res.split(' ')[0])

def is_log_format_good(filename):
	if line_number(filename) == 0 or \
		last_n_line_of_file(filename, LOG_SUMMARY_LINES) == '' or \
		last_n_line_of_file(filename, LOG_SUMMARY_LINES + 1) != '' or \
		last_n_line_of_file(filename, LOG_SUMMARY_LINES + 2) == '' :
			return False
	return True

def filenames_in_dir(dirname):
	p = Popen(['ls', dirname], shell=False, stderr=PIPE, stdout=PIPE)
	res, err = p.communicate()
	if err:
		raise Exception('ls command error')
	if sys.version_info[0] >= 3: #for python 3
		res = res.decode()
	return res.split(os.linesep)[:-1]

def get_enough_line_number(log, filepath):
	with open(filepath) as f:
		objs = map(float, f.read().splitlines())
	min_obj = min(objs)
	max_obj = max(objs)
	try:
		minimal = log.get_obj_minimal()
	except KeyError:
		sys.stderr.write("create minimal: dobj[%s][%s]: %.15g," % (log.get_dobj_key(), log.data, min_obj) + os.linesep)
		return len(objs)
	min_enough = minimal + math.fabs(minimal * TOLERANCE)

	# exception handling
	if min_enough >= max_obj and min_obj >= minimal:
		sys.stderr.write("error: tolerance %g is too large for %s" % (TOLERANCE, os.path.basename(filepath)) + os.linesep)
		return len(objs)
	if float("%.15g"%minimal) > float("%.15g"%min_obj):
		sys.stderr.write("update minimal: dobj[%s][%s]: %.15g, not %.15g"
				% (log.get_dobj_key(), log.data, min_obj, minimal) + os.linesep)
		return len(objs)

	# normal procedure
	objs_r = objs[::-1]
	min_enough_line = len(objs_r) - bisect.bisect_left(objs_r, min_enough)
	return min(min_enough_line + LINE_BUFFER, len(objs_r))

def cut_obj_column_to(filepath, outpath):
	tmp = outpath+'.cut_obj_column_to.tmp'
	with open(filepath) as f:
		line = f.readline().strip()
		tokens = line.split(' ')
		obj_index = tokens.index('obj') + 2
	with open(tmp, 'w') as f:
		column_n_of_space(filepath, obj_index, f)
	with open(outpath, 'w') as f:
		head(tmp, -LOG_SUMMARY_LINES-1, f)
	os.remove(tmp)

def exit_with_help():
	print("USAGE: " + sys.argv[0] + " filepath")
	exit(1)

def main():
	if '-h' in sys.argv:
		exit_with_help()
	filepaths = sys.argv[1:]
	for filepath in filepaths:
		try:
			filename = os.path.basename(filepath)
			filetmp = filename+'.tmp'
			cut_obj_column_to(filepath, filetmp)
			if not is_log_format_good(filepath):
				raise Exception('log not good: '+filepath)
			line_enough = get_enough_line_number(LogInfo(filepath), filetmp)
			with open(filetmp, 'w') as f:
				head(filepath, line_enough, f)
			with open(filetmp, 'a') as f:
				tail(filepath, LOG_SUMMARY_LINES+1, f)
			mv(filetmp, filepath)
		except Exception as e:
			sys.stderr.write('error: ' + filename+os.linesep)
			sys.stderr.write(str(e)+os.linesep)

if __name__ == '__main__':
	main()
