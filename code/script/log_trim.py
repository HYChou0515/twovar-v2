#!/usr/bin/env python
from subprocess import Popen, PIPE
import math, bisect
import sys, os
from itertools import islice
from liblrconf import *
from plotconst import *
from expconf import *
import traceback

LOGDIR = '../log_notime/tmp/'
TOLERANCE = 1e-9
LINE_BUFFER = 10
MAX_LINE = 10000
LOG_SUMMARY_LINES = None # assigned by find_log_summary_lines()

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

def find_log_summary_lines(filepath):
	total_line_num = line_number(filepath)
	if line_number(filepath) == 0:
		raise Exception('log not good: '+filepath)
	global LOG_SUMMARY_LINES
	LOG_SUMMARY_LINES = 1
	# find empty line (saperate empty line)
	while LOG_SUMMARY_LINES <= total_line_num:
		line = last_n_line_of_file(filepath, LOG_SUMMARY_LINES)
		if line == '':
			break
		LOG_SUMMARY_LINES += 1

	# find non-empty line (last line of log)
	while LOG_SUMMARY_LINES <= total_line_num:
		line = last_n_line_of_file(filepath, LOG_SUMMARY_LINES)
		if line != '':
			break
		LOG_SUMMARY_LINES += 1
	LOG_SUMMARY_LINES -= 1

def filenames_in_dir(dirname):
	p = Popen(['ls', dirname], shell=False, stderr=PIPE, stdout=PIPE)
	res, err = p.communicate()
	if err:
		raise Exception('ls command error')
	if sys.version_info[0] >= 3: #for python 3
		res = res.decode()
	return res.split(os.linesep)[:-1]

def get_enough_decrease_line_number(log, filepath):
	with open(filepath) as f:
		decr_rate = map(float, f.read().splitlines())
	i = 0
	buff = 0
	while i < len(decr_rate) and buff < LINE_BUFFER:
		if decr_rate[i] < 1e-15:
			buff += 1
		else:
			buff = 0
		i += 1
	return min(i, len(decr_rate))

def get_enough_optimal_line_number(log, filepath):
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

def cut_key_column_to(filepath, outpath, key):
	tmp = outpath+'.cut_' + key + '_column_to.tmp'
	with open(filepath) as f:
		line = f.readline().strip()
		tokens = line.split(' ')
		obj_index = tokens.index(key) + 2
	with open(tmp, 'w') as f:
		column_n_of_space(filepath, obj_index, f)
	with open(outpath, 'w') as f:
		head(tmp, -LOG_SUMMARY_LINES, f)
	os.remove(tmp)

def compress(filepath, filetmp, max_line):
	cmpr_line = int(math.ceil((line_number(filepath)-LOG_SUMMARY_LINES)*1.0/max_line))
        MAX_FILE_LINES = 4000000
	BUFF_READ_LINE = 4000
	SUM='sum'
	EXACT='exact'
	DECR_RATE='DR'
	SUCS_RATE='SR'
	class Aggregator(object):
		def __init__(self, colnames, cmpr_line):
			self.colnames = colnames
			self.column_type = {
					'iter': (EXACT, int, lambda x: "%d "%x),
					't': (EXACT, float, lambda x: "%f "%x),
					'obj': (EXACT, float, lambda x: "%.16g "%x),
					'decr_rate': (DECR_RATE, float, lambda x: "%.3e "%x),
					'actsize': (SUM, int, lambda x: "%d "%x),
					'updsize': (SUM, int, lambda x: "%d "%x),
					'sucpair': (SUM, int, lambda x: "%d "%x),
					'sucs_rate': (SUCS_RATE, lambda s: float(s[:-1])/100, lambda x: "%.2f%% "%(100*x)),
					'nSV': (EXACT, int, lambda x: "%d "%x),
					'nBSV': (EXACT, int, lambda x: "%d "%x),
					'nFree': (EXACT, int, lambda x: "%d "%x),
					'nNonSV': (EXACT, int, lambda x: "%d "%x),
					'PGmax': (EXACT, float, lambda x: "%.16g "%x),
					'PGmin': (EXACT, float, lambda x: "%.16g "%x),
					'PGdiff': (EXACT, float, lambda x: "%.3f "%x),
					'Gmax': (EXACT, float, lambda x: "%.16g "%x),
					'Gmin': (EXACT, float, lambda x: "%.16g "%x),
					'Gdiff': (EXACT, float, lambda x: "%.3f "%x),
					'n_exchange': (SUM, int, lambda x: '-1 ' if x<0 else "%d "%x),
					'alpha_diff': (SUM, float, lambda x: "%.16g "%x),
					'nr_pos_y': (EXACT, int, lambda x: "%d "%x),
					'nr_neg_y': (EXACT, int, lambda x: "%d "%x),
			}
			assert len(self.colnames) == len(self.column_type)
			self.cols = [None] * len(self.colnames)
			self.counter = 0
			self.cmpr_line = cmpr_line
			self.last_obj = float('nan')
			assert self.cmpr_line >= 1

		def append(self, vals):
			self.counter += 1
			for i in range(len(self.colnames)):
				colname = self.colnames[i]
				t, to_val, to_str = self.column_type[colname]
				val = to_val(vals[i])
				if self.cols[i] is None:
					self.cols[i] = val
				else:
					if t == SUM:
						self.cols[i] += val
					if t == EXACT:
						self.cols[i] = val
			if self.counter == self.cmpr_line:
				out_str = ''
				for i in range(len(self.colnames)):
					colname = self.colnames[i]
					out_str += "%s "%colname
					t, to_val, to_str = self.column_type[colname]
					if t == SUM:
						out_str += to_str(self.cols[i])
					if t == EXACT:
						out_str += to_str(self.cols[i])
					if t == DECR_RATE:
						this_obj = self.cols[self.colnames.index('obj')]
						out_str += to_str((self.last_obj-this_obj)/ math.fabs(this_obj))
						self.last_obj = this_obj
					if t == SUCS_RATE:
						sucpair = self.cols[self.colnames.index('sucpair')]
						updsize = self.cols[self.colnames.index('updsize')]
						out_str += to_str(sucpair/updsize)
				self.cols = [None] * len(self.colnames)
				self.counter = 0
				return out_str
			else:
				return None

	with open(filepath) as f:
		colnames = f.readline().strip().split(' ')[::2]
		aggregator = Aggregator(colnames, cmpr_line)

	def next_n_lines(f, n):
	    return [x.strip() for x in islice(f, n)]
	with open(filepath) as fin, open(filetmp, 'w') as fout:
                line_count = 0
		while line_count < MAX_FILE_LINES:
                        line_count += BUFF_READ_LINE
			lines = next_n_lines(fin, BUFF_READ_LINE)
			if len(lines) == 0:
				break
			for line in lines:
				if line == '':
					break
				vals = line.split(' ')[1::2]
				if len(vals) != len(aggregator.colnames):
					break
				out_str = aggregator.append(vals)
				if out_str is not None:
					fout.write(out_str + os.linesep)
                assert line_count < MAX_FILE_LINES, 'file too large'

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
			find_log_summary_lines(filepath)
			before_line_num = line_number(filepath)
			loginfo = LogInfo(filepath)
			filetmp = filename+'.tmp'
			#trim
			cut_key_column_to(filepath, filetmp, 'obj')
			line_opt_enough = get_enough_optimal_line_number(loginfo, filetmp)
			cut_key_column_to(filepath, filetmp, 'decr_rate')
			line_decr_enough = get_enough_decrease_line_number(loginfo, filetmp)
			line_enough = min(line_opt_enough, line_decr_enough)
			with open(filetmp, 'w') as f:
				head(filepath, line_enough, f)
			with open(filetmp, 'a') as f:
				tail(filepath, LOG_SUMMARY_LINES, f)
			mv(filetmp, filepath)
			#compress
			compress(filepath, filetmp, MAX_LINE)
			with open(filetmp, 'a') as f:
				tail(filepath, LOG_SUMMARY_LINES, f)
			mv(filetmp, filepath)

			after_line_num = line_number(filepath)
			print("%s\tbefore\t%d\tafter\t%d\ttrim\t%d" % (filename, before_line_num, after_line_num, before_line_num-after_line_num))
		except Exception as e:
			sys.stderr.write('error: ' + filename+os.linesep)
			_, _, tb = sys.exc_info()
			traceback.print_tb(tb) # Fixed format

if __name__ == '__main__':
	main()
