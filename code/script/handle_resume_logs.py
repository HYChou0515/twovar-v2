#!/usr/bin/env python

import sys
from enum import Enum
import shutil

class State(Enum):
	FIND_FIRST_SUMMARY = 1
	FIND_SECOND_LOG = 2
	MERGE = 3

def merge_two(fname):
	# merge two logs (in one file) to one log
	state = State.FIND_FIRST_SUMMARY
	second_log_linenum = None
	ignore_start_linenum = None
	with open(fname) as f:
		for linenum, line in enumerate(f):
			if state is State.FIND_FIRST_SUMMARY and line.startswith('iter '):
				continue
			elif (state is State.FIND_FIRST_SUMMARY and not line.startswith('iter ')) or \
				(state is State.FIND_SECOND_LOG and not line.startswith('iter ')):
				state = State.FIND_SECOND_LOG
			elif state is State.FIND_SECOND_LOG and line.startswith('iter '):
				state = State.MERGE
				second_log_linenum = linenum
				break
			else:
				sys.stderr.write("State machine error\n")
				return -2

	if state is State.FIND_FIRST_SUMMARY:
		sys.stderr.write("Error: broken log\n")
		return -1
	if state is State.FIND_SECOND_LOG:
		return 0
	if state is not State.MERGE:
		sys.stderr.write("State machine error\n")
		return -2

	#state is State.MERGE
	target_iter = int((line.split(' ')[1]))
	with open(fname) as f:
		at_least_oneline = False
		for linenum, line in enumerate(f):
			if line.startswith('iter '):
				at_least_oneline = True
				check_iter = int((line.split(' ')[1]))
				if check_iter >= target_iter:
					ignore_start_linenum = linenum
					break
			elif at_least_oneline and not line.startswith('iter '):
				ignore_start_linenum = linenum-1
			else:
				sys.stderr.write("Error: broken log\n")
				return -1

	print('second_log_linenum = %d' % second_log_linenum)
	print('ignore_start_linenum = %d' % ignore_start_linenum)
	# ignore line such that ignore_start_linenum <= linenum < second_log_linenum
	if ignore_start_linenum is None or second_log_linenum is None:
		sys.stderr.write("State machine error\n")
		return -2
	fname_tmp = fname+'.tmp'
	with open(fname) as f, open(fname_tmp, 'w') as ftmp:
		for linenum, line in enumerate(f):
			if ignore_start_linenum <= linenum < second_log_linenum:
				# ignoring
				pass
			else:
				ftmp.write(line)
	shutil.move(fname_tmp, fname)
	return 1

def main():
	if len(sys.argv) <= 1:
		sys.stderr.write("Usage: %s logpath ..\n" % sys.argv[0])
		sys.exit(-1)
	for fname in sys.argv[1:]:
		trial = 1
		while True:
			print('%s with trial=%d' % (fname, trial))
			trial += 1
			exit_code = merge_two(fname)
			if exit_code == 0:
				break
			elif exit_code == 1:
				continue
			else:
				sys.stderr.write('Error Occur in file: %s' % fname)
				break

if __name__=='__main__':
	main()
