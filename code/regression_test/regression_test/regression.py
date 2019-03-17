#!/usr/bin/env python

import sys, os, errno, re, filecmp, shutil, platform, struct, time
from threading import Thread, Event
from subprocess import *
from Queue import Queue
try:
	from tabulate import tabulate
	TABULATE_FOUND = True
except ImportError:
	TABULATE_FOUND = False

class Option:

	def __init__(self, argv):
		if len(argv) < 3:
			self.exit_with_help()
		elif 3 <= len(argv) <= 5:
			self.dir_a = argv[1]
			self.dir_b = argv[2]
			self.cmd_path = 'cmds'
			self.result_dir = 'results/'
			if len(argv) >= 4:
				self.cmd_path = argv[3]
			if len(argv) >= 5:
				self.result_dir = argv[4]
		else:
			self.exit_with_help()
		self.nr_workers = 8
		self.log_path = os.path.join(self.result_dir, 'log.txt')
		self.output_dir = os.path.join(self.result_dir, 'output')
		self.check_param()
		self.check_size()
		self.check_files()

	def check_param(self):
		for path in [self.dir_a, self.dir_b]:
			if not os.path.exists(path):
				sys.stderr.write('Directory {0} does not exist.\n'.format(path))
				self.exit_with_help()

	def check_size(self):
		size_a = os.path.getsize(self.dir_a)
		size_b = os.path.getsize(self.dir_b)
		if abs(size_a - size_b) > size_a/3:
			print('A big difference between sizes of two directories!')

	def check_files(self):
		cmp_rst = filecmp.dircmp(self.dir_a, self.dir_b);
		if len(cmp_rst.left_only) != 0:
			print("Only in raw dir:")
			print(cmp_rst.left_only)
			self.exit_with_help()
		if len(cmp_rst.right_only) != 0:
			print("Only in new dir:")
			print(cmp_rst.right_only)
			self.exit_with_help()

	def clean(self):
		for linear_dir in [self.dir_a, self.dir_b]:
			if os.path.exists(linear_dir):
				shutil.rmtree(linear_dir)

	def exit_with_help(self):
		print('''\
Usage: {0} old_dir new_dir cmd_path[="cmds"] result_dir[="results/"]'''.format(sys.argv[0]))
		sys.exit(1)

class Worker(Thread):
	global compiled

	def __init__(self, wid, job_queue, rst_queue, options, compiled):
		Thread.__init__(self)
		self.wid = wid
		self.tid = None
		self.job_queue = job_queue
		self.rst_queue = rst_queue
		self.linear_dirs = [options.dir_a, options.dir_b]
		self.nr_workers = options.nr_workers
		self.out_dirs = None
		self.compiled = compiled

	def mkdirs(self):
		for out_dir in self.out_dirs:
			try:
				os.makedirs(out_dir)
			except OSError as e:
				if e.errno != errno.EEXIST:
					raise

	def run(self):
		while True:
			item = self.job_queue.get(timeout=1000000)
			if isinstance(item,int):
				nr_dead_workers = item + 1
				if nr_dead_workers == self.nr_workers:
					self.rst_queue.put(None)
				else:
					self.job_queue.put(nr_dead_workers)
				break

			tid, setting = item[0], item[1]
			job_type = re.findall(r'type="(.*?)"', setting)[0]
			if job_type != 'linear_compile':
				self.compiled.wait()
			self.out_dirs =  [os.path.join(options.output_dir, '{0}.{1}'.format(tid,x)) for x in self.linear_dirs]
			self.tid = tid
			self.mkdirs()
			if job_type == 'linear_compile':
				self.linear_compile(setting)
				self.compiled.set()
			elif job_type == 'linear_basic':
				self.linear_basic(setting)
			elif job_type == 'linear_python':
				self.linear_python(setting)
			elif job_type == 'linear_matlab':
				self.linear_matlab(setting)
			elif job_type == 'linear_nomodel':
				self.linear_nomodel(setting)
			elif job_type == 'linear_win':
				self.linear_win(setting)
			elif job_type == 'linear_win_check':
				self.linear_win_check(setting)
			else:
				sys.stderr.write('Unkown job type {0}\n'.format(job_type))

	def run_one(self, cmd, setting):
		std_redirect = '>> {out_dir}/stdout 2>> {out_dir}/stderr'
		error = False

		for i in range(2):
			if platform.system() == 'Windows':
				std_redirect = '>> {out_dir}\stdout 2>> {out_dir}\stderr'
				sub_cmd = '{0} {1}'.format(cmd, std_redirect)
			else:
				sub_cmd = 'sh -c "{0}" {1}'.format(cmd, std_redirect)

			sub_cmd = sub_cmd.format(linear_dir=self.linear_dirs[i], out_dir=self.out_dirs[i])
			p = Popen(sub_cmd, shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
			p.communicate()
			if p.returncode:
				error = True

			logfile='{out_dir}/logfile'.format(out_dir=self.out_dirs[i])
			if os.path.isfile(logfile):
				tmp='{out_dir}/tmp'.format(out_dir=self.out_dirs[i])
				viewlog_cmd = 'sh -c "viewlog.bash -f 1,3- {logfile}" > {tmp} && mv {tmp} {logfile}'.format(logfile=logfile, tmp=tmp)
				p = Popen(viewlog_cmd, shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
				p.communicate()
				if p.returncode: error = True

		diffs=filecmp.dircmp(self.out_dirs[0], self.out_dirs[1]).diff_files
		msg = ''
		if error:
			msg = 'error '
		if len(diffs) != 0:
			msg += 'diff=[' + ','.join(str(diff) for diff in diffs) + ']'
		else:
			msg += 'pass'
		self.rst_queue.put((self.tid, msg, setting))

	def linear_compile(self, setting):
		if platform.system() == 'Windows':
			self.rst_queue.put((self.tid, 'skip', setting))
			return
		cmd = 'make -C {linear_dir}'
		self.run_one(cmd, setting)

	def linear_proto(self, trainer, predictor, setting):
		resume = os.path.join('{out_dir}', 'resume')
		logfile = os.path.join('{out_dir}', 'logfile')
		dataset = re.findall(r'dataset="(.*?)"', setting)[0]
		train_options = re.findall(r'train_options="(.*?)"', setting)[0]
		predict_options = re.findall(r'predict_options="(.*?)"', setting)[0]
		if not os.path.exists(dataset):
			raise Exception("dataset not exist")

		cmd = ' '.join([trainer, train_options, dataset, logfile, resume])
		self.run_one(cmd, setting)

	def linear_basic(self, setting):
		if platform.system() == 'Windows':
			self.rst_queue.put((self.tid, 'skip', setting))
			return os.path.join('{linear_dir}\windows', 'train.exe')
		self.linear_proto('{linear_dir}/train', '{linear_dir}/predict', setting)

	def linear_python(self, setting):
		if platform.system() == 'Windows':
			cmd = 'linear_python.py '+ os.path.join('{linear_dir}', 'python') + ' {out_dir}'
		else:
			cmd = './linear_python.py '+ os.path.join('{linear_dir}', 'python') + ' {out_dir}'
		self.run_one(cmd, setting)

	def linear_nomodel(self, setting):
		dataset = re.findall(r'dataset="(.*?)"', setting)[0]
		options = re.findall(r'train_options="(.*?)"', setting)[0]
		if platform.system() == 'Windows':
			cmd = ' '.join([os.path.join('{linear_dir}',os.path.join('windows','train.exe')), options, dataset])
		else:
			cmd = ' '.join([os.path.join('{linear_dir}','train'), options, dataset])
		self.run_one(cmd, setting)

	def linear_matlab(self, setting):
		cmd = 'matlab -nodesktop -nosplash -r "linear_matlab(\'{linear_dir}/matlab\', \'{out_dir}\'); exit"'.replace('"', r'\"')
		if platform.system() == 'Windows':
			cmd = 'matlab -nodesktop -nosplash -r "linear_matlab(\'{linear_dir}/windows\', \'{out_dir}\'); exit"'.replace('"', r'\"')
		self.run_one(cmd, setting)

	def linear_win(self, setting):
		if platform.system() != 'Windows':
			self.rst_queue.put((self.tid, 'skip', setting))
			return
		self.linear_proto(os.path.join('{linear_dir}\windows', 'train'), '{linear_dir}\windows\predict', setting)

	def linear_win_check(self, setting):
		if platform.system() != 'Windows':
			self.rst_queue.put((self.tid, 'skip', setting))
			return
		dirs = self.linear_dirs
		items = ['liblinear.dll', 'train.exe', 'predict.exe']
		items_to_check = []
		infos = []
		msg = 'pass'

		for linear_dir in dirs:
			for linear_item in items:
				items_to_check.append(open(linear_dir + '\windows\\' + linear_item, 'rb'))

		for index, item in enumerate(items_to_check):
			item.seek(60)
			header_offset = struct.unpack("<L", item.read(4))[0]
			item.seek(header_offset + 4)
			info = struct.unpack("<H", item.read(2))[0]
			infos.append(info)

		for i in xrange(0, len(items)):
			if infos[i] != infos[i + len(items)]:
				msg = 'diff'
		self.rst_queue.put((self.tid, msg, setting))

if __name__ == '__main__':

	options = Option(sys.argv)

	if(os.path.exists(options.result_dir)):
		shutil.rmtree(options.result_dir)
	os.makedirs(options.result_dir)
	log_file = open(options.log_path, "w")
	job_queue, rst_queue = Queue(0), Queue(0)

	job_queue.put((0,'type="linear_compile"'))
	for i, line in enumerate(open(options.cmd_path)):
		line = line.strip()
		if line.startswith('#') or line == '':
			continue
		job_queue.put((i+1, line))
	job_queue.put(0)

	compiled = Event()
	for i in range(options.nr_workers):
		Worker(i, job_queue, rst_queue, options, compiled).start()

	formatted_rst = []
	while True:
		rst = rst_queue.get(timeout=1000000)
		if rst is None: break
		formatted_rst.append(rst)
		tid, rst, setting = rst[0], rst[1], rst[2]
		str_out = '{0} {1} {2}'.format(str(tid).ljust(4), rst.ljust(15), setting)
		print(str_out)
	if TABULATE_FOUND:
		log_file.write(tabulate(formatted_rst, tablefmt="plain"))
	else:
		for rst in formatted_rst:
			tid, rst, setting = rst[0], rst[1], rst[2]
			str_out = '{0} {1} {2}'.format(str(tid).ljust(4), rst.ljust(15), setting)
			log_file.write(str_out+os.linesep)
	log_file.close()
