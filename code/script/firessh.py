#!/usr/bin/env python
#
# Rong-En Fan, 2006
#
# This is derived from grid.py in LIBSVM package. This script was first
# written by Chih-Jen Lin. Rong-En Fan modified it to be more general.
#
# This file is distributed under the same license as LIBSVM package.
#
# ChangeLogs
#  2009.01.04 (rafan): use list and execvp for call() to avoid meta char problem in shell

# vim:ts=4:sw=4:expandtab

import os, sys, traceback
import Queue
import getpass
import re
from threading import Thread
from string import find, split, join
from time import time, sleep
from subprocess import *
from random import shuffle

class Job:
    def __init__(self, no, cmd, fout_name):
        self.no = str(no)
        self.cmd = str(cmd)
		self.fout_name = fout_name
    def __str__(self):
        return "[%s] %s" % (self.no, self.cmd)
    def get(self):
        return (self.cmd, self.fout_name)

class WorkerStopToken:  # used to notify the worker to stop
        pass

class Worker(Thread):
    def __init__(self,name,job_queue,result_queue):
        Thread.__init__(self)
        self.name = name
        self.job_queue = job_queue
        self.result_queue = result_queue
    def run(self):
        while 1:
            job = self.job_queue.get()
            if job is WorkerStopToken:
                self.job_queue.put(job)
                #print 'worker %s stop.' % self.name
                break
            try:
                (job1, during) = self.run_one(job.get())
            except:
                # we failed, let others do that and we just quit
                traceback.print_exception(sys.exc_type, sys.exc_value, sys.exc_traceback)
                self.job_queue.put(job)
                print 'worker %s quit.' % self.name
                break
            else:
                self.result_queue.put((self.name, job, during))

class SSHWorker(Worker):
    def __init__(self,name,job_queue,result_queue,host, relcwd):
        Worker.__init__(self,name,job_queue,result_queue)
        self.host = host
        self.relcwd = relcwd
    def run_one(self,job):
		cmd, fout_name = job
        start = time()
		run_cmds = ["ssh", "-x", "%s" % self.host, "cd ~/%s; %s" % (self.relcwd, cmd)]
		if fout_name is not None and len(fout_name) != 0:
			with open(fout_name, 'w') as fout:
				call(run_cmds, shell = False, stdout=fout)
		else:
			call(run_cmds, shell = False)
        during = time() - start
        return (job, during)

class Grid:
    def __init__(self, workers, jobs):
        self.workers = workers
        self.jobs = []
        for i in range(len(jobs)):
			if '>' in jobs[i]:
				cmd,fout_name = jobs[i].split('>')
	            self.jobs.append(Job(i, cmd, fout_name.strip()))
			else:
	            self.jobs.append(Job(i, jobs[i], None))

        # put jobs in queue
        self.job_queue = Queue.Queue(0)
        self.result_queue = Queue.Queue(0)

        for i in range(len(self.jobs)):
            self.job_queue.put(self.jobs[i])

        # hack the queue to become a stack --
        # this is important when some thread
        # failed and re-put a job. If we still
        # use FIFO, the job will be put
        # into the end of the queue, and the graph
        # will only be updated in the end

        def _put(self,item):
            if sys.hexversion >= 0x020400A1:
                self.queue.appendleft(item)
            else:
                self.queue.insert(0,item)

        import new
        self.job_queue._put = new.instancemethod(_put, self.job_queue, self.job_queue.__class__)

    def go(self):
        elapsed_time = 0

        cwd = os.getenv('PWD')
        if not os.path.exists(cwd):
            cwd = os.getcwd()
		home = os.getenv('HOME')
		if os.path.commonprefix([cwd, home]) != home:
			sys.stderr.write('current not in home dir\n')
			exit(-1)
		relpath = os.path.relpath(cwd, home)

        # fire ssh workers
        for host in self.workers:
			sleep(0.1)
            SSHWorker(host, self.job_queue, self.result_queue, host, relpath).start()

        # gather results

        done_jobs = {}

        for job in self.jobs:
            while not done_jobs.has_key(job):
                (worker, job1, during) = self.result_queue.get()
                elapsed_time = elapsed_time + during
                done_jobs[job1] = job1
                print "[%s] (%s) %s" % (worker, during, job1.get())

        print "Elapsed Time (all workers): %s" % elapsed_time

        self.job_queue.put(WorkerStopToken)

hosts = [
	"r05922082@linux1",
	"r05922082@linux2",
	"r05922082@linux3",
	"r05922082@linux4",
	"r05922082@linux5",
	"r05922082@linux6",
	"r05922082@linux7",
	"r05922082@linux8",
	"r05922082@linux9",
	"r05922082@linux10",
	"r05922082@linux11",
	"r05922082@linux12",
	"r05922082@linux13",
	"r05922082@linux14",
	"r05922082@linux15",
]
host_cores = [
	32,32,32,32,24,24,24,24,24,24,24,24,24,24,16
#	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
]

if __name__ == '__main__':
	if len(sys.argv) != 2:
		sys.stderr.write('Usage: ' + sys.argv[0] + ' jobs\n')
		exit(-1)
	jobname = sys.argv[1]
	with open(jobname) as f:
		jobs = f.read().strip().split('\n')
	ssh_workers = []
	assert len(host_cores) == len(hosts)
	for i in range(len(hosts)):
		for j in range(host_cores[i]):
			ssh_workers.append(hosts[i])
	shuffle(jobs)
	shuffle(ssh_workers)
	grid = Grid(ssh_workers, jobs)
	grid.go()
