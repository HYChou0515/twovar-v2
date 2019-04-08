#!/usr/bin/env python
import argparse
import matplotlib
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import math

class Plotter(object):
	def __init__(self, args):
		self.filename = args.filename
		self.xlim = args.xlim
                if args.outdir is not None:
                    self.outdir = args.outdir
                else:
                    self.outdir = os.path.basename(self.filename)

	def draw(self):
		# Create target Directory
		try:
			os.mkdir(self.outdir)
		except IOError:
			sys.exit('Directory: \''+self.outdir+'\' already existed')

		with open(self.filename) as f:
			line = f.readline().strip()
			start=False
			self.pts = []
                        self.avg_rand = []
                        self.idxs = []
			i = 0
			while True:
				line = f.readline().strip()
				if not start and line!='===':
					continue
				if not start and line == '===':
					start = True
					continue
				if start and line == '===':
					break
                                self.idxs.append(int(line))
				line = f.readline().strip()
                                self.avg_rand.append(float(line))
				line = f.readline().strip()
                                ys = map(lambda smgd: float(smgd)/self.avg_rand[-1], line.split(' '))
				if self.xlim < float('Inf'):
					ys = ys[:int(math.ceil(self.xlim))]
				self.pts.append(ys)
				i+=1
                print('start cauculating summary')
                with open(os.path.join(self.outdir, 'summary'), "w") as flog:
                    for i in range(len(self.pts)):
                        flog.write("%d %.3g %.3g\n" % (self.idxs[i], self.avg_rand[i], np.average(self.pts[i][i:])))
                print('start plotting')
                self.yticks=None
                self.ymax = None
                ymax_i = np.argmax([max(ys) for ys in self.pts])
                self.draw_one(ymax_i) # ymax_i first to set up the ymax and yticks
                for i in range(len(self.pts)):
                        if i == ymax_i: continue
                        sys.stdout.write("\r%d" % i)
                        sys.stdout.flush()
                        print(i)
                        self.draw_one(i)

	def draw_one(self, i):
		fig = plt.figure()
		ys = self.pts[i]
		xs = [ii for ii in range(len(ys))]
		plt.axvline(self.idxs[i], alpha=0.7, color='k', linestyle='-', linewidth=0.5)
		plt.axhline(1, alpha=0.7, color='r', linestyle='-', linewidth=0.5)
		plt.plot(xs, ys, 'o', alpha=0.7, markersize=1, figure=fig)
		plt.yscale('symlog')
		if self.yticks is None:
			self.yticks, _ = plt.yticks()
			_, self.ymax = plt.ylim()
		plt.ylim(-0.25, self.ymax)
		plt.yticks(self.yticks)
                basename = ('%d'%(self.idxs[i])).zfill(int(math.ceil(math.log10(self.idxs[-1]))+3))
		fig.savefig(os.path.join(self.outdir, basename+'.png'),
				format='png', dpi=100)
		plt.close(fig)

class Parser:
	def __init__(self):
		self.parser = argparse.ArgumentParser()

	def parse_option(self):
                #optional arguments
                self.parser.add_argument('--xlim', dest='xlim',
                                type=float, action='store',
                                default=float('Inf'),
                                help='x limit')
                self.parser.add_argument('--outdir', dest='outdir',
                                type=str, action='store',
                                default=None,
                                help='output directory')
		#positional arguments
		self.parser.add_argument('filename', type=str,
				action='store',
				help='input file')
		return self.parser.parse_args()

def main():
	parser = Parser()
	args = parser.parse_option()

	plotter = Plotter(args)
	plotter.draw()

if __name__ == '__main__':
	main()
