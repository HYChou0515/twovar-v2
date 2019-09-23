#!/usr/bin/env python
import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os
import numpy as np
import math

class Plotter(object):
	def __init__(self, args):
		self.file_names = [f.name for f in args.files]
		self.figname = args.figname

		self.xaxis = args.xaxis*2 -1
		self.xlabel = args.xaxis*2 -2
		self.xlim = args.xlim

		self.yaxis = args.yaxis*2 -1
		self.ylabel = args.yaxis*2 -2

		self.yrange = [None, None]
		if args.yrange is not None:
			try:
				self.yrange[0] = float(args.yrange[0])
			except Exception:
				pass
			try:
				self.yrange[1] = float(args.yrange[1])
			except Exception:
				pass

		self.args = args

	def draw(self):
		# plot setting
		self.fig = plt.figure()
		cmap = plt.get_cmap('hsv')
		NUM_COLOR = max(8, len(self.file_names))
		self.colors = [cmap(j) for j in np.linspace(0, 1, NUM_COLOR+1)]
		self.makr = ["--", "-.", "--", "o-.", "-", "o-"]
		self.makr = self.makr * int(math.ceil(float(len(self.colors))/len(self.makr)))

		# plotting
		for idx,fname in enumerate(self.file_names):
			self.draw_one_curve(idx, fname)

		# figure setup and save
		with open(self.file_names[0], 'r') as f:
			fields = f.readline().strip().split(' ')
			plt.xlabel(fields[self.xlabel])
			plt.ylabel(fields[self.ylabel])

		if self.args.xtype == 'int':
			ax = self.fig.gca()
			ax.xaxis.set_major_locator(MaxNLocator(integer=True))

		yrange_min, yrange_max = plt.ylim()
		if self.yrange[0] is not None:
			yrange_min = self.yrange[0]
		if self.yrange[1] is not None:
			yrange_max = self.yrange[1]
		plt.ylim(yrange_min, yrange_max)

		plt.legend(loc=0)
		self.fig.savefig(self.figname,format='png',dpi=100)

	def draw_one_curve(self, idx, fname):
		xs = self.get_column_of(fname, self.xaxis)
		ys = [y/1024/1024-0.5 for y in self.get_column_of(fname, self.yaxis)]
		#ys = self.get_column_of(fname, self.yaxis)
		xs, ys = self.condition_on_xy(xs, ys)
		plt.plot(xs, ys, self.makr[idx],
				color=self.colors[idx],
				figure=self.fig, label=os.path.basename(fname),
				linewidth=3, markersize=6, markevery=0.1)

	def get_column_of(self, fname, fieldid):
		nums = []
		with open(fname, 'r') as f:
			for line in f:
				fields = line.strip().split(' ')
				try:
					field = fields[fieldid]
					if(field[-1] == '%'):
						nums.append(float(field[:-1])/100)
					else:
						nums.append(float(field))
				except:
					break
		return nums

	def condition_on_xy(self, xs, ys):
		#xlim
		xlim_idx = 0
		while xlim_idx < len(xs) and xs[xlim_idx] < self.xlim:
			xlim_idx += 1
		xs = xs[:xlim_idx]
		ys = ys[:xlim_idx]

		return (xs, ys)

class Parser:
	def __init__(self):
		self.parser = argparse.ArgumentParser()

	def parse_option(self):
		self.parser.add_argument('--fig-name', dest='figname',
				type=str, action='store',
				default='figure.png',
				help='output figure name')
		self.parser.add_argument('--xlim', dest='xlim',
				type=float, action='store',
				default=float('Inf'),
				help='x limit')
		self.parser.add_argument('--xtype', dest='xtype',
				type=str, choices=['float', 'int'],
				default='float',
				help='type of x (shown in xticks)')
		self.parser.add_argument('--yrange', dest='yrange',
				nargs=2, metavar=('ymin_range', 'ymax_range'),
				default=None,
				help='range of y axis, \'None\' for not specified')
		#positional arguments
		self.parser.add_argument('xaxis', type=int,
				help='xaxis position(see viewlog)')
		self.parser.add_argument('yaxis', type=int,
				help='yaxis position(see viewlog)')
		self.parser.add_argument('files', type=argparse.FileType('r'),
				nargs='+', help='input files')
		return self.parser.parse_args()

def main():
	parser = Parser()
	args = parser.parse_option()

	plotter = Plotter(args)
	plotter.draw()

if __name__ == '__main__':
	main()
