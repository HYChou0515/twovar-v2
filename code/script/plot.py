#!/usr/bin/env python
import math
import sys
import os
import numpy as np
from config import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import argparse
import itertools
def format_label(x,pos):
#	return "10$^{%.5f}$" % np.log10(x)
	return "%1.e" %x
def format_labelx(x,pos):
	return "%.1f" %x

def expo_slice_avg(l,base):
	assert base>1
	head = 0
	length = 1
	ret = []
	while head+length < len(l):
		part = l[head:head+int(length)]
		ret.append(sum(part)/float(len(part)))
		head += int(length)
		length = length*base
	return ret


class Plotter(object):
	def __init__(self, stype, loss, dataset, clist, eps, args):
		self.dataset = dataset
		self.clist = clist
		self.stype = stype
		self.loss = loss
		self.eps = eps
		self.logpath = args.logpath
		self.figpath = args.figpath

		matplotlib.rc('xtick', labelsize=20)
		matplotlib.rc('ytick', labelsize=20)

	def draw_all(self):
		for dstr,c in itertools.product(self.dataset, self.clist):
			self.dstr = dstr
			self.c = c
			self.draw()

	def draw(self):
		self.setup_color_marker()
		self.draw_fig()

	def setup_color_marker(self):
		self.fig = plt.figure()
		cmap = plt.get_cmap('hsv')
		totnum = len(self.stype)
		for st in self.stype:
			if st in semigd:
				totnum = totnum + len(rlist) -1
		self.colors = [cmap(j) for j in np.linspace(0, 1, totnum+1)]
		self.makr = ["--", "-.", "o-.", "-", "o-"]
		self.makr = self.makr * int(math.ceil(float(len(self.colors))/len(self.makr)))

	def init_new_fig(self):
		pass

	def draw_fig(self):
		self.init_new_fig()
		clridx = 0
		for tp in self.stype:
			getrlist = self.get_rlist(tp)
			e = self.get_eps(tp, self.eps)
			for r in getrlist:
				try:
					ys, xs = self.get_xy(tp, e, r)
				except IOError:
					continue

				lb = uselabel[tp]
				if tp in semigd:
					lb = "%s_%s" % (lb, r)
				plt.plot(xs, ys, self.makr[clridx],
						color=self.colors[clridx],
						figure=self.fig, label=lb,
						linewidth=3, markersize=6, markevery=0.1)
				clridx += 1
		self.setup_fig()

	def get_rlist(self, tp):
		if tp not in semigd:
			return [1]
		else:
			return rlist

	def get_eps(self, tp, orig_e):
		if tp not in shrink:
			return 0.1
		else:
			return orig_e

	def get_logfile(self, tp, e, r):
		if tp not in semigd:
			if tp in oneclass:
				logname = "%s_s%d_c1_e%g_n%g"% (self.dstr, tp, e, self.c)
			else:
				logname = "%s_s%d_c%g_e%g"% (self.dstr, tp, self.c, e)
		else:
			if tp in oneclass:
				logname = "%s_s%d_c1_e%g_n%g_r%g"% (self.dstr, tp, e, self.c, r)
			else:
				logname = "%s_s%d_c%g_e%g_r%g"% (self.dstr, tp, self.c, e, r)
		print(logname)
		try:
			return open(self.logpath+logname,"r")
		except IOError:
			print("cannot find file: "+self.logpath+logname)
			raise IOError

	def get_minimal(self, tp):
		if tp in biasobj:
			return dobj["bias{}c{}".format(self.loss,self.c)][self.dstr]
		elif tp in oneclass:
			return dobj["one{}c{}".format(self.loss,self.c)][self.dstr]
		else:
			return dobj["{}c{}".format(self.loss,self.c)][self.dstr]

	def setup_fig(self):
		plt.legend(loc=0)
		self.fig.savefig(self.figpath+self.get_figname_fmt()%(self.dstr, self.c),format='png',dpi=100)

class CdPlotter(Plotter):
	PLOTTYPE = "cd"
	YLIM = (1e-9, 1)
	def init_new_fig(self):
		self.min_y = 10
	def get_xlim(self, tp):
		key = "s%d_c%g_iter" % (tp, self.c)
		if key in dlim.keys():
			if self.dstr in dlim[key]:
				return dlim[key][self.dstr]
			else:
				return 100
		else:
			return 100
	def get_figname_fmt(self):
		sname = ''
		for s in self.stype:
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_cdsteps.png"
	def get_xy(self, tp, e, r):
		try:
			logfile = self.get_logfile(tp, e, r)
		except IOError:
			raise IOError
		ys = []
		xs = []
		CDsteps = 0
		minimal = self.get_minimal(tp)

		for line in logfile:
			line = line.split(' ')
			if 'iter' in line and 'obj' in line:
				val = float(line[line.index('obj')+1])
				relval = math.fabs((val - minimal)/minimal)
				if not (CdPlotter.YLIM[0] < relval < CdPlotter.YLIM[1]):
					break
				if CDsteps > 1e6:
					break
				ys.append(relval)
				self.min_y = min(relval,self.min_y)
				CDsteps = CDsteps + (float(line[line.index('updsize')+1]))
				xs.append(CDsteps)
				#print('CDsteps: %.16g\t relval: %.16g' % (CDsteps, relval))
		logfile.close()
		return ys, xs
	def setup_fig(self):
		plt.ylim(self.min_y,1)
		if(self.min_y > 1.0e-2):
			subsyy = [2,3,4,5,7]
		else:
			subsyy = []
		plt.yscale('log', subsy=subsyy, figure=self.fig)
		if(self.min_y > 1.0e-2):
			plt.tick_params(axis='y', which='minor', labelsize=14)
			plt.gca().yaxis.set_minor_formatter(FuncFormatter(format_label))
		else:
			plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=7))
		plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		plt.tick_params(axis='x', which='major', labelsize=20)
		Plotter.setup_fig(self)

class TimePlotter(Plotter):
	PLOTTYPE = "time"
	YLIM = (1e-9, 1)
	XLIM = (0, 100)
	def init_new_fig(self):
		self.min_y = 10
	def get_xlim(self, tp):
		key = "s%d_c%g_shrink" % (tp, self.c)
		if key in dlim.keys() and self.dstr in dlim[key]:
			return dlim[key][self.dstr]
		else:
			return 100

	def get_figname_fmt(self):
		sname = ''
		for s in self.stype:
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_time.png"
	def get_xy(self, tp, e, r):
		try:
			logfile = self.get_logfile(tp, e, r)
		except IOError:
			raise IOError
		ys = []
		xs = []
		minimal = self.get_minimal(tp)
		xlim = self.get_xlim(tp)

		for line in logfile:
			line = line.split(' ')
			if 't' in line and'obj' in line:
				val = float(line[line.index('obj')+1])
				relval = math.fabs((val - minimal)/minimal)
				t =  float(line[line.index('t')+1])
				if not (TimePlotter.YLIM[0] < relval < TimePlotter.YLIM[1]):
					break
				if not (TimePlotter.XLIM[0] < t < TimePlotter.XLIM[1]):
					break
				if t > xlim:
					break
				xs.append(t)
				ys.append(relval)
				self.min_y = min(relval,self.min_y)
				#print('t: %.16g\t relval: %.16g' % (t, relval))
		logfile.close()
		return ys, xs
	def setup_fig(self):
		plt.ylim(self.min_y,1)
		if(self.min_y > 1.0e-2):
			subsyy = [2,3,4,5,7]
		else:
			subsyy = []
		plt.yscale('log', subsy=subsyy, figure=self.fig)
		if(self.min_y > 1.0e-2):
			plt.tick_params(axis='y', which='minor', labelsize=14)
			plt.gca().yaxis.set_minor_formatter(FuncFormatter(format_label))
		else:
			plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=7))
		plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		plt.tick_params(axis='x', which='major', labelsize=20)
		Plotter.setup_fig(self)

class SucrCdPlotter(Plotter):
	PLOTTYPE = "sucrate-cd"
	YLIM = (1e-9, 1)
	def init_new_fig(self):
		self.max_y = 0
	def get_figname_fmt(self):
		sname = ''
		for s in self.stype:
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_sucrate.png"
	def get_xy(self, tp, e, r):
		try:
			logfile = self.get_logfile(tp, e, r)
		except IOError:
			raise IOError
		ys = []
		xs = []
		cdsteps = 0

		for line in logfile:
			line = line.split(' ')
			if 'iter' in line and 'obj' in line:
				val = float(line[line.index('sucpair')+1])/float(line[line.index('updsize')+1])
				ys.append(val)
				self.max_y = max(val,self.max_y)
				cdsteps = cdsteps + (float(line[line.index('updsize')+1]))
				xs.append(cdsteps)
		logfile.close()
		#slice xs and ys into average of 1,2,4,8,...
		xs = expo_slice_avg(xs, 1.05)
		ys = expo_slice_avg(ys, 1.05)
		return ys, xs
	def setup_fig(self):
		plt.ylim(0,self.max_y)
		plt.xscale('log', figure=self.fig)
		plt.tick_params(axis='x', which='major', labelsize=20)
		Plotter.setup_fig(self)
class ObjSucPlotter(Plotter):
	PLOTTYPE = "obj-suc"
	YLIM = (1e-9, 1)
	def init_new_fig(self):
		self.min_y = 10
	def get_figname_fmt(self):
		sname = ''
		for s in self.stype:
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_sucupd.png"
	def get_xy(self, tp, e, r):
		try:
			logfile = self.get_logfile(tp, e, r)
		except IOError:
			raise IOError
		ys = []
		xs = []
		sucpair = 0
		minimal = self.get_minimal(tp)

		for line in logfile:
			line = line.split(' ')
			if 'iter' in line and 'obj' in line:
				val = float(line[line.index('obj')+1])
				relval = math.fabs((val - minimal)/minimal)
				if not (ObjSucPlotter.YLIM[0] < relval < ObjSucPlotter.YLIM[1]):
					break
				ys.append(relval)
				self.min_y = min(relval,self.min_y)
				sucpair = sucpair + (float(line[line.index('sucpair')+1]))
				xs.append(sucpair)
				#print('sucpair: %.16g\t relval: %.16g' % (sucpair, relval))
		logfile.close()
		return ys, xs
	def setup_fig(self):
		plt.ylim(self.min_y,1)
		if(self.min_y > 1.0e-2):
			subsyy = [2,3,4,5,7]
		else:
			subsyy = []
		plt.yscale('log', subsy=subsyy, figure=self.fig)
		if(self.min_y > 1.0e-2):
			plt.tick_params(axis='y', which='minor', labelsize=14)
			plt.gca().yaxis.set_minor_formatter(FuncFormatter(format_label))
		else:
			plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=7))
		plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		plt.tick_params(axis='x', which='major', labelsize=20)
		Plotter.setup_fig(self)

class Parser:
	def __init__(self):
		self.parser = argparse.ArgumentParser()

	def parse_option(self):
		self.parser.add_argument('--logpath', dest='logpath',
				type=str, action='store',
				default='../logtmp/',
				help='path to load the log files')
		self.parser.add_argument('--figpath', dest='figpath',
				type=str, action='store',
				default='../figure/',
				help='path to save the figures')
		#positional arguments
		self.parser.add_argument('plottype', type=str,
				choices=[klass.PLOTTYPE for klass in Plotter.__subclasses__()],
				help='the plot type')

		return self.parser.parse_args()

def main():
	parser = Parser()
	args = parser.parse_option()

	stype = [runs[k] for k in runtype]

	#check losses are same
	losses = map(lambda st: st in L1, stype)
	assert losses.count(losses[0])==len(losses), "Need all stype have same loss"
	loss = "L1" if stype[0] in L1 else "L2"

	#check all or none in oneclass
	oneclass_types = map(lambda st: st in oneclass, stype)
	assert oneclass_types.count(oneclass_types[0])==len(oneclass_types), "Need all or none of stype are oneclass"
	if stype[0] in oneclass:
		real_clist = nlist
	else:
		real_clist = clist

	for plotter_class in Plotter.__subclasses__():
		if args.plottype == plotter_class.PLOTTYPE:
			plotter = plotter_class(stype, loss, dataset, real_clist, min(elist), args)
	plotter.draw_all()

if __name__ == '__main__':
	main()
