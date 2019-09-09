#!/usr/bin/env python
import math, bisect
import sys, os, resource
import numpy as np
from plotconst import *
from liblrconf import *
from expconf import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import argparse
import itertools
class StubException(Exception):
	pass
def format_label(x,pos):
#	return "10$^{%.5f}$" % np.log10(x)
	return "%1.e" %x
def format_labelx(x,pos):
	return "%.1f" %x

def expo_slice_avg(l,base):
	#slice l into average of base^0, base^1, base^2
	#e.g. l=[1,2,3,4,5,6,7], base=2
	#ret=[1, (2+3)/2, (4+5+6+7)/4]
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

def exponent_of(x):
	# return the exponent of x
	# input +-0.3, then return +-0.1
	# input 0.01, then return 0.01
	if x==0:
		return 0
	if x > 0:
		flag=1
	else:
		flag=-1
	return flag * 10**int(("%e"%x).split('e')[-1])

class Plotter(object):
	def __init__(self, stype, loss, dataset, clist, eps, args):
		self.dataset = dataset
		self.clist = clist
		self.stype = stype
		self.loss = loss
		self.eps = eps
		self.logpath = args.logpath+'/'
		self.figpath = args.figpath+'/'
		self.filesuffix = args.suffix
		self.relobj_denom = args.relobj_denom

		matplotlib.rc('xtick', labelsize=20)
		matplotlib.rc('ytick', labelsize=20)

		self.YLIM = YLIM # the range of Y (this value is set for relval)
		self.X_LABEL=None
		self.Y_LABEL=None

	def calculate_relval(self, measure_val, true_val):
		diff = math.fabs(measure_val-true_val)
		true_val_abs = math.fabs(true_val)
		#if true_val_abs < 1e-16:
		#	return diff
		#	#return math.fabs(2*math.atan(true_val/measure_val)-math.pi/2)
		#if diff < 1e-16:
		#	return 1e-16
		return diff/(self.relobj_denom+true_val_abs)

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
		totnum = 0
		for tp in self.stype:
			getrlist = self.get_rlist(tp)
			e = self.get_eps(tp, self.eps)
			for r in getrlist:
				try:
					self.get_logfile(tp, e, r)
					totnum += 1
				except IOError:
					continue
				except StubException:
					totnum += 1
					continue
		self.colors = [cmap(j) for j in np.linspace(0, 1, totnum+1)]
		self.makr = MARKER
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
				except StubException:
					clridx += 1
					continue
				lb = uselabel[tp]
				if is_semigd(tp):
					lb = lb % r
				plt.plot(xs, ys, self.makr[clridx],
						color=self.colors[clridx],
						figure=self.fig, label=lb,
						linewidth=3, markersize=6, markevery=0.1)
				clridx += 1

		min_xs, max_xs, min_ys, max_ys = zip(*self.xy_range) # list of tuples to tuple of lists
		plt.hlines(1e-1,0,max_xs,linestyles='-')
		plt.hlines(1e-2,0,max_xs,linestyles='--')
		plt.hlines(1e-3,0,max_xs,linestyles='-.')
		plt.subplots_adjust(bottom=0.1,left=0.1,right=0.98,top=0.98)
		self.setup_fig()
		plt.close(self.fig)

	def get_rlist(self, tp):
		if is_semigd(tp):
			return rlist
		else:
			return [1]

	def get_eps(self, tp, orig_e):
		if is_shrink(tp):
			return orig_e
		else:
			return 0.1

	def get_logfile(self, tp, e, r):
		if tp is None:
			print('Using a token')
			raise StubException
		if is_semigd(tp):
			if is_oneclass(tp) == 1:
				logname = "%s_s%d_c1_e%g_n%g_r%g"% (self.dstr, tp, e, self.c, r)
			else:
				logname = "%s_s%d_c%g_e%g_r%g"% (self.dstr, tp, self.c, e, r)
		else:
			if is_oneclass(tp) == 1:
				logname = "%s_s%d_c1_e%g_n%g"% (self.dstr, tp, e, self.c)
			else:
				logname = "%s_s%d_c%g_e%g"% (self.dstr, tp, self.c, e)
		try:
			f = open(self.logpath+logname,"r")
			print("open file: "+self.logpath+logname)
			return f
		except IOError:
			print("cannot find file: "+self.logpath+logname)
			raise IOError

	def get_obj(self, line):
		return float(line[line.index('obj')+1])

	def get_minimal(self, tp):
		if is_biasobj(tp):
			dobj_key = "bias%sc%g" % (self.loss,self.c)
		elif is_oneclass(tp) == 1: # ocsvm
			dobj_key = "one%sc%g" % (self.loss,self.c)
		elif is_oneclass(tp) == 2: # svdd
			dobj_key = "svdd%sc%g" % (self.loss,self.c)
		else:
			dobj_key = "%sc%g" % (self.loss,self.c)
		try:
			return dobj[dobj_key][self.dstr]
		except KeyError:
			print("\"" + self.dstr + "\" is not in " + dobj_key)
			raise KeyError

	def set_xy_lim(self, min_xs, max_xs, min_ys, max_ys):
		# set xlim
		print("minmax:%d maxmin:%d" % (min(max_xs), max(min_xs)))
		if max(min_xs) > min(max_xs):
			print(max(min_xs)+min(max_xs))
			plt.xlim(0, min(max(max_xs), max(min_xs)+min(max_xs)))
		else:
			print(min(max_xs) / MIN_SQUASH, max(max_xs))
			plt.xlim(0, min(min(max_xs) / MIN_SQUASH, max(max_xs)))
		# set ylim
		y_exponents = [1e1,1e3,1e7,1e12, float('Inf')]
		y_max = y_exponents[bisect.bisect_right(y_exponents, exponent_of(max(max_ys)))]
		plt.ylim(self.YLIM[0]/2, y_max)

	def setup_fig(self):
		plt.legend(loc=0)
		self.fig.savefig(self.figpath+self.get_figname_fmt()%(self.dstr, self.c)+'_'+self.filesuffix+'.png',format='png',dpi=100)

class OpPerSucs_ObjPlotter(Plotter):
	PLOTTYPE = "opspersucs_obj"
	XLIM = (9e-3, 1000)
	def init_new_fig(self):
		self.xy_range = []
		self.summary = []
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
		for s in filter(lambda s: s is not None, self.stype):
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_opspersucsobj"
	def get_xy(self, tp, e, r):
		logfile = self.get_logfile(tp, e, r)
		ys = []
		xs = []
		minimal = self.get_minimal(tp)

		min_x = 1e10
		max_x = 0
		max_y = 0
		min_y = 1e10
		iternum = 0
		nrnop = 0
		sucsize = 0
		for line in logfile:
			line = line.split(' ')
			if 'ops_per_sucs' in line and 'obj' in line:
				val = self.get_obj(line)
				relval = self.calculate_relval(val, minimal)
				ops_per_sucs = 2*float(line[line.index('ops_per_sucs')+1])
				if ops_per_sucs > self.YLIM[1] or relval < OpPerSucs_ObjPlotter.XLIM[0]:
					continue
				if ops_per_sucs < self.YLIM[0] or relval > OpPerSucs_ObjPlotter.XLIM[1]:
					continue
				iternum = int(line[line.index('iter')+1])
				sucsize += int(line[line.index('sucsize')+1])
				nrnop += ops_per_sucs * int(line[line.index('sucsize')+1])
				ys.append(ops_per_sucs)
				xs.append(relval)
				min_x = min(relval, min_x)
				max_x = max(relval, max_x)
				min_y = min(ops_per_sucs,min_y)
				max_y = max(ops_per_sucs,max_y)
		xy_range = (min_x, max_x, min_y, max_y)
		if all(xy_range):
			self.xy_range.append(xy_range)
		print("xy_range={}".format(xy_range))
		self.summary.append("iter: %d op/suc: %g" % (iternum, float(nrnop)/sucsize))
		logfile.close()
		return ys, xs
	def setup_fig(self):
		min_xs, max_xs, min_ys, max_ys = zip(*self.xy_range) # list of tuples to tuple of lists
		Plotter.set_xy_lim(self, min_xs, max_xs, min_ys, max_ys)
		plt.ylim(min(min_ys), max(max_ys))
		plt.yscale('linear', figure=self.fig)
		plt.xscale('log', figure=self.fig)
		plt.gca().xaxis.set_major_locator(plt.LogLocator(numticks=7))
		plt.xlim(min(min_xs), max(max_xs))

		summary = "\n".join(self.summary)
		plt.annotate(summary, (0,0), (0, -20), xycoords='axes fraction', textcoords='offset points', va='top')
		plt.subplots_adjust(bottom=0.15)
		Plotter.setup_fig(self)

class NrNOpPlotter(Plotter):
	PLOTTYPE = "nrnop"
	XLIM = (0, 1e12)
	def init_new_fig(self):
		self.xy_range = []
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
		for s in filter(lambda s: s is not None, self.stype):
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_nrnop"
	def get_xy(self, tp, e, r):
		logfile = self.get_logfile(tp, e, r)
		ys = []
		xs = []
		nr_n_ops = 0
		minimal = self.get_minimal(tp)

		min_x = None
		max_x = None
		max_y = 0
		min_y = 1e10
		stop = False
		for line in logfile:
			line = line.split(' ')
			if stop:
				break
			if 'iter' in line and 'obj' in line:
				val = self.get_obj(line)
				relval = self.calculate_relval(val, minimal)
				nr_n_ops = (float(line[line.index('nr_n_ops')+1]))
				if relval > self.YLIM[1] or nr_n_ops < NrNOpPlotter.XLIM[0]:
					continue
				if relval < self.YLIM[0] or nr_n_ops > NrNOpPlotter.XLIM[1]:
					#if len(ys) == 0 or len(xs) == 0:
					#	break
					#interpolate_rate = 1-max((nr_n_ops-NrNOpPlotter.XLIM[1])/(nr_n_ops-xs[-1]),
					#		(math.log(relval)-math.log(self.YLIM[0]))/(math.log(relval)-math.log(ys[-1])))
					#relval = math.exp(math.log(ys[-1]) + (math.log(relval)-math.log(ys[-1]))*interpolate_rate)
					#nr_n_ops = xs[-1] + (nr_n_ops-xs[-1])*interpolate_rate
					stop = True
				ys.append(relval)
				xs.append(nr_n_ops)
				if min_x is None:
					min_x = nr_n_ops
				max_x = nr_n_ops
				min_y = min(relval,min_y)
				max_y = max(relval,max_y)
				#print('nr_n_ops: %.16g\t relval: %.16g' % (nr_n_ops, relval))
		xy_range = (min_x, max_x, min_y, max_y)
		if all(xy_range):
			self.xy_range.append(xy_range)
		print("xy_range={}".format(xy_range))
		logfile.close()
		return ys, xs
	def setup_fig(self):
		min_xs, max_xs, min_ys, max_ys = zip(*self.xy_range) # list of tuples to tuple of lists
		Plotter.set_xy_lim(self, min_xs, max_xs, min_ys, max_ys)
		if(min(min_ys) > 1.0e-2):
			subsyy = [2,3,4,5,7]
		else:
			subsyy = []
		plt.yscale('log', subsy=subsyy, figure=self.fig)
		plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=7))
		plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		plt.tick_params(axis='x', which='major', labelsize=20)

		##add xy label and adjust the margin accordingly
		#plt.xlabel('number of $O(n)$ operations', fontsize=20)
		#plt.ylabel('$(f-f^*)/f*$', fontsize=20)
		#plt.subplots_adjust(bottom=0.13,left=0.133)

		Plotter.setup_fig(self)

class CdPlotter(Plotter):
	PLOTTYPE = "cd"
	XLIM = (0, 1e12)
	def init_new_fig(self):
		self.xy_range = []
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
		for s in filter(lambda s: s is not None, self.stype):
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_cdsteps"
	def get_xy(self, tp, e, r):
		logfile = self.get_logfile(tp, e, r)
		ys = []
		xs = []
		CDsteps = 0
		minimal = self.get_minimal(tp)

		min_x = None
		max_x = None
		max_y = 0
		min_y = 1e10
		stop = False
		for line in logfile:
			line = line.split(' ')
			if stop:
				break
			if 'iter' in line and 'obj' in line:
				val = self.get_obj(line)
				relval = self.calculate_relval(val, minimal)
				CDsteps = (float(line[line.index('iter')+1]))
				if relval > self.YLIM[1] or CDsteps < CdPlotter.XLIM[0]:
					continue
				if relval < self.YLIM[0] or CDsteps > CdPlotter.XLIM[1]:
					if len(ys) == 0 or len(xs) == 0:
						break
					interpolate_rate = 1-max((CDsteps-CdPlotter.XLIM[1])/(CDsteps-xs[-1]),
							(math.log(relval)-math.log(self.YLIM[0]))/(math.log(relval)-math.log(ys[-1])))
					relval = math.exp(math.log(ys[-1]) + (math.log(relval)-math.log(ys[-1]))*interpolate_rate)
					CDsteps = xs[-1] + (CDsteps-xs[-1])*interpolate_rate
					stop = True
				ys.append(relval)
				xs.append(CDsteps)
				if min_x is None:
					min_x = CDsteps
				max_x = CDsteps
				min_y = min(relval,min_y)
				max_y = max(relval,max_y)
				#print('CDsteps: %.16g\t relval: %.16g' % (CDsteps, relval))
		xy_range = (min_x, max_x, min_y, max_y)
		if all(xy_r is not None for xy_r in xy_range):
			print("add xy_range={}".format(xy_range))
			self.xy_range.append(xy_range)
		else:
			print("discard xy_range={}".format(xy_range))
		logfile.close()
		return ys, xs
	def setup_fig(self):
		print("xy_range={}".format(self.xy_range))
		min_xs, max_xs, min_ys, max_ys = zip(*self.xy_range) # list of tuples to tuple of lists
		Plotter.set_xy_lim(self, min_xs, max_xs, min_ys, max_ys)
		if(min(min_ys) > 1.0e-2):
			subsyy = [2,3,4,5,7]
		else:
			subsyy = []
		plt.yscale('log', subsy=subsyy, figure=self.fig)
		plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=7))
		plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		plt.tick_params(axis='x', which='major', labelsize=20)
		Plotter.setup_fig(self)

class TimePlotter(Plotter):
	PLOTTYPE = "time"
	XLIM = (0, 1e9)
	def init_new_fig(self):
		self.xy_range = []
	def get_xlim(self, tp):
		key = "s%d_c%g_shrink" % (tp, self.c)
		if key in dlim.keys() and self.dstr in dlim[key]:
			return dlim[key][self.dstr]
		else:
			return 100

	def get_figname_fmt(self):
		sname = ''
		for s in filter(lambda s: s is not None, self.stype):
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_time"
	def get_xy(self, tp, e, r):
		logfile = self.get_logfile(tp, e, r)
		ys = []
		xs = []
		minimal = self.get_minimal(tp)

		min_x = None
		max_x = None
		max_y = 0
		min_y = 1e10
		stop = False
		for line in logfile:
			line = line.split(' ')
			if stop:
				break
			if 't' in line and'obj' in line:
				val = self.get_obj(line)
				relval = self.calculate_relval(val, minimal)
				t =  float(line[line.index('t')+1])
				if relval > self.YLIM[1] or t < TimePlotter.XLIM[0]:
					continue
				if relval < self.YLIM[0] or t > TimePlotter.XLIM[1]:
					if len(ys) == 0 or len(xs) == 0:
						break
					interpolate_rate = 1-max((t-TimePlotter.XLIM[1])/(t-xs[-1]),
							(math.log(relval)-math.log(self.YLIM[0]))/(math.log(relval)-math.log(ys[-1])))
					relval = math.exp(math.log(ys[-1]) + (math.log(relval)-math.log(ys[-1]))*interpolate_rate)
					t = xs[-1] + (t-xs[-1])*interpolate_rate
					stop = True
				xs.append(t)
				ys.append(relval)
				if min_x is None:
					min_x = t
				max_x = t
				min_y = min(relval,min_y)
				max_y = max(relval,max_y)
				#print('t: %.16g\t relval: %.16g' % (t, relval))
		xy_range = (min_x, max_x, min_y, max_y)
		if all(xy_range):
			self.xy_range.append(xy_range)
		print("xy_range={}".format(xy_range))
		logfile.close()
		return ys, xs
	def setup_fig(self):
		min_xs, max_xs, min_ys, max_ys = zip(*self.xy_range) # list of tuples to tuple of lists
		Plotter.set_xy_lim(self, min_xs, max_xs, min_ys, max_ys)
		if(min(min_ys) > 1.0e-2):
			subsyy = [2,3,4,5,7]
		else:
			subsyy = []
		plt.yscale('log', subsy=subsyy, figure=self.fig)
		plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=7))
		plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		plt.tick_params(axis='x', which='major', labelsize=20)
		Plotter.setup_fig(self)

class SucrCdPlotter(Plotter):
	PLOTTYPE = "sucrate-cd"
	def init_new_fig(self):
		self.max_y = 0
	def get_figname_fmt(self):
		sname = ''
		for s in filter(lambda s: s is not None, self.stype):
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_sucrate"
	def get_xy(self, tp, e, r):
		logfile = self.get_logfile(tp, e, r)
		ys = []
		xs = []
		cdsteps = 0

		for line in logfile:
			line = line.split(' ')
			if 'iter' in line and 'obj' in line:
				val = float(line[line.index('sucpair')+1])/float(line[line.index('updsize')+1])
				cdsteps = (float(line[line.index('cdsteps')+1]))
				ys.append(val)
				self.max_y = max(val,self.max_y)
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
	def init_new_fig(self):
		self.min_y = 10
	def get_figname_fmt(self):
		sname = ''
		for s in filter(lambda s: s is not None, self.stype):
			sname += 's'+str(s)
		return "%s_"+sname+"_c%g_sucupd"
	def get_xy(self, tp, e, r):
		logfile = self.get_logfile(tp, e, r)
		ys = []
		xs = []
		sucpair = 0
		minimal = self.get_minimal(tp)

		for line in logfile:
			line = line.split(' ')
			if 'iter' in line and 'obj' in line:
				val = self.get_obj(line)
				relval = self.calculate_relval(val, minimal)
				sucpair = float(line[line.index('ttl_sucsize')+1])
				if relval > self.YLIM[1]:
					continue
				if relval < self.YLIM[0]:
					break;
				ys.append(relval)
				self.min_y = min(relval,self.min_y)
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
		self.parser.add_argument('--suffix', dest='suffix',
				type=str, action='store',
				default='',
				help='output figure suffix (before format suffix)')
		self.parser.add_argument('--relobj', dest='relobj_denom',
				type=float, action='store',
				default=0.0,
				help='relobj=(f-f*)/([this]+f*)')
		#positional arguments
		self.parser.add_argument('plottype', type=str,
				choices=[klass.PLOTTYPE for klass in Plotter.__subclasses__()],
				help='the plot type')

		return self.parser.parse_args()

def main():
	parser = Parser()
	args = parser.parse_option()

	stype = [runs[k] for k in runtype]

	#grouping l1/l2 and svc/ocsvm
	svc_l1loss_stype = []
	svc_l2loss_stype = []
	ocsvm_l1loss_stype = []
	ocsvm_l2loss_stype = []
	svdd_l1loss_stype = []
	svdd_l2loss_stype = []
	for st in stype:
		if st is None:
			# st is STUB
			ocsvm_l1loss_stype.append(st)
			svdd_l1loss_stype.append(st)
			svc_l1loss_stype.append(st)
			ocsvm_l2loss_stype.append(st)
			svdd_l2loss_stype.append(st)
			svc_l2loss_stype.append(st)
		else:
			if is_L1(st):
				if is_oneclass(st)==1:
					ocsvm_l1loss_stype.append(st)
				elif is_oneclass(st)==2:
					svdd_l1loss_stype.append(st)
				else:
					svc_l1loss_stype.append(st)
			else:
				if is_oneclass(st)==1:
					ocsvm_l2loss_stype.append(st)
				elif is_oneclass(st)==2:
					svdd_l2loss_stype.append(st)
				else:
					svc_l2loss_stype.append(st)

	#L1 ocsvm
	if len(filter(lambda s: s is not None, ocsvm_l1loss_stype)) != 0:
		for plotter_class in Plotter.__subclasses__():
			if args.plottype == plotter_class.PLOTTYPE:
				plotter = plotter_class(ocsvm_l1loss_stype, "L1", dataset, nlist, min(elist), args)
		plotter.draw_all()

	#L1 svdd
	if len(filter(lambda s: s is not None, svdd_l1loss_stype)) != 0:
		for plotter_class in Plotter.__subclasses__():
			if args.plottype == plotter_class.PLOTTYPE:
				plotter = plotter_class(svdd_l1loss_stype, "L1", dataset, clist, min(elist), args)
		plotter.draw_all()

	#L1 svc
	if len(filter(lambda s: s is not None, svc_l1loss_stype)) != 0:
		for plotter_class in Plotter.__subclasses__():
			if args.plottype == plotter_class.PLOTTYPE:
				plotter = plotter_class(svc_l1loss_stype, "L1", dataset, clist, min(elist), args)
		plotter.draw_all()

	#L2 ocsvm
	if len(filter(lambda s: s is not None, ocsvm_l2loss_stype)) != 0:
		for plotter_class in Plotter.__subclasses__():
			if args.plottype == plotter_class.PLOTTYPE:
				plotter = plotter_class(ocsvm_l2loss_stype, "L2", dataset, nlist, min(elist), args)
		plotter.draw_all()

	#L2 svdd
	if len(filter(lambda s: s is not None, svdd_l2loss_stype)) != 0:
		for plotter_class in Plotter.__subclasses__():
			if args.plottype == plotter_class.PLOTTYPE:
				plotter = plotter_class(svdd_l2loss_stype, "L2", dataset, clist, min(elist), args)
		plotter.draw_all()

	#L2 svc
	if len(filter(lambda s: s is not None, svc_l2loss_stype)) != 0:
		for plotter_class in Plotter.__subclasses__():
			if args.plottype == plotter_class.PLOTTYPE:
				plotter = plotter_class(svc_l2loss_stype, "L2", dataset, clist, min(elist), args)
		plotter.draw_all()

def memory_limit():
	soft, hard = resource.getrlimit(resource.RLIMIT_AS)
	resource.setrlimit(resource.RLIMIT_AS, (get_memory() * 1024 / 2, hard))

def get_memory():
	with open('/proc/meminfo', 'r') as mem:
		free_memory = 0
		for i in mem:
			sline = i.split()
			if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
				free_memory += int(sline[1])
	return free_memory

if __name__ == '__main__':
	assert sys.version_info[:3] == (2,7,12)
	assert matplotlib.__version__ == '1.5.3'
	memory_limit()
	try:
		main()
	except MemoryError as e:
		sys.stderr.write('Memory exceed limit\n')
		sys.stderr.write(str(e))

