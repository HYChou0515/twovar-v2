#!/usr/bin/env python
import math
import sys
import os
import numpy as np
from config import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
def format_label(x,pos):
#	return "10$^{%.5f}$" % np.log10(x)
	return "%1.e" %x
def format_labelx(x,pos):
	return "%.1f" %x
def draw_by_s(dstr, s, c, eps, LOSS):
	fig = plt.figure()
	cmap = plt.get_cmap('hsv')
	makr = ["--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,]
	totnum = len(s)
	for st in s:
		if st in semigd:
			totnum = totnum + len(rlist)

	colors = [cmap(j) for j in np.linspace(0, 1, totnum+1)]
	clridx = 0
	minirelval = 10
	for tp in s:
		if tp not in semigd:
			getrlist = [1]
		else:
			getrlist = rlist

		if tp not in shrink:
			e = 0.1
		else:
			e = eps
		for r in getrlist:
			relvals = []
			tms = []
			upd = []
			if tp not in semigd:
				logname = "%s_s%d_c%g_e%g"% (dstr, tp, c, e)
				if tp in oneclass:
					logname = "%s_s%d_c1_e%g_n%g"% (dstr, tp, e, c)
			else:
				logname = "%s_s%d_c%g_e%g_r%g"% (dstr, tp, c, e, r)
				if tp in oneclass:
					logname = "%s_s%d_c1_e%g_n%g_r%g"% (dstr, tp, e, c, r)
			CDsteps = 0
			if tp in biasobj:
				minimal = dobj["bias{}c{}".format(LOSS,c)][dstr]
			elif tp in oneclass:
				minimal = dobj["one{}c{}".format(LOSS,c)][dstr]
			else:
				minimal = dobj["{}c{}".format(LOSS,c)][dstr]
			print(logname)
			try:
				fl = open(LOGPATH+logname,"r")
			except IOError:
				print("error")
				continue
			if "s%d_c%g_iter"%(tp, c) in dlim.keys():
				if dstr in dlim["s%d_c%g_iter"%(tp, c)]:
					xlim = dlim["s%d_c%g_iter"%(tp, c)][dstr]
				else:
					xlim = 10000
			else:
				xlim = 10000
			for line in fl:
				line = line.split(' ')
				if 'iter' in line and 'obj' in line:
					#if CDsteps > xlim:
					#	break
					#if int(line[line.index('iter')+1]) > xlim:
					#	break
					val = float(line[line.index('obj')+1])
					relval = math.fabs((val - minimal)/minimal)
					if relval < 1e-8:
						break
					relvals.append(relval)
					CDsteps = CDsteps + (float(line[line.index('updsize')+1]))
					#CDsteps = CDsteps + (float(line[line.index('sucpair')+1]))
					tms.append(CDsteps)
					#tms.append(line[line.index('iter')+1])
					print('CDsteps: %.16g\t relval: %.16g' % (CDsteps, relval))
			#label setting
			lb = uselabel[tp]
			if tp in semigd:
				lb = "%s_%s" % (lb, r)
			minirelval = min(min(relvals),minirelval)
			if(option['gt']=='plot'):
				plt.plot(tms, relvals, makr[clridx], color=colors[clridx], figure=fig, label=lb, linewidth=3, markersize=6, markevery=0.1)
			if(option['gt']=='loglog'):
				plt.loglog(tms, relvals, makr[clridx], color=colors[clridx], figure=fig, label=lb, linewidth=3, markersize=6, markevery=0.1)
			clridx += 1


	plt.axhline(y=0.5, linestyle='--')
	plt.ylim(minirelval,1)
	if(minirelval > 1.0e-2):
		subsyy = [2,3,4,5,7]
	else:
		subsyy = []
	plt.yscale('log', subsy=subsyy, figure=fig)
	plt.tick_params(axis='y', which='major', labelsize=30)
	plt.tick_params(axis='x', which='major', labelsize=20)
	if(minirelval > 1.0e-2):
		plt.tick_params(axis='y', which='minor', labelsize=14)
		plt.gca().yaxis.set_minor_formatter(FuncFormatter(format_label))
	else:
		plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=7))
	if(option['gt']=='plot'):
		plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
	#plt.tight_layout()
	plt.legend(loc=0)
	tmps = ''
	for k in s:
		tmps += 's'+str(k)
	fig.savefig(FIGPATH+dstr+'_'+tmps+'_c'+str(c)+'_e'+str(e)+'iter_obj.png',format='png',dpi=200)

def draw_shrink(dstr, s, c, elist, LOSS):
	fig = plt.figure()
	cmap = plt.get_cmap('hsv')
	makr = ["--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,]
	totnum = len(s)
	for st in s:
		if st in semigd:
			totnum = totnum + len(rlist)
	colors = [cmap(j) for j in np.linspace(0, 1, totnum+1)]
	#colors = [(1,0,0,1),(1,0,0,1),(0,0,1,1),(0,0,1,1)]
	clridx = 0
	if c>=1:
		c = int(c)
	lasttime = 0
	minirelval = 10
	for tp in s:
		if tp in semigd:
			getrlist = rlist
		else:
			getrlist = [1]
		for r in getrlist:
			relvals = []
			tms = []
			lasttime = 0
			if tp not in shrink:
				eps = [0.1]
			else:
				eps = elist
			if tp in biasobj:
				minimal = dobj["bias{}c{}".format(LOSS,c)][dstr]
			elif tp in oneclass:
				minimal = dobj["one{}c{}".format(LOSS,c)][dstr]
			else:
				minimal = dobj["{}c{}".format(LOSS,c)][dstr]
			if "s%d_c%g_shrink"%(tp, c) in dlim.keys():
				if dstr in dlim["s%d_c%g_shrink"%(tp, c)]:
					xlim = dlim["s%d_c%g_shrink"%(tp, c)][dstr]
				else:
					xlim = 100
			else:
				xlim = 100

			for e in eps:
				if tp not in semigd:
					logname = "%s_s%d_c%g_e%g"% (dstr, tp, c, e)
					if tp in oneclass:
						logname = "%s_s%d_c1_e%g_n%g"% (dstr, tp, e, c)
				else:
					logname = "%s_s%d_c%g_e%g_r%g"% (dstr, tp, c, e, r)
					if tp in oneclass:
						logname = "%s_s%d_c1_e%g_n%g_r%g"% (dstr, tp, e, c, r)
				print(logname)
				try:
					fl = open(LOGPATH+logname,"r")
				except IOError:
					print("error")
					continue
				for line in fl:
					line = line.split(' ')
					if 't' in line and'obj' in line:
						val = float(line[line.index('obj')+1])
						relval = math.fabs((val - minimal)/minimal)
						if relval < 1e-4:
							break
						t =  float(line[line.index('t')+1])
						if t < lasttime:
							continue
						else:
							lasttime = t
						if t > xlim:
							break;
						tms.append(t)
						relvals.append(relval)
						print('t: %.16g\t relval: %.16g' % (t, relval))

			minirelval = min(min(relvals),minirelval)
			#label setting
			lb = uselabel[tp]
			if tp in semigd:
				lb = "%s_%s" % (lb, r)
			if(option['gt']=='plot'):
				plt.plot(tms, relvals, makr[clridx], color=colors[clridx], figure=fig, label=lb, linewidth=3, markersize=6, markevery=0.1)
			if(option['gt']=='loglog'):
				plt.loglog(tms, relvals, makr[clridx], color=colors[clridx], figure=fig, label=lb, linewidth=3, markersize=6, markevery=0.1)
			clridx += 1
	plt.ylim(minirelval,1)
	subsyy = []
	if(minirelval > 1.0e-2):
		subsyy = [2,3,4,5,7]
	else:
		subsyy = []
	plt.yscale('log', subsy=subsyy, figure=fig)
	plt.tick_params(axis='y', which='major', labelsize=30)
	plt.tick_params(axis='x', which='major', labelsize=20)
	if(minirelval > 1.0e-2):
		plt.tick_params(axis='y', which='minor', labelsize=14)
		plt.gca().yaxis.set_minor_formatter(FuncFormatter(format_label))
	else:
		plt.tick_params(axis='y', which='major', labelsize=30)
		plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=6))
	if(option['gt']=='plot'):
		plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
	#plt.title(dstr)
	#plt.tight_layout()
	plt.legend(loc=0)
	tmps = ''
	for k in s:
		tmps += 's'+str(k)
	fig.savefig(FIGPATH+dstr+'_'+tmps+'_c'+str(c)+'_e'+str(e)+'shrink_time.png',format='png',dpi=200)

def draw_shrink_new(dstr, s, c, eps, LOSS):
	fig = plt.figure()
	cmap = plt.get_cmap('hsv')
	makr = ["--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,"--", "ko-.", "-", "ko-." ,]
	totnum = len(s)
	for st in s:
		if st in semigd:
			totnum = totnum + len(rlist)
	colors = [cmap(j) for j in np.linspace(0, 1, totnum+1)]
	#colors = [(1,0,0,1),(1,0,0,1),(0,0,1,1),(0,0,1,1)]
	clridx = 0
	if c>=1:
		c = int(c)
	lasttime = 0
	minirelval = 10
	for tp in s:
		if tp in semigd:
			getrlist = rlist
		else:
			getrlist = [1]
		for r in getrlist:
			relvals = []
			tms = []
			lasttime = 0
			if tp not in shrink:
				e = 0.01
			if tp in biasobj:
				minimal = dobj["bias{}c{}".format(LOSS,c)][dstr]
			elif tp in oneclass:
				minimal = dobj["one{}c{}".format(LOSS,c)][dstr]
			else:
				minimal = dobj["{}c{}".format(LOSS,c)][dstr]
			if "s%d_c%g_shrink"%(tp, c) in dlim.keys():
				if dstr in dlim["s%d_c%g_shrink"%(tp, c)]:
					xlim = dlim["s%d_c%g_shrink"%(tp, c)][dstr]
				else:
					xlim = 5000
			else:
				xlim = 5000
			for e in eps:
				if tp not in semigd:
					logname = "%s_s%d_c%g_e%g"% (dstr, tp, c, e)
					if tp in oneclass:
						logname = "%s_s%d_c1_e%g_n%g"% (dstr, tp, e, c)
				else:
					logname = "%s_s%d_c%g_e%g_r%g"% (dstr, tp, c, e, r)
					if tp in oneclass:
						logname = "%s_s%d_c1_e%g_n%g_r%g"% (dstr, tp, e, c, r)
				try:
					fl = open(LOGPATH+logname,"r")
				except IOError:
					print("error")
					continue
				for line in fl:
					line = line.split(' ')
					if 't' in line and'obj' in line:
						val = float(line[line.index('obj')+1])
						relval = math.fabs((val - minimal)/minimal)
						if relval < 1e-10:
							break
						t =  float(line[line.index('t')+1])
						if t < lasttime:
							continue
						else:
							lasttime = t
						if t > xlim:
							break;
						tms.append(t)
						relvals.append(relval)

			minirelval = min(min(relvals),minirelval)
			#label setting
			lb = uselabel[tp]
			if tp in semigd:
				lb = "%s_%s" % (lb, r)
			if(option['gt']=='plot'):
				plt.plot(tms, relvals, makr[clridx], color=colors[clridx], figure=fig, label=lb, linewidth=3, markersize=6, markevery=0.1)
			if(option['gt']=='loglog'):
				plt.loglog(tms, relvals, makr[clridx], color=colors[clridx], figure=fig, label=lb, linewidth=3, markersize=6, markevery=0.1)
			clridx += 1
	plt.ylim(minirelval,1000)
	subsyy = []
	if(minirelval > 1.0e-2):
		subsyy = [2,3,4,5,7]
	else:
		subsyy = []
	plt.yscale('log', subsy=subsyy, figure=fig)
	plt.tick_params(axis='y', which='major', labelsize=30)
	plt.tick_params(axis='x', which='major', labelsize=20)
	if(minirelval > 1.0e-2):
		plt.tick_params(axis='y', which='minor', labelsize=14)
		plt.gca().yaxis.set_minor_formatter(FuncFormatter(format_label))
	else:
		plt.tick_params(axis='y', which='major', labelsize=30)
		plt.gca().yaxis.set_major_locator(plt.LogLocator(numticks=6))
	if(option['gt']=='plot'):
		plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
	#plt.title(dstr)
	#plt.tight_layout()
	plt.legend(loc=0)
	tmps = ''
	for k in s:
		tmps += 's'+str(k)
	fig.savefig(FIGPATH+dstr+'_'+tmps+'_c'+str(c)+'_e'+str(e)+'shrink_time.png',format='png',dpi=200)
plt.close('all')

LOGPATH = '../../logtmp/'
FIGPATH = '../../figure/'


plottype = sys.argv[1]

option = {'gt': 'plot', 'log': 'off', 't': '1'}
if '-gt' in sys.argv:
	option['gt']=sys.argv[sys.argv.index('-gt')+1]
if '-log' in sys.argv:
	option['log']= sys.argv[sys.argv.index('-log')+1]
if '-t' in sys.argv:
	option['-t'] = sys.argv[sys.argv.index('-t')+1]

stype = [runs[k] for k in runtype]
matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)
if stype[0] in L1:
	loss = "L1"
else:
	loss = "L2"
for st in stype:
	if st in oneclass:
		clist = nlist
if "cd" == plottype:
	for dstr in dataset:
		for c in clist:
			for e in elist:
				draw_by_s(dstr, stype, c, e, loss)
elif "time" == plottype:
	for dstr in dataset:
		for c in clist:
			if option['t'] == '1':
				draw_shrink(dstr, stype, c, elist, loss)
			if option['t'] == '2':
				draw_shrink_new(dstr, stype, c, elist, loss)
