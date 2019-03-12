#!/usr/bin/env python
from sklearn.datasets import *
from sklearn.preprocessing import normalize
import numpy as np
import matplotlib.pyplot as plt
from random import randint
import sys, os, math
import scipy

def a():
	d = load_svmlight_file('heart_scale')
	norm = normalize(d[0])
	dump_svmlight_file(norm, d[1], 'heart_scale_norm', False)

def b():
	data = np.loadtxt(open('a'))
	plt.scatter(data[:,0], data[:,1])
	plt.show()

def c():
	d = load_svmlight_file('data/heart_scale_norm')

	l = len(d[1])
	x=[]
	y=[]
	for s in range(20):
		i = randint(0,l-1)
		j = randint(0,l-1)
		xi = d[0][i,:]
		xj = d[0][j,:]
		x.append(xi.dot(xj.T)[0][0][0,0])
		yi = d[1][i]
		yj = d[1][j]
		y.append(yi*yj)

	plt.scatter(x,y)
	plt.show()

def d(dname):
#	d = load_svmlight_file('data/a9a')
	d = load_svmlight_file(dname)
	x = d[0]
	x = x.todense()
	y = d[1]
	_,n = x.shape
	l = len(y)
	d = np.zeros(n)
	for i in range(n):
		ci = (sum(x[:,i])/l)[0,0]
		print(ci)
		d[i] = ci
#		print(x[:,i].toarray())
		x[:,i] -= ci
	print(np.linalg.norm(d))
#	dump_svmlight_file(x, y, 'a9a_cent', False)

def e():
	d = load_svmlight_file('data/a9a')
#	d = load_svmlight_file('a9a_cent')
	x = d[0]
	x = x.todense()
	y = d[1]
	_,n = x.shape
	l = len(y)
	d = np.zeros(n)
	for i in range(n):
		ci = (np.max(x[:,i])+np.min(x[:,i]))/2
		print(ci)
		d[i] = ci
#		print(x[:,i].toarray())
		x[:,i] -= ci
	print(np.linalg.norm(d))
	dump_svmlight_file(x, y, 'a9a_cent2', False)

def mean_confidence_interval(data, confidence=0.95):
	a = 1.0 * np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
	return m, m-h, m+h

def _see_expectation(x, t2):
	l,_ = x.shape
	S = []
	PS = []
	T1 = 10
	if t2 < 1:
		T2 = int(math.floor(t2*l))
	else:
		T2 = int(t2)
	for t in range(T1):
		s = 0.0
		ps = 0.0
		for t in range(T2):
			i = randint(0,l-1)
			j = randint(0,l-1)
			dd = (x[i].dot(x[j].T))[0,0]
			s += dd
			if dd > 0:
				ps += 1.0
		S.append(s/T2)
		PS.append(ps/T2)
	a,b,c = mean_confidence_interval(S)
	d,e,f = mean_confidence_interval(PS)
	return a,b,c,d,e,f

def see_expectation(dname, t2):
	try:
		d = load_svmlight_file(dname)
	except IOError:
		sys.stderr.write("Error: cannot open data: \"%s\"\n" % dname)
		return -1
	a,b,c,d,e,f = _see_expectation(d[0], t2)
	print("%s %g (%g %g) %g (%g %g)" % (dname, a,b,c, d,e, f))

def see_many_expectation():
	for dname in sys.argv[1:]:
		see_expectation(dname)

def make_min_expect(dname):
	try:
		d = load_svmlight_file(dname)
	except IOError:
		sys.stderr.write("Error: cannot open data: \"%s\"\n" % dname)
		return -1
	x = d[0]
	l, n = x.shape
	b = x.T.dot(np.ones((l,1)))/(-l)
	dnew = os.path.join('data_minexp', os.path.basename(dname)+'.min')
	dump_svmlight_file(x+b.T, d[1], dnew, False)

def make_zero_expect(dname):
	try:
		d = load_svmlight_file(dname)
	except IOError:
		sys.stderr.write("Error: cannot open data: \"%s\"\n" % dname)
		return -1
	x = d[0]
	l, n = x.shape
	b = x.T.dot(np.ones((l,1)))/(-l)
	x = x.todense()
        print(x.nbytes)
	_,_,_,a,_,_ = _see_expectation(x+b.T, 0.1)
	if a > 0.5:
		print('besta>0.5')
		dnew = os.path.join('data_minexp', os.path.basename(dname))
		dump_svmlight_file(x+b.T, d[1], dnew, False)
		return
	print('find besta')
        print('R=1')
        print('a='+str(a))
	b = b.T
	inita = a
	besta = a
	bestR = 1.0
	t = 0
        x += b
	while abs(besta-0.5) > 1e-2 and t < 20:
		t+=1
		r = 1.0
		R = 1.0
		a = inita
		for i in range(8):
			r /= 2
			if a > 0.5:
				R += r
                                x += r*b
			else:
				R -= r
                                x -= r*b
			_,_,_,a,_,_ = _see_expectation(x, 0.1)
			if abs(a-0.5) < abs(besta-0.5):
				besta = a
				bestR = R
                        print('R='+str(R))
                        print('a='+str(a))
                x += (1-R)*b
        x -= (1-R)*b
	print(besta)
	dnew = os.path.join('data_minexp', os.path.basename(dname))
	dump_svmlight_file(x, d[1], dnew, False)

if sys.argv[1] == '0exp':
	for dname in sys.argv[2:]:
		print(dname)
		make_zero_expect(dname)
if sys.argv[1] == 'exp':
	for dname in sys.argv[2:]:
		make_min_expect(dname)
if sys.argv[1] == 'see':
	for dname in sys.argv[3:]:
		see_expectation(dname, float(sys.argv[2]))
if sys.argv[1] == 'd':
	for dname in sys.argv[2:]:
		d(dname)

