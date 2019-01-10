from liblrconf import *

logfolder='log_notime'
process_num=96

runtype = [
	"ONE_L1_SEMIGD_1000",
	"ONE_L1_SEMIGD_RAND_1000",
	"ONE_L1_SEMIGD_DUALOBJ_1000",
	"ONE_L1_SEMIGD_DUALOBJ_RAND_1000",
	"ONE_L1_CY_1000",
	"ONE_L1_RD_1000",
	"ONE_L1_CY_SH",
	"ONE_L1_RD_SH",
]
nlist = [0.1]
clist = [1.0/32, 32]
elist = [1e-12]
rlist = [
	#1,0.99,0.98,0.97,0.96,0.95,0.94,0.93,0.92,0.91,0.9,0.8,0.7,0.6,0.5,0.1,0.001
	#1,0.95,0.9,0.5,0.1,0.01,0.001
	#1,2,4,8,16,32
	0.1, 0.3,
	0.5,
	#0.9,0.99
	#1,0.95,0.9,0.5,0.1,0.001
	#,0.5
	#0.2, 0.1
	#, 0.02, 0.01
]
m = 100000
timeout = 7200
dataset = [
	"heart_scale",
	"a9a",
	"ijcnn1",
	"rcv1_train.binary",
	"real-sim",
	"news20.binary",
	"yahoojp",
	"covtype.libsvm.binary.scale",
	"yahookr",
	]


uselabel = {
		20411:"semigd-dualobj",
		20311:"semigd-pg",
		20511:"semigd-pg-rand",
		20611:"semigd-dualobj-rand",
		20211:"random",
		20212:"random-sh",
		20111:"cyclic",
		20112:"cyclic-sh",
}

# plot config

MARKER = ["--","--","o-.","o-.","-.","-."]
MIN_SQUASH = 0.2
YLIM = (1e-9, 1)
