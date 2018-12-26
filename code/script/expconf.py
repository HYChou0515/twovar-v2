from liblrconf import *

uselabel = labeltest
#runtype = ["TWO_L1_RD_SH","TWO_L1_RD_SH2"]
runtype = [
	#"ONE_L1_CY_1000",
	"ONE_L1_SEMIGD_DUALOBJ_1000",
	"TWO_L1_CY_1000",
	"BIAS_L1_SEMIGD_1000",
]
nlist = [0.1]
clist = [1]
elist = [0.1]
rlist = [
	#1,0.99,0.98,0.97,0.96,0.95,0.94,0.93,0.92,0.91,0.9,0.8,0.7,0.6,0.5,0.1,0.001
	#1,0.95,0.9,0.5,0.1,0.01,0.001
	#1,2,4,8,16,32
	10,
	#1,0.95,0.9,0.5,0.1,0.001
	#,0.5
	#0.2, 0.1
	#, 0.02, 0.01
]
m = 5000
dataset = [
	"heart_scale",
#	"a9a",
#	"ijcnn1",
#	"rcv1_train.binary",
#	"real-sim",
#	"news20.binary",
#	"yahoojp",
#	"covtype.libsvm.binary.scale",
#	"yahookr",
	]


