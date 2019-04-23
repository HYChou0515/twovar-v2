logfolder='log_notime'
resumefolder='resume'

runtype = [
		"ONE_L1_CY_1000",
		"ONE_L2_CY_1000",
		"BIAS_L1_SEMIGD_CY_FIRST_1000",
		"BIAS_L2_SEMIGD_CY_FIRST_1000",
		"ONE_L1_CY_SH",
		"ONE_L2_CY_SH",
		"BIAS_L1_SEMIGD_CY_FIRST_SH",
		"BIAS_L2_SEMIGD_CY_FIRST_SH",
]
nlist = [
		0.1,
		0.01,
		0.2
		]
clist = [1.0/32,
		1,
		32
		]
elist = [1e-12]
rlist = [
	#1,0.95,0.9,0.5,0.1,0.01,0.001
	#1,2,4,8,16,32,
	#0.001,0.002,0.004,0.008,0.016,0.032,
	#0.064,0.128,0.256,
	#0.512,0.9,0.95,
	#float('inf'),
	2,3,4
	#0.6,0.7,0.8,0.9,0.99
	#1,0.95,0.9,0.5,0.1,0.001
	#,0.5
	#0.2, 0.1
	#, 0.02, 0.01
]
m = 100000
timeout = 1200
tolerance = 1e-10
dataset = [
	"yahookr",
	"yahoojp",
	"covtype.libsvm.binary.scale",
	"a9a",
	"ijcnn1",
	"rcv1_train.binary",
	"real-sim",
	"news20.binary",
	]

PROCESS_MAX = 1

# plot config

MARKER = ["--","-.","-.","-.","-.","-."]
Y_BUFFER = 2
MIN_SQUASH = 0.2
YLIM = (1e-9, 1e20)

uselabel = {
	40411: "bias-first-l1",
	20111: "no-bias-l1",
	40412: "bias-first-l1-sh",
	20112: "no-bias-l1-sh",
	40421: "bias-first-l2",
	20121: "no-bias-l2",
	40422: "bias-first-l2-sh",
	20122: "no-bias-l2-sh",
}

