logfolder='log_notime'

runtype = [
	#"ONE_L1_CY_1000",
	#"ONE_L1_SEMIGD_DUALOBJ_1000",
	"ONE_L1_SEMIGD_DUALOBJ_YBAL_1000",
	"ONE_L2_SEMIGD_DUALOBJ_YBAL_1000",
]
nlist = [0.1]
clist = [1.0/32, 1, 32]
elist = [1e-12]
rlist = [
	#1,0.95,0.9,0.5,0.1,0.01,0.001
	#1,2,4,8,16,32,
	#0.001,0.002,0.004,0.008,0.016,0.032,
	#0.064,0.128,0.256,
	#0.512,0.9,0.95,
	#float('inf'),
	0.1, 0.2, 0.3, 0.4, 0.5,
	#0.6,0.7,0.8,0.9,0.99
	#1,0.95,0.9,0.5,0.1,0.001
	#,0.5
	#0.2, 0.1
	#, 0.02, 0.01
]
m = 100000
timeout = 7200
tolerance = 1e-10
dataset = [
	#"heart_scale",
	#"a9a",
	#"ijcnn1",
	#"rcv1_train.binary",
	#"real-sim",
	#"news20.binary",
	#"yahoojp",
	#"covtype.libsvm.binary.scale",
	#"yahookr",
	"rcv1_train.binary.bal",
	"rcv1_train.binary.min",
	"real-sim_36k",
	"real-sim_36k.min",
	"real-sim_36k.bal",
	"news20.binary_2k.bal",
	"news20.binary_2k",
	"news20.binary_2k.min",
	"yahookr_1k",
	"yahookr_1k.bal",
	"yahookr_1k.min",
	"yahoojp_2k",
	"yahoojp_2k.min",
	"yahoojp_2k.bal",
	"covtype.libsvm.binary.scale.bal",
	"covtype.libsvm.binary.scale.min",
	"a9a.min",
	"a9a.bal",
	"ijcnn1.bal",
	"ijcnn1.min",
	]

PROCESS_MAX = 96

# plot config

MARKER = ["--","--","--","--","--","o-.","o-.","o-.","o-.","o-.",":"]
Y_BUFFER = 2
MIN_SQUASH = 0.2
YLIM = (1e-9, 1e20)

uselabel = {
	20111: 'random',
	20411: 'semigd',
	20711: 'semigd-ybal'
}

