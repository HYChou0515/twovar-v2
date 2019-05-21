logfolder='log_notime'
resumefolder='resume'

runtype = [
		# FOR EXP: 2-cd vs new-4 vs nobias
		#"ONE_L1_CY_1000", #do not need to run as change of cdstep do not affect this
		"BIAS_L1_CY_1000",
		"BIAS_L1_SEMIGD_CY_FIRST_1000", #r=4
		#"STUB",
		#"STUB",
		##"ONE_L1_CY_SH", #do not need to run as change of cdstep do not affect this
		##"BIAS_L1_CY_SH",
		##"BIAS_L1_SEMIGD_CY_FIRST_SH", #r=4

		#"ONE_L2_CY_1000", #do not need to run as change of cdstep do not affect this
		"BIAS_L2_CY_1000",
		"BIAS_L2_SEMIGD_CY_FIRST_1000", #r=4
		#"STUB",
		#"STUB",
		##"ONE_L2_CY_SH", #do not need to run as change of cdstep do not affect this
		##"BIAS_L2_CY_SH",
		##"BIAS_L2_SEMIGD_CY_FIRST_SH", #r=4

		# FOR EXP: 2-cd vs new-4 vs semigd vs semigd-rd
		#"STUB",
		"ONECLASS_L1_CY_1000",
		"ONECLASS_L1_SEMIGD_CY_FIRST_1000", #r=4
		"ONECLASS_L1_SEMIGD_1000", #r=0.1
		"ONECLASS_L1_SEMIGD_RAND_1000", #r=0.1
		##"ONECLASS_L1_CY_SH",
		##"ONECLASS_L1_SEMIGD_CY_FIRST_SH", #r=4
		##"ONECLASS_L1_SEMIGD_SH", #r=0.1
		##"ONECLASS_L1_SEMIGD_RAND_SH", #r=0.1
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
	4, #0.1
]
m = 100000
timeout = 36000
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

MARKER = ["--",]
Y_BUFFER = 2
MIN_SQUASH = 0.2
YLIM = (1e-9, 1e20)

uselabel = {
		20111: "no-bias-cd",
		40611: "bias-cd",
		40411: "bias-perm-greedy",
		20121: "no-bias-cd",
		40621: "bias-cd",
		40421: "bias-perm-greedy",
		50811: "perm-cd",
		50611: "perm-greedy",
		50211: "semi-greedy",
		50511: "semi-greedy-random",
}

