logfolder='log'
resumefolder='resume'

runtype = [
		"SVDD_L1_CY_1000",
		"SVDD_L1_FIRST_1000",
		"SVDD_L1_SEMIGD_CY_FIRST_1000", #r=4
		"SVDD_L1_SEMIGD_FIRST_CY_1000", #r=0.1
		#"SVDD_L1_SEMIGD_BATCH_FIRST_CY_1000", #r=0.1

		"ONECLASS_L1_CY_1000",
		"ONECLASS_L1_FIRST_1000",
		"ONECLASS_L1_SEMIGD_CY_FIRST_1000", #r=4
		"ONECLASS_L1_SEMIGD_FIRST_CY_1000", #r=0.1
		#"ONECLASS_L1_SEMIGD_BATCH_FIRST_CY_1000", #r=0.1
]
nlist = [
		0.1,
	#	0.01,
	#	0.005
		]
clist = []
elist = [1e-2]
greedy_random_rlist = [0.1]
random_greedy_rlist = [4]
other_semigd_rlist = []
S = 6e7 # maximum #O(n) operation
m = 300000 # maximum iter
timeout = 60 # maximum seconds
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

PROCESS_MAX = 1 # max number of jobs of a `gird' command

# plot config

MARKER = ["--",]
Y_BUFFER = 2
MIN_SQUASH = 0.2
YLIM = (1e-6, 1000)

uselabel = {
		50811: "cyclic-2cd",
		50812: "cyclic-2cd-shrink",
		50311: "greedy-2cd",
		50611: "cyclic-%scd-greedy",
		50211: "greedy-%s-cyclic",
		50212: "greedy-%s-cyclic-shrink",
		51111: "cyclic-0.1-greedy-%s-cyclic",
		60811: "cyclic-2cd",
		60812: "cyclic-2cd-shrink",
		60311: "greedy-2cd",
		60611: "cyclic-%scd-greedy",
		60211: "greedy-%s-cyclic",
		61311: "greedy-%s-cyclic-sort",
		60212: "greedy-%s-cyclic-shrink",
		61111: "cyclic-0.1-greedy-%s-cyclic",
}

