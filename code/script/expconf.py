logfolder='log_time'
resumefolder='resume'

runtype = [
		## FOR EXP: 2-cd vs new-4 vs semigd vs semigd-rd
		"SVDD_L1_CY_1000",
		"SVDD_L1_FIRST_1000",
		"SVDD_L1_SEMIGD_CY_FIRST_1000", #r=4
		"SVDD_L1_SEMIGD_1000", #r=0.1
		"SVDD_L1_SEMIGD_BATCH_1000", #r=0.1

		"ONECLASS_L1_CY_1000",
		"ONECLASS_L1_FIRST_1000",
		"ONECLASS_L1_SEMIGD_CY_FIRST_1000", #r=4
		"ONECLASS_L1_SEMIGD_1000", #r=0.1
		"ONECLASS_L1_SEMIGD_BATCH_1000", #r=0.1
]
nlist = [0.1,0.01]
clist = []
elist = [1e-2]
rlist = [
	4,
	0.1
]
m = 100
S = 1e9
timeout = 600
tolerance = 1e-10
dataset = [
#	"yahookr",
#	"yahoojp",
#	"covtype.libsvm.binary.scale",
#	"a9a",
#	"ijcnn1",
#	"rcv1_train.binary",
#	"real-sim",
#	"news20.binary",
	"avazu-app"
	]

PROCESS_MAX = 1

# plot config

MARKER = ["--",]
Y_BUFFER = 2
MIN_SQUASH = 0.2
YLIM = (1e-6, 1000)

uselabel = {
		50811: "cyclic-2cd",
		50812: "cyclic-2cd-shrink",
		50611: "cyclic-%scd-greedy",
		50211: "greedy-%s-cyclic",
		50212: "greedy-%s-cyclic-shrink",
		60811: "cyclic-2cd",
		60812: "cyclic-2cd-shrink",
		60611: "cyclic-%scd-greedy",
		60211: "greedy-%s-cyclic",
		61311: "greedy-%s-cyclic-sort",
		60212: "greedy-%s-cyclic-shrink",
}

