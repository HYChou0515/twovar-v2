logfolder='log_notime'
resumefolder='resume'

runtype = [
		## FOR EXP: 2-cd vs new-4 vs semigd vs semigd-rd
		"STUB",
		"ONECLASS_L1_CY_1000",
		"ONECLASS_L1_SEMIGD_CY_FIRST_1000", #r=4
		"ONECLASS_L1_SEMIGD_CY_DUALOBJ_1000", #r=4
		"ONECLASS_L1_SEMIGD_1000", #r=0.1
		"ONECLASS_L1_SEMIGD_RAND_1000", #r=0.1
		## FOR EXP: 2-cd vs new-4 vs semigd vs semigd-rd
		"STUB",
		"SVDD_L1_CY_1000",
		"SVDD_L1_SEMIGD_CY_FIRST_1000", #r=4
		"SVDD_L1_SEMIGD_CY_DUALOBJ_1000", #r=4
		"SVDD_L1_SEMIGD_1000", #r=0.1
		"SVDD_L1_SEMIGD_RAND_1000", #r=0.1
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
elist = [1e-2]
rlist = [
	4,
	0.1
]
m = 100000
timeout = 50400
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
		50811: "perm-cd",
		50611: "perm-greedy",
		50911: "perm-greedy-2nd",
		50211: "semi-greedy",
		50511: "semi-greedy-random",
		60811: "perm-cd",
		60611: "perm-greedy",
		60911: "perm-greedy-2nd",
		60211: "semi-greedy",
		60511: "semi-greedy-random",
}

