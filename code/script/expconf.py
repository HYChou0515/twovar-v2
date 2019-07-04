logfolder='log_notime'
resumefolder='resume'

runtype = [
		## FOR EXP: 2-cd vs new-4 vs semigd vs semigd-rd
		"STUB",
		"SVDD_L1_CY_1000",
		"SVDD_L1_SEMIGD_CY_FIRST_1000", #r=4
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
		20112: "no-bias-cd-sh",
		40612: "bias-cd-sh",
		40412: "bias-perm-greedy-sh",
		20122: "no-bias-cd-sh",
		40622: "bias-cd-sh",
		40422: "bias-perm-greedy-sh",
		50812: "perm-cd-sh",
		50612: "perm-greedy-sh",
		50212: "semi-greedy-sh",
		50512: "semi-greedy-random-sh",
		60811: "perm-cd",
		60611: "perm-greedy",
		60211: "semi-greedy",
		60511: "semi-greedy-random",
}

