logfolder='log_notime'
resumefolder='resume'

runtype = [
		## FOR EXP: 2-cd vs new-4 vs semigd vs semigd-rd
		"SVDD_L1_CY_1000",
		"SVDD_L1_FIRST_1000",
		"SVDD_L1_SEMIGD_CY_FIRST_1000", #r=4
		"SVDD_L1_SEMIGD_1000", #r=0.1

		"ONECLASS_L1_CY_1000",
		"ONECLASS_L1_FIRST_1000",
		"ONECLASS_L1_SEMIGD_CY_FIRST_1000", #r=4
		"ONECLASS_L1_SEMIGD_1000", #r=0.1
]
nlist = [
		#0.2,
		0.1,
		#0.05,
		0.01,0.005,
		#0.001,0.0001
		]
clist = []
elist = [1e-2]
rlist = [
	4,
	0.1
]
m = 300000
timeout = 21600
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
YLIM = (1e-6, 1000)

uselabel = {
		50811: "cyclic-2cd",
		50311: "greedy-2cd",
		50611: "cyclic-%scd-greedy",
		50211: "greedy-%s-cyclic",
		60811: "cyclic-2cd",
		60311: "greedy-2cd",
		60611: "cyclic-%scd-greedy",
		60211: "greedy-%s-cyclic",
}

