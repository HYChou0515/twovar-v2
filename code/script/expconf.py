logfolder='log_notime'
resumefolder='resume'

runtype = [
		#"ONECLASS_L1_CY_1000",
		#"ONECLASS_L1_SEMIGD_CY_FIRST_1000",
		"BIAS_L1_CY_1000",
		#"BIAS_L2_CY_1000",
		"BIAS_L1_SEMIGD_CY_FIRST_1000",
		#"BIAS_L2_SEMIGD_CY_FIRST_1000",
]
nlist = [
		0.1,
		0.01,
		0.2
		]
clist = [#1.0/32,
		#1,
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
	3,4,0.1
	#0.6,0.7,0.8,0.9,0.99
	#1,0.95,0.9,0.5,0.1,0.001
	#,0.5
	#0.2, 0.1
	#, 0.02, 0.01
]
m = 100000
timeout = 600
tolerance = 1e-10
dataset = [
	#"yahookr",
	"yahoojp",
	#"covtype.libsvm.binary.scale",
	#"a9a",
	#"ijcnn1",
	#"rcv1_train.binary",
	#"real-sim",
	#"news20.binary",
	]

PROCESS_MAX = 1

# plot config

MARKER = ["--","o-.","o-.","o-.","-.","-."]
Y_BUFFER = 2
MIN_SQUASH = 0.2
YLIM = (1e-9, 1e20)

uselabel = {
	50211: "semigd",
	50811: "cyclic",
	50611: "semigd-first",
	50911: "semigd-second",
	50212: "semigd-sh",
	50812: "cyclic-sh",
	50612: "semigd-first-sh",
	50912: "semigd-second-sh",

	40211: "semigd",
	40411: "semigd-first",
	40611: "cyclic",
	40711: "semigd-second",
	40212: "semigd-sh",
	40412: "semigd-first-sh",
	40612: "cyclic-sh",
	40712: "semigd-second-sh",

	40221: "semigd",
	40421: "semigd-first",
	40621: "cyclic",
	40721: "semigd-second",
	40222: "semigd-sh",
	40422: "semigd-first-sh",
	40622: "cyclic-sh",
	40722: "semigd-second-sh",
}

