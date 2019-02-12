logfolder='log_notime'

runtype = [
	"ONECLASS_L1_RD_1000",
	"ONECLASS_L1_SEMIGD_1000",
	"ONECLASS_L1_FIRST_1000",
	"ONECLASS_L1_SECOND_1000",
	"ONECLASS_L2_RD_1000",
	"ONECLASS_L2_SEMIGD_1000",
	"ONECLASS_L2_FIRST_1000",
	"ONECLASS_L2_SECOND_1000",
	#"BIAS_L1_RD_1000",
	#"BIAS_L1_SEMIGD_1000",
	#"BIAS_L2_RD_1000",
	#"BIAS_L2_SEMIGD_1000",
]
nlist = [0.1]
clist = [1.0/32, 1, 32]
elist = [1e-12]
rlist = [
	#1,0.99,0.98,0.97,0.96,0.95,0.94,0.93,0.92,0.91,0.9,0.8,0.7,0.6,0.5,0.1,0.001
	#1,0.95,0.9,0.5,0.1,0.01,0.001
	#1,2,4,8,16,32
	0.1, 0.3,
	0.5,
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
	"heart_scale",
	"a9a",
	"ijcnn1",
	"rcv1_train.binary",
	"real-sim",
	"news20.binary",
	"yahoojp",
	"covtype.libsvm.binary.scale",
	"yahookr",
	]

# plot config

MARKER = ["--","--","o-.","o-.","-.","-."]
Y_BUFFER = 2
YLIM = (1e-9, 1e3)

uselabel = {
		50111:"l1-random",
		50211:"l1-semigd",
		50311:"l1-first",
		50411:"l1-second",
		50121:"l2-random",
		50221:"l2-semigd",
		50321:"l2-first",
		50421:"l2-second",
}
