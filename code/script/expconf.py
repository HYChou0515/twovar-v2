logfolder='log_notime'

runtype = [
	"ONECLASS_L1_SEMIGD_1000",
	"ONECLASS_L1_SEMIGD_RAND_1000",
	#"ONECLASS_L1_RD_1000",
	"ONECLASS_L2_SEMIGD_1000",
	"ONECLASS_L2_SEMIGD_RAND_1000",
	#"ONECLASS_L2_RD_1000",
]
nlist = [0.01, 0.1, 0.2]
clist = [1]
elist = [1e-12]
rlist = [
	#1,0.95,0.9,0.5,0.1,0.01,0.001
	1,2,4,8,16,32,
	#0.001,0.002,0.004,0.008,0.016,0.032,
	0.064,0.128,0.256,
	0.512,0.9,0.95,
	float('inf'),
	#0.1, 0.2, 0.3, 0.4, 0.5,
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

MARKER = ["--","--","--","--","--","o-.","o-.","o-.","o-.","o-.",":"]
Y_BUFFER = 2
MIN_SQUASH = 0.2
YLIM = (1e-9, 1e20)

uselabel = {
		40111:"l1-random",
		40121:"l2-random",
		40211:"l1-semigd",
		40221:"l2-semigd",
		40311:"l1-semigdrd",
		40321:"l2-semigdrd",
		50111:"l1-random",
		50211:"l1-semigd",
		50311:"l1-first",
		50411:"l1-second",
		50511:"l1-semigdrd",
		50121:"l2-random",
		50221:"l2-semigd",
		50321:"l2-first",
		50421:"l2-second",
		50521:"l2-semigdrd",
}
