logfolder='log_notime'

runtype = [
		"ONE_L1_SEMIGD_DUALOBJ_1000",
		"ONE_L1_CY_1000",
]
nlist = [1]
clist = [1.0/32,1,32]
elist = [1e-12]
rlist = [
#	0.01,0.03,
	0.1#, 0.2, 0.3, 0.4, 0.5,
]
m = 10000
timeout = 3600
tolerance = 1e-10
dataset = [
	#"heart_scale",
	"a9a.min",
	"ijcnn1.min",
	"rcv1_train.binary.min",
	"real-sim_36k.min",
	"news20.binary_2k.min",
	"yahoojp_2k.min",
	"covtype.libsvm.binary.scale.min",
	"yahookr_1k.min",
	"a9a.bal",
	"ijcnn1.bal",
	"rcv1_train.binary.bal",
	"real-sim_36k.bal",
	"news20.binary_2k.bal",
	"yahoojp_2k.bal",
	"covtype.libsvm.binary.scale.bal",
	"yahookr_1k.bal",
	"a9a",
	"ijcnn1",
	"rcv1_train.binary",
	"real-sim_36k",
	"news20.binary_2k",
	"yahoojp_2k",
	"covtype.libsvm.binary.scale",
	"yahookr_1k",
	]

# plot config

MARKER = ["--","--","--","--","--","o-.","o-.","o-.","o-.","o-.",":"]
Y_BUFFER = 2
MIN_SQUASH = 0.2
YLIM = (1e-9, 1e20)

uselabel = {
		20111:"random",
		20411:"semigd",
}
