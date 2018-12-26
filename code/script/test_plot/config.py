runs = {"ONE_L1_CY_1000iter":21, "ONE_L1_RD_1000iter":22, "ONE_L1_RD_SH":23,\
	"ONE_L2_CY_1000iter":24, "ONE_L2_RD_1000iter":25, "ONE_L2_RD_SH":26,\
        "ONE_L1_SEMIGD_1000iter":27,\
        "ONE_L1_SEMIGD_1000iter_SHUFFLE":271,\
        "ONE_L1_SEMIGD_SH":28,\
        "ONE_L2_SEMIGD_1000iter":29,\
        "ONE_L2_SEMIGD_SH":30,\
	"TWO_L1_CY_1000iter":31, "TWO_L1_RD_1000iter":32, "TWO_L1_RD_SH":33,\
	"TWO_L2_CY_1000iter":34, "TWO_L2_RD_1000iter":35, "TWO_L2_RD_SH":36,\
	"ONE_L2_CY_SH":1, "ONE_L1_CY_SH":3,\
	"TWO_L1_SEMICY_1000iter": 37,\
	"TWO_L2_SEMICY_1000iter": 38,\
	"TWO_L1_SEMIRDONE_1000iter": 371,\
	"TWO_L2_SEMIRDONE_1000iter": 381,\
	"TWO_L1_SEMIRDTWO_1000iter": 372,\
	"TWO_L2_SEMIRDTWO_1000iter": 382,\
	"TWO_L1_RD_SH2":39,\
	"TWO_L2_RD_SH2":40,\
	"BIAS_L1_RD_1000iter":41, "BIAS_L2_RD_1000iter":42,\
	"BIAS_L1_semigd_1000iter":43, "BIAS_L2_semigd_1000iter":44,\
	"BIAS_L1_RD_SH":45, "BIAS_L2_RD_SH":46,\
	"BIAS_L1_SEMIGD_SH": 47, "BIAS_L2_SEMIGD_SH": 48, \
	"ONECLASS_L1_RD_1000iter":51,\
	"ONECLASS_L1_SEMIGD_1000iter":52,\
	"ONECLASS_L1_FIRST_1000iter":53,\
	"ONECLASS_L1_SECOND_1000iter":54,\
	"ONECLASS_L1_RD_SH":55,\
	"ONECLASS_L1_SEMIGD_SH":56,\
	"ONECLASS_L2_RD_1000iter":57,\
	"ONECLASS_L2_SEMIGD_1000iter":58,\
	"ONECLASS_L2_FIRST_1000iter":59,\
	"ONECLASS_L2_SECOND_1000iter":60,\
	"ONECLASS_L2_RD_SH":61,\
	"ONECLASS_L2_SEMIGD_SH":62,\

}
alltype = runs.values()
biasobj = [41,42,43,44,45,46,47,48]
semigd = [27,271,28,29,30,43,44,47,48,52,56,58,62]
shrink = [28,29,23,26,33,36,1,3,39,40,45,46,47,48,55,56,61,62]
oneclass = [51,52,53,54,55,56,57,58,59,60,61]
L1 = [21,22,23,27,271,28,31,32,33,41,43,45,47,3,37,43,51,52,53,54,55,56,371,372]
label1 = {31:"perm", 32:"random",37:"semi-random",\
	  34:"perm", 35:"random",38:"semi-random",\
}
label2 = {21:"1-CD-perm", 22:"1-CD-random", 23:"1-CD-random",\
	  24:"1-CD-perm", 25:"1-CD-random", 26:"1-CD-random",\
	  27:"1-CD-smgd", 271:"1-CD-smgd-shf",\
	  31:"2-CD-perm", 32:"2-CD-random", 33:"2-CD-random",\
	  34:"2-CD-perm", 35:"2-CD-random", 36:"2-CD-random",\
	  3:"1-CD-perm",  1:"1-CD-perm",
          371:"2-CD-cyclic", 372:"2-CD-cyclic",\
          381:"2-CD-cyclic", 382:"2-CD-cyclic",\
}
labeltest = {21:"1-CD-perm", 22:"1-CD-random", 23:"1-CD-random",\
	  24:"1-CD-perm", 25:"1-CD-random", 26:"1-CD-random",\
	  27:"1-CD-smgd", 271:"1-CD-smgd-shf",\
	  31:"2-CD-perm", 32:"2-CD-random", 33:"2-CD-random",\
	  34:"2-CD-perm", 35:"2-CD-random", 36:"2-CD-random",\
	  3:"1-CD-perm",  1:"1-CD-perm",
          371:"2-CD-cyclic", 372:"2-CD-cyclic",\
          381:"2-CD-cyclic", 382:"2-CD-cyclic",\
	  51:"oneclass-random" ,52:"oneclass-semigd",
	  43: "semigd", 47: "semigd_shrink", 41: "random", 45: "random_shrink",
	  44: "semigd", 48: "semigd_shrink", 42: "random", 46: "random_shrink",
	  58:"oneclass-random-L2", 55:"oneclass-random-L1"
}

label3 = {41:"2-CD-bias" ,32:"2-CD-nobias", 42:"2-CD-bias" ,35:"2-CD-nobias"}
label4 = {3:"1-CD-shrink", 21:"1-CD", 32:"2-CD", 33:"2-CD-shrink",\
	  1:"1-CD-shrink", 24:"1-CD", 35:"2-CD", 36:"2-CD-shrink"}
label5 = {41:"2-CD-random", 42:"2-CD-random", 43:"2-CD-semigd", 44:"2-CD-semigd",
	51:"oneclass-random" ,52:"oneclass-semigd"
}
label6 = {33: "2-nobias-shrink", 39:"2-nobias-shrink2"}
label7 = {51: "1class-L1-random", 55: "1class-L1-random-shrink", 52: "1class-L1-semigd", 56: "1class-L1-semigd-shrink"}
label8 = {43: "semigd", 47: "semigd_shrink", 41: "random", 45: "random_shrink", 44: "semigd", 48: "semigd_shrink", 42: "random", 46: "random_shrink"}
label9 = {53:"L1-1st-order", 54:"L1-2nd-order", 51:"L1-random", 52:"L1:semigd"}
label10 = {1:"nobias_random", 3:"nobias_random", 45:"bias_random", 47:"bias_semi", 46:"bias_random", 48:"bias_semi"}
label11 = {41:"2-CD-random", 45:"2-CD-random-sh", 43:"2-CD-semigd", 47:"2-CD-semigd-sh",
	51:"oneclass-random" ,52:"oneclass-semigd"}
uselabel = labeltest
#runtype = ["TWO_L1_RD_SH","TWO_L1_RD_SH2"]
#runtype = ["ONECLASS_L1_RD_1000iter", "ONECLASS_L1_SEMIGD_1000iter", "ONECLASS_L1_RD_SH", "ONECLASS_L1_SEMIGD_SH"]
#runtype = ["BIAS_L1_RD_1000iter", "BIAS_L1_semigd_1000iter", "BIAS_L1_RD_SH", "BIAS_L1_SEMIGD_SH"]
#runtype = ["ONECLASS_L1_FIRST_1000iter", "ONECLASS_L1_SECOND_1000iter", "ONECLASS_L1_RD_1000iter", "ONECLASS_L1_SEMIGD_1000iter"]
#runtype = ["BIAS_L2_RD_1000iter", "BIAS_L2_semigd_1000iter", "BIAS_L2_RD_SH", "BIAS_L2_SEMIGD_SH"]

#runtype = ["ONE_L2_CY_SH", "ONE_L1_CY_SH"]#, "ONE_L1_CY_1000iter", "ONE_L2_CY_1000iter" ]

#runtype = ["BIAS_L2_RD_1000iter", "BIAS_L2_semigd_1000iter", "BIAS_L2_RD_SH", "BIAS_L2_SEMIGD_SH"]
#runtype = ["BIAS_L2_semigd_1000iter",  "BIAS_L2_SEMIGD_SH"]
#runtype = ["BIAS_L1_RD_1000iter", "BIAS_L1_semigd_1000iter", "BIAS_L1_RD_SH", "BIAS_L1_SEMIGD_SH", ]
#runtype = ["BIAS_L1_semigd_1000iter",  "BIAS_L1_SEMIGD_SH"]

#runtype = ["ONE_L2_CY_SH", "ONE_L1_CY_SH", "ONE_L1_CY_1000iter", "ONE_L2_CY_1000iter" ]
#runtype=["BIAS_L1_RD_1000iter", "BIAS_L2_RD_1000iter",\
#	"BIAS_L1_semigd_1000iter", "BIAS_L2_semigd_1000iter",\
#	"BIAS_L1_RD_SH", "BIAS_L2_RD_SH",\
#	"BIAS_L1_SEMIGD_SH", "BIAS_L2_SEMIGD_SH", \
#	]
#elist = [0.1, 0.01, 0.001, 0.0001]
#elist = [0.0001]
#clist = [1, 0.03125, 32, 128, 1024]
#clist = [128,1024]
#rlist = [0.2, 0.1, 0.02, 0.01, 0.002, 0.001]
#rlist = [0.1, 0.01, 0.001]
#nlist = [0.1]
#m=10000


#runtype = [
#	"TWO_L1_SEMIRDONE_1000iter",\
#	"TWO_L2_SEMIRDONE_1000iter",\
#	"TWO_L1_SEMIRDTWO_1000iter",\
#	"TWO_L2_SEMIRDTWO_1000iter",\
#"TWO_L1_CY_1000iter", "TWO_L1_RD_1000iter",\
#"TWO_L2_CY_1000iter", "TWO_L2_RD_1000iter",
#]

runtype = [
	"ONE_L1_SEMIGD_1000iter",
	"ONE_L1_CY_1000iter",
	"ONE_L1_RD_1000iter",
	"TWO_L1_RD_SH",
	"BIAS_L1_RD_1000iter",
	"BIAS_L1_SEMIGD_SH",
	"BIAS_L2_semigd_1000iter",
	#"ONECLASS_L1_RD_SH",
	#"ONECLASS_L2_SEMIGD_1000iter",
]
clist = [1]
elist = [0.01]
rlist = [#1
	#,0.5
	0.2, 0.1
	#, 0.02, 0.01
]
nlist = [0.1]
m = 1000
dataset = ["a9a",
	#"ijcnn1",
	#"news20.binary",
	#"rcv1_train.binary",
	#"real-sim",
	#"yahoojp",
#	"covtype.libsvm.binary.scale",
#	"yahookr",
	]
#dataset = ["heart_scale"]
#dataset = ["a1a",
#	"breast-cancer_scale",
#	"heart_scale",
#	"australian_scale"
#	]
#dataset = ["australian_scale"]
#dataset = ["a1a"]
#runtype=[
#   "ONECLASS_L1_RD_1000iter"\
# , "ONECLASS_L1_SEMIGD_1000iter"\
# , "ONECLASS_L1_FIRST_1000iter"\
# , "ONECLASS_L1_SECOND_1000iter"\
# , "ONECLASS_L1_RD_SH"\
# , "ONECLASS_L1_SEMIGD_SH"\
# , "ONECLASS_L2_RD_1000iter"\
# , "ONECLASS_L2_SEMIGD_1000iter"\
# , "ONECLASS_L2_FIRST_1000iter"\
# , "ONECLASS_L2_SECOND_1000iter"\
# , "ONECLASS_L2_RD_SH"\
# , "ONECLASS_L2_SEMIGD_SH"
#]

#elist = [ 0.01]
#clist = [1]
#rlist = [0.2, 0.1, 0.02, 0.01 , 0.002, 0.001]
#rlist = [0.1]
#nlist = [0.01, 0.1, 0.2]
#nlist = [0.1]

#runtype=[
#"TWO_L1_RD_SH2",\
#"TWO_L2_RD_SH2",\
#]


#runtype = ["BIAS_L1_RD_SH", "BIAS_L1_SEMIGD_SH", "ONE_L1_CY_SH"]
#runtype = ["BIAS_L2_RD_SH", "ONE_L2_CY_SH","BIAS_L2_SEMIGD_SH" ]
#elist = [ 0.1, 0.01, 0.001, 0.0001]
#clist = [1]
#rlist = [0.1]
#dataset = ["yahookr", "covtype.libsvm.binary"]
#runtype = [ "TWO_L1_CY_1000iter",  "TWO_L1_RD_1000iter","TWO_L1_SEMICY_1000iter"]
#elist = [ 0.01]
#clist = [1,8192]

#runtype = [ "TWO_L2_CY_1000iter",  "TWO_L2_RD_1000iter","TWO_L2_SEMICY_1000iter"]
#elist = [ 0.01]
#clist = [1,8192]

#dataset = ["heart_scale"]
#dataset = ["yahookr"]
#dataset = ["a9a","ijcnn1","news20.binary",\
#	   "rcv1_train.binary","real-sim","yahoojp"]


# need to be checked later
dobj = {
    "oneL2c0.1":{#THIS IS JUST A STUB, DONT USE IT AS EXP
	"a9a": 94360681.2818033695120,
	"ijcnn1": 6412124.210081208359921,
	"ijcnn1": 6412124.2100812522,
	"news20.binary": 229824.45181480217973,
#	"rcv1_train.binary": 55051.248030351707612,
	"rcv1_train.binary": 55051.248030351656,
#	"real-sim": 282831.525057278166050,
	"real-sim": 282831.52505727799,
	"yahoojp": 2867829.835220990236848,
    },
    "oneL1c0.2":{
	"a9a": 94360681.2818033695120,
	"ijcnn1": 6412124.210081208359921,
	"ijcnn1": 6412124.2100812522,
	"news20.binary": 229824.45181480217973,
#	"rcv1_train.binary": 55051.248030351707612,
	"rcv1_train.binary": 55051.248030351656,
#	"real-sim": 282831.525057278166050,
	"real-sim": 282831.52505727799,
	"yahoojp": 2867829.835220990236848,
    },
    "oneL1c0.1":{
	"a9a": 21622408.022816251963376,
#	"ijcnn1": 1569167.168492316131029,
	"ijcnn1": 1569167.1684923134,
#	"news20.binary": 41283.881972819392103,
	"news20.binary": 41283.88197281755,
	"rcv1_train.binary": 11331.28987903422420,
	"real-sim": 47167.016173025476745,
	"yahoojp": 197132.995628784934521,
    },
    "oneL1c0.01":{
	"a9a": 180123.032132618769533,
#	"ijcnn1": 15027.1238209391450,
	"ijcnn1": 15027.123809602508,
#	"news20.binary": 101.481072729227019,
	"news20.binary": 101.48107272922547,
#	"rcv1_train.binary": 97.713305289904080,
	"rcv1_train.binary": 97.713305188030290,
	"real-sim": 205.475956164127350,
	"yahoojp": 442.799337673292710,
    },
    "L1c8192":{
	"a9a": -93532332.6676734685897827,
	"ijcnn1": -69948865.0377422869205475,
	"news20.binary": -610011.4759008947294205,
	"rcv1_train.binary": -224951.4860682494472712,
	"real-sim": -2715246.1517096990719438,
	"yahoojp": -10861242.2472088448703289,
	"heart_scale": -777407.015664262231439,
	"covtype.libsvm.binary": -131408.78622,#could be small
	"yahookr": -39648244.84681542,
	"covtype.libsvm.binary.scale": -2765352655.81279516#100,000
	#-2767598655.16365861, #2561070
    },
    "L1c32":{
	"news20.binary": -5191.225594965,
	"a9a": -365535.0613,
	"ijcnn1": -273298.871245,
	"rcv1_train.binary": -4415.13148310,
	"real-sim": -26466,
	"yahoojp": -183009.81295,
    },
    "L1c1": {
	"rcv1_test": -39031.082174,
	"covtype.libsvm.binary.scale": -337953.029319151713,
	"yahookr": -144380.9835506559466012,
	"australian_scale": -200.500000,
	"splice_scale": -556.698396,
	"webspam": -12096.542471,
	"webspam_uni": -69591.024894,
	"webspam_wc_normalized_unigram": -69591.024894,
	"HIGGS": -8909820.625936,
	"SUSY": -2556563.336213,
	"ijcnn1.t": -15492.876661,
	"gisette_scale": -0.668683,
	"news20": -2562.542777,
	"epsilon_normalized": -109694.535560,
	"kdda": -821893.005048,
	"kddb": -1884911.576190,
	"url": -10736.121324,
	"url_combined": -10736.121324,
	"heart_scale": -96.4982779946964513,
	"a9a": -11433.8076970376823738,
	"ijcnn1": -8596.0038722944773326,
	"news20.binary": -2562.5427771712793401,
	"rcv1_train.binary": -1781.7303640360603367,
	"real-sim": -5402.0548794539135997,
	"yahoojp": -25932.1394146107923007,
	"covtype.libsvm.binary": -7555.482649466,
    },
    "L1c0.25":{
	"a9a": -2864.8765262514107235,
	"ijcnn1": -2186.4128780127202845,
	"news20.binary": -1525.2053920552516502,
	"rcv1_train.binary": -924.1368154180046304,
	"real-sim": -2482.5367215910819141,
    },
    "L1c0.5":{
	"a9a": -5721.6271924909697191,
	"ijcnn1": -4324.5934820310685609,
	"news20.binary": -2075.8308801364341889,
	"rcv1_train.binary": -1309.9995758980692244,
	"real-sim": -3696.1932663871161822,
    },
    "L1c2":{
	"a9a": -22856.9250076922726294,
	"ijcnn1": -17135.8241437759097607,
	"news20.binary": -2862.7117633669572569,
	"rcv1_train.binary": -2278.8731527108348018,
	"real-sim": -7690.9944430220884897,
    },
    "L1c4":{
	"a9a": -45702.2729259744082810,
	"ijcnn1": -34213.9642026758519933,
	"news20.binary": -3056.3966041606213366,
	"rcv1_train.binary": -2720.2674859238159115,
	"real-sim": -10597.8279588544974104,
    },
    "L2c1": {
	"rcv1_test": -44655.657965,
	"covtype.libsvm.binary.scale": -403423.63037590693,
	"yahookr": -140767.00986357091,
	"australian_scale": -276.689688,
	"splice_scale": None,
	"webspam": -10813.749434,
	"webspam_uni": -89938.633777,
	"webspam_wc_normalized_unigram": -89938.633777,
	"HIGGS": -9880685.418335,
	"SUSY": -3114688.507448,
	"ijcnn1.t": -21402.976954,
	"gisette_scale": -0.667967,
	"news20": -1947.384705, ## should be -1947.384704
	"epsilon_normalized": -131303.806194,
	"kdda": -652663.443215,
	"kddb": -1460362.265430,
	"url": -9139.968342,
	"url_combined": -9139.968342,
	"heart_scale": -121.134724,
	"a1a": -637.905580210829612,
	"a9a": -13742.397304382038783,
	"heart_scale": -121.134724436870442,
	"ijcnn1": -11023.190049223218011,
	"news20.binary": -1947.384704485529937,
	"phishing": -2037.970466784219525,
	"yahoojp": -25956.664097291635699,
	"news20.binary": -1947.384704488515126,
	"rcv1_train.binary": -1412.360148443366143,
	"real-sim":-4396.382163307219344,
	"covtype.libsvm.binary": -7574.7988735,
    },
    "L2c0.1":{
	"a1a": -66.886812971233,
	"a9a": -1376.324729053861,
	"heart_scale": -12.418702331243075,
	"ijcnn1": -1119.444551887607531,
	"news20.binary": -696.754382317548220,
	"phishing": -244.120822368800873,
	"rcv1_train.binary": -439.885381149034345,
	"yahoojp": -4117.569710134097477,
	"real-sim":-1140.218892969904800
    },
    "L2c10":{
    	"a1a": -6279.530312061115183,
    	"a9a": -137394.982701375702163,
    	"heart_scale": -1207.955529794298400,
    	"ijcnn1": -110050.998198572007823,
    	"news20.binary": -3383.543862333313882,
    	"phishing": -19671.482209781686834,
    	"rcv1_train.binary": -2943.364638952120004,
	"yahoojp": -104033.727589559814078,
	"real-sim": -13836.875530578598045
    },
    "L2c0.00390625":{
	"a9a": 54.898778338797818,
	"ijcnn1": 53.512353230548626,
	"news20.binary": 65.325470645634766,
	"rcv1_train.binary": 54.712015521627009,
	"real-sim": 138.481181429138871,
	"yahoojp": 239.686820890553378
    },
    "L2c0.00390625":{
	"a9a": -54.898778338797818,
	"ijcnn1": -53.512353230548626,
	"news20.binary": -65.325470645634766,
	"rcv1_train.binary": -54.712015521627009,
	"real-sim": -138.481181429138871,
	"yahoojp": -239.686820890553378
    },
    "L2c0.0078125":{
	"a9a": -108.775294204578714,
	"ijcnn1": -99.033083755195975,
	"news20.binary": -116.668204480315296,
	"rcv1_train.binary": -90.146748062907719,
	"real-sim": -222.425885862142394,
	"yahoojp": -441.563581606782179
    },
    "L2c0.015625":{
	"a9a": -216.351124429685939,
	"ijcnn1": -187.349534588965071,
	"news20.binary": -200.367017099675991,
	"rcv1_train.binary": -142.802552795678082,
	"real-sim": -351.517426993473066,
	"yahoojp": -813.971481519512167
    },
    "L2c0.03125":{
	"a9a": -431.292501258571122,
	"ijcnn1": -361.140736375657696,
	"news20.binary": -330.813983463643183,
	"rcv1_train.binary": -220.416776712844552,
	"real-sim": -549.303813357221429,
	"yahoojp": -1498.443585379727438
    },
    "L2c0.0625":{
	"a9a": -860.920227930710098,
	"ijcnn1": -706.261571875952086,
	"news20.binary": -523.835011833008252,
	"rcv1_train.binary": -334.457641102894911,
	"real-sim": -851.155668379751432,
	"yahoojp": -2745.435422664220368
    },
    "L2c0.125":{
	"a9a": -1719.890401311377218,
	"ijcnn1": -1394.724218121006970,
	"news20.binary": -791.349614888897236,
	"rcv1_train.binary": -499.569415384952492,
	"real-sim": -1308.111043933457950,
	"yahoojp": -4980.473670993681480
    },
    "L2c0.25":{
	"a9a": -3437.551392918615875,
	"ijcnn1": -2770.537714909431543,
	"news20.binary": -1133.582847857311208,
	"rcv1_train.binary": -730.881198927293099,
	"real-sim": -1989.919067061102623,
	"yahoojp": -8883.663651256469166
    },
    "L2c0.5":{
	"a9a": -6872.597884215801059,
	"ijcnn1": -5521.534140854936595,
	"news20.binary": -1531.471577739558143,
	"rcv1_train.binary": -1037.172520665836828,
	"real-sim": -2985.083459979493455,
	"yahoojp": -15449.736298034187712
    },
    "L2c2":{
	"a9a": -27481.708416428122291,
	"ijcnn1": -22026.327355073775834,
	"news20.binary": -2346.540995955459039,
	"rcv1_train.binary": -1834.611767736936827,
	"real-sim": -6333.455416994101142,
	"yahoojp": -41746.030456264306849
    },
    "L2c4":{
	"a9a": -54960.096983520641515,
	"ijcnn1": -44032.513023981635342,
	"news20.binary": -2735.986426563034911,
	"rcv1_train.binary": -2283.244225754563558,
	"real-sim": -8930.338796751411792,
	"yahoojp": -63853.477320735502872
    },
    "L2c8":{
	"a9a": -109916.701135227776831,
	"ijcnn1": -88044.839477867004462,
	"news20.binary": -3195.257905758150628,
	"rcv1_train.binary": -2769.288412378648900,
	"real-sim": -12434.688318295049612,
	"yahoojp": -92974.295920521661174
    },
    "L2c16":{
	"a9a": -219829.805222296039574,
	"ijcnn1": -176069.469845799118048,
	"news20.binary": -3900.525173059725148,
	"rcv1_train.binary": -3364.618088805269508,
	"real-sim": -17434.185508116705023,
	"yahoojp": -130738.077443597358069
    },
    "L2c32":{
	"a9a": -439655.956240487750620,
	"ijcnn1": -352118.719285823928658,
	"news20.binary": -5182.896292066204296,
	"rcv1_train.binary": -4230.336312488960175,
	"real-sim": -25259.515194111059827,
	"yahoojp": -182894.209371628297959
    },
    "L2c64":{
	"a9a": -879308.228298854199238,
	"ijcnn1": -704217.212511943303980,
	"news20.binary": -7672.547244318364392,
	"rcv1_train.binary": -5666.318037512346564,
	"real-sim": -38673.699751422020199,
	"yahoojp": -264169.677239460230339
    },
    "L2c128":{
	"a9a": -1758612.757052449276671,
	"ijcnn1": -1408414.196135769132525,
	"news20.binary": -12598.701005728493328,
	"rcv1_train.binary": -8220.630820635633427,
	"real-sim": -63110.718463130244345,
	"yahoojp": -405837.113369637809228
    },
    "L2c256":{
	"a9a": -3517221.806764540728182,
	"ijcnn1": -2816808.161968844477087,
	"news20.binary": -22388.319599819948053,
	"rcv1_train.binary": -12901.992960715813751,
	"real-sim": -109303.575283922182280,
	"yahoojp": -670914.682720793876797
    },
    "L2c512":{
	"a9a": -7034439.902311488054693,
	"ijcnn1": -5633596.092927607707679,
	"news20.binary": -41909.359917495683476,
	"rcv1_train.binary": -21672.794702927221806,
	"real-sim": -199004.633745041763177,
	"yahoojp": -1358076.966178216971457
    },
    "L2c1024":{
	"a9a": -14068876.091480575501919,
	"ijcnn1": -11267171.954491455107927,
	"news20.binary": -80789.998678168383776,
	"rcv1_train.binary": -38557.081880197503779,
	"real-sim": -375864.010833866486792,
	"yahoojp": -2287332.981667054817080
    },
    "L2c2048":{
	"a9a": -28137748.468770138919353,
	"ijcnn1": -22534323.677441917359829,
	"news20.binary": -158559.902321982255671,
	"rcv1_train.binary": -70956.610593336779857,
	"real-sim": -727087.785517628188245,
	"yahoojp": -4374941.380039221607149
    },
    "L2c4096":{
	"a9a": -56275493.222848773002625,
	"ijcnn1": -45068627.123254954814911,
	"news20.binary": -312871.115416924527381,
	"rcv1_train.binary": -134764.228345855779480,
	"real-sim": -1426363.841841514222324,
	"yahoojp": -9425674.325434481725097
    },
    "L2c8192":{
	"a9a": -112550982.730995744466782,
	"ijcnn1": -90137234.014836862683296,
	"news20.binary": -621371.247855326742865,
	"rcv1_train.binary": -261586.383875048835762,
	"real-sim": -2822565.899852479808033,
	"yahoojp": -17671666.476192910224199,
	"yahookr": -156874013.744594454,
	"heart_scale": -989244.883040017564781,
	"covtype.libsvm.binary": -403431.628650082,
	"covtype.libsvm.binary.scale": -3304653435.3253512,
    },
    "L2c16384":{
	"a9a": -225101961.746906965970993,
	"ijcnn1": -180274447.797976464033127,
	"news20.binary": -1238484.946269881678745,
	"rcv1_train.binary": -513405.296142977429554,
	"real-sim": -5607244.624500610865653,
	"yahoojp": -34485579.633781522512436
    },
    "biasL2c0.125":{
	"a9a": -1719.800997668895207,
	"ijcnn1": -1392.883689205032169,
	"news20.binary": -791.181560204815582,
	"rcv1_train.binary": -494.043381165942890,
	"real-sim": -1294.441610946594210,
    },
    "biasL2c0.25":{
	"a9a": -3437.460936031284746,
	"ijcnn1": -2768.738183005332303,
	"news20.binary": -1133.074509199425620,
	"rcv1_train.binary": -719.899668351338619,
	"real-sim": -1969.797156971597815,
    },
    "biasL2c0.5":{
	"a9a": -6872.505670026256666,
	"ijcnn1": -5519.756388105358383,
	"news20.binary": -1529.970321558319256,
	"rcv1_train.binary": -1017.546685940559200,
	"real-sim": -2953.786180694465656,
    },
    "biasL2c1":{
	"breast-cancer_scale": -59.43064203128018,
	"australian_scale": -270.2818072480945,
	"heart_scale": -114.9144550166213,
	"a1a": -637.7516571562145,
	"a9a": -13742.30343972653,
	"ijcnn1": -11021.42356208828,
	"news20.binary": -1943.053582171131,
	"rcv1_train.binary": -1380.777284961636,
	"real-sim": -4344.395126381159,
	"yahoojp": -23717.47272737786,
	"covtype.libsvm.binary":-4354.712314,
	"yahookr": -139947.251207661116496,
	"covtype.libsvm.binary.scale": -403422.79422253957#100,000
    },
    "biasL2c2":{
	"a9a": -27481.613283583359589,
	"ijcnn1": -22024.566605893283850,
	"news20.binary": -2334.075822110063200,
	"rcv1_train.binary": -1788.555839247801487,
	"real-sim": -6241.376872388566881,
    },
    "biasL2c4":{
	"a9a": -54960.001024990750011,
	"ijcnn1": -44030.755167771138076,
	"news20.binary": -2701.749896330909905,
	"rcv1_train.binary": -2222.278001758616483,
	"real-sim": -8760.226493137148282,
    },
    "biasL2c8":{
	"a9a": -109916.604693185145152,
	"ijcnn1": -88043.083073296264047,
	"news20.binary": -3107.509138781007550,
	"rcv1_train.binary": -2695.249658405296941,
	"real-sim": -12108.640672664554586,
    },
    "biasL2c16":{
	"a9a": -219829.708519455103669,
	"ijcnn1": -176067.714168953243643,
	"news20.binary": -3692.941022846368014,
	"rcv1_train.binary": -3281.055774396421384,
	"real-sim": -16805.424806433442427,
    },
    "biasL2c32":{
	"breast-cancer_scale": -1879.836599213785,
	"australian_scale": -8502.418933943414,
	"heart_scale": -3656.426734390486,
	"a1a": -20057.68582954491,
	"a9a": -439655.9494012388,
	"ijcnn1": -352116.9639732051,
	"news20.binary": -4723.837069179705,
	"rcv1_train.binary": -4140.82075436665,
	"real-sim": -24053.94159888949,
	"yahoojp":-173345.5089274935,
    },
    "biasL2c64":{
	"a9a": -879308.131385275861248,
	"ijcnn1": -704215.457381472107954,
	"news20.binary": -6701.982272093096981,
	"rcv1_train.binary": -5572.982197315027406,
	"real-sim": -36359.762257705486263,
    },
    "biasL2c128":{
	"breast-cancer_scale": -7517.100777544773,
	"a9a": -1758612.660095014143735,
	"ijcnn1": -1408412.441096334252506,
	"news20.binary": -10599.889527051478581,
	"rcv1_train.binary": -8123.214165212757507,
	"real-sim": -58661.907674048008630,
    },
    "biasL2c1024":{
	"breast-cancer_scale": -60131.54873761383,
    },
    "biasL2c0.03125":
    {
    	"breast-cancer_scale": -2.295708688115122,
    	"australian_scale": -9.024280329477792,
    	"heart_scale": -3.977967876909248,
    	"a1a": -21.76590059953335,
    	"a9a": -431.2070362905467,
	"ijcnn1": -359.1098274363673,
	"news20.binary": -330.7927242777877,
	"rcv1_train.binary": -219.4214901785117,
	"real-sim": -541.4970983808082,
	"yahoojp": -1238.511389715671,
    },
    "biasL1c0.03125":
    {
	"breast-cancer_scale": -2.276977998997026,
	"australian_scale": -6.749138721482132,
	"heart_scale": -3.720906294186185,
	"a1a": -19.86103562388033,
    	"a9a": -362.411649520055,
    	"ijcnn1": -289.1172415285611,
    	"news20.binary": -406.8329557870090,
    	"rcv1_train.binary": -286.5870788920642,
	"real-sim": -676.9317624201007,
    	"yahoojp": -973.7819315032502,
    },
    "biasL1c0.25":{
	"a9a": -2864.449193292481596,
	"ijcnn1": -2180.265734131292902,
	"news20.binary": -1525.056515180141332,
	"rcv1_train.binary": -915.792597368999282,
	"real-sim": -2451.940502592653047,
    },
    "biasL1c0.5":{
	"a9a": -5721.118992262800020,
	"ijcnn1": -4318.578282852491611,
	"news20.binary": -2075.398910304479614,
	"rcv1_train.binary": -1291.548787006288649,
	"real-sim": -3656.116728320576385,
    },
    "biasL1c1":{
	"breast-cancer_scale": -46.00399013664471,
	"australian_scale": -199.6504829745745,
	"heart_scale": -92.47337462016957,
	"a1a": -540.5750672979575,
	"a9a": -11433.387236617831149,
	"ijcnn1": -8590.159909824076749,
	"news20.binary": -2561.182971907652732,
	"rcv1_train.binary": -1745.667985918439186,
	"real-sim": -5345.122136598187353,
	"yahoojp":-23351.524127158492774,
	"yahookr":-143389.304898,
	"covtype.libsvm.binary": -4383.23183543,
	"covtype.libsvm.binary.scale":-337952.75811474106, #100,000
    },
    "biasL1c2":{
	"a9a": -22854.314487569292396,
	"ijcnn1": -17129.448134741887770,
	"news20.binary": -2857.686521048841314,
	"rcv1_train.binary": -2224.939490336275867,
	"real-sim": -7599.712544596704902,
	"covtype.libsvm.binary":-4383.2318354390 ,
    },
    "biasL1c4":{
	"a9a": -45694.182491092935379,
	"ijcnn1": -34205.631990751018748,
	"news20.binary": -3039.349619352877198,
	"rcv1_train.binary": -2642.510798770926613,
	"real-sim": -10433.096460461036258,
    },
    "biasL1c32":
    {
	"breast-cancer_scale": -1408.737046024731,
	"australian_scale": -6356.178279425147,
	"heart_scale": -3656.427734390489,
	"a1a": -16588.37269656481,
    	"a9a" : -365534.5682798115,
    	"ijcnn1" : -273293.1496785236,
    	"rcv1_train.binary": -4313.039575816452,
    	"news20.binary": -4522.361358841757,
    	"real-sim":-24858.86588074639,
    	"yahoojp":-172087.6931881657,
    },
    "biasL1c128":{
	"breast-cancer_scale": -5628.278659012924,
	"australian_scale": -25424.22358367619,
	"heart_scale": -11502.76505355036,
	"a1a": -66242.75397072639,
    },
    "biasL1c1024":{
	"breast-cancer_scale": -225628.578659012924,
	"australian_scale": -202409.9220992539,
	"heart_scale": -90911.38447171829,
	"a1a": -514633.4836791962,
    },
    "biasL1c65536":{
	"a1a": -2648889.203799085,
	"breast-cancer_scale": -399687.573817556,
    },
}

dlim = {
	"s41_c128_shrink":
	{
		"breast-cancer_scale": 0.1,
	},
	"s41_c32_shrink":
	{
		"a1a": 0.5,
		"australian_scale":0.2,
		"breast-cancer_scale": 5e-2,
		"heart_scale": 6e-2,
		"a9a": 3,
	},
	"s41_c0.03125_shrink":
	{
		"australian_scale": 0.1,
		"breast-cancer_scale": 1e-3,
		"heart_scale": 6e-4,
		"a1a": 0.015,
	},
	"s41_c1_shrink":
	{
		"a1a": 0.15,
		"breast-cancer_scale": 3e-3,
		"heart_scale": 4e-3,
		"australian_scale": 0.1,
		"a9a": 40,
		"ijcnn1": 5,
		"news20.binary": 20,
		"rcv1_train.binary": 20,
		"real-sim": 5,
		"yahoojp": 300
	},
	"s43_c1_shrink":
	{
		"a9a": 40,
		"ijcnn1": 5,
		"news20.binary": 20,
		"rcv1_train.binary": 20,
		"real-sim": 5,
		#"yahoojp": 1.5
	},
	"s43_c32_shrink":
	{
		"heart_scale": 6e-2,
		"a9a": 3,
	},
	"s42_c1_shrink":
	{
		"a1a":0.5,
	},
	"s45_c1_shrink":
	{
		"a9a": 40,
		"ijcnn1": 5,
		"news20.binary": 20,
		"rcv1_train.binary": 20,
		"real-sim": 5,
		#"yahoojp": 1.5
	},
	"s47_c1_shrink":
	{
		"a9a": 40,
		"ijcnn1": 5,
		"news20.binary": 20,
		"rcv1_train.binary": 20,
		"real-sim": 5,
		#"yahoojp": 1.5
	},
	"s51_c0.1_shrink":
	{
		"a9a": 4,
		"ijcnn1": 12,
		"news20.binary": 90,
		"rcv1_train.binary": 25,
		"real-sim": 120,
		"yahoojp": 1.5
	},
	"s52_c0.1_shrink":
	{
		"a9a": 0.2,
		"ijcnn1": 0.5,
		"news20.binary": 2,
		"rcv1_train.binary": 1,
		"real-sim": 1,
		"yahoojp": 1.5
	},
	"s55_c0.1_shrink":
	{
		"a9a": 0.2,
		"ijcnn1": 0.5,
		"news20.binary": 2,
		"rcv1_train.binary": 1,
		"real-sim": 1,
		"yahoojp": 1.5
	},
	"s56_c0.1_shrink":
	{
		"a9a": 0.2,
		"ijcnn1": 0.5,
		"news20.binary": 2,
		"rcv1_train.binary": 1,
		"real-sim": 1,
		"yahoojp": 1.5
	},
	"s53_c0.1_iter":
	{
		"a9a": 10000,
		"ijcnn1": 20000,
		"news20.binary": 10000,
		"rcv1_train.binary": 20000,
		"real-sim": 20000,
		"yahoojp": 40000
	},
	"s54_c0.1_iter":
	{
		"a9a": 10000,
		"ijcnn1": 20000,
		"news20.binary": 10000,
		"rcv1_train.binary": 20000,
		"real-sim": 20000,
		"yahoojp": 40000
	},
	"s51_c0.1_iter":
	{
		"a9a": 1000000,
		"ijcnn1": 2000000,
		"news20.binary": 1000000,
		"rcv1_train.binary": 2000000,
		"real-sim": 2000000,
		"yahoojp": 4000000
	},
	"s52_c0.1_iter":
	{
		"a9a": 10000,
		"ijcnn1": 200000,
		"news20.binary": 10000,
		"rcv1_train.binary": 20000,
		"real-sim": 20000,
		"yahoojp": 40000
	},
	"s35_c8192_iter":
	{
		"a9a":30000,
		"news20.binary":16000,
		"rcv1_train.binary":11000,
	},
	"s38_c8192_iter":
	{
		"a9a":30000,
		"news20.binary":30000,
		"rcv1_train.binary":300000,
	},
	"s34_c8192_iter":
	{
		"a9a":30000,
		"news20.binary":9500,
		"rcv1_train.binary":7500,
	},
	"s32_c8192_shrink":
	{
		"yahoojp":5000,
		"ijcnn1":150,
		"a9a":80,
		"real-sim":500,
		"rcv1_train.binary":250,
		"news20.binary":450,
	},
	"s32_c1_shrink":
	{
		"yahoojp":50,
		"ijcnn1":9,
		"a9a":10,
		"real-sim":1.5,
		"rcv1_train.binary":1,
		"news20.binary":35,
		"yahookr": 500,
	},
	"s21_c1_shrink":
	{
		"yahoojp":50,
		"ijcnn1":5,
		"a9a":5,
		"real-sim":1.5,
		"rcv1_train.binary":1,
		"news20.binary":20,
		"yahookr": 500,
	},
	"s22_c8192_shrink":
	{
		"yahoojp":4000,
		"ijcnn1":250,
		"a9a":150,
		"real-sim":600,
		"rcv1_train.binary":200,
		"news20.binary":370,
	},
	"s21_c8192_shrink":
	{
		"yahoojp":4000,
		"ijcnn1":150,
		"a9a":150,
		"real-sim":300,
		"rcv1_train.binary":200,
		"news20.binary":285,
	},
	"s22_c1_shrink":
	{
		"yahoojp":50,
		"ijcnn1":8,
		"a9a":8,
		"real-sim":1.5,
		"rcv1_train.binary":1,
		"news20.binary":34,
		"yahookr": 500,
	},
	"s26_c1_shrink":
	{
		"yahoojp":7,
		"ijcnn1":1,
		"a9a":1,
		"real-sim":1,
		"rcv1_train.binary":1,
		"news20.binary":2,
	},
	"s36_c8192_shrink":
	{
		"yahoojp":1000,
		"ijcnn1":180,
		"a9a":150,
		"real-sim":60,
		"rcv1_train.binary":15,
		"news20.binary":680,
	},
	"s1_c8192_shrink":
	{
		"yahoojp":1500,
		"ijcnn1":200,
		"a9a":150,
		"real-sim":100,
		"rcv1_train.binary":30,
		"news20.binary":800,
	},
	"s26_c8192_shrink":
	{
		"yahoojp":1500,
		"ijcnn1":200,
		"a9a":150,
		"real-sim":100,
		"rcv1_train.binary":30,
		"news20.binary":680,
	},
	"s36_c1_shrink":
	{
		"yahoojp":7,
		"ijcnn1":1,
		"a9a":1,
		"real-sim":1,
		"rcv1_train.binary":1,
		"news20.binary":2,
	},
	"s1_c1_shrink":
	{
		"yahoojp":7,
		"ijcnn1":1,
		"a9a":1,
		"real-sim":1,
		"rcv1_train.binary":1,
		"news20.binary":2,
		"covtype.libsvm.binary":1000,
	},
	"s22_c1_iter":
	{
		"yahoojp":250,
		"ijcnn1":1000,
		"a9a":1000,
		"real-sim":300,
		"rcv1_train.binary":500,
		"news20.binary":500,
	},
	"s21_c8192_iter":
	{
		"yahoojp":30000,
		"ijcnn1":30000,
		"a9a":30000,
		"real-sim":15000,
		"rcv1_train.binary":30000,
		"news20.binary":6800,
		"yahookr":25000,
	},
	"s32_c8192_iter":
	{
		"news20.binary":4240,
		"rcv1_train.binary":14300,
		"a9a":30000,
		"real-sim":9500 ,
	},
	"s22_c8192_iter":
	{
		"yahoojp":30000,
		"ijcnn1":30000,
		"a9a":30000,
		"real-sim":20000,
		"rcv1_train.binary":30000,
		"news20.binary":9000,
	},
	"s21_c1_iter":
	{
		"yahoojp":250,
		"ijcnn1":1000,
		"a9a":1000,
		"real-sim":300,
		"rcv1_train.binary":500,
		"news20.binary":500,
		"yahookr":1700
	},
	"s38_c1_iter":
	{
		"a9a":330,
		"news20.binary":20,
		"rcv1_train.binary":18,
		"heart_scale":290,
	},
	"s34_c1_iter":
	{
		"a9a":180,
		"news20.binary":20,
		"rcv1_train.binary":18,
		"heart_scale":140,
	},
	"s37_c1_iter":
	{
		"a9a":2000,
		"news20.binary":300,
		"rcv1_train.binary":2000,
		"heart_scale":460,
	},
	"s31_c1_iter":
	{
		"a9a":2000,
		"news20.binary":300,
		"rcv1_train.binary":2000,
		"heart_scale":240,
	},
	"s37_c8192_iter":
	{
		"a9a":30000,
		"news20.binary":6400,
		"rcv1_train.binary":7025,
	},
	"s31_c8192_iter":
	{
		"a9a":30000,
		"news20.binary":4150,
		"rcv1_train.binary":4325,
	},
	"s41_c32_iter":
	{
		"australian_scale":100000
	},
	"s41_c0.03125_iter":
	{
		"australian_scale":200000
	},
	"s41_c1_iter":
	{
		"australian_scale":200000,
		"ijcnn1":500,
		"a9a":500,
		"news20.binary":100,
		"real-sim":200,
		"rcv1_train.binary":200,
		"yahoojp":500,
	},
	"s32_c1_iter":
	{
		"ijcnn1":500,
		"a9a":500,
		"news20.binary":350,
		"real-sim":200,
		"rcv1_train.binary":270,
		"yahoojp":100,
		"heart_scale":200,
		"yahookr":900,
	},
	"s25_c8192_shrink":
	{
		"yahoojp":2000,
		"ijcnn1":220,
		"a9a":110,
		"real-sim":700,
		"rcv1_train.binary":200,
		"news20.binary":1000,
		"yahookr":20000,
	},
	"s35_c1_shrink":
	{
		"yahoojp":7,
		"ijcnn1":0.5,
		"a9a":1.6,
		"real-sim":1.3,
		"rcv1_train.binary":0.28,
		"news20.binary":2.5,
		"yahookr":47,
	},
	"s25_c1_shrink":
	{
		"yahoojp":7,
		"ijcnn1":0.5,
		"a9a":1.4,
		"real-sim":1,
		"rcv1_train.binary":0.25,
		"news20.binary":1.5,
		"yahookr":40,
	},
	"s35_c8192_iter":
	{
		"yahookr":10000,
	},
	"s24_c8192_shrink":
	{
		"yahoojp":2000,
		"ijcnn1":220,
		"a9a":110,
		"real-sim":700,
		"rcv1_train.binary":200,
		"news20.binary":1000,
		"yahookr":20000,
	},
	"s24_c1_shrink":
	{
		"yahoojp":3.0,
		"ijcnn1":0.4,
		"a9a":1,
		"real-sim":0.8,
		"rcv1_train.binary":0.12,
		"news20.binary":1,
		"yahookr":27,
	},
	"s35_c8192_shrink":
	{
		"yahoojp":2000,
		"ijcnn1":220,
		"a9a":110,
		"real-sim":700,
		"rcv1_train.binary":200,
		"news20.binary":1000,
		"yahookr":20000,
	},
	"s24_c0.125_iter":
	{
		"a9a":  100,
		"ijcnn1":  50,
		"news20.binary": 50 ,
		"rcv1_train.binary": 40,
		"real-sim":  50,
		"yahoojp":  50,
	},
	"s25_c0.125_iter":
	{
		"a9a":  100,
		"ijcnn1":  50,
		"news20.binary": 50 ,
		"rcv1_train.binary":  40,
		"real-sim":  50,
		"yahoojp":  50,
	},
	"s35_c0.125_iter":
	{
		"a9a":  40,
		"ijcnn1": 20 ,
		"news20.binary": 15 ,
		"rcv1_train.binary":  10,
		"real-sim":  20,
		"yahoojp":  50,
	},
	"s24_c0.25_iter":
	{
		"a9a":  150,
		"ijcnn1": 50 ,
		"news20.binary": 50 ,
		"rcv1_train.binary":  50,
		"real-sim":  50,
		"yahoojp":  50,
	},
	"s25_c0.25_iter":
	{
		"a9a":  150,
		"ijcnn1": 50 ,
		"news20.binary":  20,
		"rcv1_train.binary":  50,
		"real-sim":  50,
		"yahoojp":  50,
	},
	"s35_c0.25_iter":
	{
		"a9a":  60,
		"ijcnn1":  20,
		"news20.binary": 15 ,
		"rcv1_train.binary":  20,
		"real-sim":  20,
		"yahoojp":  50,
	},
	"s24_c0.5_iter":
	{
		"a9a":  250,
		"ijcnn1":  100,
		"news20.binary":  50,
		"rcv1_train.binary":  50,
		"real-sim": 50 ,
		"yahoojp":  50,
	},
	"s25_c0.5_iter":
	{
		"a9a": 250,
		"ijcnn1":  100,
		"news20.binary":  50,
		"rcv1_train.binary":  50,
		"real-sim":  50,
		"yahoojp":  50,
	},
	"s35_c0.5_iter":
	{
		"a9a":  100,
		"ijcnn1":  35,
		"news20.binary": 20 ,
		"rcv1_train.binary":  20,
		"real-sim":  20,
		"yahoojp":  50,
	},
	"s24_c0.1_iter":
	{
		"a9a":90,
		"a1a": 60,
		"heart_scale": 25,
		"phishing": 18,
    		"rcv1_train.binary": 13,
	},
	"s25_c0.1_iter":
	{
		"a9a":100,
		"a1a": 60,
		"heart_scale": 25,
		"phishing": 25,
    		"rcv1_train.binary": 15
	},
	"s34_c0.1_iter":
	{
		"a9a":100,
		"a1a": 50,
		"heart_scale": 25,
		"phishing": 25,
    		"rcv1_train.binary": 15
	},
	"s35_c0.1_iter":
	{
		"a9a":100,
		"a1a": 50,
		"heart_scale": 25,
		"phishing": 25,
    		"rcv1_train.binary": 15
	},

	"s24_c1_iter":
	{
		"a9a":400,
		"a1a": 500,
		"heart_scale": 360,
		"phishing": 50,
    		"rcv1_train.binary": 20,
		"ijcnn1": 80,
		"real-sim": 25,
		"yahoojp":30,
		"news20.binary":25,
		"yahookr":35,
	},

	"s25_c1_iter":
	{
		"a9a":400,
		"a1a": 500,
		"heart_scale": 360,
		"phishing": 50,
    		"rcv1_train.binary": 35,
		"real-sim": 40,
		"ijcnn1": 90,
		"yahoojp":40,
		"news20.binary":40,
		"yahookr":51,
	},
	"s35_c1_iter":
	{
		"a9a":180,
		"a1a": 500,
		"heart_scale": 200,
		"phishing": 40,
    		"rcv1_train.binary": 20,
		"news20.binary": 20,
		"yahoojp":25,
		"ijcnn1": 50,
		"real-sim": 20,
		"heart_scale": 140,
		"yahookr":25,
	},
	"s24_c10_iter":
	{
		"a9a":1000,
		"a1a": 1000,
		"heart_scale": 1000,
		"phishing": 400,
    		"rcv1_train.binary": 400
	},
	"s25_c10_iter":
	{
		"a9a":1000,
		"a1a": 1000,
		"heart_scale": 1000,
		"phishing": 400,
    		"rcv1_train.binary": 400
	},

	"s34_c10_iter":
	{
		"a9a":1000,
		"a1a": 1000,
		"heart_scale": 1000,
		"phishing": 200,
    		"rcv1_train.binary": 400
	},

	"s35_c10_iter":
	{
		"a9a":1000,
		"a1a": 1000,
		"heart_scale": 1000,
		"phishing": 200,
    		"rcv1_train.binary": 400
	},

	"s24_c2_iter":
	{
		"a9a": 800 ,
		"ijcnn1":  250,
		"news20.binary": 60 ,
		"rcv1_train.binary":  100,
		"real-sim": 75 ,
		"yahoojp":  100,
	},
	"s25_c2_iter":
	{
		"a9a": 800 ,
		"ijcnn1":  250,
		"news20.binary":  75,
		"rcv1_train.binary": 100 ,
		"real-sim": 100 ,
		"yahoojp":  75,
	},
	"s35_c2_iter":
	{
		"a9a": 350,
		"ijcnn1":  100,
		"news20.binary":  35,
		"rcv1_train.binary": 50 ,
		"real-sim":  50,
		"yahoojp":  50,
	},
	"s24_c4_iter":
	{
		"a9a": 1000 ,
		"ijcnn1":  400,
		"news20.binary":  100,
		"rcv1_train.binary": 100 ,
		"real-sim":  150,
		"yahoojp":  150,
	},
	"s25_c4_iter":
	{
		"a9a":  1000,
		"ijcnn1":  400,
		"news20.binary":  150,
		"rcv1_train.binary": 200 ,
		"real-sim":  175,
		"yahoojp": 150 ,
	},
	"s35_c4_iter":
	{
		"a9a": 750 ,
		"ijcnn1":  250,
		"news20.binary": 50,
		"rcv1_train.binary":  75,
		"real-sim":  100,
		"yahoojp": 100 ,
	},
	"s24_c8_iter":
	{
		"a9a":  1000,
		"ijcnn1":  800,
		"news20.binary":  200,
		"rcv1_train.binary": 200 ,
		"real-sim": 250 ,
		"yahoojp": 250 ,
	},
	"s25_c8_iter":
	{
		"a9a":  1000,
		"ijcnn1":  800,
		"news20.binary":  250,
		"rcv1_train.binary": 250 ,
		"real-sim": 230 ,
		"yahoojp":  250,
	},
	"s35_c8_iter":
	{
		"a9a":  1000,
		"ijcnn1":  400,
		"news20.binary": 125 ,
		"rcv1_train.binary": 150 ,
		"real-sim":  150,
		"yahoojp":  150,
	},
	"s24_c16_iter":
	{
		"a9a":  1000,
		"ijcnn1":  1000,
		"news20.binary":  400,
		"rcv1_train.binary":  400,
		"real-sim": 450 ,
		"yahoojp":  450,
	},
	"s25_c16_iter":
	{
		"a9a":  1000,
		"ijcnn1": 1000 ,
		"news20.binary":  500,
		"rcv1_train.binary": 450,
		"real-sim":  400,
		"yahoojp":  450,
	},
	"s35_c16_iter":
	{
		"a9a":  1000,
		"ijcnn1": 700 ,
		"news20.binary":  200,
		"rcv1_train.binary":  250,
		"real-sim": 250 ,
		"yahoojp":  250,
	},
	"s24_c32_iter":
	{
		"a9a":  1000,
		"ijcnn1":  1000,
		"news20.binary":  600,
		"rcv1_train.binary":  600,
		"real-sim":  800,
		"yahoojp":  600,
	},
	"s25_c32_iter":
	{
		"a9a":  1000,
		"ijcnn1": 1000 ,
		"news20.binary":  750,
		"rcv1_train.binary":  700,
		"real-sim": 750 ,
		"yahoojp":  750,
	},
	"s35_c32_iter":
	{
		"a9a":  1000,
		"ijcnn1": 1000 ,
		"news20.binary":  400,
		"rcv1_train.binary":  400,
		"real-sim":  400,
		"yahoojp":  400,
	},
	"s42_c64_iter":
	{
		"ijcnn1":2000,
		"a9a":2000,
		"news20.binary":2000,
		"real-sim":2000,
		"rcv1_train.binary":2000,
	},
	"s42_c4_iter":
	{
		"ijcnn1":500,
		"a9a":1300,
		"news20.binary":200,
		"real-sim":800,
		"rcv1_train.binary":450,
	},
	"s42_c8_iter":
	{
		"ijcnn1":800,
		"a9a":2000,
		"news20.binary":400,
		"real-sim":1600,
		"rcv1_train.binary":800,
	},
	"s35_c64_iter":
	{
		"ijcnn1":1000,
		"a9a":1000,
		"news20.binary":600,
		"real-sim":500,
		"rcv1_train.binary":500,
	},
	"s42_c128_iter":
	{
		"ijcnn1":2000,
		"a9a":2000,
		"news20.binary":2000,
		"real-sim":2000,
		"rcv1_train.binary":2000,
	},
	"s42_c0.25_iter":
	{
		"ijcnn1":90,
		"a9a":150,
		"news20.binary":40,
		"real-sim":100,
		"rcv1_train.binary":80,
	},
	"s35_c128_iter":
	{
		"ijcnn1":1000,
		"a9a":1000,
		"news20.binary":1000,
		"real-sim":1000,
		"rcv1_train.binary":1000,
	},
	"s42_c0.125_iter":
	{
		"ijcnn1":80,
		"a9a":100,
		"news20.binary":35,
		"real-sim":100,
		"rcv1_train.binary":60,
	},
	"s42_c2_iter":
	{
		"ijcnn1":250,
		"a9a":800,
		"news20.binary":100,
		"real-sim":400,
		"rcv1_train.binary":200,
	},
	"s42_c0.5_iter":
	{
		"ijcnn1":100,
		"a9a":250,
		"news20.binary":50,
		"real-sim":100,
		"rcv1_train.binary":100,
	},
	"s42_c1_iter":
	{
		"ijcnn1":140,
		"a9a":450,
		"news20.binary":50,
		"real-sim":200,
		"rcv1_train.binary":100,
		"yahoojp":100,
	},
	"s42_c16_iter":
	{
		"ijcnn1":1500,
		"a9a":2000,
		"news20.binary":700,
		"real-sim":2000,
		"rcv1_train.binary":2000,
	},
	"s42_c32_iter":
	{
		"ijcnn1":2000,
		"a9a":2000,
		"news20.binary":1500,
		"real-sim":2000,
		"rcv1_train.binary":2000,
	},

	#######for l1 shrinking#######

	"s3_c1_shrink":
	{
		"yahoojp":5,
		"ijcnn1":0.079,
		"a9a":0.20,
		"real-sim":0.7,
		"rcv1_train.binary":1,
		"news20.binary":5,
		"covtype.libsvm.binary": 287,
		"yahookr":500,
	},
	"s33_c1_shrink":
	{
		"yahoojp":5,
		"ijcnn1":0.25,
		"a9a":5,
		"real-sim":0.7,
		"rcv1_train.binary":1,
		"news20.binary":5,
		"covtype.libsvm.binary":287,
		"yahookr":500,
	},
	"s39_c1_shrink":
	{
		"yahoojp":5,
		"ijcnn1":0.25,
		"a9a":5,
		"real-sim":0.7,
		"rcv1_train.binary":1,
		"news20.binary":5,
		"covtype.libsvm.binary":287,
		"yahookr":500,
	},
	"s32_c8192_shrink":
	{
		"yahoojp":1500,
		"ijcnn1":200,
		"a9a":200,
		"real-sim":500,
		"rcv1_train.binary":150,
		"news20.binary":500,
		"yahookr":14500,
	},
	"s3_c8192_shrink":
	{
		"yahoojp":1500,
		"ijcnn1":200,
		"a9a":200,
		"real-sim":500,
		"rcv1_train.binary":50,
		"news20.binary":500,
		"yahookr":14500,
	},
	"s32_c1_shrink":
	{
		"yahoojp":5,
		"ijcnn1":1,
		"a9a":2,
		"real-sim":0.7,
		"rcv1_train.binary":1,
		"news20.binary":5,
		"covtype.libsvm.binary":287,
		"yahookr":500,
	},
	"s21_c1_shrink":
	{
		"yahoojp":5,
		"ijcnn1":1,
		"a9a":2,
		"real-sim":0.7,
		"rcv1_train.binary":1,
		"news20.binary":5,
		"covtype.libsvm.binary":287,
		"yahookr":500,
	},
	"s33_c8192_shrink":
	{
		"yahoojp":1500,
		"ijcnn1":200,
		"a9a":200,
		"real-sim":500,
		"rcv1_train.binary":50,
		"news20.binary":500,
		"yahookr":14500,
	},
	"s21_c8192_shrink":
	{
		"yahoojp":1500,
		"ijcnn1":200,
		"a9a":200,
		"real-sim":500,
		"rcv1_train.binary":150,
		"news20.binary":285,
		"yahookr":14500,
	},
}

dstart = {
    "3": {
	"rcv1_test": (1e-8, 1e-1),
	"covtype_scale": (1e-7, 1e0),
	"a9a": (1e-1, 1e0),
	"australian_scale": (1e-8, 1e-1),
	"epsilon_normalized": (1e-7, 1e-1),
	"gisette_scale": (1e-7, 1e-1),
	"ijcnn1.t": (1e-8, 1e-1),
	"HIGGS": (1e-5, 1e0),
	"kdda": (1e-6, 1e-1),
	"kddb": (1e-7, 1e-1),
	"news20": (1e-8, 1e-1),
	"SUSY": (1e-5, 1e0),
	"url_combined": (1e-4, 1e-1),
	"url": (1e-4, 1e-1),
	"webspam_wc_normalized_unigram": (1e-8, 1e-1),
	"webspam_uni": (1e-8, 1e-1),
	"webspam": (1e-7, 1e-1),
	"yahoojp": (1e-8, 1e-2),
	"yahookr": (1e-8, 1e-2),
	"splice_scale": (1e-8, 1e-1),
	"heart_scale": (1e-8, 1e-1)
    },
    "1": {
	"rcv1_test": (1e-8, 1e-1),
	"covtype_scale": (1e-7, 1e-2),
	"a9a": (1e-8, 1e-1),
	"australian_scale": (1e-8, 1e-1),
	"epsilon_normalized": (1e-8, 1e0),
	"gisette_scale": (1e-7, 1e-1),
	"ijcnn1.t": (1e-8, 1e-1),
	"HIGGS": (1e-5, 1e-1),
	"kdda": (1e-6, 1e-1),
	"kddb": (1e-8, 1e-1),
	"news20": (1e-8, 1e-1),
	"SUSY": (1e-5, 1e-1),
	"url_combined": (1e-6, 1e-1),
	"url": (1e-6, 1e-1),
	"webspam_wc_normalized_unigram": (1e-8, 1e-1),
	"webspam_uni": (1e-8, 1e-1),
	"webspam": (1e-7, 1e-1),
	"yahoojp": (1e-8, 1e-2),
	"yahookr": (1e-8, 1e-2),
	"splice_scale": (1e-8, 1e-1),
	"heart_scale": (1e-8, 1e-1)
    },
}

bestC = {
    "3": {
	"rcv1_test": 2,
	"covtype_scale": 16,
	"epsilon_normalized": 16,
	"HIGGS": 0.5,
	"kddb": 0.25,
	"news20": 2,
	"url_combined": 2,
	"url": 2,
	"webspam_wc_normalized_unigram": 128,
	"webspam_uni": 128,
	"webspam": 32,
	"yahoojp": 2,
	"yahookr": 8,
    },
    "1": {
	"rcv1_test": 0.5,
	"covtype_scale": 0.015625,
	"epsilon_normalized": 4,
	"HIGGS": 0.5,
	"kddb": 0.03125,
	"news20": 8,
	"url_combined": 2,
	"url": 2,
	"webspam_wc_normalized_unigram": 128,
	"webspam_uni": 128,
	"webspam": 32,
	"yahoojp": 0.5,
	"yahookr": 2,
	"a1a": 0.5,
	"a9a": 0.5,
	"ijcnn1": 16,
	"news20.binary": 32,
	"rcv1_train.binary": 1,
	"real-sim": 1,
    }
}

legend = {
    "eps5new": "Alg 4",
    "pass": "Async-CD",
    "ori": "LIBLINEAR"
}
