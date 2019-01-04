runs = {
	"L2R_LR" : 0,
	"ONE_L2_CY_SH" : 1,
	"L2R_L2LOSS_SVC" : 2,
	"ONE_L1_CY_SH" : 3,
	"MCSVM_CS" : 4,
	"L1R_L2LOSS_SVC" : 5,
	"L1R_LR" : 6,
	"L2R_LR_DUAL" : 7,
	"L2R_L2LOSS_SVR" : 11,
	"L2R_L2LOSS_SVR_DUAL" : 12,
	"L2R_L1LOSS_SVR_DUAL" : 13,

	"ONE_L1_CY_1000" : 20111,
	"ONE_L2_CY_1000" : 20121,

	"ONE_L1_RD_1000" : 20211,
	"ONE_L1_RD_SH" : 20212,
	"ONE_L2_RD_1000" : 20221,
	"ONE_L2_RD_SH" : 20122,

	"ONE_L1_SEMIGD_1000" : 20311,
	"ONE_L1_SEMIGD_SH" : 20312,
	"ONE_L2_SEMIGD_1000" : 20321,
	"ONE_L2_SEMIGD_SH" : 20322,

	"ONE_L1_SEMIGD_DUALOBJ_1000" : 20411,
	"ONE_L1_SEMIGD_DUALOBJ_SH" : 20412,
	"ONE_L2_SEMIGD_DUALOBJ_1000" : 20421,
	"ONE_L2_SEMIGD_DUALOBJ_SH" : 20422,

	"TWO_L1_CY_1000" : 30111,
	"TWO_L2_CY_1000" : 30121,

	"TWO_L1_RD_1000" : 30211,
	"TWO_L1_RD_SH" : 30212,
	"TWO_L1_RD_SH2" : 30213,
	"TWO_L2_RD_1000" : 30221,
	"TWO_L2_RD_SH" : 30222,
	"TWO_L2_RD_SH2" : 30223,

	"TWO_L1_SEMICY_1000" : 30311,
	"TWO_L2_SEMICY_1000" : 30321,
	"TWO_L1_SEMIRDONE_1000" : 30411,
	"TWO_L2_SEMIRDONE_1000" : 30421,
	"TWO_L1_SEMIRDTWO_1000" : 30511,
	"TWO_L2_SEMIRDTWO_1000" : 30521,
	"TWO_L1_SEMIGD_1000" : 30611,
	"TWO_L1_SEMIGD_SH" : 30612,
	"TWO_L2_SEMIGD_1000" : 30621,
	"TWO_L2_SEMIGD_SH" : 30622,

	"BIAS_L1_RD_1000" : 40111,
	"BIAS_L1_RD_SH" : 40112,
	"BIAS_L2_RD_1000" : 40121,
	"BIAS_L2_RD_SH" : 40122,

	"BIAS_L1_SEMIGD_1000" : 40211,
	"BIAS_L1_SEMIGD_SH" : 40212,
	"BIAS_L2_SEMIGD_1000" : 40221,
	"BIAS_L2_SEMIGD_SH" : 40222,

	"ONECLASS_L1_RD_1000" : 50111,
	"ONECLASS_L1_RD_SH" : 50112,
	"ONECLASS_L2_RD_1000" : 50121,
	"ONECLASS_L2_RD_SH" : 50222,

	"ONECLASS_L1_SEMIGD_1000" : 50211,
	"ONECLASS_L1_SEMIGD_SH" : 50212,
	"ONECLASS_L2_SEMIGD_1000" : 50221,
	"ONECLASS_L2_SEMIGD_SH" : 50222,

	"ONECLASS_L1_FIRST_1000" : 50311,
	"ONECLASS_L2_FIRST_1000" : 50321,
	"ONECLASS_L1_SECOND_1000" : 50411,
	"ONECLASS_L2_SECOND_1000" : 50421,
}
alltype = runs.values()
def is_biasobj(code):
	str_code = str(code)
	if len(str_code) != 5:
		return False
	else:
		return str_code[0] == "4"

def is_semigd(code):
	str_code = str(code)
	if len(str_code) != 5:
		return False
	else:
		semigd_prefix = [203, 204, 402, 502]
		return str_code[:3] in map(str, semigd_prefix)

def is_shrink(code):
	str_code = str(code)
	if len(str_code) != 5:
		return str_code in map(str, [1,3])
	else:
		return str_code[4] != "1"

def is_oneclass(code):
	str_code = str(code)
	if len(str_code) != 5:
		return False
	else:
		return str_code[0] == "5"

def is_L1(code):
	str_code = str(code)
	if len(str_code) != 5:
		return str_code in map(str, [3])
	else:
		return str_code[3] == "1"
