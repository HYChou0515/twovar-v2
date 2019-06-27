#ifndef _LIBLINEAR_H
#define _LIBLINEAR_H

#ifdef __cplusplus
extern "C" {
#endif

struct feature_node
{
	int index;
	double value;
};

struct problem
{
	int l, n;
	double *y;
	struct feature_node **x;
	double bias;            /* < 0 if no bias term */  
};

enum { 
	L2R_LR, 
	OLD_ONE_L2_CY_SH,
	L2R_L2LOSS_SVC, 
	OLD_ONE_L1_CY_SH,
	MCSVM_CS,
	L1R_L2LOSS_SVC,
	L1R_LR,
	L2R_LR_DUAL,
	//for regression
	L2R_L2LOSS_SVR = 11,
	L2R_L2LOSS_SVR_DUAL = 12,
	L2R_L1LOSS_SVR_DUAL = 13,
	// code=abcde
	// a: big category
	// bc: algo
	// d: loss 1 or 2
	// e: 1000 iter or shrink
	//for one-variable
	ONE_L1_CY_1000 = 20111,
	ONE_L1_CY_SH = 20112,
	ONE_L2_CY_1000 = 20121,
	ONE_L2_CY_SH = 20122,

	ONE_L1_RD_1000 = 20211,
	ONE_L1_RD_SH = 20212,
	ONE_L2_RD_1000 = 20221,
	ONE_L2_RD_SH = 20222,

	ONE_L1_SEMIGD_1000 = 20311,
	ONE_L1_SEMIGD_SH = 20312,
	ONE_L2_SEMIGD_1000 = 20321,
	ONE_L2_SEMIGD_SH = 20322,

	ONE_L1_SEMIGD_DUALOBJ_1000 = 20411,
	ONE_L1_SEMIGD_DUALOBJ_SH = 20412,
	ONE_L2_SEMIGD_DUALOBJ_1000 = 20421,
	ONE_L2_SEMIGD_DUALOBJ_SH = 20422,

	ONE_L1_SEMIGD_RAND_1000 = 20511,
	ONE_L1_SEMIGD_RAND_SH = 20512,
	ONE_L2_SEMIGD_RAND_1000 = 20521,
	ONE_L2_SEMIGD_RAND_SH = 20522,

	ONE_L1_SEMIGD_DUALOBJ_RAND_1000 = 20611,
	ONE_L1_SEMIGD_DUALOBJ_RAND_SH = 20612,
	ONE_L2_SEMIGD_DUALOBJ_RAND_1000 = 20621,
	ONE_L2_SEMIGD_DUALOBJ_RAND_SH = 20622,

	ONE_L1_SEMIGD_DUALOBJ_YBAL_1000 = 20711,
	ONE_L1_SEMIGD_DUALOBJ_YBAL_SH = 20712,
	ONE_L2_SEMIGD_DUALOBJ_YBAL_1000 = 20721,
	ONE_L2_SEMIGD_DUALOBJ_YBAL_SH = 20722,

	//for two-variable
	TWO_L1_CY_1000 = 30111,
	TWO_L2_CY_1000 = 30121,

	TWO_L1_RD_1000 = 30211,
	TWO_L1_RD_SH = 30212,
	TWO_L1_RD_SH2 = 30213,
	TWO_L2_RD_1000 = 30221,
	TWO_L2_RD_SH = 30222,
	TWO_L2_RD_SH2 = 30223,

	TWO_L1_SEMICY_1000 = 30311,
	TWO_L2_SEMICY_1000 = 30321,
	TWO_L1_SEMIRDONE_1000 = 30411,
	TWO_L2_SEMIRDONE_1000 = 30421,
	TWO_L1_SEMIRDTWO_1000 = 30511,
	TWO_L2_SEMIRDTWO_1000 = 30521,
	TWO_L1_SEMIGD_1000 = 30611,
	TWO_L1_SEMIGD_SH = 30612,
	TWO_L2_SEMIGD_1000 = 30621,
	TWO_L2_SEMIGD_SH = 30622,

	//for two-variable linear constraints
	BIAS_L1_RD_1000 = 40111,
	BIAS_L1_RD_SH = 40112,
	BIAS_L2_RD_1000 = 40121,
	BIAS_L2_RD_SH = 40122,

	BIAS_L1_SEMIGD_1000 = 40211,
	BIAS_L1_SEMIGD_SH = 40212,
	BIAS_L2_SEMIGD_1000 = 40221,
	BIAS_L2_SEMIGD_SH = 40222,

	BIAS_L1_SEMIGD_RAND_1000 = 40311,
	BIAS_L1_SEMIGD_RAND_SH = 40312,
	BIAS_L2_SEMIGD_RAND_1000 = 40321,
	BIAS_L2_SEMIGD_RAND_SH = 40322,

	BIAS_L1_SEMIGD_CY_FIRST_1000 = 40411,
	BIAS_L1_SEMIGD_CY_FIRST_SH = 40412,
	BIAS_L2_SEMIGD_CY_FIRST_1000 = 40421,
	BIAS_L2_SEMIGD_CY_FIRST_SH = 40422,

	BIAS_L1_SEMIGD_RD_FIRST_1000 = 40511,
	BIAS_L1_SEMIGD_RD_FIRST_SH = 40512,
	BIAS_L2_SEMIGD_RD_FIRST_1000 = 40521,
	BIAS_L2_SEMIGD_RD_FIRST_SH = 40522,

	BIAS_L1_CY_1000 = 40611,
	BIAS_L1_CY_SH = 40612,
	BIAS_L2_CY_1000 = 40621,
	BIAS_L2_CY_SH = 40622,

	BIAS_L1_SEMIGD_CY_DUALOBJ_1000 = 40711,
	BIAS_L1_SEMIGD_CY_DUALOBJ_SH = 40712,
	BIAS_L2_SEMIGD_CY_DUALOBJ_1000 = 40721,
	BIAS_L2_SEMIGD_CY_DUALOBJ_SH = 40722,

	BIAS_L1_SEMIGD_RD_DUALOBJ_1000 = 40811,
	BIAS_L1_SEMIGD_RD_DUALOBJ_SH = 40812,
	BIAS_L2_SEMIGD_RD_DUALOBJ_1000 = 40821,
	BIAS_L2_SEMIGD_RD_DUALOBJ_SH = 40822,

	//for one-class svm
	ONECLASS_L1_RD_1000 = 50111,
	ONECLASS_L1_RD_SH = 50112,

	ONECLASS_L1_SEMIGD_1000 = 50211,
	ONECLASS_L1_SEMIGD_SH = 50212,

	ONECLASS_L1_FIRST_1000 = 50311,
	ONECLASS_L1_SECOND_1000 = 50411,

	ONECLASS_L1_SEMIGD_RAND_1000 = 50511,
	ONECLASS_L1_SEMIGD_RAND_SH = 50512,

	ONECLASS_L1_SEMIGD_CY_FIRST_1000 = 50611,
	ONECLASS_L1_SEMIGD_CY_FIRST_SH = 50612,

	ONECLASS_L1_SEMIGD_RD_FIRST_1000 = 50711,
	ONECLASS_L1_SEMIGD_RD_FIRST_SH = 50712,

	ONECLASS_L1_CY_1000 = 50811,
	ONECLASS_L1_CY_SH = 50812,

	ONECLASS_L1_SEMIGD_CY_DUALOBJ_1000 = 50911,
	ONECLASS_L1_SEMIGD_CY_DUALOBJ_SH = 50912,

	ONECLASS_L1_SEMIGD_RD_DUALOBJ_1000 = 51011,
	ONECLASS_L1_SEMIGD_RD_DUALOBJ_SH = 51012,
	}; /* solver_type */

struct resume
{
	char fname[1024];
	bool read_resume;
	int iter;
	unsigned long long cdsteps;
	unsigned long long nr_n_ops;
	clock_t duration;
	unsigned long long nr_rand_calls;
	double last_obj;
	int active_size;
	double PGmax_old;
	double PGmin_old;
	double Gmax_old;
	double Gmin_old;
	int alpha_size;
	double *alpha;
	int *index;
	int w_size;
	double *w;
};

struct parameter
{
	int solver_type;
	/* these are for training only */
	double eps;	        /* stopping criteria */
	double C;
	int nr_weight;
	int *weight_label;
	double* weight;
	double p;
	double *init_sol;
	double r;       /* update ratio for semi-gd*/ 
	int max_iter;     
	int timeout; // in second     
	int max_cdstep; // this *= prob->l
	double opt_val;
	double nu;	/* for one-class formulation */

	struct resume *_resume;
	FILE* log_fp;
};

struct model
{
	struct parameter param;
	int nr_class;		/* number of classes */
	int nr_feature;
	double *w;
	int *label;		/* label of each class */
	double bias;
};

struct model* train(const struct problem *prob, const struct parameter *param);
void cross_validation(const struct problem *prob, const struct parameter *param, int nr_fold, double *target);
void find_parameter_C(const struct problem *prob, const struct parameter *param, int nr_fold, double start_C, double max_C, double *best_C, double *best_rate);

double predict_values(const struct model *model_, const struct feature_node *x, double* dec_values);
double predict(const struct model *model_, const struct feature_node *x);
double predict_probability(const struct model *model_, const struct feature_node *x, double* prob_estimates);

struct resume *load_resume(const char *resume_file_name);

int save_model(const char *model_file_name, const struct model *model_);
struct model *load_model(const char *model_file_name);

int get_nr_feature(const struct model *model_);
int get_nr_class(const struct model *model_);
void get_labels(const struct model *model_, int* label);
double get_decfun_coef(const struct model *model_, int feat_idx, int label_idx);
double get_decfun_bias(const struct model *model_, int label_idx);

void free_model_content(struct model *model_ptr);
void free_and_destroy_model(struct model **model_ptr_ptr);
void destroy_param(struct parameter *param);

const char *check_parameter(const struct problem *prob, const struct parameter *param);
int check_probability_model(const struct model *model);
int check_regression_model(const struct model *model);
void set_print_string_function(void (*print_func) (const char*));

#ifdef __cplusplus
}
#endif

#endif /* _LIBLINEAR_H */

