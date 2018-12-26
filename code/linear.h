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
	ONE_L2_CY_SH,
	L2R_L2LOSS_SVC, 
	ONE_L1_CY_SH,
	MCSVM_CS,
	L1R_L2LOSS_SVC,
	L1R_LR,
	L2R_LR_DUAL,
	//for regression
	L2R_L2LOSS_SVR = 11,
	L2R_L2LOSS_SVR_DUAL = 12,
	L2R_L1LOSS_SVR_DUAL = 13,
	//for one-variable
	ONE_L1_CY_1000 = 21,
	ONE_L1_RD_1000 = 22,
	ONE_L1_RD_SH = 23,

	ONE_L2_CY_1000 = 24,
	ONE_L2_RD_1000 = 25,
	ONE_L2_RD_SH = 26,

	ONE_L1_SEMIGD_1000 = 27,
	ONE_L1_SEMIGD_SH = 28,
	ONE_L2_SEMIGD_1000 = 29,
	ONE_L2_SEMIGD_SH = 30,

	//for two-variable
	TWO_L1_CY_1000 = 31,
	TWO_L1_RD_1000 = 32,
	TWO_L1_RD_SH = 33,
	
	TWO_L2_CY_1000 = 34,
	TWO_L2_RD_1000 = 35,
	TWO_L2_RD_SH = 36,
	TWO_L1_SEMICY_1000 = 37,
	TWO_L2_SEMICY_1000 = 38,
	TWO_L1_SEMIRDONE_1000 = 371,
	TWO_L2_SEMIRDONE_1000 = 381,
	TWO_L1_SEMIRDTWO_1000 = 372,
	TWO_L2_SEMIRDTWO_1000 = 382,
	TWO_L1_RD_SH2 = 39,
	TWO_L2_RD_SH2 = 40,
	TWO_L1_SEMIGD_1000 = 931,
	TWO_L1_SEMIGD_SH = 932,
	TWO_L2_SEMIGD_1000 = 933,
	TWO_L2_SEMIGD_SH = 934,


	//for two-variable linear constraints
	BIAS_L1_RD_1000 = 41,
	BIAS_L2_RD_1000 = 42,
	BIAS_L1_SEMIGD_1000 = 43,
	BIAS_L2_SEMIGD_1000 = 44,
	BIAS_L1_RD_SH = 45,
	BIAS_L2_RD_SH = 46,
	BIAS_L1_SEMIGD_SH = 47,
	BIAS_L2_SEMIGD_SH = 48,

	//for one-class svm
	ONECLASS_L1_RD_1000 = 51,
	ONECLASS_L1_SEMIGD_1000 = 52,
	ONECLASS_L1_FIRST_1000 = 53,
	ONECLASS_L1_SECOND_1000 = 54,
	ONECLASS_L1_RD_SH = 55,
	ONECLASS_L1_SEMIGD_SH = 56,
	ONECLASS_L2_RD_1000 = 57,
	ONECLASS_L2_SEMIGD_1000 = 58,
	ONECLASS_L2_FIRST_1000 = 59,
	ONECLASS_L2_SECOND_1000 = 60,
	ONECLASS_L2_RD_SH = 61,
	ONECLASS_L2_SEMIGD_SH = 62,
	}; /* solver_type */

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
	double nu;	/* for one-class formulation */
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

