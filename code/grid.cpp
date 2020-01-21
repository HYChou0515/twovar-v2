#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <omp.h>
#include "linear.h"
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define INF HUGE_VAL

void print_null(const char *s) {}

void exit_with_help()
{
	printf(
	"Usage: train training_set_file [log_file] [resume_file] [model_file]\n"
	"options are from stdin\n"
	"options:\n"
	"-s type : set type of solver (default 1)\n"
	"  for multi-class classification\n"
	"	 0 -- L2-regularized logistic regression (primal)\n"
	"	 1 -- L2-regularized L2-loss support vector classification (dual)\n"
	"	 2 -- L2-regularized L2-loss support vector classification (primal)\n"
	"	 3 -- L2-regularized L1-loss support vector classification (dual)\n"
	"	 4 -- support vector classification by Crammer and Singer\n"
	"	 5 -- L1-regularized L2-loss support vector classification\n"
	"	 6 -- L1-regularized logistic regression\n"
	"	 7 -- L2-regularized logistic regression (dual)\n"
	"  for regression\n"
	"	11 -- L2-regularized L2-loss support vector regression (primal)\n"
	"	12 -- L2-regularized L2-loss support vector regression (dual)\n"
	"	13 -- L2-regularized L1-loss support vector regression (dual)\n"
	"-c cost : set the parameter C (default 1)\n"
	"-p epsilon : set the epsilon in loss function of SVR (default 0.1)\n"
	"-e epsilon : set tolerance of termination criterion\n"
	"	-s 0 and 2\n"
	"		|f'(w)|_2 <= eps*min(pos,neg)/l*|f'(w0)|_2,\n"
	"		where f is the primal function and pos/neg are # of\n"
	"		positive/negative data (default 0.01)\n"
	"	-s 11\n"
	"		|f'(w)|_2 <= eps*|f'(w0)|_2 (default 0.001)\n"
	"	-s 1, 3, 4, and 7\n"
	"		Dual maximal violation <= eps; similar to libsvm (default 0.1)\n"
	"	-s 5 and 6\n"
	"		|f'(w)|_1 <= eps*min(pos,neg)/l*|f'(w0)|_1,\n"
	"		where f is the primal function (default 0.01)\n"
	"	-s 12 and 13\n"
	"		|f'(alpha)|_1 <= eps |f'(alpha0)|,\n"
	"		where f is the dual function (default 0.1)\n"
	"-B bias : if bias >= 0, instance x becomes [x; bias]; if < 0, no bias term added (default -1)\n"
	"-wi weight: weights adjust the parameter C of different classes (see README for details)\n"
	"-v n: n-fold cross validation mode\n"
	"-C : find parameter C (only for -s 0 and 2)\n"
	"-q : quiet mode (no outputs)\n"
	"-m : max iteration\n"
	"-S : max order-n operations\n"
	"-t : timeout (in second)\n"
	"-o : minimum objective value\n"
	"L2R_LR, \n"
	"OLD_ONE_L2_CY_SH,\n"
	"L2R_L2LOSS_SVC, \n"
	"OLD_ONE_L1_CY_SH,\n"
	"MCSVM_CS,\n"
	"L1R_L2LOSS_SVC,\n"
	"L1R_LR,\n"
	"L2R_LR_DUAL,\n"
	"//for regression\n"
	"L2R_L2LOSS_SVR = 11,\n"
	"L2R_L2LOSS_SVR_DUAL = 12,\n"
	"L2R_L1LOSS_SVR_DUAL = 13,\n"
	"// code=abcde\n"
	"// a: big category\n"
	"// bc: algo\n"
	"// d: loss 1 or 2\n"
	"// e: 1000 iter or shrink\n"
	"//for one-variable\n"
	"ONE_L1_CY_1000 = 20111,\n"
	"ONE_L1_CY_SH = 20112,\n"
	"ONE_L2_CY_1000 = 20121,\n"
	"ONE_L2_CY_SH = 20122,\n"
	"ONE_L1_RD_1000 = 20211,\n"
	"ONE_L1_RD_SH = 20212,\n"
	"ONE_L2_RD_1000 = 20221,\n"
	"ONE_L2_RD_SH = 20222,\n"
	"ONE_L1_SEMIGD_1000 = 20311,\n"
	"ONE_L1_SEMIGD_SH = 20312,\n"
	"ONE_L2_SEMIGD_1000 = 20321,\n"
	"ONE_L2_SEMIGD_SH = 20322,\n"
	"ONE_L1_SEMIGD_RAND_1000 = 20511,\n"
	"ONE_L1_SEMIGD_RAND_SH = 20512,\n"
	"ONE_L2_SEMIGD_RAND_1000 = 20521,\n"
	"ONE_L2_SEMIGD_RAND_SH = 20522,\n"
	"ONE_L1_SEMIGD_DUALOBJ_1000 = 20411,\n"
	"ONE_L1_SEMIGD_DUALOBJ_SH = 20412,\n"
	"ONE_L2_SEMIGD_DUALOBJ_1000 = 20421,\n"
	"ONE_L2_SEMIGD_DUALOBJ_SH = 20422,\n"
	"ONE_L1_SEMIGD_DUALOBJ_RAND_1000 = 20611,\n"
	"ONE_L1_SEMIGD_DUALOBJ_RAND_SH = 20612,\n"
	"ONE_L2_SEMIGD_DUALOBJ_RAND_1000 = 20621,\n"
	"ONE_L2_SEMIGD_DUALOBJ_RAND_SH = 20622,\n"
	"ONE_L1_SEMIGD_DUALOBJ_YBAL_1000 = 20711,\n"
	"ONE_L1_SEMIGD_DUALOBJ_YBAL_SH = 20712,\n"
	"ONE_L2_SEMIGD_DUALOBJ_YBAL_1000 = 20721,\n"
	"ONE_L2_SEMIGD_DUALOBJ_YBAL_SH = 20722,\n"
	"//for two-variable\n"
	"TWO_L1_CY_1000 = 30111,\n"
	"TWO_L2_CY_1000 = 30121,\n"
	"TWO_L1_RD_1000 = 30211,\n"
	"TWO_L1_RD_SH = 30212,\n"
	"TWO_L1_RD_SH2 = 30213,\n"
	"TWO_L2_RD_1000 = 30221,\n"
	"TWO_L2_RD_SH = 30222,\n"
	"TWO_L2_RD_SH2 = 30223,\n"
	"TWO_L1_SEMICY_1000 = 30311,\n"
	"TWO_L2_SEMICY_1000 = 30321,\n"
	"TWO_L1_SEMIRDONE_1000 = 30411,\n"
	"TWO_L2_SEMIRDONE_1000 = 30421,\n"
	"TWO_L1_SEMIRDTWO_1000 = 30511,\n"
	"TWO_L2_SEMIRDTWO_1000 = 30521,\n"
	"TWO_L1_SEMIGD_1000 = 30611,\n"
	"TWO_L1_SEMIGD_SH = 30612,\n"
	"TWO_L2_SEMIGD_1000 = 30621,\n"
	"TWO_L2_SEMIGD_SH = 30622,\n"
	"//for two-variable linear constraints\n"
	"BIAS_L1_RD_1000 = 40111,\n"
	"BIAS_L1_RD_SH = 40112,\n"
	"BIAS_L2_RD_1000 = 40121,\n"
	"BIAS_L2_RD_SH = 40122,\n"
	"BIAS_L1_SEMIGD_1000 = 40211,\n"
	"BIAS_L1_SEMIGD_SH = 40212,\n"
	"BIAS_L2_SEMIGD_1000 = 40221,\n"
	"BIAS_L2_SEMIGD_SH = 40222,\n"
	"BIAS_L1_SEMIGD_RAND_1000 = 40311,\n"
	"BIAS_L1_SEMIGD_RAND_SH = 40312,\n"
	"BIAS_L2_SEMIGD_RAND_1000 = 40321,\n"
	"BIAS_L2_SEMIGD_RAND_SH = 40322,\n"
	"BIAS_L1_SEMIGD_CY_FIRST_1000 = 40411,\n"
	"BIAS_L1_SEMIGD_CY_FIRST_SH = 40412,\n"
	"BIAS_L2_SEMIGD_CY_FIRST_1000 = 40421,\n"
	"BIAS_L2_SEMIGD_CY_FIRST_SH = 40422,\n"
	"BIAS_L1_SEMIGD_RD_FIRST_1000 = 40511,\n"
	"BIAS_L1_SEMIGD_RD_FIRST_SH = 40512,\n"
	"BIAS_L2_SEMIGD_RD_FIRST_1000 = 40521,\n"
	"BIAS_L2_SEMIGD_RD_FIRST_SH = 40522,\n"
	"BIAS_L1_CY_1000 = 40611,\n"
	"BIAS_L1_CY_SH = 40612,\n"
	"BIAS_L2_CY_1000 = 40621,\n"
	"BIAS_L2_CY_SH = 40622,\n"
	"BIAS_L1_SEMIGD_CY_DUALOBJ_1000 = 40711,\n"
	"BIAS_L1_SEMIGD_CY_DUALOBJ_SH = 40712,\n"
	"BIAS_L2_SEMIGD_CY_DUALOBJ_1000 = 40721,\n"
	"BIAS_L2_SEMIGD_CY_DUALOBJ_SH = 40722,\n"
	"BIAS_L1_SEMIGD_RD_DUALOBJ_1000 = 40811,\n"
	"BIAS_L1_SEMIGD_RD_DUALOBJ_SH = 40812,\n"
	"BIAS_L2_SEMIGD_RD_DUALOBJ_1000 = 40821,\n"
	"BIAS_L2_SEMIGD_RD_DUALOBJ_SH = 40822,\n"
	"//for one-class svm\n"
	"ONECLASS_L1_RD_1000 = 50111,\n"
	"ONECLASS_L1_RD_SH = 50112,\n"
	"ONECLASS_L1_SEMIGD_1000 = 50211,\n"
	"ONECLASS_L1_SEMIGD_SH = 50212,\n"
	"ONECLASS_L1_FIRST_1000 = 50311,\n"
	"ONECLASS_L1_SECOND_1000 = 50411,\n"
	"ONECLASS_L1_SEMIGD_RAND_1000 = 50511,\n"
	"ONECLASS_L1_SEMIGD_RAND_SH = 50512,\n"
	"ONECLASS_L1_SEMIGD_CY_FIRST_1000 = 50611,\n"
	"ONECLASS_L1_SEMIGD_CY_FIRST_SH = 50612,\n"
	"ONECLASS_L1_SEMIGD_RD_FIRST_1000 = 50711,\n"
	"ONECLASS_L1_SEMIGD_RD_FIRST_SH = 50712,\n"
	"ONECLASS_L1_CY_1000 = 50811,\n"
	"ONECLASS_L1_CY_SH = 50812,\n"
	"ONECLASS_L1_SEMIGD_CY_DUALOBJ_1000 = 50911,\n"
	"ONECLASS_L1_SEMIGD_CY_DUALOBJ_SH = 50912,\n"
	"ONECLASS_L1_SEMIGD_RD_DUALOBJ_1000 = 51011,\n"
	"ONECLASS_L1_SEMIGD_RD_DUALOBJ_SH = 51012,\n"
	"ONECLASS_L1_SEMIGD_BATCH_1000 = 51111,\n"
	"ONECLASS_L1_SEMIGD_CONV_1000 = 51211,\n"
	"ONECLASS_L1_SEMIGD_SORT_1000 = 51211,\n"
	"//for svdd\n"
	"SVDD_L1_RD_1000 = 60111,\n"
	"SVDD_L1_RD_SH = 60112,\n"
	"SVDD_L1_SEMIGD_1000 = 60211,\n"
	"SVDD_L1_SEMIGD_SH = 60212,\n"
	"SVDD_L1_FIRST_1000 = 60311,\n"
	"SVDD_L1_SECOND_1000 = 60411,\n"
	"SVDD_L1_SEMIGD_RAND_1000 = 60511,\n"
	"SVDD_L1_SEMIGD_RAND_SH = 60512,\n"
	"SVDD_L1_SEMIGD_CY_FIRST_1000 = 60611,\n"
	"SVDD_L1_SEMIGD_CY_FIRST_SH = 60612,\n"
	"SVDD_L1_SEMIGD_RD_FIRST_1000 = 60711,\n"
	"SVDD_L1_SEMIGD_RD_FIRST_SH = 60712,\n"
	"SVDD_L1_CY_1000 = 60811,\n"
	"SVDD_L1_CY_SH = 60812,\n"
	"SVDD_L1_SEMIGD_CY_DUALOBJ_1000 = 60911,\n"
	"SVDD_L1_SEMIGD_CY_DUALOBJ_SH = 60912,\n"
	"SVDD_L1_SEMIGD_RD_DUALOBJ_1000 = 61011,\n"
	"SVDD_L1_SEMIGD_RD_DUALOBJ_SH = 61012,\n"
	"SVDD_L1_SEMIGD_BATCH_1000 = 61111,\n"
	"SVDD_L1_SEMIGD_CONV_1000 = 61211,\n"
	"SVDD_L1_SEMIGD_SORT_1000 = 61311,\n"
	);
	exit(1);
}

void exit_input_error(int line_num)
{
	fprintf(stderr,"Wrong input format at line %d\n", line_num);
	exit(1);
}

static void ltrim(char* s)
{
	int slen = (int) strlen(s);
	int trim_idx = 0;
	bool end=false;
	while(trim_idx < slen && !end)
	{
		switch(s[trim_idx])
		{
			case ' ':  case '\n':
			case '\r': case '\t':
				++trim_idx;
				break;
			default:
				end=true;
		}
	}
	for(int i=0; i+trim_idx<slen+1; i++) //move \0 too
	{
		s[i] = s[i+trim_idx];
	}
}

static void rtrim(char* s)
{
	int slen = (int) strlen(s);
	int trim_idx = slen-1;
	bool end=false;
	while(trim_idx >= 0 && !end)
	{
		switch(s[trim_idx])
		{
			case ' ':  case '\n':
			case '\r': case '\t':
				--trim_idx;
				break;
			default:
				if(trim_idx+1 < slen)
					s[trim_idx+1]='\0';
				end=true;
		}
	}
}

static int trim(char* s)
{
	rtrim(s);
	ltrim(s);
	return (int) strlen(s);
}

static char *line = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
	int len;

	if(fgets(line,max_line_len,input) == NULL)
		return NULL;

	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
			break;
	}
	return line;
}

class GridItem
{
public:
	struct parameter param;
	char *model_file_name;
};

void parse_command_line(int argc, char **argv, char *input_file_name);
void parse_stdin(char *input_file_name, GridItem* grid_item, char *param_str);
void read_problem(const char *filename);
void do_cross_validation(struct parameter* param);
void do_find_parameter_C(struct parameter* param);

struct feature_node *x_space;
struct problem prob;
struct model* model_;
int flag_cross_validation;
int flag_find_C;
int flag_C_specified;
int flag_solver_specified;
int nr_fold;
double bias;

int main(int argc, char **argv)
{
	int GRID_MAX=1024;
	char input_file_name[1024];
	const char *error_msg;
	GridItem grid_items[GRID_MAX];

	parse_command_line(argc, argv, input_file_name);
	char param_str[1024];
	int param_num = 0;
	for(param_num = 0; fgets(param_str, sizeof(param_str), stdin); param_num++)
	{
		grid_items[param_num].model_file_name = Malloc(char, 1024);
		parse_stdin(input_file_name, &(grid_items[param_num]), param_str);
	}
	read_problem(input_file_name);

	#pragma omp parallel for
	for(int i = 0; i < param_num; i++)
	{
		struct parameter param = grid_items[i].param;
		error_msg = check_parameter(&prob,&param);

		if(error_msg)
		{
			fprintf(stderr,"ERROR: %s\n",error_msg);
			exit(1);
		}

		if (flag_find_C)
		{
			do_find_parameter_C(&param);
		}
		else if(flag_cross_validation)
		{
			do_cross_validation(&param);
		}
		else
		{
			model_=train(&prob, &param);
			if(save_model(grid_items[i].model_file_name, model_))
			{
				fprintf(stderr,"can't save model to file %s\n",grid_items[i].model_file_name);
				exit(1);
			}
			free_and_destroy_model(&model_);
		}
		fclose(param.log_fp);
		destroy_param(&param);
	}
	free(prob.y);
	free(prob.x);
	free(x_space);
	free(line);

	return 0;
}

void do_find_parameter_C(struct parameter* param)
{
	double start_C, best_C, best_rate;
	double max_C = 1024;
	if (flag_C_specified)
		start_C = param->C;
	else
		start_C = -1.0;
	printf("Doing parameter search with %d-fold cross validation.\n", nr_fold);
	find_parameter_C(&prob, param, nr_fold, start_C, max_C, &best_C, &best_rate);
	printf("Best C = %g  CV accuracy = %g%%\n", best_C, 100.0*best_rate);
}

void do_cross_validation(struct parameter* param)
{
	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	double *target = Malloc(double, prob.l);

	cross_validation(&prob,param,nr_fold,target);
	if(param->solver_type == L2R_L2LOSS_SVR ||
	   param->solver_type == L2R_L1LOSS_SVR_DUAL ||
	   param->solver_type == L2R_L2LOSS_SVR_DUAL)
	{
		for(i=0;i<prob.l;i++)
		{
			double y = prob.y[i];
			double v = target[i];
			total_error += (v-y)*(v-y);
			sumv += v;
			sumy += y;
			sumvv += v*v;
			sumyy += y*y;
			sumvy += v*y;
		}
		printf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
		printf("Cross Validation Squared correlation coefficient = %g\n",
				((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
				((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
			  );
	}
	else
	{
		for(i=0;i<prob.l;i++)
			if(target[i] == prob.y[i])
				++total_correct;
		printf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/prob.l);
	}

	free(target);
}

void parse_command_line(int argc, char **argv, char *input_file_name)
{
	// determine filenames
	if(argc != 2)
		exit_with_help();

	strcpy(input_file_name, argv[1]);
}
void parse_stdin(char *input_file_name, GridItem* grid_item, char *param_str)
{
	// default values
	grid_item->param.solver_type = ONE_L2_CY_SH;
	grid_item->param.C = 1;
	grid_item->param.eps = INF; // see setting below
	grid_item->param.p = 0.1;
	grid_item->param.r = 1;
	grid_item->param.max_iter = 1000;
	grid_item->param.timeout = 0;
	grid_item->param.max_cdstep = -1;
	grid_item->param.max_nr_n_ops = -1;
	grid_item->param.opt_val = -INF;
	grid_item->param.nu = 0.1;
	grid_item->param.nr_weight = 0;
	grid_item->param.weight_label = NULL;
	grid_item->param.weight = NULL;
	grid_item->param.init_sol = NULL;
	grid_item->param._resume = NULL;
	flag_cross_validation = 0;
	flag_C_specified = 0;
	flag_solver_specified = 0;
	flag_find_C = 0;
	bias = -1;

	char * key_token = strtok(param_str, " ");
	char * val_token;
	// parse options
	while(key_token != NULL)
	{
		if(key_token[0] != '-')
		{
			break;
		}
		switch(key_token[1])
		{
			case 'm':
				val_token = strtok(NULL, " ");
				grid_item->param.max_iter = atoi(val_token);
				break;

			case 'o':
				val_token = strtok(NULL, " ");
				grid_item->param.opt_val = atof(val_token);
				break;

			case 't':
				val_token = strtok(NULL, " ");
				grid_item->param.timeout = atoi(val_token);
				break;

			case 'u':
				val_token = strtok(NULL, " ");
				grid_item->param.max_cdstep = atoi(val_token);
				break;

			case 'S':
				val_token = strtok(NULL, " ");
				grid_item->param.max_nr_n_ops = atof(val_token);
				break;

			case 'r':
				val_token = strtok(NULL, " ");
				grid_item->param.r = atof(val_token);
				break;

			case 'n':
				val_token = strtok(NULL, " ");
				grid_item->param.nu = atof(val_token);
				break;

			case 's':
				val_token = strtok(NULL, " ");
				grid_item->param.solver_type = atoi(val_token);
				flag_solver_specified = 1;
				break;

			case 'c':
				val_token = strtok(NULL, " ");
				grid_item->param.C = atof(val_token);
				flag_C_specified = 1;
				break;

			case 'p':
				val_token = strtok(NULL, " ");
				grid_item->param.p = atof(val_token);
				break;

			case 'e':
				val_token = strtok(NULL, " ");
				grid_item->param.eps = atof(val_token);
				break;

			case 'B':
				val_token = strtok(NULL, " ");
				bias = atof(val_token);
				break;

			case 'w':
				val_token = strtok(NULL, " ");
				++(grid_item->param.nr_weight);
				grid_item->param.weight_label = (int *) realloc(grid_item->param.weight_label,sizeof(int)*grid_item->param.nr_weight);
				grid_item->param.weight = (double *) realloc(grid_item->param.weight,sizeof(double)*grid_item->param.nr_weight);
				grid_item->param.weight_label[grid_item->param.nr_weight-1] = atoi(&key_token[2]);
				grid_item->param.weight[grid_item->param.nr_weight-1] = atof(val_token);
				break;

			case 'v':
				val_token = strtok(NULL, " ");
				flag_cross_validation = 1;
				nr_fold = atoi(val_token);
				if(nr_fold < 2)
				{
					fprintf(stderr,"n-fold cross validation: n must >= 2\n");
					exit_with_help();
				}
				break;

			case 'q':
				set_print_string_function(print_null);
				break;

			case 'C':
				flag_find_C = 1;
				break;

			default:
				fprintf(stderr,"unknown option: -%c\n", key_token[1]);
				exit_with_help();
				break;
		}
		key_token = strtok(NULL, " ");
	}
	char log_file_name[1024]={'\0'};
	if(key_token != NULL && trim(key_token) != 0)
	{
		strcpy(log_file_name, key_token);
	}

	key_token = strtok(NULL, " ");
	if(key_token != NULL && trim(key_token) != 0)
	{
		if((grid_item->param._resume = load_resume(key_token)) == NULL)
		{
			fprintf(stderr,"can't open resume file %s\n",key_token);
			exit(1);
		}
		grid_item->param.log_fp = fopen(log_file_name, "a");
	}
	else
	{
		if(strlen(log_file_name) != 0)
			grid_item->param.log_fp = fopen(log_file_name, "w");
		else
			grid_item->param.log_fp = stdout;
	}

	key_token = strtok(NULL, " ");
	if(key_token != NULL && trim(key_token) != 0)
	{
		strcpy(grid_item->model_file_name, key_token);
	}
	else
	{
		char *p = strrchr(input_file_name,'/');
		if(p==NULL)
			p = input_file_name;
		else
			++p;
		sprintf(grid_item->model_file_name,"%s.model",p);
	}

	// default solver for parameter selection is L2R_L2LOSS_SVC
	if(flag_find_C)
	{
		if(!flag_cross_validation)
			nr_fold = 5;
		if(!flag_solver_specified)
		{
			fprintf(stderr, "Solver not specified. Using -s 2\n");
			grid_item->param.solver_type = L2R_L2LOSS_SVC;
		}
		else if(grid_item->param.solver_type != L2R_LR && grid_item->param.solver_type != L2R_L2LOSS_SVC)
		{
			fprintf(stderr, "Warm-start parameter search only available for -s 0 and -s 2\n");
			exit_with_help();
		}
	}
	if(grid_item->param.eps == INF)
	{
		switch(grid_item->param.solver_type)
		{
			case L2R_LR:
			case L2R_L2LOSS_SVC:
				grid_item->param.eps = 0.01;
				break;
			case L2R_L2LOSS_SVR:
				grid_item->param.eps = 0.001;
				break;
			case OLD_ONE_L2_CY_SH:
			case OLD_ONE_L1_CY_SH:
			case MCSVM_CS:
			case L2R_LR_DUAL:
				grid_item->param.eps = 0.1;
				break;
			case L1R_L2LOSS_SVC:
			case L1R_LR:
				grid_item->param.eps = 0.01;
				break;
			case L2R_L1LOSS_SVR_DUAL:
			case L2R_L2LOSS_SVR_DUAL:
			//for L1L2
			case ONE_L1_CY_1000:
			case ONE_L1_CY_SH:
			case ONE_L2_CY_1000:
			case ONE_L2_CY_SH:
			case ONE_L1_RD_1000:
			case ONE_L2_RD_1000:
			case ONE_L1_RD_SH:
			case ONE_L2_RD_SH:
			case ONE_L1_SEMIGD_1000:
			case ONE_L1_SEMIGD_SH:
			case ONE_L2_SEMIGD_1000:
			case ONE_L2_SEMIGD_SH:
			case ONE_L1_SEMIGD_RAND_1000:
			case ONE_L1_SEMIGD_RAND_SH:
			case ONE_L2_SEMIGD_RAND_1000:
			case ONE_L2_SEMIGD_RAND_SH:
			case ONE_L1_SEMIGD_DUALOBJ_1000:
			case ONE_L1_SEMIGD_DUALOBJ_SH:
			case ONE_L2_SEMIGD_DUALOBJ_1000:
			case ONE_L2_SEMIGD_DUALOBJ_SH:
			case ONE_L1_SEMIGD_DUALOBJ_RAND_1000:
			case ONE_L1_SEMIGD_DUALOBJ_RAND_SH:
			case ONE_L2_SEMIGD_DUALOBJ_RAND_1000:
			case ONE_L2_SEMIGD_DUALOBJ_RAND_SH:
			case ONE_L1_SEMIGD_DUALOBJ_YBAL_1000:
			case ONE_L1_SEMIGD_DUALOBJ_YBAL_SH:
			case ONE_L2_SEMIGD_DUALOBJ_YBAL_1000:
			case ONE_L2_SEMIGD_DUALOBJ_YBAL_SH:
			case TWO_L1_CY_1000:
			case TWO_L2_CY_1000:
			case TWO_L1_RD_1000:
			case TWO_L2_RD_1000:
			case TWO_L1_SEMIGD_1000:
			case TWO_L2_SEMIGD_1000:
			case TWO_L1_SEMIGD_SH:
			case TWO_L2_SEMIGD_SH:
			case TWO_L1_SEMICY_1000:
			case TWO_L2_SEMICY_1000:
			case TWO_L1_SEMIRDONE_1000:
			case TWO_L2_SEMIRDONE_1000:
			case TWO_L1_SEMIRDTWO_1000:
			case TWO_L2_SEMIRDTWO_1000:
			case BIAS_L1_RD_1000:
			case BIAS_L2_RD_1000:
			case BIAS_L1_RD_SH:
			case BIAS_L2_RD_SH:
			case BIAS_L1_CY_1000:
			case BIAS_L2_CY_1000:
			case BIAS_L1_CY_SH:
			case BIAS_L2_CY_SH:
			case BIAS_L1_SEMIGD_1000:
			case BIAS_L2_SEMIGD_1000:
			case BIAS_L1_SEMIGD_SH:
			case BIAS_L2_SEMIGD_SH:
			case BIAS_L1_SEMIGD_RAND_1000:
			case BIAS_L2_SEMIGD_RAND_1000:
			case BIAS_L1_SEMIGD_RAND_SH:
			case BIAS_L2_SEMIGD_RAND_SH:
			case BIAS_L1_SEMIGD_CY_FIRST_1000:
			case BIAS_L2_SEMIGD_CY_FIRST_1000:
			case BIAS_L1_SEMIGD_CY_FIRST_SH:
			case BIAS_L2_SEMIGD_CY_FIRST_SH:
			case BIAS_L1_SEMIGD_RD_FIRST_1000:
			case BIAS_L2_SEMIGD_RD_FIRST_1000:
			case BIAS_L1_SEMIGD_RD_FIRST_SH:
			case BIAS_L2_SEMIGD_RD_FIRST_SH:
			case BIAS_L1_SEMIGD_CY_DUALOBJ_1000:
			case BIAS_L2_SEMIGD_CY_DUALOBJ_1000:
			case BIAS_L1_SEMIGD_CY_DUALOBJ_SH:
			case BIAS_L2_SEMIGD_CY_DUALOBJ_SH:
			case BIAS_L1_SEMIGD_RD_DUALOBJ_1000:
			case BIAS_L2_SEMIGD_RD_DUALOBJ_1000:
			case BIAS_L1_SEMIGD_RD_DUALOBJ_SH:
			case BIAS_L2_SEMIGD_RD_DUALOBJ_SH:
			case TWO_L1_RD_SH:
			case TWO_L2_RD_SH:
			case TWO_L1_RD_SH2:
			case TWO_L2_RD_SH2:
			case ONECLASS_L1_RD_1000:
			case ONECLASS_L1_RD_SH:
			case ONECLASS_L1_CY_1000:
			case ONECLASS_L1_CY_SH:
			case ONECLASS_L1_SECOND_1000:
			case ONECLASS_L1_FIRST_1000:
			case ONECLASS_L1_SEMIGD_1000:
			case ONECLASS_L1_SEMIGD_SH:
			case ONECLASS_L1_SEMIGD_RAND_1000:
			case ONECLASS_L1_SEMIGD_RAND_SH:
			case ONECLASS_L1_SEMIGD_CY_FIRST_1000:
			case ONECLASS_L1_SEMIGD_CY_FIRST_SH:
			case ONECLASS_L1_SEMIGD_RD_FIRST_1000:
			case ONECLASS_L1_SEMIGD_RD_FIRST_SH:
			case ONECLASS_L1_SEMIGD_CY_DUALOBJ_1000:
			case ONECLASS_L1_SEMIGD_CY_DUALOBJ_SH:
			case ONECLASS_L1_SEMIGD_RD_DUALOBJ_1000:
			case ONECLASS_L1_SEMIGD_RD_DUALOBJ_SH:
			case ONECLASS_L1_SEMIGD_BATCH_1000:
			case ONECLASS_L1_SEMIGD_CONV_1000:
			case ONECLASS_L1_SEMIGD_SORT_1000:
			case SVDD_L1_RD_1000:
			case SVDD_L1_RD_SH:
			case SVDD_L1_CY_1000:
			case SVDD_L1_CY_SH:
			case SVDD_L1_SECOND_1000:
			case SVDD_L1_FIRST_1000:
			case SVDD_L1_SEMIGD_1000:
			case SVDD_L1_SEMIGD_SH:
			case SVDD_L1_SEMIGD_RAND_1000:
			case SVDD_L1_SEMIGD_RAND_SH:
			case SVDD_L1_SEMIGD_CY_FIRST_1000:
			case SVDD_L1_SEMIGD_CY_FIRST_SH:
			case SVDD_L1_SEMIGD_RD_FIRST_1000:
			case SVDD_L1_SEMIGD_RD_FIRST_SH:
			case SVDD_L1_SEMIGD_CY_DUALOBJ_1000:
			case SVDD_L1_SEMIGD_CY_DUALOBJ_SH:
			case SVDD_L1_SEMIGD_RD_DUALOBJ_1000:
			case SVDD_L1_SEMIGD_RD_DUALOBJ_SH:
			case SVDD_L1_SEMIGD_BATCH_1000:
			case SVDD_L1_SEMIGD_CONV_1000:
			case SVDD_L1_SEMIGD_SORT_1000:
				grid_item->param.eps = 0.01;
				break;
		}
	}
}

// read in a problem (in libsvm format)
void read_problem(const char *filename)
{
	int max_index, inst_max_index, i;
	size_t elements, j;
	FILE *fp = fopen(filename,"r");
	char *endptr;
	char *idx, *val, *label;

	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	prob.l = 0;
	elements = 0;
	max_line_len = 1024;
	line = Malloc(char,max_line_len);
	while(readline(fp)!=NULL)
	{
		char *p = strtok(line," \t"); // label

		// features
		while(1)
		{
			p = strtok(NULL," \t");
			if(p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
				break;
			elements++;
		}
		elements++; // for bias term
		prob.l++;
	}
	rewind(fp);

	prob.bias=bias;

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct feature_node *,prob.l);
	x_space = Malloc(struct feature_node,elements+prob.l);

	max_index = 0;
	j=0;
	for(i=0;i<prob.l;i++)
	{
		inst_max_index = 0; // strtol gives 0 if wrong format
		readline(fp);
		prob.x[i] = &x_space[j];
		label = strtok(line," \t\n");
		if(label == NULL) // empty line
			exit_input_error(i+1);

		prob.y[i] = strtod(label,&endptr);
		if(endptr == label || *endptr != '\0')
			exit_input_error(i+1);

		while(1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");

			if(val == NULL)
				break;

			errno = 0;
			x_space[j].index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || x_space[j].index <= inst_max_index)
				exit_input_error(i+1);
			else
				inst_max_index = x_space[j].index;

			errno = 0;
			x_space[j].value = strtod(val,&endptr);
			if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(i+1);

			++j;
		}

		if(inst_max_index > max_index)
			max_index = inst_max_index;

		if(prob.bias >= 0)
			x_space[j++].value = prob.bias;

		x_space[j++].index = -1;
	}

	if(prob.bias >= 0)
	{
		prob.n=max_index+1;
		for(i=1;i<prob.l;i++)
			(prob.x[i]-2)->index = prob.n;
		x_space[j-2].index = prob.n;
	}
	else
		prob.n=max_index;

	fclose(fp);
}
