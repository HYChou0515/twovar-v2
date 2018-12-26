#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <locale.h>
#include <algorithm>
#include "linear.h"
#include "tron.h"
//for one-class 
#include <vector>
#include <queue>
#include <time.h>
typedef signed char schar;
template <class T> static inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
#ifndef min
template <class T> static inline T min(T x,T y) { return (x<y)?x:y; }
#endif
#ifndef max
template <class T> static inline T max(T x,T y) { return (x>y)?x:y; }
#endif
template <class S, class T> static inline void clone(T*& dst, S* src, int n)
{
	dst = new T[n];
	memcpy((void *)dst,(void *)src,sizeof(T)*n);
}
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define INF HUGE_VAL

static void print_string_stdout(const char *s)
{
	fputs(s,stdout);
	//fflush(stdout);
}
static void print_null(const char *s) {}

static void (*liblinear_print_string) (const char *) = &print_string_stdout;

#if 1
static void info(const char *fmt,...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	(*liblinear_print_string)(buf);
}
#else
static void info(const char *fmt,...) {}
#endif


class sparse_operator
{
public:
	static double nrm2_sq(const feature_node *x)
	{
		double ret = 0;
		while(x->index != -1)
		{
			ret += x->value*x->value;
			x++;
		}
		return (ret);
	}

	static double dot(const double *s, const feature_node *x)
	{
		double ret = 0;
		while(x->index != -1)
		{
			ret += s[x->index-1]*x->value;
			x++;
		}
		return (ret);
	}
	static double feature_dot(const feature_node *x1, const feature_node *x2)
	{
		double ret = 0;
		while(x1->index != -1 && x2->index != -1)
		{
			if(x1->index == x2->index)
			{
				ret += x1->value * x2->value;
				++x1;
				++x2;
			}
			else
			{
				if(x1->index > x2->index)
					++x2;
				else
					++x1;
			}
		}
		return (ret);
	}

	static void axpy(const double a, const feature_node *x, double *y)
	{
		while(x->index != -1)
		{
			y[x->index-1] += a*x->value;
			x++;
		}
	}
};

class l2r_lr_fun: public function
{
public:
	l2r_lr_fun(const problem *prob, double *C);
	~l2r_lr_fun();

	double fun(double *w);
	void grad(double *w, double *g);
	void Hv(double *s, double *Hs);

	int get_nr_variable(void);

private:
	void Xv(double *v, double *Xv);
	void XTv(double *v, double *XTv);

	double *C;
	double *z;
	double *D;
	const problem *prob;
};

l2r_lr_fun::l2r_lr_fun(const problem *prob, double *C)
{
	int l=prob->l;

	this->prob = prob;

	z = new double[l];
	D = new double[l];
	this->C = C;
}

l2r_lr_fun::~l2r_lr_fun()
{
	delete[] z;
	delete[] D;
}


double l2r_lr_fun::fun(double *w)
{
	int i;
	double f=0;
	double *y=prob->y;
	int l=prob->l;
	int w_size=get_nr_variable();

	Xv(w, z);

	for(i=0;i<w_size;i++)
		f += w[i]*w[i];
	f /= 2.0;
	for(i=0;i<l;i++)
	{
		double yz = y[i]*z[i];
		if (yz >= 0)
			f += C[i]*log(1 + exp(-yz));
		else
			f += C[i]*(-yz+log(1 + exp(yz)));
	}

	return(f);
}

void l2r_lr_fun::grad(double *w, double *g)
{
	int i;
	double *y=prob->y;
	int l=prob->l;
	int w_size=get_nr_variable();

	for(i=0;i<l;i++)
	{
		z[i] = 1/(1 + exp(-y[i]*z[i]));
		D[i] = z[i]*(1-z[i]);
		z[i] = C[i]*(z[i]-1)*y[i];
	}
	XTv(z, g);

	for(i=0;i<w_size;i++)
		g[i] = w[i] + g[i];
}

int l2r_lr_fun::get_nr_variable(void)
{
	return prob->n;
}

void l2r_lr_fun::Hv(double *s, double *Hs)
{
	int i;
	int l=prob->l;
	int w_size=get_nr_variable();
	double *wa = new double[l];
	feature_node **x=prob->x;

	for(i=0;i<w_size;i++)
		Hs[i] = 0;
	for(i=0;i<l;i++)
	{
		feature_node * const xi=x[i];
		wa[i] = sparse_operator::dot(s, xi);
		
		wa[i] = C[i]*D[i]*wa[i];

		sparse_operator::axpy(wa[i], xi, Hs);
	}
	for(i=0;i<w_size;i++)
		Hs[i] = s[i] + Hs[i];
	delete[] wa;
}

void l2r_lr_fun::Xv(double *v, double *Xv)
{
	int i;
	int l=prob->l;
	feature_node **x=prob->x;

	for(i=0;i<l;i++)
		Xv[i]=sparse_operator::dot(v, x[i]);
}

void l2r_lr_fun::XTv(double *v, double *XTv)
{
	int i;
	int l=prob->l;
	int w_size=get_nr_variable();
	feature_node **x=prob->x;

	for(i=0;i<w_size;i++)
		XTv[i]=0;
	for(i=0;i<l;i++)
		sparse_operator::axpy(v[i], x[i], XTv);
}

class l2r_l2_svc_fun: public function
{
public:
	l2r_l2_svc_fun(const problem *prob, double *C);
	~l2r_l2_svc_fun();

	double fun(double *w);
	void grad(double *w, double *g);
	void Hv(double *s, double *Hs);

	int get_nr_variable(void);

protected:
	void Xv(double *v, double *Xv);
	void subXTv(double *v, double *XTv);

	double *C;
	double *z;
	double *D;
	int *I;
	int sizeI;
	const problem *prob;
};

l2r_l2_svc_fun::l2r_l2_svc_fun(const problem *prob, double *C)
{
	int l=prob->l;

	this->prob = prob;

	z = new double[l];
	D = new double[l];
	I = new int[l];
	this->C = C;
}

l2r_l2_svc_fun::~l2r_l2_svc_fun()
{
	delete[] z;
	delete[] D;
	delete[] I;
}

double l2r_l2_svc_fun::fun(double *w)
{
	int i;
	double f=0;
	double *y=prob->y;
	int l=prob->l;
	int w_size=get_nr_variable();

	Xv(w, z);

	for(i=0;i<w_size;i++)
		f += w[i]*w[i];
	f /= 2.0;
	for(i=0;i<l;i++)
	{
		z[i] = y[i]*z[i];
		double d = 1-z[i];
		if (d > 0)
			f += C[i]*d*d;
	}

	return(f);
}

void l2r_l2_svc_fun::grad(double *w, double *g)
{
	int i;
	double *y=prob->y;
	int l=prob->l;
	int w_size=get_nr_variable();

	sizeI = 0;
	for (i=0;i<l;i++)
		if (z[i] < 1)
		{
			z[sizeI] = C[i]*y[i]*(z[i]-1);
			I[sizeI] = i;
			sizeI++;
		}
	subXTv(z, g);

	for(i=0;i<w_size;i++)
		g[i] = w[i] + 2*g[i];
}

int l2r_l2_svc_fun::get_nr_variable(void)
{
	return prob->n;
}

void l2r_l2_svc_fun::Hv(double *s, double *Hs)
{
	int i;
	int w_size=get_nr_variable();
	double *wa = new double[sizeI];
	feature_node **x=prob->x;

	for(i=0;i<w_size;i++)
		Hs[i]=0;
	for(i=0;i<sizeI;i++)
	{
		feature_node * const xi=x[I[i]];
		wa[i] = sparse_operator::dot(s, xi);
		
		wa[i] = C[I[i]]*wa[i];

		sparse_operator::axpy(wa[i], xi, Hs);
	}
	for(i=0;i<w_size;i++)
		Hs[i] = s[i] + 2*Hs[i];
	delete[] wa;
}

void l2r_l2_svc_fun::Xv(double *v, double *Xv)
{
	int i;
	int l=prob->l;
	feature_node **x=prob->x;

	for(i=0;i<l;i++)
		Xv[i]=sparse_operator::dot(v, x[i]);
}

void l2r_l2_svc_fun::subXTv(double *v, double *XTv)
{
	int i;
	int w_size=get_nr_variable();
	feature_node **x=prob->x;

	for(i=0;i<w_size;i++)
		XTv[i]=0;
	for(i=0;i<sizeI;i++)
		sparse_operator::axpy(v[i], x[I[i]], XTv);
}

class l2r_l2_svr_fun: public l2r_l2_svc_fun
{
public:
	l2r_l2_svr_fun(const problem *prob, double *C, double p);

	double fun(double *w);
	void grad(double *w, double *g);

private:
	double p;
};

l2r_l2_svr_fun::l2r_l2_svr_fun(const problem *prob, double *C, double p):
	l2r_l2_svc_fun(prob, C)
{
	this->p = p;
}

double l2r_l2_svr_fun::fun(double *w)
{
	int i;
	double f=0;
	double *y=prob->y;
	int l=prob->l;
	int w_size=get_nr_variable();
	double d;

	Xv(w, z);

	for(i=0;i<w_size;i++)
		f += w[i]*w[i];
	f /= 2;
	for(i=0;i<l;i++)
	{
		d = z[i] - y[i];
		if(d < -p)
			f += C[i]*(d+p)*(d+p);
		else if(d > p)
			f += C[i]*(d-p)*(d-p);
	}

	return(f);
}

void l2r_l2_svr_fun::grad(double *w, double *g)
{
	int i;
	double *y=prob->y;
	int l=prob->l;
	int w_size=get_nr_variable();
	double d;

	sizeI = 0;
	for(i=0;i<l;i++)
	{
		d = z[i] - y[i];

		// generate index set I
		if(d < -p)
		{
			z[sizeI] = C[i]*(d+p);
			I[sizeI] = i;
			sizeI++;
		}
		else if(d > p)
		{
			z[sizeI] = C[i]*(d-p);
			I[sizeI] = i;
			sizeI++;
		}

	}
	subXTv(z, g);

	for(i=0;i<w_size;i++)
		g[i] = w[i] + 2*g[i];
}

// A coordinate descent algorithm for 
// multi-class support vector machines by Crammer and Singer
//
//  min_{\alpha}  0.5 \sum_m ||w_m(\alpha)||^2 + \sum_i \sum_m e^m_i alpha^m_i
//    s.t.     \alpha^m_i <= C^m_i \forall m,i , \sum_m \alpha^m_i=0 \forall i
// 
//  where e^m_i = 0 if y_i  = m,
//        e^m_i = 1 if y_i != m,
//  C^m_i = C if m  = y_i, 
//  C^m_i = 0 if m != y_i, 
//  and w_m(\alpha) = \sum_i \alpha^m_i x_i 
//
// Given: 
// x, y, C
// eps is the stopping tolerance
//
// solution will be put in w
//
// See Appendix of LIBLINEAR paper, Fan et al. (2008)

#define GETI(i) ((int) prob->y[i])
// To support weights for instances, use GETI(i) (i)

class Solver_MCSVM_CS
{
	public:
		Solver_MCSVM_CS(const problem *prob, int nr_class, double *C, double eps=0.1, int max_iter=100000);
		~Solver_MCSVM_CS();
		void Solve(double *w);
	private:
		void solve_sub_problem(double A_i, int yi, double C_yi, int active_i, double *alpha_new);
		bool be_shrunk(int i, int m, int yi, double alpha_i, double minG);
		double *B, *C, *G;
		int w_size, l;
		int nr_class;
		int max_iter;
		double eps;
		const problem *prob;
};

Solver_MCSVM_CS::Solver_MCSVM_CS(const problem *prob, int nr_class, double *weighted_C, double eps, int max_iter)
{
	this->w_size = prob->n;
	this->l = prob->l;
	this->nr_class = nr_class;
	this->eps = eps;
	this->max_iter = max_iter;
	this->prob = prob;
	this->B = new double[nr_class];
	this->G = new double[nr_class];
	this->C = weighted_C;
}

Solver_MCSVM_CS::~Solver_MCSVM_CS()
{
	delete[] B;
	delete[] G;
}

int compare_double(const void *a, const void *b)
{
	if(*(double *)a > *(double *)b)
		return -1;
	if(*(double *)a < *(double *)b)
		return 1;
	return 0;
}

void Solver_MCSVM_CS::solve_sub_problem(double A_i, int yi, double C_yi, int active_i, double *alpha_new)
{
	int r;
	double *D;

	clone(D, B, active_i);
	if(yi < active_i)
		D[yi] += A_i*C_yi;
	qsort(D, active_i, sizeof(double), compare_double);

	double beta = D[0] - A_i*C_yi;
	for(r=1;r<active_i && beta<r*D[r];r++)
		beta += D[r];
	beta /= r;

	for(r=0;r<active_i;r++)
	{
		if(r == yi)
			alpha_new[r] = min(C_yi, (beta-B[r])/A_i);
		else
			alpha_new[r] = min((double)0, (beta - B[r])/A_i);
	}
	delete[] D;
}

bool Solver_MCSVM_CS::be_shrunk(int i, int m, int yi, double alpha_i, double minG)
{
	double bound = 0;
	if(m == yi)
		bound = C[GETI(i)];
	if(alpha_i == bound && G[m] < minG)
		return true;
	return false;
}

void Solver_MCSVM_CS::Solve(double *w)
{
	int i, m, s;
	int iter = 0;
	double *alpha =  new double[l*nr_class];
	double *alpha_new = new double[nr_class];
	int *index = new int[l];
	double *QD = new double[l];
	int *d_ind = new int[nr_class];
	double *d_val = new double[nr_class];
	int *alpha_index = new int[nr_class*l];
	int *y_index = new int[l];
	int active_size = l;
	int *active_size_i = new int[l];
	double eps_shrink = max(10.0*eps, 1.0); // stopping tolerance for shrinking
	bool start_from_all = true;

	// Initial alpha can be set here. Note that 
	// sum_m alpha[i*nr_class+m] = 0, for all i=1,...,l-1
	// alpha[i*nr_class+m] <= C[GETI(i)] if prob->y[i] == m
	// alpha[i*nr_class+m] <= 0 if prob->y[i] != m
	// If initial alpha isn't zero, uncomment the for loop below to initialize w
	for(i=0;i<l*nr_class;i++)
		alpha[i] = 0;

	for(i=0;i<w_size*nr_class;i++)
		w[i] = 0;
	for(i=0;i<l;i++)
	{
		for(m=0;m<nr_class;m++)
			alpha_index[i*nr_class+m] = m;
		feature_node *xi = prob->x[i];
		QD[i] = 0;
		while(xi->index != -1)
		{
			double val = xi->value;
			QD[i] += val*val;

			// Uncomment the for loop if initial alpha isn't zero
			// for(m=0; m<nr_class; m++)
			//	w[(xi->index-1)*nr_class+m] += alpha[i*nr_class+m]*val;
			xi++;
		}
		active_size_i[i] = nr_class;
		y_index[i] = (int)prob->y[i];
		index[i] = i;
	}

	while(iter < max_iter)
	{
		double stopping = -INF;
		for(i=0;i<active_size;i++)
		{
			int j = i+rand()%(active_size-i);
			swap(index[i], index[j]);
		}
		for(s=0;s<active_size;s++)
		{
			i = index[s];
			double Ai = QD[i];
			double *alpha_i = &alpha[i*nr_class];
			int *alpha_index_i = &alpha_index[i*nr_class];

			if(Ai > 0)
			{
				for(m=0;m<active_size_i[i];m++)
					G[m] = 1;
				if(y_index[i] < active_size_i[i])
					G[y_index[i]] = 0;

				feature_node *xi = prob->x[i];
				while(xi->index!= -1)
				{
					double *w_i = &w[(xi->index-1)*nr_class];
					for(m=0;m<active_size_i[i];m++)
						G[m] += w_i[alpha_index_i[m]]*(xi->value);
					xi++;
				}

				double minG = INF;
				double maxG = -INF;
				for(m=0;m<active_size_i[i];m++)
				{
					if(alpha_i[alpha_index_i[m]] < 0 && G[m] < minG)
						minG = G[m];
					if(G[m] > maxG)
						maxG = G[m];
				}
				if(y_index[i] < active_size_i[i])
					if(alpha_i[(int) prob->y[i]] < C[GETI(i)] && G[y_index[i]] < minG)
						minG = G[y_index[i]];

				for(m=0;m<active_size_i[i];m++)
				{
					if(be_shrunk(i, m, y_index[i], alpha_i[alpha_index_i[m]], minG))
					{
						active_size_i[i]--;
						while(active_size_i[i]>m)
						{
							if(!be_shrunk(i, active_size_i[i], y_index[i],
											alpha_i[alpha_index_i[active_size_i[i]]], minG))
							{
								swap(alpha_index_i[m], alpha_index_i[active_size_i[i]]);
								swap(G[m], G[active_size_i[i]]);
								if(y_index[i] == active_size_i[i])
									y_index[i] = m;
								else if(y_index[i] == m)
									y_index[i] = active_size_i[i];
								break;
							}
							active_size_i[i]--;
						}
					}
				}

				if(active_size_i[i] <= 1)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}

				if(maxG-minG <= 1e-12)
					continue;
				else
					stopping = max(maxG - minG, stopping);

				for(m=0;m<active_size_i[i];m++)
					B[m] = G[m] - Ai*alpha_i[alpha_index_i[m]] ;

				solve_sub_problem(Ai, y_index[i], C[GETI(i)], active_size_i[i], alpha_new);
				int nz_d = 0;
				for(m=0;m<active_size_i[i];m++)
				{
					double d = alpha_new[m] - alpha_i[alpha_index_i[m]];
					alpha_i[alpha_index_i[m]] = alpha_new[m];
					if(fabs(d) >= 1e-12)
					{
						d_ind[nz_d] = alpha_index_i[m];
						d_val[nz_d] = d;
						nz_d++;
					}
				}

				xi = prob->x[i];
				while(xi->index != -1)
				{
					double *w_i = &w[(xi->index-1)*nr_class];
					for(m=0;m<nz_d;m++)
						w_i[d_ind[m]] += d_val[m]*xi->value;
					xi++;
				}
			}
		}

		iter++;
		if(iter % 10 == 0)
		{
			info(".");
		}

		if(stopping < eps_shrink)
		{
			if(stopping < eps && start_from_all == true)
				break;
			else
			{
				active_size = l;
				for(i=0;i<l;i++)
					active_size_i[i] = nr_class;
				info("*");
				eps_shrink = max(eps_shrink/2, eps);
				start_from_all = true;
			}
		}
		else
			start_from_all = false;
	}

	info("\noptimization finished, #iter = %d\n",iter);
	if (iter >= max_iter)
		info("\nWARNING: reaching max number of iterations\n");

	// calculate objective value
	double v = 0;
	int nSV = 0;
	for(i=0;i<w_size*nr_class;i++)
		v += w[i]*w[i];
	v = 0.5*v;
	for(i=0;i<l*nr_class;i++)
	{
		v += alpha[i];
		if(fabs(alpha[i]) > 0)
			nSV++;
	}
	for(i=0;i<l;i++)
		v -= alpha[i*nr_class+(int)prob->y[i]];
	info("Objective value = %lf\n",v);
	info("nSV = %d\n",nSV);

	delete [] alpha;
	delete [] alpha_new;
	delete [] index;
	delete [] QD;
	delete [] d_ind;
	delete [] d_val;
	delete [] alpha_index;
	delete [] y_index;
	delete [] active_size_i;
}

// A coordinate descent algorithm for 
// L1-loss and L2-loss SVM dual problems
//
//  min_\alpha  0.5(\alpha^T (Q + D)\alpha) - e^T \alpha,
//    s.t.      0 <= \alpha_i <= upper_bound_i,
// 
//  where Qij = yi yj xi^T xj and
//  D is a diagonal matrix 
//
// In L1-SVM case:
// 		upper_bound_i = Cp if y_i = 1
// 		upper_bound_i = Cn if y_i = -1
// 		D_ii = 0
// In L2-SVM case:
// 		upper_bound_i = INF
// 		D_ii = 1/(2*Cp)	if y_i = 1
// 		D_ii = 1/(2*Cn)	if y_i = -1
//
// Given: 
// x, y, Cp, Cn
// eps is the stopping tolerance
//
// solution will be put in w
// 
// See Algorithm 3 of Hsieh et al., ICML 2008

#undef GETI
#define GETI(i) (y[i]+1)
// To support weights for instances, use GETI(i) (i)

class Solver
{
public:
	Solver(int _solver_type);

	//arguments and options
	const problem *prob;
	int solver_type;	
	double *QD;
	double *diag;
	double *upper_bound;
	double ratio_update;
	schar *y;
	double *w;
	double *alpha;
	double eps;
	int *alpha_status;
	int max_set;
	int *index;
	int max_iter;

	// logged variables
	int iter;
	double PGmax_new;
	double PGmin_new;
	double Gmax;
	double Gmin;
	clock_t duration;
	int active_size;
	int success_pair;
	int update_size;
	int sv;
	int bsv;
	int freesv;
	int nonsv;
	int n_exchange;
	double alpha_diff;
	int nr_pos_y;
	int nr_neg_y;

	// local variables
	double last_obj; //the obj calculated in the last log_message
	enum SolverCate 
	{
		ONE_NOBIAS,
		TWO_NOBIAS,
		TWO_BIAS,
		ONECLASS
	};
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	SolverCate category;
	const char* solver_name;
	//other function
	double calculate_obj();
	double calculate_rho();
	void countSVs();	
	void log_message();
	void summary();
	void debug(const char *fmt, ...);

	//onetwo_nobias_update
	void one_random_shrink();
	void one_cyclic_1000();
	void one_cyclic_shrink();
	void one_random_1000();
	void one_semigd_1000();
	void one_semigd_shrink();
	void two_semicyclic_1000();
	void two_semirandom1_1000();
	void two_semirandom2_1000();
	void two_random_shrink2();
	void two_random_shrink();
	void two_cyclic_1000();
	void two_random_1000();
	void two_semigd_1000();
	void two_semigd_shrink();
	//two_bias_update
	void bias_semigd_shrink();
	void bias_semigd_1000();
	void bias_random_shrink();
	void bias_random_1000();
	//oneclass_update
	void oneclass_random_shrink();
	void oneclass_random_1000();
	void oneclass_first_1000();
	void oneclass_second_1000();
	void oneclass_semigd_1000();	
	void oneclass_semigd_shrink();
};
Solver::Solver(int _solver_type)
{
	iter = 0;
	duration = 0;
	PGmax_new = nan("");
	PGmin_new = nan("");
	Gmax = nan("");
	Gmin = nan("");
	last_obj = nan("");
	active_size = -1;
	success_pair = -1;
	n_exchange = -1;
	update_size = -1;
	solver_type = _solver_type;
	alpha_diff = nan("");
	nr_pos_y = -1;
	nr_neg_y = -1;

	switch(solver_type)
	{
		case ONE_L1_CY_SH:// ONE_L1_CY_SH
		case ONE_L2_CY_SH:// ONE_L2_CY_SH
		case ONE_L1_RD_SH:
		case ONE_L2_RD_SH:
		case ONE_L1_CY_1000:
		case ONE_L2_CY_1000:
		case ONE_L1_RD_1000:
		case ONE_L2_RD_1000:
		case ONE_L1_SEMIGD_1000:
		case ONE_L1_SEMIGD_SH:
		case ONE_L2_SEMIGD_1000:
		case ONE_L2_SEMIGD_SH:
			category = ONE_NOBIAS; 
			break;
		case TWO_L1_SEMICY_1000:
		case TWO_L2_SEMICY_1000:
		case TWO_L1_SEMIRDONE_1000:
		case TWO_L2_SEMIRDONE_1000:
		case TWO_L1_SEMIRDTWO_1000:
		case TWO_L2_SEMIRDTWO_1000:
		case TWO_L1_RD_SH:
		case TWO_L2_RD_SH:
		case TWO_L1_RD_SH2:
		case TWO_L2_RD_SH2:
		case TWO_L1_CY_1000:
		case TWO_L2_CY_1000:
		case TWO_L1_RD_1000:
		case TWO_L2_RD_1000:
		case TWO_L1_SEMIGD_1000:
		case TWO_L2_SEMIGD_1000:
		case TWO_L1_SEMIGD_SH:
		case TWO_L2_SEMIGD_SH:
			category = TWO_NOBIAS; 
			break;
		case BIAS_L1_RD_SH:
		case BIAS_L2_RD_SH:
		case BIAS_L1_RD_1000:
		case BIAS_L2_RD_1000:
		case BIAS_L1_SEMIGD_1000:
		case BIAS_L2_SEMIGD_1000:
		case BIAS_L1_SEMIGD_SH:
		case BIAS_L2_SEMIGD_SH:
			category = TWO_BIAS; 
			break;
		case ONECLASS_L1_RD_SH:
		case ONECLASS_L2_RD_SH:
		case ONECLASS_L1_RD_1000:
		case ONECLASS_L2_RD_1000:
		case ONECLASS_L1_FIRST_1000:
		case ONECLASS_L2_FIRST_1000:
		case ONECLASS_L1_SECOND_1000:
		case ONECLASS_L2_SECOND_1000:
		case ONECLASS_L1_SEMIGD_1000:
		case ONECLASS_L2_SEMIGD_1000:
		case ONECLASS_L1_SEMIGD_SH:
		case ONECLASS_L2_SEMIGD_SH:
			category = ONECLASS;
			break;
	}
#define SAVE_NAME(p) case(p): solver_name = #p; break;
	switch(solver_type)
	{
		SAVE_NAME(ONE_L1_RD_SH);
		SAVE_NAME(ONE_L2_RD_SH);
		SAVE_NAME(ONE_L1_CY_1000);
		SAVE_NAME(ONE_L2_CY_1000);
		SAVE_NAME(ONE_L1_CY_SH);
		SAVE_NAME(ONE_L2_CY_SH);
		SAVE_NAME(ONE_L1_RD_1000);
		SAVE_NAME(ONE_L2_RD_1000);
		SAVE_NAME(ONE_L1_SEMIGD_1000);
		SAVE_NAME(ONE_L1_SEMIGD_SH);
		SAVE_NAME(ONE_L2_SEMIGD_1000);
		SAVE_NAME(ONE_L2_SEMIGD_SH);
		SAVE_NAME(TWO_L1_SEMICY_1000);
		SAVE_NAME(TWO_L2_SEMICY_1000);
		SAVE_NAME(TWO_L1_SEMIRDONE_1000);
		SAVE_NAME(TWO_L2_SEMIRDONE_1000);
		SAVE_NAME(TWO_L1_SEMIRDTWO_1000);
		SAVE_NAME(TWO_L2_SEMIRDTWO_1000);
		SAVE_NAME(TWO_L1_RD_SH);
		SAVE_NAME(TWO_L2_RD_SH);
		SAVE_NAME(TWO_L1_RD_SH2);
		SAVE_NAME(TWO_L2_RD_SH2);
		SAVE_NAME(TWO_L1_CY_1000);
		SAVE_NAME(TWO_L2_CY_1000);
		SAVE_NAME(TWO_L1_RD_1000);
		SAVE_NAME(TWO_L2_RD_1000);
		SAVE_NAME(TWO_L1_SEMIGD_1000);
		SAVE_NAME(TWO_L2_SEMIGD_1000);
		SAVE_NAME(TWO_L1_SEMIGD_SH);
		SAVE_NAME(TWO_L2_SEMIGD_SH);
		SAVE_NAME(BIAS_L1_RD_SH);
		SAVE_NAME(BIAS_L2_RD_SH);
		SAVE_NAME(BIAS_L1_RD_1000);
		SAVE_NAME(BIAS_L2_RD_1000);
		SAVE_NAME(BIAS_L1_SEMIGD_1000);
		SAVE_NAME(BIAS_L2_SEMIGD_1000);
		SAVE_NAME(BIAS_L1_SEMIGD_SH);
		SAVE_NAME(BIAS_L2_SEMIGD_SH);
		SAVE_NAME(ONECLASS_L1_RD_SH);
		SAVE_NAME(ONECLASS_L2_RD_SH);
		SAVE_NAME(ONECLASS_L1_RD_1000);
		SAVE_NAME(ONECLASS_L2_RD_1000);
		SAVE_NAME(ONECLASS_L1_FIRST_1000);
		SAVE_NAME(ONECLASS_L2_FIRST_1000);
		SAVE_NAME(ONECLASS_L1_SECOND_1000);
		SAVE_NAME(ONECLASS_L2_SECOND_1000);
		SAVE_NAME(ONECLASS_L1_SEMIGD_1000);
		SAVE_NAME(ONECLASS_L2_SEMIGD_1000);
		SAVE_NAME(ONECLASS_L1_SEMIGD_SH);
		SAVE_NAME(ONECLASS_L2_SEMIGD_SH);
	}
#undef SAVE_NAME
}
double Solver::calculate_obj()
{
	int i;
	double v = 0;
	switch(category)
	{
		case ONE_NOBIAS:
		case TWO_NOBIAS:
		case TWO_BIAS:
		{
			for(i=0; i<prob->n; i++)
				v += w[i]*w[i];
			for(i=0; i<prob->l; i++)
			{
				v += alpha[i]*(alpha[i]*diag[GETI(i)] - 2);
			}
			v /= 2;
			break;
		}
		case ONECLASS:
		{
			for(i=0; i<prob->n; i++)
				v += w[i] * w[i];
			v /= 2;
			break;
		}
	}
	return v;
}
double Solver::calculate_rho()
{
	double r = 0;
	if(category == ONE_NOBIAS || category == TWO_NOBIAS)
		return 0;
	double sum = 0;
	int n_free = 0;
	double lb = -INF;
	double ub = INF;
	for(int i = 0; i < prob->l; i++)
	{
		double yy = 0;
		if(category == TWO_BIAS)
			yy = y[i];
		double yg =  sparse_operator::dot(w, prob->x[i]) - yy;
		if(alpha_status[i] == FREE)
		{
			n_free++;
			sum += yg;

		}
		else if(alpha_status[i] == UPPER_BOUND)
		{
			if(yg>lb)
				lb = yg;
		}
		else if(alpha_status[i] == LOWER_BOUND)
		{
			if(yg<ub)
				ub = yg;
		}
	}
	if(n_free>0)
		r = sum / n_free;
	else
		r = (lb + ub)/2;
	return r;
}
void Solver:: countSVs()
{
	sv = 0;
	bsv = 0;
	freesv = 0;
	nonsv = 0;
	if(category == ONE_NOBIAS || category == TWO_NOBIAS)
	{
		for(int i = 0; i < prob->l; i++)
		{
			if(alpha[i] >= upper_bound[GETI(i)])
			{
				sv++;
				bsv++;
			}
			else if(alpha[i] <= 0)
			{
				nonsv++;
			}
			else
			{
				sv++;
				freesv++;
			}
		}
	}
	else // TWO_BIAS, ONECLASS
	{
		for(int i = 0; i < prob->l; i++)
		{
			if(alpha_status[i] == FREE)
			{
				sv++;
				freesv++;
			}
			else if(alpha_status[i] == UPPER_BOUND)
			{
				sv++;
				bsv++;
			}
			else if(alpha_status[i] == LOWER_BOUND)
			{
				nonsv++;
			}
		}
	}
}
void Solver::log_message()
{
	
	info("iter %d ", iter);
	info("t %5f ", (double)(duration)/CLOCKS_PER_SEC);
	
	double new_obj = calculate_obj();
	info("obj %.16g ", new_obj);
	info("decr_rate %.3e ", (last_obj-new_obj)/fabs(new_obj));
	last_obj = new_obj;

	int pseudo_updsize;
	//default: see update_size and active_size are the same variable
	if(update_size == -1)
		pseudo_updsize = active_size;
	else
		pseudo_updsize = update_size;
	info("actsize %d ", active_size);
	info("updsize %d ", pseudo_updsize);
	info("sucpair %d ", success_pair);
	info("sucs_rate %.2f%% ", (double)success_pair/pseudo_updsize*100);
	
	countSVs();
	info("nSV %d ", sv);
	info("nBSV %d ", bsv);
	info("nFree %d ", freesv);
	info("nNonSV %d ", nonsv);

	info("PGmax %.16g ", PGmax_new);
	info("PGmin %.16g ", PGmin_new);
	info("PGdiff %.3f ", PGmax_new-PGmin_new);
	
	info("Gmax %.16g ", Gmax);
	info("Gmin %.16g ", Gmin);
	info("Gdiff %.3f ", Gmax-Gmin);

	info("n_exchange %d ", n_exchange);

	info("alpha_diff %.16g ", alpha_diff);

	info("nr_pos_y %d ", nr_pos_y);
	info("nr_neg_y %d ", nr_neg_y);

	info("\n");
}
void Solver::summary()
{
	info("\n");
	info("l = %d n = %d\n", prob->l, prob->n);
	info("eps = %f ratio_update = %f\n", eps, ratio_update);
	info("solver = %s\n", solver_name);
	info("max_iter = %d\n", max_iter);
	info("obj = %.16g rho = %.16g\n", calculate_obj(), calculate_rho());

	countSVs();
	info("nSV = %d : %.2lf%% of l\n", sv, (double)sv/prob->l*100);
	info("nBSV = %d : %.2lf%% of nSV %.2lf%% of l\n", bsv, (double)bsv/sv*100, (double)bsv/prob->l*100);
	info("nFree = %d : %.2lf%% of nSV %.2lf%% of l\n", freesv, (double)freesv/sv*100, (double)freesv/prob->l*100);
	info("nNonSV = %d : %.2lf%% of l\n", nonsv, (double)nonsv/prob->l*100);
}
void Solver::debug(const char *fmt, ...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	fprintf( stderr, "%s", buf);
	fflush(stderr);
}

static int updateAlphaStatus(double alpha, double C)
{
	if(alpha <= 0)
		return 0;
	else if (alpha >= C)
		return 1;
	else 
		return 2;
}

struct mincomp
{
	bool operator()(const feature_node& left, const feature_node& right)
	{
		return left.value > right.value;
	}
};
struct maxcomp
{
	bool operator()(const feature_node& left, const feature_node& right)
	{
		return left.value < right.value;
	}
};

void Solver::one_random_shrink()
{
	int l = max_set;
	int i, s, si;
	double C, d, G;
	// PG: projected gradient, for shrinking and stopping
	double PG;
	double PGmax_old = INF;
	double PGmin_old = -INF;
	clock_t start;
	active_size = l;

	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		PGmax_new = -INF;
		PGmin_new = INF;
		for (s=0; s<active_size; s++)
		{
			si = rand()%active_size;
			i = index[si];
			const schar yi = y[i];
			feature_node * const xi = prob->x[i];
			G = yi*sparse_operator::dot(w, xi)-1;
			C = upper_bound[GETI(i)];
			G += alpha[i]*diag[GETI(i)];

			PG = 0;
			if (alpha[i] == 0)
			{
				if (G > PGmax_old)
				{
					active_size--;
					swap(index[si], index[active_size]);
					continue;
				}
				else if (G < 0)
					PG = G;
			}
			else if (alpha[i] == C)
			{
				if (G < PGmin_old)
				{
					active_size--;
					swap(index[si], index[active_size]);
					continue;
				}
				else if (G > 0)
					PG = G;
			}
			else
				PG = G;

			PGmax_new = max(PGmax_new, PG);
			PGmin_new = min(PGmin_new, PG);

			double alpha_old = alpha[i];
			alpha[i] = min(max(alpha[i] - G/QD[i], 0.0), C);
			d = (alpha[i] - alpha_old)*yi;
			if(fabs(d) > 1.0e-16)
			{
				success_pair++;
				sparse_operator::axpy(d, xi, w);
			}
		}
		iter++;
		duration += clock() - start;
		log_message();
		if(PGmax_new - PGmin_new <= eps)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				PGmax_old = INF;
				PGmin_old = -INF;
				continue;
			}
		}
		PGmax_old = PGmax_new;
		PGmin_old = PGmin_new;
		if (PGmax_old <= 0)
			PGmax_old = INF;
		if (PGmin_old >= 0)
			PGmin_old = -INF;
	}
	summary();
}

void Solver::one_cyclic_1000()
{
	int l = max_set;
	int i, s;
	double C, d, G;
	clock_t start;
	active_size = l;
	
	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		for (i=0; i<active_size; i++)
		{
			int j = i+rand()%(active_size-i);
			swap(index[i], index[j]);
		}
		for (s=0; s<active_size; s++)
		{
			i = index[s];
			const schar yi = y[i];
			feature_node * const xi = prob->x[i];

			G = yi*sparse_operator::dot(w, xi)-1;

			C = upper_bound[GETI(i)];
			G += alpha[i]*diag[GETI(i)];

			double alpha_old = alpha[i];
			alpha[i] = min(max(alpha[i] - G/QD[i], 0.0), C);
			d = (alpha[i] - alpha_old)*yi;
			if(fabs(d) > 1.0e-16)
			{
				success_pair++;
				sparse_operator::axpy(d, xi, w);
			}
		}
		iter++;
		duration += clock() -start;
		log_message();
	}
	summary();
}

void Solver::one_semigd_shrink()
{
	info("THIS IS one semigd shrink\n");
	//TODO: implement
}

void Solver::one_semigd_1000()
{
	int l = max_set;
	int i, s;
	clock_t start;
	double G_i, C_i;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	double *G = new double[l];
	double PG;
	std::priority_queue<struct feature_node, std::vector<feature_node>, mincomp> Min_order_queue;
	int *workset = new int[l];
	int smgd_size;
	bool *need_update_ind = new bool[l];

	active_size = l;
	while (iter < max_iter)
	{
		for(i=0; i<l; i++)
		{
			need_update_ind[i]=true;
		}
		start = clock();
		success_pair = 0;
		nr_pos_y = 0;
		nr_neg_y = 0;
		for (i=0; i<active_size; i++)
		{
			int j = i+rand()%(active_size-i);
			swap(index[i], index[j]);
		}
		PGmax_new = -INF;
		PGmin_new = INF;
		
		for(i=0; i<active_size; i++)
		{
			int G_index = index[i];
			feature_node * const xi = prob->x[G_index];
			G[i] = y[G_index]*sparse_operator::dot(w, xi) -1 +alpha[G_index]*diag[GETI(G_index)];
		}
		//TODO: rlist as number not ratio
		//update_size = int(active_size*ratio_update);
		update_size = int(ratio_update);
		if(update_size < 1)
			update_size = 1;
		smgd_size = update_size;
		//update_size = l;

		for(i=0; i<active_size; i++)
		{
			PG = 0;	
			struct feature_node comp;
			comp.index = index[i];

			if (alpha_status[comp.index] == LOWER_BOUND)
			{
				if(G[i] < 0)
					PG = G[i];
			}
			else if (alpha_status[comp.index] == UPPER_BOUND)
			{
				if(G[i] > 0)
					PG = G[i];
			}
			else
				PG = G[i];

			PGmax_new = max(PGmax_new, PG);
			PGmin_new = min(PGmin_new, PG);

			comp.value = fabs(PG);
			if((int) Min_order_queue.size() < smgd_size)
			{
				need_update_ind[comp.index] = false;
				Min_order_queue.push(comp);
			}
			else if ((int) Min_order_queue.size() <= 0)
			{
				// nothing
			}
			else if (Min_order_queue.top().value < comp.value)
			{
				need_update_ind[Min_order_queue.top().index] = true;
				Min_order_queue.pop();
				need_update_ind[comp.index] = false;
				Min_order_queue.push(comp);
			}
		}
		smgd_size = (int) Min_order_queue.size();
		update_size = smgd_size;

		for(i=0; i<smgd_size; i++)
		{
			workset[smgd_size-1-i] = Min_order_queue.top().index;
			Min_order_queue.pop();
		}
		//j=0;
		//for(i=0; i<update_size; i++)
		//{
		//	if(need_update_ind[i])
		//	{
		//		workset[smgd_size+j] = i;
		//		++j;
		//	}
		//}
		//std::random_shuffle(workset, workset+update_size);
		for(s=0; s<update_size; s++)
		{
			i = workset[s];
			const schar yi = y[i];
			if(yi == +1)
				++nr_pos_y;
			if(yi == -1)
				++nr_neg_y;
			feature_node const * xi = prob->x[i];

			G_i = yi*sparse_operator::dot(w, xi)-1;
			C_i = upper_bound[GETI(i)];
			G_i += alpha[i]*diag[GETI(i)];

			double alpha_new = min(max(alpha[i] - G_i/QD[i], 0.0), C_i);
			alpha_diff = alpha_new - alpha[i];
			if (fabs(alpha_diff) > 1e-16)
			{
				success_pair++;
				alpha[i] = alpha_new;
				sparse_operator::axpy(alpha_diff*yi, xi, w);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[GETI(i)]);
			}
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
}

void Solver::one_random_1000()
{
	int l = max_set;
	int i, s;
	double C, d, G;
	clock_t start;
	active_size = l;
	
	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		nr_pos_y = 0;
		nr_neg_y = 0;
		for (s=0; s<active_size; s++)
		{
			
			i = index[rand()%active_size];
			const schar yi = y[i];
			if(yi == +1)
				++nr_pos_y;
			if(yi == -1)
				++nr_neg_y;
			feature_node * const xi = prob->x[i];
			G = yi*sparse_operator::dot(w, xi)-1;
			C = upper_bound[GETI(i)];
			G += alpha[i]*diag[GETI(i)];

			double alpha_old = alpha[i];
			alpha[i] = min(max(alpha[i] - G/QD[i], 0.0), C);
			d = (alpha[i] - alpha_old)*yi;
			if(fabs(d) > 1.0e-16)
			{
				success_pair++;
				sparse_operator::axpy(d, xi, w);
			}
		}

		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
}

void Solver::two_semicyclic_1000()
{
	int l = max_set;
	int i, j, si;
	int k=0;
	active_size = l;
	double C_i, C_j, G_i, G_j;
	clock_t start;
	int *pickindex = new int[l];
	
	for (i=0; i<active_size; i++)
	{
		pickindex[i] = index[i];
	}
	for (i=0; i<active_size; i++)
	{
		int j = i+rand()%(active_size-i);
		swap(pickindex[i], pickindex[j]);
	}
	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		for (i=0; i<active_size; i++)
		{
			int j = i+rand()%(active_size-i);
			swap(index[i], index[j]);
		}
		k = k%active_size;
		i = pickindex[k];
		for(si=0; si<active_size; si++)
		{
			j = index[si];
			if(i==j)
				continue;
			const schar yi = y[i];
			const schar yj = y[j];
				
			C_i = upper_bound[GETI(i)];
			C_j = upper_bound[GETI(j)];
			feature_node const *xi = prob->x[i];
			feature_node const *xj = prob->x[j];
			
			G_i = 0;
			G_j = 0;
			double Q_ij = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;
				
			double base = QD[i]*QD[j] - Q_ij*Q_ij;
			double newalpha_i, newalpha_j;
			if(base == 0 && (fabs(-QD[j]*G_i + Q_ij * G_j) < 1.0e-12 || fabs(-QD[i]*G_j + Q_ij * G_i <1.0e-12)))
			{
				double delta = QD[i]*alpha[i] + Q_ij*alpha[j] - G_i;
				if(Q_ij < 0)
				{
					if(QD[i]*C_i <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = 0;
					}
					else if( Q_ij*C_j >= delta )
					{
						newalpha_i = 0;
						newalpha_j = C_j;
					}
					else if(0 <= delta)
					{
						newalpha_i = delta/QD[i];
						newalpha_j = 0;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
				}
				else
				{
					if( 0 >= delta)
					{
						newalpha_i = 0;
						newalpha_j = 0;
					}
					else if(QD[i]*C_i + Q_ij*C_j <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = C_j;
					}
					else if( Q_ij*C_j <= delta)
					{
						newalpha_i = (delta-Q_ij*C_j)/QD[i];
						newalpha_j = C_j;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
					
				}
			}
			else
			{
				bool use_j = false;
				newalpha_i = min(C_i, max(0.0, alpha[i] + (-QD[j]*G_i + Q_ij * G_j)/base));
				newalpha_j = min(C_j, max(0.0, alpha[j] + (-QD[i]*G_j + Q_ij * G_i)/base));
				double d_i = newalpha_i - alpha[i];
				if(newalpha_i >= C_i)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i<=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;

				}
				else if(newalpha_i <= 0.0)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i>=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;
				}
				else{
					use_j = true;
				}
				if(use_j == true)
					newalpha_i = min(C_i, max(0.0, alpha[i] - (Q_ij*(newalpha_j - alpha[j]) + G_i)/QD[i]));
			}
			bool update = false;
			if(fabs(newalpha_i - alpha[i]) > 1.0e-16)
			{
				sparse_operator::axpy(yi*(newalpha_i - alpha[i]), prob->x[i], w);
				alpha[i] = newalpha_i;
				update = true;
			}
			if(fabs(newalpha_j - alpha[j]) > 1.0e-16)
			{
				sparse_operator::axpy(yj*(newalpha_j - alpha[j]), prob->x[j], w);
				alpha[j] = newalpha_j;
				update = true;
			}
			if(update)
				success_pair += 1;
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
}

void Solver::two_random_shrink2()
{
	int l = max_set;
	int i, j, s;
	active_size = l;
	double C_i, C_j, G_i, G_j;
	double PG_i , PG_j;
	double PGmax_old = INF;
	double PGmin_old = -INF;
	
	clock_t start;
	while (iter < max_iter)
	{
		start = clock();
		PGmax_new = -INF;
		PGmin_new = INF;
		success_pair = 0;
		for(s=0; s<active_size; s++)
		{
			int si = rand()%active_size;
			int sj = rand()%active_size;
			
			while(si==sj)
			{
				sj = rand()%active_size;
			}
			
			i = index[si];
			j = index[sj];

			const schar yi = y[i];
			const schar yj = y[j];
			C_i = upper_bound[GETI(i)];
			C_j = upper_bound[GETI(j)];
			
			feature_node const *xi = prob->x[i];
			feature_node const *xj = prob->x[j];
	
			G_i = 0;
			G_j = 0;
			double Q_ij = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;

			PG_i = 0;
			if(alpha[i] == 0)
			{
				if(G_i > PGmax_old)
				{
					active_size--;
					swap(index[si],index[active_size]);
					continue;
				}
			}
			else if (alpha[i] == C_i)
			{
				if(G_i < PGmin_old)
				{
					active_size--;
					swap(index[si], index[active_size]);
					continue;
				}
			}
			else
			{
				PG_i = G_i;
			}
			
			PGmax_new = max(PGmax_new, PG_i);
			PGmin_new = min(PGmin_new, PG_i);
			PG_j = 0;
			if(alpha[j] == 0)
			{
				if(G_j > PGmax_old)
				{
					active_size--;
					swap(index[sj],index[active_size]);
					continue;
				}
			}
			else if (alpha[j] == C_j)
			{
				if(G_j < PGmin_old)
				{
					active_size--;
					swap(index[sj], index[active_size]);
					continue;
				}
			}
			else
			{
				PG_j = G_j;
			}
			PGmax_new = max(PGmax_new, PG_j);
			PGmin_new = min(PGmin_new, PG_j);
					
			xi = prob->x[i];
			xj = prob->x[j];
			double base = QD[i]*QD[j] - Q_ij*Q_ij;
			bool use_j = false;
			double newalpha_i, newalpha_j;
			newalpha_i = min(C_i, max(0.0, alpha[i] + (-QD[j]*G_i + Q_ij * G_j)/base));
			newalpha_j = min(C_j, max(0.0, alpha[j] + (-QD[i]*G_j + Q_ij * G_i)/base));
			double d_i = newalpha_i - alpha[i];
			if(base == 0 && (fabs(-QD[j]*G_i + Q_ij * G_j) < 1.0e-12 || fabs(-QD[i]*G_j + Q_ij * G_i <1.0e-12)))
			{
				double delta = QD[i]*alpha[i] + Q_ij*alpha[j] - G_i;
				if(Q_ij < 0)
				{
					if(QD[i]*C_i <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = 0;
					}
					else if( Q_ij*C_j >= delta )
					{
						newalpha_i = 0;
						newalpha_j = C_j;
					}
					else if(0 <= delta)
					{
						newalpha_i = delta/QD[i];
						newalpha_j = 0;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
				}
				else
				{
					if( 0 >= delta)
					{
						newalpha_i = 0;
						newalpha_j = 0;
					}
					else if(QD[i]*C_i + Q_ij*C_j <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = C_j;
					}
					else if( Q_ij*C_j <= delta)
					{
						newalpha_i = (delta-Q_ij*C_j)/QD[i];
						newalpha_j = C_j;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
					
				}
			}
			else
			{
				if(newalpha_i >= C_i)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i<=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;

				}
				else if(newalpha_i <= 0.0)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i>=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;
				}
				else
					use_j = true;
				
				if(use_j == true)
					newalpha_i = min(C_i, max(0.0, alpha[i] - (Q_ij*(newalpha_j - alpha[j]) + G_i)/QD[i]));
			}
			bool update = false;
			if(fabs(newalpha_i - alpha[i]) > 1.0e-16)
			{
				sparse_operator::axpy(yi*(newalpha_i - alpha[i]), xi, w);
				alpha[i] = newalpha_i;
				update = true;
			}
			if(fabs(newalpha_j - alpha[j]) > 1.0e-16)
			{
				sparse_operator::axpy(yj*(newalpha_j - alpha[j]), xj, w);
				alpha[j] = newalpha_j;
				update = true;
			}
			if(update)
				success_pair++;
		}
		iter++;
		duration += clock() - start;
		log_message();
		if(PGmax_new - PGmin_new <= eps)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				PGmax_old = INF;
				PGmin_old = -INF;
				continue;
			}
		}
		PGmax_old = PGmax_new;
		PGmin_old = PGmin_new;
		if(PGmax_old <= 0)
			PGmax_old = INF;
		if(PGmin_old >= 0)
			PGmin_old = -INF;
	}
	summary();
}

void Solver::two_random_shrink()
{
	int l = max_set;
	int i, j, s;
	active_size = l;
	double C_i, C_j, G_i, G_j;
	double PG_i , PG_j;
	double PGmax_old = INF;
	double PGmin_old = -INF;
	
	clock_t start;
	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		PGmax_new = -INF;
		PGmin_new = INF;
		for(s=0; s<active_size; s++)
		{
			int si = rand()%active_size;
			int sj = rand()%active_size;
			
			while(si==sj)
			{
				sj = rand()%active_size;
			}
			
			i = index[si];
			j = index[sj];

			const schar yi = y[i];
			const schar yj = y[j];
			C_i = upper_bound[GETI(i)];
			C_j = upper_bound[GETI(j)];
			
			feature_node const *xi = prob->x[i];
			feature_node const *xj = prob->x[j];
	
			G_i = 0;
			G_j = 0;
			double Q_ij = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;

			PG_i = 0;
			if(alpha[i] == 0)
			{
				if(G_i > PGmax_old)
				{
					active_size--;
					swap(index[si],index[active_size]);
					continue;
				}
				else if(G_i < 0)
					PG_i = G_i;
			}
			else if (alpha[i] == C_i)
			{
				if(G_i < PGmin_old)
				{
					active_size--;
					swap(index[si], index[active_size]);
					continue;
				}
				else if(G_i > 0)
					PG_i = G_i;
			}
			else
				PG_i = G_i;
			
			PGmax_new = max(PGmax_new, PG_i);
			PGmin_new = min(PGmin_new, PG_i);
			PG_j = 0;
			if(alpha[j] == 0)
			{
				if(G_j > PGmax_old)
				{
					active_size--;
					swap(index[sj],index[active_size]);
					continue;
				}
				else if(G_j < 0)
					PG_j = G_j;
			}
			else if (alpha[j] == C_j)
			{
				if(G_j < PGmin_old)
				{
					active_size--;
					swap(index[sj], index[active_size]);
					continue;
				}
				else if(G_j > 0)
					PG_j = G_j;
			}
			else
				PG_j = G_j;

			PGmax_new = max(PGmax_new, PG_j);
			PGmin_new = min(PGmin_new, PG_j);
					
			xi = prob->x[i];
			xj = prob->x[j];
			double base = QD[i]*QD[j] - Q_ij*Q_ij;
			bool use_j = false;
			double newalpha_i, newalpha_j;
			newalpha_i = min(C_i, max(0.0, alpha[i] + (-QD[j]*G_i + Q_ij * G_j)/base));
			newalpha_j = min(C_j, max(0.0, alpha[j] + (-QD[i]*G_j + Q_ij * G_i)/base));
			double d_i = newalpha_i - alpha[i];
			if(base == 0 && (fabs(-QD[j]*G_i + Q_ij * G_j) < 1.0e-12 || fabs(-QD[i]*G_j + Q_ij * G_i <1.0e-12)))
			{
				double delta = QD[i]*alpha[i] + Q_ij*alpha[j] - G_i;
				if(Q_ij < 0)
				{
					if(QD[i]*C_i <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = 0;
					}
					else if( Q_ij*C_j >= delta )
					{
						newalpha_i = 0;
						newalpha_j = C_j;
					}
					else if(0 <= delta)
					{
						newalpha_i = delta/QD[i];
						newalpha_j = 0;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
				}
				else
				{
					if( 0 >= delta)
					{
						newalpha_i = 0;
						newalpha_j = 0;
					}
					else if(QD[i]*C_i + Q_ij*C_j <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = C_j;
					}
					else if( Q_ij*C_j <= delta)
					{
						newalpha_i = (delta-Q_ij*C_j)/QD[i];
						newalpha_j = C_j;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
					
				}
			}
			else
			{
				if(newalpha_i >= C_i)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i<=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;

				}
				else if(newalpha_i <= 0.0)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i>=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;
				}
				else
					use_j = true;
				
				if(use_j == true)
					newalpha_i = min(C_i, max(0.0, alpha[i] - (Q_ij*(newalpha_j - alpha[j]) + G_i)/QD[i]));
			}
			bool update = false;
			if(fabs(newalpha_i - alpha[i]) > 1.0e-16)
			{
				sparse_operator::axpy(yi*(newalpha_i - alpha[i]), xi, w);
				alpha[i] = newalpha_i;
				update = true;
			}
			if(fabs(newalpha_j - alpha[j]) > 1.0e-16)
			{
				sparse_operator::axpy(yj*(newalpha_j - alpha[j]), xj, w);
				alpha[j] = newalpha_j;
				update = true;
			}
			if(update)
				success_pair++;
		}
		iter++;
		duration += clock() - start;
		log_message();
		if(PGmax_new - PGmin_new <= eps)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				PGmax_old = INF;
				PGmin_old = -INF;
				continue;
			}
		}
		PGmax_old = PGmax_new;
		PGmin_old = PGmin_new;
		if(PGmax_old <= 0)
			PGmax_old = INF;
		if(PGmin_old >= 0)
			PGmin_old = -INF;

	}
	summary();
}

void Solver::two_cyclic_1000()
{
	int l = max_set;
	int full_l = prob->l;
	int i, j, s, si;
	double C_i, C_j, G_i, G_j;
	active_size = l;
	int *all_index = new int[l*(l-1)/2];
	clock_t start;
	s=0;

	for(i=0;i<l;i++)
	{
		for(j=i+1;j<l;j++)
		{
			all_index[s]=index[i]*full_l+index[j];
			s++;
		}
		
	}
	for(i=0;i<l*(l-1)/2;i++)
	{
		int j = i+rand()%(l*(l-1)/2-i);
		swap(all_index[i],all_index[j]);
	}
	s = 0;
	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		for(si=0; si<active_size; si++)
		{
			int tempidx = all_index[s];
			i = tempidx / full_l;
			j = tempidx - i*full_l;
			s++;
			s = s % (l*(l-1)/2);
			
			const schar yi = y[i];
			const schar yj = y[j];
			
			C_i = upper_bound[GETI(i)];
			C_j = upper_bound[GETI(j)];
			
			feature_node const *xi = prob->x[i];
			feature_node const *xj = prob->x[j];
			
			G_i = 0;
			G_j = 0;
			double Q_ij = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;
			
			double base = QD[i]*QD[j] - Q_ij*Q_ij;
			double newalpha_i, newalpha_j;
			if(base == 0 && (fabs(-QD[j]*G_i + Q_ij * G_j) < 1.0e-12 || fabs(-QD[i]*G_j + Q_ij * G_i <1.0e-12)))
			{
				double delta = QD[i]*alpha[i] + Q_ij*alpha[j] - G_i;
				if(Q_ij < 0)
				{
					if(QD[i]*C_i <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = 0;
					}
					else if( Q_ij*C_j >= delta )
					{
						newalpha_i = 0;
						newalpha_j = C_j;
					}
					else if(0 <= delta)
					{
						newalpha_i = delta/QD[i];
						newalpha_j = 0;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
				}
				else
				{
					if( 0 >= delta)
					{
						newalpha_i = 0;
						newalpha_j = 0;
					}
					else if(QD[i]*C_i + Q_ij*C_j <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = C_j;
					}
					else if( Q_ij*C_j <= delta)
					{
						newalpha_i = (delta-Q_ij*C_j)/QD[i];
						newalpha_j = C_j;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
					
				}
			}
			else
			{
				bool use_j = false;
				newalpha_i = min(C_i, max(0.0, alpha[i] + (-QD[j]*G_i + Q_ij * G_j)/base));
				newalpha_j = min(C_j, max(0.0, alpha[j] + (-QD[i]*G_j + Q_ij * G_i)/base));
				double d_i = newalpha_i - alpha[i];
				if(newalpha_i >= C_i)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i<=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;

				}
				else if(newalpha_i <= 0.0)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i>=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;
				}
				else{
					use_j = true;
				}
				if(use_j == true)
					newalpha_i = min(C_i, max(0.0, alpha[i] - (Q_ij*(newalpha_j - alpha[j]) + G_i)/QD[i]));
			}
			bool update = false;
			if(fabs(newalpha_i - alpha[i]) > 1.0e-16)
			{
				sparse_operator::axpy(yi*(newalpha_i - alpha[i]), prob->x[i], w);
				alpha[i] = newalpha_i;
				update = true;
			}
			if(fabs(newalpha_j - alpha[j]) > 1.0e-16)
			{
				sparse_operator::axpy(yj*(newalpha_j - alpha[j]), prob->x[j], w);
				alpha[j] = newalpha_j;
				update = true;
			}
			if(update)
				success_pair++;
			
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
}

void Solver::two_semirandom2_1000()
{
	int l = max_set;
	int i, j, si;
	active_size = l;
	double C_i, C_j, G_i, G_j;
	clock_t start;
	active_size = active_size - 1;
	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		for (i=0; i<active_size+1; i++)
		{
			int j = i+rand()%(active_size-i+1);
			swap(index[i], index[j]);
		}
		for(si=0; si<active_size; si++)
		{
			i = index[si];
			j = index[si+1];

			const schar yi = y[i];
			const schar yj = y[j];
			
			C_i = upper_bound[GETI(i)];
			C_j = upper_bound[GETI(j)];
			

			feature_node const *xi = prob->x[i];
			feature_node const *xj = prob->x[j];
			
			G_i = 0;
			G_j = 0;
			double Q_ij = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;
				
			double base = QD[i]*QD[j] - Q_ij*Q_ij;
			double newalpha_i, newalpha_j;
			if(base == 0 && (fabs(-QD[j]*G_i + Q_ij * G_j) < 1.0e-12 || fabs(-QD[i]*G_j + Q_ij * G_i <1.0e-12)))
			{
				double delta = QD[i]*alpha[i] + Q_ij*alpha[j] - G_i;
				if(Q_ij < 0)
				{
					if(QD[i]*C_i <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = 0;
					}
					else if( Q_ij*C_j >= delta )
					{
						newalpha_i = 0;
						newalpha_j = C_j;
					}
					else if(0 <= delta)
					{
						newalpha_i = delta/QD[i];
						newalpha_j = 0;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
				}
				else
				{
					if( 0 >= delta)
					{
						newalpha_i = 0;
						newalpha_j = 0;
					}
					else if(QD[i]*C_i + Q_ij*C_j <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = C_j;
					}
					else if( Q_ij*C_j <= delta)
					{
						newalpha_i = (delta-Q_ij*C_j)/QD[i];
						newalpha_j = C_j;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
					
				}
			}
			else
			{
				bool use_j = false;
				newalpha_i = min(C_i, max(0.0, alpha[i] + (-QD[j]*G_i + Q_ij * G_j)/base));
				newalpha_j = min(C_j, max(0.0, alpha[j] + (-QD[i]*G_j + Q_ij * G_i)/base));
				double d_i = newalpha_i - alpha[i];
				if(newalpha_i >= C_i)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i<=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;

				}
				else if(newalpha_i <= 0.0)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i>=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;
				}
				else{
					use_j = true;
				}
				if(use_j == true)
					newalpha_i = min(C_i, max(0.0, alpha[i] - (Q_ij*(newalpha_j - alpha[j]) + G_i)/QD[i]));
			}
			bool update = false;
			if(fabs(newalpha_i - alpha[i]) > 1.0e-16)
			{
				sparse_operator::axpy(yi*(newalpha_i - alpha[i]), prob->x[i], w);
				alpha[i] = newalpha_i;
				update = true;
			}
			if(fabs(newalpha_j - alpha[j]) > 1.0e-16)
			{
				sparse_operator::axpy(yj*(newalpha_j - alpha[j]), prob->x[j], w);
				alpha[j] = newalpha_j;
				update = true;
			}
			if(update)
				success_pair++;
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
}

void Solver::two_semirandom1_1000()
{
	int l = max_set;
	int i, j, si;
	active_size = l;
	double C_i, C_j, G_i, G_j;
	clock_t start;
	int hactive_size = active_size;
	active_size = active_size/2;
	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		for (i=0; i<hactive_size; i++)
		{
			int j = i+rand()%(hactive_size-i);
			swap(index[i], index[j]);
		}
		for(si=0; si<hactive_size-1; si+=2)
		{
			i = index[si];
			j = index[si+1];

			const schar yi = y[i];
			const schar yj = y[j];
			
			C_i = upper_bound[GETI(i)];
			C_j = upper_bound[GETI(j)];
			

			feature_node const *xi = prob->x[i];
			feature_node const *xj = prob->x[j];
			
			G_i = 0;
			G_j = 0;
			double Q_ij = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;
				
			double base = QD[i]*QD[j] - Q_ij*Q_ij;
			double newalpha_i, newalpha_j;
			if(base == 0 && (fabs(-QD[j]*G_i + Q_ij * G_j) < 1.0e-12 || fabs(-QD[i]*G_j + Q_ij * G_i <1.0e-12)))
			{
				double delta = QD[i]*alpha[i] + Q_ij*alpha[j] - G_i;
				if(Q_ij < 0)
				{
					if(QD[i]*C_i <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = 0;
					}
					else if( Q_ij*C_j >= delta )
					{
						newalpha_i = 0;
						newalpha_j = C_j;
					}
					else if(0 <= delta)
					{
						newalpha_i = delta/QD[i];
						newalpha_j = 0;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
				}
				else
				{
					if( 0 >= delta)
					{
						newalpha_i = 0;
						newalpha_j = 0;
					}
					else if(QD[i]*C_i + Q_ij*C_j <= delta)
					{
						newalpha_i = C_i;
						newalpha_j = C_j;
					}
					else if( Q_ij*C_j <= delta)
					{
						newalpha_i = (delta-Q_ij*C_j)/QD[i];
						newalpha_j = C_j;
					}
					else
					{
						newalpha_i = 0;
						newalpha_j = delta/Q_ij;
					}
					
				}
			}
			else
			{
				bool use_j = false;
				newalpha_i = min(C_i, max(0.0, alpha[i] + (-QD[j]*G_i + Q_ij * G_j)/base));
				newalpha_j = min(C_j, max(0.0, alpha[j] + (-QD[i]*G_j + Q_ij * G_i)/base));
				double d_i = newalpha_i - alpha[i];
				if(newalpha_i >= C_i)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i<=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;

				}
				else if(newalpha_i <= 0.0)
				{
					if(QD[i]*d_i + Q_ij*(newalpha_j - alpha[j]) + G_i>=0)
						newalpha_j = min(C_j, max(0.0, alpha[j] - (Q_ij*d_i + G_j)/QD[j]));
					else
						use_j = true;
				}
				else{
					use_j = true;
				}
				if(use_j == true)
					newalpha_i = min(C_i, max(0.0, alpha[i] - (Q_ij*(newalpha_j - alpha[j]) + G_i)/QD[i]));
			}
			bool update = false;
			if(fabs(newalpha_i - alpha[i]) > 1.0e-16)
			{
				sparse_operator::axpy(yi*(newalpha_i - alpha[i]), prob->x[i], w);
				alpha[i] = newalpha_i;
				update = true;
			}
			if(fabs(newalpha_j - alpha[j]) > 1.0e-16)
			{
				sparse_operator::axpy(yj*(newalpha_j - alpha[j]), prob->x[j], w);
				alpha[j] = newalpha_j;
				update = true;
			}
			if(update)
				success_pair++;
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
}

static void calculate_unbias_two_newalpha(
		std::pair<double,double> *newalpha_ij,
		double QDi, double QDj, double Q_ij, double C_i, double C_j,
		double alpha_i, double alpha_j, double G_i, double G_j)
{
	double base = QDi*QDj - Q_ij*Q_ij;
	if(base == 0 && (fabs(-QDj*G_i + Q_ij * G_j) < 1.0e-12 || fabs(-QDi*G_j + Q_ij * G_i <1.0e-12)))
	{
		double delta = QDi*alpha_i + Q_ij*alpha_j - G_i;
		if(Q_ij < 0)
		{
			if(QDi*C_i <= delta)
			{
				newalpha_ij->first = C_i;
				newalpha_ij->second = 0;
			}
			else if( Q_ij*C_j >= delta )
			{
				newalpha_ij->first = 0;
				newalpha_ij->second = C_j;
			}
			else if(0 <= delta)
			{
				newalpha_ij->first = delta/QDi;
				newalpha_ij->second = 0;
			}
			else
			{
				newalpha_ij->first = 0;
				newalpha_ij->second = delta/Q_ij;
			}
		}
		else
		{
			if( 0 >= delta)
			{
				newalpha_ij->first = 0;
				newalpha_ij->second = 0;
			}
			else if(QDi*C_i + Q_ij*C_j <= delta)
			{
				newalpha_ij->first = C_i;
				newalpha_ij->second = C_j;
			}
			else if( Q_ij*C_j <= delta)
			{
				newalpha_ij->first = (delta-Q_ij*C_j)/QDi;
				newalpha_ij->second = C_j;
			}
			else
			{
				newalpha_ij->first = 0;
				newalpha_ij->second = delta/Q_ij;
			}
			
		}
	}
	else
	{
		bool use_j = false;
		newalpha_ij->first = min(C_i, max(0.0, alpha_i + (-QDj*G_i + Q_ij * G_j)/base));
		newalpha_ij->second = min(C_j, max(0.0, alpha_j + (-QDi*G_j + Q_ij * G_i)/base));
		double d_i = newalpha_ij->first - alpha_i;
		if(newalpha_ij->first >= C_i)
		{
			if(QDi*d_i + Q_ij*(newalpha_ij->second - alpha_j) + G_i<=0)
				newalpha_ij->second = min(C_j, max(0.0, alpha_j - (Q_ij*d_i + G_j)/QDj));
			else
				use_j = true;
		}
		else if(newalpha_ij->first <= 0.0)
		{
			if(QDi*d_i + Q_ij*(newalpha_ij->second - alpha_j) + G_i>=0)
				newalpha_ij->second = min(C_j, max(0.0, alpha_j - (Q_ij*d_i + G_j)/QDj));
			else
				use_j = true;
		}
		else
		{
			use_j = true;
		}
		if(use_j == true)
			newalpha_ij->first = min(C_i, max(0.0, alpha_i - (Q_ij*(newalpha_ij->second - alpha_j) + G_i)/QDi));
	}
}

void Solver::two_random_1000()
{
	int l = max_set;
	int i, j, si;
	active_size = l;
	double C_i, C_j, G_i, G_j;
	double Q_ij;
	bool is_updated;
	std::pair<double,double> *newalpha_ij = new std::pair<double,double>;

	clock_t start;
	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		for(si=0; si<active_size; si++)
		{
			i = index[rand()%active_size];
			j = index[rand()%active_size];
			
			while(i==j)
			{
				j = index[rand()%active_size];
			}

			const schar yi = y[i];
			const schar yj = y[j];
			const feature_node *xi = prob->x[i];
			const feature_node *xj = prob->x[j];
			
			Q_ij = yi * yj * sparse_operator::feature_dot(xi, xj);
			G_i = yi*sparse_operator::dot(w, xi) -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*sparse_operator::dot(w, xj) -1 +alpha[j]*diag[GETI(j)];
				
			C_i = upper_bound[GETI(i)];
			C_j = upper_bound[GETI(j)];
			
			calculate_unbias_two_newalpha(newalpha_ij, QD[i],QD[j],Q_ij,C_i,C_j,alpha[i],alpha[j],G_i,G_j);

			is_updated = false;
			if(fabs(newalpha_ij->first - alpha[i]) > 1.0e-16)
			{
				sparse_operator::axpy(yi*(newalpha_ij->first - alpha[i]), prob->x[i], w);
				alpha[i] = newalpha_ij->first;
				is_updated = true;
			}
			if(fabs(newalpha_ij->second - alpha[j]) > 1.0e-16)
			{
				sparse_operator::axpy(yj*(newalpha_ij->second - alpha[j]), prob->x[j], w);
				alpha[j] = newalpha_ij->second;
				is_updated = true;
			}
			if(is_updated)
				success_pair++;
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
}

void Solver::two_semigd_1000()
{
	// TODO:implement
	info("Not implement yet\n");
	int l = max_set;
	int i, j, si;
	active_size = l;
	double C_i, C_j, G_i, G_j;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	double *G = new double[l];
	double PG;
	double Q_ij;
	bool is_updated;
	std::pair<double,double> *newalpha_ij = new std::pair<double,double>;

	clock_t start;

	int *Max_order_index = new int[l];
	int *Min_order_index = new int[l];
	
	std::priority_queue<struct feature_node, std::vector<feature_node>, mincomp> Min_order_queue;
	std::priority_queue<struct feature_node, std::vector<feature_node>, maxcomp> Max_order_queue;

	while (iter < max_iter)
	{
		start = clock();
		success_pair = 0;

		for(i=0; i<active_size; i++)
		{
			int G_index = index[i];
			feature_node * const xi = prob->x[G_index];
			G[i] = y[G_index]*sparse_operator::dot(w, xi) -1 +alpha[i]*diag[GETI(i)];
		}

		update_size = int(active_size*ratio_update);
		if(update_size < 1)
			update_size = 1;
		update_size = 1;

		for(i=0; i<active_size; i++)
		{
			struct feature_node comp;
			comp.index = index[i];

//			PG=0;
//			if(alpha_status[comp.index] == LOWER_BOUND)
//			{
//				if(G[i] < 0)
//					PG = G[i];
//			}
//			else if (alpha_status[comp.index] == UPPER_BOUND)
//			{
//				if(G[i] > 0)
//					PG = G[i];
//			}
//			else
//				PG = G[i];
//
//			PGmax_new = max(PGmax_new, PG);
//			PGmin_new = min(PGmin_new, PG);

			//comp.value = fabs(PG);

//			comp.value = PG;
//			if((int) Min_order_queue.size() < update_size)
//			{
//				Min_order_queue.push(comp);
//			}
//			else if(Min_order_queue.top().value < comp.value)
//			{
//				Min_order_queue.pop();
//				Min_order_queue.push(comp);
//			}
//			if((int) Max_order_queue.size() < update_size)
//			{
//				Max_order_queue.push(comp);
//			}
//			else if(Max_order_queue.top().value > comp.value)
//			{
//				Max_order_queue.pop();
//				Max_order_queue.push(comp);
//			}

//			int yc = y[comp.index];
//			comp.value = -yc * G[i];
//			if( (alpha_status[comp.index] != UPPER_BOUND && yc==+1) ||
//					(alpha_status[comp.index] != LOWER_BOUND && yc==-1) )
//			{ //Iup
//				if((int) Min_order_queue.size() < update_size)
//				{
//					Min_order_queue.push(comp);
//				}
//				else if(Min_order_queue.top().value < comp.value)
//				{
//					Min_order_queue.pop();
//					Min_order_queue.push(comp);
//				}
//			}
//			if( (alpha_status[comp.index] != UPPER_BOUND && yc==-1) ||
//					(alpha_status[comp.index] != LOWER_BOUND && yc==+1) )
//			{ //Ilow
//				if((int) Max_order_queue.size() < update_size)
//				{
//					Max_order_queue.push(comp);
//				}
//				else if(Max_order_queue.top().value > comp.value)
//				{
//					Max_order_queue.pop();
//					Max_order_queue.push(comp);
//				}
//			}
		}
		update_size = min((int)Min_order_queue.size(), (int)Max_order_queue.size());
		while((int)Max_order_queue.size() > update_size)
		{
			Max_order_queue.pop();
		}
		while((int)Min_order_queue.size() > update_size)
		{
			Min_order_queue.pop();
		}
		for(i=0; i<update_size; i++)
		{
			Max_order_index[update_size-1-i] = Min_order_queue.top().index;
			Min_order_queue.pop();

			Min_order_index[update_size-1-i] = Max_order_queue.top().index;
			Max_order_queue.pop();
		}

		for(si=0; si<update_size; si++)
		{
			i = Max_order_index[si];
			j = Min_order_index[si];

			const schar yi = y[i];
			const schar yj = y[j];
			const feature_node *xi = prob->x[i];
			const feature_node *xj = prob->x[j];
			
			Q_ij = yi * yj * sparse_operator::feature_dot(xi, xj);
			G_i = yi*sparse_operator::dot(w, xi) -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*sparse_operator::dot(w, xj) -1 +alpha[j]*diag[GETI(j)];

			C_i = upper_bound[GETI(i)];
			C_j = upper_bound[GETI(j)];
			
			calculate_unbias_two_newalpha(newalpha_ij, QD[i],QD[j],Q_ij,C_i,C_j,alpha[i],alpha[j],G_i,G_j);

			is_updated = false;
			if(fabs(newalpha_ij->first - alpha[i]) > 1.0e-16)
			{
				sparse_operator::axpy(yi*(newalpha_ij->first - alpha[i]), prob->x[i], w);
				alpha[i] = newalpha_ij->first;
				is_updated = true;
			}
			if(fabs(newalpha_ij->second - alpha[j]) > 1.0e-16)
			{
				sparse_operator::axpy(yj*(newalpha_ij->second - alpha[j]), prob->x[j], w);
				alpha[j] = newalpha_ij->second;
				is_updated = true;
			}
			if(is_updated)
				success_pair++;
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
	// TODO:test
}

void Solver::two_semigd_shrink()
{
	info("THIS IS two semigd sh\n");
	info("Not implement yet\n");
	// TODO:implement
	// TODO:test
}

void Solver::bias_semigd_shrink()
{
	int l = prob->l;
	int i, j;
	active_size = l;
	double G_i, G_j;
	int counter = min(l,3);
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	Gmax = -INF;
	Gmin = INF;
	time_t start;
	
	int *Max_order_index = new int[l];
	int *Min_order_index = new int[l];
	
	std::priority_queue<struct feature_node, std::vector<feature_node>, mincomp> Min_order_queue;
	std::priority_queue<struct feature_node, std::vector<feature_node>, maxcomp> Max_order_queue;

	double *G = new double[l];
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		Gmax = -INF;
		Gmin = INF;

		for(i=0; i<active_size; i++)
		{
			int G_index = index[i];
			feature_node * const xi = prob->x[G_index];
			G[i] = -y[G_index] * (y[G_index]*sparse_operator::dot(w, xi) -1 +alpha[G_index]*diag[GETI(G_index)]);
			
			if(y[G_index] == 1)
			{
				if(alpha_status[G_index] != UPPER_BOUND)
					Gmax = max(Gmax, G[i]);
				if(alpha_status[G_index] != LOWER_BOUND)
					Gmin = min(Gmin, G[i]);
			}
			else
			{
				if(alpha_status[G_index] != UPPER_BOUND)
					Gmin = min(Gmin, G[i]);
				if(alpha_status[G_index] != LOWER_BOUND)
					Gmax = max(Gmax, G[i]);
			}
		}

		if(Gmax - Gmin < eps && iter>1)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				counter = 0;
				continue;
			}
		}
		else if(counter-- == 0) //do shrinking
		{
			for(i=active_size-1; i>=0; i--)
			{
				bool be_shunk = false;
				int sh_idx = index[i];
				if(y[sh_idx] == 1)
				{
					if(alpha_status[sh_idx] == UPPER_BOUND)
					{
						if(G[i] > Gmax)
							be_shunk = true;
					}
					else if(alpha_status[sh_idx] == LOWER_BOUND)
					{
						if(G[i] < Gmin)
							be_shunk = true;
					}
				}
				else if(y[sh_idx] == -1)
				{
					if(alpha_status[sh_idx] == LOWER_BOUND)
					{
						if(G[i] > Gmax)
							be_shunk = true;
					}
					else if(alpha_status[sh_idx] == UPPER_BOUND)
					{
						if(G[i] < Gmin)
							be_shunk = true;
					}
				}
				if(be_shunk)
				{
					swap(index[i],index[--active_size]);
					swap(G[i], G[active_size]);
				}
				

			}
			counter = min(l,1);
		}
		update_size = int(active_size*ratio_update);
		if(update_size < 1)
		update_size = 1;

		for(i=0; i<active_size; i++)
		{
			struct feature_node comp;
			comp.index = index[i];
			comp.value = G[i];
			int yc = y[comp.index];
			if( (alpha_status[comp.index] != UPPER_BOUND && yc==+1) ||
					(alpha_status[comp.index] != LOWER_BOUND && yc==-1) )
			{
				if((int) Min_order_queue.size() < update_size)
				{
					Min_order_queue.push(comp);
				}
				else
				{
					if(Min_order_queue.top().value < comp.value)
					{
						Min_order_queue.pop();
						Min_order_queue.push(comp);
					}
				}
			}
			if( (alpha_status[comp.index] != UPPER_BOUND && yc==-1) ||
					(alpha_status[comp.index] != LOWER_BOUND && yc==+1) )
			{
				if((int) Max_order_queue.size() < update_size)
				{
					Max_order_queue.push(comp);
				}
				else
				{
					if(Max_order_queue.top().value > comp.value)
					{
						Max_order_queue.pop();
						Max_order_queue.push(comp);
					}
				}
			}
		}
		update_size = min((int)Min_order_queue.size(), (int)Max_order_queue.size());
		while((int)Max_order_queue.size() > update_size)
		{
			Max_order_queue.pop();
		}
		while((int)Min_order_queue.size() > update_size)
		{
			Min_order_queue.pop();
		}

		for(i=0; i<update_size; i++)
		{
			Max_order_index[update_size-1-i] = Min_order_queue.top().index;
			Min_order_queue.pop();

			Min_order_index[update_size-1-i] = Max_order_queue.top().index;
			Max_order_queue.pop();
		}
		for(int index_i = 0; index_i < update_size; index_i++)
		{
			i = Max_order_index[index_i];
			j = Min_order_index[index_i];

			feature_node const * xi = prob->x[i];
			feature_node const * xj = prob->x[j];
			
			const schar yi = y[i];
			const schar yj = y[j];
			double C_i = upper_bound[GETI(i)];
			double C_j = upper_bound[GETI(j)];
			
			G_i = 0;
			G_j = 0;
			double Q_ij = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;
			
			double old_alpha_i = alpha[i];
			double old_alpha_j = alpha[j];

			if(y[i] == y[j])
			{
				double quad_coef = QD[i] + QD[j] - 2*Q_ij;
				if(quad_coef <= 0)
					quad_coef = 1e-12;
				double delta = (G_i-G_j)/quad_coef;
				double sum = alpha[i] + alpha[j];
				alpha[i] -= delta;
				alpha[j] += delta;

				if(sum > C_i)
				{
					if(alpha[i] > C_i)
					{
						alpha[i] = C_i;
						alpha[j] = sum -C_i;
					}
				}
				else
				{
					if(alpha[j] < 0)
					{
						alpha[j] = 0;
						alpha[i] = sum;
					}
				}
				if(sum > C_j)
				{
					if(alpha[j] > C_j)
					{
						alpha[j] = C_j;
						alpha[i] = sum -C_j;
					}
				}
				else
				{
					if(alpha[i] < 0)
					{
						alpha[i] = 0;
						alpha[j] = sum;
					}
				}
			}
			else
			{
				double quad_coef = QD[i] + QD[j] + 2*Q_ij;
				if(quad_coef <= 0)
					quad_coef = 1e-12;
				double delta = (-G_i-G_j)/quad_coef;
				double diff = alpha[i] - alpha[j];
				alpha[i] += delta;
				alpha[j] += delta;
				if(diff > 0)
				{
					if(alpha[j] < 0)
					{
						alpha[j] = 0;
						alpha[i] = diff;
					}
				}
				else
				{
					if(alpha[i] < 0)
					{
						alpha[i] = 0;
						alpha[j] = -diff;
					}
				}
				if(diff > C_i - C_j)
				{
					if(alpha[i] > C_i)
					{
						alpha[i] = C_i;
						alpha[j] = C_i - diff;
					}
				}
				else
				{
					if(alpha[j] > C_j)
					{
						alpha[j] = C_j;
						alpha[i] = C_j + diff;
					}
				}
			}

			// update alpha status and w
			if(fabs(alpha[i]-old_alpha_i) > 1e-16)
			{
				success_pair++;
				sparse_operator::axpy(y[i]*(alpha[i]-old_alpha_i), prob->x[i], w);
				sparse_operator::axpy(y[j]*(alpha[j]-old_alpha_j), prob->x[j], w);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[GETI(i)]);
				alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[GETI(j)]);
			}
			else
			{
				alpha[i] = old_alpha_i;
				alpha[j] = old_alpha_j;
			}
		}
		iter++;
		duration += clock() - start;
		log_message();

	}
	summary();	
	delete [] alpha_status;
}

void Solver::bias_semigd_1000()
{
	int l = prob->l;
	int i, j;
	active_size = l;
	double G_i, G_j;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	time_t start;
	
	int *Max_order_index = new int[l];
	int *Min_order_index = new int[l];
	
	std::priority_queue<struct feature_node, std::vector<feature_node>, mincomp> Min_order_queue;
	std::priority_queue<struct feature_node, std::vector<feature_node>, maxcomp> Max_order_queue;

	double *G = new double[l];
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		nr_pos_y = 0;
		nr_neg_y = 0;

		for(i=0; i<active_size; i++)
		{
			int G_index = index[i];
			feature_node * const xi = prob->x[G_index];
			G[i] = -y[G_index] * (y[G_index]*sparse_operator::dot(w, xi) -1 +alpha[i]*diag[GETI(i)]);
		}

		//TODO: rlist as number not ratio
		//update_size = int(active_size*ratio_update);
		update_size = int(ratio_update/2);
		if(update_size < 1)
			update_size = 1;

		for(i=0; i<active_size; i++)
		{
			struct feature_node comp;
			comp.index = index[i];
			comp.value = G[i];
			int yc = y[comp.index];
			if( (alpha_status[comp.index] != UPPER_BOUND && yc==+1) ||
					(alpha_status[comp.index] != LOWER_BOUND && yc==-1) )
			{
				if((int) Min_order_queue.size() < update_size)
				{
					Min_order_queue.push(comp);
				}
				else
				{
					if(Min_order_queue.top().value < comp.value)
					{
						Min_order_queue.pop();
						Min_order_queue.push(comp);
					}
				}
			}
			if( (alpha_status[comp.index] != UPPER_BOUND && yc==-1) ||
					(alpha_status[comp.index] != LOWER_BOUND && yc==+1) )
			{
				if((int) Max_order_queue.size() < update_size)
				{
					Max_order_queue.push(comp);
				}
				else
				{
					if(Max_order_queue.top().value > comp.value)
					{
						Max_order_queue.pop();
						Max_order_queue.push(comp);
					}
				}
			}
		}
		update_size = min((int)Min_order_queue.size(), (int)Max_order_queue.size());
		while((int)Max_order_queue.size() > update_size)
		{
			Max_order_queue.pop();
		}
		while((int)Min_order_queue.size() > update_size)
		{
			Min_order_queue.pop();
		}

		for(i=0; i<update_size; i++)
		{
			
			Max_order_index[update_size-1-i] = Min_order_queue.top().index;
		//info("%d %g\n",Min_order_queue.top().index, Min_order_queue.top().value);
			Min_order_queue.pop();

			Min_order_index[update_size-1-i] = Max_order_queue.top().index;
		//info("%d %g\n",Max_order_queue.top().index, Max_order_queue.top().value);
			Max_order_queue.pop();
		}
		// std::random_shuffle(Max_order_index, Max_order_index+update_size);
		// std::random_shuffle(Min_order_index, Min_order_index+update_size);
		for(int index_i = 0; index_i < update_size; index_i++)
		{
			i = Max_order_index[index_i];
			j = Min_order_index[index_i];
			//i = rand()%active_size;
			//j = rand()%active_size;

			feature_node const * xi = prob->x[i];
			feature_node const * xj = prob->x[j];
			
			const schar yi = y[i];
			if(yi == +1)
				++nr_pos_y;
			if(yi == -1)
				++nr_neg_y;
			const schar yj = y[j];
			if(yj == +1)
				++nr_pos_y;
			if(yj == -1)
				++nr_neg_y;
			double C_i = upper_bound[GETI(i)];
			double C_j = upper_bound[GETI(j)];
			
			G_i = 0;
			G_j = 0;
			double Q_ij = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;
			
			double old_alpha_i = alpha[i];
			double old_alpha_j = alpha[j];

			if(y[i] == y[j])
			{
				double quad_coef = QD[i] + QD[j] - 2*Q_ij;
				if(quad_coef <= 0)
					quad_coef = 1e-12;
				double delta = (G_i-G_j)/quad_coef;
				double sum = alpha[i] + alpha[j];
				alpha[i] -= delta;
				alpha[j] += delta;

				if(sum > C_i)
				{
					if(alpha[i] > C_i)
					{
						alpha[i] = C_i;
						alpha[j] = sum -C_i;
					}
				}
				else
				{
					if(alpha[j] < 0)
					{
						alpha[j] = 0;
						alpha[i] = sum;
					}
				}
				if(sum > C_j)
				{
					if(alpha[j] > C_j)
					{
						alpha[j] = C_j;
						alpha[i] = sum -C_j;
					}
				}
				else
				{
					if(alpha[i] < 0)
					{
						alpha[i] = 0;
						alpha[j] = sum;
					}
				}
			}
			else
			{
				double quad_coef = QD[i] + QD[j] + 2*Q_ij;
				if(quad_coef <= 0)
					quad_coef = 1e-12;
				double delta = (-G_i-G_j)/quad_coef;
				double diff = alpha[i] - alpha[j];
				alpha[i] += delta;
				alpha[j] += delta;
				if(diff > 0)
				{
					if(alpha[j] < 0)
					{
						alpha[j] = 0;
						alpha[i] = diff;
					}
				}
				else
				{
					if(alpha[i] < 0)
					{
						alpha[i] = 0;
						alpha[j] = -diff;
					}
				}
				if(diff > C_i - C_j)
				{
					if(alpha[i] > C_i)
					{
						alpha[i] = C_i;
						alpha[j] = C_i - diff;
					}
				}
				else
				{
					if(alpha[j] > C_j)
					{
						alpha[j] = C_j;
						alpha[i] = C_j + diff;
					}
				}
			}

			// update alpha status and w
			if(fabs(alpha[i]-old_alpha_i) > 1e-16)
			{
				success_pair++;
				sparse_operator::axpy(y[i]*(alpha[i]-old_alpha_i), prob->x[i], w);
				sparse_operator::axpy(y[j]*(alpha[j]-old_alpha_j), prob->x[j], w);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[GETI(i)]);
				alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[GETI(j)]);
			}
			else
			{
				alpha[i] = old_alpha_i;
				alpha[j] = old_alpha_j;
			}
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();	
	delete [] alpha_status;
}

void Solver::bias_random_shrink()
{
	int l = prob->l;
	int i, j;
	active_size = l;
	double G_i, G_j;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	double PGmax_old = INF;
	double PGmin_old = -INF;
	time_t start;
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		PGmax_new = -INF;
		PGmin_new = INF;
		for(int s = 0; s < active_size; s++)
		{
			int si = rand()%active_size;
			int sj = rand()%active_size;
			
			while(si==sj)
			{
				sj = rand()%active_size;
			}

			i = index[si];
			j = index[sj];

			feature_node const * xi = prob->x[i];
			feature_node const * xj = prob->x[j];
			
			const schar yi = y[i];
			const schar yj = y[j];
			double C_i = upper_bound[GETI(i)];
			double C_j = upper_bound[GETI(j)];
		
			double Q_ij = 0;
			G_i = 0;
			G_j = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			xi = prob->x[i];
			xj = prob->x[j];
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;
			double yG_i = -yi * G_i;
			double yG_j = -yj * G_j;
			bool to_be_shrunk = false;
			//shrink
			if( (alpha_status[i] == UPPER_BOUND && yi == +1) ||
				(alpha_status[i] == LOWER_BOUND && yi == -1) )
			{
				if(yG_i > PGmax_old)
				{
					active_size--;
					swap(index[si], index[active_size]);
					to_be_shrunk = true;
					continue;
				}
				else
					PGmin_new = min(PGmin_new, yG_i);
			}
			else if( (alpha_status[i] == UPPER_BOUND && yi == -1) ||
				(alpha_status[i] == LOWER_BOUND && yi == +1) )
			{
				if(yG_i < PGmin_old)
				{
					active_size--;
					swap(index[si], index[active_size]);
					to_be_shrunk = true;
					continue;
				}
				else
					PGmax_new = max(PGmax_new, yG_i);
			}
			else
			{
				PGmax_new = max(PGmax_new, yG_i);
				PGmin_new = min(PGmin_new, yG_i);
			}
			if( (alpha_status[j] == UPPER_BOUND && yj == +1) ||
				(alpha_status[j] == LOWER_BOUND && yj == -1) )
			{
				if(yG_j > PGmax_old)
				{
					active_size--;
					swap(index[sj], index[active_size]);
					to_be_shrunk = true;
					continue;
				}
				else
					PGmin_new = min(PGmin_new, yG_j);
			}
			else if( (alpha_status[j] == UPPER_BOUND && yj == -1) ||
				(alpha_status[j] == LOWER_BOUND && yj == +1) )
			{
				if(yG_j < PGmin_old)
				{
					active_size--;
					swap(index[sj], index[active_size]);
					to_be_shrunk = true;
					continue;
				}
				else
					PGmax_new = max(PGmax_new, yG_j);
			}
			else
			{
				PGmax_new = max(PGmax_new, yG_j);
				PGmin_new = min(PGmin_new, yG_j);
			}
			if(to_be_shrunk)
				continue;
			
			double old_alpha_i = alpha[i];
			double old_alpha_j = alpha[j];

			if(y[i] == y[j])
			{
				double quad_coef = QD[i] + QD[j] - 2*Q_ij;
				if(quad_coef <= 0)
					quad_coef = 1e-12;
				double delta = (G_i-G_j)/quad_coef;
				double sum = alpha[i] + alpha[j];
				alpha[i] -= delta;
				alpha[j] += delta;

				if(sum > C_i)
				{
					if(alpha[i] > C_i)
					{
						alpha[i] = C_i;
						alpha[j] = sum -C_i;
					}
				}
				else
				{
					if(alpha[j] < 0)
					{
						alpha[j] = 0;
						alpha[i] = sum;
					}
				}
				if(sum > C_j)
				{
					if(alpha[j] > C_j)
					{
						alpha[j] = C_j;
						alpha[i] = sum -C_j;
					}
				}
				else
				{
					if(alpha[i] < 0)
					{
						alpha[i] = 0;
						alpha[j] = sum;
					}
				}
			}
			else
			{
				double quad_coef = QD[i] + QD[j] + 2*Q_ij;
				if(quad_coef <= 0)
					quad_coef = 1e-12;
				double delta = (-G_i-G_j)/quad_coef;
				double diff = alpha[i] - alpha[j];
				alpha[i] += delta;
				alpha[j] += delta;
				if(diff > 0)
				{
					if(alpha[j] < 0)
					{
						alpha[j] = 0;
						alpha[i] = diff;
					}
				}
				else
				{
					if(alpha[i] < 0)
					{
						alpha[i] = 0;
						alpha[j] = -diff;
					}
				}
				if(diff > C_i - C_j)
				{
					if(alpha[i] > C_i)
					{
						alpha[i] = C_i;
						alpha[j] = C_i - diff;
					}
				}
				else
				{
					if(alpha[j] > C_j)
					{
						alpha[j] = C_j;
						alpha[i] = C_j + diff;
					}
				}
			}

			// update alpha status and w
			if(fabs(alpha[i]-old_alpha_i) > 1e-16)
			{
				success_pair++;
				sparse_operator::axpy(y[i]*(alpha[i]-old_alpha_i), prob->x[i], w);
				sparse_operator::axpy(y[j]*(alpha[j]-old_alpha_j), prob->x[j], w);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[GETI(i)]);
				alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[GETI(j)]);
			}
			else
			{
				alpha[i] = old_alpha_i;
				alpha[j] = old_alpha_j;
			}
		}
		if(PGmax_new - PGmin_new <= eps)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				PGmax_old = INF;
				PGmin_old = -INF;
				continue;
			}
		}
		PGmax_old = PGmax_new;
		PGmin_old = PGmin_new;
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
	delete [] alpha_status;
}

void Solver::bias_random_1000()
{
	int l = prob->l;
	int i, j;
	active_size = l;
	double G_i, G_j;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	time_t start;
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		for(int si = 0; si < active_size; si++)
		{
			i = rand()%active_size;
			j = rand()%active_size;
			
			feature_node const * xi = prob->x[i];
			feature_node const * xj = prob->x[j];
			
			const schar yi = y[i];
			const schar yj = y[j];
			double C_i = upper_bound[GETI(i)];
			double C_j = upper_bound[GETI(j)];
			
			double Q_ij = 0;
			G_i = 0;
			G_j = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			G_i = yi*G_i -1 +alpha[i]*diag[GETI(i)];
			G_j = yj*G_j -1 +alpha[j]*diag[GETI(j)];
			Q_ij = yi * yj * Q_ij;

			double old_alpha_i = alpha[i];
			double old_alpha_j = alpha[j];

			if(y[i] == y[j])
			{
				double quad_coef = QD[i] + QD[j] - 2*Q_ij;
				if(quad_coef <= 0)
					quad_coef = 1e-12;
				double delta = (G_i-G_j)/quad_coef;
				double sum = alpha[i] + alpha[j];
				alpha[i] -= delta;
				alpha[j] += delta;

				if(sum > C_i)
				{
					if(alpha[i] > C_i)
					{
						alpha[i] = C_i;
						alpha[j] = sum -C_i;
					}
				}
				else
				{
					if(alpha[j] < 0)
					{
						alpha[j] = 0;
						alpha[i] = sum;
					}
				}
				if(sum > C_j)
				{
					if(alpha[j] > C_j)
					{
						alpha[j] = C_j;
						alpha[i] = sum -C_j;
					}
				}
				else
				{
					if(alpha[i] < 0)
					{
						alpha[i] = 0;
						alpha[j] = sum;
					}
				}
			}
			else
			{
				double quad_coef = QD[i] + QD[j] + 2*Q_ij;
				if(quad_coef <= 0)
					quad_coef = 1e-12;
				double delta = (-G_i-G_j)/quad_coef;
				double diff = alpha[i] - alpha[j];
				alpha[i] += delta;
				alpha[j] += delta;
				if(diff > 0)
				{
					if(alpha[j] < 0)
					{
						alpha[j] = 0;
						alpha[i] = diff;
					}
				}
				else
				{
					if(alpha[i] < 0)
					{
						alpha[i] = 0;
						alpha[j] = -diff;
					}
				}
				if(diff > C_i - C_j)
				{
					if(alpha[i] > C_i)
					{
						alpha[i] = C_i;
						alpha[j] = C_i - diff;
					}
				}
				else
				{
					if(alpha[j] > C_j)
					{
						alpha[j] = C_j;
						alpha[i] = C_j + diff;
					}
				}
			}
			// update alpha status and w
			if(fabs(alpha[i]-old_alpha_i) > 1e-16)
			{
				success_pair++;
				sparse_operator::axpy(y[i]*(alpha[i]-old_alpha_i), prob->x[i], w);
				sparse_operator::axpy(y[j]*(alpha[j]-old_alpha_j), prob->x[j], w);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[GETI(i)]);
				alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[GETI(j)]);
			}
			else
			{
				alpha[i] = old_alpha_i;
				alpha[j] = old_alpha_j;
			}
		}
		iter++;
		duration += clock() - start;
		log_message();

	}
	summary();
	delete [] alpha_status;
}

void Solver::oneclass_random_shrink()
{
	clock_t start;
	int l = prob->l;
	int i, j;
	active_size = l;
	double G_i, G_j;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	double Gmax_old = INF, Gmin_old = -INF;
	
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		n_exchange = 0;
		Gmax = -INF;
		Gmin = INF;
		for(int index_i = 0; index_i<active_size; index_i++)
		{
			int si = rand()%active_size;
			int sj = rand()%active_size;
			
			i = index[si];
			j = index[sj];

			feature_node const * xi = prob->x[i];
			feature_node const * xj = prob->x[j];
			
			double Q_ij = 0;
			G_i = 0;
			G_j = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			
			int alp_st = alpha_status[i];
			bool shrink = false;
			if(alp_st == LOWER_BOUND)
			{
				if(-G_i < Gmin_old)
				{
					//shrink
					active_size--;
					swap(index[si], index[active_size]);
					shrink = true;
				}
				else
					Gmax = max(-G_i, Gmax);
			}
			else if (alp_st == UPPER_BOUND)
			{
				if(-G_i > Gmax_old)
				{
					//shrink
					active_size--;
					swap(index[si], index[active_size]);
					shrink = true;
				}
				else
					Gmin = min(-G_i, Gmin);
			}
			else
			{
				Gmax = max(-G_i, Gmax);
				Gmin = min(-G_i, Gmin);
			}

			alp_st = alpha_status[j];
			if(alp_st == LOWER_BOUND)
			{
				if(-G_j < Gmin_old)
				{
					//shrink
					active_size--;
					swap(index[sj], index[active_size]);
					shrink = true;
				}
				Gmax = max(-G_j, Gmax);
			}
			else if (alp_st == UPPER_BOUND)
			{
				if(-G_j > Gmax_old)
				{
					//shrink
					active_size--;
					swap(index[sj], index[active_size]);
					shrink = true;
				}
				Gmin = min(-G_j, Gmin);
			}
			else
			{
				Gmax = max(-G_j, Gmax);
				Gmin = min(-G_j, Gmin);
			}
			if(shrink)
				continue;
			xi = prob->x[i];
			xj = prob->x[j];

			double C_i = upper_bound[2];
			double C_j = upper_bound[2];

			double old_alpha_i = alpha[i];
			double old_alpha_j = alpha[j];

			double quad_coef = QD[i] + QD[j] - 2*Q_ij;
			if(quad_coef <= 0)
				quad_coef = 1e-12;
			double delta = (G_i-G_j)/quad_coef;
			double sum = alpha[i] +alpha[j];
			alpha[i] -= delta;
			alpha[j] += delta;
			if(sum > C_i)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = sum -C_i;
				}
			}
			else
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = sum -C_j;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = sum;
				}
			}
			// update alpha status and w
			if(fabs(alpha[i]-old_alpha_i) > 1e-16)
			{
				success_pair++;
				if(fabs(alpha[i]-old_alpha_i) == upper_bound[2])
					n_exchange++;
				sparse_operator::axpy(alpha[i]-old_alpha_i, prob->x[i], w);
				sparse_operator::axpy(alpha[j]-old_alpha_j, prob->x[j], w);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[2]);
				alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[2]);
			}
			else
			{
				alpha[i] = old_alpha_i;
				alpha[j] = old_alpha_j;
			}
		}
		iter++;
		duration += clock() - start;
		log_message();
		if(Gmax-Gmin< eps)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				Gmax_old = INF;
				Gmin_old = -INF;
			}
		}
		else
		{
			Gmax_old = Gmax;
			Gmin_old = Gmin;
		}
	}
	summary();
}

void Solver::oneclass_random_1000()
{
	clock_t start;
	int l = prob->l;
	int i, j;
	active_size = l;
	double G_i, G_j;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		n_exchange = 0;
		for(int index_i = 0; index_i<active_size; index_i++)
		{
			i = rand()%active_size;
			j = rand()%active_size;
			
			feature_node const * xi = prob->x[i];
			feature_node const * xj = prob->x[j];
			
			double Q_ij = 0;
			G_i = 0;
			G_j = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			
			xi = prob->x[i];
			xj = prob->x[j];

			double C_i = upper_bound[2];
			double C_j = upper_bound[2];

			double old_alpha_i = alpha[i];
			double old_alpha_j = alpha[j];

			double quad_coef = QD[i] + QD[j] - 2*Q_ij;
			if(quad_coef <= 0)
				quad_coef = 1e-12;
			double delta = (G_i-G_j)/quad_coef;
			double sum = alpha[i] +alpha[j];
			alpha[i] -= delta;
			alpha[j] += delta;
			if(sum > C_i)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = sum -C_i;
				}
			}
			else
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = sum -C_j;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = sum;
				}
			}
			// update alpha status and w
			if(fabs(alpha[i]-old_alpha_i) > 1e-16)
			{
				success_pair++;
				if(fabs(alpha[i]-old_alpha_i) == upper_bound[2])
					n_exchange++;
				sparse_operator::axpy(alpha[i]-old_alpha_i, prob->x[i], w);
				sparse_operator::axpy(alpha[j]-old_alpha_j, prob->x[j], w);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[2]);
				alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[2]);
			}
			else
			{
				alpha[i] = old_alpha_i;
				alpha[j] = old_alpha_j;
			}
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
}

void Solver::oneclass_first_1000()
{
	clock_t start;
	int l = prob->l;
	int i, j;
	double G_i, G_j;
	double *G = new double[l];
	active_size = l;
	update_size = 1;
	int Gmax_index= -1, Gmin_index = -1;
	double Gmax=-INF, Gmin = INF;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		Gmax = -INF;
		Gmin = INF;
		Gmax_index =-1;
		Gmin_index =-1;
		for(i=0; i<active_size; i++)
		{
			int G_index = index[i];
			feature_node * const xi = prob->x[G_index];
			G[i] = -sparse_operator::dot(w, xi);
			if(alpha_status[G_index] != UPPER_BOUND)
			{
				if(Gmax <= G[i])
				{
					Gmax = G[i];
					Gmax_index = G_index;
				}
				
			}
			if(alpha_status[G_index] != LOWER_BOUND)
			{
				if(Gmin >= G[i])
				{
					Gmin = G[i];
					Gmin_index = G_index;
				}	
			}
			
		}

		i = Gmax_index;
		j = Gmin_index;
		G_i = -G[Gmax_index];
		G_j = -G[Gmin_index];
		feature_node const * xi = prob->x[i];
		feature_node const * xj = prob->x[j];
		double Q_ij = 0;
		while(xi->index != -1 && xj->index != -1)
		{
			if(xi->index == xj->index)
			{
				Q_ij += xi->value * xj->value;
				++xi;
				++xj;
			}
			else
			{
				if(xi->index > xj->index)
					++xj;
				else
					++xi;
	
			}
		}

		double C_i = upper_bound[2];
		double C_j = upper_bound[2];

		double old_alpha_i = alpha[i];
		double old_alpha_j = alpha[j];

		double quad_coef = QD[i] + QD[j] - 2*Q_ij;
		if(quad_coef <= 0)
			quad_coef = 1e-12;
		double delta = (G_i-G_j)/quad_coef;
		double sum = alpha[i] + alpha[j];
		alpha[i] -= delta;
		alpha[j] += delta;
		if(sum > C_i)
		{
			if(alpha[i] > C_i)
			{
				alpha[i] = C_i;
				alpha[j] = sum -C_i;
			}
		}
		else
		{
			if(alpha[j] < 0)
			{
				alpha[j] = 0;
				alpha[i] = sum;
			}
		}
		if(sum > C_j)
		{
			if(alpha[j] > C_j)
			{
				alpha[j] = C_j;
				alpha[i] = sum -C_j;
			}
		}
		else
		{
			if(alpha[i] < 0)
			{
				alpha[i] = 0;
				alpha[j] = sum;
			}
		}
		// update alpha status and w
		if(fabs(alpha[i]-old_alpha_i) > 1e-16)
		{
			success_pair++;
			sparse_operator::axpy(alpha[i]-old_alpha_i, prob->x[i], w);
			sparse_operator::axpy(alpha[j]-old_alpha_j, prob->x[j], w);
			alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[2]);
			alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[2]);
		}
		else
		{
			alpha[i] = old_alpha_i;
			alpha[j] = old_alpha_j;
		}
		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
	delete [] G;
}

void Solver::oneclass_second_1000()
{
	clock_t start;
	int l = prob->l;
	int i, j;
	active_size = l;
	update_size = 1;
	double G_i, G_j;
	double *G = new double[l];
	int Gmax_index= -1, Gmin_index = -1;
	double Gmax=-INF, smin = INF;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		Gmax = -INF;
		for(i=0; i<active_size; i++)
		{
			int G_index = index[i];
			feature_node * const xi = prob->x[G_index];
			G[i] = -sparse_operator::dot(w, xi);
			if(alpha_status[G_index] != UPPER_BOUND)
			{
				if(Gmax <= G[i])
				{
					Gmax = G[i];
					Gmax_index = G_index;
				}
			}
		}

		smin = INF;
		double Q_ij = 0;
		for(j = 0; j<active_size; j++)
		{
			int G_index = index[j];
			if(alpha_status[G_index] != LOWER_BOUND)
			{
				feature_node const * xj = prob->x[G_index];
				feature_node const * xi = prob->x[Gmax_index];
				if(G[G_index] < G[Gmax_index])
				{
					double Q_ijpair = 0;
					while(xi->index != -1 && xj->index != -1)
					{
						if(xi->index == xj->index)
						{
							Q_ijpair += xi->value * xj->value;
							++xi;
							++xj;
						}
						else
						{
							if(xi->index > xj->index)
								++xj;
							else
								++xi;
						}
					}
					xj = prob->x[G_index];
					xi = prob->x[Gmax_index];
					double b = G[Gmax_index] - G[G_index];
					double a = QD[Gmax_index] + QD[G_index] - 2*Q_ijpair;
					double obj = -b*b / (2*a);
					if(obj < smin)
					{
						smin = obj;
						Gmin_index = G_index;
						Q_ij = Q_ijpair;
					}
				}
			}
		}
		i = Gmax_index;
		j = Gmin_index;
		G_i = -G[Gmax_index];
		G_j = -G[Gmin_index];
		
		double C_i = upper_bound[2];
		double C_j = upper_bound[2];

		double old_alpha_i = alpha[i];
		double old_alpha_j = alpha[j];

		double quad_coef = QD[i] + QD[j] - 2*Q_ij;
		if(quad_coef <= 0)
			quad_coef = 1e-12;
		double delta = (G_i-G_j)/quad_coef;
		double sum = alpha[i] +alpha[j];
		alpha[i] -= delta;
		alpha[j] += delta;
		if(sum > C_i)
		{
			if(alpha[i] > C_i)
			{
				alpha[i] = C_i;
				alpha[j] = sum -C_i;
			}
		}
		else
		{
			if(alpha[j] < 0)
			{
				alpha[j] = 0;
				alpha[i] = sum;
			}
		}
		if(sum > C_j)
		{
			if(alpha[j] > C_j)
			{
				alpha[j] = C_j;
				alpha[i] = sum -C_j;
			}
		}
		else
		{
			if(alpha[i] < 0)
			{
				alpha[i] = 0;
				alpha[j] = sum;
			}
		}
		// update alpha status and w
		if(fabs(alpha[i]-old_alpha_i) > 1e-16)
		{
			success_pair++;
			sparse_operator::axpy(alpha[i]-old_alpha_i, prob->x[i], w);
			sparse_operator::axpy(alpha[j]-old_alpha_j, prob->x[j], w);
			alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[2]);
			alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[2]);
		}
		else
		{
			alpha[i] = old_alpha_i;
			alpha[j] = old_alpha_j;
		}

		iter++;
		duration += clock() - start;
		log_message();
	}
	summary();
	delete [] G;
}

void Solver::oneclass_semigd_1000()
{
	clock_t start;
	int l = prob->l;
	int i, j;
	active_size = l;
	double G_i, G_j;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	
	int *Max_order_index = new int[l];
	int *Min_order_index = new int[l];
				
	std::priority_queue<struct feature_node, std::vector<feature_node>, maxcomp> Max_order_queue;
	std::priority_queue<struct feature_node, std::vector<feature_node>, mincomp> Min_order_queue;
	
	int success_all = 0;
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		n_exchange = 0;
		//shrinking
		update_size = int(active_size*ratio_update);
		if(update_size < 1)
			update_size = 1;

		for(i=0; i<active_size; i++)
		{
			struct feature_node comp;
			feature_node const * xi = prob->x[index[i]];
			comp.index = index[i];
			comp.value = -sparse_operator::dot(w, xi);
			if(alpha_status[comp.index] != UPPER_BOUND)
			{
				if((int) Min_order_queue.size() <  update_size)
					Min_order_queue.push(comp);
				else
				{
					if(Min_order_queue.top().value < comp.value)
					{
						Min_order_queue.pop();
						Min_order_queue.push(comp);
					}
				}
			}
			if(alpha_status[comp.index] != LOWER_BOUND)
			{
				if((int) Max_order_queue.size() < update_size)
					Max_order_queue.push(comp);
				else
				{
					if(Max_order_queue.top().value > comp.value)
					{
						Max_order_queue.pop();
						Max_order_queue.push(comp);
					}
				}
			}
		}
		
		update_size = min((int)Min_order_queue.size(), (int)Max_order_queue.size());
		while((int)Max_order_queue.size() > update_size)
			Max_order_queue.pop();
		while((int)Min_order_queue.size() > update_size)
			Min_order_queue.pop();
		
		for(i=0; i<update_size; i++)
		{
			Max_order_index[update_size-1-i] = Min_order_queue.top().index;
			Min_order_index[update_size-1-i] = Max_order_queue.top().index;
			Min_order_queue.pop();
			Max_order_queue.pop();
		}
		for(int index_i = 0; index_i<update_size; index_i++)
		{
			i = Max_order_index[index_i];
			j = Min_order_index[index_i];
			
			feature_node const * xi = prob->x[i];
			feature_node const * xj = prob->x[j];
			
			double Q_ij = 0;
			G_i = 0;
			G_j = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}
			
			xi = prob->x[i];
			xj = prob->x[j];

			double C_i = upper_bound[2];
			double C_j = upper_bound[2];

			double old_alpha_i = alpha[i];
			double old_alpha_j = alpha[j];

			double quad_coef = QD[i] + QD[j] - 2*Q_ij;
			if(quad_coef <= 0)
				quad_coef = 1e-12;
			double delta = (G_i-G_j)/quad_coef;
			double sum = alpha[i] +alpha[j];
			alpha[i] -= delta;
			alpha[j] += delta;
			if(sum > C_i)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = sum -C_i;
				}
			}
			else
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = sum -C_j;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = sum;
				}
			}
			// update alpha status and w
			if(fabs(alpha[i]-old_alpha_i) > 1e-16)
			{
				success_pair++;
				if(fabs(alpha[i]-old_alpha_i) == upper_bound[2])
					n_exchange++;
				sparse_operator::axpy(alpha[i]-old_alpha_i, prob->x[i], w);
				sparse_operator::axpy(alpha[j]-old_alpha_j, prob->x[j], w);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[2]);
				alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[2]);
			}
			else
			{
				alpha[i] = old_alpha_i;
				alpha[j] = old_alpha_j;
			}
		}

		iter++;
		duration += clock() - start;
		log_message();
		success_all += success_pair;
	}
	summary();
	delete [] Max_order_index;
	delete [] Min_order_index;
}

void Solver::oneclass_semigd_shrink()
{
	clock_t start;
	int l = prob->l;
	int i, j;
	active_size = l;
	double G_i, G_j;
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	Gmax = -INF;
	Gmin = INF;
	
	int *Max_order_index = new int[l];
	int *Min_order_index = new int[l];
				
	std::priority_queue<struct feature_node, std::vector<feature_node>, maxcomp> Max_order_queue;
	std::priority_queue<struct feature_node, std::vector<feature_node>, mincomp> Min_order_queue;

	
	int success_all = 0;
	double *G = new double[l];
	while(iter < max_iter)
	{
		start = clock();
		success_pair = 0;
		n_exchange = 0;
		Gmax = -INF;
		Gmin = INF;
		for(i=0; i<active_size; i++)
		{
			int G_index = index[i];
			feature_node * const xi = prob->x[G_index];
			G[i] = -sparse_operator::dot(w, xi);
			if(alpha_status[G_index] != UPPER_BOUND)
			{
				if(Gmax <= G[i])
					Gmax = G[i];
			}
			if(alpha_status[G_index] != LOWER_BOUND)
			{
				if(Gmin >= G[i])
					Gmin = G[i];
			}
		}
		//shrinking
		if(Gmax - Gmin < eps)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				continue;
			}
		}
		else
		{
			for(i=active_size-1; i>=0; i--)
			{
				bool be_shunk = false;
				int shrinking_index = index[i];
				if(alpha_status[shrinking_index] == UPPER_BOUND)
				{
					if(G[i] > Gmax)
						be_shunk = true;
				}
				else if(alpha_status[shrinking_index] == LOWER_BOUND)
				{
					if(G[i] < Gmin)
						be_shunk = true;
				}

				if(be_shunk)
				{
					swap(index[i],index[--active_size]);
					swap(G[i],G[active_size]);
				}

			}
		}
		update_size = int(active_size*ratio_update);
		if(update_size < 1)
			update_size = 1;

		for(i=0; i<active_size; i++)
		{
			struct feature_node comp;
			comp.index = index[i];
			comp.value = G[i];
			if(alpha_status[comp.index] != UPPER_BOUND)
			{
				if((int) Min_order_queue.size() <  update_size)
					Min_order_queue.push(comp);
				else
				{
					if(Min_order_queue.top().value < comp.value)
					{
						Min_order_queue.pop();
						Min_order_queue.push(comp);
					}
				}
			}
			if(alpha_status[comp.index] != LOWER_BOUND)
			{
				if((int) Max_order_queue.size() < update_size)
					Max_order_queue.push(comp);
				else
				{
					if(Max_order_queue.top().value > comp.value)
					{
						Max_order_queue.pop();
						Max_order_queue.push(comp);
					}
				}
			}
		}
		
		update_size = min((int)Min_order_queue.size(), (int)Max_order_queue.size());
		while((int)Max_order_queue.size() > update_size)
			Max_order_queue.pop();
		while((int)Min_order_queue.size() > update_size)
			Min_order_queue.pop();
		
		for(i=0; i<update_size; i++)
		{
			Max_order_index[update_size-1-i] = Min_order_queue.top().index;
			Min_order_index[update_size-1-i] = Max_order_queue.top().index;
			Min_order_queue.pop();
			Max_order_queue.pop();
		}
		for(int index_i = 0; index_i<update_size; index_i++)
		{
			i = Max_order_index[index_i];
			j = Min_order_index[index_i];
			
			feature_node const * xi = prob->x[i];
			feature_node const * xj = prob->x[j];
			
			double Q_ij = 0;
			G_i = 0;
			G_j = 0;
			while(xi->index != -1 && xj->index != -1)
			{
				if(xi->index == xj->index)
				{
					Q_ij += xi->value * xj->value;
					G_i += xi->value * w[xi->index-1];
					G_j += xj->value * w[xj->index-1];
					++xi;
					++xj;
				}
				else
				{
					if(xi->index > xj->index)
					{
						G_j += xj->value * w[xj->index-1];
						++xj;
					}
					else
					{
						G_i += xi->value * w[xi->index-1];
						++xi;
					}
				}
			}
			while(xi->index != -1)
			{
				G_i += xi->value * w[xi->index-1];
				++xi;
			}
			while(xj->index != -1)
			{
				G_j += xj->value * w[xj->index-1];
				++xj;
			}

			double C_i = upper_bound[2];
			double C_j = upper_bound[2];

			double old_alpha_i = alpha[i];
			double old_alpha_j = alpha[j];

			double quad_coef = QD[i] + QD[j] - 2*Q_ij;
			if(quad_coef <= 0)
				quad_coef = 1e-12;
			double delta = (G_i-G_j)/quad_coef;
			double sum = alpha[i] +alpha[j];
			alpha[i] -= delta;
			alpha[j] += delta;
			if(sum > C_i)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = sum -C_i;
				}
			}
			else
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = sum -C_j;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = sum;
				}
			}
			// update alpha status and w
			if(fabs(alpha[i]-old_alpha_i) > 1e-16)
			{
				success_pair++;
				if(fabs(alpha[i]-old_alpha_i) == upper_bound[2])
					n_exchange++;
				sparse_operator::axpy(alpha[i]-old_alpha_i, prob->x[i], w);
				sparse_operator::axpy(alpha[j]-old_alpha_j, prob->x[j], w);
				double nG_i = sparse_operator::dot(w, prob->x[i]);
				double nG_j = sparse_operator::dot(w, prob->x[j]);
				info("%g %g %g %g | %g %g\n", alpha[i], alpha[j], nG_i, nG_j, G_i, G_j);
				alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[2]);
				alpha_status[j] = updateAlphaStatus(alpha[j],upper_bound[2]);
			}
			else
			{
				alpha[i] = old_alpha_i;
				alpha[j] = old_alpha_j;
			}
		}

		iter++;
		duration += clock() - start;
		log_message();
		success_all += success_pair;
	}
	summary();
	delete [] Max_order_index;
	delete [] Min_order_index;
}

static void oneclass_update(
	const problem *prob, double *w,
	const parameter *param)
{
	int solver_type = param->solver_type;
	int l = prob->l;
	double nu = param->nu;
	int i;
	double *QD = new double[l];
	double *alpha =  new double[l];
	int *alpha_status = new int[l];
	int *index = new int[l];
	// default solver_type: ONE_L2_CY_SH
	double diag[3] = {1, 0, 1};//need to check
	double upper_bound[3] = {INF, 0, INF};
	enum {LOWER_BOUND, UPPER_BOUND, FREE};

	if(param->solver_type == ONECLASS_L1_RD_1000
	|| param->solver_type == ONECLASS_L1_SEMIGD_1000
	|| param->solver_type == ONECLASS_L1_RD_SH
	|| param->solver_type == ONECLASS_L1_SECOND_1000
	|| param->solver_type == ONECLASS_L1_FIRST_1000
	|| param->solver_type == ONECLASS_L1_SEMIGD_SH)
	{
		diag[0] = 0;
		diag[2] = 0;
		upper_bound[0] = 1;
		upper_bound[2] = 1;
	}

	// Initial alpha can be set here. Note that
	// 0 <= alpha[i] <= upper_bound[GETI(i)]
	int n = (int)(nu*l); 
 	for(i=0; i<n; i++)
	{
		alpha[i] = 1;
	}
	if(n<l)
	{
		alpha[i] = nu*l-n;
	}
	for(i=n+1; i<l; i++)
	{
		alpha[i] = 0;
	}
	for(i=0; i<l; i++)
	{
		alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[2]);
		feature_node * const xi = prob->x[i];
		QD[i] = sparse_operator::nrm2_sq(xi);
		sparse_operator::axpy(alpha[i], xi, w);
		index[i] = i;
	}
	Solver solver = Solver(solver_type);
	solver.prob = prob;
	solver.w = w;
	solver.alpha = alpha;
	solver.eps = param->eps;
	solver.QD = QD;
	solver.alpha_status = alpha_status;
	solver.index = index;
	solver.diag = diag;
	solver.upper_bound = upper_bound;
	solver.ratio_update = param->r;
	solver.max_iter = param->max_iter;
	switch(solver_type)
	{
		
		case ONECLASS_L1_RD_SH:
		case ONECLASS_L2_RD_SH:
			solver.oneclass_random_shrink();
			break;
		case ONECLASS_L1_RD_1000:
		case ONECLASS_L2_RD_1000:
			solver.oneclass_random_1000();
			break;
		case ONECLASS_L1_FIRST_1000:
		case ONECLASS_L2_FIRST_1000:
			solver.oneclass_first_1000();
			break;
		
		case ONECLASS_L1_SECOND_1000:
		case ONECLASS_L2_SECOND_1000:
			solver.oneclass_second_1000();
			break;
		case ONECLASS_L1_SEMIGD_1000:
		case ONECLASS_L2_SEMIGD_1000:
			solver.oneclass_semigd_1000();
			break;
		
		case ONECLASS_L1_SEMIGD_SH:
		case ONECLASS_L2_SEMIGD_SH:
			solver.oneclass_semigd_shrink();
			break;
	}

	delete [] QD;
	delete [] alpha;
	delete [] alpha_status;
	delete [] index;
}
static void two_bias_update(
	const problem *prob, double *w, double eps,
	const double Cp, const double Cn, const parameter *param)
{
	int solver_type = param->solver_type;
	int l = prob->l;
	int w_size = prob->n;
	int i;
	double *QD = new double[l];
	double *alpha =  new double[l];
	schar *y = new schar[l];
	enum {LOWER_BOUND, UPPER_BOUND, FREE};
	// default solver_type: ONE_L2_CY_SH
	double diag[3] = {0.5/Cn, 0, 0.5/Cp};
	double upper_bound[3] = {INF, 0, INF};
	int *alpha_status = new int[l];
	int *index = new int[l];

	if( solver_type == BIAS_L1_RD_1000 ||
		solver_type == BIAS_L1_SEMIGD_1000 ||
		solver_type == BIAS_L1_SEMIGD_SH ||
		solver_type == BIAS_L1_RD_SH )
	{
		diag[0] = 0;
		diag[2] = 0;
		upper_bound[0] = Cn;
		upper_bound[2] = Cp;
	}

	for(i=0; i<l; i++)
	{
		if(prob->y[i] > 0)
		{
			y[i] = +1;
		}
		else
		{
			y[i] = -1;
		}
	}

	// Initial alpha can be set here. Note that
	// 0 <= alpha[i] <= upper_bound[GETI(i)]
	for(i=0; i<l; i++)
		alpha[i] = 0;

	for(i=0; i<w_size; i++)
		w[i] = 0;
	for(i=0; i<l; i++)
	{
		alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[GETI(i)]);
		QD[i] = diag[GETI(i)];
		feature_node * const xi = prob->x[i];
		QD[i] += sparse_operator::nrm2_sq(xi);
		sparse_operator::axpy(y[i]*alpha[i], xi, w);
		index[i] = i;
	}
	Solver solver = Solver(solver_type);
	solver.prob = prob;
	solver.w = w;
	solver.alpha = alpha;
	solver.eps = eps;
	solver.QD = QD;
	solver.y = y;
	solver.diag = diag;
	solver.index = index;
	solver.upper_bound = upper_bound;
	solver.alpha_status = alpha_status;
	solver.ratio_update = param->r;
	solver.max_iter = param->max_iter;
	switch(solver_type)
	{
		case BIAS_L1_RD_SH:
		case BIAS_L2_RD_SH:
		{
			solver.bias_random_shrink();
			break;
		}
		case BIAS_L1_RD_1000:
		case BIAS_L2_RD_1000:
		{
			solver.bias_random_1000();
			break;
		}
		case BIAS_L1_SEMIGD_1000:
		case BIAS_L2_SEMIGD_1000:
		{
			solver.bias_semigd_1000();
			break;
		}
		case BIAS_L1_SEMIGD_SH:
		case BIAS_L2_SEMIGD_SH:
		{
			solver.bias_semigd_shrink();
			break;
		}
	}

	delete [] QD;
	delete [] alpha;
	delete [] y;
}
static void onetwo_nobias_update(
	const problem *prob, double *w, double eps,
	const double Cp, const double Cn, const parameter *param)
{
	int solver_type = param->solver_type;
	int l = prob->l;
	int max_set = l;
	int w_size = prob->n;
	double *QD = new double[l];
	double *alpha =  new double[l];
	int *alpha_status = new int[l];
	int *index = new int[l];
	schar *y = new schar[l];
	int i = 0;
	
	double diag[3] = {0.5/Cn, 0, 0.5/Cp};
	double upper_bound[3] = {INF, 0, INF};
	if(solver_type == ONE_L1_CY_1000 
	 ||solver_type == ONE_L1_RD_1000 
	 ||solver_type == ONE_L1_CY_SH
	 ||solver_type == ONE_L1_RD_SH
	 ||solver_type == ONE_L1_SEMIGD_1000
	 ||solver_type == ONE_L1_SEMIGD_SH
	 ||solver_type == TWO_L1_CY_1000 
	 ||solver_type == TWO_L1_RD_1000 
	 ||solver_type == TWO_L1_SEMIGD_1000 
	 ||solver_type == TWO_L1_SEMIGD_SH 
	 ||solver_type == TWO_L1_SEMICY_1000 
	 ||solver_type == TWO_L1_SEMIRDONE_1000 
	 ||solver_type == TWO_L1_SEMIRDTWO_1000 
	 ||solver_type == TWO_L1_RD_SH
	 ||solver_type == TWO_L1_RD_SH2
	 )
	{
		diag[0] = 0;
		diag[2] = 0;
		upper_bound[0] = Cn;
		upper_bound[2] = Cp;
	}

	for(i=0; i<l; i++)
	{
		if(prob->y[i] > 0)
		{
			y[i] = +1;
		}
		else
		{
			y[i] = -1;
		}
	}

	// Initial alpha can be set here. Note that
	// 0 <= alpha[i] <= upper_bound[GETI(i)]
	for(i=0; i<w_size; i++)
		w[i] = 0;

	for(i=0; i<l; i++)
	{
		QD[i] = diag[GETI(i)];
		feature_node * const xi = prob->x[i];
		QD[i] += sparse_operator::nrm2_sq(xi);
		if(QD[i] == 0)
		{
			alpha[i] = upper_bound[GETI(i)];
			index[--max_set] = i;
		}
		else
		{
			alpha[i] = 0;
			index[i-l+max_set] = i;
		}
		sparse_operator::axpy(y[i]*alpha[i], xi, w);
		alpha_status[i] = updateAlphaStatus(alpha[i],upper_bound[GETI(i)]);
	}
	Solver solver = Solver(solver_type);
	solver.prob = prob;
	solver.w = w;
	solver.alpha = alpha;
	solver.eps = eps;
	solver.QD = QD;
	solver.y = y;
	solver.diag = diag;
	solver.max_set = max_set;
	solver.alpha_status = alpha_status;
	solver.ratio_update = param->r;
	solver.index = index;
	solver.upper_bound = upper_bound;
	solver.max_iter = param->max_iter;
	switch(solver_type)
	{
		case ONE_L1_RD_SH:
		case ONE_L2_RD_SH:
		{
			solver.one_random_shrink();
			break;
		}
		case ONE_L1_CY_1000:
		case ONE_L2_CY_1000:
		{
			solver.one_cyclic_1000();
			break;	
		}
		case ONE_L1_CY_SH:
		case ONE_L2_CY_SH:
		{
			solver.one_cyclic_shrink();
			break;
		}
		case ONE_L1_RD_1000:
		case ONE_L2_RD_1000:
		{
			solver.one_random_1000();
			break;
		}
		case ONE_L1_SEMIGD_1000:
		case ONE_L2_SEMIGD_1000:
		{
			solver.one_semigd_1000();
			break;
		}
		case ONE_L1_SEMIGD_SH:
		case ONE_L2_SEMIGD_SH:
		{
			solver.one_semigd_shrink();
			break;
		}
		case TWO_L1_RD_SH:
		case TWO_L2_RD_SH:
		{
			solver.two_random_shrink();
			break;
		}
		case TWO_L1_RD_SH2:
		case TWO_L2_RD_SH2:
		{
			solver.two_random_shrink2();
			break;
		}
		case TWO_L1_CY_1000:
		case TWO_L2_CY_1000:
		{
			solver.two_cyclic_1000();
			break;	
		}
		case TWO_L1_RD_1000:
		case TWO_L2_RD_1000:
		{
			solver.two_random_1000();
			break;
		}
		case TWO_L1_SEMIGD_1000:
		case TWO_L2_SEMIGD_1000:
		{
			solver.two_semigd_1000();
			break;
		}
		case TWO_L1_SEMIGD_SH:
		case TWO_L2_SEMIGD_SH:
		{
			solver.two_semigd_shrink();
			break;
		}
		case TWO_L1_SEMICY_1000:
		case TWO_L2_SEMICY_1000:
		{
			solver.two_semicyclic_1000();
			break;
		}
		case TWO_L1_SEMIRDONE_1000:
		case TWO_L2_SEMIRDONE_1000:
		{
			solver.two_semirandom1_1000();
			break;
		}
		case TWO_L1_SEMIRDTWO_1000:
		case TWO_L2_SEMIRDTWO_1000:
		{
			solver.two_semirandom2_1000();
			break;
		}
	}
	delete [] index;
	delete [] QD;
	delete [] alpha;
	delete [] y;
}

void Solver::one_cyclic_shrink()
{
	int l = max_set;
	int i, s;
	double C, d, G;
	clock_t start;
	// PG: projected gradient, for shrinking and stopping
	double PG;
	double PGmax_old = INF;
	double PGmin_old = -INF;
	
	active_size = l;
	start = clock();
	while (iter < max_iter)
	{
		PGmax_new = -INF;
		PGmin_new = INF;
		start = clock();
		success_pair = 0;
		for (i=0; i<active_size; i++)
		{
			int j = i+rand()%(active_size-i);
			swap(index[i], index[j]);
		}
		for (s=0; s<active_size; s++)
		{
			i = index[s];
			const schar yi = y[i];
			feature_node * const xi = prob->x[i];

			G = yi*sparse_operator::dot(w, xi)-1;

			C = upper_bound[GETI(i)];
			G += alpha[i]*diag[GETI(i)];

			PG = 0;
			if (alpha[i] == 0)
			{
				if (G > PGmax_old)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}
				else if (G < 0)
					PG = G;
			}
			else if (alpha[i] == C)
			{
				if (G < PGmin_old)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}
				else if (G > 0)
					PG = G;
			}
			else
				PG = G;

			PGmax_new = max(PGmax_new, PG);
			PGmin_new = min(PGmin_new, PG);

			if(fabs(PG) > 1.0e-12)
			{
				double alpha_old = alpha[i];
				alpha[i] = min(max(alpha[i] - G/QD[i], 0.0), C);
				d = (alpha[i] - alpha_old)*yi;
				sparse_operator::axpy(d, xi, w);
				success_pair++;
			}
		}
		iter++;
		duration += clock()- start;
		log_message();
		if(PGmax_new - PGmin_new <= eps)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				PGmax_old = INF;
				PGmin_old = -INF;
				continue;
			}
		}
		PGmax_old = PGmax_new;
		PGmin_old = PGmin_new;
		if (PGmax_old <= 0)
			PGmax_old = INF;
		if (PGmin_old >= 0)
			PGmin_old = -INF;
	}
	summary();
}

// A coordinate descent algorithm for 
// L1-loss and L2-loss epsilon-SVR dual problem
//
//  min_\beta  0.5\beta^T (Q + diag(lambda)) \beta - p \sum_{i=1}^l|\beta_i| + \sum_{i=1}^l yi\beta_i,
//    s.t.      -upper_bound_i <= \beta_i <= upper_bound_i,
// 
//  where Qij = xi^T xj and
//  D is a diagonal matrix 
//
// In L1-SVM case:
// 		upper_bound_i = C
// 		lambda_i = 0
// In L2-SVM case:
// 		upper_bound_i = INF
// 		lambda_i = 1/(2*C)
//
// Given: 
// x, y, p, C
// eps is the stopping tolerance
//
// solution will be put in w
//
// See Algorithm 4 of Ho and Lin, 2012   

#undef GETI
#define GETI(i) (0)
// To support weights for instances, use GETI(i) (i)

static void solve_l2r_l1l2_svr(
	const problem *prob, double *w, const parameter *param,
	int solver_type)
{
	int l = prob->l;
	double C = param->C;
	double p = param->p;
	int w_size = prob->n;
	double eps = param->eps;
	int i, s, iter = 0;
	int max_iter = 1000;
	int active_size = l;
	int *index = new int[l];

	double d, G, H;
	double Gmax_old = INF;
	double Gmax_new, Gnorm1_new;
	double Gnorm1_init = -1.0; // Gnorm1_init is initialized at the first iteration
	double *beta = new double[l];
	double *QD = new double[l];
	double *y = prob->y;

	// L2R_L2LOSS_SVR_DUAL
	double lambda[1], upper_bound[1];
	lambda[0] = 0.5/C;
	upper_bound[0] = INF;

	if(solver_type == L2R_L1LOSS_SVR_DUAL)
	{
		lambda[0] = 0;
		upper_bound[0] = C;
	}

	// Initial beta can be set here. Note that
	// -upper_bound <= beta[i] <= upper_bound
	for(i=0; i<l; i++)
		beta[i] = 0;

	for(i=0; i<w_size; i++)
		w[i] = 0;
	for(i=0; i<l; i++)
	{
		feature_node * const xi = prob->x[i];
		QD[i] = sparse_operator::nrm2_sq(xi);
		sparse_operator::axpy(beta[i], xi, w);

		index[i] = i;
	}


	while(iter < max_iter)
	{
		Gmax_new = 0;
		Gnorm1_new = 0;

		for(i=0; i<active_size; i++)
		{
			int j = i+rand()%(active_size-i);
			swap(index[i], index[j]);
		}

		for(s=0; s<active_size; s++)
		{
			i = index[s];
			G = -y[i] + lambda[GETI(i)]*beta[i];
			H = QD[i] + lambda[GETI(i)];

			feature_node * const xi = prob->x[i];
			G += sparse_operator::dot(w, xi);

			double Gp = G+p;
			double Gn = G-p;
			double violation = 0;
			if(beta[i] == 0)
			{
				if(Gp < 0)
					violation = -Gp;
				else if(Gn > 0)
					violation = Gn;
				else if(Gp>Gmax_old && Gn<-Gmax_old)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}
			}
			else if(beta[i] >= upper_bound[GETI(i)])
			{
				if(Gp > 0)
					violation = Gp;
				else if(Gp < -Gmax_old)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}
			}
			else if(beta[i] <= -upper_bound[GETI(i)])
			{
				if(Gn < 0)
					violation = -Gn;
				else if(Gn > Gmax_old)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}
			}
			else if(beta[i] > 0)
				violation = fabs(Gp);
			else
				violation = fabs(Gn);

			Gmax_new = max(Gmax_new, violation);
			Gnorm1_new += violation;

			// obtain Newton direction d
			if(Gp < H*beta[i])
				d = -Gp/H;
			else if(Gn > H*beta[i])
				d = -Gn/H;
			else
				d = -beta[i];

			if(fabs(d) < 1.0e-12)
				continue;

			double beta_old = beta[i];
			beta[i] = min(max(beta[i]+d, -upper_bound[GETI(i)]), upper_bound[GETI(i)]);
			d = beta[i]-beta_old;

			if(d != 0)
				sparse_operator::axpy(d, xi, w);
		}

		if(iter == 0)
			Gnorm1_init = Gnorm1_new;
		iter++;
		if(iter % 10 == 0)
			info(".");

		if(Gnorm1_new <= eps*Gnorm1_init)
		{
			if(active_size == l)
				break;
			else
			{
				active_size = l;
				info("*");
				Gmax_old = INF;
				continue;
			}
		}

		Gmax_old = Gmax_new;
	}

	info("\noptimization finished, #iter = %d\n", iter);
	if(iter >= max_iter)
		info("\nWARNING: reaching max number of iterations\nUsing -s 11 may be faster\n\n");

	// calculate objective value
	double v = 0;
	int nSV = 0;
	for(i=0; i<w_size; i++)
		v += w[i]*w[i];
	v = 0.5*v;
	for(i=0; i<l; i++)
	{
		v += p*fabs(beta[i]) - y[i]*beta[i] + 0.5*lambda[GETI(i)]*beta[i]*beta[i];
		if(beta[i] != 0)
			nSV++;
	}

	info("Objective value = %lf\n", v);
	info("nSV = %d\n",nSV);

	delete [] beta;
	delete [] QD;
	delete [] index;
}


// A coordinate descent algorithm for 
// the dual of L2-regularized logistic regression problems
//
//  min_\alpha  0.5(\alpha^T Q \alpha) + \sum \alpha_i log (\alpha_i) + (upper_bound_i - \alpha_i) log (upper_bound_i - \alpha_i),
//    s.t.      0 <= \alpha_i <= upper_bound_i,
// 
//  where Qij = yi yj xi^T xj and 
//  upper_bound_i = Cp if y_i = 1
//  upper_bound_i = Cn if y_i = -1
//
// Given: 
// x, y, Cp, Cn
// eps is the stopping tolerance
//
// solution will be put in w
//
// See Algorithm 5 of Yu et al., MLJ 2010

#undef GETI
#define GETI(i) (y[i]+1)
// To support weights for instances, use GETI(i) (i)

void solve_l2r_lr_dual(const problem *prob, double *w, double eps, double Cp, double Cn)
{
	int l = prob->l;
	int w_size = prob->n;
	int i, s, iter = 0;
	double *xTx = new double[l];
	int max_iter = 1000;
	int *index = new int[l];	
	double *alpha = new double[2*l]; // store alpha and C - alpha
	schar *y = new schar[l];
	int max_inner_iter = 100; // for inner Newton
	double innereps = 1e-2;
	double innereps_min = min(1e-8, eps);
	double upper_bound[3] = {Cn, 0, Cp};

	for(i=0; i<l; i++)
	{
		if(prob->y[i] > 0)
		{
			y[i] = +1;
		}
		else
		{
			y[i] = -1;
		}
	}
	
	// Initial alpha can be set here. Note that
	// 0 < alpha[i] < upper_bound[GETI(i)]
	// alpha[2*i] + alpha[2*i+1] = upper_bound[GETI(i)]
	for(i=0; i<l; i++)
	{
		alpha[2*i] = min(0.001*upper_bound[GETI(i)], 1e-8);
		alpha[2*i+1] = upper_bound[GETI(i)] - alpha[2*i];
	}

	for(i=0; i<w_size; i++)
		w[i] = 0;
	for(i=0; i<l; i++)
	{
		feature_node * const xi = prob->x[i];
		xTx[i] = sparse_operator::nrm2_sq(xi);
		sparse_operator::axpy(y[i]*alpha[2*i], xi, w);
		index[i] = i;
	}

	while (iter < max_iter)
	{
		for (i=0; i<l; i++)
		{
			int j = i+rand()%(l-i);
			swap(index[i], index[j]);
		}
		int newton_iter = 0;
		double Gmax = 0;
		for (s=0; s<l; s++)
		{
			i = index[s];
			const schar yi = y[i];
			double C = upper_bound[GETI(i)];
			double ywTx = 0, xisq = xTx[i];
			feature_node * const xi = prob->x[i];
			ywTx = yi*sparse_operator::dot(w, xi);
			double a = xisq, b = ywTx;

			// Decide to minimize g_1(z) or g_2(z)
			int ind1 = 2*i, ind2 = 2*i+1, sign = 1;
			if(0.5*a*(alpha[ind2]-alpha[ind1])+b < 0)
			{
				ind1 = 2*i+1;
				ind2 = 2*i;
				sign = -1;
			}

			//  g_t(z) = z*log(z) + (C-z)*log(C-z) + 0.5a(z-alpha_old)^2 + sign*b(z-alpha_old)
			double alpha_old = alpha[ind1];
			double z = alpha_old;
			if(C - z < 0.5 * C)
				z = 0.1*z;
			double gp = a*(z-alpha_old)+sign*b+log(z/(C-z));
			Gmax = max(Gmax, fabs(gp));

			// Newton method on the sub-problem
			const double eta = 0.1; // xi in the paper
			int inner_iter = 0;
			while (inner_iter <= max_inner_iter)
			{
				if(fabs(gp) < innereps)
					break;
				double gpp = a + C/(C-z)/z;
				double tmpz = z - gp/gpp;
				if(tmpz <= 0)
					z *= eta;
				else // tmpz in (0, C)
					z = tmpz;
				gp = a*(z-alpha_old)+sign*b+log(z/(C-z));
				newton_iter++;
				inner_iter++;
			}

			if(inner_iter > 0) // update w
			{
				alpha[ind1] = z;
				alpha[ind2] = C-z;
				sparse_operator::axpy(sign*(z-alpha_old)*yi, xi, w);
			}
		}

		iter++;
		if(iter % 10 == 0)
			info(".");

		if(Gmax < eps)
			break;

		if(newton_iter <= l/10)
			innereps = max(innereps_min, 0.1*innereps);

	}

	info("\noptimization finished, #iter = %d\n",iter);
	if (iter >= max_iter)
		info("\nWARNING: reaching max number of iterations\nUsing -s 0 may be faster (also see FAQ)\n\n");

	// calculate objective value

	double v = 0;
	for(i=0; i<w_size; i++)
		v += w[i] * w[i];
	v *= 0.5;
	for(i=0; i<l; i++)
		v += alpha[2*i] * log(alpha[2*i]) + alpha[2*i+1] * log(alpha[2*i+1])
			- upper_bound[GETI(i)] * log(upper_bound[GETI(i)]);
	info("Objective value = %lf\n", v);

	delete [] xTx;
	delete [] alpha;
	delete [] y;
	delete [] index;
}

// A coordinate descent algorithm for 
// L1-regularized L2-loss support vector classification
//
//  min_w \sum |wj| + C \sum max(0, 1-yi w^T xi)^2,
//
// Given: 
// x, y, Cp, Cn
// eps is the stopping tolerance
//
// solution will be put in w
//
// See Yuan et al. (2010) and appendix of LIBLINEAR paper, Fan et al. (2008)

#undef GETI
#define GETI(i) (y[i]+1)
// To support weights for instances, use GETI(i) (i)

static void solve_l1r_l2_svc(
	problem *prob_col, double *w, double eps,
	double Cp, double Cn)
{
	int l = prob_col->l;
	int w_size = prob_col->n;
	int j, s, iter = 0;
	int max_iter = 1000;
	int active_size = w_size;
	int max_num_linesearch = 20;

	double sigma = 0.01;
	double d, G_loss, G, H;
	double Gmax_old = INF;
	double Gmax_new, Gnorm1_new;
	double Gnorm1_init = -1.0; // Gnorm1_init is initialized at the first iteration
	double d_old, d_diff;
	double loss_old, loss_new;
	double appxcond, cond;

	int *index = new int[w_size];
	schar *y = new schar[l];
	double *b = new double[l]; // b = 1-ywTx
	double *xj_sq = new double[w_size];
	feature_node *x;

	double C[3] = {Cn,0,Cp};

	// Initial w can be set here.
	for(j=0; j<w_size; j++)
		w[j] = 0;

	for(j=0; j<l; j++)
	{
		b[j] = 1;
		if(prob_col->y[j] > 0)
			y[j] = 1;
		else
			y[j] = -1;
	}
	for(j=0; j<w_size; j++)
	{
		index[j] = j;
		xj_sq[j] = 0;
		x = prob_col->x[j];
		while(x->index != -1)
		{
			int ind = x->index-1;
			x->value *= y[ind]; // x->value stores yi*xij
			double val = x->value;
			b[ind] -= w[j]*val;
			xj_sq[j] += C[GETI(ind)]*val*val;
			x++;
		}
	}

	while(iter < max_iter)
	{
		Gmax_new = 0;
		Gnorm1_new = 0;

		for(j=0; j<active_size; j++)
		{
			int i = j+rand()%(active_size-j);
			swap(index[i], index[j]);
		}

		for(s=0; s<active_size; s++)
		{
			j = index[s];
			G_loss = 0;
			H = 0;

			x = prob_col->x[j];
			while(x->index != -1)
			{
				int ind = x->index-1;
				if(b[ind] > 0)
				{
					double val = x->value;
					double tmp = C[GETI(ind)]*val;
					G_loss -= tmp*b[ind];
					H += tmp*val;
				}
				x++;
			}
			G_loss *= 2;

			G = G_loss;
			H *= 2;
			H = max(H, 1e-12);

			double Gp = G+1;
			double Gn = G-1;
			double violation = 0;
			if(w[j] == 0)
			{
				if(Gp < 0)
					violation = -Gp;
				else if(Gn > 0)
					violation = Gn;
				else if(Gp>Gmax_old/l && Gn<-Gmax_old/l)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}
			}
			else if(w[j] > 0)
				violation = fabs(Gp);
			else
				violation = fabs(Gn);

			Gmax_new = max(Gmax_new, violation);
			Gnorm1_new += violation;

			// obtain Newton direction d
			if(Gp < H*w[j])
				d = -Gp/H;
			else if(Gn > H*w[j])
				d = -Gn/H;
			else
				d = -w[j];

			if(fabs(d) < 1.0e-12)
				continue;

			double delta = fabs(w[j]+d)-fabs(w[j]) + G*d;
			d_old = 0;
			int num_linesearch;
			for(num_linesearch=0; num_linesearch < max_num_linesearch; num_linesearch++)
			{
				d_diff = d_old - d;
				cond = fabs(w[j]+d)-fabs(w[j]) - sigma*delta;

				appxcond = xj_sq[j]*d*d + G_loss*d + cond;
				if(appxcond <= 0)
				{
					x = prob_col->x[j];
					sparse_operator::axpy(d_diff, x, b);
					break;
				}

				if(num_linesearch == 0)
				{
					loss_old = 0;
					loss_new = 0;
					x = prob_col->x[j];
					while(x->index != -1)
					{
						int ind = x->index-1;
						if(b[ind] > 0)
							loss_old += C[GETI(ind)]*b[ind]*b[ind];
						double b_new = b[ind] + d_diff*x->value;
						b[ind] = b_new;
						if(b_new > 0)
							loss_new += C[GETI(ind)]*b_new*b_new;
						x++;
					}
				}
				else
				{
					loss_new = 0;
					x = prob_col->x[j];
					while(x->index != -1)
					{
						int ind = x->index-1;
						double b_new = b[ind] + d_diff*x->value;
						b[ind] = b_new;
						if(b_new > 0)
							loss_new += C[GETI(ind)]*b_new*b_new;
						x++;
					}
				}

				cond = cond + loss_new - loss_old;
				if(cond <= 0)
					break;
				else
				{
					d_old = d;
					d *= 0.5;
					delta *= 0.5;
				}
			}

			w[j] += d;

			// recompute b[] if line search takes too many steps
			if(num_linesearch >= max_num_linesearch)
			{
				info("#");
				for(int i=0; i<l; i++)
					b[i] = 1;

				for(int i=0; i<w_size; i++)
				{
					if(w[i]==0) continue;
					x = prob_col->x[i];
					sparse_operator::axpy(-w[i], x, b);
				}
			}
		}

		if(iter == 0)
			Gnorm1_init = Gnorm1_new;
		iter++;
		if(iter % 10 == 0)
			info(".");

		if(Gnorm1_new <= eps*Gnorm1_init)
		{
			if(active_size == w_size)
				break;
			else
			{
				active_size = w_size;
				info("*");
				Gmax_old = INF;
				continue;
			}
		}

		Gmax_old = Gmax_new;
	}

	info("\noptimization finished, #iter = %d\n", iter);
	if(iter >= max_iter)
		info("\nWARNING: reaching max number of iterations\n");

	// calculate objective value

	double v = 0;
	int nnz = 0;
	for(j=0; j<w_size; j++)
	{
		x = prob_col->x[j];
		while(x->index != -1)
		{
			x->value *= prob_col->y[x->index-1]; // restore x->value
			x++;
		}
		if(w[j] != 0)
		{
			v += fabs(w[j]);
			nnz++;
		}
	}
	for(j=0; j<l; j++)
		if(b[j] > 0)
			v += C[GETI(j)]*b[j]*b[j];

	info("Objective value = %lf\n", v);
	info("#nonzeros/#features = %d/%d\n", nnz, w_size);

	delete [] index;
	delete [] y;
	delete [] b;
	delete [] xj_sq;
}

// A coordinate descent algorithm for 
// L1-regularized logistic regression problems
//
//  min_w \sum |wj| + C \sum log(1+exp(-yi w^T xi)),
//
// Given: 
// x, y, Cp, Cn
// eps is the stopping tolerance
//
// solution will be put in w
//
// See Yuan et al. (2011) and appendix of LIBLINEAR paper, Fan et al. (2008)

#undef GETI
#define GETI(i) (y[i]+1)
// To support weights for instances, use GETI(i) (i)

static void solve_l1r_lr(
	const problem *prob_col, double *w, double eps,
	double Cp, double Cn)
{
	int l = prob_col->l;
	int w_size = prob_col->n;
	int j, s, newton_iter=0, iter=0;
	int max_newton_iter = 100;
	int max_iter = 1000;
	int max_num_linesearch = 20;
	int active_size;
	int QP_active_size;

	double nu = 1e-12;
	double inner_eps = 1;
	double sigma = 0.01;
	double w_norm, w_norm_new;
	double z, G, H;
	double Gnorm1_init = -1.0; // Gnorm1_init is initialized at the first iteration
	double Gmax_old = INF;
	double Gmax_new, Gnorm1_new;
	double QP_Gmax_old = INF;
	double QP_Gmax_new, QP_Gnorm1_new;
	double delta, negsum_xTd, cond;

	int *index = new int[w_size];
	schar *y = new schar[l];
	double *Hdiag = new double[w_size];
	double *Grad = new double[w_size];
	double *wpd = new double[w_size];
	double *xjneg_sum = new double[w_size];
	double *xTd = new double[l];
	double *exp_wTx = new double[l];
	double *exp_wTx_new = new double[l];
	double *tau = new double[l];
	double *D = new double[l];
	feature_node *x;

	double C[3] = {Cn,0,Cp};

	// Initial w can be set here.
	for(j=0; j<w_size; j++)
		w[j] = 0;

	for(j=0; j<l; j++)
	{
		if(prob_col->y[j] > 0)
			y[j] = 1;
		else
			y[j] = -1;

		exp_wTx[j] = 0;
	}

	w_norm = 0;
	for(j=0; j<w_size; j++)
	{
		w_norm += fabs(w[j]);
		wpd[j] = w[j];
		index[j] = j;
		xjneg_sum[j] = 0;
		x = prob_col->x[j];
		while(x->index != -1)
		{
			int ind = x->index-1;
			double val = x->value;
			exp_wTx[ind] += w[j]*val;
			if(y[ind] == -1)
				xjneg_sum[j] += C[GETI(ind)]*val;
			x++;
		}
	}
	for(j=0; j<l; j++)
	{
		exp_wTx[j] = exp(exp_wTx[j]);
		double tau_tmp = 1/(1+exp_wTx[j]);
		tau[j] = C[GETI(j)]*tau_tmp;
		D[j] = C[GETI(j)]*exp_wTx[j]*tau_tmp*tau_tmp;
	}

	while(newton_iter < max_newton_iter)
	{
		Gmax_new = 0;
		Gnorm1_new = 0;
		active_size = w_size;

		for(s=0; s<active_size; s++)
		{
			j = index[s];
			Hdiag[j] = nu;
			Grad[j] = 0;

			double tmp = 0;
			x = prob_col->x[j];
			while(x->index != -1)
			{
				int ind = x->index-1;
				Hdiag[j] += x->value*x->value*D[ind];
				tmp += x->value*tau[ind];
				x++;
			}
			Grad[j] = -tmp + xjneg_sum[j];

			double Gp = Grad[j]+1;
			double Gn = Grad[j]-1;
			double violation = 0;
			if(w[j] == 0)
			{
				if(Gp < 0)
					violation = -Gp;
				else if(Gn > 0)
					violation = Gn;
				//outer-level shrinking
				else if(Gp>Gmax_old/l && Gn<-Gmax_old/l)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}
			}
			else if(w[j] > 0)
				violation = fabs(Gp);
			else
				violation = fabs(Gn);

			Gmax_new = max(Gmax_new, violation);
			Gnorm1_new += violation;
		}

		if(newton_iter == 0)
			Gnorm1_init = Gnorm1_new;

		if(Gnorm1_new <= eps*Gnorm1_init)
			break;

		iter = 0;
		QP_Gmax_old = INF;
		QP_active_size = active_size;

		for(int i=0; i<l; i++)
			xTd[i] = 0;

		// optimize QP over wpd
		while(iter < max_iter)
		{
			QP_Gmax_new = 0;
			QP_Gnorm1_new = 0;

			for(j=0; j<QP_active_size; j++)
			{
				int i = j+rand()%(QP_active_size-j);
				swap(index[i], index[j]);
			}

			for(s=0; s<QP_active_size; s++)
			{
				j = index[s];
				H = Hdiag[j];

				x = prob_col->x[j];
				G = Grad[j] + (wpd[j]-w[j])*nu;
				while(x->index != -1)
				{
					int ind = x->index-1;
					G += x->value*D[ind]*xTd[ind];
					x++;
				}

				double Gp = G+1;
				double Gn = G-1;
				double violation = 0;
				if(wpd[j] == 0)
				{
					if(Gp < 0)
						violation = -Gp;
					else if(Gn > 0)
						violation = Gn;
					//inner-level shrinking
					else if(Gp>QP_Gmax_old/l && Gn<-QP_Gmax_old/l)
					{
						QP_active_size--;
						swap(index[s], index[QP_active_size]);
						s--;
						continue;
					}
				}
				else if(wpd[j] > 0)
					violation = fabs(Gp);
				else
					violation = fabs(Gn);

				QP_Gmax_new = max(QP_Gmax_new, violation);
				QP_Gnorm1_new += violation;

				// obtain solution of one-variable problem
				if(Gp < H*wpd[j])
					z = -Gp/H;
				else if(Gn > H*wpd[j])
					z = -Gn/H;
				else
					z = -wpd[j];

				if(fabs(z) < 1.0e-12)
					continue;
				z = min(max(z,-10.0),10.0);

				wpd[j] += z;

				x = prob_col->x[j];
				sparse_operator::axpy(z, x, xTd);
			}

			iter++;

			if(QP_Gnorm1_new <= inner_eps*Gnorm1_init)
			{
				//inner stopping
				if(QP_active_size == active_size)
					break;
				//active set reactivation
				else
				{
					QP_active_size = active_size;
					QP_Gmax_old = INF;
					continue;
				}
			}

			QP_Gmax_old = QP_Gmax_new;
		}

		if(iter >= max_iter)
			info("WARNING: reaching max number of inner iterations\n");

		delta = 0;
		w_norm_new = 0;
		for(j=0; j<w_size; j++)
		{
			delta += Grad[j]*(wpd[j]-w[j]);
			if(wpd[j] != 0)
				w_norm_new += fabs(wpd[j]);
		}
		delta += (w_norm_new-w_norm);

		negsum_xTd = 0;
		for(int i=0; i<l; i++)
			if(y[i] == -1)
				negsum_xTd += C[GETI(i)]*xTd[i];

		int num_linesearch;
		for(num_linesearch=0; num_linesearch < max_num_linesearch; num_linesearch++)
		{
			cond = w_norm_new - w_norm + negsum_xTd - sigma*delta;

			for(int i=0; i<l; i++)
			{
				double exp_xTd = exp(xTd[i]);
				exp_wTx_new[i] = exp_wTx[i]*exp_xTd;
				cond += C[GETI(i)]*log((1+exp_wTx_new[i])/(exp_xTd+exp_wTx_new[i]));
			}

			if(cond <= 0)
			{
				w_norm = w_norm_new;
				for(j=0; j<w_size; j++)
					w[j] = wpd[j];
				for(int i=0; i<l; i++)
				{
					exp_wTx[i] = exp_wTx_new[i];
					double tau_tmp = 1/(1+exp_wTx[i]);
					tau[i] = C[GETI(i)]*tau_tmp;
					D[i] = C[GETI(i)]*exp_wTx[i]*tau_tmp*tau_tmp;
				}
				break;
			}
			else
			{
				w_norm_new = 0;
				for(j=0; j<w_size; j++)
				{
					wpd[j] = (w[j]+wpd[j])*0.5;
					if(wpd[j] != 0)
						w_norm_new += fabs(wpd[j]);
				}
				delta *= 0.5;
				negsum_xTd *= 0.5;
				for(int i=0; i<l; i++)
					xTd[i] *= 0.5;
			}
		}

		// Recompute some info due to too many line search steps
		if(num_linesearch >= max_num_linesearch)
		{
			for(int i=0; i<l; i++)
				exp_wTx[i] = 0;

			for(int i=0; i<w_size; i++)
			{
				if(w[i]==0) continue;
				x = prob_col->x[i];
				sparse_operator::axpy(w[i], x, exp_wTx);
			}

			for(int i=0; i<l; i++)
				exp_wTx[i] = exp(exp_wTx[i]);
		}

		if(iter == 1)
			inner_eps *= 0.25;

		newton_iter++;
		Gmax_old = Gmax_new;

		info("iter %3d  #CD cycles %d\n", newton_iter, iter);
	}

	info("=========================\n");
	info("optimization finished, #iter = %d\n", newton_iter);
	if(newton_iter >= max_newton_iter)
		info("WARNING: reaching max number of iterations\n");

	// calculate objective value

	double v = 0;
	int nnz = 0;
	for(j=0; j<w_size; j++)
		if(w[j] != 0)
		{
			v += fabs(w[j]);
			nnz++;
		}
	for(j=0; j<l; j++)
		if(y[j] == 1)
			v += C[GETI(j)]*log(1+1/exp_wTx[j]);
		else
			v += C[GETI(j)]*log(1+exp_wTx[j]);

	info("Objective value = %lf\n", v);
	info("#nonzeros/#features = %d/%d\n", nnz, w_size);

	delete [] index;
	delete [] y;
	delete [] Hdiag;
	delete [] Grad;
	delete [] wpd;
	delete [] xjneg_sum;
	delete [] xTd;
	delete [] exp_wTx;
	delete [] exp_wTx_new;
	delete [] tau;
	delete [] D;
}

// transpose matrix X from row format to column format
static void transpose(const problem *prob, feature_node **x_space_ret, problem *prob_col)
{
	int i;
	int l = prob->l;
	int n = prob->n;
	size_t nnz = 0;
	size_t *col_ptr = new size_t [n+1];
	feature_node *x_space;
	prob_col->l = l;
	prob_col->n = n;
	prob_col->y = new double[l];
	prob_col->x = new feature_node*[n];

	for(i=0; i<l; i++)
		prob_col->y[i] = prob->y[i];

	for(i=0; i<n+1; i++)
		col_ptr[i] = 0;
	for(i=0; i<l; i++)
	{
		feature_node *x = prob->x[i];
		while(x->index != -1)
		{
			nnz++;
			col_ptr[x->index]++;
			x++;
		}
	}
	for(i=1; i<n+1; i++)
		col_ptr[i] += col_ptr[i-1] + 1;

	x_space = new feature_node[nnz+n];
	for(i=0; i<n; i++)
		prob_col->x[i] = &x_space[col_ptr[i]];

	for(i=0; i<l; i++)
	{
		feature_node *x = prob->x[i];
		while(x->index != -1)
		{
			int ind = x->index-1;
			x_space[col_ptr[ind]].index = i+1; // starts from 1
			x_space[col_ptr[ind]].value = x->value;
			col_ptr[ind]++;
			x++;
		}
	}
	for(i=0; i<n; i++)
		x_space[col_ptr[i]].index = -1;

	*x_space_ret = x_space;

	delete [] col_ptr;
}

// label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
static void group_classes(const problem *prob, int *nr_class_ret, int **label_ret, int **start_ret, int **count_ret, int *perm)
{
	int l = prob->l;
	int max_nr_class = 16;
	int nr_class = 0;
	int *label = Malloc(int,max_nr_class);
	int *count = Malloc(int,max_nr_class);
	int *data_label = Malloc(int,l);
	int i;

	for(i=0;i<l;i++)
	{
		int this_label = (int)prob->y[i];
		int j;
		for(j=0;j<nr_class;j++)
		{
			if(this_label == label[j])
			{
				++count[j];
				break;
			}
		}
		data_label[i] = j;
		if(j == nr_class)
		{
			if(nr_class == max_nr_class)
			{
				max_nr_class *= 2;
				label = (int *)realloc(label,max_nr_class*sizeof(int));
				count = (int *)realloc(count,max_nr_class*sizeof(int));
			}
			label[nr_class] = this_label;
			count[nr_class] = 1;
			++nr_class;
		}
	}

	//
	// Labels are ordered by their first occurrence in the training set. 
	// However, for two-class sets with -1/+1 labels and -1 appears first, 
	// we swap labels to ensure that internally the binary SVM has positive data corresponding to the +1 instances.
	//
	if (nr_class == 2 && label[0] == -1 && label[1] == 1)
	{
		swap(label[0],label[1]);
		swap(count[0],count[1]);
		for(i=0;i<l;i++)
		{
			if(data_label[i] == 0)
				data_label[i] = 1;
			else
				data_label[i] = 0;
		}
	}

	int *start = Malloc(int,nr_class);
	start[0] = 0;
	for(i=1;i<nr_class;i++)
		start[i] = start[i-1]+count[i-1];
	for(i=0;i<l;i++)
	{
		perm[start[data_label[i]]] = i;
		++start[data_label[i]];
	}
	start[0] = 0;
	for(i=1;i<nr_class;i++)
		start[i] = start[i-1]+count[i-1];

	*nr_class_ret = nr_class;
	*label_ret = label;
	*start_ret = start;
	*count_ret = count;
	free(data_label);
}

static void train_one(const problem *prob, const parameter *param, double *w, double Cp, double Cn)
{
	//inner and outer tolerances for TRON
	double eps = param->eps;
	double eps_cg = 0.1;
	if(param->init_sol != NULL)
		eps_cg = 0.5;

	int pos = 0;
	int neg = 0;
	for(int i=0;i<prob->l;i++)
		if(prob->y[i] > 0)
			pos++;
	neg = prob->l - pos;
	double primal_solver_tol = eps*max(min(pos,neg), 1)/prob->l;
	function *fun_obj=NULL;
	switch(param->solver_type)
	{
		case L2R_LR:
		{
			double *C = new double[prob->l];
			for(int i = 0; i < prob->l; i++)
			{
				if(prob->y[i] > 0)
					C[i] = Cp;
				else
					C[i] = Cn;
			}
			fun_obj=new l2r_lr_fun(prob, C);
			TRON tron_obj(fun_obj, primal_solver_tol, eps_cg);
			tron_obj.set_print_string(liblinear_print_string);
			tron_obj.tron(w);
			delete fun_obj;
			delete[] C;
			break;
		}
		case L2R_L2LOSS_SVC:
		{
			double *C = new double[prob->l];
			for(int i = 0; i < prob->l; i++)
			{
				if(prob->y[i] > 0)
					C[i] = Cp;
				else
					C[i] = Cn;
			}
			fun_obj=new l2r_l2_svc_fun(prob, C);
			TRON tron_obj(fun_obj, primal_solver_tol, eps_cg);
			tron_obj.set_print_string(liblinear_print_string);
			tron_obj.tron(w);
			info("\n%.15lf\n",fun_obj->fun(w));
			delete fun_obj;
			delete[] C;
			break;
		}
		case L1R_L2LOSS_SVC:
		{
			problem prob_col;
			feature_node *x_space = NULL;
			transpose(prob, &x_space ,&prob_col);
			solve_l1r_l2_svc(&prob_col, w, primal_solver_tol, Cp, Cn);
			delete [] prob_col.y;
			delete [] prob_col.x;
			delete [] x_space;
			break;
		}
		case L1R_LR:
		{
			problem prob_col;
			feature_node *x_space = NULL;
			transpose(prob, &x_space ,&prob_col);
			solve_l1r_lr(&prob_col, w, primal_solver_tol, Cp, Cn);
			delete [] prob_col.y;
			delete [] prob_col.x;
			delete [] x_space;
			break;
		}
		case L2R_LR_DUAL:
			solve_l2r_lr_dual(prob, w, eps, Cp, Cn);
			break;
		case L2R_L2LOSS_SVR:
		{
			double *C = new double[prob->l];
			for(int i = 0; i < prob->l; i++)
				C[i] = param->C;

			fun_obj=new l2r_l2_svr_fun(prob, C, param->p);
			TRON tron_obj(fun_obj, param->eps);
			tron_obj.set_print_string(liblinear_print_string);
			tron_obj.tron(w);
			delete fun_obj;
			delete[] C;
			break;

		}
		case L2R_L1LOSS_SVR_DUAL:
			solve_l2r_l1l2_svr(prob, w, param, L2R_L1LOSS_SVR_DUAL);
			break;
		case L2R_L2LOSS_SVR_DUAL:
			solve_l2r_l1l2_svr(prob, w, param, L2R_L2LOSS_SVR_DUAL);
			break;
		case ONECLASS_L1_RD_1000:
		case ONECLASS_L1_SEMIGD_1000:
		case ONECLASS_L1_RD_SH:
		case ONECLASS_L1_SEMIGD_SH:
		case ONECLASS_L1_FIRST_1000:
		case ONECLASS_L1_SECOND_1000:
		case ONECLASS_L2_SECOND_1000:
		case ONECLASS_L2_FIRST_1000:
		case ONECLASS_L2_RD_1000:
		case ONECLASS_L2_SEMIGD_1000:
		case ONECLASS_L2_RD_SH:
		case ONECLASS_L2_SEMIGD_SH:
			oneclass_update(prob, w, param);
			break;
		case ONE_L2_CY_SH: //ONE_L1_CY_SH
		case ONE_L1_CY_SH: //ONE_L2_CY_SH
		case ONE_L1_CY_1000:
		case ONE_L2_CY_1000:
		case ONE_L1_RD_1000:
		case ONE_L2_RD_1000:
		case ONE_L1_RD_SH:
		case ONE_L2_RD_SH:
		case ONE_L1_SEMIGD_1000:
		case ONE_L1_SEMIGD_SH:
		case ONE_L2_SEMIGD_1000:
		case ONE_L2_SEMIGD_SH:
		case TWO_L1_CY_1000:
		case TWO_L2_CY_1000:
		case TWO_L1_SEMICY_1000:
		case TWO_L2_SEMICY_1000:
		case TWO_L1_SEMIRDONE_1000:
		case TWO_L2_SEMIRDONE_1000:
		case TWO_L1_SEMIRDTWO_1000:
		case TWO_L2_SEMIRDTWO_1000:
		case TWO_L1_RD_1000:
		case TWO_L2_RD_1000:
		case TWO_L1_SEMIGD_1000:
		case TWO_L2_SEMIGD_1000:
		case TWO_L1_SEMIGD_SH:
		case TWO_L2_SEMIGD_SH:
		case TWO_L1_RD_SH:
		case TWO_L2_RD_SH:
		case TWO_L1_RD_SH2:
		case TWO_L2_RD_SH2:
			onetwo_nobias_update(prob, w, eps, Cp, Cn, param);
			break;
		case BIAS_L1_RD_1000:
		case BIAS_L2_RD_1000:
		case BIAS_L1_SEMIGD_1000:
		case BIAS_L2_SEMIGD_1000:
		case BIAS_L1_SEMIGD_SH:
		case BIAS_L2_SEMIGD_SH:
		case BIAS_L1_RD_SH:
		case BIAS_L2_RD_SH:
			two_bias_update(prob, w, eps, Cp, Cn, param);
			break;
		default:
			fprintf(stderr, "ERROR: unknown solver_type\n");
			break;
	}
}

// Calculate the initial C for parameter selection
static double calc_start_C(const problem *prob, const parameter *param)
{
	int i;
	double xTx,max_xTx;
	max_xTx = 0;
	for(i=0; i<prob->l; i++)
	{
		xTx = 0;
		feature_node *xi=prob->x[i];
		while(xi->index != -1)
		{
			double val = xi->value;
			xTx += val*val;
			xi++;
		}
		if(xTx > max_xTx)
			max_xTx = xTx;
	}

	double min_C = 1.0;
	if(param->solver_type == L2R_LR)
		min_C = 1.0 / (prob->l * max_xTx);
	else if(param->solver_type == L2R_L2LOSS_SVC)
		min_C = 1.0 / (2 * prob->l * max_xTx);

	return pow( 2, floor(log(min_C) / log(2.0)) );
}


//
// Interface functions
//
model* train(const problem *prob, const parameter *param)
{
	int i,j;
	int l = prob->l;
	int n = prob->n;
	int w_size = prob->n;
	model *model_ = Malloc(model,1);

	if(prob->bias>=0)
		model_->nr_feature=n-1;
	else
		model_->nr_feature=n;
	model_->param = *param;
	model_->bias = prob->bias;

	if(check_regression_model(model_) 
	|| param->solver_type == ONECLASS_L1_RD_1000
	|| param->solver_type == ONECLASS_L1_SEMIGD_1000
	|| param->solver_type == ONECLASS_L1_RD_SH
	|| param->solver_type == ONECLASS_L1_SEMIGD_SH
	|| param->solver_type == ONECLASS_L1_FIRST_1000
	|| param->solver_type == ONECLASS_L1_SECOND_1000
	|| param->solver_type == ONECLASS_L2_SECOND_1000
	|| param->solver_type == ONECLASS_L2_FIRST_1000
	|| param->solver_type == ONECLASS_L2_RD_1000
	|| param->solver_type == ONECLASS_L2_SEMIGD_1000
	|| param->solver_type == ONECLASS_L2_RD_SH
	|| param->solver_type == ONECLASS_L2_SEMIGD_SH)
	{
		if(param->solver_type == ONECLASS_L1_RD_1000
		|| param->solver_type == ONECLASS_L1_SEMIGD_1000
		|| param->solver_type == ONECLASS_L1_RD_SH
		|| param->solver_type == ONECLASS_L1_SEMIGD_SH
		|| param->solver_type == ONECLASS_L1_FIRST_1000
		|| param->solver_type == ONECLASS_L1_SECOND_1000
		|| param->solver_type == ONECLASS_L2_SECOND_1000
		|| param->solver_type == ONECLASS_L2_FIRST_1000
		|| param->solver_type == ONECLASS_L2_RD_1000
		|| param->solver_type == ONECLASS_L2_SEMIGD_1000
		|| param->solver_type == ONECLASS_L2_RD_SH
		|| param->solver_type == ONECLASS_L2_SEMIGD_SH)
		{
			w_size = w_size + 1;
			model_->bias = 1;
		}
		model_->w = Malloc(double, w_size);
		for(i=0; i<w_size; i++)
			model_->w[i] = 0;
		model_->nr_class = 2;
		model_->label = NULL;
		train_one(prob, param, model_->w, 0, 0);
	}
	else
	{
		int nr_class;
		int *label = NULL;
		int *start = NULL;
		int *count = NULL;
		int *perm = Malloc(int,l);

		// group training data of the same class
		group_classes(prob,&nr_class,&label,&start,&count,perm);

		model_->nr_class=nr_class;
		model_->label = Malloc(int,nr_class);
		for(i=0;i<nr_class;i++)
			model_->label[i] = label[i];

		// calculate weighted C
		double *weighted_C = Malloc(double, nr_class);
		for(i=0;i<nr_class;i++)
			weighted_C[i] = param->C;
		for(i=0;i<param->nr_weight;i++)
		{
			for(j=0;j<nr_class;j++)
				if(param->weight_label[i] == label[j])
					break;
			if(j == nr_class)
				fprintf(stderr,"WARNING: class label %d specified in weight is not found\n", param->weight_label[i]);
			else
				weighted_C[j] *= param->weight[i];
		}

		// constructing the subproblem
		feature_node **x = Malloc(feature_node *,l);
		for(i=0;i<l;i++)
			x[i] = prob->x[perm[i]];

		int k;
		problem sub_prob;
		sub_prob.l = l;
		sub_prob.n = n;
		sub_prob.x = Malloc(feature_node *,sub_prob.l);
		sub_prob.y = Malloc(double,sub_prob.l);

		for(k=0; k<sub_prob.l; k++)
			sub_prob.x[k] = x[k];

		// multi-class svm by Crammer and Singer
		if(param->solver_type == MCSVM_CS)
		{
			model_->w=Malloc(double, n*nr_class);
			for(i=0;i<nr_class;i++)
				for(j=start[i];j<start[i]+count[i];j++)
					sub_prob.y[j] = i;
			Solver_MCSVM_CS Solver(&sub_prob, nr_class, weighted_C, param->eps);
			Solver.Solve(model_->w);
		}
		else
		{
			if(nr_class == 2)
			{
				model_->w=Malloc(double, w_size);

				int e0 = start[0]+count[0];
				k=0;
				for(; k<e0; k++)
					sub_prob.y[k] = +1;
				for(; k<sub_prob.l; k++)
					sub_prob.y[k] = -1;
				
				if(param->init_sol != NULL)
					for(i=0;i<w_size;i++)
						model_->w[i] = param->init_sol[i];
				else
					for(i=0;i<w_size;i++)
						model_->w[i] = 0;

				train_one(&sub_prob, param, model_->w, weighted_C[0], weighted_C[1]);
			}
			else
			{
				model_->w=Malloc(double, w_size*nr_class);
				double *w=Malloc(double, w_size);
				for(i=0;i<nr_class;i++)
				{
					int si = start[i];
					int ei = si+count[i];

					k=0;
					for(; k<si; k++)
						sub_prob.y[k] = -1;
					for(; k<ei; k++)
						sub_prob.y[k] = +1;
					for(; k<sub_prob.l; k++)
						sub_prob.y[k] = -1;

					if(param->init_sol != NULL)
						for(j=0;j<w_size;j++)
							w[j] = param->init_sol[j*nr_class+i];
					else
						for(j=0;j<w_size;j++)
							w[j] = 0;

					train_one(&sub_prob, param, w, weighted_C[i], param->C);

					for(int j=0;j<w_size;j++)
						model_->w[j*nr_class+i] = w[j];
				}
				free(w);
			}

		}

		free(x);
		free(label);
		free(start);
		free(count);
		free(perm);
		free(sub_prob.x);
		free(sub_prob.y);
		free(weighted_C);
	}
	return model_;
}

void cross_validation(const problem *prob, const parameter *param, int nr_fold, double *target)
{
	int i;
	int *fold_start;
	int l = prob->l;
	int *perm = Malloc(int,l);
	if (nr_fold > l)
	{
		nr_fold = l;
		fprintf(stderr,"WARNING: # folds > # data. Will use # folds = # data instead (i.e., leave-one-out cross validation)\n");
	}
	fold_start = Malloc(int,nr_fold+1);
	for(i=0;i<l;i++) perm[i]=i;
	for(i=0;i<l;i++)
	{
		int j = i+rand()%(l-i);
		swap(perm[i],perm[j]);
	}
	for(i=0;i<=nr_fold;i++)
		fold_start[i]=i*l/nr_fold;

	for(i=0;i<nr_fold;i++)
	{
		int begin = fold_start[i];
		int end = fold_start[i+1];
		int j,k;
		struct problem subprob;

		subprob.bias = prob->bias;
		subprob.n = prob->n;
		subprob.l = l-(end-begin);
		subprob.x = Malloc(struct feature_node*,subprob.l);
		subprob.y = Malloc(double,subprob.l);

		k=0;
		for(j=0;j<begin;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			++k;
		}
		for(j=end;j<l;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			++k;
		}
		struct model *submodel = train(&subprob,param);
		for(j=begin;j<end;j++)
			target[perm[j]] = predict(submodel,prob->x[perm[j]]);
		free_and_destroy_model(&submodel);
		free(subprob.x);
		free(subprob.y);
	}
	free(fold_start);
	free(perm);
}

void find_parameter_C(const problem *prob, const parameter *param, int nr_fold, double start_C, double max_C, double *best_C, double *best_rate)
{
	// variables for CV
	int i;
	int *fold_start;
	int l = prob->l;
	int *perm = Malloc(int, l);
	double *target = Malloc(double, prob->l);
	struct problem *subprob = Malloc(problem,nr_fold);

	// variables for warm start
	double ratio = 2;
	double **prev_w = Malloc(double*, nr_fold);
	for(i = 0; i < nr_fold; i++)
		prev_w[i] = NULL;
	int num_unchanged_w = 0;
	struct parameter param1 = *param;
	void (*default_print_string) (const char *) = liblinear_print_string;

	if (nr_fold > l)
	{
		nr_fold = l;
		fprintf(stderr,"WARNING: # folds > # data. Will use # folds = # data instead (i.e., leave-one-out cross validation)\n");
	}
	fold_start = Malloc(int,nr_fold+1);
	for(i=0;i<l;i++) perm[i]=i;
	for(i=0;i<l;i++)
	{
		int j = i+rand()%(l-i);
		swap(perm[i],perm[j]);
	}
	for(i=0;i<=nr_fold;i++)
		fold_start[i]=i*l/nr_fold;

	for(i=0;i<nr_fold;i++)
	{
		int begin = fold_start[i];
		int end = fold_start[i+1];
		int j,k;

		subprob[i].bias = prob->bias;
		subprob[i].n = prob->n;
		subprob[i].l = l-(end-begin);
		subprob[i].x = Malloc(struct feature_node*,subprob[i].l);
		subprob[i].y = Malloc(double,subprob[i].l);

		k=0;
		for(j=0;j<begin;j++)
		{
			subprob[i].x[k] = prob->x[perm[j]];
			subprob[i].y[k] = prob->y[perm[j]];
			++k;
		}
		for(j=end;j<l;j++)
		{
			subprob[i].x[k] = prob->x[perm[j]];
			subprob[i].y[k] = prob->y[perm[j]];
			++k;
		}

	}

	*best_rate = 0;
	if(start_C <= 0)
		start_C = calc_start_C(prob,param);
	param1.C = start_C;

	while(param1.C <= max_C)
	{
		//Output disabled for running CV at a particular C
		set_print_string_function(&print_null);

		for(i=0; i<nr_fold; i++)
		{
			int j;
			int begin = fold_start[i];
			int end = fold_start[i+1];

			param1.init_sol = prev_w[i];
			struct model *submodel = train(&subprob[i],&param1);

			int total_w_size;
			if(submodel->nr_class == 2)
				total_w_size = subprob[i].n;
			else
				total_w_size = subprob[i].n * submodel->nr_class;

			if(prev_w[i] == NULL)
			{
				prev_w[i] = Malloc(double, total_w_size);
				for(j=0; j<total_w_size; j++)
					prev_w[i][j] = submodel->w[j];
			}
			else if(num_unchanged_w >= 0)
			{
				double norm_w_diff = 0;
				for(j=0; j<total_w_size; j++)
				{
					norm_w_diff += (submodel->w[j] - prev_w[i][j])*(submodel->w[j] - prev_w[i][j]);
					prev_w[i][j] = submodel->w[j];
				}
				norm_w_diff = sqrt(norm_w_diff);

				if(norm_w_diff > 1e-15)
					num_unchanged_w = -1;
			}
			else
			{
				for(j=0; j<total_w_size; j++)
					prev_w[i][j] = submodel->w[j];
			}

			for(j=begin; j<end; j++)
				target[perm[j]] = predict(submodel,prob->x[perm[j]]);

			free_and_destroy_model(&submodel);
		}
		set_print_string_function(default_print_string);

		int total_correct = 0;
		for(i=0; i<prob->l; i++)
			if(target[i] == prob->y[i])
				++total_correct;
		double current_rate = (double)total_correct/prob->l;
		if(current_rate > *best_rate)
		{
			*best_C = param1.C;
			*best_rate = current_rate;
		}

		info("log2c=%7.2f\trate=%g\n",log(param1.C)/log(2.0),100.0*current_rate);
		num_unchanged_w++;
		if(num_unchanged_w == 3)
			break;
		param1.C = param1.C*ratio;
	}

	if(param1.C > max_C && max_C > start_C) 
		info("warning: maximum C reached.\n");
	free(fold_start);
	free(perm);
	free(target);
	for(i=0; i<nr_fold; i++)
	{
		free(subprob[i].x);
		free(subprob[i].y);
		free(prev_w[i]);
	}
	free(prev_w);
	free(subprob);
}

double predict_values(const struct model *model_, const struct feature_node *x, double *dec_values)
{
	int idx;
	int n;
	if(model_->bias>=0)
		n=model_->nr_feature+1;
	else
		n=model_->nr_feature;
	double *w=model_->w;
	int nr_class=model_->nr_class;
	int i;
	int nr_w;
	if(nr_class==2 && model_->param.solver_type != MCSVM_CS)
		nr_w = 1;
	else
		nr_w = nr_class;

	const feature_node *lx=x;
	for(i=0;i<nr_w;i++)
		dec_values[i] = 0;
	for(; (idx=lx->index)!=-1; lx++)
	{
		// the dimension of testing data may exceed that of training
		if(idx<=n)
			for(i=0;i<nr_w;i++)
				dec_values[i] += w[(idx-1)*nr_w+i]*lx->value;
	}

	if(nr_class==2)
	{
		if(check_regression_model(model_))
			return dec_values[0];
		else if(model_->param.solver_type != ONECLASS_L1_RD_1000
		 	||model_->param.solver_type != ONECLASS_L1_SEMIGD_1000
		 	||model_->param.solver_type != ONECLASS_L1_RD_SH
		 	||model_->param.solver_type != ONECLASS_L1_SEMIGD_SH
			||model_->param.solver_type != ONECLASS_L1_FIRST_1000
			||model_->param.solver_type != ONECLASS_L1_SECOND_1000
			||model_->param.solver_type != ONECLASS_L2_SECOND_1000
			||model_->param.solver_type != ONECLASS_L2_FIRST_1000
	 		||model_->param.solver_type != ONECLASS_L2_RD_1000
		 	||model_->param.solver_type != ONECLASS_L2_SEMIGD_1000
		 	||model_->param.solver_type != ONECLASS_L2_RD_SH
	 		||model_->param.solver_type != ONECLASS_L2_SEMIGD_SH)
			return (dec_values[0]>0)?1:-1;
		else
			return (dec_values[0]>0)?model_->label[0]:model_->label[1];
	}
	else
	{
		int dec_max_idx = 0;
		for(i=1;i<nr_class;i++)
		{
			if(dec_values[i] > dec_values[dec_max_idx])
				dec_max_idx = i;
		}
		return model_->label[dec_max_idx];
	}
}

double predict(const model *model_, const feature_node *x)
{
	double *dec_values = Malloc(double, model_->nr_class);
	double label=predict_values(model_, x, dec_values);
	free(dec_values);
	return label;
}

double predict_probability(const struct model *model_, const struct feature_node *x, double* prob_estimates)
{
	if(check_probability_model(model_))
	{
		int i;
		int nr_class=model_->nr_class;
		int nr_w;
		if(nr_class==2)
			nr_w = 1;
		else
			nr_w = nr_class;

		double label=predict_values(model_, x, prob_estimates);
		for(i=0;i<nr_w;i++)
			prob_estimates[i]=1/(1+exp(-prob_estimates[i]));

		if(nr_class==2) // for binary classification
			prob_estimates[1]=1.-prob_estimates[0];
		else
		{
			double sum=0;
			for(i=0; i<nr_class; i++)
				sum+=prob_estimates[i];

			for(i=0; i<nr_class; i++)
				prob_estimates[i]=prob_estimates[i]/sum;
		}

		return label;
	}
	else
		return 0;
}

class Solver_type_table
{
	public:
	Solver_type_table()
	{
		for(size_t i=0; i<init_size; i++)
		{
			array[i]="UNDEFINED";
		}
#define SAVE_NAME(p) array[p] = #p;
		SAVE_NAME(L2R_LR);
		SAVE_NAME(ONE_L2_CY_SH);
		SAVE_NAME(L2R_L2LOSS_SVC);
		SAVE_NAME(ONE_L1_CY_SH);
		SAVE_NAME(MCSVM_CS);
		SAVE_NAME(L1R_L2LOSS_SVC);
		SAVE_NAME(L1R_LR);
		SAVE_NAME(L2R_LR_DUAL);
		SAVE_NAME(L2R_L2LOSS_SVR);
		SAVE_NAME(L2R_L2LOSS_SVR_DUAL);
		SAVE_NAME(L2R_L1LOSS_SVR_DUAL);
		//for one-variable
		SAVE_NAME(ONE_L1_CY_1000);
		SAVE_NAME(ONE_L1_RD_1000);
		SAVE_NAME(ONE_L1_RD_SH);
		SAVE_NAME(ONE_L2_CY_1000);
		SAVE_NAME(ONE_L2_RD_1000);
		SAVE_NAME(ONE_L2_RD_SH);
		SAVE_NAME(ONE_L1_SEMIGD_1000);
		SAVE_NAME(ONE_L1_SEMIGD_SH);
		SAVE_NAME(ONE_L2_SEMIGD_1000);
		SAVE_NAME(ONE_L2_SEMIGD_SH);
		//for two-variable
		SAVE_NAME(TWO_L1_CY_1000);
		SAVE_NAME(TWO_L1_RD_1000);
		SAVE_NAME(TWO_L1_RD_SH);
		SAVE_NAME(TWO_L2_CY_1000);
		SAVE_NAME(TWO_L2_RD_1000);
		SAVE_NAME(TWO_L2_RD_SH);
		SAVE_NAME(TWO_L1_SEMICY_1000);
		SAVE_NAME(TWO_L2_SEMICY_1000);
		SAVE_NAME(TWO_L1_SEMIRDONE_1000);
		SAVE_NAME(TWO_L2_SEMIRDONE_1000);
		SAVE_NAME(TWO_L1_SEMIRDTWO_1000);
		SAVE_NAME(TWO_L2_SEMIRDTWO_1000);
		SAVE_NAME(TWO_L1_RD_SH2);
		SAVE_NAME(TWO_L2_RD_SH2);
		SAVE_NAME(TWO_L1_SEMIGD_1000);
		SAVE_NAME(TWO_L1_SEMIGD_SH);
		SAVE_NAME(TWO_L2_SEMIGD_1000);
		SAVE_NAME(TWO_L2_SEMIGD_SH);
		//for two-variable linear constraints
		SAVE_NAME(BIAS_L1_RD_1000);
		SAVE_NAME(BIAS_L2_RD_1000);
		SAVE_NAME(BIAS_L1_SEMIGD_1000);
		SAVE_NAME(BIAS_L2_SEMIGD_1000);
		SAVE_NAME(BIAS_L1_RD_SH);
		SAVE_NAME(BIAS_L2_RD_SH);
		SAVE_NAME(BIAS_L1_SEMIGD_SH);
		SAVE_NAME(BIAS_L2_SEMIGD_SH);
		//for one-class svm
		SAVE_NAME(ONECLASS_L1_RD_1000);
		SAVE_NAME(ONECLASS_L1_SEMIGD_1000);
		SAVE_NAME(ONECLASS_L1_FIRST_1000);
		SAVE_NAME(ONECLASS_L1_SECOND_1000);
		SAVE_NAME(ONECLASS_L1_RD_SH);
		SAVE_NAME(ONECLASS_L1_SEMIGD_SH);
		SAVE_NAME(ONECLASS_L2_RD_1000);
		SAVE_NAME(ONECLASS_L2_SEMIGD_1000);
		SAVE_NAME(ONECLASS_L2_FIRST_1000);
		SAVE_NAME(ONECLASS_L2_SECOND_1000);
		SAVE_NAME(ONECLASS_L2_RD_SH);
		SAVE_NAME(ONECLASS_L2_SEMIGD_SH);
#undef SAVE_NAME
	}

	const char *operator[](size_t index) const
	{
		if(index > init_size)
		{
			fprintf(stderr, "ERROR: index exceeded array size\n");
			return NULL;
		}
		return array[index];
	}
	static const char * UNDEFINED;

	private:
		static const size_t init_size=1024;
		const char * array[init_size];
}solver_type_table;

int save_model(const char *model_file_name, const struct model *model_)
{
	int i;
	int nr_feature=model_->nr_feature;
	int n;
	const parameter& param = model_->param;

	if(model_->bias>=0)
		n=nr_feature+1;
	else
		n=nr_feature;
	int w_size = n;
	FILE *fp = fopen(model_file_name,"w");
	if(fp==NULL) return -1;

	char *old_locale = setlocale(LC_ALL, NULL);
	if (old_locale)
	{
		old_locale = strdup(old_locale);
	}
	setlocale(LC_ALL, "C");

	int nr_w;
	if(model_->nr_class==2 && model_->param.solver_type != MCSVM_CS)
		nr_w=1;
	else
		nr_w=model_->nr_class;

	fprintf(fp, "solver_type %s\n", solver_type_table[param.solver_type]);
	fprintf(fp, "nr_class %d\n", model_->nr_class);

	if(model_->label)
	{
		fprintf(fp, "label");
		for(i=0; i<model_->nr_class; i++)
			fprintf(fp, " %d", model_->label[i]);
		fprintf(fp, "\n");
	}

	fprintf(fp, "nr_feature %d\n", nr_feature);

	fprintf(fp, "bias %.16g\n", model_->bias);

	fprintf(fp, "w\n");
	for(i=0; i<w_size; i++)
	{
		int j;
		for(j=0; j<nr_w; j++)
			fprintf(fp, "%.16g ", model_->w[i*nr_w+j]);
		fprintf(fp, "\n");
	}

	setlocale(LC_ALL, old_locale);
	free(old_locale);

	if (ferror(fp) != 0 || fclose(fp) != 0) return -1;
	else return 0;
}

//
// FSCANF helps to handle fscanf failures.
// Its do-while block avoids the ambiguity when
// if (...)
//    FSCANF();
// is used
//
#define FSCANF(_stream, _format, _var)do\
{\
	if (fscanf(_stream, _format, _var) != 1)\
	{\
		fprintf(stderr, "ERROR: fscanf failed to read the model\n");\
		EXIT_LOAD_MODEL()\
	}\
}while(0)
// EXIT_LOAD_MODEL should NOT end with a semicolon.
#define EXIT_LOAD_MODEL()\
{\
	setlocale(LC_ALL, old_locale);\
	free(model_->label);\
	free(model_);\
	free(old_locale);\
	return NULL;\
}
struct model *load_model(const char *model_file_name)
{
	FILE *fp = fopen(model_file_name,"r");
	if(fp==NULL) return NULL;

	int i;
	int nr_feature;
	int n;
	int nr_class;
	double bias;
	model *model_ = Malloc(model,1);
	parameter& param = model_->param;

	model_->label = NULL;

	char *old_locale = setlocale(LC_ALL, NULL);
	if (old_locale)
	{
		old_locale = strdup(old_locale);
	}
	setlocale(LC_ALL, "C");

	char cmd[81];
	while(1)
	{
		FSCANF(fp,"%80s",cmd);
		if(strcmp(cmd,"solver_type")==0)
		{
			FSCANF(fp,"%80s",cmd);
			int i;
			for(i=0;solver_type_table[i];i++)
			{
				if(strcmp(solver_type_table[i],cmd)==0)
				{
					param.solver_type=i;
					break;
				}
			}
			if(solver_type_table[i] == NULL)
			{
				fprintf(stderr,"unknown solver type.\n");
				EXIT_LOAD_MODEL()
			}
		}
		else if(strcmp(cmd,"nr_class")==0)
		{
			FSCANF(fp,"%d",&nr_class);
			model_->nr_class=nr_class;
		}
		else if(strcmp(cmd,"nr_feature")==0)
		{
			FSCANF(fp,"%d",&nr_feature);
			model_->nr_feature=nr_feature;
		}
		else if(strcmp(cmd,"bias")==0)
		{
			FSCANF(fp,"%lf",&bias);
			model_->bias=bias;
		}
		else if(strcmp(cmd,"w")==0)
		{
			break;
		}
		else if(strcmp(cmd,"label")==0)
		{
			int nr_class = model_->nr_class;
			model_->label = Malloc(int,nr_class);
			for(int i=0;i<nr_class;i++)
				FSCANF(fp,"%d",&model_->label[i]);
		}
		else
		{
			fprintf(stderr,"unknown text in model file: [%s]\n",cmd);
			EXIT_LOAD_MODEL()
		}
	}

	nr_feature=model_->nr_feature;
	if(model_->bias>=0)
		n=nr_feature+1;
	else
		n=nr_feature;
	int w_size = n;
	int nr_w;
	if(nr_class==2 && param.solver_type != MCSVM_CS)
		nr_w = 1;
	else
		nr_w = nr_class;

	model_->w=Malloc(double, w_size*nr_w);
	for(i=0; i<w_size; i++)
	{
		int j;
		for(j=0; j<nr_w; j++)
			FSCANF(fp, "%lf ", &model_->w[i*nr_w+j]);
		if (fscanf(fp, "\n") !=0)
		{
			fprintf(stderr, "ERROR: fscanf failed to read the model\n");
			EXIT_LOAD_MODEL()
		}
	}

	setlocale(LC_ALL, old_locale);
	free(old_locale);

	if (ferror(fp) != 0 || fclose(fp) != 0) return NULL;

	return model_;
}

int get_nr_feature(const model *model_)
{
	return model_->nr_feature;
}

int get_nr_class(const model *model_)
{
	return model_->nr_class;
}

void get_labels(const model *model_, int* label)
{
	if (model_->label != NULL)
		for(int i=0;i<model_->nr_class;i++)
			label[i] = model_->label[i];
}

// use inline here for better performance (around 20% faster than the non-inline one)
static inline double get_w_value(const struct model *model_, int idx, int label_idx) 
{
	int nr_class = model_->nr_class;
	int solver_type = model_->param.solver_type;
	const double *w = model_->w;

	if(idx < 0 || idx > model_->nr_feature)
		return 0;
	if(check_regression_model(model_))
		return w[idx];
	else 
	{
		if(label_idx < 0 || label_idx >= nr_class)
			return 0;
		if(nr_class == 2 && solver_type != MCSVM_CS)
		{
			if(label_idx == 0)
				return w[idx];
			else
				return -w[idx];
		}
		else
			return w[idx*nr_class+label_idx];
	}
}

// feat_idx: starting from 1 to nr_feature
// label_idx: starting from 0 to nr_class-1 for classification models;
//            for regression models, label_idx is ignored.
double get_decfun_coef(const struct model *model_, int feat_idx, int label_idx)
{
	if(feat_idx > model_->nr_feature)
		return 0;
	return get_w_value(model_, feat_idx-1, label_idx);
}

double get_decfun_bias(const struct model *model_, int label_idx)
{
	int bias_idx = model_->nr_feature;
	double bias = model_->bias;
	if(bias <= 0)
		return 0;
	else
		return bias*get_w_value(model_, bias_idx, label_idx);
}

void free_model_content(struct model *model_ptr)
{
	if(model_ptr->w != NULL)
		free(model_ptr->w);
	if(model_ptr->label != NULL)
		free(model_ptr->label);
}

void free_and_destroy_model(struct model **model_ptr_ptr)
{
	struct model *model_ptr = *model_ptr_ptr;
	if(model_ptr != NULL)
	{
		free_model_content(model_ptr);
		free(model_ptr);
	}
}

void destroy_param(parameter* param)
{
	if(param->weight_label != NULL)
		free(param->weight_label);
	if(param->weight != NULL)
		free(param->weight);
	if(param->init_sol != NULL)
		free(param->init_sol);
}

const char *check_parameter(const problem *prob, const parameter *param)
{
	if(param->eps <= 0)
		return "eps <= 0";

	if(param->C <= 0)
		return "C <= 0";

	if(param->p < 0)
		return "p < 0";

	if(param->solver_type != L2R_LR
		&& param->solver_type != ONE_L2_CY_SH
		&& param->solver_type != L2R_L2LOSS_SVC
		&& param->solver_type != ONE_L1_CY_SH
		&& param->solver_type != MCSVM_CS
		&& param->solver_type != L1R_L2LOSS_SVC
		&& param->solver_type != L1R_LR
		&& param->solver_type != L2R_LR_DUAL
		&& param->solver_type != L2R_L2LOSS_SVR
		&& param->solver_type != L2R_L2LOSS_SVR_DUAL
		&& param->solver_type != L2R_L1LOSS_SVR_DUAL
		&& param->solver_type != ONE_L1_CY_1000
		&& param->solver_type != ONE_L1_RD_1000
		&& param->solver_type != ONE_L1_RD_SH
		&& param->solver_type != ONE_L2_CY_1000
		&& param->solver_type != ONE_L2_RD_1000
		&& param->solver_type != ONE_L2_RD_SH
		&& param->solver_type != ONE_L1_SEMIGD_1000
		&& param->solver_type != ONE_L1_SEMIGD_SH
		&& param->solver_type != ONE_L2_SEMIGD_1000
		&& param->solver_type != ONE_L2_SEMIGD_SH
		&& param->solver_type != TWO_L1_CY_1000
		&& param->solver_type != TWO_L1_RD_1000
		&& param->solver_type != TWO_L1_RD_SH
		&& param->solver_type != TWO_L2_CY_1000
		&& param->solver_type != TWO_L2_RD_1000
		&& param->solver_type != TWO_L2_RD_SH
		&& param->solver_type != TWO_L1_SEMICY_1000
		&& param->solver_type != TWO_L2_SEMICY_1000
		&& param->solver_type != TWO_L1_SEMIRDONE_1000
		&& param->solver_type != TWO_L2_SEMIRDONE_1000
		&& param->solver_type != TWO_L1_SEMIRDTWO_1000
		&& param->solver_type != TWO_L2_SEMIRDTWO_1000
		&& param->solver_type != TWO_L1_RD_SH2
		&& param->solver_type != TWO_L2_RD_SH2
		&& param->solver_type != TWO_L1_SEMIGD_1000
		&& param->solver_type != TWO_L1_SEMIGD_SH
		&& param->solver_type != TWO_L2_SEMIGD_1000
		&& param->solver_type != TWO_L2_SEMIGD_SH
		&& param->solver_type != BIAS_L1_RD_1000
		&& param->solver_type != BIAS_L2_RD_1000
		&& param->solver_type != BIAS_L1_SEMIGD_1000
		&& param->solver_type != BIAS_L2_SEMIGD_1000
		&& param->solver_type != BIAS_L1_SEMIGD_SH
		&& param->solver_type != BIAS_L2_SEMIGD_SH
		&& param->solver_type != BIAS_L1_RD_SH
		&& param->solver_type != BIAS_L2_RD_SH
	 	&& param->solver_type != ONECLASS_L1_RD_1000
	 	&& param->solver_type != ONECLASS_L1_SEMIGD_1000
	 	&& param->solver_type != ONECLASS_L1_RD_SH
	 	&& param->solver_type != ONECLASS_L1_SEMIGD_SH
		&& param->solver_type != ONECLASS_L1_FIRST_1000
		&& param->solver_type != ONECLASS_L1_SECOND_1000
		&& param->solver_type != ONECLASS_L2_SECOND_1000
		&& param->solver_type != ONECLASS_L2_FIRST_1000
	 	&& param->solver_type != ONECLASS_L2_RD_1000
	 	&& param->solver_type != ONECLASS_L2_SEMIGD_1000
	 	&& param->solver_type != ONECLASS_L2_RD_SH
	 	&& param->solver_type != ONECLASS_L2_SEMIGD_SH)
		return "unknown solver type";

	if(param->init_sol != NULL 
		&& param->solver_type != L2R_LR && param->solver_type != L2R_L2LOSS_SVC)
		return "Initial-solution specification supported only for solver L2R_LR and L2R_L2LOSS_SVC";

	return NULL;
}

int check_probability_model(const struct model *model_)
{
	return (model_->param.solver_type==L2R_LR ||
			model_->param.solver_type==L2R_LR_DUAL ||
			model_->param.solver_type==L1R_LR);
}

int check_regression_model(const struct model *model_)
{
	return (model_->param.solver_type==L2R_L2LOSS_SVR ||
			model_->param.solver_type==L2R_L1LOSS_SVR_DUAL ||
			model_->param.solver_type==L2R_L2LOSS_SVR_DUAL);
}

void set_print_string_function(void (*print_func)(const char*))
{
	if (print_func == NULL)
		liblinear_print_string = &print_string_stdout;
	else
		liblinear_print_string = print_func;
}

