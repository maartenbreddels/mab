#pragma once
#include <math.h>
#include <iostream>
#include <gsl/gsl_min.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <functional>
#include <nlopt.h>

using namespace std;


namespace gd {


static double nlopt_func(unsigned n, const double *x, double *grad, void *data);
static double nlopt_func_grad(unsigned n, const double *x, double *grad, void *data);

template<typename F=std::function<void(double)>, typename T=double>
class MinimizerNLopt {
	public:
		F f;
		MinimizerNLopt(F f) : f(f) {
		}
	
		T findMinimum(T a, T b) {
			nlopt_opt opt;
			opt = nlopt_create(NLOPT_LD_MMA, 2);
			//nlopt_set_lower_bounds(opt, lb);
			//nlopt_set_min_objective(opt, myfunc, NULL);
			double lower_bounds[2] = { 0, 0 };
			nlopt_set_lower_bounds(opt, lower_bounds);
			nlopt_set_xtol_rel(opt, 1e-4);
			double x[2] = { 1.234, 5.678 };  /* some initial guess */
			double minf; /* the minimum objective value, upon return */
			
			if (nlopt_optimize(opt, x, &minf) < 0) {
				printf("nlopt failed!\n");
			}
			else {
				printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
			}
			nlopt_destroy(opt);
		}

};
template<typename F=std::function<double(double)>, typename T=double>
class MinimizerNLopt2 {
	public:
		F f;
		nlopt_opt opt;
		double x;
		MinimizerNLopt2(F f, nlopt_algorithm algorithm=NLOPT_LN_NELDERMEAD) : f(f) {
			opt = nlopt_create(algorithm, 1);
			x = 0;
			nlopt_set_min_objective(opt, nlopt_func, this);
		}
	
		T optimize() {
			//nlopt_set_lower_bounds(opt, lb);
			//double lower_bounds[2] = { 0, 0 };
			//nlopt_set_lower_bounds(opt, lower_bounds);
			//nlopt_set_xtol_rel(opt, 1e-4);
			//double x[2] = { 1.234, 5.678 };  /* some initial guess */
			double minf; /* the minimum objective value, upon return */
			
			if (nlopt_optimize(opt, &x, &minf) < 0) {
				printf("nlopt failed!\n");
			}
			else {
				//printf("found minimum at f(%g) = %0.10g\n", x, minf);
			}
			return x;
			//nlopt_destroy(opt);
		}

};

static double nlopt_func(unsigned n, const double *x, double *grad, void *data) {
	//
	MinimizerNLopt2<> *opt = (MinimizerNLopt2<>*)data;
	/*data.f( 
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }*/
    return opt->f(*x);
}


template<typename F=std::function<double(double)>, typename T=double>
class MinimizerNLopt2grad {
	public:
		F f;
		F g;
		nlopt_opt opt;
		double x;
		MinimizerNLopt2grad(F f, F g, nlopt_algorithm algorithm=NLOPT_LD_LBFGS) : f(f), g(g) {
			opt = nlopt_create(algorithm, 1);
			x = 0;
			nlopt_set_min_objective(opt, nlopt_func_grad, this);
		}
	
		T optimize() {
			//nlopt_set_lower_bounds(opt, lb);
			//double lower_bounds[2] = { 0, 0 };
			//nlopt_set_lower_bounds(opt, lower_bounds);
			//nlopt_set_xtol_rel(opt, 1e-4);
			//double x[2] = { 1.234, 5.678 };  /* some initial guess */
			double minf; /* the minimum objective value, upon return */
			int status = nlopt_optimize(opt, &x, &minf);
			if (status < 0) {
				printf("nlopt failed! status %d\n", status);
			}
			else {
				//printf("found minimum at f(%g) = %0.10g\n", x, minf);
			}
			return x;
			//nlopt_destroy(opt);
		}

};

static double nlopt_func_grad(unsigned n, const double *x, double *grad, void *data) {
	//
	MinimizerNLopt2grad<> *opt = (MinimizerNLopt2grad<>*)data;
	/*data.f( 
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }*/
    double v = opt->f(*x);
	double g = opt->g(*x);
	if(grad)
		grad[0] = g;
	return v;
}




template<typename C>
double gslWrapper(double x, void * params);

static void check_gsl_minimizer_result(int result) {
	if(result == GSL_EMAXITER) {
		throw std::range_error("the maximum number of subdivisions was exceeded\n");
	}
	if(result == GSL_EROUND) {
		throw std::range_error("cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table. \n");
	}
	if(result == GSL_ESING) {
		throw std::range_error("a non-integrable singularity or other bad integrand behavior was found in the integration interval. \n");
	}
	if(result == GSL_EDIVERGE) {
		throw std::range_error("the integral is divergent, or too slowly convergent to be integrated numerically. \n");
	}
}

template<typename C>
double gslWrapper3(double x, void * params);

template<typename F=std::function<double(double)>, typename T=double>
class MinimizerGSL {
	public:
		F f;
		MinimizerGSL(F f) : f(f) {
		}
	
		T findMinimum(T a, T b) {
			int status;
			int iter = 0, max_iter = 100;
			const gsl_min_fminimizer_type *rT;
			gsl_min_fminimizer *s;
			double r = a+1e-3; //(a+b)/2;
			//, r_expected = 0;
			//double x_lo = 0.0, x_hi = 5.0;


			gsl_function gsl_f;
			T (* function_pointer)(T, void * params) = &(gslWrapper3<MinimizerGSL<F,T> >);
			gsl_f.function = function_pointer;
			gsl_f.params = this;


			rT = gsl_min_fminimizer_brent;
			s = gsl_min_fminimizer_alloc (rT);
			gsl_min_fminimizer_set (s, &gsl_f, r, a, b);
			
			/*printf ("using %s method\n", 
					gsl_root_fsolver_name (s));*/
			
			/*printf ("%5s [%9s, %9s] %9s %10s %9s\n",
					"iter", "lower", "upper", "root", 
					"err", "err(est)");*/
			printf("starting..\n");
			do
				{
				iter++;
				status = gsl_min_fminimizer_iterate (s);
				r = gsl_min_fminimizer_x_minimum (s);
				a = gsl_min_fminimizer_x_lower (s);
				b = gsl_min_fminimizer_x_upper (s);
				status = gsl_min_test_interval (a, b, 0.001, 0);
			
				/*if (status == GSL_SUCCESS)
					printf ("Converged:\n");
				*/
				printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
						iter, a, b,
						r,  a - b);
				check_gsl_minimizer_result(status);
			} while (status == GSL_CONTINUE && iter < max_iter);
			gsl_min_fminimizer_free (s);
			//return status;
			return r;

		}

};

template<typename C>
double gslWrapper3(double x, void * params) {
	C* object = (C*) params;
	return object->f(x);
}


}
