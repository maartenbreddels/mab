#pragma once
#include "math.h"
#include<iostream>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <functional>
#include <stdexcept>

using namespace std;


namespace gd {


template<typename C>
double gslWrapper(double x, void * params);

static void check_gsl_rootfinder_result(int result) {
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
double gslWrapper2(double x, void * params);

template<typename C>
double gslWrapper3a(double x, void * params);
template<typename C>
double gslWrapper3b(double x, void * params);
template<typename C>
void gslWrapper3c(double x, void * params, double * y, double *dy);

template<typename F=std::function<double(double)>, typename T=double>
	class RootFinderGSL {
	public:
		F f;
	RootFinderGSL(F f) : f(f) {
	}
		
	T findRoot(T a, T b, double epsabs=0., double epsrel=1e-6) {
			int status;
			int iter = 0, max_iter = 100;
			const gsl_root_fsolver_type *rT;
			gsl_root_fsolver *s;
			double r = 0; //, r_expected = 0;
			//double x_lo = 0.0, x_hi = 5.0;
			
			
			gsl_function gsl_f;
			T (* function_pointer)(T, void * params) = &(gslWrapper2<RootFinderGSL<F,T> >);
			gsl_f.function = function_pointer;
			gsl_f.params = this;
			
			
			rT = gsl_root_fsolver_brent;
			s = gsl_root_fsolver_alloc (rT);
			gsl_root_fsolver_set (s, &gsl_f, a, b);
			
			/*printf ("using %s method\n", 
					gsl_root_fsolver_name (s));*/
			
			/*printf ("%5s [%9s, %9s] %9s %10s %9s\n",
					"iter", "lower", "upper", "root", 
					"err", "err(est)");*/
			
			do
			{
				iter++;
				status = gsl_root_fsolver_iterate (s);
				r = gsl_root_fsolver_root (s);
				a = gsl_root_fsolver_x_lower (s);
				b = gsl_root_fsolver_x_upper (s);
				status = gsl_root_test_interval (a, b,
				                                 epsabs, epsrel);
				
				/*if (status == GSL_SUCCESS)
					printf ("Converged:\n");
				
				printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
						iter, a, b,
						r, r - r_expected, 
						a - b);*/
				check_gsl_rootfinder_result(status);
			} while (status == GSL_CONTINUE && iter < max_iter);
			gsl_root_fsolver_free (s);
			//return status;
			return r;
			
		}
		T dI(T x) {
			return f(x);
		}
		
	};

template<typename F=std::function<void(double, double*, double*)>, typename T=double>
	class RootFinderGSLDerivative {
	public:
		F f;
		const gsl_root_fdfsolver_type *rT;
	RootFinderGSLDerivative(F f, const gsl_root_fdfsolver_type* rT=gsl_root_fdfsolver_secant) : f(f), rT(rT) {
	}
		
		T findRoot(T x_begin, double epsabs=0., double epsrel=1e-6) {
			int status;
			int iter = 0, max_iter = 100;
			gsl_root_fdfsolver *s;
			double x0 = x_begin;
			double x = x0;
			//double x_lo = 0.0, x_hi = 5.0;
			
			
			gsl_function_fdf gsl_f;
			double (* function_pointer_a)(T, void * params) = &(gslWrapper3a<RootFinderGSLDerivative<F,T> >);
			double (* function_pointer_b)(T, void * params) = &(gslWrapper3b<RootFinderGSLDerivative<F,T> >);
			void (* function_pointer_c)(T, void * params, T*, T*) = &(gslWrapper3c<RootFinderGSLDerivative<F,T> >);
			gsl_f.f = function_pointer_a;
			gsl_f.df = function_pointer_b;
			gsl_f.fdf = function_pointer_c;
			gsl_f.params = this;
			
			
			s = gsl_root_fdfsolver_alloc (rT);
			gsl_root_fdfsolver_set (s, &gsl_f, x);
			
			/*printf ("using %s method\n", 
					gsl_root_fsolver_name (s));*/
			
			/*printf ("%5s [%9s, %9s] %9s %10s %9s\n",
					"iter", "lower", "upper", "root", 
					"err", "err(est)");*/
			
			do
			{
				iter++;
				status = gsl_root_fdfsolver_iterate (s);
				check_gsl_rootfinder_result(status);
				x0 = x;
				x = gsl_root_fdfsolver_root (s);
				status = gsl_root_test_delta (x, x0, epsabs, epsrel);
				check_gsl_rootfinder_result(status);
				
				//a = gsl_root_fdfsolver_x_lower (s);
				//b = gsl_root_fdfsolver_x_upper (s);
				//status = gsl_root_test_interval (a, b,
				//                               0, 0.001);
				
				/*
				if (status == GSL_SUCCESS)
					printf ("Converged:\n");
				
				printf ("%5d %10.7f %10.7f\n",
				        iter, x, x - x0);
				
				*/
			} while (status == GSL_CONTINUE && iter < max_iter);
			gsl_root_fdfsolver_free (s);
			//return status;
			return x;
			
		}
		
	};

template<typename C>
double gslWrapper2(double x, void * params) {
	C* object = (C*) params;
	return object->dI(x);
}

template<class C>
double gslWrapper3a(double x, void * params) {
	C* object = (C*) params;
	double y, dy;
	object->f(x, &y, &dy);
	return y;
}

template<class C>
double gslWrapper3b(double x, void * params) {
	C* object = (C*) params;
	double y, dy;
	object->f(x, &y, &dy);
	return dy;
}

template<class C>
void gslWrapper3c(double x, void * params, double * y, double *dy) {
	C* object = (C*) params;
	object->f(x, y, dy);
}


}
