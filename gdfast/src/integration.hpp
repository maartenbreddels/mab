#pragma once
#include "math.h"
#include<iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <functional>
#include <stdexcept>

using namespace std;


namespace gd {

template<typename F=std::function<double(double)>, typename T=double>
class IntegratorGrid {
public:
	F f;
	T x1, x2, dx;
	T *ygrid;
	int N;
	IntegratorGrid(F f, T x1, T x2, T dx) : f(f), x1(x1), x2(x2), dx(dx) {
		N = int((x2-x1)/dx);
		ygrid = new T[N];
		int i = 0;
		double x = 0;
		while(i < N) {
			ygrid[i] = f(x); 
			x += dx;
			i++;
		}
	}
	
	T integrate(T a, T b) {
		int istart = int((a-x1)/dx);
		int iend = int((b-x1)/dx);
		double xoffset = a - istart*dx;
		//cout << "i - " << istart << " - " << iend << endl;
		//cout << istart << " " << xoffset << endl;
		double dy = (ygrid[istart+1]-ygrid[istart]);
		double y1 = ygrid[istart+1];
		double ystart =  ygrid[istart] + dy/dx * xoffset;
		//cout << ystart << endl;
		double I = 0;
		double h = dx-xoffset;
		I += h/2 * (ystart+y1);
		int i = istart;
		while(i < iend) {
			double dI = dx/2 * (ygrid[i] +ygrid[i+1]);
			I += dI;  
			i++;
		}
		
		xoffset = b - iend*dx;
		dy = (ygrid[iend+1]-ygrid[iend]);
		y1 = ygrid[iend];
		double yend =  ygrid[iend] + dy/dx * xoffset;
		//cout << "yend " << yend << endl;
		h = xoffset;
		I += h/2 * (y1+yend);
		
		return I;
	}
};

template<typename C>
double gslWrapper(double x, void * params);

static void check_gsl_integration_result(int result) {
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

template<typename F=std::function<double(double)>, typename T=double>
class IntegratorGSL {
	public:
		F f;
		double epsrel, epsabs;
		int limit;
		gsl_integration_workspace * w; 
		gsl_function gsl_f;
		
		IntegratorGSL(F f, int limit=10000, double epsrel=1e-6, double epsabs=1e-6) : f(f), epsrel(epsrel), epsabs(epsabs), limit(limit) {
			T (* function_pointer)(T, void * params) = &(gslWrapper<IntegratorGSL<F,T> >);
			gsl_f.function = function_pointer;
			gsl_f.params = this;
			w = gsl_integration_workspace_alloc (limit);
		}
		~IntegratorGSL() {
			gsl_integration_workspace_free(w);
		}
		T integrate(T a, T b) {
			double integral = 0;
			double err = 0;
			int result = gsl_integration_qags(&gsl_f, a, b, epsabs, epsrel, limit, w, &integral, &err);
			check_gsl_integration_result(result);
			return integral;
		}
		T integrate_non_adaptive(T a, T b) {
			double integral = 0;
			double err = 0;
			size_t nevals;
			int result = gsl_integration_qng(&gsl_f, a, b, epsabs, epsrel, &integral, &err, &nevals);
			check_gsl_integration_result(result);
			return integral;
		}
		T integrate_no_singularities(T a, T b, int key=GSL_INTEG_GAUSS15) {
			double integral = 0;
			double err = 0;
			size_t nevals;
			int result = gsl_integration_qag(&gsl_f, a, b, epsabs, epsrel, limit, key, w, &integral, &err);
			check_gsl_integration_result(result);
			return integral;
		}
		T integrate_to_sigularity(T a, T b) {
			double integral = 0;
			double err = 0;
			int result;
			double points[2];
			points[0] = a;
			points[1] = b;
			result = gsl_integration_qagp(&gsl_f, points, 2, epsabs, epsrel, limit, w, &integral, &err);
			check_gsl_integration_result(result);
			return integral;
		}
		T integrate_to_inf(T a, int limit=1000) {
			double integral = 0;
			double err = 0;
			int result;
			result = gsl_integration_qagiu(&gsl_f, a, 0, 1e-4, limit, w, &integral, &err);
			check_gsl_integration_result(result);
			return integral;
		}
		T integrate_inf_to_inf(int limit=1000) {
			double integral = 0;
			double err = 0;
			int result;
			result = gsl_integration_qagi(&gsl_f, 0, 1e-2, limit, w, &integral, &err);
			check_gsl_integration_result(result);
			return integral;
		}
		T dI(T x) {
			return f(x);
		}
};

template<typename C>
double gslWrapper(double x, void * params) {
	C* object = (C*) params;
	return object->dI(x);
}


}

