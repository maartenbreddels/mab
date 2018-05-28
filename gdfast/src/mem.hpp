#pragma once
#include <nlopt.h>
#include <cmath>
#include <assert.h>
#include <stdio.h>
#include "integration.hpp"

double nlopt_func_mem(unsigned n, const double *x, double *grad, void *my_func_data);
double myconstraint_mem(unsigned n, const double *x, double *grad, void *data);

namespace gd {
	class MEM {
	public:
		double sigma, gamma, xmax;
		double norm, l2, l4;
		MEM(double sigma, double gamma, double xmax) : sigma(sigma), gamma(gamma), xmax(xmax) {
			optimize();
		}
		double fopt(double norm, double l2, double l4) {
			this->norm = fabs(norm);
			this->l2 = l2;
			this->l4 = l4;
			double m0 = this->moment(0);
			double m2 = this->moment(2);
			double m4 = this->moment(4);
			double k = m4/(m2*m2);
			double b = this->operator()(xmax);
			printf("%f %f %f\n", m0, m2, m4);
			printf(" %f %f %f ", m0, m2/(sigma*sigma), k/gamma);
			printf(" %f %f %f ", norm, l2, l4);
			double error = pow((m2/(sigma*sigma) - 1), 2) + pow(k/gamma - 1, 2) + pow(m0 - 1, 2) + b;
			printf(" b = %f error = %f\n", b, error);
			return error;
		}
		
		double operator()(double x) {
			double xsq = x * x;
			return norm * exp(l2 * xsq  + l4 * xsq * xsq);
		}
		
		double moment(int k) {
			auto f = [this,&k](double x) {
				return this->operator()(x) * pow(x, k);
			};
			IntegratorGSL<> integratorGSL(f);
			double integral = 2 * integratorGSL.integrate(0,xmax);
			return integral;
		}
		
	private:
		void optimize() {
			nlopt_opt opt;
			nlopt_algorithm algorithm=NLOPT_LN_NELDERMEAD;
			opt = nlopt_create(algorithm, 3);
			double x[3] = {1,0.5, 0};
			nlopt_set_min_objective(opt, nlopt_func_mem, this);
			//nlopt_set_xtol_abs(opt, 1e-5);
			nlopt_set_stopval(opt, 1e-5);
			//nlopt_add_inequality_constraint(opt, myconstraint_mem, this, 1e-8);
			//nlopt_set_min_objective(opt, myfunc, NULL);

			double minf = 0;
			if (nlopt_optimize(opt, &x[0], &minf) < 0) {
				printf("nlopt failed!\n");
			}
			else {
				printf("found minimum at f(%g,%g,%g) = %0.10g\n", x[0],x[1],x[2], minf);
			}
		}
	};
};

double nlopt_func_mem(unsigned n, const double *x, double *grad, void *my_func_data)
{
	double norm = x[0];
	double l2 = x[1];
	double l4 = x[2];
	assert(grad == 0);
    return ((gd::MEM*)my_func_data)->fopt(norm, l2, l4);
}

double myconstraint_mem(unsigned n, const double *x, double *grad, void *data)
{
	//double norm = x[0];
	double l2 = x[1];
	double l4 = x[2];
	gd::MEM* mem = (gd::MEM*)data;
	// 2 * xmax * l4 > -l2
	assert(grad == 0);
    return -l2 - 2 * mem->xmax * l4 ;
}
