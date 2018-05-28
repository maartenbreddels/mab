#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <nlopt.h>
#include "datatypes.hpp"

namespace mab {
using namespace std;
using namespace gd;

class Dataset {
	
};

class ModelParameters {
public:
	virtual void length() = 0;
};

class Parameter {
public:
	Parameter(double* ptr, string name) : ptr(ptr), name(name) {
	}
	virtual void set(double value) { *ptr = value; }
	virtual double get() { return *ptr; }
	virtual string getName() { return name; }
    bool operator==(Parameter const& x) const { return ptr == x.ptr; }
    bool operator!=(Parameter const& x) const { return ptr != x.ptr; }
private:
	double* ptr;
	string name;
};

class Model {
public:
	virtual int parameter_count() { return parameters.size(); }
	virtual void add_paramaters(vector<Parameter> parameters) {
		this->parameters.insert(this->parameters.end(), parameters.begin(), parameters.end() );
	}
	virtual vector<Parameter> get_parameters() { return parameters; }
	virtual double logL(double x, double sigma_x) = 0;
	virtual double likelihood(double x, double sigma_x)  {
		return exp(logL(x, sigma_x));
	}
	virtual double sample_twister() {
		return this->sample(engine);
	};
	virtual double sample(std::mt19937 &engine) = 0;
	virtual ~Model() {}
protected:
	vector<Parameter> parameters;
	std::mt19937 engine;
};

class Gaussian : public Model {
public:
	Gaussian(double mu, double sigma) : mu(mu), sigma(sigma), param_mu(&this->mu, "mu"), param_sigma(&this->sigma, "sigma") {
		parameters.push_back(param_mu);
		parameters.push_back(param_sigma);
	}
	virtual double logL(double x, double sigma_x) {
		double xp = (mu-x);
		double sigmasq = this->sigma*this->sigma + sigma_x*sigma_x;
		return (-(xp*xp)/(2*sigmasq)) - log(sqrt(sigmasq)*sqrt(2*M_PI));
	}
	Parameter parameter_mu() { return param_mu; }
	Parameter parameter_sigma() { return param_sigma; }
	virtual double sample(std::mt19937 &engine) {
		normal_distribution<double> d(mu, sigma);
		return d(engine);
	}
	double get_mu() { return mu; }
	double get_sigma() { return sigma; }
private:
	double mu, sigma;
	Parameter param_mu;
	Parameter param_sigma;
};

template<class T, class V>
int find_index(T begin, T end, V value) {
	int index = 0;
	for(T i = begin; i != end; i++) {
		if(*i == value)
			return index;
		index++;
	}
}

class ModelSum : public Model {
	public:
	ModelSum(Model* model1, Model* model2, double ratio=1) : model1(model1), model2(model2), ratio(ratio), param_ratio(&this->ratio, "ratio") {
		add_paramaters(model1->get_parameters());
		add_paramaters(model2->get_parameters());
		parameters.push_back(param_ratio);
	}
	virtual double logL(double x, double sigma_x) {
		// w1 + w2 == 1
		// ratio = w1/w2
		ratio = fabs(ratio);
		double w1 =  ratio / (1 + ratio);
		double w2 = 1 / (1 + ratio);
		return log(w1*exp(model1->logL(x, sigma_x)) + w2*exp(model2->logL(x, sigma_x)) );
	}
	virtual double sample(std::mt19937 &engine) {
		uniform_real_distribution<double> d(0., 1.);
		double u = d(engine);
		ratio = fabs(ratio);
		double w1 =  ratio / (1 + ratio);
		if(u < w1) {
			return model1->sample(engine);
		} else {
			return model2->sample(engine);
		}
	}
	double get_ratio() { return ratio; }
	Model *model1, *model2;
	double ratio;
	Parameter param_ratio;
};

template<class T, class V>
int find_index(T list, V value) {
	return find_index(list.begin(), list.end(), value);
}

static double nlopt_ModelFitterWrapper(unsigned n, const double *x, double *grad, void *data);


class ModelFitter {
public:
	struct data_t {
		ModelFitter* fitter;
		double* x;
		double* sigma_x;
		int N;
	};
	ModelFitter(Model* model, bool debug) : model(model), debug(debug) {
		dimension = model->parameter_count();
		opt = nlopt_create(NLOPT_LN_NELDERMEAD, dimension);
		lower_bounds = new double[dimension];
		upper_bounds = new double[dimension];
		for(int i = 0; i < dimension; i++) {
			lower_bounds[i] = -HUGE_VAL;
			upper_bounds[i] = HUGE_VAL;
		}
		values = new double[dimension];
	}
	virtual ~ModelFitter() {
		nlopt_destroy(opt);
		delete lower_bounds;
		delete upper_bounds;
	}
	
	Model* get_model() { 
		return this->model;
	}
	void set_lower_bound(Parameter parameter, double value) {
		int index = find_index(model->get_parameters(), parameter);
		lower_bounds[index] = value;
		nlopt_set_lower_bounds(opt, lower_bounds);
	}
	
	void set_upper_bound(Parameter parameter, double value) {
		int index = find_index(model->get_parameters(), parameter);
		upper_bounds[index] = value;
		nlopt_set_upper_bounds(opt, lower_bounds);
	}
	
	void fit(double_vector x, double_vector sigma_x) {
		fit_(x.data().begin(), sigma_x.data().begin(), x.size());
	}
	void fit_(double* x, double* sigma_x, int N) {
		data_t d = {this, x, sigma_x, N};
		nlopt_set_min_objective(opt, nlopt_ModelFitterWrapper, &d);
		double* parameter_values = new double[dimension];
		auto parameters = model->get_parameters();
		int index = 0;
		for(auto i = parameters.begin(); i != parameters.end(); i++) {
			parameter_values[index] = (*i).get();
			index++;
		}
		double min_logL;
		if (nlopt_optimize(opt, parameter_values, &min_logL) < 0) {
			printf("nlopt failed!\n");
		}
		else {
			printf("nlopt ok!\n");
			set_parameters(parameter_values);
			//printf("found minimum at f(%g) = %0.10g\n", x[0], min_logL);
		}
		fflush(stdout);
		delete[] parameter_values;
	}
	
	void set_parameters(const double *x) {
		auto parameters = model->get_parameters();
		int index = 0;
		for(auto i = parameters.begin(); i != parameters.end(); i++) {
			(*i).set(x[index]);
			index++;
		}
	}
	double f(const double *x, double *grad, data_t* data) {
		double logL = 0;
		set_parameters(x);
		for(int i = 0; i < data->N; i++) {
			logL += model->logL(data->x[i], data->sigma_x[i]);
		}
		if(debug) {
			auto parameters = model->get_parameters();
			printf("parameters: [");
			for(auto i = parameters.begin(); i != parameters.end(); i++) {
				printf("%s=% 8f ", (*i).getName().c_str(), (*i).get());
			}
			printf("] ");
			printf("logL=%g ", logL);
			printf("\n");
		}
		return -logL;
	}
	double* lower_bounds;
	double* upper_bounds;
	double* values;
	int dimension;
	Model* model;
	bool debug;
	nlopt_opt opt;
};

static double nlopt_ModelFitterWrapper(unsigned n, const double *x, double *grad, void *data) {
	ModelFitter::data_t* d = (ModelFitter::data_t*)data;
    return d->fitter->f(x, grad, d);
}


}