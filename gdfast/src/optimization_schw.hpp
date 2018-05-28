#include "datatypes.hpp"
#include <Eigen/Core>
#include <Eigen/Eigen>
#include "parameters.hpp"

#include <nlopt.h>

namespace gd {

using namespace Eigen;
using namespace std;

template<class F>
double nlopt_wrapper(unsigned n, const double *x, double *grad, void *data);
template<typename F=std::function<double(unsigned n, const double *x, double *grad)>, typename T=double>
class MinimizerNLoptTest {
	typedef MinimizerNLoptTest<F, T> type;
public:
	F f;
	int N;
	MinimizerNLoptTest(F f, int N) : f(f), N(N) {
	}
	//template<int N, class... Ts>
	void optimize(int n_eval, double* x) {
	//const int N = sizeof...(Ts);
		nlopt_opt opt;
		opt = nlopt_create(NLOPT_LD_MMA, N);
		opt = nlopt_create(NLOPT_LD_LBFGS, N);
		//opt = nlopt_create(NLOPT_AUGLAG, N);
	//opt = nlopt_create(NLOPT_LD_VAR1, N);
	//opt = nlopt_create(NLOPT_LN_COBYLA, N);
		
	//nlopt_set_lower_bounds(opt, lb);
		nlopt_set_min_objective(opt, nlopt_wrapper<type>, this);
		//double lower_bounds[N] = { 0 };
		double* lower_bounds = new double[N];
		double* upper_bounds = new double[N];
//lower_bounds[0] = 0.1;
		//lower_bounds[1] = 1;
		for(int i = 0; i < N; i++) {
			lower_bounds[i] = -HUGE_VAL;
			upper_bounds[i] =  HUGE_VAL;
		}
		
		//double upper_bounds[N];
		//upper_bounds[0] = 10;
		//upper_bounds[1] = 40;
		//for(int i = 2; i < N; i++) {
		//	lower_bounds[i] = -1;
		//	upper_bounds[i] = 1;
		//}
		//nlopt_set_lower_bounds(opt, lower_bounds);
		//nlopt_set_upper_bounds(opt, upper_bounds);
		//for(int i = 0; i < N; i++)
		//	printf("lower_bounds[%d] = %f\n", i, lower_bounds[i]);
		//for(int i = 0; i < N; i++)
		//	printf("upper [%d] = %f\n", i, upper_bounds[i]);
	//double steps[N] = {1e-8, 1e-10};
	//nlopt_set_initial_step1(opt, 1e-5);
		
		nlopt_set_xtol_rel(opt, 1e-8);
		//nlopt_set_ftol_abs(opt, 1e-9);
		nlopt_set_maxeval(opt, n_eval);
		//double x[N] = { initial_values... };
		double minf;
		int ret = nlopt_optimize(opt, x, &minf);
		if (ret < 0) {
			printf("nlopt failed!: %d\n", ret);
		}
		else {
			printf("found minimum at f(...");
			//for(int i = 0; i < N; i++)
			//	printf("%g, ", x[i]);
			printf(") = %0.10g\n", minf);
		}
		nlopt_destroy(opt);
	}
};

template<class F>
MinimizerNLoptTest<F> nnlopt_optimize(F f, int N)
{
	MinimizerNLoptTest<F> opt(f, N);
	return opt;
}

namespace schw {

};



class OptimizationMatrixChiSquare {
public:
	OptimizationMatrixChiSquare(double_matrix model_matrix, double_vector observed) {
		this->model_matrix = MatrixXd::Map(model_matrix.data().begin(), model_matrix.size2(), model_matrix.size1());
		this->observed = VectorXd::Map(observed.data().begin(), observed.size());
	}
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		//double total = x.sum();
		//VectorXd xn = x/total;
		//double norm = totalmass.dot(x);
		//int N = pmatrix.rows();
		//printf("[N=%d f=%f %f]\n", N, norm, log(norm));
		//double value = ((pmatrix * x).cwise().log()).sum() - N*log(norm);// - pmatrix.rows() * log(x.sum());
		double chisq = (model_matrix * x - observed).cwise().square().sum();
		return -0.5 * chisq;
		//return value;
	}
	void _dlogpdx(double_vector _x, double_vector _gradient) {
		VectorXd x = VectorXd::Map(_x.data().begin(), _x.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dlogpdx(x, gradient);
	}
	
	template<class T1, class T2>
	void dlogpdx(T1& x, T2& gradient) {
		//double total = x.sum();
		//VectorXd xn = x/total;
		VectorXd model_values = model_matrix * x;
		//int N = model_matrix.rows();
		//VectorXd Pi_inv = Pi.cwise().inverse();
		//VectorXd t1 = -pmatrix.transpose() * Pi_inv;
		//VectorXd g = t1.cwise() + x.size()/(x.sum());
		//gradient.setZero();
		int Nx = x.size();
		//printf("matrix: %d %d observed: %d\n", model_matrix.cols(), model_matrix.rows(), observed.size());
		//double norm = totalmass.dot(x);
		//int N = pmatrix.rows();
		int Nj = model_matrix.rows();
		for(int k = 0; k < Nx; k++) {
			for(int j = 0; j < Nj; j++) {
				gradient(k) += -(model_values(j) - observed(j))*model_matrix(j,k);
			}
			//gradient(k) += -N/norm*totalmass(k);
		}
		//gradient += g;//.cwise();//*x;
		
	}
	MatrixXd model_matrix;
	VectorXd observed;
};



class OptimizationMatrix {
public:
	OptimizationMatrix(double_matrix pmatrix, double_vector totalmass) {
		this->pmatrix = MatrixXd::Map(pmatrix.data().begin(), pmatrix.size2(), pmatrix.size1());
		this->totalmass = VectorXd::Map(totalmass.data().begin(), totalmass.size());
	}
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		//double total = x.sum();
		//VectorXd xn = x/total;
		double norm = totalmass.dot(x);
		int N = pmatrix.rows();
		//printf("[N=%d f=%f %f]\n", N, norm, log(norm));
		double value = ((pmatrix * x).cwise().log()).sum() - N*log(norm);// - pmatrix.rows() * log(x.sum());
		return value;
	}
	void _dlogpdx(double_vector _x, double_vector _gradient) {
		VectorXd x = VectorXd::Map(_x.data().begin(), _x.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dlogpdx(x, gradient);
	}
	
	template<class T1, class T2>
	void dlogpdx(T1& x, T2& gradient) {
		//double total = x.sum();
		//VectorXd xn = x/total;
		VectorXd Pi = (pmatrix * x);
		int N = pmatrix.rows();
		//VectorXd Pi_inv = Pi.cwise().inverse();
		//VectorXd t1 = -pmatrix.transpose() * Pi_inv;
		//VectorXd g = t1.cwise() + x.size()/(x.sum());
		//gradient.setZero();
		int Nx = x.size();
		//printf("matrix: %d %d Pi: %d\n", pmatrix.cols(), pmatrix.rows(), Pi.size());
		double norm = totalmass.dot(x);
		//int N = pmatrix.rows();
		int Nj = pmatrix.rows();
		for(int k = 0; k < Nx; k++) {
			for(int j = 0; j < Nj; j++) {
				gradient(k) += pmatrix(j,k)/Pi(j);
			}
			gradient(k) += -N/norm*totalmass(k);
		}
		//gradient += g;//.cwise();//*x;
		
	}
	MatrixXd pmatrix;
	VectorXd totalmass;
};

class OptimizationMatrixForegroundConditional {
public:
	OptimizationMatrixForegroundConditional(double_matrix pmatrix, double_matrix pcmatrix, double_vector ratios, double_vector p_v_non_member) {
		this->pmatrix = MatrixXd::Map(pmatrix.data().begin(), pmatrix.size2(), pmatrix.size1());
		this->pcmatrix = MatrixXd::Map(pcmatrix.data().begin(), pcmatrix.size2(), pcmatrix.size1());
		//this->totalmass = VectorXd::Map(totalmass.data().begin(), totalmass.size());
		this->ratios = VectorXd::Map(ratios.data().begin(), ratios.size());
		this->p_v_non_member = VectorXd::Map(p_v_non_member.data().begin(), p_v_non_member.size());
		p_member = this->ratios.cwise() / (1 + this->ratios.cwise());
	}
	
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		//double total = x.sum();
		//VectorXd xn = x/total;
		//double norm = totalmass.dot(x);
		//int N = pmatrix.rows();
		//printf("[N=%d f=%f %f]\n", N, norm, log(norm));
		// sum( p(v|R) = sum( p(v,R)/p(R) )
		VectorXd p_v_member = (pmatrix * x).cwise() / (pcmatrix * x);
		double value = (p_v_member.cwise()*p_member +p_v_non_member.cwise() * ((p_member*-1).cwise() + 1)).cwise().log().sum(); // - N*log(norm);// - pmatrix.rows() * log(x.sum());
		return value;
	}
	void _dlogpdx(double_vector _x, double_vector _gradient) {
		VectorXd x = VectorXd::Map(_x.data().begin(), _x.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dlogpdx(x, gradient);
	}
	
	template<class T1, class T2>
	void dlogpdx(T1& x, T2& gradient) {
		//double total = x.sum();
		//VectorXd xn = x/total;
		VectorXd Pi = (pmatrix * x);
		VectorXd PCi = (pcmatrix * x);
		//VectorXd p = Pi.cwise() / PCi;
		
		VectorXd p_v_member = (pmatrix * x).cwise() / (pcmatrix * x);
		VectorXd p = (p_v_member.cwise()*p_member +p_v_non_member.cwise() * ((p_member*-1).cwise() + 1));
		
		//int N = pmatrix.rows();
		//VectorXd Pi_inv = Pi.cwise().inverse();
		//VectorXd t1 = -pmatrix.transpose() * Pi_inv;
		//VectorXd g = t1.cwise() + x.size()/(x.sum());
		//gradient.setZero();
		int Nx = x.size();
		//printf("matrix: %d %d Pi: %d\n", pmatrix.cols(), pmatrix.rows(), Pi.size());
		//double norm = totalmass.dot(x);
		//int N = pmatrix.rows();
		int Nj = pmatrix.rows();
		//printf(" %d %d %d %d %d %d \n", Nj, Nx, Pi.size(), PCi.size(), p.size(), pcmatrix.cols());
		for(int k = 0; k < Nx; k++) {
			for(int j = 0; j < Nj; j++) {
				//gradient(k) += pmatrix(j,k)/Pi(j);
				//gradient(k) += 1./p(j) * (pmatrix(j,k) / PCi(j)  - Pi(j)/pow(PCi(j), 2) * pcmatrix(j,k) );
				gradient(k) += 1./p(j) * (pmatrix(j,k) * p_member(j) / PCi(j)  - Pi(j) * p_member(j) /pow(PCi(j), 2) * pcmatrix(j,k) );
			}
			//gradient(k) += -N/norm*totalmass(k);
		}
		//exit(0);
		//gradient += g;//.cwise();//*x;
		
	}
	MatrixXd pmatrix, pcmatrix;
	VectorXd totalmass;
	VectorXd ratios;
	VectorXd p_member;
	VectorXd p_v_non_member;
	
};



class OptimizationQP {
public:
	MatrixXd P;
	VectorXd q;
	OptimizationQP(double_matrix P, double_vector q) {
		this->P = MatrixXd::Map(P.data().begin(), P.size2(), P.size1());
		this->q = VectorXd::Map(q.data().begin(), q.size());
	}
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		//double total = x.sum();
		VectorXd Px = P * x;
		double value = x.dot(Px)/2 + q.dot(x);
		//double norm = totalmass.dot(x);
		//double value = ((pmatrix * x).cwise().log().cwise() * counts).sum() - this->N*log(norm);// - pmatrix.rows() * log(x.sum());
		return value;
	}
	void _dlogpdx(double_vector _x, double_vector _gradient) {
		VectorXd x = VectorXd::Map(_x.data().begin(), _x.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dlogpdx(x, gradient);
	}
	
	template<class T1, class T2>
	void dlogpdx(T1& x, T2& gradient) {
		gradient = gradient + P * x + q;
	}
};


class OptimizationMatrixN {
public:
	OptimizationMatrixN(double_matrix pmatrix, double_vector counts, double_vector totalmass) {
		this->pmatrix = MatrixXd::Map(pmatrix.data().begin(), pmatrix.size2(), pmatrix.size1());
		this->counts = VectorXd::Map(counts.data().begin(), counts.size());
		this->totalmass = VectorXd::Map(totalmass.data().begin(), totalmass.size());
		this->N = this->counts.sum();
	}
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		//double total = x.sum();
		//VectorXd xn = x/total;
		double norm = totalmass.dot(x);
		double value = ((pmatrix * x).cwise().log().cwise() * counts).sum() - this->N*log(norm);// - pmatrix.rows() * log(x.sum());
		return value;
	}
	void _dlogpdx(double_vector _x, double_vector _gradient) {
		VectorXd x = VectorXd::Map(_x.data().begin(), _x.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dlogpdx(x, gradient);
	}
	
	template<class T1, class T2>
	void dlogpdx(T1& x, T2& gradient) {
		/*double total = x.sum();
		VectorXd xn = x/total;
		VectorXd pi = (pmatrix * xn);
		VectorXd t = pi.cwise().inverse();
		MatrixXd m = pmatrix;
		for(int i = 0; i < pmatrix.cols(); i++) {
			m.col(i) = m.col(i).cwise() * counts;
		}
		//cout << counts;
		//printf("matrix: %d %d vector: %d\n", pmatrix.cols(), pmatrix.rows(), counts.size());
		//MatrixXd m = pmatrix.transpose().cwise() * counts;
		//VectorXd g1 = -(pmatrix.transpose().colwise() * counts) * t;
		//g1 = g1.cwise() + pmatrix.rows() * 1/(x.sum());
		//g1 = (g1.cwise() * counts)  + (counts*(1/(x.sum())));
		//g1 = (g1.cwise() * counts)  + (counts*(1/(x.sum())));
		//g1 = (g1) ;//  + (counts*(1/(x.sum())));
		//gradient.setZero();
		VectorXd ex = counts*(1./x.sum());
		double bla = counts.sum()/(x.sum());
		gradient += (-(m.transpose() * t)).cwise() + bla;
		//gradient = g1;//.cwise();//*x;*/
		//double total = x.sum();
		//VectorXd xn = x/total;
		VectorXd Pi = (pmatrix * x);
		//VectorXd Pi_inv = Pi.cwise().inverse();
		//VectorXd t1 = -pmatrix.transpose() * Pi_inv;
		//VectorXd g = t1.cwise() + x.size()/(x.sum());
		//gradient.setZero();
		int Nx = x.size();
		double norm = totalmass.dot(x);
		//printf("matrixN: %d %d Pi: %d\n", pmatrix.cols(), pmatrix.rows(), Pi.size());
		int Nj = pmatrix.rows();
		for(int k = 0; k < Nx; k++) {
			for(int j = 0; j < Nj; j++) {
				gradient(k) += counts(j) * pmatrix(j,k)/Pi(j) - counts(j)/norm*totalmass(k);
			}
		}
		//gradient += g;//.cwise();//*x;
		
	}
	MatrixXd pmatrix;
	VectorXd counts;
	VectorXd totalmass;
	double N;
};




class OptimizationNormalize {
public:
	double value, error;
	OptimizationNormalize(double value, double error) : value(value), error(error) {
	}
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		double value = -pow((this->value-x.sum())/error, 2);
		return value;
	}
	void _dlogpdx(double_vector _x, double_vector _gradient) {
		VectorXd x = VectorXd::Map(_x.data().begin(), _x.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dlogpdx(x, gradient);
	}
	
	template<class T1, class T2>
	void dlogpdx(T1& x, T2& gradient) {
		gradient = gradient.cwise() + 2*(value-x.sum())/pow(error,2);
	}	
};

class OptimizationNormalizeMass {
public:
	VectorXd mass_vector;
	double totalmass, error;
	OptimizationNormalizeMass(double_vector mass_vector, double totalmass, double error) : totalmass(totalmass), error(error) {
		this->mass_vector = VectorXd::Map(mass_vector.data().begin(), mass_vector.size());
	}
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		double mass = mass_vector.dot(x);
		double value = -pow((this->totalmass-mass)/error, 2);
		return value;
	}
	void _dlogpdx(double_vector _x, double_vector _gradient) {
		VectorXd x = VectorXd::Map(_x.data().begin(), _x.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dlogpdx(x, gradient);
	}
	
	template<class T1, class T2>
	void dlogpdx(T1& x, T2& gradient) {
		//gradient = gradient.cwise() + 2*(((mass_vector.cwise() - this->totalmass).cwise()) /pow(error,2)).cwise();
		int Nx = x.size();
		//double norm = totalmass.dot(x);
		//printf("matrixN: %d %d Pi: %d\n", pmatrix.cols(), pmatrix.rows(), Pi.size());
		//int Nj = pmatrix.rows();
		double mass = mass_vector.dot(x);
		double value = (this->totalmass-mass);
		double err = pow(error, 2);
		for(int k = 0; k < Nx; k++) {
			//for(int j = 0; j < Nj; j++) {
			gradient(k) += 2 * value/err * mass_vector(k);
			//}
		}
	}	
};

class OptimizationEntropy {
public:
	double scale;
	OptimizationEntropy(double scale) : scale(scale)  {
	}
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		double value = -scale * (x.cwise() * x.cwise().log()).sum();
		return value;
	}
	void _dlogpdx(double_vector _x, double_vector _gradient) {
		VectorXd x = VectorXd::Map(_x.data().begin(), _x.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dlogpdx(x, gradient);
	}
	
	template<class T1, class T2>
	void dlogpdx(T1& x, T2& gradient) {
		//gradient = gradient.cwise() + 2*(value-x.sum())/pow(error,2);
		gradient += -scale * ((x.cwise().log().cwise()+1));
	}	
};



/*
class OptimizationMatrixN {
public:
	OptimizationMatrixN(double_matrix pmatrix, double_vector counts) : pmatrix(pmatrix), counts(counts) {
		this->pmatrix = MatrixXd::Map(pmatrix.data().begin(), pmatrix.size2(), pmatrix.size1());
		this->counts = VectorXd::Map(counts.data().begin(), counts.size());
	}
	double _logp(double_vector _x) {
		Map<VectorXd> x = VectorXd::Map(_x.data().begin(), _x.size());
		return logp(x);
	}
	template<class T>
	double logp(T& x) {
		VectorXd x = u.cwise().exp();
		double total = x.sum();
		VectorXd xn = x/total;
		//VectorXd rho2d = rho2dmatrix * x;
		//cout << "cols " << pmatrix.cols() << endl;
		//double logLkinematics = (pmatrix * x).cwise().log().sum() - pmatrix.rows() * log(x.sum());
		double value = ((pmatrix * xn).cwise().log().cwise() * count.cwise()).sum();// - pmatrix.rows() * log(x.sum());
		return value;
	}
	MatrixXd pmatrix;
	VectorXd counts;
};*/

class OptimizationProblemSchw {
public:
	OptimizationProblemSchw(double_matrix _pmatrix, double_matrix _rho2dmatrix, double_vector _orbitweights, double_vector _rho2d_target, double_vector _rho2d_error, double error_x, double entropy_scale, bool kin, bool light, bool norm) : error_x(error_x), entropy_scale(entropy_scale), kin(kin), light(light), norm(norm) {
		pmatrix = MatrixXd::Map(_pmatrix.data().begin(), _pmatrix.size2(), _pmatrix.size1());
		rho2dmatrix = MatrixXd::Map(_rho2dmatrix.data().begin(), _rho2dmatrix.size2(), _rho2dmatrix.size1());
		orbitweights = VectorXd::Map(_orbitweights.data().begin(), _orbitweights.size());
		rho2d_target = VectorXd::Map(_rho2d_target.data().begin(), _rho2d_target.size());
		rho2d_error = VectorXd::Map(_rho2d_error.data().begin(), _rho2d_error.size());
	}
	
	void optimize(int n_eval, int steps, double_vector _u) {
		int iteration = 0;
		//int steps = 1;
		auto g = [&](int n, const double *x, double *grad) -> double {
			//toy.scale = x[0];
			//toy.M = exp(x[1]);
			Map<VectorXd> u(x, n);
			double minlogL = -this->likelihood(u);
			if(grad) {
				Map<VectorXd> gradvector(grad, n);
				this->dfdx(u, gradvector);
			}
			if((iteration % steps) == 0) { 
				printf("current point(...");
				//for(int i = 0; i < n; i++)
				//	printf("%g, ", x[i]);
				printf(") = %20f\n", minlogL);
			}
			/*
			double chisq = fit.optimize_function(grad);
			if((iteration % steps) == 0) { 
				if(grad) {
					printf(" gradient(");
					for(int i = 0; i < n; i++)
						printf("%g, ", grad[i]);
				}
				printf(")");
				printf(" chisq = %f iteration = %d\n", chisq, iteration);
			}
		//grad[0] = fit.denergyChisqdb(J1, J2);
		//return fit.energyChisq(J1, J2);*/
			iteration++;
			return minlogL;
		};
		auto o = nnlopt_optimize(g, orbitweights.rows());
		//int n = atoi(argv[1]);
		o.optimize(n_eval, _u.data().begin());
		cout << "interations: " << iteration << endl;
	}
	
	void _hessian(double_matrix _h, double_vector _u, bool du=true) {
		Map<MatrixXd> h = MatrixXd::Map(_h.data().begin(), _h.size2(), _h.size1());
		VectorXd u = VectorXd::Map(_u.data().begin(), _u.size());
		hessian(h, u, du);
	}
	
	template<class T1, class T2>
	double hessian(T1& h, T2& u, bool du=true) {
		VectorXd x = u.cwise().exp();
		
		assert(h.rows() == h.cols());
		int N = h.rows();
		VectorXd gradient = VectorXd::Zero(N);
		this->dfdx(u, gradient);
		h.setZero();
		//cout << "g " << gradient(0) << endl;
		
		
		
		
		cout << "rows = " << pmatrix.rows() << endl;
		cout << "cols = " << pmatrix.cols() << endl;
		if(kin) {
			VectorXd w = (pmatrix * x).cwise().inverse().cwise().square();
			cout << "els = " << w.rows() << endl;
			for(int k = 0; k < N; k++) {
				for(int l = 0; l < N; l++) {
					for(int i = 0; i < pmatrix.rows(); i++) {
						if(du) 
							h(l,k) += w(i) * pmatrix(i,l) * pmatrix(i,k) * x(k) * x(l);
						else
							h(l,k) += w(i) * pmatrix(i,l) * pmatrix(i,k);
					}
				}
			}
			
		}
		for(int k = 0; k < N; k++) {
			for(int l = 0; l < N; l++) {
				if((k == l) and du)
					h(k,l) += gradient(k);
				if(light) {
					for(int i = 0; i < rho2dmatrix.rows(); i++) {
						if(du)
							h(k,l) += 2 * rho2dmatrix(i,l) * rho2dmatrix(i,k) / pow(rho2d_error(i), 2) * x(k) * x(l);
						else
							h(k,l) += 2 * rho2dmatrix(i,l) * rho2dmatrix(i,k) / pow(rho2d_error(i), 2);
					}
				}
				if(norm) {
					if(du)
						h(k, l) += 2/pow(error_x, 2) * x(k) * x(l);
					else
						h(k, l) += 2/pow(error_x, 2);
				}
				//cout << "h " << h(k,l) << endl;
			} 
		} 
		//MatrixXd::Map(_h.data().begin(), _h.size2(), _h.size1()) = h;
	}

	double  _likelihood(double_vector _u) {
		Map<VectorXd> u = VectorXd::Map(_u.data().begin(), _u.size());
		return likelihood(u);
	}
	// return f(x) to optimize
	template<class T>
	double likelihood(T& u) {
		//VectorXd u = VectorXd::Map(_u.data().begin(), _u.size());
		VectorXd x = u.cwise().exp();
		double total = x.sum();
		VectorXd xn = x/total;
		VectorXd rho2d = rho2dmatrix * x;
		//cout << "cols " << pmatrix.cols() << endl;
		//double logLkinematics = (pmatrix * x).cwise().log().sum() - pmatrix.rows() * log(x.sum());
		double logLkinematics = (pmatrix * xn).cwise().log().sum();// - pmatrix.rows() * log(x.sum());
		// + 52082.8 + 1500;
		double logLnorm = -pow((1.000001-x.sum())/error_x, 2);
		double logLdensity = -((rho2d_target-rho2d).cwise()/rho2d_error).cwise().pow(2).sum();
		//cout << "u = " << u << endl;
		//cout << "x = " << x << endl;
		//cout << ">" << logLkinematics << " " << logLnorm << " " << logLdensity << endl;
		//assert(0);
		//double entropy = -k * (x.cwise() * x.cwise().log()).sum();
		//return logLkinematics + logLnorm + logLdensity;
		double logL = 0;
		logL += -entropy_scale * (x.cwise() * x.cwise().log()).sum();
		if(kin)
			logL += logLkinematics;
		if(light)
			logL += logLdensity;
		if(norm)
			logL += logLnorm ;
		return logL;
		
		//return logLdensity;// - entropy; // logL \propto -entropy.. ?
		//return 1; //logLnorm;// - entropy; // logL \propto -entropy.. ?
		//return logLdensity;
		//return 0;
		//return logLkinematics;
		//-sum((rho2d_true-rho2d)**2/rho2d_error**2)*fudgefactor
	}
	
	// gradient
	void _dfdx(double_vector _u, double_vector _gradient) {
		VectorXd u = VectorXd::Map(_u.data().begin(), _u.size());
		Map<VectorXd> gradient = VectorXd::Map(_gradient.data().begin(), _gradient.size());
		dfdx(u, gradient);
	}
	
	template<class T1, class T2>
	void dfdx(T1& u, T2& gradient) {
		VectorXd x = u.cwise().exp();
		VectorXd rho2d = rho2dmatrix * x;
		/*g1 =  -sum(dpidxk / pi, axis=1)
		gu = -(2*(1.001-sum(u))/error_x**2)
		g2 = -sum(2*(rho2d_true-rho2d)/rho2d_error**2 *rho2ds, axis=1) * fudgefactor
		g = g1*u+gu*u# -g2 #-(-g1+g2)
		g += g2*u*/
		VectorXd pi = (pmatrix * x);
		//pmatrix.cwise() / pi;
		VectorXd t = pi.cwise().inverse();
		VectorXd g1 = -pmatrix.transpose() * t;
		g1 = g1.cwise() + pmatrix.rows() * 1/(x.sum()); // * x;
		VectorXd gu = -2*(1.000001-x.sum())/pow(error_x,2) * x;
		t = (rho2d_target-rho2d).cwise()/rho2d_error.cwise().pow(2) * -2;
		VectorXd g2 = rho2dmatrix.transpose() * t;
		//VectorXd g(g2.size());
		gradient.setZero();
		gradient = entropy_scale * ((x.cwise().log().cwise()+1).cwise() * x);
		if(kin)
			gradient += g1.cwise()*x;
		if(light)
			gradient += g2.cwise()*x;
		if(norm)
			gradient += gu;

		//VectorXd::Map(_gradient.data().begin(), _gradient.size()) = g;
		//VectorXd::Map(_gradient.data().begin(), _gradient.size()) = g2.cwise()*x;
		
		/*VectorXd::Map(_g1.data().begin(), _g1.size()) = g1.cwise()*x;
		VectorXd::Map(_gu.data().begin(), _gu.size()) = gu;
		VectorXd::Map(_g2.data().begin(), _g2.size()) = g2.cwise()*x;*/
	}
	
	double f_and_gradient(double_vector x, double &fx, double_vector &gradientvector) {
		return 0.;
	}
	MatrixXd pmatrix;
	MatrixXd rho2dmatrix;
	VectorXd orbitweights;
	VectorXd rho2d_target;
	VectorXd rho2d_error;
	double error_x;
	double entropy_scale;
	double k;
	bool kin, light, norm;
};
template<class F>
double nlopt_wrapper(unsigned n, const double *x, double *grad, void *data) {
	F *opt = (F*)data;
	return opt->f(n, x, grad);
}


};