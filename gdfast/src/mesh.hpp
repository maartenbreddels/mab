#pragma once
#include <Eigen/Core>
#include "datatypes.hpp"
#include <Eigen/Eigen>
#include "profile.hpp"
#include "integration.hpp"
#include "transformation.hpp"

namespace gd {

using namespace std;
USING_PART_OF_NAMESPACE_EIGEN

class Mesh1d {
public:
	virtual int findindex(double) { return 0; }
	virtual double indexto_x(int) { return 0; }
	//virtual void indexto_x(int) { return 0; }
	virtual bool inrange(double) { return false; }
	virtual int length() {return 0;}
};

template<class Base=Mesh1d>
class Mesh1dRegular : public Base {
public:
	Mesh1dRegular(double x1, double x2, int length) : x1(x1), x2(x2), _length(length) {}
	int findindex(double x) {
		return x < x1 ? 0 :
			(x >= x2 ? _length-1 : (int)((x-x1)/(x2-x1)*(_length)));
	}
	double indexto_x(int i) {
		return x1 + (x2-x1)/length()*i;
	}
	bool inrange(double r) {
		return (r >= x1) && (r < x2); 
	}
	int length() { return _length; }
	double x1, x2;
	int _length;
};

template<class T=double>
class BasisTriLeft {
public:
	T operator()(T x) {
		return x >= -1 ? (x < 0 ? (x+1) : 0) : 0;
	}
};

template<class T=double>
class BasisTriRight {
public:
	T operator()(T x) {
		return x >= 0 ? (x < 1 ? (1-x) : 0) : 0;
	}
};

	

template<class T=double>
class BasisTriMesh1dRegular {
public:
	BasisTriMesh1dRegular(T _x1, T _x2, int _n_nodes) : x1(_x1), x2(_x2), n_nodes(_n_nodes), mesh(_x1, _x2, _n_nodes), ctrans(_n_nodes+1,_n_nodes+1) {
		MatrixXd m = MatrixXd::Zero(_n_nodes+1,_n_nodes+1);
		//m.setZero();
		m(0,0) = 1./3;
		m(_n_nodes,_n_nodes) = 1./3;
		for(int i = 1; i < (_n_nodes+1); i++) {
			m(i, i-1) = 1./6;
			m(i-1, i) = 1./6;
		}  
		for(int i = 1; i < (_n_nodes); i++) {
			m(i, i) = 2./3;
		}  
		T scale = (x2-x1)/n_nodes;
		m = m * scale;
		//cout << m << endl;
		ctrans = m.inverse();
		//cout << ctrans << endl;
	}
	void testprofile(Profile* profile, double_vector v, bool dotrans) {
		double* vp = v.data().begin();
		double* array = vp;
		//int size = x.size();
		//double scale = 1./size;
		//cout << "size = " << size << endl;
		for(int i = 0; i < n_nodes; i++) {
			double integral = 0;
			T xleft = mesh.indexto_x(i);
			T xright = mesh.indexto_x(i+1);
			T dx = (xright-xleft);
			{
				BasisTriRight<T> triright;
				auto f = [&](double x) { return triright((x-xleft)/dx) * profile->densityr(x); };
				IntegratorGSL<> integratorGSL(f); // the integrator
				integral = integratorGSL.integrate(xleft, xright);
				array[i] +=  integral;
			}
			{
				BasisTriLeft<T> trileft;
				auto f = [&](double x) { return trileft((x-xleft)/dx-1) * profile->densityr(x); };
				IntegratorGSL<> integratorGSL(f); // the integrator
				integral = integratorGSL.integrate(xleft, xright);
				array[i+1] +=  integral;
			}
		}
		if(dotrans) {
			VectorXd v_alias = VectorXd::Map(v.data().begin(), v.size());
			VectorXd vtrans = ctrans * v_alias;
			VectorXd::Map(v.data().begin(), v.size()) = vtrans;
		}
	}

	void test(double_vector x, double_vector y, double_vector v, bool dotrans) {
		double* xp = x.data().begin();
		double* yp = y.data().begin();
		double* vp = v.data().begin();
		int size = x.size();
		double scale = 1./size * (x2-x1);
		//cout << "size = " << size << endl;
		for(int i = 0; i < size; i++) {
			//cout << "i = " << i << endl;
			this->operator()(xp[i], yp[i]*scale, vp);
			
		}
		if(dotrans) {
			VectorXd v_alias = VectorXd::Map(v.data().begin(), v.size());
			VectorXd vtrans = ctrans * v_alias;
			VectorXd::Map(v.data().begin(), v.size()) = vtrans;
		}
	}


	template<class A>
	T operator()(T x, T y, A array) {
		if(mesh.inrange(x)) {
			int indexleft = mesh.findindex(x);
			T xleft = mesh.indexto_x(indexleft);
			T xright = mesh.indexto_x(indexleft+1);
			T fraction = (x-xleft)/(xright-xleft);
			//cout << "x = " << x << " xleft = " << xleft << " xright = " << xright << " fraction = " << fraction << endl;
			BasisTriLeft<T> trileft;
			BasisTriRight<T> triright;
			//cout << "y = " << y << " triright(fraction) = " << triright(fraction) << " trileft (fraction-1) = " << trileft (fraction-1) << endl;
			array[indexleft]   +=  triright(fraction) * y;
			array[indexleft+1] +=  trileft (fraction-1) * y;
			// TODO: use trileft(fraction) + triright(fraction) == 1
		}
		return 0; // TODO: not finished.., remove code? 
	}
	
	T x1, x2;
	int n_nodes;
	Mesh1dRegular<> mesh;
	MatrixXd ctrans;
};

template<class Basis, class T=double>
class MeshRegularNodal1d {
public:
	typedef Basis basis_type;
	T x1, x2;
	int n_cells;
	Mesh1dRegular<> mesh;
	MatrixXd ctrans;
	Transformation1d_in_3d* transformation;
	int dof;
	enum { dof_per_cell = Basis::degree+1 };

	int get_dof() { return dof;}
	int get_n_cells() { return n_cells;}
	int dof_index(int cell_index, int local_index) {
		return  cell_index*(dof_per_cell-1)+local_index;
	}


	template<int I, class B=Basis>
	struct util {
		typedef util<I-1, typename B::next_type> next_type;
		next_type next;
		T integrate(int i, double xleft, double xright, double dx, Profile* profile) {
			if(i == I) {  
				B basis;
				auto f = [&](double x) { return basis((x-xleft)/dx) * profile->densityr(x); };
				IntegratorGSL<> integratorGSL(f); // the integrator
				return integratorGSL.integrate(xleft, xright);
			} else {
				return next.integrate(i, xleft, xright, dx, profile);
			}
		}
		template<class F>
		T integrate2(int i, double xleft, double xright, double dx, F f) {
			if(i == I) {  
				B basis;
				auto f2 = [&](double x) { return basis((x-xleft)/dx) * f(x); };
				IntegratorGSL<> integratorGSL(f2); // the integrator
				return integratorGSL.integrate(xleft, xright);
			} else {
				return next.integrate2(i, xleft, xright, dx, f);
			}
		}
		template<class Array>
		T eval(Array& array, int index, double xleft, double dx, double x) {
			B basis;
			//cout << "eval: " << index << " " << xleft << " " << dx << " " << x << " u=" << ((x-xleft)/dx) << " " << basis((x-xleft)/dx) << " " << array(index) << endl;
			return basis((x-xleft)/dx) * array(index) + next.eval(array, index-1, xleft, dx, x);
		}
		template<class Array>
		T gradient(Array& array, int index, double xleft, double dx, double x) {
			B basis;
			//cout << "grad: " << index << " " << xleft << " " << dx << " " << x << " u=" << ((x-xleft)/dx) << " " << basis((x-xleft)/dx) << " " << basis.dfdx((x-xleft)/dx) << " " << array(index) << endl;
			return basis.dfdx((x-xleft)/dx)/dx * array(index) + next.gradient(array, index-1, xleft, dx, x);
		}
	};
	template<class B>
	struct util<-1, B> {
		T integrate(int, double, double, double, Profile*) {
			return 0;
		}
		template<class F>
		T integrate2(int, double, double, double, F) {
			return 0;
		}
		template<class Array>
		T eval(Array&, int, double, double, double) {
			return 0;
		}
		template<class Array>
		T gradient(Array&, int, double, double, double) {
			return 0;
		}
	};

	template<int I, int J, class B1=Basis, class B2=Basis>
	struct selfintegrator {
		typedef selfintegrator<I-1, J, typename B1::next_type, B2> next_typeI;
		typedef selfintegrator<I, J-1, B1, typename B2::next_type> next_typeJ;
		next_typeI nextI;
		next_typeJ nextJ;
		T integrate(int i, int j) {
			if((i == I)) { // first search for right i
				if(j == J) { // then right j
					B1 basis1;
					B2 basis2;
					auto f = [&](double x) { return basis1(x) * basis2(x); };
					IntegratorGSL<> integratorGSL(f);
					return integratorGSL.integrate(0, 1);
				}  else {
					return nextJ.integrate(i, j);
				}
			} else {
				return nextI.integrate(i, j);
			}
		}
		T integrate_grad(int i, int j, double xleft, double xright, double dx, Transformation1d_in_3d* transformation) {
			if((i == I)) { // first search for right i
				if(j == J) { // then right j
					B1 basis1;
					B2 basis2;
					//auto f = [&](double x) { return basis1.dfdx((x-xleft)/dx) / dx * basis2.dfdx((x-xleft)/dx) / dx *  x * x;}; //*/; };
					auto f = [&](double x) -> double {
						//double r = tan(u*M_PI/2);
						double u = x;
						//double t = tan(u*M_PI/2);
						//double c = cos(u*M_PI/2);
						//return 16./2 * basis1.dfdx((x-xleft)/dx) / dx * basis2.dfdx((x-xleft)/dx) / dx * t * t * c * c;
						return basis1.dfdx((x-xleft)/dx) / dx * basis2.dfdx((x-xleft)/dx) / dx * transformation->laplace_u1_1_times_d3xdu(u) * transformation->laplace_u1_2(u);
					}; //*/; };
					IntegratorGSL<> integratorGSL(f);
					return integratorGSL.integrate(xleft, xright);
				}  else {
					return nextJ.integrate_grad(i, j, xleft, xright, dx, transformation);
				}
			} else {
				return nextI.integrate_grad(i, j, xleft, xright, dx, transformation);
			}
		}
	};
	template<int I, class B1, class B2>
	struct selfintegrator<I, -1, B1, B2> {
		T integrate(int, int) {
			return 0;
		}
		T integrate_grad(int, int, double, double, double, Transformation1d_in_3d*) {
			return 0;
		}
	};
	template<int J, class B1, class B2>
	struct selfintegrator<-1, J, B1, B2> {
		T integrate(int, int) {
			return 0;
		}
		T integrate_grad(int, int, double, double, double, Transformation1d_in_3d* ) {
			return 0;
		}
	};

	double integrate_gradshape(int cell_index, int i, int j) {
		typedef selfintegrator<Basis::degree, Basis::degree> selfintegrator_type;
		selfintegrator_type si; 
		T xleft = mesh.indexto_x(cell_index);
		T xright = mesh.indexto_x(cell_index+1);
		T dx = (xright-xleft);
		//cout << "integrate: " << cell_index << " i " << i << " "  << xleft << " to " << xright << endl;
		return si.integrate_grad(i, j, xleft, xright, dx, transformation);
	}

	template<class F>
	double integrate_shape(int cell_index, int i, F f) {
		T xleft = mesh.indexto_x(cell_index);
		T xright = mesh.indexto_x(cell_index+1);
		T dx = (xright-xleft);
		util<Basis::degree, Basis> integrator;
		//cout << "integrate shape: " << cell_index << " i " << i << " " << xleft << " to " << xright << endl;
		return integrator.integrate2(i, xleft, xright, dx, f);
		//return -integrator.integrate2(dof_per_cell-1-i, xright, xleft, dx, f);
	}

	double eval(VectorXd& solution, double x) {
		int cell_index = mesh.findindex(x);
		T xleft = mesh.indexto_x(cell_index);
		T xright = mesh.indexto_x(cell_index+1);
		T dx = xright-xleft;
		double v = 0;
		util<Basis::degree, Basis> util;
		if(mesh.inrange(x))
			v = util.eval(solution, dof_index(cell_index, dof_per_cell-1), xleft, dx, x);
		//cout << "eval: " << v << endl;
		return v;
		/*for(int i = 0; i < dof_per_cell; i++) {
			ev
		}*/
	}

	double gradient(VectorXd& solution, double x) {
		int cell_index = mesh.findindex(x);
		T xleft = mesh.indexto_x(cell_index);
		T xright = mesh.indexto_x(cell_index+1);
		T dx = xright-xleft;
		double v = 0;
		util<Basis::degree, Basis> util;
		if(mesh.inrange(x))
			v = util.gradient(solution, dof_index(cell_index, dof_per_cell-1), xleft, dx, x);
		//cout << "grad: " << " " << cell_index << " " << n_cells << " " << v << endl;
		return v;
		/*for(int i = 0; i < dof_per_cell; i++) {
			ev
		}*/
	}


	MeshRegularNodal1d(T _x1, T _x2, int n_cells, Transformation1d_in_3d* transformation) : x1(_x1), x2(_x2), n_cells(n_cells), mesh(_x1, _x2, n_cells), ctrans(1, 1), transformation(transformation) {
		if(Basis::degree == 0) {
			//dof_per_cell = 1;
			dof = n_cells;
		} else {
			dof = 1 + n_cells + (dof_per_cell-2)*n_cells; // 1 dof per border + dofs inside the cel
		}
		//cout << "n_cells = " << n_cells << " dof = " << dof << " dof_per_cell = " << dof_per_cell << endl;
		MatrixXd m = MatrixXd::Zero(dof, dof);
		ctrans.resize(dof, dof);
		T scale = (x2-x1)/n_cells;

		T integrals[dof_per_cell][dof_per_cell];
		for(int j = 0; j < dof_per_cell; j++) {
			for(int k = 0; k < (j+1); k++) {
				typedef selfintegrator<Basis::degree, Basis::degree> selfintegrator_type;
				selfintegrator_type si; 
				double integral = si.integrate(j,k);
				//cout << j << " " << k << " " << integral << endl;
				integrals[j][k] = integral;
				integrals[k][j] = integral;
			}
		}
		for(int i = 0; i < n_cells; i++) {
			for(int j = 0; j < dof_per_cell; j++) {
				for(int k = 0; k < dof_per_cell; k++) {
					int i1 = i*(dof_per_cell-1)+j;
					int i2 = i*(dof_per_cell-1)+k;
					m(i1, i2) = m(i1, i2) + integrals[j][k];
				}
			}
		}
		//cout << m << endl;
		m = m * scale;
		ctrans = m.inverse();
	}
	void testprofile(Profile* profile, double_vector v, bool dotrans) {
		double* vp = v.data().begin();
		assert((int)v.size() == dof);
		double* array = vp;
		//int size = x.size();
		//double scale = 1./size;
		//cout << "size = " << size << endl;
		for(int i = 0; i < n_cells; i++) {
			for(int j = 0; j < dof_per_cell; j++) {
				T xleft = mesh.indexto_x(i);
				T xright = mesh.indexto_x(i+1);
				T dx = (xright-xleft);
				util<Basis::degree, Basis> integrator;
				double integral = integrator.integrate(j, xleft, xright, dx, profile);
				/**/
				cout << "> " << i << " " << j << " " << (i*(dof_per_cell-1)+j) << " " << integral << endl; 
				//array[i*(dof_per_cell-1)+(dof_per_cell-1-j)] += integral;
				array[i*(dof_per_cell-1)+j] += integral;
				
			} 
		}
		if(dotrans) {
			VectorXd v_alias = VectorXd::Map(v.data().begin(), v.size());
			VectorXd vtrans = ctrans * v_alias;
			VectorXd::Map(v.data().begin(), v.size()) = vtrans;
		}
	}

	void test(double_vector x, double_vector y, double_vector v, bool dotrans) {
		/*double* xp = x.data().begin();
		double* yp = y.data().begin();
		double* vp = v.data().begin();
		int size = x.size();
		double scale = 1./size * (x2-x1);
		//cout << "size = " << size << endl;
		for(int i = 0; i < size; i++) {
			//cout << "i = " << i << endl;
			this->operator()(xp[i], yp[i]*scale, vp);
			
		}
		if(dotrans) {
			VectorXd v_alias = VectorXd::Map(v.data().begin(), v.size());
			VectorXd vtrans = ctrans * v_alias;
			VectorXd::Map(v.data().begin(), v.size()) = vtrans;
		}*/
	}


	template<class A>
	T operator()(T x, T y, A array) {
		/*if(mesh.inrange(x)) {
			int indexleft = mesh.findindex(x);
			T xleft = mesh.indexto_x(indexleft);
			T xright = mesh.indexto_x(indexleft+1);
			T fraction = (x-xleft)/(xright-xleft);
			//cout << "x = " << x << " xleft = " << xleft << " xright = " << xright << " fraction = " << fraction << endl;
			BasisTriLeft<T> trileft;
			BasisTriRight<T> triright;
			//cout << "y = " << y << " triright(fraction) = " << triright(fraction) << " trileft (fraction-1) = " << trileft (fraction-1) << endl;
			array[indexleft]   +=  triright(fraction) * y;
			array[indexleft+1] +=  trileft (fraction-1) * y;
			// TODO: use trileft(fraction) + triright(fraction) == 1
		}*/
	}
	
};


}