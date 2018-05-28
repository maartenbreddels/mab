#include "polynomial.hpp"
#include "vectorize.hpp"
#include <boost/python.hpp>
#include <sstream>

namespace gd {
using namespace boost::python;


template<int degree, class Ctor>
void py_export_polynomial_lagrange() {
	typedef Polynomial<1, degree, LagrangeBasis<degree>> P;
	stringstream name;
	name << "PolynomialLagrange" << degree;
	object c = class_< P, boost::noncopyable >(name.str().c_str(), Ctor())
		.def("__call__", &P::operator())
		;
	c.attr("degree") = degree; 
}

template<int N>
void py_export_hermite_basis() {
	typedef HermiteBasis<N> B;
	stringstream name;
	name << "HermiteBasis" << N;
	object c = class_< B, boost::noncopyable >(name.str().c_str(), init<>())
		.def("__call__", vectorize(&B::operator()))
		;
	//c.attr("degree") = degree; 
}

void py_export_polynomial() {
	//py_export_polynomial_lagrange<0>();
	py_export_polynomial_lagrange<0, init<double>>();
	py_export_polynomial_lagrange<1, init<double, double>>();
	py_export_polynomial_lagrange<2, init<double, double, double>>();
	py_export_polynomial_lagrange<3, init<double, double, double, double>>();
	py_export_polynomial_lagrange<4, init<double, double, double, double, double>>();
	py_export_polynomial_lagrange<5, init<double, double, double, double, double, double>>();
	py_export_polynomial_lagrange<6, init<double, double, double, double, double, double, double>>();
	py_export_hermite_basis<0>();
	py_export_hermite_basis<1>();
	py_export_hermite_basis<2>();
	py_export_hermite_basis<3>();
	py_export_hermite_basis<4>();
	py_export_hermite_basis<5>();
	py_export_hermite_basis<6>();
}

}