#include <boost/python.hpp>
#include "sparse_grid.hpp"
#include <string>
#include "vectorize.hpp"

namespace gd {
using namespace boost::python;
using namespace std;


template<class T, class Ctor>
void py_export_tri_basis_level(string name, Ctor ctor) {
	class_<T>(name.c_str(), ctor)
		.def("__call__", &T::operator())
		;
}

template<class T, typename... X>
void wrap_eval(T& t, object callback) {
	auto f = [&](X... x) -> double {
		object result = callback(x...);
		extract<double> value(result);
		return value();
	};
	t.eval(f); 
}

void py_export_sparse_grid() {
	/*py_export_tri_basis_level<BasisTriHierLevel<0>>("BasisTriHierLevel0", init<>());
	py_export_tri_basis_level<BasisTriHierLevel<1>>("BasisTriHierLevel1", init<>());
	py_export_tri_basis_level<BasisTriHierLevel<2>>("BasisTriHierLevel2", init<>());
	py_export_tri_basis_level<BasisTriHierLevel<3>>("BasisTriHierLevel3", init<>());*/
	/*{
		typedef Vector<BasisSetTriHierLevel, 3> T;
		class_<T>("v", init<>())
			//.def("__call__", vectorize(&T::operator()))
			//.def("test", &T::test)
			;
	}*/
	{
		typedef Product<2, BasisSetTriHierLevel, 0> T;
		class_<T>("P2", init<>())
			.def("print_", &T::print)
			.def("eval", &wrap_eval<T, double, double>)
//.def("test", &T::test)
			;
	}
	{
		typedef Product<3, BasisSetTriHierLevel, 3> T;
		class_<T>("P3", init<>())
			.def("print_", &T::print)
			//.def("test", &T::test)
			;
	}
	{
		typedef SparseGrid<2, BasisSetTriHierLevel, 3> T;
		class_<T>("SparseGrid2", init<>())
			.def("eval", &wrap_eval<T, double, double>)
			.def("test", &T::test2d)
			.def("__call__", &T::operator()<double, double>)
//.def("test", &T::test)
			;
	}
	{
		typedef SparseGrid<3, BasisSetTriHierLevel, 3> T;
		class_<T>("SparseGrid3", init<>())
			.def("eval", &wrap_eval<T, double, double, double>)
			//.def("test", &T::test2d)
			.def("__call__", &T::operator()<double, double, double>)
//.def("test", &T::test)
			;
	}
	{
		typedef SparseGrid<1, BasisSetTriHierLevel, 5> T;
		class_<T>("SparseGrid1", init<>())
			.def("eval", &wrap_eval<T, double>)
			.def("test", &T::test1d)
			.def("__call__", vectorize(&T::operator()<double>))
//.def("test", &T::test)
			;
	}
	/*
	typedef SparseGrid<1, BasisTriHierLevel, 3> T;
	class_<T>("SparseTest", init<>())
		.def("__call__", vectorize(&T::operator()))
		.def("test", &T::test)
		;
	typedef SparseGrid<2, BasisTriHierLevel, 2> T2;
	class_<T2>("SparseTest2", init<>())
		.def("__call__", (&T2::operator()))
		.def("test", &T2::test)
		;
	*/
}

} // namespace gd
