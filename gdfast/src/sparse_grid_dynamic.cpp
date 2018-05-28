#include "sparse_grid_dynamic.hpp"
#include <boost/python.hpp>
#include "datatypes.hpp"

namespace gd { 
using namespace boost::python;

template<class T>
	void wrap_eval(T& t, object callback) {
		auto f = [&](int dim, double* point) -> double {
			object result;
			if(dim == 1)
				result = callback(point[0]);
			if(dim == 2)
				result = callback(point[0], point[1]);
			if(dim == 3)
				result = callback(point[0], point[1], point[2]);
			if(dim == 4)
				result = callback(point[0], point[1], point[2], point[3]);
			if(dim == 5)
				result = callback(point[0], point[1], point[2], point[3], point[4]);
			if(dim == 6)
				result = callback(point[0], point[1], point[2], point[3], point[4], point[5]);
			extract<double> value(result);
			return value();
		};
		t.eval(f); 
	}

template<class T, typename... X>
	double wrap_call(T& t, X... x) {
		/*auto f = [&](X... x) -> double {
			object result = callback(x...);
			extract<double> value(result);
			return value();
		};
		t.eval(f);*/
		double xarray[sizeof...(X)] = {x...};
		return t(xarray);
		//for(int i = 0; i << 
	}

template<class T, typename X0, typename... X>
	void wrap_call_vec(T& t, X0 x0, X... x) {
		/*auto f = [&](X... x) -> double {
			object result = callback(x...);
			extract<double> value(result);
			return value();
		};
		t.eval(f);*/
		const int dim = 1+sizeof...(X);
		X0* xarray[dim] = {&x0, &x...};
		int indices[dim] = {0};
		int count[dim];
		double point[dim];
		for(int d = 0; d < dim; d++) {
			count[d] = xarray[d]->size();
			printf("count[%d] = %d\n", d, count[d]);
			point[d] = xarray[d]->operator[](0);
		}
//for(int i = 0; i << 
		bool done = false;
		while(!done) { // loop over indices
			t(point);
			for(int d = dim-1; d >= 0; d--) {
				indices[d]++;
				
				if(indices[d]>=count[d]) {
					indices[d] = 0;
					if(d == 0) done = true;
					else
					point[d] = xarray[d]->operator[](indices[d]);
				} else {
					point[d] = xarray[d]->operator[](indices[d]);
					break;
				}
			}
		}
	}



void py_export_sparse_grid_dynamic() {
{
	/*class_<BasisSet, boost::noncopyable>("BasisSet", no_init);
	class_<BasisSetHierTriangular, bases<BasisSet>, boost::noncopyable >("BasisSetHierTriangular", init<>())
		.def("count", &BasisSetHierTriangular::count) 
		.def("__call__", &BasisSetHierTriangular::operator()) 
		;
	*/
	class_<BasisSetHierTriangular, boost::noncopyable >("BasisSetHierTriangulara", init<>())
		.def("count", &BasisSetHierTriangular::count) 
		.def("__call__", &BasisSetHierTriangular::operator()) 
		;
	class_<BasisSetHierTriangularM, boost::noncopyable >("BasisSetHierTriangular", init<>())
		.def("count", &BasisSetHierTriangularM::count) 
		.def("center", &BasisSetHierTriangularM::center) 
		.def("__call__", &BasisSetHierTriangularM::operator()) 
		;
	
	{
		typedef SparseGrid T;
		class_<T>("SparseGridD", init<int,int>())
			//.def("test", &T::test);
			.def("marginalize1d", &T::marginalize1d)
			.def("marginalize2d", &T::marginalize2d)
			.def("integrate", &T::integrate)
			.def("eval", &wrap_eval<T>)
			.def("__call__", &wrap_call<T, double>)
			.def("__call__", &wrap_call<T, double, double>)
			.def("__call__", &wrap_call<T, double, double, double>)
			.def("__call__", &wrap_call<T, double, double, double, double>)
			.def("__call__", &wrap_call<T, double, double, double, double, double>)
			.def("__call__", &wrap_call<T, double, double, double, double, double, double>)
			.def("__call__", &wrap_call_vec<T, double_vector>)
			.def("__call__", &wrap_call_vec<T, double_vector, double_vector>)
			.def("__call__", &wrap_call_vec<T, double_vector, double_vector, double_vector>)
			.def("__call__", &wrap_call_vec<T, double_vector, double_vector, double_vector, double_vector>)
			.def("__call__", &wrap_call_vec<T, double_vector, double_vector, double_vector, double_vector, double_vector>)
			.def("__call__", &wrap_call_vec<T, double_vector, double_vector, double_vector, double_vector, double_vector, double_vector>)
//.def("print_", &T::print)
			//.def("eval", &wrap_eval<T, double, double>)
	//.def("test", &T::test)
			;
	}
	{
		typedef RegularGrid T;
		class_<T>("RegularGridD", init<int,int>())
			//.def("test", &T::test);
			.def("eval", &wrap_eval<T>)
			.def("__call__", &wrap_call<T, double>)
			.def("__call__", &wrap_call<T, double, double>)
			.def("__call__", &wrap_call<T, double, double, double>)
			.def("__call__", &wrap_call<T, double, double, double, double>)
			.def("__call__", &wrap_call<T, double, double, double, double, double>)
			.def("__call__", &wrap_call<T, double, double, double, double, double, double>)
//.def("print_", &T::print)
			//.def("eval", &wrap_eval<T, double, double>)
	//.def("test", &T::test)
			;
	}
}
}
}