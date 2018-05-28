#include "transformation.hpp"
#include <boost/python.hpp>
#include <pyublas/numpy.hpp>

namespace gd {
using namespace boost::python;
/*

class Transformation1d_in_3d {
public:
	virtual double transform(double r) = 0;
	virtual double inverse_transform(double r) = 0;
	virtual double d3xdu(double u) = 0;
	virtual double drdu(double u) = 0;
	virtual double laplace_u1_1(double u) = 0;
	virtual double laplace_u1_2(double u) = 0;
	virtual double laplace_u1_1_times_d3xdu(double u) {
		return this->laplace_u1_1(u) * this->d3xdu(u); 
	}
};

*/
template<class T, class Ctor>
void py_export_transformation_type(const char* name) {
	class_<T, bases<Transformation1d_in_3d>, boost::noncopyable  >(name, Ctor())
		.def("transform", (&T::transform))
		.def("inverse_transform", (&T::inverse_transform))
		.def("d3xdu", (&T::d3xdu))
		.def("drdu", (&T::drdu))
		.def("laplace_u1_1", (&T::laplace_u1_1))
		.def("laplace_u1_2", (&T::laplace_u1_2))
		.def("laplace_u1_1_times_d3xdu", (&T::laplace_u1_1_times_d3xdu))
		;
}
void py_export_transformation() {
	class_< Transformation1d_in_3d, boost::noncopyable >("Transformation1d_in_3d", no_init)
		;
	py_export_transformation_type<TransformationSphericalLog_in_3d, init<double>>("TransformationSphericalLog_in_3d");
	py_export_transformation_type<TransformationSphericalTan_in_3d, init<>>("TransformationSphericalTan_in_3d");
	py_export_transformation_type<TransformationSpherical_in_3d, init<>>("TransformationSpherical_in_3d");
}


}
