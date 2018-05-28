#include "profile3d.hpp"
#include <boost/python.hpp>
#include <pyublas/numpy.hpp>
#include "optimization.hpp"
#include "rootfinder.hpp"

namespace gd {
using namespace boost::python;




void py_export_profile_3d() {
	class_< ProfileModel3d, boost::noncopyable >("ProfileModel3d", no_init)
		;
	class_<ProfileModel3d1C, bases<ProfileModel3d>, boost::noncopyable  >("ProfileModel3d1C", init<Profile3d*>())
		.def("density",   (&ProfileModel3d1C::density))
		.def("potential", (&ProfileModel3d1C::potential))
		.def("dphidx", (&ProfileModel3d1C::dphidx))
		.def("dphidy", (&ProfileModel3d1C::dphidy))
		.def("dphidz", (&ProfileModel3d1C::dphidz))
		;
	class_<ProfileModel3d2C, bases<ProfileModel3d>, boost::noncopyable  >("ProfileModel3d2C", init<Profile3d*,Profile3d*>())
		.def("density",   (&ProfileModel3d2C::density))
		.def("potential", (&ProfileModel3d2C::potential))
		.def("dphidx", (&ProfileModel3d2C::dphidx))
		.def("dphidy", (&ProfileModel3d2C::dphidy))
		.def("dphidz", (&ProfileModel3d2C::dphidz))
		;
	class_< Profile3d, boost::noncopyable >("Profile3d", no_init)
		;
	class_< Density3d, boost::noncopyable >("Density3d", no_init)
		;

	class_<NFW3d, bases<Profile3d, Density3d> >("NFW3d", init<double, double, double, double>())
		.def("density",   (&NFW3d::density))
		.def("potential", (&NFW3d::potential))
		.def("dphidx", (&NFW3d::dphidx))
		.def("dphidy", (&NFW3d::dphidy))
		.def("dphidz", (&NFW3d::dphidz))
		.def_readonly("c", &NFW3d::c)
		.def_readonly("mass200", &NFW3d::mass200)
		.def_readonly("r200", &NFW3d::r200)
		.def_readonly("rs", &NFW3d::rs)
		.def_readonly("rho0", &NFW3d::rho0);
	;
	class_<NFW3dTriax, bases<Profile3d, Density3d> >("NFW3dTriax", init<double, double, double, double, double, double, double, double>())
		.def("density",   (&NFW3dTriax::density))
		.def("potential", (&NFW3dTriax::potential))
		.def("dphidx", (&NFW3dTriax::dphidx))
		.def("dphidy", (&NFW3dTriax::dphidy))
		.def("dphidz", (&NFW3dTriax::dphidz))
		.def_readonly("c", &NFW3dTriax::c)
		.def_readonly("mass200", &NFW3dTriax::mass200)
		.def_readonly("r200", &NFW3dTriax::r200)
		.def_readonly("rs", &NFW3dTriax::rs)
		.def_readonly("rho0", &NFW3dTriax::rho0);
	;
	class_<MiyamotoNagai3d, bases<Profile3d, Density3d> >("MiyamotoNagai3d", init<double, double, double, double>())
		.def("density",   (&MiyamotoNagai3d::density))
		.def("potential", (&MiyamotoNagai3d::potential))
		.def("dphidx", (&MiyamotoNagai3d::dphidx))
		.def("dphidy", (&MiyamotoNagai3d::dphidy))
		.def("dphidz", (&MiyamotoNagai3d::dphidz))
		.def_readonly("a", &MiyamotoNagai3d::a)
		.def_readonly("b", &MiyamotoNagai3d::b)
		.def_readonly("mass", &MiyamotoNagai3d::mass)
	;
	
}


}
