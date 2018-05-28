// Boost Includes ==============================================================
//#define PYUBLAS_HAVE_BOOST_BINDINGS
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/cstdint.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Includes ====================================================================
#include "profile.hpp"
#include "jeans.hpp"
#include "orbit_integrator.hpp"


//#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
//#include <numpy/arrayobject.h>
//#include <numpyconv.hpp>

// Using =======================================================================
using namespace boost::python;
using namespace gd;

// Declarations ================================================================
namespace  {


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(gdfast)
{
	//import_array();
	//boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

	class_< Profile, boost::noncopyable >("Profile", no_init)
		.def("densityr", pure_virtual(&Profile::densityr))
	;
	class_<Plummer, bases<Profile> >("Plummer", init<double, double, double>())
		.def("densityr", (&Profile::densityr))
		.def("dphidr", (&Profile::dphidr))
		.def("dphidr2", (&Plummer::dphidr2))
	;
	class_<Hernquist, bases<Profile> >("Hernquist", init<double, double, double>())
		.def("densityr", (&Profile::densityr))
		.def("dphidr", (&Profile::dphidr))
	;

	class_< OrbitIntegrator, boost::noncopyable >("OrbitIntegrator", no_init)
	;
	class_<OrbitIntegratorEuler, bases<OrbitIntegrator> >("OrbitIntegratorEuler", init<double, double, double, double>())
		.def("integrate", (&OrbitIntegratorEuler::integrate))
	;
	class_<OrbitIntegratorLeapFrog, bases<OrbitIntegrator> >("OrbitIntegratorLeapFrog", init<double, double, double, double>())
		.def("integrate", (&OrbitIntegratorLeapFrog::integrate))
	;

	//mab::numpy::init();

}
