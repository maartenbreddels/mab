#include "galaxy.hpp"
#include <boost/python.hpp>
//#include "quadrature.hpp"

namespace gd {
using namespace boost::python;


void py_export_galaxy() {
	class_< Galaxy, boost::noncopyable >("Galaxy", no_init)
		//.def("densityr", pure_virtual(&Profile::densityr))
	;
	class_<Galaxy1C_constant_anisotropy, bases<Galaxy> >("Galaxy1C_constant_anisotropy", init<Profile*, double, double>())
		.def("findzmax", (&Galaxy1C_constant_anisotropy::findzmax))
		.def("rmax_at_E", (&Galaxy1C_constant_anisotropy::rmax_at_E))
		.def("Tr", (&Galaxy1C_constant_anisotropy::Tr))
	;

	class_<Galaxy2C_constant_anisotropy, bases<Galaxy> >("Galaxy2C_constant_anisotropy", init<Profile*, Profile*, double, double>())
		.def("potentialr", (&Galaxy2C_constant_anisotropy::potentialr))
		.def("findzmax", (&Galaxy2C_constant_anisotropy::findzmax))
		.def("rmax_at_E", (&Galaxy2C_constant_anisotropy::rmax_at_E))
		.def("Tr", (&Galaxy2C_constant_anisotropy::Tr))
	;

}

double Galaxy1C_constant_anisotropy::Tr(double E, double L, double ra, double rp) {
	double Lsq = L*L;
	auto f = [&](double r) {
		return 2/sqrt(-Lsq/r*r-2*this->potentialr(r)+2*E);
	};
	IntegratorGSL<> integratorGSL(f); // the integrator
	//return integratorGSL.integrate(ra, rp, 100000, 0, 1e-4);
	double center = (ra+rp)/2;
	return integratorGSL.integrate_to_sigularity(rp, center) + integratorGSL.integrate_to_sigularity(center, ra);
	
}



}
