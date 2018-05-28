#include "poisson_fem2.hpp"
#include "polynomial.hpp"
#include <boost/python.hpp>
#include "profile_python.hpp"
//#include <sstream>

namespace gd {
using namespace boost::python;

/*
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
*/
template<class MeshType>
void py_export_PoissonSolver1d(string name) {
	//typedef Polynomial<1, degree, LagrangeBasis<degree>> P;
	{
		typedef PoissonSolver1d<MeshType> PoissonSolver1dType;
		stringstream name;
		name << "PoissonSolver1dL" << MeshType::basis_type::degree;
		/*object c =*/
		class_< PoissonSolver1dType, boost::noncopyable >(name.str().c_str(), init<Density*, MeshType*, double>())
			.def("solve", &PoissonSolver1dType::solve)
			.def("__call__", &PoissonSolver1dType::operator())
			.def("gradient", &PoissonSolver1dType::gradient)
			;
	}
	
	{
	
		typedef ProfileNumerical1d<MeshType> ProfileNumerical1dType; 
		//stringstream name;
		//name << "ProfileNumerical1dL" << MeshType::basis_type::degree;
		py_export_profile_profile<ProfileNumerical1dType, init<int, Density*, Transformation1d_in_3d*, double, double, double>>(name);
		
	}
	//c.attr("degree") = degree; 
}

void py_export_poisson_fem2() {
	{
		//typedef MeshRegularNodal1d<LagrangeBasis<1>> MeshType;
		//typedef PoissonSolver1d<MeshType> PoissonSolver1dType; 
		py_export_PoissonSolver1d<MeshRegularNodal1d<LagrangeBasis<1>>>("ProfileNumerical1dL1");
		py_export_PoissonSolver1d<MeshRegularNodal1d<LagrangeBasis<2>>>("ProfileNumerical1dL2");
		py_export_PoissonSolver1d<MeshRegularNodal1d<LagrangeBasis<3>>>("ProfileNumerical1dL3");
		
	}
	/*{
		typedef MeshRegularNodal1d<LagrangeBasis<2>> MeshType;
		typedef PoissonSolver1d<MeshType> PoissonSolver1dType; 
		
		class_< PoissonSolver1dType, boost::noncopyable >("PoissonSolver1dL2", init<Profile*, MeshType*>())
			.def("solve", &PoissonSolver1dType::solve)
			.def("__call__", &PoissonSolver1dType::operator())
			.def("gradient", &PoissonSolver1dType::gradient)
			;
	}
	{
		typedef MeshRegularNodal1d<LagrangeBasis<3>> MeshType;
		typedef PoissonSolver1d<MeshType> PoissonSolver1dType; 
		
		class_< PoissonSolver1dType, boost::noncopyable >("PoissonSolver1dL3", init<Profile*, MeshType*>())
			.def("solve", &PoissonSolver1dType::solve)
			.def("__call__", &PoissonSolver1dType::operator())
			.def("gradient", &PoissonSolver1dType::gradient)
			;
	}*/
}

}
using namespace gd;

extern "C" int main() {
}