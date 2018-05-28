#include "mesh.hpp"
#include <boost/python.hpp>
#include <stdio.h>
#include <algorithm>
#include "polynomial.hpp"

namespace gd {

using namespace boost::python;
using namespace boost;
using namespace std;

void py_export_mesh() {
	class_<BasisTriMesh1dRegular<> >("BasisTriMesh1dRegular", init<double, double, int>())
		.def("test", (&BasisTriMesh1dRegular<>::test))
		.def("testprofile", (&BasisTriMesh1dRegular<>::testprofile))
		;
	{
		typedef MeshRegularNodal1d<LagrangeBasis<1>> T;
		class_<T>("MeshRegularNodal1dLagrange1", init<double, double, int, Transformation1d_in_3d*>())
			//.def("test", (&BasisTriMesh1dRegular<>::test))
			.def("testprofile", (&T::testprofile))
			.def("get_dof", (&T::get_dof))
			;
	}
	{
		typedef MeshRegularNodal1d<LagrangeBasis<2>> T;
		class_<T>("MeshRegularNodal1dLagrange2", init<double, double, int, Transformation1d_in_3d*>())
			//.def("test", (&BasisTriMesh1dRegular<>::test))
			.def("testprofile", (&T::testprofile))
			.def("get_dof", (&T::get_dof))
			;
	}
	{
		typedef MeshRegularNodal1d<LagrangeBasis<3>> T;
		class_<T>("MeshRegularNodal1dLagrange3", init<double, double, int, Transformation1d_in_3d*>())
			//.def("test", (&BasisTriMesh1dRegular<>::test))
			.def("testprofile", (&T::testprofile))
			.def("get_dof", (&T::get_dof))
			;
	}
}

}