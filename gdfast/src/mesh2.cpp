#include "mesh2.hpp"
#include <boost/python.hpp>
#include <stdio.h>
#include <algorithm>
#include "polynomial.hpp"

namespace gd {

using namespace boost::python;
using namespace boost;
using namespace std;

void py_export_mesh2() {
	//class_<BasisTriMesh1dRegular<> >("BasisTriMesh1dRegular", init<double, double, int>())
	//.def("test", (&BasisTriMesh1dRegular<>::test))
	//.def("testprofile", (&BasisTriMesh1dRegular<>::testprofile))
	//;
	{
		typedef MeshRegularNodal<3, LagrangeBasis<0>> T;
		class_<T>("MeshRegularNodal3dLagrange0", init<double, double, double, double, double, double, int, int, int>())
			//.def("test", (&T::test))
			.def("get_dof", (&T::get_dof))
			.def("get_dof_per_cell", (&T::get_dof_per_cell))
			.def("eval", (&T::eval_))
			.def("dof_index", &T::dof_index)
			.def("basis_uv", &T::basis_uv)
			.def("solve_coordinates", &T::solve_coordinates)
			;
	}
	{
		typedef MeshRegularNodal<3, LagrangeBasis<1>> T;
		class_<T>("MeshRegularNodal3dLagrange1", init<double, double, double, double, double, double, int, int, int>())
			//.def("test", (&T::test))
			.def("get_dof", (&T::get_dof))
			.def("get_dof_per_cell", (&T::get_dof_per_cell))
			.def("eval", (&T::eval_))
			.def("dof_index", &T::dof_index)
			.def("basis_uv", &T::basis_uv)
			.def("solve_coordinates", &T::solve_coordinates)
			;
	}
	{
		typedef MeshRegularNodal<2, LagrangeBasis<0>> T;
		class_<T>("MeshRegularNodal2dLagrange0", init<double, double, double, double, int, int>())
			.def("test", (&T::test))
			.def("get_dof", (&T::get_dof))
			.def("get_dof_per_cell", (&T::get_dof_per_cell))
			.def("eval", (&T::eval_))
			.def("dof_index", &T::dof_index)
			.def("basis_uv", &T::basis_uv)
			.def("solve_coordinates", &T::solve_coordinates)
			;
	}
	{
		typedef MeshRegularNodal<2, LagrangeBasis<1>> T;
		class_<T>("MeshRegularNodal2dLagrange1", init<double, double, double, double, int, int>())
			.def("test", (&T::test))
			.def("get_dof", (&T::get_dof))
			.def("get_dof_per_cell", (&T::get_dof_per_cell))
			.def("eval", (&T::eval_))
			.def("dof_index", &T::dof_index)
			.def("basis_uv", &T::basis_uv)
			.def("solve_coordinates", &T::solve_coordinates)
			;
	}
	{
		typedef MeshRegularNodal<2, LagrangeBasis<2>> T;
		class_<T>("MeshRegularNodal2dLagrange2", init<double, double, double, double, int, int>())
			.def("test", (&T::test))
			.def("get_dof", (&T::get_dof))
			.def("get_dof_per_cell", (&T::get_dof_per_cell))
			.def("eval", (&T::eval_))
			.def("dof_index", &T::dof_index)
			.def("basis_uv", &T::basis_uv)
			.def("solve_coordinates", &T::solve_coordinates)
			;
	}
	
	{
		typedef MeshRegularNodal<1, LagrangeBasis<0>> T;
		class_<T>("MeshRegularNodal1dLagrange0", init<double, double,  int>())
			.def("get_dof", (&T::get_dof))
			.def("get_dof_per_cell", (&T::get_dof_per_cell))
			.def("eval", (&T::eval_))
			.def("dof_index", &T::dof_index)
			.def("basis_uv", &T::basis_uv)
			.def("solve_coordinates", &T::solve_coordinates)
			;
	}
	{
		typedef MeshRegularNodal<1, LagrangeBasis<1>> T;
		class_<T>("MeshRegularNodal1dLagrange1", init<double, double, int>())
			.def("get_dof", (&T::get_dof))
			.def("get_dof_per_cell", (&T::get_dof_per_cell))
			.def("eval", (&T::eval_))
			.def("dof_index", &T::dof_index)
			.def("basis_uv", &T::basis_uv)
			.def("solve_coordinates", &T::solve_coordinates)
			;
	}
	{
		typedef MeshRegularNodal<1, LagrangeBasis<2>> T;
		class_<T>("MeshRegularNodal1dLagrange2", init<double, double, int>())
			.def("get_dof", (&T::get_dof))
			.def("get_dof_per_cell", (&T::get_dof_per_cell))
			.def("eval", (&T::eval_))
			.def("dof_index", &T::dof_index)
			.def("basis_uv", &T::basis_uv)
			.def("solve_coordinates", &T::solve_coordinates)
			;
	}	
/*{
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
	}*/
}

}