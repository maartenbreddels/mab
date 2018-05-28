#include "df_grid.hpp"
#include "datatypes.hpp"
#include <boost/python.hpp>

namespace gd {
using namespace boost::python;

void py_export_schw_df_grid() {
	class_<DFGrid, boost::noncopyable >("DFGrid", init<int, int, double, double, Galaxy*>())
		.def("output_grid", &DFGrid::output_grid) 
		.def("n_dofs", &DFGrid::n_dofs) 
		.def("refine_global", &DFGrid::refine_global)
		.def("shape_value", &DFGrid::shape_value)
		.def("print_dof_indices_per_face", &DFGrid::print_dof_indices_per_face)
		.def("print_dof_indices_per_face_on_level", &DFGrid::print_dof_indices_per_face_on_level)
		.def("print_vertices", &DFGrid::print_vertices)
			;
}
	

};