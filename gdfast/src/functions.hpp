#pragma once
#include "profile.hpp"
#include "datatypes.hpp"

namespace gd {
	
void py_export_functions();

void matmul_2d_3d(double_matrix M, double_vector x, double_vector y, double_vector xt, double_vector yt, double_vector zt);
void matmul_3d_3d(double_matrix M, double_vector x, double_vector y, double_vector z, double_vector xt, double_vector yt, double_vector zt);
/** Converts cartesial velocities to spherical
 @param xv lala
*/
void velocity_cartesian_to_spherical(double_vector xv, double_vector yv, double_vector zv, double_vector vxv, double_vector vyv, double_vector vzv, double_vector vrv, double_vector vphiv, double_vector vthetav);

void add_to_matrix_indexed(double_matrix, int_vector, int_vector, double_vector);
}
