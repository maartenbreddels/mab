#include "functions.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace gd {

using namespace boost::numeric::ublas;

void py_export_functions() {
	using namespace boost::python;
	def("matmul_2d_3d", matmul_2d_3d);
	def("matmul_3d_3d", matmul_3d_3d);
	def("velocity_cartesian_to_spherical", velocity_cartesian_to_spherical);
	def("add_to_matrix_indexed", add_to_matrix_indexed);
}

void add_to_matrix_indexed(double_matrix M, int_vector _i, int_vector _j, double_vector weight) {
	double* Mp = M.data().begin();
	double* wp = weight.data().begin();
	int* ip = _i.data().begin();
	int* jp = _j.data().begin();
	int size1 = M.size1();
	int size2 = M.size2();
	int length = min(min(_i.size(), _j.size()), weight.size());
	for(int k = 0; k < length; k++) {
		int i = *ip++;
		int j = *jp++;
		double w = *wp++;
		if((i < size1) && (j < size2) && (i >= 0) && (j >= 0)) {
			int index = j + i * size2;
			Mp[index] += w;
		}
	}
}

void matmul_2d_3d(double_matrix M, double_vector x, double_vector y, double_vector xt, double_vector yt, double_vector zt) {
	double* Mp = M.data().begin();
	int Mstride = M.size2();
	double* xp = x.data().begin();
	double* yp = y.data().begin();
	double* xtp = xt.data().begin();
	double* ytp = yt.data().begin();
	double* ztp = zt.data().begin();
	int size = x.size();
	for(int i = 0; i < size; i++) {
		//for(int j = 0; j < 3; j++) {
		*xtp++ = (*xp) * Mp[0+Mstride*0] + (*yp) * Mp[0+Mstride*1];
		*ytp++ = (*xp) * Mp[1+Mstride*0] + (*yp) * Mp[1+Mstride*1];
		*ztp++ = (*xp) * Mp[2+Mstride*0] + (*yp) * Mp[2+Mstride*1];
		xp++; yp++;
		//}
	}
} 
void matmul_3d_3d(double_matrix M, double_vector x, double_vector y, double_vector z, double_vector xt, double_vector yt, double_vector zt) {
	double* Mp = M.data().begin();
	int Mstride = M.size2();
	double* xp = x.data().begin();
	double* yp = y.data().begin();
	double* zp = z.data().begin();
	double* xtp = xt.data().begin();
	double* ytp = yt.data().begin();
	double* ztp = zt.data().begin();
	int size = x.size();
	for(int i = 0; i < size; i++) {
		*xtp++ = (*xp) * Mp[0+Mstride*0] + (*yp) * Mp[0+Mstride*1] + (*zp) * Mp[0+Mstride*2];
		*ytp++ = (*xp) * Mp[1+Mstride*0] + (*yp) * Mp[1+Mstride*1] + (*zp) * Mp[1+Mstride*2];
		*ztp++ = (*xp) * Mp[2+Mstride*0] + (*yp) * Mp[2+Mstride*1] + (*zp) * Mp[2+Mstride*2];
		xp++; yp++; zp++;
	}
} 

void velocity_cartesian_to_spherical(double_vector xv, double_vector yv, double_vector zv, double_vector vxv, double_vector vyv, double_vector vzv, double_vector vrv, double_vector vphiv, double_vector vthetav) {
	double *xp = xv.data().begin();
	double *yp = yv.data().begin();
	double *zp = zv.data().begin();
	double *vxp = vxv.data().begin();
	double *vyp = vyv.data().begin();
	double *vzp = vzv.data().begin();
	double *vrp = vrv.data().begin();
	double *vphip = vphiv.data().begin();
	double *vthetap = vthetav.data().begin();
	double *xpend = xv.data().end();
	while(xp != xpend) {
		double x = *xp++;
		double y = *yp++;
		double z = *zp++;
		double vx = *vxp++;
		double vy = *vyp++;
		double vz = *vzp++;
		double r = sqrt(x*x+y*y+z*z);
		double rhosq = (x*x+y*y);
		double rho = sqrt(rhosq);
		//double cosatan2xy = x / rho;
		//double sinatan2xy = y / rho;
		*vrp++ = (vx*x + vy*y + vz*z)/r;
		*vphip++ = (vy*x - vx*y)/rho;
		*vthetap++ = (vx*z*x/rho + vy*z*y/rho - vz*rho)/r;
		
	}
}


}
