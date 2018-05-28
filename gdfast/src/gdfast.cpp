// Boost Includes ==============================================================
//#define PYUBLAS_HAVE_BOOST_BINDINGS
#include <boost/python.hpp>
#include <pyublas/numpy.hpp>
#include <gsl/gsl_errno.h>

#include <numpy/arrayobject.h>
//#include <ndarray/boost/python/ndarray/Array.hpp>
//#include <ndarray/boost/python/ndarray/Vector.hpp>
#include <stdexcept>
/*
ndarray::Array<double,1,1> gaussian(int size, double sigma, double mu) {
	ndarray::Array<double,1,1> result = ndarray::allocate(ndarray::makeVector(size));
	ndarray::Array<double,1,1>::Iterator i = result.begin();
	for (int n = 0; n < size; ++n, ++i) {
		double x = (n - mu) / sigma;
		*i = std::exp(-0.5 * x * x);
	}
	return result;
}*/

//#include <boost/python/numeric.hpp>
//#include <boost/cstdint.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
namespace mab {
	void py_export_mab_parameters();
}
namespace gd {
	void py_export_mem();
	void py_export_profile();
	void py_export_profile_axi();
	void py_export_profile_2d();
	void py_export_profile_3d();
	void py_export_jeans();
	void py_export_orbit_integrator();
	void py_export_galaxy();
	void py_export_aperture();
	void py_export_functions();
	//void py_export_poisson_fem();
	void py_export_mesh();
	void py_export_polynomial();
	void py_export_poisson_fem2();
	void py_export_transformation();
	void py_export_mesh2();
//void py_export_google_profiler();
	void py_export_coordinate_systems();
	void py_export_sparse_grid();
	void py_export_sparse_grid_dynamic();
	void py_export_df();
	void py_export_torus_spherical();
	void py_export_tri_orbitreader();
}
// Using =======================================================================
using namespace boost::python;
using namespace gd;

// Declarations ================================================================
namespace  {


}// namespace 

void bla()
{
namespace bp = boost::python;
Py_Initialize();
bp::object main = bp::import("__main__");
bp::object global(main.attr("__dict__"));
global["np"] = bp::import("numpy");
//bp::exec("print np.array([1,2,3])", global, global);
//import_array1(false);
}
// Module ======================================================================
BOOST_PYTHON_MODULE(gdfast)
{
	gsl_set_error_handler_off ();
	import_array(); // need to call this at the beginning of the module
	bla();
	//py_export_mem();
	py_export_profile();
	py_export_profile_axi();
	py_export_profile_2d();
	py_export_profile_3d();
	py_export_jeans();
	py_export_orbit_integrator();
	py_export_galaxy();
	py_export_aperture();
	py_export_functions();
	py_export_df();
	//py_export_poisson_fem();
	py_export_mesh();	
	py_export_polynomial();
	py_export_poisson_fem2();
	py_export_transformation();
	py_export_mesh2();
	//py_export_google_profiler();
	py_export_coordinate_systems();
	py_export_sparse_grid();
	//py_export_sparse_grid_dynamic();
	py_export_torus_spherical();
	mab::py_export_mab_parameters();
	//py_export_tri_orbitreader();
	//boost::python::def("gaussian", &gaussian);
	//py_export_df();
	/*boost::python::def(
	                    "gaussian", &gaussian, boost::python::args("size", "sigma", "mu"),
	                    "array = gaussian(size, sigma, mu)\n\n"
	                    "Create a 1-d Gaussian function with scale sigma and offset mu.\n"
	                  );
	*/
	
}
