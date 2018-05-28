// Boost Includes ==============================================================
//#define PYUBLAS_HAVE_BOOST_BINDINGS
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/cstdint.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Includes ====================================================================
void py_export_profile();
void py_export_jeans();


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
	py_export_profile();
	py_export_jeans();


	//mab::numpy::init();

}
