#include "google_profiler.hpp"
#include <boost/python.hpp>


namespace gd {
using namespace boost::python;

void py_export_google_profiler() {
	def("ProfilerStart", ProfilerStart);
	def("ProfilerStop", ProfilerStop);
}

}