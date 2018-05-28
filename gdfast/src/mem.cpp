#include "mem.hpp"
#include <boost/python.hpp>

namespace gd {
using namespace boost::python;

void py_export_mem() {
	
	class_< MEM, boost::noncopyable >("MEM", init<double, double, double>( (boost::python::arg("sigma"), boost::python::arg("gamma"), boost::python::arg("xmax"))   ))
		.def("__call__", (&MEM::operator()))
		.def("moment", (&MEM::moment))
		;
}

};