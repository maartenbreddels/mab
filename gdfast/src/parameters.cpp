#include "parameters.hpp"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

namespace mab {

using namespace boost::python;

void py_export_mab_parameters() {
	class_<std::vector<Parameter> >("parameter_vector")
        .def(vector_indexing_suite<std::vector<Parameter> >())
    ;	
	class_<Parameter, boost::noncopyable  >("Parameter", 
		init<double*, std::string>((boost::python::arg("ptr"), boost::python::arg("name")))
											)
		.def("set", &Parameter::set)
		.def("get", &Parameter::get)
		.add_property("name", &Parameter::getName)
		;
	
	
	class_<Model, boost::noncopyable  >("Model", no_init)
		.add_property("parameters", &Model::get_parameters)
		.def("sample", &Model::sample_twister)
		;
	class_<Gaussian, bases<Model> >("Gaussian", 
		init<double, double>((boost::python::arg("mu"), boost::python::arg("sigma"))) 
		)
		.def("logL", &Gaussian::logL)
		.def("__call__", &Model::likelihood)
		.add_property("mu", &Gaussian::get_mu)
		.add_property("sigma", &Gaussian::get_sigma)
		;
	class_<ModelSum, bases<Model> >("ModelSum", 
		init<Model*, Model*, double>((boost::python::arg("model1"), boost::python::arg("model2"), boost::python::arg("ratio"))) 
		)
		.def("logL", &ModelSum::logL)
		.def("__call__", &Model::likelihood)
		.add_property("ratio", &ModelSum::get_ratio)
		;


	class_<ModelFitter, boost::noncopyable  >("ModelFitter", 
		init<Model*, bool>((boost::python::arg("model"), boost::python::arg("debug")))
		)
		.def("fit", &ModelFitter::fit)
		.def("get_model", &ModelFitter::get_model, return_internal_reference<1>())
		//.add_property("model", &ModelFitter::get_model)
		;

}


}