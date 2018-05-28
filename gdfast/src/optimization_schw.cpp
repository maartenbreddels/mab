#include "optimization_schw.hpp"
#include "datatypes.hpp"
#include <boost/python.hpp>

namespace gd {
using namespace boost::python;

void py_export_schw_opt2() {
	class_<OptimizationProblemSchw >("OptimizationProblemSchw", init<double_matrix, double_matrix, double_vector, double_vector, double_vector, double, double, bool, bool, bool>())
		.def("optimize", (&OptimizationProblemSchw::optimize))
			.def("likelihood", (&OptimizationProblemSchw::_likelihood))
			.def("dfdx", (&OptimizationProblemSchw::_dfdx))
			.def("hessian", (&OptimizationProblemSchw::_hessian))
			;
			
	class_<OptimizationMatrixChiSquare>("OptimizationMatrixChiSquare", init<double_matrix, double_vector>())
			.def("logp", (&OptimizationMatrixChiSquare::_logp))
			.def("dlogpdx", (&OptimizationMatrixChiSquare::_dlogpdx))
			;
	class_<OptimizationMatrix>("OptimizationMatrix", init<double_matrix, double_vector>())
			.def("logp", (&OptimizationMatrix::_logp))
			.def("dlogpdx", (&OptimizationMatrix::_dlogpdx))
			;
	class_<OptimizationMatrixForegroundConditional>("OptimizationMatrixForegroundConditional", init<double_matrix, double_matrix, double_vector, double_vector>())
			.def("logp", (&OptimizationMatrixForegroundConditional::_logp))
			.def("dlogpdx", (&OptimizationMatrixForegroundConditional::_dlogpdx))
			;
	class_<OptimizationMatrixN>("OptimizationMatrixN", init<double_matrix, double_vector, double_vector>())
			.def("logp", (&OptimizationMatrixN::_logp))
			.def("dlogpdx", (&OptimizationMatrixN::_dlogpdx))
			;
	class_<OptimizationQP>("OptimizationQP", init<double_matrix, double_vector>())
			.def("logp", (&OptimizationQP::_logp))
			.def("dlogpdx", (&OptimizationQP::_dlogpdx))
			;
	class_<OptimizationNormalize>("OptimizationNormalize", init<double, double>())
			.def("logp", (&OptimizationNormalize::_logp))
			.def("dlogpdx", (&OptimizationNormalize::_dlogpdx))
			;
	class_<OptimizationNormalizeMass>("OptimizationNormalizeMass", init<double_vector, double, double>())
			.def("logp", (&OptimizationNormalizeMass::_logp))
			.def("dlogpdx", (&OptimizationNormalizeMass::_dlogpdx))
			;
	class_<OptimizationEntropy>("OptimizationEntropy", init<double>())
			.def("logp", (&OptimizationEntropy::_logp))
			.def("dlogpdx", (&OptimizationEntropy::_dlogpdx))
			;
}
	

};