// Boost Includes ==============================================================
//#define PYUBLAS_HAVE_BOOST_BINDINGS
#include "datatypes.hpp"
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/cstdint.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "profile.hpp"
#include <Eigen/Core>
USING_PART_OF_NAMESPACE_EIGEN

using namespace boost::python;
using namespace std;

namespace gd {
	class SchwOpt1 {
		public:
			void test() {
				printf("hello\n");
			}
			void test2(Profile * profile) {
				printf("hello: %f\n", profile->densityr(1));
			}
			void test3() {
				Matrix3f m3;
				m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
				Matrix4f m4 = Matrix4f::Identity();
				Vector4i v4(1, 2, 3, 4);

				std::cout << "m3\n" << m3 << "\nm4:\n"
						<< m4 << "\nv4:\n" << v4 << std::endl;

			}
			void test4(double_matrix _A, double_vector _b, double_vector _c) {
				MatrixXd A = MatrixXd::Map(_A.data().begin(), _A.size2(), _A.size1());
				std::cout << A << std::endl;
				std::cout << _A.size1() << std::endl;
				std::cout << _A.size2() << std::endl;
				std::cout << _b.size() << std::endl;
				//VectorXd b = VectorXd::Map(_b.data().begin(), b.size());
				//VectorXd b = Eigen::Map<VectorXd>(_b.data().begin(), 1, b.size());
				VectorXd b = VectorXd::Map(_b.data().begin(), _b.size());
				std::cout << b << std::endl;
				VectorXd c = VectorXd::Map(_c.data().begin(), _c.size());
				//c = A * b;
				VectorXd::Map(_c.data().begin(), _c.size()) = A.transpose() * b;
				
				/*double array[3];
				//VectorXd v = Eigen::Map<Vector3d>(array, 2);
				VectorXd v = VectorXd::Map(array, 2);
				v.fill(10.);
				double data[6] = {1., 2., 3., 4., 5., 6.};
				//Matrix2i mat2x2(data);
				//MatrixXi mat2x2 = Map<Matrix2i>(data);
				MatrixXd mat2x2 = Eigen::Map<MatrixXd>(data, 3, 2);
				cout << v << endl;
				cout << mat2x2 << endl;
				cout << (mat2x2 *v )<< endl;*/
				
			}
	};
	void py_export_schw_opt() {
		class_<SchwOpt1 >("SchwOpt1", init<>())
				.def("test", (&SchwOpt1::test))
				.def("test2", (&SchwOpt1::test2))
				.def("test3", (&SchwOpt1::test3))
				.def("test4", (&SchwOpt1::test4))
				
				;
	}
	void py_export_schw_opt2();
//void py_export_schw_df_grid();
}
// Using =======================================================================
using namespace gd;

BOOST_PYTHON_MODULE(gdfast_schw)
{
	py_export_schw_opt();
	py_export_schw_opt2();
	//py_export_schw_df_grid();
}
