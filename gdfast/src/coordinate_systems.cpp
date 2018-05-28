#include "coordinate_systems.hpp"
#include <boost/python.hpp>
using namespace gd;
using namespace boost::python;

namespace gd {
struct custom_vector_to_python_tuple
{
	static PyObject* convert(Vector3 const& v)
	{
		boost::python::tuple t = make_tuple(v(0), v(1), v(2));
		return boost::python::incref(t.ptr());
	}
};


void py_export_coordinate_systems() {
	
	boost::python::to_python_converter<
		Vector3,
		custom_vector_to_python_tuple>();
	
	class_<CoordinateSystem, boost::noncopyable >("CoordinateSystem", no_init)
	                         ;
	class_<Cartesian, bases<CoordinateSystem> >("Cartesian", init<>())
		;
	class_<Cylindrical, bases<CoordinateSystem> >("Cylindrical", init<>())
		;
	class_<SphericalGalactic, bases<CoordinateSystem> >("SphericalGalactic", init<>())
		.def("ex", &SphericalGalactic::ex)
		.def("ey", &SphericalGalactic::ey)
		.def("ez", &SphericalGalactic::ez)
	;
	
	class_<Coordinate>("Coordinate", init<double,double,double, CoordinateSystem*>())
		.def("to_cartesian", &Coordinate::to_cartesian)
		;
	class_<VelocityCoordinate>("VelocityCoordinate", init<double,double,double, CoordinateSystem*>())
		.def("to_cartesian", &VelocityCoordinate::to_cartesian)
		;
//class_<Position>("Position", init<double,double,double, CoordinateSystem*>())
	//	.def("to_cartesian", &Coordinate::to_cartesian)
	//	;
	
	
	class_<ReferenceFrameBase, boost::noncopyable  >("ReferenceFrameBase", no_init)
		;
	
	class_<Position>("Position", init<Coordinate*, ReferenceFrameBase*>())
		.def("to", &Position::to)
		.def("to_coordinate_system", &Position::to_coordinate_system)
		.def("to_global", &Position::to_global)
		;
	class_<Velocity>("Velocity", init<VelocityCoordinate*, ReferenceFrameBase*, Position*, ReferenceFrameBase*>())
		.def("to", &Velocity::to)
		.def("to_coordinate_system", &Velocity::to_coordinate_system)
		.def("to_global", &Velocity::to_global)
		;
	
	
	class_<ZeroReferenceFrame, bases<ReferenceFrameBase> >("ZeroReferenceFrame", init<>())
		;
	class_<ReferenceFrame, bases<ReferenceFrameBase> >("ReferenceFrame", init<Position*, ReferenceFrameBase*>())
		;
	class_<RotatedReferenceFrame, bases<ReferenceFrameBase> >("RotatedReferenceFrame", init<double, ReferenceFrameBase*>())
		;
	class_<MovingReferenceFrame, bases<ReferenceFrameBase> >("MovingReferenceFrame", init<Velocity*, ReferenceFrameBase*>())
		;
	class_<EqReferenceFrame, bases<ReferenceFrameBase> >("EqReferenceFrame", init<double, double, double, ReferenceFrameBase*>())
		;
//class_<TranslatedCoordinateSystem, bases<CoordinateSystem> >("TranslatedCoordinateSystem", init<double, double, double, CoordinateSystem*>())
	//	;
	//class_<Cylindrical , bases<CoordinateSystem> >("Cylindrical", init<Position*>())
	//	;
	//class_<RotatedCoordinateSystem, bases<CoordinateSystem> >("RotatedCoordinateSystem", init<CoordinateSystem*>())
	//;

}

}
