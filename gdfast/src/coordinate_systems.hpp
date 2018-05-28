
#include <Eigen/Dense>
#include <Eigen/Geometry>
//#include <s

namespace gd {
using namespace Eigen;
//using namespace std;
typedef Matrix<double, 3, 1> Vector3;
typedef Matrix<double, 3, 3> Matrix3;


/*

 System3(System2(System1(..)))

*/
class CoordinateSystem {
public:
	CoordinateSystem() {
	}/*
	virtual Matrix3 matrix(double x, double y, double z) {
		return this->matrix_local(x, y, z) * parent->matrix(x, y, z);
	}*/
	virtual Vector3 ex(double x, double y, double z) = 0;
	virtual Vector3 ey(double x, double y, double z) = 0;
	virtual Vector3 ez(double x, double y, double z) = 0;
	//virtual Vector3 get_origin_cartesian() { return Vector3(0,0,0); }
	virtual Vector3 to_cartesian(Vector3d v) = 0;
	virtual Vector3 from_cartesian(Vector3d v) = 0;

	Matrix3 basis_matrix(double x, double y, double z) {
		Matrix3 m;
		m.block<3,1>(0,0) = ex(x, y, z);
		m.block<3,1>(0,1) = ey(x, y, z);
		m.block<3,1>(0,2) = ez(x, y, z);
		return m;
	}
	/*Matrix3 matrix_local(double x, double y, double z) {
		Matrix3 m;
		m.block<3,1>(0,0) = ex(x, y, z);
		m.block<3,1>(0,1) = ey(x, y, z);
		m.block<3,1>(0,2) = ez(x, y, z);
		return m;
	}
	Matrix3 matrix(Vector3 v) {
		return this->matrix(v(0),v(1),v(2));
	}

	Matrix3 transformation_matrix(double x, double y, double z, CoordinateSystem* prime_system) {
		Matrix3 m1 = this->matrix(x, y, z);
		Matrix3 m2 = prime_system->matrix(x, y, z);
		return m1*m2.transpose();
	}*/
};

class Cartesian : public CoordinateSystem {
public:
	Vector3 ex(double x, double y, double z) {
		return Vector3(1, 0, 0);
	}
	Vector3 ey(double x, double y, double z) {
		return Vector3(0, 1, 0);
	}
	Vector3 ez(double x, double y, double z) {
		return Vector3(0, 0, 1);
	}
	Vector3 to_cartesian(Vector3d v) {
		return v;
	}
	Vector3 from_cartesian(Vector3d v) {
		return v;
	}
};

class Cylindrical : public CoordinateSystem {
public:
	Vector3 ex(double x, double y, double z) {
		double rho = sqrt(x*x+y*y);
		return Vector3(x/rho, y/rho, 0);
	}
	Vector3 ey(double x, double y, double z) {
		double rho = sqrt(x*x+y*y);
		return Vector3(y/rho, -x/rho, 0);
	}
	Vector3 ez(double x, double y, double z) {
		return Vector3(0, 0, 1);
	}
	Vector3 to_cartesian(Vector3d v) {
		double rho = v(0);
		double theta = v(1);
		double z = v(2);
		double x = rho * cos(theta);
		double y = rho * sin(theta);
		Vector3 cartesian(x,y,z);
		return cartesian;
	}
	Vector3 from_cartesian(Vector3d v) {
		double x = v(0);
		double y = v(1);
		double z = v(2);
		double rho = sqrt(x*x+y*y);
		double theta = atan2(y,x);
		Vector3 cylindrical(rho,theta,z);
		return cylindrical;
	}
};


class SphericalGalactic : public CoordinateSystem {
public:
	SphericalGalactic() {
		//rotation = AngleAxisd(M_PI, Vector3::UnitZ());
		//std::cout << "rot[" << rotation << "]" << std::endl; 
	}
	Vector3 ex(double x, double y, double z) {
		double r = sqrt(x*x+y*y+z*z);
		return Vector3(x/r, y/r, z/r);
	}
	Vector3 ey(double x, double y, double z) {
		double rho = sqrt(x*x+y*y);
		return Vector3(-y/rho, x/rho, 0);
	}
	Vector3 ez(double x, double y, double z) {
		double r = sqrt(x*x+y*y+z*z);
		double rho = sqrt(x*x+y*y);
		Vector3 v = Vector3(z*x/rho/r, z*y/rho/r, -rho/r);
		return v;
	}
	Vector3 to_cartesian(Vector3d v) {
		double r = v(0);
		double phi = v(1);
		double theta = v(2);
		double x = r * sin(theta) * cos(phi);
		double y = r * sin(theta) * sin(phi);
		double z = r * cos(theta); 
		Vector3 cartesian(x,y,z);
		return cartesian;
	}
	Vector3 from_cartesian(Vector3d cartesian) {
		double x = cartesian(0);
		double y = cartesian(1);
		double z = cartesian(2);
		double r = sqrt(x*x+y*y+z*z);
		double phi = atan2(y,x); 
		double theta = acos(z/r);
		Vector3 spherical(r,phi,theta);
		return spherical;
	}
};



class Coordinate {
public:
	CoordinateSystem* coordinate_system;
	Vector3 x;
	Coordinate(double x1, double x2, double x3, CoordinateSystem* coordinate_system) : coordinate_system(coordinate_system), x(x1, x2, x3) {
	}
	/*virtual Matrix3 matrix(double x, double y, double z) {
		return this->matrix_local(x, y, z) * parent->matrix(x, y, z);
	}
	virtual Vector3 ex(double x, double y, double z) = 0;
	virtual Vector3 ey(double x, double y, double z) = 0;
	virtual Vector3 ez(double x, double y, double z) = 0;
	//virtual Vector3 get_origin_cartesian() { return Vector3(0,0,0); }
	virtual Vector3 cartesian_transformation(Vector3 v) { return v;}*/
	virtual Vector3 to_cartesian() {
		return coordinate_system->to_cartesian(x);
	}
	//virtual Vector3 to_cartesian_local(Vector3d v) = 0;
	//virtual Vector3 from_cartesian_local(Vector3d v) = 0;
	/*Matrix3 matrix_local(double x, double y, double z) {
		Matrix3 m;
		m.block<3,1>(0,0) = ex(x, y, z);
		m.block<3,1>(0,1) = ey(x, y, z);
		m.block<3,1>(0,2) = ez(x, y, z);
		return m;
	}
	Matrix3 matrix(Vector3 v) {
		return this->matrix(v(0),v(1),v(2));
	}
	
	Matrix3 transformation_matrix(double x, double y, double z, CoordinateSystem* prime_system) {
		Matrix3 m1 = this->matrix(x, y, z);
		Matrix3 m2 = prime_system->matrix(x, y, z);
		return m1*m2.transpose();
	}*/
};

class VelocityCoordinate {
public:
	CoordinateSystem* coordinate_system;
	Vector3 v;
	VelocityCoordinate(double v1, double v2, double v3, CoordinateSystem* coordinate_system) : coordinate_system(coordinate_system), v(v1, v2, v3) {
	}
	/*virtual Matrix3 matrix(double x, double y, double z) {
		return this->matrix_local(x, y, z) * parent->matrix(x, y, z);
	}
	virtual Vector3 ex(double x, double y, double z) = 0;
	virtual Vector3 ey(double x, double y, double z) = 0;
	virtual Vector3 ez(double x, double y, double z) = 0;
	//virtual Vector3 get_origin_cartesian() { return Vector3(0,0,0); }
	virtual Vector3 cartesian_transformation(Vector3 v) { return v;}*/
	virtual Vector3 to_cartesian(double x, double y, double z) {
		return \
			v(0) * coordinate_system->ex(x, y, z) +\
			v(1) * coordinate_system->ey(x, y, z) +\
			v(2) * coordinate_system->ez(x, y, z);
	}
	virtual Vector3 to_cartesian_vec(Vector3 p) {
		return to_cartesian(p(0), p(1), p(2));
	}

};

class ReferenceFrameBase {
public:
	Matrix3 basis_matrix(double x, double y, double z) {
		Matrix3 m;
		m.block<3,1>(0,0) = ex(x, y, z);
		m.block<3,1>(0,1) = ey(x, y, z);
		m.block<3,1>(0,2) = ez(x, y, z);
		return m;
	}
	virtual Vector3 ex(double x, double y, double z) = 0;
	virtual Vector3 ey(double x, double y, double z) = 0;
	virtual Vector3 ez(double x, double y, double z) = 0;
	virtual Vector3 to_global(Vector3 local) = 0;
	virtual Vector3 to_local(Vector3 global) = 0;
	virtual Vector3 to_global_velocity(Vector3 v_local) { return v_local; }
	virtual Vector3 to_local_velocity(Vector3 v_global) { return v_global; }
};


class Position {
public:
	Coordinate* c;
	ReferenceFrameBase* f;
	Position(Coordinate* c, ReferenceFrameBase* f) : c(c), f(f) {
	}
	Vector3 to(ReferenceFrameBase* target_frame) {
		return target_frame->to_local( f->to_global(c->to_cartesian()) );;
	}
	Vector3 to_global() {
		return f->to_global(c->to_cartesian());
	}
	Vector3 to_coordinate_system(ReferenceFrameBase* target_frame, CoordinateSystem* target_coordinate_system) {
		Vector3 v = to(target_frame);
		return target_coordinate_system->from_cartesian(v);
	}
};

class Velocity {
public:
	VelocityCoordinate* vc;
	ReferenceFrameBase* f;
	Position* p;
	ReferenceFrameBase* pos_frame;
	Velocity(VelocityCoordinate* vc, ReferenceFrameBase* f, Position* p, ReferenceFrameBase* pos_frame) : vc(vc), f(f), p(p), pos_frame(pos_frame) {
	}
	Vector3 to_coordinate_system(ReferenceFrameBase* target_frame, CoordinateSystem* target_coordinate_system) {
		Vector3 v = to(target_frame);
		Vector3 local_coordinate = p->to(target_frame); //f->to_global(Vector3(0, 0, 0));
		double x = local_coordinate(0);
		double y = local_coordinate(1);
		double z = local_coordinate(2);
		//std::cout << "[x,y,z = ( " << x << ", " << y << ", " << z <<")]" << std::endl;
		Matrix3 m = target_coordinate_system->basis_matrix(x, y, z);
		return m.inverse() * v; 
	}
	Vector3 to(ReferenceFrameBase* target_frame) {
		Vector3 global_velocity = to_global();
		return target_frame->to_local_velocity(global_velocity);
	}
	Vector3 to_global() {
		Vector3 local_coordinate = p->to(pos_frame); //f->to_global(Vector3(0, 0, 0));
		Vector3 local_velocity = vc->to_cartesian_vec(local_coordinate);
		Vector3 global_velocity = f->to_global_velocity(local_velocity);
		return global_velocity;
	}	
};

class ZeroReferenceFrame : public ReferenceFrameBase {
public:
	ZeroReferenceFrame()  {
	}
	Vector3 ex(double x, double y, double z) {
		return Vector3(1, 0, 0);
	}
	Vector3 ey(double x, double y, double z) {
		return Vector3(0, 1, 0);
	}
	Vector3 ez(double x, double y, double z) {
		return Vector3(0, 0, 1);
	}
	virtual Vector3 to_global(Vector3 local) { return local; }
	virtual Vector3 to_local(Vector3 global) { return global; }
};

class ReferenceFrame : public ReferenceFrameBase {
public:
	Vector3 x0;
	ReferenceFrameBase* parent;
	ReferenceFrame(Position* origin, ReferenceFrameBase* parent) : x0(origin->to(parent)), parent(parent) {
	}
	virtual Vector3 to_global(Vector3 local) { return parent->to_global(local + x0); }
	virtual Vector3 to_local(Vector3 global) { return parent->to_local(global) - x0; }
	Vector3 ex(double x, double y, double z) {
		return Vector3(1, 0, 0);
	}
	Vector3 ey(double x, double y, double z) {
		return Vector3(0, 1, 0);
	}
	Vector3 ez(double x, double y, double z) {
		return Vector3(0, 0, 1);
	}
};

class RotatedReferenceFrame : public ReferenceFrameBase {
public:
	Matrix3 rotation;
	Matrix3 rotation_inverse;
	ReferenceFrameBase *parent;
	RotatedReferenceFrame(double angle, ReferenceFrameBase* parent) : parent(parent) {
		rotation = AngleAxisd(angle, Vector3::UnitZ());
		rotation_inverse = rotation.inverse();
	}
	virtual Vector3 to_global_velocity(Vector3 v_local) { return parent->to_global_velocity(rotation*v_local); }
	virtual Vector3 to_local_velocity(Vector3 v_global) { return rotation_inverse*parent->to_local_velocity(v_global); }
	virtual Vector3 to_global(Vector3 local) { return parent->to_global(rotation*local); }
	virtual Vector3 to_local(Vector3 global) { return rotation_inverse*parent->to_local(global); }
	//Vector3 get_origin_cartesian() { return origin->to_cartesian(); }
	Vector3 ex(double x, double y, double z) {
		return rotation * Vector3::UnitX();
	}
	Vector3 ey(double x, double y, double z) {
		return rotation * Vector3::UnitY();
	}
	Vector3 ez(double x, double y, double z) {
		return rotation * Vector3::UnitZ();
	}
};

class EqReferenceFrame : public ReferenceFrameBase {
public:
	Matrix3 rotation;
	Matrix3 rotation_inverse;
	ReferenceFrameBase *parent;
	EqReferenceFrame(double theta0, double a_NGP, double d_NGP, ReferenceFrameBase* parent) : parent(parent) {
		Matrix3 a,b,c;
		c << 	cos(a_NGP),  sin(a_NGP), 0,
				sin(a_NGP), -cos(a_NGP), 0,
				0, 0, 1;
		b << 	-sin(d_NGP), 0, cos(d_NGP),
				0, -1, 0,
				cos(d_NGP), 0, sin(d_NGP);
		a << 	cos(theta0),  sin(theta0), 0,
				sin(theta0), -cos(theta0), 0,
				0, 0, 1;
		rotation = a*b*c;
		rotation_inverse = rotation.inverse();
		//std::cout << "[" << rotation << "]" << std::endl;
		//std::cout << rotation_inverse << std::endl;
	}
	virtual Vector3 to_global_velocity(Vector3 v_local) { return parent->to_global_velocity(rotation*v_local); }
	virtual Vector3 to_local_velocity(Vector3 v_global) { return rotation_inverse*parent->to_local_velocity(v_global); }
	virtual Vector3 to_global(Vector3 local) { return parent->to_global(rotation*local); }
	virtual Vector3 to_local(Vector3 global) { return rotation_inverse*parent->to_local(global); }
	//Vector3 get_origin_cartesian() { return origin->to_cartesian(); }
	Vector3 ex(double x, double y, double z) {
		return rotation * Vector3::UnitX();
	}
	Vector3 ey(double x, double y, double z) {
		return rotation * Vector3::UnitY();
	}
	Vector3 ez(double x, double y, double z) {
		return rotation * Vector3::UnitZ();
	}
};

class MovingReferenceFrame : public ReferenceFrameBase {
public:
	ReferenceFrameBase* parent;
	Vector3 velocity_frame;
	MovingReferenceFrame(Velocity* velocity , ReferenceFrameBase* parent) : parent(parent) {
		velocity_frame = velocity->to(parent);
	}
	virtual Vector3 to_global_velocity(Vector3 v_local) { return parent->to_global_velocity(v_local+velocity_frame); }
	virtual Vector3 to_local_velocity(Vector3 v_global) { return parent->to_local_velocity(v_global)-velocity_frame; }
	virtual Vector3 to_global(Vector3 local) { return parent->to_global(local); }
	virtual Vector3 to_local(Vector3 global) { return parent->to_local(global); }
	Vector3 ex(double x, double y, double z) {
		return Vector3::UnitX();
	}
	Vector3 ey(double x, double y, double z) {
		return Vector3::UnitY();
	}
	Vector3 ez(double x, double y, double z) {
		return Vector3::UnitZ();
	}
};



/*
class TranslatedReferenceFrame : public CoordinateSystem {
public:
	Vector3 translation;
	TranslatedReferenceFrame(double x, double y, double z, CoordinateSystem* parent) : CoordinateSystem(parent), translation(x, y, z) {
	}
	Vector3 ex(double x, double y, double z) {
		return Vector3(1, 0, 0);
	}
	Vector3 ey(double x, double y, double z) {
		return Vector3(0, 1, 0);
	}
	Vector3 ez(double x, double y, double z) {
		return Vector3(0, 0, 1);
	}
	Vector3 to_cartesian_local(Vector3d v) {
		return v+translation;
	}
	Vector3 from_cartesian_local(Vector3d cartesian) {
		return cartesian-translation;
	}
};*/
/*
class Cylindrical : public CoordinateSystem {
public:
	Position* origin;
	Cylindrical(Position* origin) : origin(origin) {
	}
	Vector3 get_origin_cartesian() { return origin->to_cartesian(); }
	Vector3 ex(double x, double y, double z) {
		double rho = sqrt(x*x+y*y);
		return Vector3(x/rho, y/rho, 0);
	}
	Vector3 ey(double x, double y, double z) {
		double rho = sqrt(x*x+y*y);
		return Vector3(y/rho, -x/rho, 0);
	}
	Vector3 ez(double x, double y, double z) {
		return Vector3(0, 0, 1);
	}
	Vector3 to_cartesian_local(Vector3d v) {
		double rho = v(0);
		double theta = v(1);
		double z = v(2);
		double x = rho * cos(theta);
		double y = rho * sin(theta);
		Vector3 cartesian(x,y,z);
		return cartesian;
	}
	Vector3 from_cartesian_local(Vector3d v) {
		double x = v(0);
		double y = v(1);
		double z = v(2);
		double rho = sqrt(x*x+y*y);
		double theta = atan2(y,x);
		Vector3 cylindrical(rho,theta,z);
		return cylindrical;
	}
};
*/
/*
class RotatedCoordinateSystem : public CoordinateSystem {
public:
	Matrix3 rotation;
	Matrix3 rotation_inverse;
	RotatedCoordinateSystem(CoordinateSystem* parent) : CoordinateSystem(parent) {
		rotation = AngleAxisd(M_PI, Vector3::UnitZ());
		rotation_inverse = rotation.inverse();
	}
	//Vector3 get_origin_cartesian() { return origin->to_cartesian(); }
	Vector3 ex(double x, double y, double z) {
		return rotation * parent->ex(x, y, z);
	}
	Vector3 ey(double x, double y, double z) {
		return rotation * parent->ey(x, y, z);
	}
	Vector3 ez(double x, double y, double z) {
		return rotation * parent->ez(x, y, z);
	}
	Vector3 to_cartesian_local(Vector3d v) {
		//return rotation*parent->to_cartesian(v);
		return rotation*v;
	}
	Vector3 from_cartesian_local(Vector3d cartesian) {
		//return parent->from_cartesian(rotation_inverse * cartesian);
		return rotation_inverse*cartesian;
	}
};

*/


/*
class Spherical {
	
	Vector3 ex(double x, double y, double z) {
		
	}
	Vector3 ey(double x, double y, double z) {
	}
	Vector3 ez(double x, double y, double z) {
	}
};*/

}