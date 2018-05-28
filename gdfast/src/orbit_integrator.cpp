#include "orbit_integrator.hpp"
#include "profile.hpp"

#include <boost/python.hpp>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
# include <boost/python/return_opaque_pointer.hpp>
# include <boost/python/def.hpp>
# include <boost/python/module.hpp>
# include <boost/python/return_value_policy.hpp>

//struct void
BOOST_PYTHON_OPAQUE_SPECIALIZED_TYPE_ID(gsl_odeiv_step_type);

namespace gd {
using namespace boost::python;
/*
struct gsl_odeiv_step_type_to_pyint
{
	gsl_odeiv_step_type_to_pyint() {
		boost::python::converter::registry::push_back(
			&convertible,
			&construct,
			boost::python::type_id<gsl_odeiv_step_type_to_pyint>()
		);
	}
	static PyObject* convert(const gsl_odeiv_step_type* x)
	{        
		return PyInt_FromLong((long int)x);
	}
	static const gsl_odeiv_step_type* execute(PyObject* x)
	{
		return (const gsl_odeiv_step_type*)PyInt_AsLong(x);
	}
};
*/



void py_export_orbit_integrator() {
	class_< OrbitIntegrator, boost::noncopyable >("OrbitIntegrator", no_init);
	class_< OrbitIntegratorSet, boost::noncopyable >("OrbitIntegratorSet", no_init);

	class_<OrbitIntegratorEuler, bases<OrbitIntegrator> >("OrbitIntegratorEuler", init<double, double, double, double>())
		.def("integrate", (&OrbitIntegratorEuler::integrate))
	;
	class_<OrbitIntegratorLeapFrog, bases<OrbitIntegrator> >("OrbitIntegratorLeapFrog", init<double, double, double, double>())
		.def("integrate", (&OrbitIntegratorLeapFrog::integrate))
	;
	class_<OrbitIntegratorSetLeapFrog, bases<OrbitIntegratorSet> >("OrbitIntegratorSetLeapFrog", init<ProfileModel*, int>())
			.def("integrate", (&OrbitIntegratorSetLeapFrog::integrate))
					;
	class_<OrbitIntegratorSetLeapFrog3d, bases<OrbitIntegratorSet> >("OrbitIntegratorSetLeapFrog3d", init<ProfileModel3d*, int>())
			.def("integrate", (&OrbitIntegratorSetLeapFrog3d::integrate))
					;
	class_<OrbitIntegratorGSL>("OrbitIntegratorGSL", init<ProfileModel*, optional<double, double, const gsl_odeiv_step_type*> >())
		.def("integrate", (&OrbitIntegratorGSL::integrate))
		;
	class_<OrbitIntegratorAxiGSL>("OrbitIntegratorAxiGSL", init<ProfileModelAxi*, optional<double, double, const gsl_odeiv_step_type*> >())
		.def("integrate", (&OrbitIntegratorAxiGSL::integrate))
		;
	class_<OrbitIntegrator2dGSL>("OrbitIntegrator2dGSL", init<ProfileModel2d*, optional<double, double, const gsl_odeiv_step_type*> >())
		.def("integrate", (&OrbitIntegrator2dGSL::integrate))
		;
	class_<OrbitIntegratorGSL3d>("OrbitIntegratorGSL3d", init<ProfileModel*, optional<double, double, const gsl_odeiv_step_type*> >())
		.def("integrate", (&OrbitIntegratorGSL3d::integrate))
		;
	opaque<gsl_odeiv_step_type>();
	scope().attr("gsl_odeiv_step_rk2") = (gsl_odeiv_step_type*)&gsl_odeiv_step_rk2;
	scope().attr("gsl_odeiv_step_rk4") = (gsl_odeiv_step_type*)&gsl_odeiv_step_rk4;
	scope().attr("gsl_odeiv_step_rkf45") = (gsl_odeiv_step_type*)&gsl_odeiv_step_rkf45;
	scope().attr("gsl_odeiv_step_rkck") = (gsl_odeiv_step_type*)&gsl_odeiv_step_rkck;
	scope().attr("gsl_odeiv_step_rk8pd") = (gsl_odeiv_step_type*)&gsl_odeiv_step_rk8pd;
	scope().attr("gsl_odeiv_step_rk2imp") = (gsl_odeiv_step_type*)&gsl_odeiv_step_rk2imp;
	scope().attr("gsl_odeiv_step_rk4imp") = (gsl_odeiv_step_type*)&gsl_odeiv_step_rk4imp;
	//to_python_converter<const gsl_odeiv_step_type*, gsl_odeiv_step_type_to_pyint> ();
	//lvalue_from_pytype<pyfile_to_FILE,&PyFile_Type> ();
	//lvalue_from_pytype<gsl_odeiv_step_type_to_pyint, &PyInt_Type> ();
	//opaque

	//handle<> h(gsl_odeiv_step_rkck);
	//scope().attr("gsl_odeiv_step_rkck") = h;
}

int f(double t, const double y[], double f[], void *params) {
	//OrbitIntegratorGSL *oi = ( OrbitIntegratorGSL*)params;
	ProfileModel* profile_model = (ProfileModel*)params;
	double r = sqrt(y[0]*y[0] + y[1]*y[1]);
	double dphidr = profile_model->dphidr(r)* 3.2407764868054621e-17;
	//double dphidr = G * mass * r / pow((r*r + scale*scale), (3./2)) * 3.2407764868054621e-17; 
	double Fx = -dphidr*y[0]/r;
	double Fy = -dphidr*y[1]/r;
	//printf ("%.5e %.5e %.5e %.5e\n", r, dphidr, Fx, Fy);
	f[0] = y[2]*3.2407764868054621e-17;
	f[1] = y[3]*3.2407764868054621e-17;
	f[2] = Fx;
	f[3] = Fy;
	return GSL_SUCCESS;
}

struct f_axi_params_type {
	ProfileModelAxi* profile_model;
	double Lz;
};
int f_axi(double t, const double y[], double f[], void *params) {
	//OrbitIntegratorGSL *oi = ( OrbitIntegratorGSL*)params;
	ProfileModelAxi* profile_model = ((f_axi_params_type*)params)->profile_model;
	double Lz = ((f_axi_params_type*)params)->Lz;
	double x = y[0];
	double z = y[1];
	double dphidx_eff = profile_model->dphidx_eff(x, z, Lz)* 3.2407764868054621e-17;
	double dphidz_eff = profile_model->dphidz_eff(x, z, Lz)* 3.2407764868054621e-17;
//double dphidr = G * mass * r / pow((r*r + scale*scale), (3./2)) * 3.2407764868054621e-17; 
	double Fx = -dphidx_eff;
	double Fz = -dphidz_eff;
	//printf ("%.5e %.5e %.5e %.5e\n", r, dphidr, Fx, Fy);
	f[0] = y[2]*3.2407764868054621e-17;
	f[1] = y[3]*3.2407764868054621e-17;
	f[2] = Fx;
	f[3] = Fz;
	return GSL_SUCCESS;
}
struct orbit2d_params_type {
	ProfileModel2d* profile_model;
	double angular_velocity;
};

int f2d(double t, const double y[], double f[], void *params) {
	//OrbitIntegratorGSL3d *oi = ( OrbitIntegratorGSL3d*)params;
	orbit2d_params_type* orbit2d_params = (orbit2d_params_type*)params;
	//ProfileModel2d* profile_model = (ProfileModel2d*)params;
	//double r = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
	double angular_velocity = orbit2d_params->angular_velocity * 3.2407764868054621e-17;
	double angular_velocity_sq = angular_velocity * angular_velocity;
	//printf("%f ", angular_velocity);
	double dphidx = orbit2d_params->profile_model->dphidx(y[0], y[1])* 3.2407764868054621e-17 - angular_velocity_sq*y[0];
	double dphidy = orbit2d_params->profile_model->dphidy(y[0], y[1])* 3.2407764868054621e-17 - angular_velocity_sq*y[1];
//double dphidr = G * mass * r / pow((r*r + scale*scale), (3./2)) * 3.2407764868054621e-17; 
	double Fx = -dphidx + 2 * angular_velocity * y[3];
	double Fy = -dphidy - 2 * angular_velocity * y[2];
	//printf ("%.5e %.5e %.5e %.5e\n", r, dphidr, Fx, Fy);
	f[0] = y[2]*3.2407764868054621e-17;
	f[1] = y[3]*3.2407764868054621e-17;
	f[2] = Fx;
	f[3] = Fy;
	return GSL_SUCCESS;
}

int f3d(double t, const double y[], double f[], void *params) {
	//OrbitIntegratorGSL3d *oi = ( OrbitIntegratorGSL3d*)params;
	ProfileModel* profile_model = (ProfileModel*)params;
	double r = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
	double dphidr = profile_model->dphidr(r)* 3.2407764868054621e-17;
	//double dphidr = G * mass * r / pow((r*r + scale*scale), (3./2)) * 3.2407764868054621e-17; 
	double Fx = -dphidr*y[0]/r;
	double Fy = -dphidr*y[1]/r;
	double Fz = -dphidr*y[2]/r;
	//printf ("%.5e %.5e %.5e %.5e\n", r, dphidr, Fx, Fy);
	f[0] = y[3]*3.2407764868054621e-17;
	f[1] = y[4]*3.2407764868054621e-17;
	f[2] = y[5]*3.2407764868054621e-17;
	f[3] = Fx;
	f[4] = Fy;
	f[5] = Fz;
	return GSL_SUCCESS;
}
/*
int  (hEadjust) (void * state, size_t dim, unsigned int ord, const double y[], const double yerr[], const double yp[], double * h) {
	return GSL_ODEIV_HADJ_NIL;
}

const gsl_odeiv_control_type Econtrol_type = {"E control", NULL, NULL, &hEadjust, NULL};
		
struct Econtrol_state {
	Profile* profile;
	double maxrelE, time;
};

gsl_odeiv_control* gsl_odeiv_control_E(Profile* profile, double maxrelE, double time) {
	gsl_odeiv_control* c = (gsl_odeiv_control*)malloc(sizeof(gsl_odeiv_control));
	//gsl_odeiv_control * c = gsl_odeiv_control_alloc (&Econtrol_type);
	
	Econtrol_state* econtrol_state = (Econtrol_state*)malloc(sizeof(Econtrol_state));
	econtrol_state->profile = profile;
	econtrol_state->maxrelE = maxrelE;
	econtrol_state->time = time;
	
	c->state = econtrol_state;
	c->type = &Econtrol_type;
	return c;
	 
		
}
*/
void OrbitIntegrator2dGSL::integrate(int N, double angular_velocity, double_vector dt, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix yout , double_matrix vxout, double_matrix vyout) {
	double* xoutp = xout.data().begin();
	double* youtp = yout.data().begin();
	double* vxoutp = vxout.data().begin();
	double* vyoutp = vyout.data().begin();
	int no_particles = dt.size();
	
	
	for(int i = 0; i < no_particles; i++) {
		//double t = 0;
		//double Fx, Fy;
		double  x = q0(i,0);
		double  y = q0(i,1);
		double vx = v0(i,0);
		double vy = v0(i,1);
		//const gsl_odeiv_step_type * T = gsl_odeiv_step_rkck;
		double t1 = 0.0, t2 = dt(i);
		double stepsize = dt(i);
		double phasespace[4];
		//printf("%p %p\n", gsl_odeiv_step_type_, gsl_odeiv_step_rkck);
		gsl_odeiv_step * s		= gsl_odeiv_step_alloc (gsl_odeiv_step_type_, 4);
		gsl_odeiv_control * c	= gsl_odeiv_control_y_new (error_abs, error_rel);
		gsl_odeiv_evolve * e	= gsl_odeiv_evolve_alloc (4);
		
		orbit2d_params_type params = {profile_model, angular_velocity};
		gsl_odeiv_system sys = {f2d, NULL, 4, &params};
		phasespace[0] = x;
		phasespace[1] = y;
		phasespace[2] = vx;
		phasespace[3] = vy;
		
		for(int j = 0; j < N; j++) {
			
			int status = GSL_SUCCESS;
			while(t1 < t2) {
				status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t1, t2, &stepsize, phasespace);
				if (status != GSL_SUCCESS)
					break;
				//printf ("%.5e %.5e %.5e %.5e\n", t, phasespace[0], phasespace[1], stepsize);
			}
			t2 += dt(i);
			if (status != GSL_SUCCESS) {
				printf ("error integrating\n");
			}
			xoutp[i*N+j] = phasespace[0];
			youtp[i*N+j] = phasespace[1];
			vxoutp[i*N+j] = phasespace[2];
			vyoutp[i*N+j] = phasespace[3];
		}
	}
}

void OrbitIntegratorGSL3d::integrate(int N, double_vector dt, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix yout, double_matrix zout, double_matrix vxout, double_matrix vyout, double_matrix vzout) {
	double* xoutp = xout.data().begin();
	double* youtp = yout.data().begin();
	double* zoutp = zout.data().begin();
	double* vxoutp = vxout.data().begin();
	double* vyoutp = vyout.data().begin();
	double* vzoutp = vzout.data().begin();
	int no_particles = dt.size();

	for(int i = 0; i < no_particles; i++) {
		//double t = 0;
		//double Fx, Fy;
		double  x = q0(i,0);
		double  y = q0(i,1);
		double  z = q0(i,2);
		double vx = v0(i,0);
		double vy = v0(i,1);
		double vz = v0(i,2);
		//const gsl_odeiv_step_type * T = gsl_odeiv_step_rkck;
		double t1 = 0.0, t2 = dt(i);
		double stepsize = dt(i);
		double phasespace[6];
		//printf("%p %p\n", gsl_odeiv_step_type_, gsl_odeiv_step_rkck);
		gsl_odeiv_step * s		= gsl_odeiv_step_alloc (gsl_odeiv_step_type_, 6);
		gsl_odeiv_control * c	= gsl_odeiv_control_y_new (error_abs, error_rel);
		gsl_odeiv_evolve * e	= gsl_odeiv_evolve_alloc (6);
		
		gsl_odeiv_system sys = {f3d, NULL, 6, profile_model};
		phasespace[0] = x;
		phasespace[1] = y;
		phasespace[2] = z;
		phasespace[3] = vx;
		phasespace[4] = vy;
		phasespace[5] = vz;

		for(int j = 0; j < N; j++) {

			int status = GSL_SUCCESS;
			while(t1 < t2) {
				status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t1, t2, &stepsize, phasespace);
				if (status != GSL_SUCCESS)
					break;
				//printf ("%.5e %.5e %.5e %.5e\n", t, phasespace[0], phasespace[1], stepsize);
			}
			t2 += dt(i);
			if (status != GSL_SUCCESS) {
				printf ("error integrating\n");
			}
			xoutp[i*N+j] = phasespace[0];
			youtp[i*N+j] = phasespace[1];
			zoutp[i*N+j] = phasespace[2];
			vxoutp[i*N+j] = phasespace[3];
			vyoutp[i*N+j] = phasespace[4];
			vzoutp[i*N+j] = phasespace[5];
		}
	}
}

void OrbitIntegratorGSL::integrate(int N, double_vector dt, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix yout, double_matrix vxout, double_matrix vyout) {
	double* xoutp = xout.data().begin();
	double* youtp = yout.data().begin();
	double* vxoutp = vxout.data().begin();
	double* vyoutp = vyout.data().begin();
	int no_particles = dt.size();

	for(int i = 0; i < no_particles; i++) {
		//double t = 0;
		//double Fx, Fy;
		double  x = q0(i,0);
		double  y = q0(i,1);
		double vx = v0(i,0);
		double vy = v0(i,1);
		//const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
		double t1 = 0.0, t2 = dt(i);
		double stepsize = dt(i);
		double phasespace[4];

		gsl_odeiv_step * s		= gsl_odeiv_step_alloc (gsl_odeiv_step_type_, 4);
		gsl_odeiv_control * c	= gsl_odeiv_control_y_new (error_abs, error_rel);
		gsl_odeiv_evolve * e	= gsl_odeiv_evolve_alloc (4);
		
		gsl_odeiv_system sys = {f, NULL, 4, profile_model};
		phasespace[0] = x;
		phasespace[1] = y;
		phasespace[2] = vx;
		phasespace[3] = vy;

		for(int j = 0; j < N; j++) {

			int status = GSL_SUCCESS;
			while(t1 < t2) {
				status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t1, t2, &stepsize, phasespace);
				if (status != GSL_SUCCESS)
					break;
				//printf ("%.5e %.5e %.5e %.5e\n", t, phasespace[0], phasespace[1], stepsize);
			}
			t2 += dt(i);
			if (status != GSL_SUCCESS) {
				printf ("error integrating\n");
			}
			xoutp[i*N+j] = phasespace[0];
			youtp[i*N+j] = phasespace[1];
			vxoutp[i*N+j] = phasespace[2];
			vyoutp[i*N+j] = phasespace[3];
		}
	}
}
	

void OrbitIntegratorAxiGSL::integrate(int N, double_vector dt, double_vector Lz, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix zout, double_matrix vxout, double_matrix vzout) {
	double* xoutp = xout.data().begin();
	double* zoutp = zout.data().begin();
	double* vxoutp = vxout.data().begin();
	double* vzoutp = vzout.data().begin();
	int no_particles = dt.size();
	
	for(int i = 0; i < no_particles; i++) {
		//double t = 0;
		//double Fx, Fy;
		double  x = q0(i,0);
		double  z = q0(i,1);
		double vx = v0(i,0);
		double vz = v0(i,1);
		double Lzi = Lz(i);
//const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
		double t1 = 0.0, t2 = dt(i);
		double stepsize = dt(i);
		double phasespace[4];
		
		gsl_odeiv_step * s		= gsl_odeiv_step_alloc (gsl_odeiv_step_type_, 4);
		gsl_odeiv_control * c	= gsl_odeiv_control_y_new (error_abs, error_rel);
		gsl_odeiv_evolve * e	= gsl_odeiv_evolve_alloc (4);
		
		f_axi_params_type params = {profile_model, Lzi};
		gsl_odeiv_system sys = {f_axi, NULL, 4, &params};
		phasespace[0] = x;
		phasespace[1] = z;
		phasespace[2] = vx;
		phasespace[3] = vz;
		
		for(int j = 0; j < N; j++) {
			
			int status = GSL_SUCCESS+1;
			do {
				status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t1, t2, &stepsize, phasespace);
			/*if (status != GSL_SUCCESS) {
				printf ("error integrating\n");
				break;*/
			} while(status != GSL_SUCCESS); 
				//printf ("%.5e %.5e %.5e %.5e\n", t, phasespace[0], phasespace[1], stepsize);
			//}
			t2 += dt(i);
			if (status != GSL_SUCCESS) {
				printf ("error integrating\n");
				exit(-1);
			}
			xoutp[i*N+j] = phasespace[0];
			zoutp[i*N+j] = phasespace[1];
			vxoutp[i*N+j] = phasespace[2];
			vzoutp[i*N+j] = phasespace[3];
		}
	}
}

}
	
	
