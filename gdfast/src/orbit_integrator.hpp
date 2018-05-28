#pragma once

#include <gsl/gsl_odeiv.h>
#include "datatypes.hpp"
#include <tuple>
#include "profile.hpp"
#include "profile3d.hpp"
#include "profile_axi.hpp"
#include "profile_2d.hpp"

namespace gd {
	
class ProfileModel;
class Profile;
using namespace std;

void py_export_orbit_integrator();
	
class OrbitIntegrator {
};
class OrbitIntegratorSet {
};

class OrbitIntegratorEuler : public OrbitIntegrator {
public:
	double x, y;
	double vx, vy;
	OrbitIntegratorEuler(double x0, double y0, double vx0, double vy0) : x(x0), y(y0), vx(vy0), vy(vy0) {
		
	}
	void integrate(Profile* profile, int N, double_vector x_out, double_vector y_out, double_vector vx_out, double_vector vy_out, double dt=0.1) {
		//double t = 0;
		//double i = 0;
		double Fx, Fy;
		
		for(int i = 0; i < N; i++) {
			tie(Fx, Fy) = profile->dphidxy(x, y);
			Fx = -Fx;
			Fy = -Fy;
			
			vx += Fx * dt;
			x += vx * dt;
			
			vy += Fy * dt;
			y += vy * dt;
			
			x_out[i] = x;
			y_out[i] = y;
			vx_out[i] = vx;
			vy_out[i] = vy;
			
		
		}
	}
	
};

class OrbitIntegratorLeapFrog : public OrbitIntegrator {
public:
	double x, y;
	double vx, vy;
	OrbitIntegratorLeapFrog(double x0, double y0, double vx0, double vy0) : x(x0), y(y0), vx(vy0), vy(vy0) {
	}
	void integrate(Profile* profile, int N, double_vector x_out, double_vector y_out, double_vector vx_out, double_vector vy_out, double dt=0.1) {
		//double t = 0;
		//double i = 0;
		double Fx, Fy;
		
		double hdt = dt / 2;
		tie(Fx, Fy) = profile->dphidxy(x, y);
		Fx = -Fx;
		Fy = -Fy;
			
		vx += Fx * hdt;
		vy += Fy * hdt;
		
		for(int i = 0; i < N; i++) {
			x += vx * dt;
			y += vy * dt;

			tie(Fx, Fy) = profile->dphidxy(x, y);
			Fx = -Fx;
			Fy = -Fy;
			
			vx += Fx * dt;
			vy += Fy * dt;
			
			x_out[i] = x;
			y_out[i] = y;
			vx_out[i] = vx;
			vy_out[i] = vy;
			
		
		}
	}
	
};

class OrbitIntegratorSetLeapFrog : public OrbitIntegratorSet {
public:
OrbitIntegratorSetLeapFrog(ProfileModel* profile_model, int M) : profile_model(profile_model), M(M) {}
	void integrate(int N, double_vector dt, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix yout, double_matrix vxout, double_matrix vyout) {
		//float_matrix xout(N, x0.size());
		//float_matrix yout(N, x0.size());
		//float_matrix vxout(N, x0.size());
// 		//float_matrix vyout(N, x0.size());
		double* xoutp = xout.data().begin();
		double* youtp = yout.data().begin();
		double* vxoutp = vxout.data().begin();
		double* vyoutp = vyout.data().begin();
		int no_particles = dt.size();
		//float G = 4.30200406461e-06;
		double r;
		//double dphidr; 
		for(int i = 0; i < no_particles; i++) {
			//double t = 0;
			double Fx, Fy;
			double  x = q0(i,0);
			double  y = q0(i,1);
			double vx = v0(i,0);
			double vy = v0(i,1);
			r = sqrt(x*x+y*y);
			/*double  x = x0(i);
			double  y = y0(i);
			double vx = vx0(i);
			double vy = vy0(i);*/
			double dti = dt(i)/M;
			
			double hdti = dti / 2;
			//tie(Fx, Fy) = profile->dphidxy(x,x y);
			double dphidr;

			for(int j = 0; j < N; j++) {
				// half a kick
				dphidr = profile_model->dphidr(r)* 3.2407764868054621e-17;
				r = sqrt(x*x+y*y);
				Fx = -dphidr*x/r;
				Fy = -dphidr*y/r;
				vx += Fx * hdti;
				vy += Fy * hdti;
				
				// in between steps can be full kicks and steps
				for(int k = 0; k < (M-1); k++) {
					// step
					x += vx * dti* 3.2407764868054621e-17;
					y += vy * dti* 3.2407764868054621e-17;
					// kick
					r = sqrt(x*x+y*y);
					dphidr = profile_model->dphidr(r)* 3.2407764868054621e-17;
					Fx = -dphidr*x/r;
					Fy = -dphidr*y/r;
					vx += Fx * dti;
					vy += Fy * dti;
				}
				//  step
				x += vx * dti* 3.2407764868054621e-17;
				y += vy * dti* 3.2407764868054621e-17;
				// half a kick
				r = sqrt(x*x+y*y);
				dphidr = profile_model->dphidr(r)* 3.2407764868054621e-17;
				Fx = -dphidr*x/r;
				Fy = -dphidr*y/r;
				vx += Fx * hdti;
				vy += Fy * hdti;
				
				/*
				xout (i,j) = x;
				yout (i,j) = y;
				vxout(i,j) = vx;
				vyout(i,j) = vy;
				/*/
				xoutp[i*N+j] = x;
				youtp[i*N+j] = y;
				vxoutp[i*N+j] = vx;
				vyoutp[i*N+j] = vy;
				/**/
			}
		}
	}
	ProfileModel* profile_model;
	int M;
};

class OrbitIntegratorSetLeapFrog3d : public OrbitIntegratorSet {
public:
OrbitIntegratorSetLeapFrog3d(ProfileModel3d* profile_model, int M) : profile_model(profile_model), M(M) {}
	void integrate(int N, double_vector dt, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix yout, double_matrix zout, double_matrix vxout, double_matrix vyout, double_matrix vzout) {
		//float_matrix xout(N, x0.size());
		//float_matrix yout(N, x0.size());
		//float_matrix vxout(N, x0.size());
// 		//float_matrix vyout(N, x0.size());
		double* xoutp = xout.data().begin();
		double* youtp = yout.data().begin();
		double* zoutp = zout.data().begin();
		double* vxoutp = vxout.data().begin();
		double* vyoutp = vyout.data().begin();
		double* vzoutp = vzout.data().begin();
		int no_particles = dt.size();
		//float G = 4.30200406461e-06;
		//double dphidr; 
		for(int i = 0; i < no_particles; i++) {
			//double t = 0;
			double  x = q0(i,0);
			double  y = q0(i,1);
			double  z = q0(i,2);
			double vx = v0(i,0);
			double vy = v0(i,1);
			double vz = v0(i,2);
			double dti = dt(i)/M;
			
			double hdti = dti / 2;
			//tie(Fx, Fy) = profile->dphidxy(x,x y);
			//double dphidr;

			for(int j = 0; j < N; j++) {
				// half a kick
				vx += hdti * profile_model->dphidx(x, y, z)* 3.2407764868054621e-17;
				vy += hdti * profile_model->dphidy(x, y, z)* 3.2407764868054621e-17;
				vz += hdti * profile_model->dphidz(x, y, z)* 3.2407764868054621e-17;
				
				// in between steps can be full kicks and steps
				for(int k = 0; k < (M-1); k++) {
					// step
					x += vx * dti* 3.2407764868054621e-17;
					y += vy * dti* 3.2407764868054621e-17;
					z += vz * dti* 3.2407764868054621e-17;
					// kick
					vx += dti * profile_model->dphidx(x, y, z)* 3.2407764868054621e-17;
					vy += dti * profile_model->dphidy(x, y, z)* 3.2407764868054621e-17;
					vz += dti * profile_model->dphidz(x, y, z)* 3.2407764868054621e-17;
				}
				//  step
				x += vx * dti* 3.2407764868054621e-17;
				y += vy * dti* 3.2407764868054621e-17;
				z += vz * dti* 3.2407764868054621e-17;
				// half a kick
				vx += hdti * profile_model->dphidx(x, y, z)* 3.2407764868054621e-17;
				vy += hdti * profile_model->dphidy(x, y, z)* 3.2407764868054621e-17;
				vz += hdti * profile_model->dphidz(x, y, z)* 3.2407764868054621e-17;
				xoutp[i*N+j] = x;
				youtp[i*N+j] = y;
				zoutp[i*N+j] = z;
				vxoutp[i*N+j] = vx;
				vyoutp[i*N+j] = vy;
				vzoutp[i*N+j] = vz;
				/**/
			}
		}
	}
	ProfileModel3d* profile_model;
	int M;
};



class OrbitIntegratorGSL : public OrbitIntegratorSet {
public:
OrbitIntegratorGSL(ProfileModel* profile_model, double error_abs=1e-5, double error_rel=0, const gsl_odeiv_step_type* gsl_odeiv_step_type_ = gsl_odeiv_step_rk8pd) : profile_model(profile_model), error_abs(error_abs), error_rel(error_rel), gsl_odeiv_step_type_(gsl_odeiv_step_type_) {}
	void integrate(int N, double_vector dt, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix yout, double_matrix vxout, double_matrix vyout);
	ProfileModel* profile_model;
	double error_abs;
	//void integrate(int N, double dt, tuple<double, double, double> q0, tuple<double, double, double> p0, ...)
	double error_rel;
	const gsl_odeiv_step_type* gsl_odeiv_step_type_;
};

class OrbitIntegratorAxiGSL : public OrbitIntegratorSet {
public:
OrbitIntegratorAxiGSL(ProfileModelAxi* profile_model, double error_abs=1e-5, double error_rel=0, const gsl_odeiv_step_type* gsl_odeiv_step_type_ = gsl_odeiv_step_rk8pd) : profile_model(profile_model), error_abs(error_abs), error_rel(error_rel), gsl_odeiv_step_type_(gsl_odeiv_step_type_) {}
	void integrate(int N, double_vector dt, double_vector Lz, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix out, double_matrix vxout, double_matrix vzout);
	ProfileModelAxi* profile_model;
	double error_abs;
	double error_rel;
	const gsl_odeiv_step_type* gsl_odeiv_step_type_;
};

class OrbitIntegratorGSL3d  {
public:
OrbitIntegratorGSL3d(ProfileModel* profile_model, double error_abs=1e-5, double error_rel=0, const gsl_odeiv_step_type* gsl_odeiv_step_type_ = gsl_odeiv_step_rk8pd) : profile_model(profile_model), error_abs(error_abs), error_rel(error_rel), gsl_odeiv_step_type_(gsl_odeiv_step_type_) {} 
	void integrate(int N, double_vector dt, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix yout, double_matrix zout, double_matrix vxout, double_matrix vyout, double_matrix vzout);
	ProfileModel* profile_model;
	double error_abs;
	double error_rel;
	const gsl_odeiv_step_type* gsl_odeiv_step_type_;
};

class OrbitIntegrator2dGSL {
public:
	OrbitIntegrator2dGSL(ProfileModel2d* profile_model, double error_abs=1e-5, double error_rel=0, const gsl_odeiv_step_type* gsl_odeiv_step_type_ = gsl_odeiv_step_rk8pd) : profile_model(profile_model), error_abs(error_abs), error_rel(error_rel), gsl_odeiv_step_type_(gsl_odeiv_step_type_) {} 
	void integrate(int N, double angular_velocity, double_vector dt, double_matrix q0, double_matrix v0, double_matrix xout, double_matrix yout, double_matrix vxout, double_matrix vyout);
	ProfileModel2d* profile_model;
	double error_abs;
	double error_rel;
	const gsl_odeiv_step_type* gsl_odeiv_step_type_;
};


}


