# -*- coding: utf-8 -*-
#return
import numpy
from numpy import *
try:
	import pycuda.autoinit
	import pycuda.driver as drv
	import numpy
	from pycuda.compiler import SourceModule
except:
	pass
from mab.astrounits import *

orbit_integrator_cuda_source_template = """
/*const float G = %%(G)f;
const float stellarmass = %%(stellarmass)f;
const float dmmass = %%(stellarmass)f;
const float stellar_scale = %%(stellar_scale)f;
const float dm_scale = %%(dm_scale)f;*/
	
__global__ void integrate_orbit_ME(int N, int no_particles, float* dts, float* x0, float* y0, float* z0, float* vx0, float* vy0, float* vz0, float* xout, float* yout, float* zout, float* vxout, float* vyout, float* vzout)
{
	int nr = threadIdx.x + blockDim.x*(blockIdx.x+ gridDim.x*blockIdx.y);
	if(nr >= no_particles)
		return;
	
	float x = x0[nr];
	float y = y0[nr];
	float z = z0[nr];
	float vx = vx0[nr];
	float vy = vy0[nr];
	float vz = vz0[nr];
	int index = nr * N;
	float dt = dts[nr];
	
	float r = sqrt(x*x+y*y+z*z);
	float dphidr = %(dphidr)s;
	dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
	//G * mass * r / pow((r*r + scale*scale), (3.f/2));
	float Fx = -dphidr*x/r;
	float Fy = -dphidr*y/r;
	float Fz = -dphidr*z/r;
	//vx += Fx * dt/2;
	//vy += Fy * dt/2;
	//vz += Fz * dt/2;
	
	for(int j = 0; j < N; j++) {
		for(int k = 0; k < %(M)d; k++) {
			r = sqrt(x*x+y*y+z*z);
			dphidr = %(dphidr)s;
			dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
			//G * mass * r / pow((r*r + scale*scale), (3.f/2));
			Fx = -dphidr*x/r;
			Fy = -dphidr*y/r;
			Fz = -dphidr*z/r;

			vx += Fx * dt;
			vy += Fy * dt;
			vz += Fz * dt;
			
			x += vx * dt * %(kpc_per_km)r;
			y += vy * dt * %(kpc_per_km)r;
			z += vz * dt * %(kpc_per_km)r;
			
			
		}
		xout[index] = x;
		yout[index] = y;
		zout[index] = z;
		vxout[index] = vx;// - Fx * dt/2;
		vyout[index] = vy;// - Fy * dt/2;
		vzout[index] = vz;// - Fz * dt/2;
		index++;
	}
	
}

__global__ void integrate_orbit_LF(int N, int no_particles, float* dts, float* x0, float* y0, float* z0, float* vx0, float* vy0, float* vz0, float* xout, float* yout, float* zout, float* vxout, float* vyout, float* vzout)
{
	int nr = threadIdx.x + blockDim.x*(blockIdx.x+ gridDim.x*blockIdx.y);
	if(nr >= no_particles)
		return;
	
	float x = x0[nr];
	float y = y0[nr];
	float z = z0[nr];
	float vx = vx0[nr];
	float vy = vy0[nr];
	float vz = vz0[nr];
	int index = nr * N;
	float dt = dts[nr];
	
	float r = sqrt(x*x+y*y+z*z);
	float dphidr = %(dphidr)s;
	dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
	//G * mass * r / pow((r*r + scale*scale), (3.f/2));
	float Fx = -dphidr*x/r;
	float Fy = -dphidr*y/r;
	float Fz = -dphidr*z/r;
	vx += Fx * dt/4;
	vy += Fy * dt/4;
	vz += Fz * dt/4;
	//x += vx * dt/2 * %(kpc_per_km)r;
	//y += vy * dt/2 * %(kpc_per_km)r;
	//z += vz * dt/2 * %(kpc_per_km)r;
	
	for(int j = 0; j < N; j++) {
		for(int k = 0; k < %(M)d; k++) {
			x += vx * dt/2 * %(kpc_per_km)r;
			y += vy * dt/2 * %(kpc_per_km)r;
			z += vz * dt/2 * %(kpc_per_km)r;
			
			r = sqrt(x*x+y*y+z*z);
			dphidr = %(dphidr)s;
			dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
			//G * mass * r / pow((r*r + scale*scale), (3.f/2));
			Fx = -dphidr*x/r;
			Fy = -dphidr*y/r;
			Fz = -dphidr*z/r;

			vx += Fx * dt;
			vy += Fy * dt;
			vz += Fz * dt;
			
			x += vx * dt/2 * %(kpc_per_km)r;
			y += vy * dt/2 * %(kpc_per_km)r;
			z += vz * dt/2 * %(kpc_per_km)r;
			
			
		}
		xout[index] = x;// - vx * dt/4 * %(kpc_per_km)r;
		yout[index] = y;// - vx * dt/4 * %(kpc_per_km)r;
		zout[index] = z;// - vx * dt/4 * %(kpc_per_km)r;
		vxout[index] = vx;// - Fx * dt/2;
		vyout[index] = vy;// - Fy * dt/2;
		vzout[index] = vz;// - Fz * dt/2;
		index++;
	}
	
}

__global__ void integrate_orbit_LF3(int N, int no_particles, float dt, float* x0, float* y0, float* z0, float* vx0, float* vy0, float* vz0, float* xout, float* yout, float* zout, float* vxout, float* vyout, float* vzout)
{
	int nr = threadIdx.x + blockDim.x*(blockIdx.x+ gridDim.x*blockIdx.y);
	if(nr >= no_particles)
		return;
	
	float x = x0[nr];
	float y = y0[nr];
	float z = z0[nr];
	float vx = vx0[nr];
	float vy = vy0[nr];
	float vz = vz0[nr];
	int index = nr * N;
	
	float r = sqrt(x*x+y*y+z*z);
	float dphidr = %(dphidr)s;
	dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
	//G * mass * r / pow((r*r + scale*scale), (3.f/2));
	float Fx = -dphidr*x/r;
	float Fy = -dphidr*y/r;
	float Fz = -dphidr*z/r;
	//vx += Fx * dt/2;
	//vy += Fy * dt/2;
	//vz += Fz * dt/2;
	//x += vx * dt/2 * %(kpc_per_km)r;
	//y += vy * dt/2 * %(kpc_per_km)r;
	//z += vz * dt/2 * %(kpc_per_km)r;
	
	for(int j = 0; j < N; j++) {
		for(int k = 0; k < %(M)d; k++) {
		
			// step 1
			r = sqrt(x*x+y*y+z*z);
			dphidr = %(dphidr)s;
			dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
			//G * mass * r / pow((r*r + scale*scale), (3.f/2));
			Fx = -dphidr*x/r;
			Fy = -dphidr*y/r;
			Fz = -dphidr*z/r;
			vx += Fx * dt * 7./24;
			vy += Fy * dt * 7./24;
			vz += Fz * dt * 7./24;

			x += vx * dt*2/3 * %(kpc_per_km)r;
			y += vy * dt*2/3 * %(kpc_per_km)r;
			z += vz * dt*2/3 * %(kpc_per_km)r;
			
		
			// step 2
			r = sqrt(x*x+y*y+z*z);
			dphidr = %(dphidr)s;
			dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
			//G * mass * r / pow((r*r + scale*scale), (3.f/2));
			Fx = -dphidr*x/r;
			Fy = -dphidr*y/r;
			Fz = -dphidr*z/r;
			vx += Fx * dt * 3./4;
			vy += Fy * dt * 3./4;
			vz += Fz * dt * 3./4;

			x += vx * dt*-2/3 * %(kpc_per_km)r;
			y += vy * dt*-2/3 * %(kpc_per_km)r;
			z += vz * dt*-2/3 * %(kpc_per_km)r;
			
		
			// step 3
			r = sqrt(x*x+y*y+z*z);
			dphidr = %(dphidr)s;
			dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
			//G * mass * r / pow((r*r + scale*scale), (3.f/2));
			Fx = -dphidr*x/r;
			Fy = -dphidr*y/r;
			Fz = -dphidr*z/r;
			vx += Fx * dt * -1./24;
			vy += Fy * dt * -1./24;
			vz += Fz * dt * -1./24;

			x += vx * dt* %(kpc_per_km)r;
			y += vy * dt* %(kpc_per_km)r;
			z += vz * dt* %(kpc_per_km)r;
			
		}
		xout[index] = x;// - vx * dt/2 * %(kpc_per_km)r;
		yout[index] = y;// - vx * dt/2 * %(kpc_per_km)r;
		zout[index] = z;// - vx * dt/2 * %(kpc_per_km)r;
		vxout[index] = vx;// - Fx * dt/2;
		vyout[index] = vy;// - Fy * dt/2;
		vzout[index] = vz;// - Fz * dt/2;
		index++;
	}
	
}
const float LF4_a1 = 0.67560359597982889;
const float LF4_a2 = -0.17560359597982886;
const float LF4_a3 = LF4_a2;
const float LF4_a4 = LF4_a1;

//const float LF4_b1 = 0;
const float LF4_b2 = 1.3512071919596578;
const float LF4_b3 = -1.7024143839193153;
const float LF4_b4 = LF4_b2;

__global__ void integrate_orbit_LF4(int N, int no_particles, float* dts, float* x0, float* y0, float* z0, float* vx0, float* vy0, float* vz0, float* xout, float* yout, float* zout, float* vxout, float* vyout, float* vzout, float* maxEdiff)
{
	int nr = threadIdx.x + blockDim.x*(blockIdx.x+ gridDim.x*blockIdx.y);
	if(nr >= no_particles)
		return;
	
	float x = x0[nr];
	float y = y0[nr];
	float z = z0[nr];
	float vx = vx0[nr];
	float vy = vy0[nr];
	float vz = vz0[nr];
	float dt = dts[nr];
	int index = nr * N;
	
	float r = sqrt(x*x+y*y+z*z);
	float dphidr = %(dphidr)s;
	dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
	float Fx = -dphidr*x/r;
	float Fy = -dphidr*y/r;
	float Fz = -dphidr*z/r;
	
	float E0 = 0.5 * (vx*vx + vy*vy + vz*vz) + %(Epot)s;
	float E, Ediff;
	
	for(int j = 0; j < N; j++) {
		for(int k = 0; k < %(M)d; k++) {
		
			// step 1
			x += vx * dt*LF4_a1 * %(kpc_per_km)r;
			y += vy * dt*LF4_a1 * %(kpc_per_km)r;
			z += vz * dt*LF4_a1 * %(kpc_per_km)r;
			
		
			// step 2
			r = sqrt(x*x+y*y+z*z);
			dphidr = %(dphidr)s;
			dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
			//G * mass * r / pow((r*r + scale*scale), (3.f/2));
			Fx = -dphidr*x/r;
			Fy = -dphidr*y/r;
			Fz = -dphidr*z/r;
			vx += Fx * dt * LF4_b2;
			vy += Fy * dt * LF4_b2;
			vz += Fz * dt * LF4_b2;

			x += vx * dt*LF4_a2 * %(kpc_per_km)r;
			y += vy * dt*LF4_a2 * %(kpc_per_km)r;
			z += vz * dt*LF4_a2 * %(kpc_per_km)r;
			
			// step 3
			r = sqrt(x*x+y*y+z*z);
			dphidr = %(dphidr)s;
			dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
			//G * mass * r / pow((r*r + scale*scale), (3.f/2));
			Fx = -dphidr*x/r;
			Fy = -dphidr*y/r;
			Fz = -dphidr*z/r;
			vx += Fx * dt * LF4_b3;
			vy += Fy * dt * LF4_b3;
			vz += Fz * dt * LF4_b3;

			x += vx * dt*LF4_a3 * %(kpc_per_km)r;
			y += vy * dt*LF4_a3 * %(kpc_per_km)r;
			z += vz * dt*LF4_a3 * %(kpc_per_km)r;
			
			// step 4
			r = sqrt(x*x+y*y+z*z);
			dphidr = %(dphidr)s;
			dphidr *= %(kpc_per_km)r; // go from km/kpc km/s**2 to km/s**2 
			//G * mass * r / pow((r*r + scale*scale), (3.f/2));
			Fx = -dphidr*x/r;
			Fy = -dphidr*y/r;
			Fz = -dphidr*z/r;
			vx += Fx * dt * LF4_b4;
			vy += Fy * dt * LF4_b4;
			vz += Fz * dt * LF4_b4;

			x += vx * dt*LF4_a4 * %(kpc_per_km)r;
			y += vy * dt*LF4_a4 * %(kpc_per_km)r;
			z += vz * dt*LF4_a4 * %(kpc_per_km)r;
			
			float r = sqrt(x*x+y*y+z*z);
			E = 0.5 * (vx*vx + vy*vy + vz*vz) + %(Epot)s;
			Ediff = abs(E-E0);
			
			if(((j == 0) && (k == 0)) || (Ediff > maxEdiff[nr]))
				maxEdiff[nr] = Ediff;
			
		}
		xout[index] = x;// - vx * dt/2 * %(kpc_per_km)r;
		yout[index] = y;// - vx * dt/2 * %(kpc_per_km)r;
		zout[index] = z;// - vx * dt/2 * %(kpc_per_km)r;
		vxout[index] = vx;// - Fx * dt/2;
		vyout[index] = vy;// - Fy * dt/2;
		vzout[index] = vz;// - Fz * dt/2;
		index++;
	}
	
}
"""

#class OrbitIntegrator(object):
#	def
 
def integrate_orbit_cuda_LF(stellar_profile, dm_profile, x, y, z, vx, vy, vz, dt, N, M, splits=1, debug=False, filename="cudacode.c", usefloat=False):
	variables = {}
	variables["dphidr"] = stellar_profile.c_code_dphidr()
	variables["Epot"] = stellar_profile.c_code_potential()
	if dm_profile.M > 1:
		if debug:
			print "including DM halo"
		variables["dphidr"] += "+" + dm_profile.c_code_dphidr()
		variables["Epot"] += " + " + dm_profile.c_code_potential()
	variables["M"] = M
	variables["kpc_per_km"] = (KM/KPC).asNumber()
	if debug:
		print "kpc per km =", variables["kpc_per_km"] 
	orbit_integrator_cuda_source = orbit_integrator_cuda_source_template % variables
	if debug:
		f = open(filename, "w")
		print >>f, orbit_integrator_cuda_source
		f.close() 
	module = SourceModule(orbit_integrator_cuda_source)
	integrate_orbit = module.get_function("integrate_orbit_LF4")
	#numpy.int32(no_particles), numpy.float64(dt/M), drv.In(x0), drv.In(y0), drv.In(vx0), drv.In(vy0), drv.Out(xout), drv.Out(yout), drv.Out(vxout)
	no_particles = len(x)
	if usefloat:
		x = x.astype(float32)
		y = y.astype(float32)
		z = z.astype(float32)
		vx = vx.astype(float32)
		vy = vy.astype(float32)
		vz = vz.astype(float32)
		dt = dt.astype(float32)
		xout = zeros((len(x), N), dtype=float32)
	else:
		xout = zeros((len(x), N), dtype=float64)
	yout = zeros_like(xout)
	zout = zeros_like(xout)
	vxout = zeros_like(xout)
	vyout = zeros_like(xout)
	vzout = zeros_like(xout)
	maxEdiff = zeros_like(dt)
	particles_per_split = int(math.ceil(no_particles/double(splits)))
	if debug:
		print particles_per_split, "particles per split"
	for i in range(splits):
		no_particles_i = particles_per_split
		index1 = i * particles_per_split
		index2 = (i+1) * particles_per_split
		if index2 > no_particles:
			index2 = no_particles
		no_particles_i = index2 - index1 
		gridsizex = 32
		gridsizey = no_particles_i/(32*gridsizex)
		if gridsizey * gridsizex * 32 < no_particles_i:
			gridsizey += 1 
		assert gridsizey * gridsizex * 32 >= no_particles_i
		grid = (gridsizex, gridsizey)
		if debug:
			print "step", i, "from", index1, "to", index2
			#print z
			print grid, gridsizey * gridsizex * 32, no_particles_i
		integrate_orbit(numpy.int32(N), numpy.int32(no_particles_i), drv.In(dt[index1:index2]), drv.In(x[index1:index2]), drv.In(y[index1:index2]), drv.In(z[index1:index2]), drv.In(vx[index1:index2]), drv.In(vy[index1:index2]), drv.In(vz[index1:index2]), drv.Out(xout[index1:index2,:]), drv.Out(yout[index1:index2,:]), drv.Out(zout[index1:index2,:]), drv.Out(vxout[index1:index2,:]), drv.Out(vyout[index1:index2,:]), drv.Out(vzout[index1:index2,:]), drv.InOut(maxEdiff[index1:index2]), block=(32,1,1), grid=grid)
	return xout, yout, zout, vxout, vyout, vzout, maxEdiff
 
 
 
	
