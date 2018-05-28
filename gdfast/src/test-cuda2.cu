#include <stdio.h>
#include "orbit_integrator_cuda.cu"

#define N 256
#define N_TOT N * J
float x_h[N_TOT], y_h[N_TOT], vx_h[N_TOT], vy_h[N_TOT];
float *x_d, *y_d, *vx_d, *vy_d;
cudaError_t err;

int main(int argc, char** argv) {
	for(int i = 0; i < N; i++) {
		x_h[i*J] = 1;
		y_h[i*J] = 0;
		vx_h[i*J] = 0;
		vy_h[i*J] = 0.1;
	}
	//printf("x[%d] = %f\n", N-1, x[N-1]);
	
	if( argc > 1) {
		cudaMalloc((void**) &x_d, sizeof(float)*N_TOT);
		cudaMalloc((void**) &y_d, sizeof(float)*N_TOT);
		cudaMalloc((void**) &vx_d, sizeof(float)*N_TOT);
		cudaMalloc((void**) &vy_d, sizeof(float)*N_TOT);
		err = cudaGetLastError ();
		printf("malloc: %s\n", cudaGetErrorString(err));
		
		
		cudaMemcpy(x_d, x_h, sizeof(float)*N_TOT, cudaMemcpyHostToDevice);
		err = cudaGetLastError ();
		printf("copy: %s\n", cudaGetErrorString(err));
		cudaMemcpy(y_d, y_h, sizeof(float)*N_TOT, cudaMemcpyHostToDevice);
		err = cudaGetLastError ();
		printf("copy: %s\n", cudaGetErrorString(err));
		cudaMemcpy(vx_d, vx_h, sizeof(float)*N_TOT, cudaMemcpyHostToDevice);
		err = cudaGetLastError ();
		printf("copy: %s\n", cudaGetErrorString(err));
		cudaMemcpy(vy_d, vy_h, sizeof(float)*N_TOT, cudaMemcpyHostToDevice);
		err = cudaGetLastError ();
		printf("copy: %s\n", cudaGetErrorString(err));
		
		dim3 dimBlock(512,64);
		//square<<<dimBlock, 512>>>(x_device);
		integrate_orbit_euler<<<32, N/32>>>(x_d, y_d, vx_d, vy_d);
		err = cudaGetLastError ();
		printf("call: %s\n", cudaGetErrorString(err));
		//*
		cudaMemcpy(x_h, x_d, sizeof(float)*N_TOT, cudaMemcpyDeviceToHost);
		cudaMemcpy(y_h, y_d, sizeof(float)*N_TOT, cudaMemcpyDeviceToHost);
		cudaMemcpy(vx_h, vx_d, sizeof(float)*N_TOT, cudaMemcpyDeviceToHost);
		cudaMemcpy(vy_h, vy_d, sizeof(float)*N_TOT, cudaMemcpyDeviceToHost);
		/**/
		err = cudaGetLastError ();
		printf("cpy: %s\n", cudaGetErrorString(err));
	} else {
		for(int i = 0; i < N; i++) {
			int nr = i;
			float x = x_h[nr * J];
			float y = y_h[nr * J];
			float vx = vx_h[nr * J];
			float vy = vy_h[nr * J];
			float dt = 0.01;
			for(int j = 0; j < J; j++) {
				for(int k = 0; k < K; k++) {
					float r = sqrt(x*x+y*y);
					float dphidr = G * mass * r / pow((r*r + scale*scale), (3./2));
					float Fx = -dphidr*x/r;
					float Fy = -dphidr*y/r;
					
					vx += Fx * dt;
					x += vx * dt;
					
					vy += Fy * dt;
					y += vy * dt;
					
				}
				x_h[nr*J+j] = x;
				y_h[nr*J+j] = y;
				vx_h[nr*J+j] = vx;
				vy_h[nr*J+j] = vy;
			}
		}
	}
	
	for(int i = 0; i < 1; i++) {
		printf(" x[%d] = %f\n", i, x_h[J*i]);
		printf(" y[%d] = %f\n", i, y_h[J*i]);
		printf("vx[%d] = %f\n", i, vx_h[J*i]);
		printf("vy[%d] = %f\n", i, vy_h[J*i]);
	}
	//printf("x[%d] = %f\n", N-1, x[N-1]);
	return 0;
}
	

