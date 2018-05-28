#include <stdio.h>
#include "orbit_integrator_cuda.cu"

#define N 512 * 512 * 64
float x[N];
float *x_device;
cudaError_t err;

int main(int argc, char** argv) {
	for(int i = 0; i < N; i++) {
		x[i] = i;
	}
	printf("x[%d] = %f\n", N-1, x[N-1]);
	
	if( argc > 1) {
		cudaMalloc((void**) &x_device, sizeof(float)*N);
		err = cudaGetLastError ();
		printf("malloc: %s\n", cudaGetErrorString(err));
		
		
		cudaMemcpy(x_device, x, sizeof(float)*N, cudaMemcpyHostToDevice);
		err = cudaGetLastError ();
		printf("copy: %s\n", cudaGetErrorString(err));
		
		dim3 dimBlock(512,64);
		square<<<dimBlock, 512>>>(x_device);
		err = cudaGetLastError ();
		printf("call: %s\n", cudaGetErrorString(err));
		
		cudaMemcpy(x, x_device, sizeof(float)*N, cudaMemcpyDeviceToHost);
		err = cudaGetLastError ();
		printf("cpy: %s\n", cudaGetErrorString(err));
	} else {
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < M; j++)
				x[i] = x[i] + M;
		}
	}
	
	for(int i = 0; i < 10; i++) {
		printf("x[%d] = %f\n", i, x[i]);
	}
	printf("x[%d] = %f\n", N-1, x[N-1]);
	return 0;
}
	
