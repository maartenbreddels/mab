#define M 1

#define J 500
#define K 100
#define G 4.30200406461e-06
#define mass 1e6
#define scale 1.


__global__ void integrate_orbit_euler(float* x_out, float* y_out, float* vx_out, float *vy_out)
{
	int nr = threadIdx.x + blockDim.x*(blockIdx.x + blockIdx.y*gridDim.x);
	float x = x_out[nr * J];
	float y = y_out[nr * J];
	float vx = vx_out[nr * J];
	float vy = vy_out[nr * J];
	float dt = 0.0001;
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
		x_out[nr*J+j] = x;
		y_out[nr*J+j] = y;
		vx_out[nr*J+j] = vx;
		vy_out[nr*J+j] = vy;
	}
	
}


__global__ void square(float* x)
{
	int i = threadIdx.x + blockDim.x*(blockIdx.x + blockIdx.y*gridDim.x);
	for(int j = 0; j < M; j++)
		x[i] = x[i] + M;
	//x[i] = log(x[i]+1) + x[i] * x[i] + sqrt(x[i]/(x[i]+100)); // * 2; //x[i];
}

