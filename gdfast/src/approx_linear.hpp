#pragma once
#include "nlinear.hpp"

template<typename Grid, typename T=double>
class ApproxLinear {
public:
	ApproxLinear(Grid* grid, double x0, double x1) : grid(grid), nodes(grid->length()+1), x0(x0), x1(x1) {
		nodevalues = new T[nodes];
	}

	void bin(double_vector x) {
	}
	void bin(double_vector x, double_vector data, double_vector binned) {
		double* xp = x.data().begin();
		int size = x.size();
		//NLinear<1, double> nlinear
		for(int i=0; i < size; i++) {
			double xi = xp[i];
			if(grid->inrange(xi)) {
				int index = grid->findindex(xi);
				weight1 = (1 - abs(x));
				
			}
		}
		
	}

	T* nodevalues;
	int nodes;
	Grid* grid;
};