#include <stdio.h>
#include <cmath>
#include <assert.h>

#undef SG_DEBUG

class BasisTriangular {// : public BasisSet {
public:
	BasisTriangular() {}
	int count(int level) {
		return level;
	}
	double center(int level, int i) {
		return 1./(level-1) * i;
	}
	bool inrange(int level, int i, double x) {
		return true;
	}
	double operator()(int level, int i, double x) {
		int count = this->count(level);
		const double hl = 1./(count-1);// pow(2, -level);
		double center = this->center(level, i);
		if((x > 1) || (x < 0)) // outside range?
			return 0;
		else {
			return x >= (center+hl) ? 0 : (
					x <= (center-hl) ? 0 : 1-fabs((x-center)/hl)
				);
			/*return x >= (j+1)*hl ? 0 : (
					x < (j-1)*hl ? 0 : 1-fabs(x/hl-j)
				);*/
		}
	}
};

class BasisSet {
public:
	virtual int count(int level) = 0;
	virtual double operator()(int level, int i, double x) = 0;
	virtual double center(int level, int i) = 0;
};

class BasisSetHierTriangular {// : public BasisSet {
public:
	BasisSetHierTriangular() {}
	int count(int level) {
		if(level == 0) return 2;
		else return 1<<(level-1);
	}
	double weight(int level) {
		if(level == 0) return 0.5;
		else return 0.5/(1<<(level-1));
	}
	double center(int level, int i) {
		if(level == 0) {
			return i == 0 ? 0 : 1;
		} else {
			int count = 1<<(level-1);;
			//int count = this->count(level);
			const double hl = 1./count;// pow(2, -level);
			//int j = i * 2 + 1;
			return (i+0.5) * hl;
		}
	}
	bool inrange(int level, int i, double x) {
		if(level == 0) {
			return ((x <= 1) && (x >= 0));
		} else {
			int count = 1<<(level-1);;
			const double hl = 1./count;// pow(2, -level);
			double center = this->center(level, i);
			if((x > 1) || (x < 0)) // outside range?
				return false;
			else {
				return x >= (center+hl/2) ? false : 
					(
						x <= (center-hl/2) ? false : true
					);
			}
		}
	}
	double operator()(int level, int i, double x) {
		if(level > 0) {
			int count = 1<<(level-1);;
			const double hl = 1./count;// pow(2, -level);
			double center = this->center(level, i);
			//int j = i * 2 + 1;
			//printf("[l=%d i=%d x=%f hl=%f c=%f]", level, i, x, hl, center);
			if((x > 1) || (x < 0)) // outside range?
				return 0;
			else {
			return x >= (center+hl/2) ? 0 : (
			                                  x <= (center-hl/2) ? 0 : 1-fabs((x-center)/hl*2)
			                                );
				/*return x >= (j+1)*hl ? 0 : (
						x < (j-1)*hl ? 0 : 1-fabs(x/hl-j)
					);*/
			}
		} else {
			if((x > 1) || (x < 0)) // outside range or odd, skip
				return 0;
			else {
				if(i == 0)
					return x > 1 ? 0 : (x < 0 ? 0 : 1-x);
				else if(i == 1)
					return x > 1 ? 0 : (x < 0 ? 0 : x);
				else
					return 0;
			}
		}
	}
};


class BasisSetHierTriangularM {// : public BasisSet {
public:
	BasisSetHierTriangularM() {}
	int count(int level) {
		if(level < 2) return 1<<level;
		else return 1<<(level-1);
	}
	double weight(int level) {
		return 0.5 / count(level);
	}
	double center(int level, int i) {
		if(level == 0) {
			return 0.5;
		} else if(level == 1) {
			return i == 0 ? 0. : 1.;
		} else {
			int count = 1<<(level-1);
			//int count = this->count(level);
			const double hl = 1./count;// pow(2, -level);
			//int j = i * 2 + 1;
			return (i+0.5) * hl;
		}
	}
	/*bool inrange(int level, int i, double x) {
		if(level == 0) {
			return ((x <= 1) && (x >= 0));
		} else {
			int count = 1<<(level-1);;
			const double hl = 1./count;// pow(2, -level);
			double center = this->center(level, i);
			if((x > 1) || (x < 0)) // outside range?
				return false;
			else {
			return x >= (center+hl/2) ? false : 
				(
				  x <= (center-hl/2) ? false : true
				);
			}
		}
	}*/
	double operator()(int level, int i, double x) {
		if(level != 1) {
			int count = this->count(level); //1<<(level-1);
			const double hl = 1./count;// pow(2, -level);
			double center = this->center(level, i);
			//int j = i * 2 + 1;
			//printf("[l=%d i=%d x=%f hl=%f c=%f]", level, i, x, hl, center);
			if((x > 1) || (x < 0)) // outside range?
				return 0;
			else {
				return x >= (center+hl/2) ? 0 : (
			                                  x <= (center-hl/2) ? 0 : 1-fabs((x-center)/hl*2)
			                                );
				/*return x >= (j+1)*hl ? 0 : (
						x < (j-1)*hl ? 0 : 1-fabs(x/hl-j)
					);*/
			}
		} else { // level == 1
			if((x > 1) || (x < 0)) // outside range or odd, skip
				return 0;
			else {
				if(i == 0)
					return x > 0.5 ? 0 : (x < 0 ? 0 : 1-x*2);
				else if(i == 1)
					return x > 1 ? 0 : (x < 0.5 ? 0 : x*2-1);
				else
					return 0;
			}
		}
	}
	};

class SparseGrid {
public:
	double *values;
	int count;
	int dim, max_level;
	BasisSetHierTriangularM basis_set;
	SparseGrid(int dim, int max_level) : count(0), dim(dim), max_level(max_level) {
		auto f = [&](int level, int* levels) {
			int levelcount = 1;
			for(int d = 0; d < dim; d++) {
				levelcount *= this->basis_set.count(levels[d]);
			}
			/*
			int norm1;
			for(int d = 0; d < dim; d++) {
				printf("%d ", levels[d]);
				norm1 += levels[d];
			}
			printf(" norm1=%d", norm1);
			printf(" levelcount=%d\n", levelcount);
			/**/
			this->count += levelcount;
			
		};
		loop(f);
		printf("count=%d\n", count);
		values = new double[count];
		for(int i = 0; i < count; i++) values[i] = 0;
	}
	SparseGrid marginalize1d(int dim_m) {
		SparseGrid marginalized = SparseGrid(1, max_level);
		//double value = 0;
		int value_index = 0;
		int value_index1d = 0;
		//int level1d = 0;
		//int index1d = 0;
		auto f = [&](int level, int* levels, int* indices, double* point) {
			value_index1d = 0;
			auto g = [&](int level1d, int* levels1d, int* indices1d, double* point1d) {
				if((levels1d[0] == levels[dim_m]) && (indices1d[0] == indices[dim_m])) {
					double w = 1;
					for(int d = 0; d < dim; d++) {
						if(d != dim_m) {
							w *= basis_set.weight(levels[d]);
						}
					}
					marginalized.values[value_index1d] += values[value_index] * w;
				}
				value_index1d++;
			};
			marginalized.loop_indices(g);
			value_index++;
		};
		loop_indices(f);
		return marginalized;
	}
	
	SparseGrid marginalize2d(int dim1_m, int dim2_m) {
		SparseGrid marginalized = SparseGrid(2, max_level);
		//double value = 0;
		int value_index = 0;
		int value_index2d = 0;
		//int level1d = 0;
		//int index1d = 0;
		auto f = [&](int level, int* levels, int* indices, double* point) {
			value_index2d = 0;
			auto g = [&](int level2d, int* levels2d, int* indices2d, double* point2d) {
				if((levels2d[0] == levels[dim1_m]) && (indices2d[0] == indices[dim1_m]) && (levels2d[1] == levels[dim2_m]) && (indices2d[1] == indices[dim2_m])) {
					double w = 1;
					for(int d = 0; d < dim; d++) {
						if((d != dim1_m) && (d != dim2_m)) {
							w *= basis_set.weight(levels[d]);
						}
					}
					marginalized.values[value_index2d] += values[value_index] * w;
				}
				value_index2d++;
			};
			marginalized.loop_indices(g);
			value_index++;
		};
		loop_indices(f);
		return marginalized;
	}
	
	double integrate(int max_level) {
		double value = 0;
		int value_index = 0;
		auto f = [&](int level, int* levels, int* indices, double *point) {
			double weight = 1;
			for(int d = 0; d < dim; d++) {
				/*if(!basis_set.inrange(levels[d], indices[d], point[d])) {
					weight = 0;
					break;
				}*/
				weight *= basis_set.weight(levels[d]);
				//printf("l=%d i=%d p=%f v=%f\n", levels[d], indices[d], point[d], basis_set.operator()(levels[d], indices[d], point[d]));
			}
			//printf("v = % 10.4f w = % 10.4f\n", values[value_index], weight);
			value += values[value_index] * weight; 
			value_index++;
		};
		loop_indices(f, max_level);
		return value;
	}
	
	double operator()(double* point) {
		return call_level(max_level, point);
	}
	double call_level(int max_level, double* point) {
		double value = 0;
		//int value_index = 0;
		auto f = [&](int level, int* levels, int* indices, int value_index) {
			double weight = 1;
			for(int d = 0; d < dim; d++) {
				/*if(!basis_set.inrange(levels[d], indices[d], point[d])) {
					weight = 0;
					break;
				}*/
				weight *= basis_set.operator()(levels[d], indices[d], point[d]);
				//printf("l=%d i=%d p=%f v=%f\n", levels[d], indices[d], point[d], basis_set.operator()(levels[d], indices[d], point[d]));
			}
			if(weight != 0) {
#ifdef SG_DEBUG
				printf("use l [");
				for(int d = 0; d < dim; d++) {
					printf("%d ", indices[d]);
				}
				printf("]\n");
#endif
			}
				
			//printf("v = % 10.4f w = % 10.4f\n", values[value_index], weight);
			value += values[value_index] * weight; 
			//value_index++;
		};
		loop_indices2(f, point, max_level);
		return value;
	}

	template<class F>
	void eval(F function) {
		int value_index = 0;
		auto f = [&](int level, int* levels, int* indices, double* point) {
			values[value_index] = function(dim, point);
			//printf(" f(...) = % 10.4f", values[value_index]);
			if(level > 0) {
				double pvalue = this->call_level(level, point); 
				double surplus = values[value_index] - pvalue;
				//printf(" surplus = % 10.4f pvalue = % 10.4f", surplus, pvalue);
				values[value_index] = surplus;
			}
			//printf("\n");
			value_index++;
		};
		loop_indices(f);
		//printf("count=%d\n", count);
	}

	template<class F>
	inline void loop_indices2(F f, double* point) {
		loop_indices(f, point, max_level);
	}
	template<class F>
	inline void loop_indices2(F f, double* point, int max_level) {
		int *cumcounts = new int[dim];
		int *counts = new int[dim];
		int *indices = new int[dim];
		int *target_indices = new int[dim];
//double *point = new double[dim];
		int value_index = 0;
#ifdef SG_DEBUG
		printf("p[");
		for(int d = 0; d < dim; d++) {
			printf("%f ", point[d]);
		}
		printf("]\n");
#endif
		
		auto g = [&](int level, int* levels) {
#ifdef SG_DEBUG
			printf("level [");
			for(int d = 0; d < dim; d++) {
				printf("%d ", levels[d]);
			}
			printf("]\n");
#endif
			for(int d = 0; d < dim; d++) {
				counts[d] = this->basis_set.count(levels[d]);
				indices[d] = 0;
				target_indices[d] = (int)(counts[d] * point[d]);
				if(point[d] == 1.0)
					target_indices[d] = counts[d] - 1;
				//if(levels[d] == 0)
				//target_indices[d] = -1;
			}
			cumcounts[dim-1] = 1;
			for(int d = dim-2; d >= 0; d--) {
				cumcounts[d] = cumcounts[d+1] * counts[d+1];
			}
			/*printf("\t\tcounts[");
			for(int d = 0; d < dim; d++) {
				printf("%d ", counts[d]);
			}
			printf("]\n");
			printf("\t\tcumcounts[");
			for(int d = 0; d < dim; d++) {
				printf("%d ", cumcounts[d]);
			}
			printf("]\n");*/
#ifdef SG_DEBUG
			printf("\t\tchk[");
			for(int d = 0; d < dim; d++) {
				printf("%f %d ", point[d], counts[d]);
			}
			printf("]\n");
#endif
#ifdef SG_DEBUG
			printf("target index[");
			for(int d = 0; d < dim; d++) {
				printf("%d ", target_indices[d]);
			}
			printf("]\n");
#endif
			/*
			bool has_no_zero_level = true;
			for(int d = 0; d < dim; d++) {
				if(levels[d] == 0) {
					has_no_zero_level = false;
				}
			}*/
			//if(has_no_zero_level)
 			if(1) {
	 			//int level_index = 0; 
	 			int level_count = 1; //counts[0];
	 			int level_index = 0; //indices[0];
	 			/*for(int d = 0; d < dim; d++) {
					indices[d] = (int)(counts[d] * point[d]);
					if(indices[d] >= counts[d])
						indices[d] = counts[d]-1;
					level_count *= counts[d];
				}*/
				for(int d = 0; d < dim; d++) {
					//level_index *= counts[d-1];
					level_index += target_indices[d] * cumcounts[d];
					level_count *= counts[d];
				}
				value_index += level_index;
	 			f(level, levels, target_indices, value_index);
				value_index -= level_index;
				value_index += level_count;
				
			} else /**/
			{
	
				bool done = false;
				int overlaps = 0;
				while(!done) { // loop over indices
					//for(int d = 0; d < dim; d++) {
					//	point[d] = this->basis_set.center(levels[d], indices[d]);
					//}
					bool overlap = true;
					//* this codes does not work with the M basis
					for(int d = 0; d < dim ; d++) {
						//overlap = overlap && ((levels[d] == 0) || (indices[d] == target_indices[d]));
						//*
						//if(levels[d] != 0)
						{
							if(indices[d] != target_indices[d]) {
								overlap = false;
								break;
							}
						}/**/
					}
					/**/
					if(overlap) {
						f(level, levels, indices, value_index);
						overlaps++;
					}
					value_index++;
					for(int d = dim-1; d >= 0; d--) {
						indices[d]++;
						if(indices[d]>=counts[d]) {
							indices[d] = 0;
							if(d == 0) done = true;
						} else {
							break;
						}
					}
				}
				assert(overlaps == 1);
				//printf("%d, ", overlaps);
			}
		};
		loop(g, max_level);
		/*if(max_level == this->max_level) {
			printf("valueindex = %d count = %d\n", value_index, count);
			assert(value_index == count);
		}*/
		delete[] cumcounts;
		delete[] counts;
		delete[] indices;
		delete[] target_indices;
		//delete[] point;
	}









	template<class F>
	inline void loop_indices(F f) {
		loop_indices(f, max_level);
	}
	template<class F>
	inline void loop_indices(F f, int max_level) {
		int *counts = new int[dim];
		int *indices = new int[dim];
		double *point = new double[dim];
		//int value_index = 0;
		auto g = [&](int level, int* levels) {
			for(int d = 0; d < dim; d++) {
				counts[d] = this->basis_set.count(levels[d]);
				indices[d] = 0;
			}
			bool done = false;
			while(!done) { // loop over indices
				for(int d = 0; d < dim; d++) {
					point[d] = this->basis_set.center(levels[d], indices[d]);
				}
				f(level, levels, indices, point);
				for(int d = dim-1; d >= 0; d--) {
					indices[d]++;
					if(indices[d]>=counts[d]) {
						indices[d] = 0;
						if(d == 0) done = true;
					} else {
						break;
					}
				}
			}
		};
		loop(g, max_level);
		delete[] counts;
		delete[] indices;
		delete[] point;
	}
	template<class F>
	inline void loop(F f) {
		loop(f, this->max_level);
	}
	template<class F>
	inline void loop(F f, int max_level) {
		int* levels = new int[dim];
		if(dim == 1) {
			for(int l = 0; l < max_level; l++) {
				levels[0] = l;
				f(l, levels);
			}
		} else {
			int* maxs = new int[dim];
			for(int l = 0; l < max_level; l++) {
				//printf("level: %d\n", l);
				for(int d = 0; d < dim; d++) {
					levels[d] = 0;
					maxs[d] = 0;
				}
				
				int norm_left = l;
				bool done = false;
				while(!done) { // loop l1=0..l, l2=0..l-l1, l3=0..l-l1-l2
					norm_left = l;
					for(int d = 0; d < dim-1; d++) {
						//printf(" norm left: %d, sub: %d ", norm_left, levels[d]);
						norm_left -= levels[d];
						//printf("norm left: %d\n", norm_left);
					}

					//if(norm_left < 0)
					//	break;
					levels[dim-1] = norm_left;
					int total = l;
					//int levelcount = 1;
					for(int d = 0; d < dim; d++) {
						//printf(" %d ", levels[d]);
						maxs[d] = total;
						total -= levels[d];
						//levelcount *= basis_set.count(levels[d]);
					}
					f(l, levels);
					//printf(" levelcount=%d ", levelcount);
					//count += levelcount;
					//fflush(stdout);
					//printf("\n max ");
					for(int d = 0; d < dim; d++) {
						//printf(" %d ", maxs[d]);
					}
					//fflush(stdout);
					//printf("\n");
					//int total = 0;//levels[0];
					for(int d = dim-2; d >= 0; d--) {
						//printf("level[%d]++ ", d);
						levels[d]++;
						if(levels[d]>maxs[d]) {
							levels[d] = 0;
							if(d == 0) done = true;
							//printf("level[%d]=0 ", d);
							//total += levels[d];
						} else {
							break;
						}
					}
					//printf("\n");
					//if(levels[dim-1] > 0) break;
					
				}
				
				/*for(int d = 0; d < dim-1; dim++) {
					while(norm_left > 0) {
						levels[d]
					}
				}*/
			}
		}
		delete[] levels;
	}

	template<class F>
		inline void loop_(F f, int max_level) {
			int* levels = new int[dim];
			if(dim == 1) {
				for(int l = 0; l < max_level; l++) {
					levels[0] = l;
					f(l, levels);
				}
			} else {
				int* maxs = new int[dim];
				for(int l = 0; l < max_level; l++) {
				//printf("level: %d\n", l);
					for(int d = 0; d < dim; d++) {
						levels[d] = 0;
						maxs[d] = l;
					}
					
					bool done = false;
					while(!done) { // loop l1=0..l, l2=0..l-l1, l3=0..l-l1-l2
						int norm = 0;
						double a = 0;
						for(int d = 0; d < dim; d++) {
							norm += levels[d];
							a += pow(4, levels[d]);
						}
						int n = max_level;//s- 1;
						bool condition = norm - 1./5 * log2(a) <= (n + dim -1) -1./5 * log2(pow(4, n)+4*dim- 4);
						//condition = norm == l;
						if(norm == l)
							f(l, levels);
					//fflush(stdout);
					//printf("\n");
					//int total = 0;//levels[0];
						for(int d = dim-1; d >= 0; d--) {
						//printf("level[%d]++ ", d);
							levels[d]++;
							if(levels[d]>maxs[d]) {
								levels[d] = 0;
								if(d == 0) done = true;
							//printf("level[%d]=0 ", d);
							//total += levels[d];
							} else {
								break;
							}
						}
					//printf("\n");
					//if(levels[dim-1] > 0) break;
						
					}
					
				/*for(int d = 0; d < dim-1; dim++) {
					while(norm_left > 0) {
						levels[d]
					}
				}*/
				}
			}
			delete[] levels;
		}
	
	/*double call_level(int max_level, double* point) {
		int *counts = new int[dim];
		int *indices = new int[dim];
		double value = 0;
		int value_index = 0;


		auto f = [&](int level, int* levels) {
			for(int d = 0; d < dim; d++) {
				counts[d] = this->basis_set.count(levels[d]);
				indices[d] = 0;
			}
			bool done = false;
			while(!done) {
				double weight = 1;
				for(int d = 0; d < dim; d++) {
					weight *= basis_set.operator()(levels[d], indices[d], point[d]);
					//printf("l=%d i=%d p=%f v=%f\n", levels[d], indices[d], point[d], basis_set.operator()(levels[d], indices[d], point[d]));
				}
				//printf("v = % 10.4f w = % 10.4f\n", values[value_index], weight);
				value += values[value_index] * weight; 
				value_index++;
				for(int d = dim-1; d >= 0; d--) {
					indices[d]++;
					if(indices[d]>=counts[d]) {
						indices[d] = 0;
						if(d == 0) done = true;
					} else {
						break;
					}
				}
			}
		};
		loop(f, max_level);
		return value;
	}*/
	/*template<class F>
	void eval_old(F function) {
		int *counts = new int[dim];
		int *indices = new int[dim];
		double *point = new double[dim];
		int value_index = 0;
		auto f = [&](int level, int* levels) {
			for(int d = 0; d < dim; d++) {
				counts[d] = this->basis_set.count(levels[d]);
				indices[d] = 0;
			}
			bool done = false;
			while(!done) {
				printf("\ni:");
				for(int d = 0; d < dim; d++) {
					printf(" %d", indices[d]);
				}
				for(int d = 0; d < dim; d++) {
					point[d] = this->basis_set.center(levels[d], indices[d]);
				}
				printf(" p:");
				for(int d = 0; d < dim; d++) {
					printf(" % 6.2f", point[d]);
				}
				for(int d = 0; d < dim; d++) {
					//printf("w = %basis_set.operator()(levels[d], indices[d], point[d]);
					printf("l=%d i=%d p=%f v=%f\n", levels[d], indices[d], point[d], basis_set.operator()(levels[d], indices[d], point[d]));
				}

				values[value_index] = function(dim, point);
				printf(" f(...) = % 10.4f", values[value_index]);
				if(level > 0) {
					double pvalue = this->call_level(level, point); 
					double surplus = values[value_index] - pvalue;
					printf(" surplus = % 10.4f pvalue = % 10.4f", surplus, pvalue);
					values[value_index] = surplus;
				} 
				value_index++;
				printf("\n");
				for(int d = dim-1; d >= 0; d--) {
					indices[d]++;
					if(indices[d]>=counts[d]) {
						indices[d] = 0;
						if(d == 0) done = true;
					} else {
						break;
					}
				}
			}
		};
		loop(f);
		printf("count=%d\n", count);
	}*/

};



class RegularGrid {
public:
	double *values;
	int count;
	int dim, max_level;
	BasisTriangular basis;
	RegularGrid(int dim, int max_level) : dim(dim), max_level(max_level) {
		count = pow(basis.count(max_level), dim);
		values = new double[count];
	}
	double operator()(double* point) {
		return call_level(max_level, point);
	}
	double call_level(int max_level, double* point) {
		int *counts = new int[dim];
		int *indices = new int[dim];
		int value_index = 0;
		double value = 0;
		for(int d = 0; d < dim; d++) {
			counts[d] = basis.count(max_level);
			indices[d] = 0;
		}
		bool done = false;
		while(!done) { // loop over indices
			double weight = 1;
			for(int d = 0; d < dim; d++) {
				weight *= basis(max_level, indices[d], point[d]);
			}
			value += values[value_index++] * weight; 
			for(int d = dim-1; d >= 0; d--) {
				indices[d]++;
				if(indices[d]>=counts[d]) {
					indices[d] = 0;
					if(d == 0) done = true;
				} else {
					break;
				}
			}
		}
		return value;
	}

	template<class F>
	void eval(F function) {
		int *counts = new int[dim];
		int *indices = new int[dim];
		double *point = new double[dim];
		int value_index = 0;
		for(int d = 0; d < dim; d++) {
			counts[d] = basis.count(max_level);
			indices[d] = 0;
		}
		bool done = false;
		while(!done) { // loop over indices
			//printf("f(");
			for(int d = 0; d < dim; d++) {
				point[d] = basis.center(max_level, indices[d]);
				//printf("% 5.3f, ", point[d]);
			}
			//printf(") = ");
			values[value_index++] = function(dim, point);
			//printf("%f", values[value_index-1]);
			//for(int d = 0; d < dim; d++)
			//	printf("% 3d ", indices[d]);
			//printf(" = \n");
			for(int d = dim-1; d >= 0; d--) {
				indices[d]++;
				if(indices[d]>=counts[d]) {
					indices[d] = 0;
					if(d == 0) done = true;
				} else {
					break;
				}
			}
		}
	}

};

