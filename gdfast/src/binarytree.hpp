

namespace mab {


class BinaryNode {
public:
	BinaryNode(BinaryTree* tree) : tree(tree) {
		
	}
private:
	int get_dimension() { return tree->get_dimension(); }
	int max_child_nodes() { return tree->getDimension(); }
	BinaryTree* tree;
	BinaryNode childNodes;
}

class BinaryTree {
public:
	BinaryTree(int dimension);
	void optimize();
	template<class F>
	void eval(F function) {
		int value_index = 0;
		auto f = [&](int level, int* levels, int* indices, double* point) {
			values[value_index] = function(dim, point);
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
	
	int get_dimension() { return dimension; }
private:
	int dimension;
};


}