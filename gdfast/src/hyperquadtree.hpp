#include <iostream>
#include <tuple>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include "nlinear.hpp"
#include <assert.h>
#include <cmath>
#include <list>
#include <vector>
#include <iostream>
#include <stdexcept>

using namespace std;

/*template<typename... Elements>
ostream & operator << (ostream &stream, tuple<Elements...> &t)  {
	
	return stream << get<0>(t);
}*/

/*template<class Tuple, size_t N>
ostream & operator << (ostream &stream, tuple<Elements...>) {
	
	return stream << " .. ";
}

template<class Tuple, size_t N>
ostream & operator << (ostream &stream, tuple<Elements...>) {
	
	return stream << " .. ";
}*/

template<class T, size_t N>
struct print_tuple_helper {
	static void print(ostream &stream, T t) {
		print_tuple_helper<T, N-1>::print(stream, t);
		stream << get<N>(t) << "\t";
	}
};

template<class T>
struct print_tuple_helper<T, 0> {
	static void print(ostream &stream, T t) {
		stream << get<0>(t) << "\t";
	}
};

/*template<typename... Elements>
ostream & operator << (ostream &stream, tuple<Elements...> &t)  {
	typename tuple<Elements...> Tuple;
	print_tuple_helper<Tuple, tuple_size<Tuple>::value>(t);
	return stream << " ";
}
*/

/* returns 0 is v=-inf, otherwise return exp(x)  */
/*template<class T>
T save_exp(T v) {
	cout << v << " ";
	return ((v < 0) && isinf(v)) ? 0 : exp(v);
}*/

template<class T>
ostream & test_pr (ostream &stream, T &t)  {
	//typename tuple<Elements...> Tuple;
	print_tuple_helper<T, tuple_size<T>::value-1>::print(stream, t);
	//return stream << " ";
	return stream;
}

template<class T, size_t I>
struct do_tuple_helper {
	template<class F>
		static void do_tuple(T& t, size_t index, F f) {
			if(index == I)
				f(get<I>(t), index == tuple_size<T>::value-1);
			else
				do_tuple_helper<T, I-1>::do_tuple(t, index, f);
		}
};

template<class T>
struct do_tuple_helper<T, 0> {
	template<class F>
		static void do_tuple(T& t, size_t index, F f) {
			if(index == 0)
				f(get<0>(t), index == tuple_size<T>::value-1);
		}
};

template<class T, class F>
void do_tuple(T& t, size_t index, F f) {
	do_tuple_helper<T, tuple_size<T>::value-1>::do_tuple(t, index, f);
}

template<class T, size_t I>
struct do_tuple_pair_helper {
	template<class F>
		static void do_tuple_pair(const T& t1, const T& t2, size_t index, F f) {
			if(index == I)
				f(get<I>(t1), get<I>(t2), index == tuple_size<T>::value-1);
			else
				do_tuple_pair_helper<T, I-1>::do_tuple_pair(t1, t2, index, f);
		}
};

template<class T>
struct do_tuple_pair_helper<T, 0> {
	template<class F>
		static void do_tuple_pair(const T& t1, const T& t2, size_t index, F f) {
			if(index == 0)
				f(get<0>(t1), get<0>(t2), index == tuple_size<T>::value-1);
		}
};


template<class T, class F>
void do_tuple_pair(const T& t1, const T& t2, size_t index, F f) {
	do_tuple_pair_helper<T, tuple_size<T>::value-1>::do_tuple_pair(t1, t2, index, f);
}

template<typename... Elements>
ostream & operator << (ostream &stream, tuple<Elements...> t)  {
	for(int i = 0; i < (int)tuple_size<tuple<Elements...>>::value; i++) {
		do_tuple(t, i, [&](double v, bool last) {
			stream << v << (last ? "" : " ");
		});
	}
	return stream;
}

template<typename... Elements>
istream & operator >> (istream &stream, tuple<Elements...>& t)  {
	for(int i = 0; i < (int)tuple_size<tuple<Elements...>>::value; i++) {
		do_tuple(t, i, [&](double& v, bool last) {
			stream >> v;
		});
	}
	return stream;
}



template<int N>
class NTreeDim {
}; 

template<>
class NTreeDim<2> {
public:
	enum { nodes = 4 };
}; 

template<int N, typename T, typename... Elements>
class NTree;
	
template<int D, typename T>
T copy_coordinate(T p1, T p2) {
	T p = p1;
	get<D>(p) = get<D>(p2);
	return p;
}
template<int N, typename T, int Nparent, typename... Elements>
class NTreeNode {
	public:
	typedef NTree<Nparent, T, Elements...> _NTree;
	typedef typename _NTree::Point Point;
	typedef ::NLinear<N, Nparent, T> NLinear;
	Point p1, p2;
	double x1, x2;
	_NTree* tree;
	NTreeNode<N-1, T, Nparent, Elements...> left, right;
	enum { DIM = Nparent-N };

	template<typename... Tail>
	NTreeNode(_NTree *tree, typename _NTree::Point p1, typename _NTree::Point p2, double volume=1.) : tree(tree), left(tree, copy_coordinate<Nparent-N>(p1, p1), copy_coordinate<Nparent-N>(p2, p1), volume/2), right(tree, copy_coordinate<Nparent-N>(p1, p2), copy_coordinate<Nparent-N>(p2, p2), volume/2) {
		x1 = get<DIM>(p1);
		x2 = get<DIM>(p2);
		//updatepoints();
	}

	/*template<typename... Tail>
	NTreeNode(_NTree *tree, typename _NTree::Point p1, typename _NTree::Point p2, vector<bool>::iterator *split_iter) : tree(tree), left(tree, copy_coordinate<Nparent-N>(p1, p1), copy_coordinate<Nparent-N>(p2, p1), split_iter), right(tree, copy_coordinate<Nparent-N>(p1, p2), copy_coordinate<Nparent-N>(p2, p2), split_iter) {
		x1 = get<DIM>(p1);
		x2 = get<DIM>(p2);
		//updatepoints();
	}*/


	/*typedef vector<tuple<int,int,bool>>::iterator test_iter;
	NTreeNode(_NTree *tree, vector<Point> points, test_iter *connection_iter, int level) : tree(tree), p1(points[get<0>(**connection_iter)]), p2(points[get<1>(**connection_iter)]), x1(get<DIM>(points[get<0>(**connection_iter)])), x2(get<DIM>(points[get<1>(**connection_iter)])), left(tree, points, connection_iter, level+1), right(tree, points, connection_iter, level+1) {
		
		/* this would do the same as the x1(..) and x2(..) but it needs to be done before the iterator has moved beyond its current point
		int i1 = get<0>(*connection_iter);
		int i2 = get<1>(*connection_iter);
		x1 = get<DIM>(points[i1]);
		x2 = get<DIM>(points[i2]);
		*/
		/*int l = level;
		while(l--) cout << " ";
		cout << "x1, x2: " << x1 << " " << x2 << endl;
		cout << "p1, p2: " << p1 << " " << p2 << endl;
		//cout << "x1, x2: " << x1 << " " << x2 << endl;
	}*/

	void get_vertex_levels(int level) {
		left.get_vertex_levels(level);
		right.get_vertex_levels(level);
	}

	void check_maximum_level_difference(int level, int max_level_difference) {
		if((N == Nparent) && !haschildren()) {
			if(this->max_level_difference(level) > max_level_difference) {
				split();
			}
		} else {
			left.check_maximum_level_difference(level, max_level_difference);
			right.check_maximum_level_difference(level, max_level_difference);
		}
	}

	int max_level_difference(int level) {
		assert(!haschildren()); // this doesn't make sense to ask if we don't have children
		return max(left.max_level_difference(level), right.max_level_difference(level));
	}

	template<class C>
	void output_state(C& out_stream) {
		if(N == Nparent) {
			if(issplit())
				out_stream << true << " ";
			else
				out_stream << false << " ";
		}
		left.output_state(out_stream);
		right.output_state(out_stream);
	}

	void split(vector<bool>::iterator *split_iter) {
		if(N == Nparent) {
			bool splitme = **split_iter;
			//cout << "splitme? :" << splitme << endl;
			(*split_iter)++;
			if(splitme) {
				split();
				left.split(split_iter);
				right.split(split_iter);
			}
		} else {
			left.split(split_iter);
			right.split(split_iter);
		}
	}
	
	//T* subnodeValues[NTreeDim<N>::nodes];
	/*void printinfo() {
		cout << "from " << x1 << " to " << x2 << endl;
		cout << "\t";
		left.printinfo();
		cout << "left:" << endl;
		
	}*/
	template<typename F>
	void walkrange(F f) {
		left.walkrange(f);
		right.walkrange(f);
	}
	template<typename F>
	void visit(F f) {
		left.visit(f);
		right.visit(f);
	}
	template<typename F>
	void walk(F f) {
		left.walk(f);
		right.walk(f);
	}
	/*template<typename F>
	void eval(F f) {
		left.eval(f);
		right.eval(f);
	}*/
	T eval(Point p) {
		if(haschildren()) {
			T xmid = (x1+x2)/2;
			T x = get<DIM>(p);
			assert(x >= x1);
			assert(x <= x2);
			if(x < xmid)
				return left.eval(p);
			else
				return right.eval(p);
		} else {
			return nlinear().eval(p);
			//return 0;
		}
	}


	void split() {
		double xmid = (x1+x2)/2;
		Point p1;
		Point mid;
		Point p2;
		get<DIM>(p1) = x1;
		get<DIM>(mid) = xmid;
		get<DIM>(p2) = x2;
		left.split(p1, mid);
		right.split(mid, p2);
		//updatepoints();
	}
	void split(Point p1, Point p2) {
		double xmid = (x1+x2)/2;
		Point mid;
		
		mid = p2;
		get<DIM>(p1) = x1;
		get<DIM>(mid) = xmid;
		left.split(p1, mid);

		mid = p1;
		get<DIM>(mid) = xmid;
		get<DIM>(p2) = x2;
		right.split(mid, p2);
		//updatepoints();
	}
	T integrate(double &error) {
		if(!haschildren()) {
			return current_level_integrate(x1, x2);
		} else {
			if(hasnograndchildren() && (N==Nparent)) {
				double I1 = current_level_integrate(x1, x2);
				double error2 = 0;
				double I2 = left.integrate(error2) + right.integrate(error2);
				assert(error2 == 0);
				assert(error2 == 0);
				error += fabs(I1-I2);
				return I2;
			} else {
				return left.integrate(error) + right.integrate(error);
			}
		}
	}

	void geterrors(int maxlevel, int level=0) {
		if(haschildren() && (level <= maxlevel)) {
			if(hasnograndchildren() && (N==Nparent)) {
				// we are at a point where we have just 1 level below us
				Point p;
				double maxerror = get_maxerror(nlinear(), p);
				//cout << "maxerror = " << maxerror << " " << p << endl;
				double volume = this->volume();
				typename _NTree::ErrorTuple t(this, maxerror, maxerror * volume );
				tree->errornodes.push_back( t );
			}
			else {
				// pass on to lower levels
				left._geterrors(maxlevel, level);
				right._geterrors(maxlevel, level);
			}
		}
		
	}
	void _geterrors(int maxlevel, int level) {
		left._geterrors(maxlevel, level);
		right._geterrors(maxlevel, level);
	}

	template<class NLinear>
	double get_maxerror(NLinear current_level_nlinear, Point& p) {
		double maxerror = 0;
		double xmid = (x1+x2)/2;
		get<DIM>(p) = x1;
		maxerror = left.get_maxerror(current_level_nlinear, p);
		get<DIM>(p) = xmid;
		maxerror = max(maxerror, left.get_maxerror(current_level_nlinear, p));
		get<DIM>(p) = x2;
		maxerror = max(maxerror, right.get_maxerror(current_level_nlinear, p));
		return maxerror;
	}
	
	T volume() {
		return (x2-x1)*left.volume();
	}

	T error() {
		return 0;
	}

	template<typename... Coordinates>
	T current_level_integrate(T x1, T x2, Coordinates... coordinates) {
		//cout << "current_level_integrate: " << x1 << ", " << x2 << endl;
		return (nlinear().integrate(x1, x2, coordinates...));
	}
	NLinear nlinear() {
		auto l1 = left.nlinear();
		auto l2 = right.nlinear();
		return NLinear(l1, l2, x1, x2);
	}
	/*template<typename... Coordinates>
	T current_level_integrate(Coordinates... coordinates) {
		current_level_integrate(x1, x2, coordinates...);
	}*/

	void optimize() {
		double xmid = (x1+x2)/2;
		double volume = (xmid-x1);
		left.optimize(this, volume, x1, xmid);
		right.optimize(this, volume, xmid, x2);
	}
	
	template<typename... Coordinates>
	void optimize(NTreeNode<Nparent, T, Nparent, Elements...> *root, double volume, Coordinates... coordinates) {
		double xmid = (x1+x2)/2;
		volume *= (xmid-x1);
		left.optimize(root, volume, coordinates..., x1, xmid);
		right.optimize(root, volume, coordinates..., xmid, x2);
	}
	template<class F>
	void visit_leafs(F f) {
		left.visit_leafs(f);
		right.visit_leafs(f);
	}
	
	bool issplit() {
		return haschildren(); // is has children, means split
	}
	bool haschildren() {
		return left.haschildren(); // if left has, right has it as well
	}
	bool hasnochildren() {
		return left.hasnochildren() && right.hasnochildren();
	}
	bool hasnograndchildren() {
		return left.hasnograndchildren() && right.hasnograndchildren();
	}
	bool allhavegrandchildren() {
		return left.allhavegrandchildren() && right.allhavegrandchildren();
	}
private:
	void updatepoints() {
		struct F {
			_NTree *tree;
			F(_NTree *tree) : tree(tree) {}
			void operator() (NTreeNode<1, T, Nparent, Elements...>  *node, Elements... elements) {
				tree->add_point(elements...);
			}
		};
		F f(tree);
		visit(f);
	}
};

template<typename T, int Nparent, typename... Elements>
class NTreeNode<1, T, Nparent, Elements...> {
	public:
	typedef NTree<Nparent, T, Elements...> _NTree;
	typedef typename _NTree::Point Point;
	_NTree* tree;
	double x1, x2;
	enum { DIM = Nparent-1 };
	const typename _NTree::Point *p1, *p2;
	T v1, v2;
	typedef NTreeNode<Nparent, T, Nparent, Elements...> NSubTreeNode;
	NSubTreeNode *subleft, *subright;
	typedef ::NLinear<1, Nparent, T> NLinear;
	double point_volume;
	enum { subvertices = 1 };
		

	NTreeNode(_NTree *tree, typename _NTree::Point p1, typename _NTree::Point p2, double point_volume) : tree(tree), subleft(NULL), subright(NULL), point_volume(point_volume) {
		this->p1 = tree->add_point(p1); 
		this->p2 = tree->add_point(p2); 
		x1 = get<DIM>(p1);
		x2 = get<DIM>(p2);
	}
	/*NTreeNode(_NTree *tree, typename _NTree::Point p1, typename _NTree::Point p2, vector<bool>::iterator *split_iter) : tree(tree), subleft(NULL), subright(NULL) {
		this->p1 = tree->add_point(p1); 
		this->p2 = tree->add_point(p2); 
		x1 = get<DIM>(p1);
		x2 = get<DIM>(p2);
	}*/
	NTreeNode(_NTree *tree, vector<Point> points, vector<tuple<int,int,bool>>::iterator* connection_iter, int level, double point_volume) : tree(tree), subleft(NULL), subright(NULL), point_volume(point_volume) {
		int i1 = get<0>(**connection_iter);
		int i2 = get<1>(**connection_iter);
		bool last = get<2>(**connection_iter);
		p1 = tree->add_point(points[i1]);
		p2 = tree->add_point(points[i2]);
		x1 = get<DIM>(*p1);
		x2 = get<DIM>(*p2);
		int l = level;
		while(l--) cout << " ";
		cout << "x1, x2: " << x1 << " " << x2 << endl;
		(*connection_iter)++;
		if(!last) {
			subleft = new NSubTreeNode(tree, points, connection_iter, level+1, volume/2);
			subright = new NSubTreeNode(tree, points, connection_iter, level+1, volume/2);
		}

	}
	
	~NTreeNode() {
		if(subleft)
			delete subleft;
		if(subright)
			delete subright;
	}

	void get_vertex_levels(int level) {
		tree->has_minimum_level(p1, level);
		tree->has_minimum_level(p2, level);
		if(haschildren()) {
			subleft->get_vertex_levels(level+1);
			subright->get_vertex_levels(level+1);
		}
	}

	int max_level_difference(int level) {
		assert(!haschildren()); // this doesn't make sense to ask if we don't have children
		return max(abs(tree->levels[*p1]-level), abs(tree->levels[*p2]-level));
	}

	void check_maximum_level_difference(int level, int max_level_difference) {
		assert(haschildren()); // we shouldn't end up here if we don't
		subleft->check_maximum_level_difference(level+1, max_level_difference);
		subright->check_maximum_level_difference(level+1, max_level_difference);
	}

	template<class C>
	void output_state(C& out_stream) {
		if(haschildren()) {
			subleft->output_state(out_stream);
			subright->output_state(out_stream);
		} else {
		}
	}

	T eval(Point p) {
		assert(haschildren());
		T xmid = (x1+x2)/2;
		T x = get<DIM>(p);
		assert(x >= x1);
		assert(x <= x2);
		if(x < xmid)
			return subleft->eval(p);
		else
			return subright->eval(p);
	}

	template<typename F, typename... Tail>
	void walkrange(F f, Tail... tail) {
		f(tail..., *p1);
		f(tail..., *p2);
	}
	template<typename... Tail>
	void split(Point p1, Point p2) {
		double xmid = (x1+x2)/2;
		Point mid;

		get<DIM>(p1) = x1;
		get<DIM>(p2) = x2;

		//cout << "x1: " << x1 <<", " << xmid << ", " << x2  << endl;

		mid = p2;
		get<DIM>(mid) = xmid;
		if(subleft)
			subleft->split();
		else
			subleft = new NSubTreeNode(tree, p1, mid, point_volume/2);
		mid = p1;
		get<DIM>(mid) = xmid;
		if(subright)
			subright->split();
		else
			subright = new NSubTreeNode(tree, mid, p2, point_volume/2);
	}
	void split(vector<bool>::iterator *split_iter) {
		assert(haschildren());
		subleft->split(split_iter);
		subright->split(split_iter);
	}

	template<typename F, typename... Tail>
	void visit(F f, Tail... tail) {
		//f(this, tail..., x1);
		//f(this, tail..., x2);
	}
	template<typename F>
	void walk(F f) {
		if(haschildren()) {
			f(p1, p2, false);
			subleft->walk(f);
			subright->walk(f);
		}
		else {
			f(p1, p2, true);
		}
	}
	/*template<typename F, typename... Tail>
	void eval(F f, Tail... tail) {
		v1 = f(tail..., x1);
		v2 = f(tail..., x2);
		typename _NTree::Point p1(tail..., x1);
		typename _NTree::Point p2(tail..., x2);
		tree->values[p1] = v1;
		tree->values[p2] = v2;
	}*/
	
	NLinear nlinear() {
		return nlinear(x1, x2);
	}
	NLinear nlinear(T x1, T x2) {
		//cout << "nlinear<1>: " << x1 << ", " << x2 << "| " << tree->values[*p1] << ", " << tree->values[*p2] << endl;
		tree->assert_has_value(p1);
		tree->assert_has_value(p2);
		return NLinear(tree->values[*p1], tree->values[*p2], x1, x2);
	}
	T integrate(double &error) {
		assert(haschildren());
		return subleft->integrate(error) + subright->integrate(error);
	}
	/*double subintegrate() {
		assert(haschildren);
		return subleft->integrate() + subright->integrate();
	}*/
	bool haschildren() {
		return (subleft != NULL) && (subright != NULL);
	}
	bool hasnochildren() {
		return (subleft == NULL) && (subright == NULL);
	}
	bool allhavegrandchildren() {
		return (subleft != NULL) && (subright != NULL) && subleft->haschildren() && subright->haschildren(); 
	}
	bool hasnograndchildren() {
		assert(haschildren());
		return subleft->hasnochildren() && subright->hasnochildren();
	}
	void _geterrors(int maxlevel, int level) {
		subleft->geterrors(maxlevel, level+1);
		subright->geterrors(maxlevel, level+1);
	}

	template<class NLinear>
	double get_maxerror(NLinear current_level_nlinear, Point& p) {
		double maxerror = 0;
		double xmid = (x1+x2) / 2;
		get<DIM>(p) = x1;
		maxerror = abs(exp(current_level_nlinear.eval(p)) - exp(subleft->eval(p)));
		//cout << "  " << p << " " << maxerror  << " " << current_level_nlinear.eval(p) << " " << subleft->eval(p) << endl; 
		get<DIM>(p) = xmid;
		maxerror = max(maxerror, abs(exp(current_level_nlinear.eval(p)) - exp(subleft->eval(p))));
		//cout << "  " << p << " " << maxerror << " " << current_level_nlinear.eval(p) << " " << subleft->eval(p) << endl; 
		get<DIM>(p) = x2;
		maxerror = max(maxerror, abs(exp(current_level_nlinear.eval(p)) - exp(subright->eval(p))));
		//cout << "  " << p << " " << maxerror << " " << current_level_nlinear.eval(p) << " " << subright->eval(p) << endl; 
		//cout << "  " << endl; 
		return maxerror;
	}

	T volume() {
		return (x2-x1);
	}

	template<typename... Coordinates>
	void optimize(NTreeNode<Nparent, T, Nparent, Elements...> *root, double volume, Coordinates... coordinates) {
		assert(haschildren());
		double xmid = (x1+x2)/2;
 		volume *= (xmid-x1);
		double error;
		if(subleft->hasnochildren()) {
			//cout << "left" << endl;
			//cout << "from " << I1 << " to " << I2 << endl;
			//cout << "extra" << x1 << ", " << xmid << endl;
			double I1 = root->current_level_integrate(coordinates..., x1, xmid);
			double I2 = subleft->integrate(error);
			typename _NTree::ErrorTuple t(subleft, fabs(I1-I2) * pow(volume,0), fabs((I1-I2)/I2) );
			tree->errornodes.push_back( t );
			
			/*if( (fabs((I1-I2)) > abserror) || (fabs((I1-I2)/I2) > relerror)) {
				subleft->split();
			}*/
		} else {
			subleft->optimize();
		}
		if(subright->hasnochildren()) {
			//cout << "right" << endl;
			double I1 = root->current_level_integrate(coordinates..., xmid, x2);
			double I2 = subright->integrate(error);
			//cout << "from " << I1 << " to " << I2 << endl;
			typename _NTree::ErrorTuple t(subright, fabs(I1-I2) * pow(volume,0), fabs((I1-I2)/I2) );
			tree->errornodes.push_back( t );
			/*if( (fabs((I1-I2)) > abserror) || (fabs((I1-I2)/I2) > relerror)) {
				subright->split();
			}*/
		} else {
			subright->optimize();
		}
	}
	template<class F>
	void visit_leafs(F f) {
		if(haschildren()) {
			subleft->visit_leafs(f);
			subright->visit_leafs(f);
		} else {
			f(*this);
		}
	}
	/*void printinfo() {
		cout << "from " << x1 << " to " << x2 << endl;
	}*/
};


/*template<typename... Coordinates>
class Point : public tuple<Coordinates...> {

};*/

template<class Point>
struct PointCompare {
	bool operator()(const Point& p1, const Point& p2) const {
		bool precedes = true;
		for(int i = 0; i < tuple_size<Point>::value; i++) {
			do_tuple_pair(p1, p2, i, [&](const double v1, const double v2, bool last) {
				precedes = precedes & (v1 < v2); 
			});
		}
		return precedes;
	}
};

template<int N, typename T, typename... Elements>
class NTree {
public:
	typedef tuple<Elements...> Point;
	typedef NTreeNode<N, T, N, Elements...> RootNode;
	typedef NTreeNode<1, T, N, Elements...> LeafNode;
	typedef tuple<RootNode*, double, double> ErrorTuple;
	//map<Point, T, PointCompare<Point> > values;
	Point p1, p2;
	map<Point, int> indices;
	set<Point> points;
	map<Point, T> values;
	int point_index;
	
	map<Point, int> levels;
	map<Point, T> volumes;
	map<int, Point> pointmap;
	list< ErrorTuple > errornodes;
	int new_points;

	RootNode rootNode;

	
	

	NTree(Point p1, Point p2) : p1(p1), p2(p2), indices(), points(), values(), point_index(0), new_points(0), rootNode(this, p1, p2) {
		;
		//add_point(p1);
		//add_point(p2);
	}

	/*NTree(vector<Point> points, vector<tuple<int,int,bool>> connectionlist) : indices(), points(), values(), point_index(0), rootNode(this, points, &connectionlist.begin(), 0) {
		;
		//add_point(p1);
		//add_point(p2);
	}*/
	
	NTree(Point p1, Point p2, vector<bool>::iterator *split_iter) : p1(p1), p2(p2), indices(), points(), values(), point_index(0), rootNode(this, p1, p2, split_iter) {
		;
		//add_point(p1);
		//add_point(p2);
	}	
	~NTree() {
	}
	
	
	void output(string name) {
		_accumulate_volumes();
		ofstream f;
		f.open((name + ".state").c_str());
		/*map<Point, int> indices;
		int index = 0;
		//for_each(points.begin(), points.end(), [&](Point p) {
		for(int i = 0; i < pointmap.size(); i++) {
			Point p = pointmap[i];
			//indices[p] = index;
			//index++;
			for(int i = 0; i < N; i++) {
				//f << (get<0>(p)) << " ";
				//print_tuple(
				//f << p;
				//test_pr(f, p);
				do_tuple<Point>(p, i, [&f](const T v, bool last) {
					f << v << (last ? "" : " ");
				});
			}
			f << endl;
		}
		f.close();

		f.open(filename_connect);


		rootNode.walk([&](const Point* p1, const Point *p2, bool end) {
			//f << "(" << (*p1) << ") , (" << (*p2) << ") ";
			f << this->indices[*p1] << " " << this->indices[*p2] << " " << end << endl; 
			//f << endl;
			//cout << indices[*p1] << " " << indices[*p2] << " " << end << endl;
		});*/
		f << p1 << " "  << p2 << endl;
		rootNode.output_state(f);
		f.close();

		ofstream f_known;
		ofstream f_unknown;
		f_known.open((name + ".known").c_str());
		f_unknown.open((name + ".unknown").c_str());
		f_known.precision(30);
		f_unknown.precision(30);
		int count_known = 0;
		int count_unknown = 0;
		for_each(points.begin(), points.end(), [&](Point p) {
			//cout << this->values.count(p) << endl;
			if(this->values.count(p) > 0) {
				f_known << scientific << p << " " << this->volumes[p] << " " << this->values[p] << endl;
				count_known++;
			} else {
				f_unknown << scientific << p << " " << this->volumes[p]  << endl;
				count_unknown++;
			}
		});
		printf("known: %d unknown: %d\n", count_known, count_unknown);
	}

	void readpoints(string name) {
		ifstream f;
		f.open(name);
		cout << "reading from filename: " << name << endl;
		if(!f.is_open())
			throw runtime_error((string("file doesn't exists: ") + name));

		while(!f.eof()) {
			T value, volume;
			Point p;
			f >> p >> volume >> value;
			if(f.good()) {
				//cout << p << " = " << value << " " << f.good() << endl;
				values[p] = value;
			}
		}

	}
	
	/*template< class tuple<T> >
	void _output_pt(ofstream& f, ) {
		f << (get<i>(p)) << " ";
	}*/
	
	
	/*const Point* add_point(Elements... elements) {
		Point p(elements...);
		if(points.
		const Point* pstored = &(*(get<0>(points.insert(p))));
		cout << "point_index = " << point_index << " check: ";
		indices[*pstored] = point_index++;
		cout << "point_index = " << indices[*pstored] << endl;
		assert(pstored != NULL);
		return pstored;
	}*/
	const Point* add_point(Point& p) {
		const Point* pstored;
		if(points.count(p) == 0) {
			auto it = points.insert(p);
			pstored = &(*(it.first));
			//cout << "point_index = " << point_index << " check: ";
			pointmap[point_index] = *pstored;
			indices[*pstored] = point_index++;
			new_points++;
			//cout << "point_index = " << indices[*pstored] << " " << p << endl;
		} else {
			auto it = points.insert(p);
			pstored = &(*(it.first));
		}
		assert(pstored != NULL);
		return pstored;
	}
	template<int Da=0, int Db=1>
	void output_ps(const char *filename, double scalex=1., double scaley=1., double r=1., double x0=50, double y0=50, double border=10) {
		ofstream myfile;
		myfile.open(filename);
		//myfile << "!PS" << endl;
		T x1 = get<Da>(p1); 
		T x2 = get<Da>(p2); 
		T y1 = get<Db>(p1); 
		T y2 = get<Db>(p2); 
		
		myfile << "%!PS-Adobe-3.0" << endl;
		myfile << "%%BoundingBox: " << (x0-border) << " " << (y0-border) << " " << (x0+scalex+border) << " " << (y0+scaley+border) << endl;
		for_each(points.begin(), points.end(), [&](Point p) {
			double x = get<Da>(p);
			double y = get<Db>(p);
			x = (x-x1)/(x2-x1);
			y = (y-y1)/(y2-y1);
			myfile << (x0+x*scalex) << " "  << (y0+y*scaley) << " " << r << " 0 360 arc closepath fill" << endl;
		});
		myfile.close();
	}
	template<int Da=0, int Db=1>
	void output_gpl(const char *filename, double scalex=1., double scaley=1., double r=1.) {
		ofstream myfile;
		myfile.open(filename);
		myfile << "# GNUPLOT format: x y value" << endl;
		for_each(values.begin(), values.end(), [&](tuple<Point, T> t) {
			Point p = get<0>(t);
			T value = get<1>(t);
			myfile << (get<Da>(p)*scalex) << " "  << (get<Db>(p)*scaley) << " " << value << endl;
		});
		myfile.close();
	}
	template<typename F>
	void eval(F f) {
		for_each(points.begin(), points.end(), [&](Point p) {
			T value = exp(f(p));
			this->values[p] = value;
		});
	}

	void assert_has_value(const Point *p) {
		if(values.count(*p) == 0) {
			stringstream ss;
			ss << "missing value for point: " << (*p);
			throw runtime_error(ss.str());
		}
	}
	
	T integrate(double &error) {
		return rootNode.integrate(error);
	}

	void optimize(double absfraction=0.2, double relfraction=0.05, int Nabsmin=1, int Nrelmin=0, int maxlevel=8, int max_level_difference=0) {
		errornodes.clear();
		//rootNode.optimize();
		rootNode.geterrors(maxlevel);
		errornodes.sort([](ErrorTuple& first, ErrorTuple& second) { return get<1>(first) > get<1>(second); } );
		//cout << "size: " << errornodes.size() << endl;
		int Nmax = max(Nabsmin, (int)(errornodes.size() *  absfraction));
		auto it = errornodes.begin();
		new_points = 0;
		while((new_points < Nmax) && (it != errornodes.end())) {
			(get<0>(*it))->split();
			cout << "new points (1) = " << new_points << " (" << Nmax << ")" << endl;
			it++;
			errornodes.pop_front();
			if(max_level_difference > 0) {
				rootNode.get_vertex_levels(0);
				rootNode.check_maximum_level_difference(0, max_level_difference);
			}
		}
		errornodes.sort([](ErrorTuple& first, ErrorTuple& second) { return get<2>(first) > get<2>(second); } );
		Nmax = max(Nrelmin, (int)(errornodes.size() *  relfraction));
		it = errornodes.begin();
		new_points = 0;
		while((new_points < Nmax) && (it != errornodes.end())) {
			(get<0>(*it))->split();
			cout << "new points (2) = " << new_points << " (" << Nmax << ")" << endl;
			it++;
			errornodes.pop_front();
			if(max_level_difference > 0) {
				rootNode.get_vertex_levels(0);
				rootNode.check_maximum_level_difference(0, max_level_difference);
			}
		}
		//cout << "size: " << errornodes.size() << " left" << endl;
		/*for(auto it = errornodes.begin(); it != errornodes.end(); it++) {
			cout << "errors: " << get<1>(*it) << endl;
		}*/

	}

	void limit_level_difference(int max_level_difference) {
		rootNode.get_vertex_levels(0);
		rootNode.check_maximum_level_difference(0, max_level_difference);
	}

	void has_minimum_level(const Point *p1, int level) {
		levels[*p1] = max(level, levels[*p1]);
		
	}
	template<typename... C>
	T eval(C... coordinates) {
		//rootNode.eval(coordinates...);
		return exp(rootNode.eval(Point(coordinates...)));;
	}

	T error() {
		return rootNode.error();
	}
protected:
	void _accumulate_volumes() {
		volumes.clear();
		rootNode.visit_leafs([&](LeafNode& leafnode) {
			// half of the volume associated to one vertex.
			this->volumes[*leafnode.p1] += leafnode.point_volume/2;
			this->volumes[*leafnode.p2] += leafnode.point_volume/2;
		});
	}

	//T* subnodeValues[NTreeDim<N>::nodes];
};