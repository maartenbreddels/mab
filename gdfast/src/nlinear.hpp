#include <fstream>
#include <cmath>
#include <tuple>
#include <iostream>

using namespace std;




template<int N, int Nparent, typename T>
class NLinear;

template<int Nparent, typename T>
class NLinear<1, Nparent, T>;

template<int N, int Nparent, typename T>
class NLinear {
public:
	typedef NLinear<N-1, Nparent, T> Sublinear;
	enum { DIM = Nparent-N };
	Sublinear sublinear1, sublinear2 ;
	T x0, x1;
	NLinear(Sublinear& sublinear1, Sublinear& sublinear2, T x0=0.0, T x1=1.0) : sublinear1(sublinear1), sublinear2(sublinear2), x0(x0), x1(x1) {
	}
	
	template<typename... Coordinates>
	T operator()(T x, Coordinates... coordinates) {
		return eval(x, coordinates...);
	}
	template<typename... Coordinates>
	T eval(T x, Coordinates... coordinates)
	{
		T t0 = sublinear1.eval(coordinates...);
		T t1 = sublinear2.eval(coordinates...);
		//return (isinf(t0) && isinf(t1) && (t0 < 0) && (t1 < 0)) ? log(0) : (x1-x)/(x1-x0) * t0 + (x-x0)/(x1-x0) * t1; 
		return (x1-x)/(x1-x0) * t0 + (x-x0)/(x1-x0) * t1; 
		//return (x1-x)/(x1-x0) * t0 + (x-x0)/(x1-x0) * t1; 
	}
	
	template<class Tuple>
	T operator()(Tuple t) {
		return eval(t);
	}
	template<class Tuple>
	T eval(Tuple t) {
		T x = get<DIM>(t);
		T t1 = sublinear1.eval(t);
		T t2 = sublinear2.eval(t);
		//return (isinf(t1) && isinf(t2) && (t1 < 0) && (t2 < 0)) ? log(0) : (x1-x)/(x1-x0) * t1 + (x-x0)/(x1-x0) * t2; 
		return (x1-x)/(x1-x0) * t1 + (x-x0)/(x1-x0) * t2;
	}
	/*template<class Tuple>
	T eval(Tuple )
	{
		T t1 = sublinear1.eval(coordinates...);
		T t2 = sublinear2.eval(coordinates...);
		return (x-x1)/(x2-x1) * t1 + (x2-x)/(x2-x1) * t2; 
	}*/
	
	T integrate() {
		return integrate(x0, x1);
	}
	template<typename... Coordinates>
	T integrate(T a, T b, Coordinates... coordinates) {
		return integrate2(a, b, x0, x1, coordinates...);
	}
	template<typename... Coordinates>
	T integrate2(T a, T b, T x0, T x1, Coordinates... coordinates) {
		double t1 = (sublinear1.integrate(coordinates...));
		double t2 = (sublinear2.integrate(coordinates...));
		//double f1 = (-(a*t2+b*t2-t2*x0+t1*x1)/(x0-x1)); 
		/*double in = (exp((b*t1+a*t2)/(x0-x1)+f1)-exp((a*t1+b*t2)/(x0-x1)+f1))*(x0-x1)/(t1-t2);
		if(t1 == t2)
			return -a * exp(t1) + b * exp(t2);
		else
			return in;*/
		return -(a-b)*(a*(t1-t2) + b*(t1-t2) + 2*t2*x0 - 2*t1*x1)/(2*(x0-x1));
	}
};
const double epsilon = 1e-10;

template<int Nparent, typename T>
class NLinear<1, Nparent, T> {
public:
	NLinear(T t1, T t2, T x0=0.0, T x1=1.0) : t1(t1), t2(t2), x0(x0), x1(x1) {
	}
	T eval(T x) {
		return (x1-x)/(x1-x0) * t1 + (x-x0)/(x1-x0) * t2;
		//return (x1-x)/(x1-x0) * log(t1+epsilon) + (x-x0)/(x1-x0) * log(t2+epsilon);
		//cout << "[" << t1  << " " << t2 << ", " << log(t1 < 0 ? 0 : t1) << ", " << log(t2 < 0 ? 0 : t2) << "]";
		//return (x1-x)/(x1-x0) * log(t1 < 0 ? 0 : t1) + (x-x0)/(x1-x0) * log(t2 < 0 ? 0 : t2);
		// if both endpoints are 0, return -inf
		//cout << "[" << ((t1 == 0) && (t2 == 0) ? log(0) : (x1-x)/(x1-x0) * log(t1) + (x-x0)/(x1-x0) * log(t2)) << " " << t1 << " " << t2 << "]";
		//return (t1 == 0) && (t2 == 0) ? log(0) : (x1-x)/(x1-x0) * log(t1) + (x-x0)/(x1-x0) * log(t2);
	}
	template<class Tuple>
	T eval(Tuple t) {
		T x = get<Nparent-1>(t);
		return eval(x);
	}
	T integrate() {
		//return (t1+t2)/2;
		return integrate(x0, x1);
	}
	T integrate(T a, T b) {
		//return (a-b)*(t1*(a+b-2)-t2*(a+b))/2;
		return integrate(a, b, x0, x1);
	}
	/*T integrate(T a, T b, T x0, T x1) {
		double f1 = (-(a*t2+b*t2-t2*x0+t1*x1)/(x0-x1)); 
		double in = (exp((b*t1+a*t2)/(x0-x1)+f1)-exp((a*t1+b*t2)/(x0-x1)+f1))*(x0-x1)/(t1-t2);
		//cout << f1 << " " << ((b*t1+a*t2)/(x0-x1)+f1) << " " << ((a*t1+b*t2)/(x0-x1)+f1) << " " <<  (x0-x1)/(t1-t2) << " | ";
		//cout << t1 << " " << t2 << " " << in << endl;
		if(t1 == t2)
			return -a * exp(t1) + b * exp(t2);
		else
			return in;
			//return -(a-b)*(a*(t1-t2) + b*(t1-t2) + 2*t2*x0 - 2*t1*x1)/(2*(x0-x1));
	}*/
	T integrate(T a, T b, T x0, T x1) {
		double t1 = (this->t1); // exp 
		double t2 = (this->t2);
		//return exp(-(a*t2+b*t2-t2*x0+t1*x1)/(x0-x1))*(exp((b*t1+a*t2)/(x0-x1))-exp((a*t1+b*t2)/(x0-x1)))*(x0-x1)/(t1-t2);
		return -(a-b)*(a*(t1-t2) + b*(t1-t2) + 2*t2*x0 - 2*t1*x1)/(2*(x0-x1));
	}
	T t1, t2, x0, x1;
};
