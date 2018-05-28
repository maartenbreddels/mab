#include "cpptools.hpp"
#include <cmath>
#include <stdio.h>
#include <tuple>

double mysum() {
	return 0;
}

template<typename Head, typename ...Tail>
double mysum(Head head, Tail... tail) {
	return head + mysum(tail...);
}

template<int N>
struct lagrange_basis {
	//typedef monomial<N-1> subtype;
	typedef typename genseq<N-1>::type seq;
	enum { size = N };
	template<typename... Tail>
	//double operator()(double x, double a, Tail... tail) {
	double operator()(double x, Tail... tail) {
		//int N = sizeof...(Tail);
		//return a * lj<sizeof...(Tail), N>(x) + this->operator()(x, tail...);
		//dot(a, tail..., 
		std::tuple<Tail...> t(tail...);
		return test(x, t, typename genseq<N-1>::type());
		//return 0;
	}
	template<int I>
	double unit_vector(double x) {
		//int N = sizeof...(Tail);
		//return a * lj<sizeof...(Tail), N>(x) + this->operator()(x, tail...);
		//dot(a, tail..., 
		//std::tuple<Tail...> t(tail...);
		return I;
		//return  lj<I, N>(x);
		//return 0;
	}
	template<typename T, int ...S>
	double test(double x, T t, ::seq<S...>) {
		return mysum(std::get<S>(t) * lj<S, N>(x)...);
		return 0;
	}
	double operator()(double x) {
		return 0;
	}
	template<int J, int K>
	double lj(double x) {
		printf("x = %.2f J=%d K=%d\n", x, J, K);
		typename genseq<K-1>::type s;
		return ljm<J, K>(x, s);
	}
	template<int J, int K, int M, int ...S>
	double ljm(double x, ::seq<M, S...> s) {
		//typename genseq<M-1>::type snext;
		double xm = (1.*M)/(K-1.);
		double xj = (1.*J)/(K-1.);
		::seq<S...> snext;
		if(J == M) {
			return ljm<J,K>(x, snext);
		} else {
			double y = (x-xm) / (xj-xm);
			//int j = J;
			//int k = K;
			//int m = M;
			printf("x = %.2f J=%d K=%d M=%d xm=%.2f xj=%.2f y=%.2f\n", x, J, K, M, xm, xj, y);
			return y * ljm<J,K>(x, snext);
		}
	}
	template<int J, int K>
	double ljm(double x, ::seq<> s) {
		return 1;
	}
};


template<int N>
struct monomial_basis {
	//typedef monomial<N-1> subtype;
	typedef typename genseq<N-1>::type seq;
	enum { size = N };
	template<typename... Tail>
	double operator()(double x, double a, Tail... tail) {
		int N = sizeof...(Tail);
		return a * pow(x, N) + this->operator()(x, tail...);
	}
	double operator()(double x) {
		return 0;
	}
};

template<typename B>
struct vector {
	typedef B basis;
	double a[B::size];
	template<typename... Args>
	vector(Args... args) : a{args...} {
	}
	
	double operator()(double x) {
		typename basis::seq s;
		return call(x, s);
	}
	template<int ...S>
	double call(double x, seq<S...>) {
		basis b;
		//double t = b.template unit_vector<1>(0.3);
		//return  b(x, this->x[S]...);
		return mysum((a[S] * (b.template unit_vector<S>(x)))...);
		//return mysum(a[S]...);
	}
	void print() {
	}
};


