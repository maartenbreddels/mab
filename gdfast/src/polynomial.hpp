#include <fstream>
#include <cmath>
#include <tuple>
#include <iostream>

using namespace std;


template<size_t N, class T=double>
struct HermiteBasis {
	typedef HermiteBasis<N, T> type;
	typedef HermiteBasis<N-1, T> next_type;
	next_type next;
	T operator()(T x) { return x * next(x) - (N-1) * next.next(x); }
	//T dfdx(T x) { return Degree == 0 ? 0 : pow(x, Degree-1)/Degree; }
};

/* sentinel */
template<class T>
struct HermiteBasis<1, T> {
	typedef HermiteBasis<1, T> type;
	typedef HermiteBasis<0, T> next_type;
	next_type next;
	T operator()(T x) { return x; }
};

/* sentinel */
template<class T>
struct HermiteBasis<0, T> {
	T operator()(T) { return 1; }
};

template<class T>
struct HermiteBasis<-1, T> {
	T operator()(T) { return 0; }
};












template<int K, int J, int F=K, class T=double>
struct LagrangePolynomialKJF {
	typedef LagrangePolynomialKJF<K, J, F-1> next_type;
	next_type next;
	T operator()(T x) { return next(x) * (x-T(F)/K)/(T(J)/K-T(F)/K); }
	//T dfdx(T x) { return next.dfdx(x) * (1)/(T(J)/K-T(F)/K); }
	T dfdx(T x) { return next(x) * (1)/(T(J)/K-T(F)/K) + next.dfdx(x) * (x-T(F)/K)/(T(J)/K-T(F)/K); }
	/*cout << "nup K=" << K << " J=" << J << " F=" << F << " " << x << " " <<  ((x-T(F)/K)/(T(J)/K-T(F)/K))  << " ---- " << (T(J)/K-T(F)/K) << endl;*/
};

/* special for J==F */
template<int K, int J, class T>
struct LagrangePolynomialKJF<K,J,J,T> {
	typedef LagrangePolynomialKJF<K, J, J-1> next_type;
	next_type next;
	T operator()(T  x) { return next(x) * 1; }
	T dfdx(T x) { return next.dfdx(x) * 1; }
};

/* sentinel */
template<int K, int J, class T>
struct LagrangePolynomialKJF<K,J,-1,T> {
	T operator()(T) { return 1; }
	T dfdx(T) { return 0; }
};

template<int K, int J=K, class T=double>
struct LagrangeBasis {
	enum { degree = K };
	typedef LagrangePolynomialKJF<K,J> l_j_type;
	typedef LagrangeBasis<K, J, T> type;
	typedef LagrangeBasis<K,J-1,T> next_type;
	typedef LagrangeBasis<K-1,K-1> derivative_type;
	next_type next;
	l_j_type l_j;	
	T operator()(T x) { return l_j(x); }
	T dfdx(T x) { return l_j.dfdx(x); }
};

/* sentinel */
template<int K, class T>
struct LagrangeBasis<K, -1, T> {
	T operator()(T) { return 0; }
};



template<size_t Degree, class T=double>
struct MonomialBasis {
	typedef MonomialBasis<Degree, T> type;
	typedef MonomialBasis<Degree-1, T> next_type;
	typedef MonomialBasis<Degree-1, T> derivative_type;
	T operator()(T x) { return pow(x, Degree); }
	T dfdx(T x) { return Degree == 0 ? 0 : pow(x, Degree-1)/Degree; }
};

/* sentinel */
template<class T>
struct MonomialBasis<-1, T> {
	T operator()(T) { return 0; }
};

template<size_t Dim, size_t Degree, class Basis=MonomialBasis<Degree, double>, class RootBasis=Basis, class T=double, class Tstore=double>
struct Polynomial {
	typedef Polynomial<Dim, Degree, Basis, RootBasis, T, Tstore> type;
	typedef Polynomial<Dim, Degree-1, typename Basis::next_type, RootBasis, T, Tstore> next_type;
	typedef Polynomial<Dim-1, RootBasis::degree, RootBasis, RootBasis, T, Tstore> sub_type;
	typedef typename Basis::type basis_type;
//typedef typename Basis::next next_basis_type;
	enum { degree = Degree }; 
	basis_type basis;
	sub_type a;
	next_type next;
	
	template<class... Ts>
	Polynomial(sub_type a, Ts... ts) : a(a), next(ts...) {}
	//Polynomial(T& a, next_type& next) : a(a), next(next) {}
	
	template<class... Ts>
	T operator()(T x, Ts... rest) { return a(rest...) * basis(x) + next(x, rest...); }
	//type operator+(type& rhs) { return type(a + rhs.a, next + rhs.next); }
};

template<size_t Dim, class Basis, class RootBasis, class T, class Tstore>
struct Polynomial<Dim, -1, Basis, RootBasis, T, Tstore> {
	typedef Polynomial<Dim, -1, Basis, RootBasis, T> type;
	template<class... Ts>
	T operator()(T x, Ts... rest) { return 0; }
	type operator+(type& rhs) { return type(); }
};



template<size_t Degree, class Basis, class RootBasis, class T, class Tstore>
struct Polynomial<1, Degree, Basis, RootBasis, T, Tstore> {
	typedef Polynomial<1, Degree, Basis, RootBasis, T, Tstore> type;
	typedef Polynomial<1, Degree-1, typename Basis::next_type, RootBasis, T, Tstore> next_type;
	typedef typename Basis::type basis_type;
//typedef typename Basis::next next_basis_type;
	typedef Polynomial<1, Degree-1, typename Basis::derivative_type, typename RootBasis::derivative_type, T, Tstore> derivative_type;
	enum { degree = Degree }; 
	basis_type basis;
	Tstore a;
	next_type next;
	
	template<class... Ts>
	Polynomial(Tstore a, Ts... ts) : a(a), next(ts...) {}
	Polynomial(T& a, next_type& next) : a(a), next(next) {}
	
	T operator()(T x) { return a * basis(x) + next(x); }
	T dfdx(T x) { return a * basis.dfdx(x) + next.dfdx(x); }
	type operator+(type& rhs) { return type(a + rhs.a, next + rhs.next); }
};

template<class Basis, class RootBasis, class T, class Tstore>
struct Polynomial<1, -1, Basis, RootBasis, T, Tstore> {
	typedef Polynomial<1, -1, Basis, RootBasis, T, Tstore> type;
	T operator()(T) { return 0; }
	T dfdx(T) { return 0; }
	type operator+(type& rhs) { return type(); }
};


/*
template<size_t Dim, size_t Degree, class Basis=MonomialBasis<Degree, double>, class T=double, class Tstore=double>
struct NPolynomial {
	typedef NPolynomial<Dim, Degree, Basis, T, Tstore> type;
	typedef NPolynomial<Dim-1, Degree, Basis, T, Tstore> next_type;
	typedef Polynomial<Degree, Basis, T, Tstore> polynomial_type;
	//typedef Polynomial<Degree-1, typename Basis::next_type, T, Tstore> next_type;
	typedef typename Basis::type basis_type;
//typedef typename Basis::next next_basis_type;
	enum { degree = Degree }; 
	basis_type basis;
	next_type next;
	//Tstore a;
	polynomial_type p;
	
	template<class... Ts>
	NPolynomial(polynomial_type p, Ts&&... ts) : p(p), next(ts...) {}
	//Polynomial(T& a, next_type& next) : a(a), next(next) {}
	
	T operator()(T x) { return p(x) * basis(x) + next(x); }
	//type operator+(type& rhs) { return type(a + rhs.a, next + rhs.next); }
};

template<size_t Degree, class Basis, class T, class Tstore>
struct NPolynomial<-1, Degree, Basis, T, Tstore> {
	//typedef Polynomial<-1, Basis, T> type;
	T operator()(T x) { return 0; }
	//type operator+(type& rhs) { return type(); }
};


*/