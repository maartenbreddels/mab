#pragma once
//#include "fadiff.h"
#include <memory>
#include <iostream>
#include <cmath>
namespace ad {

#undef CT_AD

using namespace std;

template<class T>
class Neg;
template<class T1, class T2>
class Mul;
template<class F, class G>
class Div;
template<class T1, class T2>
class Sum;


template<class T>
struct Zero {
	typedef T type;
	type value() const { return 0; }
	type evaluate() const { return 0; }
	type reevaluate() const { return 0; }
	void reset(){};
	type d(int i) const { return 0; }
	template<class Tother>
	type d(Tother& x) const { return 0; }
	template<int i>
	type d() const { return  0; }
	template<class X>
	shared_ptr<Zero<type>> dx(X &x)  {
		return shared_ptr<Zero<type>>(new Zero<type>);
	}
};
template<class T>
struct One {
	typedef T type;
	template<class X>
	struct result_of_partial {
		typedef Zero<type> partial_type;
	};

	type value() const { return 1; }
	type evaluate() const { return 1; }
	type reevaluate() const { return 1; }
	void reset(){};
	type d(int i) const { return 0; }
	template<class Tother>
	type d(Tother& x) const { return 0; }
	template<int i>
	type d() const { return  0; }
	template<class X>
	shared_ptr<Zero<type>> dx(X &x)  {
		return shared_ptr<Zero<type>>(new Zero<type>);
	}
};

template<class T>
struct C {
	typedef T type;
	typedef C<T> vartype;
	template<class X>
	struct result_of_partial {
		typedef Zero<type> partial_type;
	};
	C(T c) : c(c) { }
	T c;
	type value() const { return c; }
	type evaluate() const { return c; }
	type reevaluate() const { return c; }
	void reset(){};
	type d(int i) const { return 0; }
	template<class Tother>
	type d(Tother& x) const { return 0; }
	template<int i>
	type d() const { return  0; }
	template<class X>
	shared_ptr<Zero<type>> dx(X &x)  {
		return shared_ptr<Zero<type>>(new Zero<type>);
	}

#ifdef CT_AD
	dtype dx(int i) const { return C<T>(0); }
#endif
};

template<class T1>
int ddx(T1 const*const du, shared_ptr<T1> x) {
	T1* xptr = x.get();
	return xptr == du ? 1 : 0;
}

/*template<class T1, class T2>
int ddx(T1 const*const du, shared_ptr<T2> x) {
	return  0;
}*/

template<class T, int index>
struct V;

template<class T, int index>
shared_ptr<C<T>> ddxx(V<T, index> const*const du, shared_ptr<V<T, index>> x) {
	V<T, index>* xptr = x.get();
	return shared_ptr<C<T>>(new C<T>(xptr == du ? 1 : 0));
}/*
template<class T1, class T2>
shared_ptr<C<T1>> ddxx(V<T2> const*const du, shared_ptr<V<T2>> x) {
	return shared_ptr<C<T1>>(new C<T1>(0));
}*/

template<class T1, class T2>
struct One_or_Zero {
	typedef Zero<double> type; 
	type value() { return type(); }
};

template<class T>
struct One_or_Zero<T, T> {
	typedef One<double> type; 
	type value() { return type(); }
};

template<class T1, class T2>
void testme(T1 &t1, T2 &t2) {
}

template<class T1, class T2>
auto partial(shared_ptr<T1> f, shared_ptr<T2> x) -> decltype(f->dx(x)) {
	return f->dx(x);
}

template<class T, int index>
struct V {
	enum { index_ = index };
	typedef T type;
	typedef V<T, index> vartype;
	
	template<class X>
	struct result_of_partial {
		typedef typename One_or_Zero<X, vartype>::type partial_type;
	};
	

	V(T value) : v(value) { }
	type value() const { return v; }
	type evaluate() const { return v; }
	type reevaluate() const { return v; }
	void reset(){};
	template<class Other>
	type d(shared_ptr<Other> &x) const {
		return index == Other::index_ ? 1 : 0;
	}
	/*template<class OtherType>
	shared_ptr<typename One_or_Zero<type, OtherType>::type> dx(OtherType &x) {
		One_or_Zero<type, OtherType> T;
		return T.value();
	}*/
	/*template<int indexdsa>
	shared_ptr<One<double>> dx(shared_ptr<vartype> &x) {
		//cout << "1 == 1" << endl;
		return shared_ptr<One<double>>(new One<double>());
	}

	template<class OtherType>
	shared_ptr<Zero<double>> dx(shared_ptr<OtherType> &x) {
		//cout << "1 != 0" << endl;
		return shared_ptr<Zero<double>>(new Zero<double>());
	}*/
	template<class OtherType>
	shared_ptr<typename One_or_Zero<vartype, OtherType>::type > dx(shared_ptr<OtherType> &other) {
		typedef typename One_or_Zero<vartype, OtherType>::type T;
		return shared_ptr<T>(new T());
	}
	/*template<class OtherType>
	shared_ptr<One<type>> dx(shared_ptr<OtherType> &other) {
		typedef One<type> T;
		return shared_ptr<T>(new T());
	}*/
	void operator=(T value) {
		v = value;
	}
protected:
	T v;
};

template<int i>
struct varindex {
	enum { index = i };
};

template<class T, class index>
shared_ptr<V<T, index::index>> var(T v, index i) {
	typedef V<T, index::index> Vt;
	return shared_ptr<Vt>(new Vt(v));
}

template<class T, int index>
T get(shared_ptr<V<T, index>> v) {
	return v->value();
}

template<class T>
auto get(T t) -> decltype(t->value()) {
	return t->value();
}
template<class T, int index>
void set(shared_ptr<V<T, index>> v, T value) {
	return v->operator=(value);
}




template<class T1, class T2, class Child>
class BinaryExpr {
public:
	typedef typename T1::type type; 
	shared_ptr<T1> t1;
	shared_ptr<T2> t2;
private:
	type v;
public:
	bool evaluated;
	BinaryExpr(shared_ptr<T1> t1, shared_ptr<T2> t2) : t1(t1), t2(t2), v(-2), evaluated(false) { }

	type reevaluate() {
		reset();
		return this->evaluate();
	}
	void reset() {
		evaluated = false;
		t1->reset();
		t2->reset();
	}
	type evaluate() {
		if(!evaluated) {
			t1->evaluate();
			t2->evaluate();
			v = static_cast<Child*>(this)->_evaluate();
			evaluated = true;
		}
		return value(); 
	}
	type value() const {
		return v;
	}
	//virtual type value() const { return -1; }
};

template<class T1, class T2, class X>
auto psum(T1& t1, T2& t2, X& x) -> decltype(partial(t1, x)+partial(t2, x)) {
	return t1->dx(x) + t2->dx(x);
}


template<class C>
struct simplify {
	typedef C type;
};



template<class T1, class T2>
class Sum : public BinaryExpr<T1, T2, Sum<T1, T2>> {
public:
	typedef typename T1::type type;
	typedef BinaryExpr<T1, T2, Sum<T1, T2>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename simplify<Sum<	typename T1::template result_of_partial<X>::partial_type,
						typename T2::template result_of_partial<X>::partial_type>>::type partial_type;
	};

	Sum(shared_ptr<T1> t1, shared_ptr<T2> t2) : Base(t1, t2) { }
	type _evaluate() const { return this->t1->value() + this->t2->value(); }
	template<class X>
	type d(X &x) const {
		return this->t1->d(x) + this->t2->d(x);
	}

	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type> dx(shared_ptr<X> &x)  {
		return partial(this->t1, x) + partial(this->t2, x);
	}

};


template<class T1, class T2>
shared_ptr<Sum<T1, T2>> operator+(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	typedef Sum<T1, T2> Tt;
	return shared_ptr<Tt>(new Tt(t1, t2));
}

// 0 + 0 = 0
template<class T1, class T2>
struct simplify<Sum<Zero<T1>, Zero<T2>>> {
	typedef Zero<T1> type;
};
template<class T1, class T2>
shared_ptr<Zero<T1>>operator+(shared_ptr<Zero<T1>> t1, shared_ptr<Zero<T2>> t2) {
	return t1;
}
// 0 + x = x
template<class T1, class T2>
struct simplify<Sum<Zero<T1>, T2>> {
	typedef T2 type;
};
template<class T1, class T2>
shared_ptr<T2> operator+(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return t2;
}
// x + 0 = x
template<class T1, class T2>
struct simplify<Sum<T1, Zero<T2>>> {
	typedef T1 type;
};
template<class T1, class T2>
shared_ptr<T2> operator+(shared_ptr<T2> t2, shared_ptr<Zero<T1>> t1) {
	return t2;
}
// x + x = 2x
template<class T>
struct simplify<Sum<T, T>> {
	typedef Mul<C<double>, T> type;
};
template<class T>
shared_ptr<Mul<C<double>, T>> operator+(shared_ptr<T> t1, shared_ptr<T> t2) {
	return shared_ptr<C<double>>(new C<double>(2)) * t1;
}


// adding up constants
template<class T2>
shared_ptr<Sum<C<double>, T2>> operator+(double t1, shared_ptr<T2> t2) {
	return shared_ptr<C<double>>(new C<double>(t1)) + t2;
}
template<class T2>
shared_ptr<Sum<T2, C<double>>> operator+(shared_ptr<T2> t2, double t1) {
	return t2 + shared_ptr<C<double>>(new C<double>(t1));
}


template<class T1, class T2>
class Sub : public BinaryExpr<T1, T2, Sub<T1, T2>> {
public:
	typedef typename T1::type type;
	typedef BinaryExpr<T1, T2, Sub<T1, T2>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename simplify<Sub<	typename T1::template result_of_partial<X>::partial_type,
						typename T2::template result_of_partial<X>::partial_type>>::type partial_type;
	};

	Sub(shared_ptr<T1> t1, shared_ptr<T2> t2) : Base(t1, t2) { }
	type _evaluate() const { return this->t1->value() - this->t2->value(); }
	template<class X>
	type d(X &x) const {
		return this->t1->d(x) - this->t2->d(x);
	}

	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type> dx(shared_ptr<X> &x)  {
		return partial(this->t1, x) - partial(this->t2, x);
	}

};


template<class T1, class T2>
shared_ptr<Sub<T1, T2>> operator-(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	typedef Sub<T1, T2> Tt;
	return shared_ptr<Tt>(new Tt(t1, t2));
}

// 0 - 0 = 0
template<class T1, class T2>
struct simplify<Sub<Zero<T1>, Zero<T2>>> {
	typedef Zero<T1> type;
};
template<class T1, class T2>
shared_ptr<Zero<T1>>operator-(shared_ptr<Zero<T1>> t1, shared_ptr<Zero<T2>> t2) {
	return t1;
}
// 0 - x = -x
template<class T1, class T2>
struct simplify<Sub<Zero<T1>, T2>> {
	typedef Neg<T2> type;
};
template<class T1, class T2>
shared_ptr<Neg<T2>> operator-(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return shared_ptr<Neg<T2>>(new Neg<T2>(t2));
}
// x - 0 = x
template<class T1, class T2>
struct simplify<Sub<T1, Zero<T2>>> {
	typedef T1 type;
};
template<class T1, class T2>
shared_ptr<T2> operator-(shared_ptr<T2> t2, shared_ptr<Zero<T1>> t1) {
	return t2;
}
// x - x = 0
template<class T>
struct simplify<Sub<T, T>> {
	typedef Zero<double> type;
};
template<class T>
shared_ptr<Zero<double>> operator-(shared_ptr<T> t1, shared_ptr<T> t2) {
	return shared_ptr<Zero<double>>(new Zero<double>());
}

// subtracting constants
template<class T2>
shared_ptr<Sub<C<double>, T2>> operator-(double t1, shared_ptr<T2> t2) {
	return shared_ptr<C<double>>(new C<double>(t1)) - t2;
}
template<class T2>
shared_ptr<Sub<T2, C<double>>> operator-(shared_ptr<T2> t2, double t1) {
	return t2 - shared_ptr<C<double>>(new C<double>(t1));
}


template<class F, class G, class X>
struct result_of_dF_G {
	typedef typename simplify<Mul<	typename F::template result_of_partial<X>::partial_type, G>>::type type;
};

template<class F, class G, class X>
struct result_of_F_dG {
	typedef typename simplify<Mul<F, typename G::template result_of_partial<X>::partial_type>>::type type;
};

template<class F, class G>
struct result_of_F_over_G {
	typedef typename simplify<Div<F, G>>::type type;
};

template<class F, class G>
class Mul : public BinaryExpr<F, G, Mul<F, G>> {
public:
	typedef typename F::type type;
	typedef BinaryExpr<F, G, Mul<F, G>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename	simplify<	Sum<
										typename result_of_dF_G<F, G, X>::type,
										typename result_of_dF_G<G, F, X>::type
									>
							>::type partial_type;
	};
	Mul(shared_ptr<F> t1, shared_ptr<G> t2) : Base(t1, t2) { }
	type _evaluate() const { return this->t1->value() * this->t2->value(); }
	template<class X>
	type d(X &x) const {
		return this->t1->d(x) * this->t2->value() + this->t2->d(x) * this->t1->value();
	}
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type>  dx(shared_ptr<X> &x) const {
		return this->t1->dx(x) * this->t2 + this->t2->dx(x) * this->t1;
	}
};

template<class T1, class T2>
shared_ptr<Mul<T1, T2>> operator*(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	return shared_ptr<Mul<T1, T2>>(new Mul<T1, T2>(t1, t2));
}

template<class T1, class T2>
shared_ptr<Mul<C<double>, T2>> operator*(T1 t1, shared_ptr<T2> t2) {
	return shared_ptr<C<double>>(new C<double>(t1)) * t2; 
}
template<class T2, class T1>
shared_ptr<Mul<C<T1>, T2>> operator*(shared_ptr<T2> t2, T1 t1) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) * t2; 
}

// 0 * x = 0
template<class T1, class T2>
struct simplify<Mul<Zero<T1>, T2>> {
	typedef Zero<T1> type;
};
template<class T1, class T2>
shared_ptr<Zero<T1>> operator*(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return t1;
}

// x * 0  = 0
template<class T1, class T2>
struct simplify<Mul<T1, Zero<T2>>> {
	typedef Zero<T2> type;
};
template<class T1, class T2>
shared_ptr<Zero<T1>> operator*(shared_ptr<T2> t2, shared_ptr<Zero<T1>> t1) {
	return t1;
}
// 0 * x = 0
template<class T1, class T2>
shared_ptr<Zero<T1>> operator*(shared_ptr<Zero<T1>> t1, shared_ptr<Zero<T2>> t2) {
	return t1;
}
// 1 * x = x
template<class T1, class T2>
struct simplify<Mul<One<T1>, T2>> {
	typedef T2 type;
};
template<class T1, class T2>
shared_ptr<T2> operator*(shared_ptr<One<T1>> t1, shared_ptr<T2> t2) {
	return shared_ptr<T2>(t2);
}
// x * 1 = x
template<class T1, class T2>
struct simplify<Mul<T1, One<T2>>> {
	typedef T1 type;
};
template<class T1, class T2>
shared_ptr<T2> operator*(shared_ptr<T2> t2, shared_ptr<One<T1>> t1) {
	return shared_ptr<T2>(t2);
}



template<class F, class G>
class Div : public BinaryExpr<F, G, Div<F, G>> {
public:
	typedef typename F::type type;
	typedef BinaryExpr<F, G, Div<F, G>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename	simplify<	Sub<
										typename simplify<
											Div<typename F::template result_of_partial<X>::partial_type, G>
										>::type,
										typename simplify<
											Div<typename simplify<
													Mul<F, typename G::template result_of_partial<X>::partial_type>
												>::type,
												typename simplify<
													Mul<G, G>
												>::type
											>
										>::type
									>
							>::type partial_type;
	};
	Div(shared_ptr<F> t1, shared_ptr<G> t2) : Base(t1, t2) { }
	type _evaluate() const { return this->t1->value() / this->t2->value(); }
	template<class X>
	type d(X &x) const {
		return this->t1->d(x)/ this->t2->value() - this->t1->value() * this->t2->d(x)/(this->t2->value()*this->t2->value());
	}
	//auto dx(X &x) const -> decltype(Base::t1->dx(x)/ Base::t2 - Base::t1 * Base::t2->dx(x)/(Base::t2*Base::t2)) {
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type>  dx(shared_ptr<X> &x) const {
		return this->t1->dx(x)/this->t2 - this->t1*this->t2->dx(x)/(this->t2*this->t2);
	}


	/*template<class X>
	type d(X &x) const {
		return this->t1->d(x) * this->t2->value() + this->t2->d(x) * this->t1->value();
	}
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type>  dx(shared_ptr<X> &x) const {
		return this->t1->dx(x) * this->t2 + this->t2->dx(x) * this->t1;
	}*/
};

template<class T1, class T2>
shared_ptr<Div<T1, T2>> operator/(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	return shared_ptr<Div<T1, T2>>(new Div<T1, T2>(t1, t2));
}

template<class T1, class T2>
shared_ptr<Div<C<T1>, T2>> operator/(T1 t1, shared_ptr<T2> t2) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) * t2; 
}
template<class T2, class T1>
shared_ptr<Div<C<T1>, T2>> operator/(shared_ptr<T2> t2, T1 t1) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) * t2; 
}

// 0 / x = 0
template<class T1, class T2>
struct simplify<Div<Zero<T1>, T2>> {
	typedef Zero<T1> type;
};
template<class T1, class T2>
shared_ptr<Zero<T1>> operator/(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return t1;
}

// x / 1 = x
template<class T1, class T2>
struct simplify<Div<T1, One<T2>>> {
	typedef T1 type;
};
template<class T1, class T2>
shared_ptr<T2> operator/(shared_ptr<T2> t2, shared_ptr<One<T1>> t1) {
	return shared_ptr<T2>(t2);
}

template<class T2>
shared_ptr<Div<C<double>, T2>> operator/(double t1, shared_ptr<T2> t2) {
	return shared_ptr<C<double>>(new C<double>(t1))/t2;
}
template<class T2>
shared_ptr<Div<T2, C<double>>> operator/(shared_ptr<T2> t2, double t1) {
	return t2/shared_ptr<C<double>>(new C<double>(t1));
}


/*
template<class T1, class T2>
shared_ptr<Div<C<T1>, T2>> operator/(T1 t1, shared_ptr<T2> t2) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) / t2; 
}
template<class T2, class T1>
shared_ptr<Div<C<T1>, T2>> operator/(shared_ptr<T2> t2, T1 t1) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) / t2; 
}*/


/*template<class T1, class T2>
shared_ptr<Mul<C<T1>, T2>> operator/(T1 t1, shared_ptr<T2> t2) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) * t2; 
}*/


template<class T, class Child>
class UnaryExpr {
public:
	typedef typename T::type type; 
	shared_ptr<T> t;
private:
	type v;
public:
	bool evaluated;
	UnaryExpr(shared_ptr<T> t) : t(t), v(-2), evaluated(false) { }

	type reevaluate() {
		reset();
		return this->evaluate();
	}
	void reset() {
		evaluated = false;
		t->reset();
	}
	type evaluate() {
		if(!evaluated) {
			t->evaluate();
			v = static_cast<Child*>(this)->_evaluate();
			evaluated = true;
		}
		return value();
	}
	type value() const {
		return v;
	}
	//virtual type value() const { return -1; }
};
template<class T>
class Log : public UnaryExpr<T, Log<T>> {
public:
	typedef typename T::type type;
	typedef UnaryExpr<T, Log<T>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename simplify<
									Div<
										typename T::template result_of_partial<X>::partial_type,
										Mul<C<double>, T>
									> 
								>::type partial_type;
	};
	double logbase;
	Log(shared_ptr<T> t, double logbase=M_E) : Base(t), logbase(logbase) { }
	type _evaluate() const { return ::log(this->t->value()) / ::log(logbase); }
	template<class X>
	type d(X &x) const {
		return this->t->d(x) / (this->t->value()*::log(logbase));
	}
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type>  dx(shared_ptr<X> &x) const {
		return this->t->dx(x) / (::log(logbase)*this->t);
	}

};

template<class T>
shared_ptr<Log<T>> log(shared_ptr<T> t) {
	return shared_ptr<Log<T>>(new Log<T>(t));
}

template<class T>
shared_ptr<Log<T>> log10(shared_ptr<T> t) {
	return shared_ptr<Log<T>>(new Log<T>(t,10));
}

template<class T>
shared_ptr<Log<T>> log2(shared_ptr<T> t) {
	return shared_ptr<Log<T>>(new Log<T>(t,10));
}

template<class T>
class Neg;

template<class T>
shared_ptr<Neg<T>> negate(shared_ptr<T> t) {
	return shared_ptr<Neg<T>>(new Neg<T>(t));
}

template<class T>
class Neg : public UnaryExpr<T, Neg<T>> {
public:
	typedef typename T::type type;
	typedef UnaryExpr<T, Neg<T>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename simplify<Neg<	typename T::template result_of_partial<X>::partial_type> >::type partial_type;
	};
	Neg(shared_ptr<T> t) : Base(t) { }
	type _evaluate() const { return -this->t->value(); }
	template<class X>
	type d(X &x) const {
		return - this->t->d(x);
	}
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type>  dx(shared_ptr<X> &x) const {
		return -this->t->dx(x);
	}

};

template<class T>
shared_ptr<Neg<T>> operator-(shared_ptr<T> t) {
	return shared_ptr<Neg<T>>(new Neg<T>(t));
}

template<class T>
shared_ptr<Zero<T>> operator-(shared_ptr<Zero<T>> z) {
	return z;
}



/*
template<class T1, class T2>
shared_ptr<Mul<C<T1>, T2>> operator*(T1 t1, shared_ptr<T2> t2) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) * t2; 
}*/



/*
template<class T>
struct Scale {
	typedef typename T::type type;
	typedef typename T::vartype vartype;
	Scale(double c, T t) : c(c), t(t) { }
	T t;
	double c;

	type eval() const {
		double value = c*t.eval();
		return value;
	}
	type d(int i)  const {
		return c * t.d(i);
	}
	type d(vartype &x)  const {
		return c * t.d(x);
	}
	template<int i>
	type d() const {
		return c * t.d<i>();
	}
};
template<class T>
Scale<T> operator*(double c, T t) {
	return Scale<T>(c, t);
}
template<class T>
Scale<T> operator*(int c, T t) {
	return Scale<T>(c, t);
}*/

template<class T>
class Exp : public UnaryExpr<T, Exp<T>> {
public:
	typedef typename T::type type;
	typedef UnaryExpr<T, Exp<T>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename simplify<
									Mul<Exp<T>, typename T::template result_of_partial<X>::partial_type>
								>::type partial_type;
	};
	Exp(shared_ptr<T> t) : Base(t) { }
	type _evaluate() const { return ::exp(this->t->value()); }
	template<class X>
	type d(X &x) const {
		return this->value() * this->t->d(x);
	}
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type> dx(shared_ptr<X> &x) {
		return shared_ptr<Exp<T>>(this) * this->t->dx(x);
	}

};

template<class T>
shared_ptr<Exp<T>> exp(shared_ptr<T> t) {
	return shared_ptr<Exp<T>>(new Exp<T>(t));
}



template<class F, class G>
class Pow : public BinaryExpr<F, G, Pow<F, G>> {
public:
	typedef typename F::type type;
	typedef BinaryExpr<F, G, Pow<F, G>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename	simplify<	Mul<
											Pow<F, G>,
											typename simplify<
												Sum<
													typename result_of_F_over_G<typename result_of_F_dG<G, F, X>::type, F>::type,
													typename simplify<Mul<Log<F>, typename G::template result_of_partial<X>::partial_type>>::type
												>
											>::type
										>
									>::type partial_type;
	};
	Pow(shared_ptr<F> t1, shared_ptr<G> t2) : Base(t1, t2) { }
	type _evaluate() const { return ::pow(this->t1->value(), this->t2->value()); }
	template<class X>
	type d(X &x) const {
		return this->value() * (this->t2->value() * this->t1->d(x) / this->t1->value() + ::log(this->t1->value())*this->t2->d(x)); 
	}
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type>  dx(shared_ptr<X> &x) {
		//return shared_ptr<Pow<F,G>>(this) * (this->t2 * this->t1->dx(x) / this->t1);
		return shared_ptr<Pow<F,G>>(this) * (this->t2 * this->t1->dx(x) / this->t1 + log(this->t1)*this->t2->dx(x));
		//  + ::log(this->t1->value())*this->t2->d(x)
	}
};

template<class T1, class T2>
shared_ptr<Pow<T1, T2>> pow(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	return shared_ptr<Pow<T1, T2>>(new Pow<T1, T2>(t1, t2));
}

template<class T1, class T2>
shared_ptr<Pow<C<T1>, T2>> pow(T1 t1, shared_ptr<T2> t2) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) * t2; 
}
template<class T2, class T1>
shared_ptr<Pow<T2, C<T1>>> pow(shared_ptr<T2> t2, T1 t1) {
	return pow(t2, shared_ptr<C<T1>>(new C<T1>(t1))); 
}


/*template<class T1, class T2>
class PowOld2 : public BinaryExpr<T1, T2, PowOld2<T1, T2>> {
public:
	typedef typename T1::type type;
	typedef BinaryExpr<T1, T2, PowOld2<T1, T2>> Base;
	PowOld2(shared_ptr<T1> t1, shared_ptr<T2> t2) : Base(t1, t2) { }
	type _evaluate() const { return ::pow(this->t1->value(), this->t2->value()); }
	template<class X>
	type d(X &x) const {
		return this->value() * (this->t2->value() * this->t1->d(x) / this->t1->value() + ::log(this->t1->value())*this->t2->d(x)); 
	}
};



template<class T1, class T2>
shared_ptr<Pow<T1, T2>> pow(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	return shared_ptr<Pow<T1, T2>>(new Pow<T1, T2>(t1, t2));
}

template<class T1, class T2>
shared_ptr<Pow<C<T1>, T2>> pow(T1 t1, shared_ptr<T2> t2) {
	return pow(shared_ptr<C<T1>>(new C<T1>(t1)), t2); 
}
template<class T2, class T1>
shared_ptr<Pow<T2, C<T1>>> pow(shared_ptr<T2> t2, T1 t1) {
	return pow(t2, shared_ptr<C<T1>>(new C<T1>(t1))); 
}*/

template<class T1, class T2> /* 0 ^ x == 0 (if x != 0) */
shared_ptr<Zero<T1>> pow(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return t1;
}

template<class T1, class T2> /* x^0 == 1 */
shared_ptr<Zero<T1>> pow(shared_ptr<T2> t2, shared_ptr<Zero<T1>> t1) {
	return shared_ptr<C<T1>>(new One<T1>());
}

template<class T1, class T2> /* 0 ^ 0 == 1 */
shared_ptr<Zero<T1>> pow(shared_ptr<Zero<T1>> t1, shared_ptr<Zero<T2>> t2) {
	return shared_ptr<C<T1>>(new One<T1>());
}

template<class T>
auto sqrt(shared_ptr<T> t) -> decltype(pow(t, 0.5)) {
	return pow(shared_ptr<T>(t), shared_ptr<C<double>>(new C<double>(0.5)));
}

template<class T>
auto sqr(shared_ptr<T> t) -> decltype(pow(t, 2)) {
	return pow(t, 2);
}


template<class T1, class T2>
struct PowOld {
	typedef typename T1::type type; 
	typedef typename T1::vartype vartype;
#ifdef CT_AD
	typedef Pow<T1, T2> f1;
	typedef Mul<T2, typename T1::dtype> fa1;
	typedef Div<fa1, T1> fa;
	typedef Mul<Log<T1>, typename T2::dtype> fb;
	typedef Sum<fa, fb> f2;
	typedef Mul<f1, f2> dtype;
#endif
	PowOld(T1 t1, T2 t2) : t1(t1), t2(t2) { }
	T1 t1;
	T2 t2;

	type eval() const {
		return ::pow(t1.eval(), t2.eval());
	}
	type d(int i) const {
		return eval() * (t2.eval() * t1.d(i) / t1.eval() + ::log(t1.eval())*t2.d(i)); 
	}
#ifdef CT_AD
	dtype dx(int i) const {
		return pow(t1, t2) * ((t2 * t1.dx(i)) / t1  + log(t1)*t2.dx(i)); 
	}
#endif
	type d(vartype &x) const {
		return t1.d(x) + t2.d(x);
	}
	template<int i>
	type d() const {
		return eval() * (t2.eval() * t1.d<i>() / t1.eval() + ::log(t1.eval())*t2.d<i>());
	}
};

/*
template<class T1>
struct Pow<T1, C<double>> {
	typedef typename T1::type type; 
	typedef typename T1::vartype vartype;
	Pow(T1& t1, C<double> t2) : t1(t1), t2(t2) { }
	T1& t1;
	C<double> t2;

	type eval() const {
		return ::pow(t1.eval(), t2.eval());
	}
	type d(int i) const {
		return eval() * (t2.eval() * t1.d(i) / t1.eval() + ::log(t1.eval())*t2.d(i)); 
	}
	type d(vartype &x) const {
		return t1.d(x) + t2.d(x);
	}
	template<int i>
	type d() const {
		return eval() * (t2.eval() * t1.d<i>() / t1.eval() + ::log(t1.eval())*t2.d<i>());
	}
};

template<class T2>
struct Pow<C<double>, T2> {
	typedef typename T2::type type; 
	typedef typename T2::vartype vartype;
	Pow(C<double> t1, T2& t2) : t1(t1), t2(t2) { }
	C<double> t1;
	T2& t2;

	type eval() const {
		return ::pow(t1.eval(), t2.eval());
	}
	type d(int i) const {
		return eval() * (t2.eval() * t1.d(i) / t1.eval() + ::log(t1.eval())*t2.d(i)); 
	}
	type d(vartype &x) const {
		return t1.d(x) + t2.d(x);
	}
	template<int i>
	type d() const {
		return eval() * (t2.eval() * t1.d<i>() / t1.eval() + ::log(t1.eval())*t2.d<i>());
	}
};*/




template<class T>
struct Cos;

template<class T>
shared_ptr<Cos<T>> cos(shared_ptr<T> t);


template<class T>
class Sin : public UnaryExpr<T, Sin<T>> {
public:
	typedef typename T::type type;
	typedef UnaryExpr<T, Sin<T>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename simplify<
									Mul< Cos<T>, typename T::template result_of_partial<X>::partial_type>
								>::type partial_type;
	};
	Sin(shared_ptr<T> t) : Base(t) { }
	type _evaluate() const { return ::sin(this->t->value()); }
	template<class X>
	type d(X &x) const {
		return ::cos(this->t->value()) * this->t->d(x);
	}
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type> dx(shared_ptr<X> &x) {
		return cos(this->t) * this->t->dx(x);
	}

};

template<class T>
shared_ptr<Sin<T>> sin(shared_ptr<T> t) {
	return shared_ptr<Sin<T>>(new Sin<T>(t));
}


template<class T>
class Cos : public UnaryExpr<T, Cos<T>> {
public:
	typedef typename T::type type;
	typedef UnaryExpr<T, Cos<T>> Base;
	template<class X>
	struct result_of_partial {
		typedef typename simplify<
									Mul< Neg<Sin<T>>, typename T::template result_of_partial<X>::partial_type>
								>::type partial_type;
	};
	Cos(shared_ptr<T> t) : Base(t) { }
	type _evaluate() const { return ::cos(this->t->value()); }
	template<class X>
	type d(X &x) const {
		return -::sin(this->t->value()) * this->t->d(x);
	}
	template<class X>
	shared_ptr<typename result_of_partial<X>::partial_type> dx(shared_ptr<X> &x) {
		return -sin(this->t) * this->t->dx(x);
	}

};

template<class T>
shared_ptr<Cos<T>> cos(shared_ptr<T> t) {
	return shared_ptr<Cos<T>>(new Cos<T>(t));
}


}

