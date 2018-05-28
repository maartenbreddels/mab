#pragma once
//#include "fadiff.h"
#include <memory>
#include <iostream>
#include <cmath>
namespace ad {

#undef CT_AD

using namespace std;

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
#ifdef CT_AD
	typedef C<T> dtype;
#endif
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

template<class T1, class T2>
class Sum : public BinaryExpr<T1, T2, Sum<T1, T2>> {
public:
	typedef typename T1::type type;
	typedef BinaryExpr<T1, T2, Sum<T1, T2>> Base; 

	Sum(shared_ptr<T1> t1, shared_ptr<T2> t2) : Base(t1, t2) { }
	type _evaluate() const { return this->t1->value() + this->t2->value(); }
	template<class X>
	type d(X &x) const {
		return this->t1->d(x) + this->t2->d(x);
	}

	template<class X>
	auto dx(X &x) -> decltype(psum(Base::t1, Base::t2, x)) {
		return psum(this->t1, this->t2, x);
	}

};
template<class T1, class T2>
shared_ptr<Sum<T1, T2>> operator+(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	typedef Sum<T1, T2> Tt;
	return shared_ptr<Tt>(new Tt(t1, t2));
}
template<class T1, class T2>
shared_ptr<Zero<T1>>operator+(shared_ptr<Zero<T1>> t1, shared_ptr<Zero<T2>> t2) {
	return t1;
}
template<class T1, class T2>
shared_ptr<T2> operator+(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return t2;
}
template<class T1, class T2>
shared_ptr<T2> operator+(shared_ptr<T2> t2, shared_ptr<Zero<T1>> t1) {
	return t2;
}

template<class T>
Sum<T, C<double>> operator+(T t, double c) {
	return Sum<T, C<double>>(t, C<double>(c));
}
template<class T>
Sum<C<double>, T> operator+(double c, T t) {
	return Sum<C<double>, T>(C<double>(c), t);
}
template<class T>
Sum<T, C<double>> operator+(T t, int c) {
	return Sum<T, C<double>>(t, C<double>(c));
}
template<class T>
Sum<C<double>, T> operator+(int c, T t) {
	return Sum<C<double>, T>(C<double>(c), t);
}

template<class T1, class T2>
class Sub : public BinaryExpr<T1, T2, Sub<T1, T2>> {
public:
	typedef typename T1::type type;
	typedef BinaryExpr<T1, T2, Sub<T1, T2>> Base; 

	Sub(shared_ptr<T1> t1, shared_ptr<T2> t2) : Base(t1, t2) { }
	type _evaluate() const { return this->t1->value() - this->t2->value(); }
	template<class X>
	type d(X &x) const {
		return this->t1->d(x) - this->t2->d(x);
	}

	template<class X>
	auto dx(X &x) -> decltype(Base::t1->dx(x) - Base::t2->dx(x)) {
		return this->t1->dx(x) - this->t2->dx(x);
	}

};
template<class T1, class T2>
shared_ptr<Sub<T1, T2>> operator-(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	typedef Sub<T1, T2> Tt;
	return shared_ptr<Tt>(new Tt(t1, t2));
}
template<class T1, class T2>
shared_ptr<Zero<T1>>operator-(shared_ptr<Zero<T1>> t1, shared_ptr<Zero<T2>> t2) {
	return t1;
}
template<class T1, class T2>
shared_ptr<T2> operator-(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return t2;
}
template<class T1, class T2>
shared_ptr<T2> operator-(shared_ptr<T2> t2, shared_ptr<Zero<T1>> t1) {
	return t2;
}



template<class T1, class T2>
class Mul : public BinaryExpr<T1, T2, Mul<T1, T2>> {
public:
	typedef typename T1::type type;
	typedef BinaryExpr<T1, T2, Mul<T1, T2>> Base;
	Mul(shared_ptr<T1> t1, shared_ptr<T2> t2) : Base(t1, t2) { }
	type _evaluate() const { return this->t1->value() * this->t2->value(); }
	template<class X>
	type d(X &x) const {
		return this->t1->d(x) * this->t2->value() + this->t2->d(x) * this->t1->value();
	}
	template<class X>
	auto dx(X &x) const -> decltype(partial(Base::t1, x) * Base::t2 + partial(Base::t2, x) * Base::t1) {
		return this->t1->dx(x) * this->t2 + this->t2->dx(x) * this->t1;
	}
};

template<class T1, class T2>
shared_ptr<Mul<T1, T2>> operator*(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	return shared_ptr<Mul<T1, T2>>(new Mul<T1, T2>(t1, t2));
}

template<class T1, class T2>
shared_ptr<Mul<C<T1>, T2>> operator*(T1 t1, shared_ptr<T2> t2) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) * t2; 
}
template<class T2, class T1>
shared_ptr<Mul<C<T1>, T2>> operator*(shared_ptr<T2> t2, T1 t1) {
	return shared_ptr<C<T1>>(new C<T1>(t1)) * t2; 
}

template<class T1, class T2>
shared_ptr<Zero<T1>> operator*(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return t1;
}

template<class T1, class T2>
shared_ptr<Zero<T1>> operator*(shared_ptr<T2> t2, shared_ptr<Zero<T1>> t1) {
	return t1;
}

template<class T1, class T2>
shared_ptr<Zero<T1>> operator*(shared_ptr<Zero<T1>> t1, shared_ptr<Zero<T2>> t2) {
	return t1;
}

template<class T1, class T2>
class Div : public BinaryExpr<T1, T2, Div<T1, T2>> {
public:
	typedef typename T1::type type;
	typedef BinaryExpr<T1, T2, Div<T1, T2>> Base;
	Div(shared_ptr<T1> t1, shared_ptr<T2> t2) : Base(t1, t2) { }
	type _evaluate() const { return this->t1->value() / this->t2->value(); }
	template<class X>
	type d(X &x) const {
		return this->t1->d(x)/ this->t2->value() - this->t1->value() * this->t2->d(x)/(this->t2->value()*this->t2->value());
	}
	//auto dx(X &x) const -> decltype(Base::t1->dx(x)/ Base::t2 - Base::t1 * Base::t2->dx(x)/(Base::t2*Base::t2)) {
	template<class X>
	auto dx(X &x) const -> decltype(partial(Base::t1, x)/Base::t2 - Base::t1 * partial(Base::t2, x) / (Base::t2*Base::t2)) {
		return this->t1->dx(x)/this->t2 - this->t1*this->t2->dx(x)/(this->t2*this->t2);
	}

};

template<class T1, class T2>
shared_ptr<Div<T1, T2>> operator/(shared_ptr<T1> t1, shared_ptr<T2> t2) {
	return shared_ptr<Div<T1, T2>>(new Div<T1, T2>(t1, t2));
}

template<class T1, class T2>
shared_ptr<Zero<T1>> operator/(shared_ptr<Zero<T1>> t1, shared_ptr<T2> t2) {
	return t1;
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
	type value() {
		return v;
	}
	//virtual type value() const { return -1; }
};
template<class T>
class Log : public UnaryExpr<T, Log<T>> {
public:
	typedef typename T::type type;
	typedef UnaryExpr<T, Log<T>> Base;
	double logbase;
	Log(shared_ptr<T> t, double logbase=M_E) : Base(t), logbase(logbase) { }
	type _evaluate() const { return ::log(this->t->value()) / ::log(logbase); }
	template<class X>
	type d(X &x) const {
		return this->t->d(x) / (this->v);
	}
	template<class X>
	auto dx(X &x) const -> decltype(Base::t->dx(x) / Base::t) {
		return this->t->dx(x) / this->t;
	}

};

template<class T>
shared_ptr<Log<T>> log(shared_ptr<T> t) {
	return shared_ptr<Log<T>>(new Log<T>(t));
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
	Neg(shared_ptr<T> t) : Base(t) { }
	type _evaluate() const { return -this->t->value(); }
	template<class X>
	type d(X &x) const {
		return - this->t->d(x);
	}
	template<class X>
	auto dx(X &x) -> decltype(-partial(Base::t, x)) {
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
struct Exp {
	typedef typename T::type type;
	typedef typename T::vartype vartype;
#ifdef CT_AD
	typedef Mul<T, typename T::dtype> dtype;
#endif
	Exp(T t) : t(t) { }
	T t;

	type eval() const {
		return ::exp(t.eval());
	}
	type d(int i) const {
		return eval() * t.d(i);
	}
	type d(vartype &x) const {
		return eval() * t.d(x);
	}
#ifdef CT_AD
	dtype dx(int i) const {
		return t * t.dx(i);
	}
#endif
	template<int i>
	type d() const {
		return eval() * t.d<i>();
	}
};
template<class T>
Exp<T> exp(T t) {
	return Exp<T>(t);
}

template<class T1, class T2>
class Pow : public BinaryExpr<T1, T2, Pow<T1, T2>> {
public:
	typedef typename T1::type type;
	typedef BinaryExpr<T1, T2, Pow<T1, T2>> Base;
	Pow(shared_ptr<T1> t1, shared_ptr<T2> t2) : Base(t1, t2) { }
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
}

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
	return pow(t, 0.5);
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
struct Sin {
	typedef typename T::type type;
	typedef typename T::vartype vartype;
#ifdef CT_AD
	typedef Mul<Cos<T>, typename T::dtype> dtype;
#endif
	Sin(T t) : t(t) { }
	T t;

	type eval() const {
		return ::sin(t.eval());
	}
	type d(int i) const {
		return ::cos(t.eval()) * t.d(i);
	}
#ifdef CT_AD
	dtype dx(int i) const {
		return cos(t) * t.dx(i);
	}
#endif
	type d(vartype &x) const {
		return ::cos(t.eval()) * t.d(x);
	}
	template<int i>
	type d() const {
		return ::cos(t.eval()) * t.d<i>();
	}
};
template<class T>
Sin<T> sin(T t) {
	return Sin<T>(t);
}

template<class T>
struct Cos {
	typedef typename T::type type;
	typedef typename T::vartype vartype;
#ifdef CT_AD
	typedef Neg<Mul<Sin<T>, typename T::dtype>> dtype;
#endif
	Cos(T t) : t(t) { }
	T t;

	type eval() const {
		return ::cos(t.eval());
	}
	type d(int i) const {
		return -::sin(t.eval()) * t.d(i);
	}
#ifdef CT_AD
	dtype dx(int i) const {
		return -(sin(t) * t.dx(i));
	}
#endif
	type d(vartype &x) const {
		return -::sin(t.eval()) * t.d(x);
	}
	template<int i>
	type d() const {
		return -::sin(t.eval()) * t.d<i>();
	}
	/*Sin<T> dx(int i) {
		return -Sin<T>(t.eval()) * t.dx<i>();
	}*/
};
template<class T>
Cos<T> cos(T t) {
	return Cos<T>(t);
}

}

