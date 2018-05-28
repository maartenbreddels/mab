#pragma once
//#include "fadiff.h"
#include <memory>
#include <iostream>
#include <cmath>
namespace ad {

#undef CT_AD

template<class T>
struct Zero {
	typedef T type;
	type eval(int i) { return 0; }
	
};
using namespace std;

template<class T>
struct One {
	typedef T type;
	typedef Zero<type> derivative;
	type eval(int i) { return 1; }
	
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
	type eval() const { return c; }
	type d(int i) const { return 0; }
	template<class Tother>
	type d(Tother& x) const { return 0; }
	template<int i>
	type d() const { return  0; }
#ifdef CT_AD
	dtype dx(int i) const { return C<T>(0); }
#endif
};

template<class T, int index_>
struct V {
	enum { index = index_ };
	typedef T type;
	typedef V<T, index_> vartype;
#ifdef CT_AD
	typedef C<T> dtype;
#endif
	V(T value) : v(new double) { *this->v = value; }
	//V(const V<T, index_>& other) v(other.v) { cout << "c1 "; this->v = other.v;}
	//V(V<T, index_>& other) { cout << "c2 "; this->v = other.v;}

	type eval() const { return *v; }
	//type d(int i) { derivative der; return i == index ? der.eval(i) : 0; }
	type d(int i) const { /*cout << "[1: " << index << "," << i << "," << eval() << "]" << endl;*/ return i == index_ ? 1 : 0; }
	template<class Tother>
	type d(Tother& x) const { return (void*)&x == (void*)this? 1 : 0; }
	template<int i>
	type d() const { /*cout << "[2: " << index << "," << i << "," << eval() << "]" << endl;*/ return i == index ? 1 : 0; }

#ifdef CT_AD
	dtype dx(int i) const { return C<T>(i == index ? 1 : 0); }
#endif

	void operator=(T value) {
		*v = value;
	}
protected:
	std::shared_ptr<T> v;

};



template<class T1, class T2>
struct Sum;

template<class T>
struct Neg {
	typedef typename T::type type; 
	typedef typename T::vartype vartype;
#ifdef CT_AD
 	typedef Neg<typename T::dtype> dtype;
#endif
	Neg(T& t) : t(t) { }
	//Sum(T1 t1, T2 t2) : t1(t1), t2(t2) { }
	T t;

	type eval() const {
		return- t.eval();
	}
#ifdef CT_AD
	dtype dx(int i) const {
		return -t.dx(i);
	}
#endif
	type d(int i) const {
		return -t.d(i);
	}
	type d(vartype &x) const {
		return -t.d(x);
	}
	template<int i>
	type d() const {
		return -t.d<i>();
	}
};
template<class T>
Neg<T> operator-(T t) {
	return Neg<T>(t);
}

template<class T1, class T2>
struct Sum {
	typedef typename T1::type type; 
	typedef typename T1::vartype vartype;
#ifdef CT_AD
 	typedef Sum<typename T1::dtype, typename T2::dtype> dtype;
#endif
	Sum(T1 t1, T2 t2) : t1(t1), t2(t2) { }
	//Sum(T1 t1, T2 t2) : t1(t1), t2(t2) { }
	T1 t1;
	T2 t2;

	type eval() const {
		return t1.eval() + t2.eval();
	}
#ifdef CT_AD
	dtype dx(int i) const {
		return t1.dx(i) + t2.dx(i);
	}
#endif
	type d(int i) const {
		return t1.d(i) + t2.d(i);
	}
	type d(vartype &x) const {
		return t1.d(x) + t2.d(x);
	}
	template<int i>
	type d() const {
		return t1.d<i>() + t2.d<i>();
	}
};
template<class T1, class T2>
Sum<T1, T2> operator+(T1 t1, T2 t2) {
	return Sum<T1, T2>(t1, t2);
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
struct Sub {
	typedef typename T1::type type; 
	typedef typename T1::vartype vartype;
#ifdef CT_AD
 	typedef Sub<typename T1::dtype, typename T2::dtype> dtype;
#endif
	Sub(T1 t1, T2 t2) : t1(t1), t2(t2) { }
	//Sum(T1 t1, T2 t2) : t1(t1), t2(t2) { }
	T1 t1;
	T2 t2;

	type eval() const {
		return t1.eval() - t2.eval();
	}
#ifdef CT_AD
	dtype dx(int i) const {
		return t1.dx(i) - t2.dx(i);
	}
#endif
	type d(int i) const {
		return t1.d(i) - t2.d(i);
	}
	type d(vartype &x) const {
		return t1.d(x) - t2.d(x);
	}
	template<int i>
	type d() const {
		return t1.d<i>() - t2.d<i>();
	}
};
template<class T1, class T2>
Sub<T1, T2> operator-(T1 t1, T2 t2) {
	return Sub<T1, T2>(t1, t2);
}
template<class T>
Sub<T, C<double>> operator-(T t, double c) {
	return Sub<T, C<double>>(t, C<double>(c));
}
template<class T>
Sub<C<double>, T> operator-(double c, T t) {
	return Sub<C<double>, T>(C<double>(c), t);
}
template<class T>
Sub<T, C<double>> operator-(T t, int c) {
	return Sub<T, C<double>>(t, C<double>(c));
}
template<class T>
Sub<C<double>, T> operator-(int c, T t) {
	return Sub<C<double>, T>(C<double>(c), t);
}


template<class T1, class T2>
struct Mul {
	typedef typename T1::type type; 
	typedef typename T1::vartype vartype;
#ifdef CT_AD
	typedef Sum<Mul<typename T1::dtype, T2>, Mul<typename T2::dtype, T1>> dtype;
#endif
	Mul(T1 t1, T2 t2) : t1(t1), t2(t2) { }
	T1 t1;
	T2 t2;

	type eval() const {
		return t1.eval() * t2.eval();
	}
	type d(int i) const {
		//cout << "t" << t1.d(i) << endl;
		return t1.d(i) * t2.eval() + t2.d(i) * t1.eval(); 
	}
#ifdef CT_AD
	dtype dx(int i) const {
		//cout << "t" << t1.d(i) << endl;
		return t1.dx(i) * t2 + t2.dx(i) * t1; 
	}
#endif
	type d(vartype &x) const {
		return t1.d(x) * t2.eval() + t2.d(x) * t1.eval();
	}
	template<int i>
	type d() const {
		return t1.d<i>() * t2.eval() + t2.d<i>() * t1.eval(); 
	}
};
template<class T1, class T2>
const Mul<T1, T2> operator*(T1 t1, T2 t2) {
	return Mul<T1, T2>(t1, t2);
}
template<class T>
Mul<C<double>, T> operator*(double c, T t) {
	return Mul<C<double>, T>(C<double>(c), t);
}
template<class T>
Mul<C<double>, T> operator*(int c, T t) {
	return Mul<C<double>, T>(C<double>(c), t);
}
template<class T>
Mul<C<double>, T> operator*(T t, double c) {
	return Mul<C<double>, T>(C<double>(c), t);
}
template<class T>
Mul<C<double>, T> operator*(T t, int c) {
	return Mul<C<double>, T>(C<double>(c), t);
}


template<class T1, class T2>
struct Div {
	typedef typename T1::type type; 
	typedef typename T1::vartype vartype;
#ifdef CT_AD
	typedef Div<typename T1::dtype, T2> fa;
	typedef Div<T1, Mul<T2, T2>> fb1;
	typedef Mul<fb1, typename T2::dtype> fb;
	typedef Sub<fa, fb> dtype;
#endif
	Div(T1 t1, T2 t2) : t1(t1), t2(t2) { }
	T1 t1;
	T2 t2;

	type eval() const {
		return t1.eval() / t2.eval();
	}
	type d(int i) const {
		//cout << "t" << t1.d(i) << endl;
		return t1.d(i) / t2.eval() - t1.eval()/(t2.eval() *t2.eval()) * t2.d(i);
	}
#ifdef CT_AD
	dtype dx(int i) const {
		//cout << "t" << t1.d(i) << endl;
		return t1.dx(i) / t2 - (t1/(t2*t2)) * t2.dx(i); 
	}
#endif
	type d(vartype &x) const {
		return t1.d(x) / t2.eval() - t1.eval()/(t2.eval() *t2.eval()) * t2.d(x);
	}
	template<int i>
	type d() const {
		return t1.d<i>() / t2.eval() - t1.eval()/(t2.eval() *t2.eval()) * t2.d<i>(); 
	}
};

template<class T1, class T2>
const Div<T1, T2> operator/(T1 t1, T2 t2) {
	return Div<T1, T2>(t1, t2);
}
template<class T>
Div<C<double>, T> operator/(double c, T t) {
	return Div<C<double>, T>(C<double>(c), t);
}
template<class T>
Div<C<double>, T> operator/(int c, T t) {
	return Div<C<double>, T>(C<double>(c), t);
}
template<class T>
Div<T, C<double>> operator/(T t, double c) {
	return Div<T, C<double>>(t, C<double>(c));
}
template<class T>
Div<T, C<double>> operator/(T t, int c) {
	return Div<T, C<double>>(t, C<double>(c));
}


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


template<class T>
struct Log {
	typedef typename T::type type;
	typedef typename T::vartype vartype;
#ifdef CT_AD
	typedef Div<Div<typename T::dtype, T>, C<double>> dtype; 
#endif
	Log(T t, double base) : t(t), base(base) { }
	T t;
	double base;

	type eval() const {
		return ::log(t.eval()) / ::log(base);
	}
	type d(int i) const {
		return 1/t.eval() * t.d(i) / ::log(base);
	}
#ifdef CT_AD
	dtype dx(int i) const {
		return t.dx(i)/t / ::log(base);
	}
#endif
	type d(vartype &x) const {
		return 1/t.eval() * t.d(x) / ::log(base);
	}
	template<int i>
	type d() const {
		return 1/t.eval() * t.d<i>() / ::log(base);
	}
};
template<class T>
Log<T> log(T t) {
	return Log<T>(t, M_E);
}
template<class T>
Log<T> log10(T t) {
	return Log<T>(t, 10);
}
template<class T>
Log<T> log2(T t) {
	return Log<T>(t,2);
}

template<class T1, class T2>
struct Pow {
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
	Pow(T1 t1, T2 t2) : t1(t1), t2(t2) { }
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

template<class T1, class T2>
const Pow<T1, T2> pow(T1 t1, T2 t2) {
	return Pow<T1, T2>(t1, t2);
}

/* special case of variable^constant */

template<class T>
const Pow<T, C<double>> pow(T t, int c) {
	return Pow<T, C<double>>(t, C<double>(c));
}
template<class T>
const Pow<T, C<double>> pow(T t, double c) {
	return Pow<T, C<double>>(t, C<double>(c));
}
template<class T>
const Pow<T, C<double>> sqrt(T t) {
	return Pow<T, C<double>>(t, C<double>(0.5));
}
template<class T>
const Pow<T, C<double>> sqr(T t) {
	return Pow<T, C<double>>(t, C<double>(2));
}

/* special case of constant^variable */
template<class T>
const Pow<C<double>, T> pow(int c, T t) {
	return Pow<C<double>, T>(C<double>(c), t);
}
template<class T>
const Pow<C<double>, T> pow(double c, T t) {
	return Pow<C<double>, T>(C<double>(c), t);
}


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

