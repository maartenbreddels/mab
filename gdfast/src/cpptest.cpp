#include "cpptools.hpp"
#include <iostream>
#include <tuple>

using namespace std;

template<typename ...Args>
void call(Args... args) {
}
tuple<double, int> t;

template<int Head, int ...Tail>
void printseq(seq<Head, Tail...> q) {
	//cout << "bla" << endl;
	cout << Head << ", ";
	seq<Tail...> b;
	printseq(b);
}
template<int X>
void printseq(seq<X> q) {
	//cout << "last" << endl;
	cout << X << endl;
}

//template<int ...Seq>
//template< seq<typename X, typename Y> >
template<int ...X>
void test(seq<X...> q) {
	printseq(q);
	//cout << (Seq) << endl;
	//call(get<Seq>(t)...);
}

int main() {
	genseq<5>::type a;
	test(a);
	seq<0,1,2,3> s;
	test(s);
	seq_prepend<9, seq<0,1,2,3>> ::type s2;
	test(s2);
	seq_append<11, seq<0,1,2,3>> ::type s3;
	test(s3);
}