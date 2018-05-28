#include "sparse_grid2.hpp"


struct cat {
	bool operator()(int i) {
		println("i = ", i);
		return true;
	};
};
struct cat2 {
	bool operator()(int i, int j) {
		println("i = ", i, " j = ", j, " norm = ", i + j);
		return true;
	};
};
struct cat3 {
	bool operator()(int i, int j, int k) {
		println("i = ", i, " j = ", j, " k = ", k, " norm = ", i + j + k);
		return true;
	};
};

int main() {
	typedef point<2, 3, seq<0, 1, 2>> ptype;
	ptype p;
	p.print();
	println(">>", get<0, ptype::Seq>::value);
	println(">>", get<1, ptype::Seq>::value);
	println(">>", get<2, ptype::Seq>::value);
	//ptype::next pnext;
	//pnext.print()
	println(">>", norm1_test<2, ptype::Seq>::value);
	println(">>", norm1_test<3, ptype::Seq>::value);
	println(">>", norm1_test<4, ptype::Seq>::value);
	seq_apply<seq<0,1,2>> a;
	a(cat());
	typedef seq_product<seq<0,1,2,3>, seq<0,1,2,3>> prtype;
	//prtype pr;
	//pr.test();
	seq_product_apply_sparse<3, 3, prtype> a2;
	
	a2(cat2());

	typedef seq_product<seq<0,1,2,3>, seq<0,1,2,3>, seq<0,1,2,3>> prtype3;
	seq_product_apply_sparse<3, 2, prtype3> a3;
	a3(cat3());
	//iter<norm1_test, 3>
	//p.next().print();
}