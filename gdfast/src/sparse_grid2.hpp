#include "cpptools.hpp"


/*
template<int I, typename Condition, bool test, typename Seq>
struct seq_inter;

template<int I, int ...S>
struct seq_inter<I, Condition, true, seq<S...>> {
	seq_inter<
};
*/

template<int el, typename Seq>
struct get;
template<int el, int Head, int... Tail>
struct get<el, seq<Head, Tail...>> {
	enum { value = get<el-1, seq<Tail...>>::value };
};
template<int Head, int... Tail>
struct get<0, seq<Head, Tail...>> {
	enum { value = Head };
};



template<typename Seq>
struct seq_sum;
template<int Head, int ...Tail>
struct seq_sum<seq<Head, Tail...>> {
	enum { value = Head + seq_sum<seq<Tail...>>::value };
};
template<>
struct seq_sum<seq<>> {
	enum { value = 0 };
};

template<int D, int Level, typename Seq>
struct point;
template<int D, int Level, int ...S>
struct point<D, Level, seq<S...>> {
	typedef seq<S...> Seq;
	//typedef next_point<D, Level, seq_iter<seq<S...>>::type> next;
	void print() {
		println(S...);
		println(seq_sum<seq<S...>>::value);
	}
};


template<typename... Seqs>
struct seq_product {
};

/*
template<int... S, typename... SeqsTail>
struct seq_product<seq<S...>, SeqsTail...> {
	void test() {
		println(">>>>>>>>", S...);
	}
};
*/

template<typename P>
struct seq_product_apply;

template<int Head, int... Tail, typename... SeqsTail>
struct seq_product_apply<  seq_product<seq<Head, Tail...>, SeqsTail...>  > {
	template<typename F, typename... Ts>
	void operator()(F f, Ts... ts) {
		seq_product_apply<seq_product<SeqsTail...>> next2;
		next2(f, Head, ts...);
		seq_product_apply<seq_product<seq<Tail...>, SeqsTail...>> next1;
		next1(f, ts...);
	}
};

template<typename... SeqsTail>
struct seq_product_apply<  seq_product<seq<>, SeqsTail...>  > {
	template<typename F, typename... Ts>
	void operator()(F f, Ts... ts) {
	}
};


template<int Head, int... Tail>
struct seq_product_apply< seq_product<seq<Head, Tail...>> > {
	template<typename F, typename... Ts>
	void operator()(F f, Ts... ts) {
		f(ts..., Head);
		seq_product_apply<seq_product<seq<Tail...>>> next1;
		next1(f, ts...);
	}
};
template<>
struct seq_product_apply< seq_product<seq<>> > {
	template<typename F, typename... Ts>
	void operator()(F f, Ts... ts) {
	}
};


template<int Level, int Left, typename P>
struct seq_product_apply_sparse;

template<int Level, int Left, int Head, int... Tail, typename... SeqsTail>
struct seq_product_apply_sparse<Level, Left, seq_product<seq<Head, Tail...>, SeqsTail...>  > {
	template<typename F, typename... Ts>
	void operator()(F f, Ts... ts) {
		if(Left-Head >= 0) {
			seq_product_apply_sparse<Level, Left-Head, seq_product<SeqsTail...>> next2;
			next2(f, Head, ts...);
			seq_product_apply_sparse<Level, Left, seq_product<seq<Tail...>, SeqsTail...>> next1;
			next1(f, ts...);
		}
	}
};

template<int Level, int Left, typename... SeqsTail>
struct seq_product_apply_sparse<Level, Left, seq_product<seq<>, SeqsTail...>  > {
	template<typename F, typename... Ts>
	void operator()(F f, Ts... ts) {
	}
};


template<int Level, int Left, int Head, int... Tail>
struct seq_product_apply_sparse<Level, Left, seq_product<seq<Head, Tail...>> > {
	template<typename F, typename... Ts>
	void operator()(F f, Ts... ts) {
		//println(Left, Head, Left+Head);
		//if(Left-Head==0)
		f(ts..., Left);
		//seq_product_apply_sparse<Level, Left, seq_product<seq<Tail...>>> next1;
		//next1(f, ts...);
	}
};
template<int Level, int Left>
struct seq_product_apply_sparse<Level, Left, seq_product<seq<>> > {
	template<typename F, typename... Ts>
	void operator()(F f, Ts... ts) {
	}
};


template<int Norm, typename Seq>
struct norm1_test;

template<int Norm, int... S>
struct norm1_test<Norm, seq<S...>> {
	enum { value = Norm == seq_sum<seq<S...>>::value };
};

