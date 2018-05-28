#include <iostream>

void println() {
	std::cout << std::endl;
}

template<typename Head, typename ...Tail>
void println(Head head, Tail... tail) {
	std::cout << head << " ";
	println(tail...);
}


template<int ...T>
struct seq {
};

template<typename S>
struct seq_apply;

template<int Head, int... Tail>
struct seq_apply<seq<Head, Tail...>> {
	template<typename F>
	void operator()(F f) {
		f(Head);
		seq_apply<seq<Tail...>> next;
		next(f);
	}
};

template<>
struct seq_apply<seq<>> {
	template<typename F>
	void operator()(F f) {
	}
};

template<int Head, typename S>
struct seq_prepend;

template<int Head, int ...Tail>
struct seq_prepend<Head, seq<Tail...>>
{
	typedef seq<Head, Tail...> type;
};



template<int Last, typename S>
struct seq_append;

template<int Last, int ...Head>
struct seq_append<Last, seq<Head...>>
{
	typedef seq<Head..., Last> type;
};





template<int N, int ...Seq>
struct genseq;

template<>
struct genseq<0>
{
	typedef seq<0> type;
};

template<int N, int ...Seq>
struct genseq
{
	typedef typename seq_append<N, typename genseq<N-1>::type>::type type;
};






