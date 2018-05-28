

namespace ad {

template<class F, class G>
auto operator+(F f, G g) {
	auto result = f + g;
	auto derivative = f.dx(u) + g.dx(u);
	return make_expr(
}

};