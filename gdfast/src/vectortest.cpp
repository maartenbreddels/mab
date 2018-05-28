#include "vector.hpp"
#include <iostream>

using namespace std;



int main() {
	typedef lagrange_basis<2> basis;
	vector<basis> v(1.,1.);
	for(int i = 0; i < 5; i++) {
		double x = 0.25 * i;
		double y = v(x);
		cout << x << " = " << y << endl;
	}
	//cout << (1 + 2 * pow(x,1) + 3 * pow(x,2))  << endl;
}