#pragma once

#include <iostream>
 

namespace gd {
using namespace std;

void cout_many() {
}
template<typename Head, typename... Tail>
void cout_many(Head head, Tail... tail...) {
	cout << " " << head << " ";
	//cout_many(head);
	cout_many(tail...);
}

template<int Level>
class BasisSetTriHierLevel {
public:
	static constexpr int size = (1<<(Level-1)); //+1;
	static constexpr double hl = 1./size; //pow(2, -(Level-1));
	template<int I>
	class BasisFunction{
	public:
		static constexpr int size = (1<<(Level-1)); //+1;
		static constexpr double hl = 1./size; //pow(2, -(Level-1));
		static constexpr double center = (1./size) * (I+0.5);
		double operator()(double x) {
			const double hl = pow(2, -Level);
			int j = I * 2 + 1;
			if((x > 1) || (x < 0)) // outside range?
				return 0;
			else {
				return x > (j+1)*hl ? 0 : (x < (j-1)*hl ? 0 : 1-fabs(x/hl-j));
			}
		}
	};
};


template<>
class BasisSetTriHierLevel<0> {
public:
	static const int size = 2;//
	static constexpr double hl = 1./size; //pow(2, -(Level-1));
	template<int I>
	class BasisFunction{
	public:
		static const int size = 2;//
		static constexpr double hl = 1./size; //pow(2, -(Level-1));
		static constexpr double center = I;
		double operator()(double x) {
			if((x > 1) || (x < 0)) // outside range or odd, skip
				return 0;
			else {
				if(I == 0)
					return x > 1 ? 0 : (x < 0 ? 0 : 1-x);
				else if(I == 1)
					return x > 1 ? 0 : (x < 0 ? 0 : x);
				else
					return 0;
			}
		}
	};
};
/*
template<template<int> class BasisSetTemplate, int Level>
class Vector {
public:
	typedef BasisSetTemplate<Level> BasisSet;
	static const int size = BasisSet::size;
	template<int I=0, int Count=BasisSet::size>
	class Coordinate {
	public:
		typedef Coordinate<I+1, Count> NextCoordinate;
		typedef typename BasisSet::template BasisFunction<I> BasisFunction;
		NextCoordinate nextCoordinate;
		BasisFunction basisFunction;
		double operator()(double x, double* coefficients) {
			return basisFunction(x) * coefficients[0] + nextCoordinate(x, &coefficients[1]); 
		}
		double operator()(double x) {
			return basisFunction(x) + nextCoordinate(x); 
		}
	};
	template<int Count>
	class Coordinate<Count,Count> {
	public:
		double operator()(double x, double* coordinates) {
			return 0;
		}
		double operator()(double x) {
			return 0;
		}
	};
	Coordinate<> coordinates;
	//double coordinate[size];
	double operator()(double x, double *coefficients) {
		return coordinates(x, &coefficients[0]);
	}
	double operator()(double x) {
		return coordinates(x);
	}
};*/

template<int Dim, template<int> class BasisSetTemplate, int Norm, int NormLeft=Norm, int I=0>
class Product;

template<int Dim, template<int> class BasisSetTemplate, int Norm, int NormLeft, int I>
class Product {
public:
	Product<Dim, BasisSetTemplate, Norm, NormLeft-1, I+1> next;
	static const int Level = I;
	//Vector<BasisSetTemplate, Level> vector;
	typedef BasisSetTemplate<Level> BasisSet;
	static const int size = BasisSet::size;
	template<int J=0, int Count=BasisSet::size>
	class Coordinate {
	public:
		typedef Coordinate<J+1, Count> NextCoordinate;
		typedef typename BasisSet::template BasisFunction<J> BasisFunction;
		Product<Dim-1, BasisSetTemplate, Norm, NormLeft, 0> dimnext;
		NextCoordinate nextCoordinate;
		BasisFunction basisFunction;
		/*double operator()(double x, double* coefficients) {
			return basisFunction(x) * coefficients[0] + nextCoordinate(x, &coefficients[1]); 
		}*/
		template<typename... X>
		double operator()(double xhead, X... xtail) {
			//cout << "Dim="<<Dim <<" Level="<<Level <<" J="<<J << " x="<<xhead <<" w"<<basisFunction(xhead) << endl;
			return basisFunction(xhead) * dimnext(xtail...) + nextCoordinate(xhead, xtail...); 
		}
		template<typename F, typename... X>
		void eval(F f, X... x) {
			double xtail = basisFunction.center;
			//cout << "Dim="<<Dim <<" Level="<<Level <<" J="<<J << " x[i]="<<xtail << endl;
			dimnext.eval(f, x..., xtail);
			nextCoordinate.eval(f, x...);
		}
		template<typename LowerLevel, typename... X>
		void convert_to_surplus(LowerLevel& lower_level, X... x) {
			double xtail = basisFunction.center;
			dimnext.convert_to_surplus(lower_level, x..., xtail);
			nextCoordinate.convert_to_surplus(lower_level, x...);
		}
	};
	template<int Count>
	class Coordinate<Count,Count> {
	public:
		double operator()(double x, double* coordinates) {
			return 0;
		}
		double operator()(double x) {
			return 0;
		}
		template<typename F, typename... X>
		void eval(F f, X... x) {
		}
		template<typename... X>
		double operator()(double xhead, X... xtail) {
			return 0; 
		}
		template<typename LowerLevel, typename... X>
		void convert_to_surplus(LowerLevel& lower_level, X... x) {
		}
	};
	Coordinate<> coordinates;
	//double coordinate[size];


	void print() {
		cout << "I = " << I << " ";
		//dimnext.print();
		cout << endl;
		next.print();
	}
	template<typename F, typename... X>
	void eval(F f, X... x) {
		coordinates.eval(f, x...);
		next.eval(f, x...);
	}
	template<typename... X>
	double operator()(X... x) {
		return coordinates(x...) + next(x...);
	}
	template<typename LowerLevel, typename... X>
	void convert_to_surplus(LowerLevel& lower_level, X... x) {
		coordinates.convert_to_surplus(lower_level, x...);
		next.convert_to_surplus(lower_level, x...);
	}
};
template<int Dim, template<int> class BasisSetTemplate, int Norm, int I>
class Product<Dim, BasisSetTemplate, Norm, -1, I>  {
public:
	void print() {
		cout << "end" << endl;
	}
	template<typename F, typename... X>
	void eval(F f, X... x) {
	}
	template<typename... X>
	double operator()(X... x) {
		return 0;
	}
	template<typename LowerLevel, typename... X>
	void convert_to_surplus(LowerLevel& lower_level, X... x) {
	}
};
template<template<int> class BasisSetTemplate, int Norm, int NormLeft>
class Product<1, BasisSetTemplate, Norm, NormLeft, 0> {
public:
	//Vector<BasisSetTemplate, Norm-I> vector;
	static const int Level = NormLeft;
	//Vector<BasisSetTemplate, Level> vector;

	typedef BasisSetTemplate<Level> BasisSet;
	static const int size = BasisSet::size;
	template<int J=0, int Count=BasisSet::size>
	class Coordinate {
	public:
		typedef Coordinate<J+1, Count> NextCoordinate;
		typedef typename BasisSet::template BasisFunction<J> BasisFunction;
		NextCoordinate nextCoordinate;
		BasisFunction basisFunction;
		double coordinate;
		/*double operator()(double x, double* coefficients) {
			return basisFunction(x) * coefficients[0] + nextCoordinate(x, &coefficients[1]); 
		}*/
		double operator()(double x) {
			//cout << "Dim="<<1 <<" Level="<<Level <<" J="<<J << " x="<<x <<" w"<<basisFunction(x) <<" coordinate="<<coordinate << endl;
			return coordinate * basisFunction(x) + nextCoordinate(x); 
		}
		template<typename F, typename... X>
		void eval(F f, X... x) {
			double xtail = basisFunction.center;
			coordinate = f(x..., xtail);
			//cout << "Dim="<<1 <<" Level="<<Level <<" J="<<J << " x[i]="<<xtail <<" coordinate="<<coordinate <<" this="<<this <<" data=" <<(&coordinate) <<endl;
			nextCoordinate.eval(f, x...);
		}
		template<typename LowerLevel, typename... X>
		void convert_to_surplus(LowerLevel& lower_level, X... x) {
			double xtail = basisFunction.center;
			double value = lower_level(Norm-1, x..., xtail);
			//cout << "surplus Dim="<<1 <<" Level="<<Level <<" J="<<J << " x[i]="<<xtail <<" coordinate="<<coordinate <<" value="<<value <<" surplus=" <<(coordinate-value) <<endl;
			coordinate -= value;
			nextCoordinate.convert_to_surplus(lower_level, x...);
		}
	};
	template<int Count>
	class Coordinate<Count,Count> {
	public:
		double operator()(double x, double* coordinates) {
			return 0;
		}
		double operator()(double x) {
			return 0;
		}
		template<typename F, typename... X>
		void eval(F f, X... x) {
		}
		template<typename LowerLevel, typename... X>
		void convert_to_surplus(LowerLevel& lower_level, X... x) {
		}
	};
	Coordinate<0, BasisSet::size> coordinates;


	void print() {
		cout << "Left = " << (NormLeft) << " " << endl;
	}
	template<typename F, typename... X>
	void eval(F f, X... x) {
		//cout << "this=" <<this << " sizeof=" << sizeof(Coordinate<0, BasisSet::size>) <<endl;
		//cout << " p1=" <<((void*)&coordinates) << " p2=" <<((void*)&coordinates.nextCoordinate) << endl;  
		coordinates.eval(f, x...);
		//next.eval(f, x...);
	}
	template<typename... X>
	double operator()(X... x) {
		return coordinates(x...);
	}
	template<typename LowerLevel, typename... X>
	void convert_to_surplus(LowerLevel& lower_level, X... x) {
		coordinates.convert_to_surplus(lower_level, x...);
	}
};


template<int Dim, template<int> class BasisSetTemplate, int MaxLevel>
class SparseGrid {
public:
	template<int I=0, int Count=MaxLevel+1>
	class SparseGridLevel {
	public:
		Product<Dim, BasisSetTemplate, I> product;
		SparseGridLevel<I+1,Count> next;
		template<typename F, typename... X>
		void eval(F f, X... x) {
			product.eval(f, x...);
			next.eval(f, x...);
			if(I == 0) {
				next.convert_to_surplus(*this);
			}
		}
		template<typename... X>
		double operator()(X... x) {
			return product(x...) + next(x...);
		}
		template<typename... X>
		double operator()(int max_level, X... x) {
			//cout << "op(..) ";
			double y = product(x...);
			//cout << "max_level=" <<max_level;
			//cout_many(x...);
			//cout << " y=" << y;
			//cout <<"I=" <<I;
			if(I < max_level)
			{
				//cout << " next...";
				y += next(max_level, x...);
			}
			//cout << endl;
			return y;
		}
		template<typename... X>
		double test(X... x) {
			return product(x...);// + next(x...);
		}
		template<typename LowerLevel>
		void convert_to_surplus(LowerLevel& top_level) {
			//cout << "convert to surplus" <<I << endl;
			product.convert_to_surplus(top_level);
			next.convert_to_surplus(top_level);
		}
	};
	template<int Count>
	class SparseGridLevel<Count, Count> {
	public:
		template<typename F, typename... X>
		void eval(F f, X... x) {
		}
		template<typename... X>
		double operator()(X... x) {
			return 0;
		}
		template<typename LowerLevel>
		void convert_to_surplus(LowerLevel& lower_level) {
		}
	};
	SparseGridLevel<> levels;
	template<typename F, typename... X>
	void eval(F f, X... x) {
		levels.eval(f, x...);
		//next.eval(f, x...);
	}
	template<typename... X>
	double operator()(X... x) {
		return levels(x...);
	}
	void test2d() {
		auto f = [](double x, double y) {
			return x;
		};
		levels.eval(f);
	}
	void test1d() {
		auto f = [](double x) {
			return x;
		};
		levels.eval(f);
	}
};

/*
template<int Level>
class BasisTriHierLevel;

template<int Level>
class BasisTriHierLevel {
public:
	static const int dof = (1<<(Level-1)); //+1; 
	static const int total_dof = dof + BasisTriHierLevel<Level-1>::total_dof;
	typedef BasisTriHierLevel<Level-1> sub;
	//enum { size = (2<<Level)+1; }
	double x(int i) {
		const double hl = pow(2, -(Level-1));
		return hl * (i+0.5);
	}
	double operator()(int i, double x) {
		const double hl = pow(2, -Level);
		int j = i * 2 + 1;
		if((x > 1) || (x < 0)) // outside range?
			return 0;
		else {
			return x > (j+1)*hl ? 0 : (x < (j-1)*hl ? 0 : 1-fabs(x/hl-j));
		}
	}
};

template<>
class BasisTriHierLevel<0> {
public:
	static const int dof = 2; 
	static const int total_dof = dof; 
	//enum { size = (2<<Level)+1; }
	double operator()(int i, double x) {
		if((x > 1) || (x < 0)) // outside range or odd, skip
			return 0;
		else {
			if(i == 0)
				return x > 1 ? 0 : (x < 0 ? 0 : 1-x);
			else if(i == 1)
				return x > 1 ? 0 : (x < 0 ? 0 : x);
			else
				return 0;
		}
	}
	double x(int i) {
		return i == 0 ? 0 : 1;
	}
};

template<int Level>
class BasisTriHier {
};

namespace internal {
	template<int... Indices>
	class sum;


	template<int  I, int ... Indices>
	class sum<I, Indices...> {
		//enum { value = I+(sum<Indices...>::value) };
		static const int value = I + sum<Indices...>::value;
	};
	template<int I>
	class sum<I> {
		//enum { value = 0 };
		static const int value = I;
	};

	template<int... Indices>
	class MultiIndex {
		enum { norm = sum<Indices...>::value };
	};


}

template<int Dim, template<int> class Basis, int Level, int CurrentLevel, int I=Basis<CurrentLevel>::dof-1>
class SparseGridLevelFunction;
template<int Dim, template<int> class Basis, int Level, int CurrentLevel, int J=CurrentLevel>
class SparseGridLevel;


// DIM == 1

template<template<int> class Basis, int Level, int CurrentLevel, int I>
class SparseGridLevelFunction<1, Basis, Level, CurrentLevel, I> {
public:
	//static const int dof = Basis<Level, I>::dof;
	SparseGridLevelFunction<1, Basis, Level, CurrentLevel, I-1> nextfunction;
	double a;
	Basis<CurrentLevel> basis;
	template<class F, typename... X>
	void eval(F f, X... xtail) {
		nextfunction.eval(f, xtail...);
		double x = basis.x(I);
		double y = f(xtail..., x);
		a = y;// - g(xtail..., x);
		//a[i] = y - subgrid(x); // surplus
		cout << "  at level: " << CurrentLevel << " i = " << I << " y = " << y << " surplus = " << a << " x = ";
		cout_many(xtail..., x);
		cout << endl;
	}
	double operator()(double x) {
		double w = basis(I, x);
		return w * a + nextfunction(x);
	};
};
template<template<int> class Basis, int Level, int CurrentLevel>
class SparseGridLevelFunction<1, Basis, Level, CurrentLevel, -1> {
public:
	template<class F,  typename... X>
	void eval(F f, X... xtail) {
		return;
	}
	double operator()(double x) { return 0; }
};

template<int Dim, template<int> class Basis, int Level, int CurrentLevel, int I>
class SparseGridLevelFunction {
public:
	//static const int dof = Basis<Level, I>::dof;
	SparseGridLevelFunction<Dim, Basis, Level, CurrentLevel, I-1> nextfunction;
	SparseGridLevel<Dim-1, Basis, Level, Level> otherdimension;
	Basis<CurrentLevel> basis;
	template<class F>
	void eval(F f) {
		//sublevel.eval(f, g);
		double x1 = basis.x(I);
		cout << "at level: " << CurrentLevel << " i = " << I << " xi = " << x1 << endl;
		otherdimension.eval(f,x1);
		nextfunction.eval(f);
		//double x1 = basis.x(I);
		//double y = f(x1, x2);
		//a = y - g(x1, x2);
		//a[i] = y - subgrid(x); // surplus
		//cout << "at level: " << Level << " i = " << I << " y = " << y << " surplus = " << a << " x1 = " << x1 << " x2 = " << x2 << endl;
	}
	double operator()(double x1, double x2) {
		double w = basis(I, x1);
		return w * (otherdimension(x2))  + nextfunction(x1, x2);
	};
};
template<int Dim, template<int> class Basis, int Level, int CurrentLevel>
class SparseGridLevelFunction<Dim, Basis, Level, CurrentLevel, -1> {
public:
	//static const int dof = Basis<Level, I>::dof;
	double a;
	//Basis<Level> basis;
	template<class F, typename... X>
	void eval(F f, X... x) { return; }
	double operator()(double x, double y) { return 0; }
};


template<template<int> class Basis, int Level, int CurrentLevel, int J>
class SparseGridLevel<1, Basis, Level, CurrentLevel, J> {
public:
	//static const int dof = Basis<Level, I>::dof;
	SparseGridLevelFunction<1, Basis, Level, CurrentLevel> levelfunction;
	SparseGridLevel<1, Basis, Level, CurrentLevel-1> lowerlevel;
	double a;
	//Basis<Level> basis;
	template<class F, typename... X>
	void eval(F f, X... x) {
		lowerlevel.eval(f, x...);
		levelfunction.eval(f, x...);
	}
	double operator()(double x) {
		return lowerlevel(x) + levelfunction(x);
	}

};
template<template<int> class Basis, int Level, int CurrentLevel>
class SparseGridLevel<1, Basis, Level, CurrentLevel, -1> {
public:
	//static const int dof = Basis<Level, I>::dof;
	double a;
	//Basis<Level> basis;
	template<class F, typename... X>
	void eval(F f, X... x) {
		return;
	}
	double operator()(double x) { return 0; }
};

template<int Dim, template<int> class Basis, int Level, int CurrentLevel, int J>
class SparseGridLevel {
public:
	//static const int dof = Basis<Level, I>::dof;
	SparseGridLevelFunction<Dim, Basis, Level, CurrentLevel> levelfunction;
	SparseGridLevel<Dim, Basis, Level, CurrentLevel-1> lowerlevel;
	double a;
	//Basis<Level> basis;
	template<class F, typename... X>
	void eval(F f, X... x) {
		lowerlevel.eval(f, x...);
		levelfunction.eval(f, x...);
	}
	template<typename... X>
	double operator()(X... x) {
		return lowerlevel(x...) + levelfunction(x...);
	}

};
template<int Dim, template<int> class Basis, int Level, int CurrentLevel>
class SparseGridLevel<Dim, Basis, Level, CurrentLevel, -1> {
public:
	//static const int dof = Basis<Level, I>::dof;
	double a;
	Basis<Level> basis;
	template<class F, typename... X>
	void eval(F f, X... x) {
		return;
	}
	template<typename... X>
	double operator()(X... x) { return 0; }
};


template<int DIM, template<int> class Basis, int Level>
class SparseGridLevels<1, Basis, Level> {
	SparseGridLevel<1, Basis, Level-1> sub;
	Basis<Level, I> basis;
};

template<int DIM, template<int> class Basis, int Level>
class SparseGrid;

template< template<int> class Basis, int Level>
class SparseGrid<1, Basis, Level> {
public:
	SparseGridLevel<1, Basis, Level, Level> gridlevels;
	void test() {
		auto f = [](double x) { return x*x; };
		//auto g = [](double x) { return 0; };
		gridlevels.eval(f);
	}
	double operator()(double x) {
		return gridlevels(x);
		
		return 0;
	};
};
template< template<int> class Basis, int Level>
class SparseGrid<2, Basis, Level> {
public:
	SparseGridLevel<2, Basis, Level, Level> gridlevels;
	void test() {
		//auto f = [](double x, double y) { return sin(x) + cos(y); };
		auto f = [](double x, double y) { return x*x; };
		//auto g = [](double x, double y) { return 0; };
		gridlevels.eval(f);
		//gridlevels.to_
		//gridlevels.eval(f);
	}
	double operator()(double x1, double x2) {
		return gridlevels(x1, x2);
		
		return 0;
	};
};
*/
/*
template< template<int> class Basis>
class SparseGrid<1, Basis, 0> {
public:
	Basis<0> basis;
	static const int dof = Basis<0>::dof;
	double a[dof];
	double operator()(double x) {
		double y = 0;
		for(int i = 0; i < dof; i++) {
			double w = basis(i, x);
			y += w * a[i];
		}
		return y;
	};
	template<class F>
	void eval(F f) {
		for(int i = 0; i < dof; i++) {
			double x = basis.x(i);
			a[i] = f(x);
			cout << "at level: " << 0 << " i = " << i << " a = " << a[i] << " x = " << x << endl;
		}
	}
};

template<int Dim, template<int> class Basis, int Level>
class SparseGrid2 {
public:

};
*/
/*
template<int Dim, template<int> class Basis, int Level>
class SparseGrid {
public:
	SparseGrid<Dim, Basis, Level-1> subgrid;
	SparseGrid<Dim-1, Basis, Level> subgriddim;
	static const int dof = Basis<Level>::dof;
	Basis<Level> basis;

	void tes2t() {
		auto f = [](double x, double y) { return x*x + y*y ; };
		eval(f);
	}
	template<class F>
	void eval(F f) {
		//auto g[]
		subgrid.eval(f);
		for(int i = 0; i < dof; i++) {
			double x = basis.x(i);
			double y = f(x);
			a[i] = y - subgrid(x); // surplus
			cout << "at level: " << Level << " i = " << i << " y = " << y << " surplus = " << a[i] << " x = " << x << endl;
		}
	}
	template<typename... Xtail>
	double operator()(double x, Xtail... xtail) {
		double y = subgrid(x);
		for(int i = 0; i < dof; i++) {
			double w = basis(i, x);
			y += w * subgriddim(xtail...);
		}
		return y;
	};
};

template<int Dim, template<int> class Basis, int Level>
class SparseGrid<Dim, Basis, 0> {
public:
	//SparseGrid<Dim, Basis, Level-1> subgrid;
	SparseGrid<Dim-1, Basis, Level> subgriddim;
	static const int dof = Basis<Level>::dof;
	Basis<Level> basis;

	void tes2t() {
		auto f = [](double x, double y) { return x*x + y*y ; };
		eval(f);
	}
	template<class F>
	void eval(F f) {
		//auto g[]
		subgrid.eval(f);
		for(int i = 0; i < dof; i++) {
			double x = basis.x(i);
			double y = f(x);
			a[i] = y - subgrid(x); // surplus
			cout << "at level: " << Level << " i = " << i << " y = " << y << " surplus = " << a[i] << " x = " << x << endl;
		}
	}
	template<typename... Xtail>
	double operator()(double x, Xtail... xtail) {
		double y = subgrid(x);
		for(int i = 0; i < dof; i++) {
			double w = basis(i, x);
			y += w * subgriddim(xtail...);
		}
		return y;
	};
};
*/

/*
template<int DIM, class Basis>
class SparseGrid {
};

template<>
class SparseGrid<1> {
public:
	double x1, x2;
	SparseGrid(double x1, double x2) : x1(x1), x2(x2) {
	}
};

template<>
class SparseGrid<0> {
};

template<int DIM, int Level>
class SparseGrid {
	//typename Grid<DIM-1,LevelTail...> subgrid;
	
	template<int point_dim>
	double operator()(int j) {
		if(point_dim == DIM) {
			return pow(2, -LevelHead) * j;
		} else {
			return subgrid<point_dim>(j);
		}
	}
	template<int point_dim>
	double operator()(double x, int j) {
		double hl = pow(2, -LevelHead);
		if(point_dim == DIM) {
			return 1 - abs(x/hl - j);
		} else {
			return subgrid<point_dim>(j);
		}
	}
};
template<int DIM, int Levels, int MaxLevels>
class SparseGridLevel {
	template<int point_dim>
	double operator()(int j) {
		if(point_dim == DIM) {
			return pow(2, -LevelHead) * j;
		} else {
			return subgrid<point_dim>(j);
		}
	}
	template<int point_dim>
	double operator()(double x, int j) {
		double hl = pow(2, -LevelHead);
		if(point_dim == DIM) {
			return 1 - abs(x/hl - j);
		} else {
			return subgrid<point_dim>(j);
		}
	}
};


class GridBasisSpace {
	template<typename... Ts>
	operator()(double xi, Ts... ) {
		return grid.basis(xi, 
	}
}
*/
}