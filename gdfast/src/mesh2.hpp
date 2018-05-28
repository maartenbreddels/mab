#pragma once
#include <Eigen/Core>
#include "datatypes.hpp"
#include <Eigen/Eigen>
#include "profile.hpp"
#include "integration.hpp"
#include "transformation.hpp"

namespace gd {

using namespace std;
USING_PART_OF_NAMESPACE_EIGEN

class Grid1d {
public:
	virtual int findindex(double) { return 0; }
	virtual double indexto_x(int) { return 0; }
	//virtual void indexto_x(int) { return 0; }
	virtual bool inrange(double) { return false; }
	virtual int length() {return 0;}
};

template<class Base=Grid1d>
class Grid1dRegular : public Base {
public:
	Grid1dRegular(double x1, double x2, int gridpoints) : x1(x1), x2(x2), _gridpoints(gridpoints) {}
	int findindex(double x) {
		return x < x1 ? 0 :
			(x >= x2 ? _gridpoints-1 : (int)((x-x1)/(x2-x1)*(_gridpoints-1)));
	}
	double indexto_x(int i) {
		return x1 + (x2-x1)/(length()-1)*i;
	}
	bool inrange(double r) {
		return (r >= x1) && (r < x2); 
	}
	int length() { return _gridpoints; }
	double x1, x2;
	int _gridpoints;
};

/*class Grid2d {
public:
	virtual int findindex(double, double) { return 0; }
	virtual bool inrange(double, double) { return false; }
	virtual int length() {return 0;}
};*/

template<class GridX=Grid1d, class GridY=Grid1d>
class Grid2d {
public:
	typedef GridX GridXType;
	typedef GridY GridYType;
	Grid2d(GridX* gridx, GridY* gridy) : gridx(gridx), gridy(gridy) {}
	void findindex2d(double x, double y, int& xi, int &yi) {
		//double r = sqrt(x*x+y*y);
		//double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI); 
		xi = gridx->findindex(x);
		yi = gridy->findindex(y);
	} 
	int findindex(double x, double y) {
		//double r = sqrt(x*x+y*y);
		//double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI); 
		return gridx->findindex(x) * gridy->length() + gridy->findindex(y);;
	} 
	void indexto_xy(int i, int j, double &x, double &y) {
		x = gridx->indexto_x(i);
		y = gridy->indexto_x(j);
	}
	bool inrange(double x, double y) {
		//double r = sqrt(x*x+y*y);
		//double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI);
		return gridx->inrange(x) && gridy->inrange(y);
	}
	int length() {
		return gridx->length() * gridy->length();
	}
	GridX* gridx;
	GridY* gridy;
};

template<class GridX=Grid1d, class GridY=Grid1d, class GridZ=Grid1d>
class Grid3d {
public:
	typedef GridX GridXType;
	typedef GridY GridYType;
	typedef GridZ GridZType;
	Grid3d(GridX* gridx, GridY* gridy, GridZ* gridz) : gridx(gridx), gridy(gridy), gridz(gridz) {}
	void findindex3d(double x, double y, double z, int& xi, int &yi, int &zi) {
		//double r = sqrt(x*x+y*y);
		//double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI); 
		xi = gridx->findindex(x);
		yi = gridy->findindex(y);
		zi = gridz->findindex(z);
	} 
	int findindex(double x, double y, double z) {
		//double r = sqrt(x*x+y*y);
		//double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI); 
		return (gridx->findindex(x) * gridy->length() + gridy->findindex(y)) * gridz->length() + gridz->findindex(z);
	} 
	void indexto_xyz(int i, int j, int k, double &x, double &y, double &z) {
		x = gridx->indexto_x(i);
		y = gridy->indexto_x(j);
		z = gridz->indexto_x(k);
	}
	bool inrange(double x, double y, double z) {
		//double r = sqrt(x*x+y*y);
		//double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI);
		return gridx->inrange(x) && gridy->inrange(y) && gridz->inrange(z);
	}
	int length() {
		return gridx->length() * gridy->length() * gridz->length();
	}
	GridX* gridx;
	GridY* gridy;
	GridZ* gridz;
};

/*
		basis_evaluator<2, Basis::degree, Basis> basis_evaluator;
		if(mesh.inrange(x))
			value = basis_evaluator.eval(solution, this, xi, u, yi, v);
*/

template<int DIM, class B, int I=B::degree, class BI=B, class T=double>
struct basis_evaluator;

template<int DIM, class B, int I, class BI, class T>
struct basis_evaluator
{
	typedef basis_evaluator<DIM, B, I-1, typename BI::next_type, T> next_typeI;
	typedef basis_evaluator<DIM-1, B, B::degree, B, T> next_typeDIM;
	next_typeI nextI;
	next_typeDIM nextDIM;

	template<class... Ts>
	T eval(int i, double u, Ts... ts) {
		if((i % (B::degree+1)) == I) { // first search for right i
			BI basis;
			return basis(u) * nextDIM.eval(i/(B::degree+1), ts...);
		} else {
			return nextI.eval(i, u, ts...); 
		}
	}
};

template<int DIM, class B,  class BI, class T>
struct basis_evaluator<DIM, B, -1, BI, T>
{
	template<class... Ts>
	T eval(int, double, Ts...) {
		return 1;
	}
};

template<class B, int I, class BI, class T>
struct basis_evaluator<0, B, I, BI, T>
{
	T eval(int) {
		return 1;
	}
};


template<int DIM, int I, class B, class BI, class T>
struct basis_function_integrator;

//template<int DIM, int I, class B, class BI, class T>
//struct basis_function_integrator<2,I,B,BI,T>
template<int DIM, int I, class B, class BI, class T>
struct basis_function_integrator
{
	typedef basis_function_integrator<DIM,   I-1, B, typename BI::next_type, T> next_typeI;
	typedef basis_function_integrator<DIM-1, B::degree, B, B, T> next_typeDIM;
	T x1, x2;
	next_typeI nextI;
	next_typeDIM nextDIM;

	template<class... Ts>
	basis_function_integrator(T x1, T x2, Ts... ts) : x1(x1), x2(x2), nextI(x1, x2, ts...), nextDIM(ts...) {
		//cout << "DIM " << DIM << " " << x1 << " " << x2 << endl;
	} 
	
	template<class F>
	T integrate_function2d(F f, int i) {
		if((i % (B::degree+1)) == I) { // first search for right i
			BI basis;
			double dx = (x2-x1);
			auto fi = [&](double x) -> double { return basis((x-this->x1)/dx) / dx  * this->nextDIM.integrate_function2d_y(f, i / (B::degree+1), x); };
			IntegratorGSL<> integratorGSL(fi);
			//cout << "integrate: " << x1 << " to " << x2 << endl;
			return integratorGSL.integrate(x1, x2);
		} else {
			return nextI.integrate_function2d(f, i);
		}
	}
	template<class F>
	T integrate_function2d_y(F f, int i, double x) {
		if((i % (B::degree+1)) == I) { // first search for right i
			BI basis;
			double dx = (x2-x1);
			auto fi = [&](double y) -> double { return basis((y-this->x1)/dx)/dx * f(x, y); };
			IntegratorGSL<> integratorGSL(fi);
			//cout << "integrate: " << x1 << " to " << x2 << endl;
			return integratorGSL.integrate(x1, x2);
		} else {
			return nextI.integrate_function2d_y(f, i, x);
		}
	}
};

template<int DIM, class B, class BI, class T>
struct basis_function_integrator<DIM, -1, B, BI, T>
{
	template<class... Ts>
	basis_function_integrator(T, T, Ts...) {} 
	template<class F>
	T integrate_function2d(F, int) { return 1; }
	template<class F>
	T integrate_function2d_y(F, int, double) { return 1; }
};

template<int I, class B, class BI, class T>
struct basis_function_integrator<0, I, B, BI, T>
{
	template<class F>
	T integrate_function(F, int) { return 1; }
};


template<int DIM, int I, int J, class B, class BI, class BJ, class T=double>
struct basis_integrator;

template<int DIM, int I, int J, class B, class BI, class BJ, class T>
struct basis_integrator
{
	typedef basis_integrator<DIM, I-1, J, B, typename BI::next_type, BJ, T> next_typeI;
	typedef basis_integrator<DIM, I, J-1, B, BI, typename BJ::next_type, T> next_typeJ;
	typedef basis_integrator<DIM-1, B::degree, B::degree, B, B, B, T> next_typeDIM;
	next_typeI nextI;
	next_typeJ nextJ;
	next_typeDIM nextDIM;


	//template<class... Ts>
	T integrate(int i, int j) {
		if((i % (B::degree+1)) == I) { // first search for right i
			if((j  % (B::degree+1)) == J) { // then right j
				BI basis1;
				BJ basis2;
				auto f = [&](double x) ->double { return basis1(x) * basis2(x) * this->nextDIM.integrate(i / (B::degree+1), j / (B::degree+1)); };
				IntegratorGSL<> integratorGSL(f);
				return integratorGSL.integrate(0, 1);
			}  else {
				return nextJ.integrate(i, j);
			}
		} else {
			return nextI.integrate(i, j);
		}
	}
};

// sentinels
template<int I, int J, class B, class BI, class BJ, class T>
struct basis_integrator<0, I, J, B, BI, BJ, T>
{
	template<class... Ts>
	T integrate(int, int) {
		return 1;
	}
};

template<int DIM, int J, class B, class BI, class BJ, class T>
struct basis_integrator<DIM, -1, J, B, BI, BJ, T>
{
	template<class... Ts>
	T integrate(int, int) {
		return 1;
	}
};

template<int DIM, int I, class B, class BI, class BJ, class T>
struct basis_integrator<DIM, I, -1, B, BI, BJ, T>
{
	template<class... Ts>
	T integrate(int, int) {
		return 1;
	}
};




template<int DIM, class Basis, class T=double>
struct MeshRegularNodalHelper;

template<class Basis, class T>
struct MeshRegularNodalHelper<0, Basis, T> {
	typedef Basis basis_type;
	enum { dof_per_cell = 1 };
};


template<int DIM, class Basis, class T>
struct MeshRegularNodalHelper {
	typedef Basis basis_type;
	typedef MeshRegularNodalHelper<DIM, Basis, T> type;
	typedef MeshRegularNodalHelper<DIM-1, Basis, T> sub_type;
	enum { dof_per_cell = (Basis::degree+1) * sub_type::dof_per_cell };

	/*int get_dof() { return dof;}
	int get_n_cells() { return n_cells;}
	int dof_index(int cell_index, int local_index) {
		return  cell_index*(dof_per_cell-1)+local_index;
	}*/
};





template<int DIM, class Basis, class T=double>
class MeshRegularNodal;

template<class Basis, class T>
class MeshRegularNodal<2, Basis, T> {
public:
	typedef Basis basis_type;
	T x1, x2;
	int n_cells_x;
	int n_cells_y;
	Grid1dRegular<> xgrid;
	Grid1dRegular<> ygrid;
	Grid2d<Grid1dRegular<>, Grid1dRegular<>> grid;
	MatrixXd M;
	//Transformation1d_in_3d* transformation;
	int dof, dofx, dofy;
	int dof1d;
	enum { dof_per_cell = MeshRegularNodalHelper<2, Basis, T>::dof_per_cell };
	enum { dof_per_cell1d = MeshRegularNodalHelper<1, Basis, T>::dof_per_cell };

	MeshRegularNodal(T x1, T y1, T x2, T y2, int n_cells_x, int n_cells_y) : n_cells_x(n_cells_x), n_cells_y(n_cells_y), xgrid(x1, x2, n_cells_x+1), ygrid(y1, y2, n_cells_y+1), grid(&xgrid, &ygrid), M(1, 1) {
		
		if(Basis::degree == 0) {
			dofx = n_cells_x;
			dofy = n_cells_y;
		} else {
			dofx = (1 + n_cells_x) + (dof_per_cell1d-2)*n_cells_x; // 1 dof per border + dofs inside the cel
			dofy = (1 + n_cells_y) + (dof_per_cell1d-2)*n_cells_y; // 1 dof per border + dofs inside the cel
		}
		dof = dofx * dofy;
		//cout << "n_cells_x = " << n_cells_x << " dofx = " << dofx << endl; 
		//cout << "n_cells_y = " << n_cells_y << " dofy = " << dofy << endl; 
		//cout << "dof = " << dof << " dof_per_cell = " << dof_per_cell << endl;
		//MatrixXd m = MatrixXd::Zero(dof, dof);
		M.resize(dof, dof);
		M =  MatrixXd::Zero(dof, dof);
		//T scale = 1; //TODO: (x2-x1)/n_cells;

		//basis_integrator<1, Basis::degree, Basis::degree, Basis, Basis, Basis, T> bi;
		//cout  << "test 00 " << bi.integrate(0, 0) << endl;
		//cout  << "test 01 " << bi.integrate(0, 1) << endl;
		//cout  << "test 11 " << bi.integrate(1, 1) << endl;

		T integrals[dof_per_cell][dof_per_cell];
		for(int j = 0; j < dof_per_cell; j++) {
			for(int k = 0; k < (j+1); k++) {
				basis_integrator<2, Basis::degree, Basis::degree, Basis, Basis, Basis, T> bi;
				T integral = bi.integrate(j, k);
				//typedef selfintegrator<Basis::degree, Basis::degree> selfintegrator_type;
				//selfintegrator_type si; 
				//double integral = si.integrate(j,k);
				//cout << j << " " << k << " " << integral << endl;
				integrals[j][k] = integral;
				integrals[k][j] = integral;
			}
		}
		for(int xi = 0; xi < n_cells_x; xi++) {
			for(int yi = 0; yi < n_cells_y; yi++) {
				for(int j = 0; j < dof_per_cell; j++) {
				for(int k = 0; k < dof_per_cell; k++) {
					int i1 = this->dof_index(xi, yi, j);
					int i2 = this->dof_index(xi, yi, k);
					//cout << "xi = " << xi << " yi = " << yi << " j = " << j << " k = " << k;
					//cout << "            i1 = " << i1 << " i2 = " << i2 << endl;
					M(i1, i2) += integrals[j][k];
				}}
			}
		}
		/*for(int i = 0; i < n_cells; i++) {
			for(int j = 0; j < dof_per_cell; j++) {
				for(int k = 0; k < dof_per_cell; k++) {
					int i1 = i*(dof_per_cell-1)+j;
					int i2 = i*(dof_per_cell-1)+k;
					m(i1, i2) = m(i1, i2) + integrals[j][k];
				}
			}
		}*/
		//cout << M << endl;
		//m = m * scale;
		//ctrans = m.inverse();
	}

	double eval_(double_vector solution_, double x, double y) {
		VectorXd solution = VectorXd::Map(solution_.data().begin(), solution_.size());
		return eval(solution, x, y);
	}

	double eval(VectorXd& solution, double x, double y) {
		/*int cell_index = mesh.findindex(x);
		T xleft = mesh.indexto_x(cell_index);
		T xright = mesh.indexto_x(cell_index+1);
		*/
		int xi = xgrid.findindex(x);
		int yi = ygrid.findindex(y);
		T x1 = xgrid.indexto_x(xi);
		T x2 = xgrid.indexto_x(xi+1);
		T y1 = ygrid.indexto_x(yi);
		T y2 = ygrid.indexto_x(yi+1);
		double u = (x-x1)/(x2-x1);
		double v = (y-y1)/(y2-y1);
		//grid.index
		//T dx = xright-xleft;
		double value = 0;
		//cout << "xi = " << xi;
		//cout << " yi = " << yi;
		if(grid.inrange(x, y)) {
			basis_evaluator<2, Basis> basis_evaluator;
			for(int i = 0; i < dof_per_cell; i++) {
				//cout << " i = " << i;
				//cout << " index = " << dof_index(xi, yi, i) << endl;
				double a = solution(dof_index(xi, yi, i));
				value += a * basis_evaluator.eval(i, u, v);
			}
		} else {
			cout << "x and y not in range: (" << x << "," << y << ")" << endl;
		}
		return value;
	}

	double basis_uv(int i, double u, double v) {
		basis_evaluator<2, Basis> basis_evaluator;
		return basis_evaluator.eval(i, u, v);
	}

	void solve_coordinates(double_vector inner_products, double_vector coordinates) {
		assert((int)inner_products.size() == dof);
		assert((int)coordinates.size() == dof);
		VectorXd x = VectorXd::Map(inner_products.data().begin(), inner_products.size());
		VectorXd a(dof);
		M.llt().solve(x, &a);
		VectorXd::Map(coordinates.data().begin(), coordinates.size()) = a;
	}

	void test(double_vector result_, double scale1, double scale2) {
		auto f = [&](double x, double y) -> double { return cos(x * scale1 + y * scale2); };
		VectorXd x = VectorXd::Zero(dof);
		assert((int)result_.size() == dof);
		for(int yi = 0; yi < n_cells_y; yi++) {
			for(int xi = 0; xi < n_cells_x; xi++) {
				for(int i = 0; i < dof_per_cell; i++) {
					double x1, x2;
					double y1, y2;
					grid.indexto_xy(xi, yi, x1, y1); 
					grid.indexto_xy(xi+1, yi+1, x2, y2);
					//cout << "integrate x[" << x1 << " " << x2 << "] y[" << y1 << " " << y2 << "]";
					basis_function_integrator<2, Basis::degree, Basis, Basis, T> bfi(x1, x2, y1, y2);
					double a = bfi.integrate_function2d(f, i);
					x(dof_index(xi, yi, i)) += a;
					//cout << " " << a;
 					//cout << endl;
				}
			}
		}
		VectorXd a(dof);
		M.llt().solve(x, &a);
		//VectorXd result = VectorXd::Map(result_.data().begin(), result_.size());
		VectorXd::Map(result_.data().begin(), result_.size()) = a;
		//result = a;
	}

	//int get_dof() { return dof; }
	int dof_index(int x_index, int y_index, int basis_index) {
		if(Basis::degree == 0) {
			return x_index + y_index * dofx;
		}
		int xi = x_index*(dof_per_cell1d-1) + basis_index % dof_per_cell1d;
		int yi = y_index*(dof_per_cell1d-1) + basis_index / dof_per_cell1d;
		return xi + yi * dofx;
	}
	int get_dof() { return dof; }
	int get_dof_per_cell() { return dof_per_cell; } 
	/*int get_n_cells() { return grid->;}
	int dof_index(int cell_index, int local_index) {
		return  cell_index*(dof_per_cell-1)+local_index;
	}*/

}; // class MeshRegularNodal<2...>


template<class Basis, class T>
class MeshRegularNodal<3, Basis, T> {
	public:
		typedef Basis basis_type;
		T x1, x2;
		int n_cells_x;
		int n_cells_y;
		int n_cells_z;
		Grid1dRegular<> xgrid;
		Grid1dRegular<> ygrid;
		Grid1dRegular<> zgrid;
		Grid3d<Grid1dRegular<>, Grid1dRegular<>, Grid1dRegular<>> grid;
		MatrixXd M;
	//Transformation1d_in_3d* transformation;
		int dof, dofx, dofy, dofz;
		int dof1d;
		enum { dof_per_cell = MeshRegularNodalHelper<3, Basis, T>::dof_per_cell };
		enum { dof_per_cell1d = MeshRegularNodalHelper<1, Basis, T>::dof_per_cell };
		
	MeshRegularNodal(T x1, T y1, T z1, T x2, T y2, T z2, int n_cells_x, int n_cells_y, int n_cells_z) : n_cells_x(n_cells_x), n_cells_y(n_cells_y), n_cells_z(n_cells_z), xgrid(x1, x2, n_cells_x+1), ygrid(y1, y2, n_cells_y+1), zgrid(z1, z2, n_cells_z+1), grid(&xgrid, &ygrid, &zgrid), M(1, 1) {
		
		if(Basis::degree == 0) {
			dofx = n_cells_x;
			dofy = n_cells_y;
			dofz = n_cells_z;
		} else {
			dofx = (1 + n_cells_x) + (dof_per_cell1d-2)*n_cells_x; // 1 dof per border + dofs inside the cel
			dofy = (1 + n_cells_y) + (dof_per_cell1d-2)*n_cells_y; // 1 dof per border + dofs inside the cel
			dofz = (1 + n_cells_z) + (dof_per_cell1d-2)*n_cells_z; // 1 dof per border + dofs inside the cel
		}
		dof = dofx * dofy * dofz;
		cout << "n_cells_x = " << n_cells_x << " dofx = " << dofx << endl; 
		cout << "n_cells_y = " << n_cells_y << " dofy = " << dofy << endl; 
		cout << "n_cells_z = " << n_cells_z << " dofz = " << dofz << endl; 
		cout << "dof = " << dof << " dof_per_cell = " << dof_per_cell << endl;
		//MatrixXd m = MatrixXd::Zero(dof, dof);
		M.resize(dof, dof);
		M =  MatrixXd::Zero(dof, dof);
		//T scale = 1; //TODO: (x2-x1)/n_cells;
		
		//basis_integrator<1, Basis::degree, Basis::degree, Basis, Basis, Basis, T> bi;
		//cout  << "test 00 " << bi.integrate(0, 0) << endl;
		//cout  << "test 01 " << bi.integrate(0, 1) << endl;
		//cout  << "test 11 " << bi.integrate(1, 1) << endl;
		
		T integrals[dof_per_cell][dof_per_cell];
		for(int j = 0; j < dof_per_cell; j++) {
			for(int k = 0; k < (j+1); k++) {
				basis_integrator<2, Basis::degree, Basis::degree, Basis, Basis, Basis, T> bi;
				T integral = bi.integrate(j, k);
				//typedef selfintegrator<Basis::degree, Basis::degree> selfintegrator_type;
				//selfintegrator_type si; 
				//double integral = si.integrate(j,k);
				//cout << j << " " << k << " " << integral << endl;
				integrals[j][k] = integral;
				integrals[k][j] = integral;
			}
		}
		for(int xi = 0; xi < n_cells_x; xi++) {
			for(int yi = 0; yi < n_cells_y; yi++) {
				for(int zi = 0; zi < n_cells_z; zi++) {
					for(int j = 0; j < dof_per_cell; j++) {
						for(int k = 0; k < dof_per_cell; k++) {
							int i1 = this->dof_index(xi, yi, zi, j);
							int i2 = this->dof_index(xi, yi, zi, k);
							//cout << "xi = " << xi << " yi = " << yi << " zi = " << zi << " j = " << j << " k = " << k;
							//cout << "            i1 = " << i1 << " i2 = " << i2 << endl;
							M(i1, i2) += integrals[j][k];
					}}
				}
			}
		}
		/*for(int i = 0; i < n_cells; i++) {
			for(int j = 0; j < dof_per_cell; j++) {
				for(int k = 0; k < dof_per_cell; k++) {
					int i1 = i*(dof_per_cell-1)+j;
					int i2 = i*(dof_per_cell-1)+k;
					m(i1, i2) = m(i1, i2) + integrals[j][k];
				}
			}
		}*/
		//cout << M << endl;
		//m = m * scale;
		//ctrans = m.inverse();
	}
		
		double eval_(double_vector solution_, double x, double y, double z) {
			VectorXd solution = VectorXd::Map(solution_.data().begin(), solution_.size());
			return eval(solution, x, y, z);
		}
		
		double eval(VectorXd& solution, double x, double y, double z) {
		/*int cell_index = mesh.findindex(x);
		T xleft = mesh.indexto_x(cell_index);
		T xright = mesh.indexto_x(cell_index+1);
		*/
			int xi = xgrid.findindex(x);
			int yi = ygrid.findindex(y);
			int zi = zgrid.findindex(z);
			T x1 = xgrid.indexto_x(xi);
			T x2 = xgrid.indexto_x(xi+1);
			T y1 = ygrid.indexto_x(yi);
			T y2 = ygrid.indexto_x(yi+1);
			T z1 = zgrid.indexto_x(zi);
			T z2 = zgrid.indexto_x(zi+1);
			double u = (x-x1)/(x2-x1);
			double v = (y-y1)/(y2-y1);
			double w = (z-z1)/(z2-z1);
		//grid.index
		//T dx = xright-xleft;
			double value = 0;
		//cout << "xi = " << xi;
		//cout << " yi = " << yi;
			if(grid.inrange(x, y, z)) {
				basis_evaluator<3, Basis> basis_evaluator;
				for(int i = 0; i < dof_per_cell; i++) {
				//cout << " i = " << i;
				//cout << " index = " << dof_index(xi, yi, i) << endl;
					double a = solution(dof_index(xi, yi, zi, i));
					value += a * basis_evaluator.eval(i, u, v, w);
				}
			} else {
				cout << "x, y and z not in range: (" << x << "," << y << "," << z << ")" << endl;
			}
			return value;
		}
		
		double basis_uv(int i, double u, double v, double w) {
			basis_evaluator<3, Basis> basis_evaluator;
			return basis_evaluator.eval(i, u, v, w);
		}
		
		void solve_coordinates(double_vector inner_products, double_vector coordinates) {
			assert((int)inner_products.size() == dof);
			assert((int)coordinates.size() == dof);
			VectorXd x = VectorXd::Map(inner_products.data().begin(), inner_products.size());
			VectorXd a(dof);
			M.llt().solve(x, &a);
			VectorXd::Map(coordinates.data().begin(), coordinates.size()) = a;
		}
		
		/*void test(double_vector result_, double scale1, double scale2) {
			auto f = [&](double x, double y) -> double { return cos(x * scale1 + y * scale2); };
			VectorXd x = VectorXd::Zero(dof);
			assert(result_.size() == dof);
			for(int yi = 0; yi < n_cells_y; yi++) {
				for(int xi = 0; xi < n_cells_x; xi++) {
					for(int i = 0; i < dof_per_cell; i++) {
						double x1, x2;
						double y1, y2;
						grid.indexto_xy(xi, yi, x1, y1); 
						grid.indexto_xy(xi+1, yi+1, x2, y2);
					//cout << "integrate x[" << x1 << " " << x2 << "] y[" << y1 << " " << y2 << "]";
						basis_function_integrator<2, Basis::degree, Basis, Basis, T> bfi(x1, x2, y1, y2);
						double a = bfi.integrate_function2d(f, i);
						x(dof_index(xi, yi, i)) += a;
					//cout << " " << a;
 					//cout << endl;
					}
				}
			}
			VectorXd a(dof);
			M.llt().solve(x, &a);
		//VectorXd result = VectorXd::Map(result_.data().begin(), result_.size());
			VectorXd::Map(result_.data().begin(), result_.size()) = a;
		//result = a;
		}*/
		
	//int get_dof() { return dof; }
		int dof_index(int x_index, int y_index, int z_index, int basis_index) {
			if(Basis::degree == 0) {
				return x_index + y_index * dofx + z_index * dofx * dofy;
			}
			int xi = x_index*(dof_per_cell1d-1) + basis_index % dof_per_cell1d;
			int yi = y_index*(dof_per_cell1d-1) + (basis_index / dof_per_cell1d)  % dof_per_cell1d ;
			int zi = z_index*(dof_per_cell1d-1) + basis_index / (dof_per_cell1d * dof_per_cell1d);
			return xi + yi * dofx + zi * dofx * dofy;
		}
		int get_dof() { return dof; }
		int get_dof_per_cell() { return dof_per_cell; } 
	/*int get_n_cells() { return grid->;}
	int dof_index(int cell_index, int local_index) {
		return  cell_index*(dof_per_cell-1)+local_index;
	}*/
		
	};

template<class Basis, class T>
class MeshRegularNodal<1, Basis, T> {
public:
	typedef Basis basis_type;
	T x1, x2;
	int n_cells_x;
	Grid1dRegular<> xgrid;
	Grid1dRegular<>& grid;
	MatrixXd M;
//Transformation1d_in_3d* transformation;
	int dof, dofx;
	int dof1d;
	enum { dof_per_cell = MeshRegularNodalHelper<1, Basis, T>::dof_per_cell };
	enum { dof_per_cell1d = MeshRegularNodalHelper<1, Basis, T>::dof_per_cell };
	
	MeshRegularNodal(T x1, T x2, int n_cells_x) : n_cells_x(n_cells_x), xgrid(x1, x2, n_cells_x+1), grid(xgrid), M(1, 1) {
	
	if(Basis::degree == 0) {
		dofx = n_cells_x;
	} else {
		dofx = (1 + n_cells_x) + (dof_per_cell1d-2)*n_cells_x; // 1 dof per border + dofs inside the cel
	}
	dof = dofx;
	cout << "n_cells_x = " << n_cells_x << " dofx = " << dofx << endl; 
	cout << "dof = " << dof << " dof_per_cell = " << dof_per_cell << endl;
	//MatrixXd m = MatrixXd::Zero(dof, dof);
	M.resize(dof, dof);
	M =  MatrixXd::Zero(dof, dof);
	//T scale = 1; //TODO: (x2-x1)/n_cells;
	
	//basis_integrator<1, Basis::degree, Basis::degree, Basis, Basis, Basis, T> bi;
	//cout  << "test 00 " << bi.integrate(0, 0) << endl;
	//cout  << "test 01 " << bi.integrate(0, 1) << endl;
	//cout  << "test 11 " << bi.integrate(1, 1) << endl;
	
	T integrals[dof_per_cell][dof_per_cell];
	for(int j = 0; j < dof_per_cell; j++) {
		for(int k = 0; k < (j+1); k++) {
			basis_integrator<2, Basis::degree, Basis::degree, Basis, Basis, Basis, T> bi;
			T integral = bi.integrate(j, k);
			//typedef selfintegrator<Basis::degree, Basis::degree> selfintegrator_type;
			//selfintegrator_type si; 
			//double integral = si.integrate(j,k);
			//cout << j << " " << k << " " << integral << endl;
			integrals[j][k] = integral;
			integrals[k][j] = integral;
		}
	}
	for(int xi = 0; xi < n_cells_x; xi++) {
		for(int j = 0; j < dof_per_cell; j++) {
			for(int k = 0; k < dof_per_cell; k++) {
				int i1 = this->dof_index(xi, j);
				int i2 = this->dof_index(xi, k);
				//cout << "xi = " << xi << " yi = " << yi << " zi = " << zi << " j = " << j << " k = " << k;
				//cout << "            i1 = " << i1 << " i2 = " << i2 << endl;
				M(i1, i2) += integrals[j][k];
			}}
	}
	/*for(int i = 0; i < n_cells; i++) {
		for(int j = 0; j < dof_per_cell; j++) {
			for(int k = 0; k < dof_per_cell; k++) {
				int i1 = i*(dof_per_cell-1)+j;
				int i2 = i*(dof_per_cell-1)+k;
				m(i1, i2) = m(i1, i2) + integrals[j][k];
			}
		}
	}*/
	//cout << M << endl;
	//m = m * scale;
	//ctrans = m.inverse();
}
	
	double eval_(double_vector solution_, double x) {
		VectorXd solution = VectorXd::Map(solution_.data().begin(), solution_.size());
		return eval(solution, x);
	}
	
	template<class Array>
	double eval(Array& solution, double x) {
	/*int cell_index = mesh.findindex(x);
	T xleft = mesh.indexto_x(cell_index);
	T xright = mesh.indexto_x(cell_index+1);
	*/
		int xi = xgrid.findindex(x);
		T x1 = xgrid.indexto_x(xi);
		T x2 = xgrid.indexto_x(xi+1);
		double u = (x-x1)/(x2-x1);
	//grid.index
	//T dx = xright-xleft;
		double value = 0;
	//cout << "xi = " << xi;
	//cout << " yi = " << yi;
		if(grid.inrange(x)) {
			basis_evaluator<1, Basis> basis_evaluator;
			for(int i = 0; i < dof_per_cell; i++) {
			//cout << " i = " << i;
			//cout << " index = " << dof_index(xi, yi, i) << endl;
				double a = solution[dof_index(xi, i)];
				value += a * basis_evaluator.eval(i, u);
			}
		} else {
		cout << "x not in range: (" << x << endl;
		}
		return value;
	}	
	
	double eval(VectorXd& solution, double x) {
	/*int cell_index = mesh.findindex(x);
	T xleft = mesh.indexto_x(cell_index);
	T xright = mesh.indexto_x(cell_index+1);
	*/
		int xi = xgrid.findindex(x);
		T x1 = xgrid.indexto_x(xi);
		T x2 = xgrid.indexto_x(xi+1);
		double u = (x-x1)/(x2-x1);
	//grid.index
	//T dx = xright-xleft;
		double value = 0;
	//cout << "xi = " << xi;
	//cout << " yi = " << yi;
		if(grid.inrange(x)) {
			basis_evaluator<1, Basis> basis_evaluator;
			for(int i = 0; i < dof_per_cell; i++) {
			//cout << " i = " << i;
			//cout << " index = " << dof_index(xi, yi, i) << endl;
				double a = solution(dof_index(xi, i));
				value += a * basis_evaluator.eval(i, u);
			}
		} else {
			cout << "x not in range: (" << x << endl;
		}
		return value;
	}
	
	double basis_uv(int i, double u) {
		basis_evaluator<1, Basis> basis_evaluator;
		return basis_evaluator.eval(i, u);
	}
	
	void solve_coordinates(double_vector inner_products, double_vector coordinates) {
		assert((int)inner_products.size() == dof);
		assert((int)coordinates.size() == dof);
		VectorXd x = VectorXd::Map(inner_products.data().begin(), inner_products.size());
		VectorXd a(dof);
		M.llt().solve(x, &a);
		VectorXd::Map(coordinates.data().begin(), coordinates.size()) = a;
	}
	
	/*void test(double_vector result_, double scale1, double scale2) {
		auto f = [&](double x, double y) -> double { return cos(x * scale1 + y * scale2); };
		VectorXd x = VectorXd::Zero(dof);
		assert(result_.size() == dof);
		for(int yi = 0; yi < n_cells_y; yi++) {
			for(int xi = 0; xi < n_cells_x; xi++) {
				for(int i = 0; i < dof_per_cell; i++) {
					double x1, x2;
					double y1, y2;
					grid.indexto_xy(xi, yi, x1, y1); 
					grid.indexto_xy(xi+1, yi+1, x2, y2);
				//cout << "integrate x[" << x1 << " " << x2 << "] y[" << y1 << " " << y2 << "]";
					basis_function_integrator<2, Basis::degree, Basis, Basis, T> bfi(x1, x2, y1, y2);
					double a = bfi.integrate_function2d(f, i);
					x(dof_index(xi, yi, i)) += a;
				//cout << " " << a;
				//cout << endl;
				}
			}
		}
		VectorXd a(dof);
		M.llt().solve(x, &a);
	//VectorXd result = VectorXd::Map(result_.data().begin(), result_.size());
		VectorXd::Map(result_.data().begin(), result_.size()) = a;
	//result = a;
	}*/
	
//int get_dof() { return dof; }
	int dof_index(int x_index, int basis_index) {
		if(Basis::degree == 0) {
			return x_index;
		}
		int xi = x_index*(dof_per_cell1d-1) + basis_index % dof_per_cell1d;
		return xi;
	}
	int get_dof() { return dof; }
	int get_dof_per_cell() { return dof_per_cell; } 
/*int get_n_cells() { return grid->;}
int dof_index(int cell_index, int local_index) {
	return  cell_index*(dof_per_cell-1)+local_index;
}*/
	
};

}