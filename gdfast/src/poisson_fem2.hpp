#pragma once
#include "profile.hpp"
#include "mesh.hpp"
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include "transformation.hpp"
#undef TEST_SPARSE


namespace gd {

USING_PART_OF_NAMESPACE_EIGEN;

template<class Mesh>
class PoissonSolver1d;
		
const double ProfileNumerical1d_epsilon = 1e-5;
template<class Mesh>
	class ProfileNumerical1d : public Profile {
	int n;
	Mesh mesh;
	Density* density;
	PoissonSolver1d<Mesh> solver;
	VectorXd solution;
	Transformation1d_in_3d* transformation;
	double G;
public:
	
	ProfileNumerical1d(int n, Density* density, Transformation1d_in_3d* transformation, double G, double u1, double u2) : mesh(u1, u2, n, transformation), density(density), solver(density, &mesh, G), solution(mesh.get_dof()), transformation(transformation), G(G) {
		solver.solve_(solution, 1., 1., 0.);
	}
	virtual double densityr(double r) {
		return density->densityr(r);
	}
	virtual double densityR(double) {
		return 0;
	}
	virtual double I(double, double) {
		return 0;
	}  
	virtual double dphidr(double r) {
		double u = transformation->inverse_transform(r);
		//double u = atan(r) * 2 / M_PI;
		//double jacobian = M_PI/2 / pow(cos(u*M_PI/2), 2);
		//double du
		//return mesh.gradient(solution, u) / jacobian;
		return mesh.gradient(solution, u) / transformation->drdu(u); //jacobian;
	}
	virtual double potentialr(double r) {
		//double u = atan(r) * 2 / M_PI;
		double u = transformation->inverse_transform(r);
		return mesh.eval(solution, u);
	}
};

template<class Mesh>
class PoissonSolver1d {
public:
	typedef Mesh mesh_type;
	Density* density;
	Mesh* mesh;
	VectorXd lasta;
	double G;
	
	PoissonSolver1d(Density* density, Mesh* mesh, double G) : density(density), mesh(mesh), lasta(mesh->get_dof()), G(G) {
	}
	
	double operator()(double x) {
		return mesh->eval(lasta, x);
	}
	
	double gradient(double x) {
		return mesh->gradient(lasta, x);
	}
	
	void solve(double_vector v, double scale1, double scale2, double boundary_value=0) {
		VectorXd v_copy = VectorXd::Map(v.data().begin(), v.size());
		solve_(v_copy, scale1, scale2, boundary_value);
	}
	void solve_(VectorXd& v, double scale1, double scale2, double boundary_value=0) {
		
		/*
		solve 'a' from the linear system Ma=x using FEM (Galerkin method)
		such that \Phi(r) = \sum_i a_i \phi_i(r) is the solution to the 
		poisson eq: \delta^2 \Phi(r) = 4 pi G rho(r) 
		*/
		int dof = mesh->get_dof();
#ifdef TEST_SPARSE
		Eigen::DynamicSparseMatrix<double> Ms(dof, dof);
#else
		MatrixXd M = MatrixXd::Zero(dof, dof);
#endif
		VectorXd x = VectorXd::Zero(dof);
		
		int dof_per_cell = Mesh::dof_per_cell;
		for(int cell_index = 0; cell_index < mesh->get_n_cells(); cell_index++) {
			for(int i = 0; i < dof_per_cell; i++) {
				for(int j = 0; j < dof_per_cell; j++) {
					//M(cell_index*(dof_per_cell-1)+i, cell_index*(dof_per_cell-1)+j) = integrate dphi_i * dphi_j
#ifdef TEST_SPARSE
					Ms.coeffRef(mesh->dof_index(cell_index, i), mesh->dof_index(cell_index, j)) +=
#else
					M(mesh->dof_index(cell_index, i), mesh->dof_index(cell_index, j)) +=
#endif
						mesh->integrate_gradshape(cell_index, i, j) * scale1; // * r * r;
						/*integrate dphi_i * dphi_j*/
				}
				//x(mesh->index(cell_index, i)) += integrate phi_i * 4 * M_PI * density->densityr(r) dr;
				//auto f = [&](double r) { return this->density->densityr(r) * r * r; }; //*/ };
				auto f = [&](double u) {
					double r = this->mesh->transformation->transform(u);
					//double r = tan(u*M_PI/2);
					//double s = sin(u*M_PI/2);
					//double c = cos(u*M_PI/2);
					//return this->density->densityr(r) * 2 * M_PI*M_PI * s*s/pow(c,4);
					return this->density->densityr(r) * this->mesh->transformation->d3xdu(u);
				};
				//double G = 1;
				//cout << mesh->dof_index(cell_index, i) << " = " << (-4 * M_PI * G * mesh->integrate_shape(cell_index, i, f)) << endl;
				x(mesh->dof_index(cell_index, i)) += -(4 * M_PI) * G * mesh->integrate_shape(cell_index, i, f) * scale2;
			}
		}
		// set boundary condition, Phi(r_end) = 0
		for(int i = 1; i < Mesh::dof_per_cell; i++) {
#ifdef TEST_SPARSE
			Ms.coeffRef(dof-1-i,dof-1) = boundary_value;
			Ms.coeffRef(dof-1,dof-1-i) = boundary_value;
#else
			M(dof-1-i,dof-1) = boundary_value;
			M(dof-1,dof-1-i) = boundary_value;
#endif
		}
		
#ifdef TEST_SPARSE
		Ms.coeffRef(dof-1,dof-1) = 1;
#else
		M(dof-1,dof-1) = 1;
#endif

#ifdef TEST_SPARSE
		typedef Eigen::SparseMatrix<double> SparseMatrixType;
		SparseMatrixType M(Ms);
#endif
		x(dof-1) = 0;
		
		//cout << M << endl;
		//cout << "next" << endl << x << endl;
		// solve 'a'
#ifdef TEST_SPARSE
		VectorXd a = x;
		Eigen::SparseLLT<SparseMatrixType,Eigen::Cholmod> sparseLLT(M);
		sparseLLT.solveInPlace(a);
#else
		//VectorXd a = M.inverse() * x;
		VectorXd a(dof);
		M.llt().solve(x, &a);
#endif
		lasta = a;
		// copy to v 
		//VectorXd::Map(v.data().begin(), v.size()) = a;
		v = a;
		//cout << "x: " << endl <<  (x) << endl;
		
		//cout << "solution: " << endl <<  (a) << endl;
		//cout << "test" << endl <<  (M * a) << endl;
	}
};

}