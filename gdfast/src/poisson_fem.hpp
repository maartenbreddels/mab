#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_refinement.h>
#include <grid/grid_out.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>

#include <dofs/dof_tools.h>

#include <fe/fe_values.h>
#include <base/quadrature_lib.h>

#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/fe_field_function.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/constraint_matrix.h>

#include <numerics/data_out.h>
#include <numerics/error_estimator.h>
#include <fstream>
#include <iostream>
#include "datatypes.hpp"

namespace gd {

using namespace dealii;

template<int dim>
class PoissonFEM {
public:
	PoissonFEM(int order);
	void make_grid();
	void setup_system();
	void refine_grid();
	void assemble_system ();
	void solve ();
	void output_results () const;
	void output_grid(char* filename) const;
	double eval(double r);
	double grad(double r);

	int n_active_cells() const;

	void getpoints(double_vector p); // for 1d only
	void global_refine(int refine);


	Triangulation<dim>     triangulation;
	FE_Q<dim>              fe;
	DoFHandler<dim>        dof_handler;
	
	ConstraintMatrix     hanging_node_constraints;
	
	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;
	
	Vector<double>       solution;
	Vector<double>       system_rhs;

};

}