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
 
 using namespace dealii;

 template <int dim>
 class RightHandSide : public Function<dim> 
 {
   public:
     RightHandSide () : Function<dim>() {}
     
     virtual double value (const Point<dim>   &p,
                           const unsigned int  component = 0) const;
 };
 template <int dim>
 double RightHandSide<dim>::value (const Point<dim> &p,
                                   const unsigned int /* component */) const 
{
	double rsq = 0;
	for (unsigned int i=0; i<dim; ++i)
		rsq += std::pow(p(i), 2);
	double r = sqrt(rsq);
	double x = p(0);
	double y = p(1);
	double phi = fmod(atan2(y, x)+2*M_PI, 2*M_PI);
	return -exp(-pow(r/4-phi, 2));
	//return ((r > 10) && (r < 15.0)) * 1.; //3/(4*M_PI) * pow(1+pow(r,2), -5./2); 
	//return -1; //-3/(4*M_PI) * pow(1+pow(r,2), -5./2);
}

template<int dim>
class LaplaceProblem 
{
public:
	LaplaceProblem ();
     void run ();
   private:
	void make_grid();
	void setup_system();
	void refine_grid();
	void assemble_system ();
	void solve ();
	void output_results () const;

	Triangulation<dim>     triangulation;
	FE_Q<dim>              fe;
	DoFHandler<dim>        dof_handler;
	
	ConstraintMatrix     hanging_node_constraints;
	
	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;
	
	Vector<double>       solution;
	Vector<double>       system_rhs;
 };

template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                 fe (1),
                dof_handler (triangulation)
 {}

template <int dim>
void LaplaceProblem<dim>::make_grid()
{
	const Point<dim> center (0,0);

   GridGenerator::hyper_ball(triangulation, center, 40) ;//, 10);
   const HyperShellBoundary<dim> boundary_description(center);
   //const HyperShellBoundary<2> boundary_description2(center);
   triangulation.set_boundary (0, boundary_description);
   //triangulation.set_boundary (1, boundary_description2);

   triangulation.refine_global (2);
	triangulation.set_boundary (0);

   std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
   std::cout << "Total number of cells: "
            << triangulation.n_cells()
            << std::endl;

}

template <int dim>
void  LaplaceProblem<dim>::setup_system() {
   dof_handler.distribute_dofs (fe);
   std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
 
   sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
   solution.reinit (dof_handler.n_dofs());
   system_rhs.reinit (dof_handler.n_dofs());


	hanging_node_constraints.clear ();
	DoFTools::make_hanging_node_constraints (dof_handler,
                                            hanging_node_constraints);
hanging_node_constraints.condense (sparsity_pattern);

	hanging_node_constraints.close ();

   sparsity_pattern.compress();
 
   system_matrix.reinit (sparsity_pattern);
 
}

template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{
   QGauss<dim>  quadrature_formula(2);
   FEValues<dim> fe_values (fe, quadrature_formula, 
                         update_values | update_gradients | update_JxW_values | update_quadrature_points);

	const RightHandSide<dim> right_hand_side;

   
   const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   const unsigned int   n_q_points    = quadrature_formula.size();
 
   FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
   Vector<double>       cell_rhs (dofs_per_cell);
 
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
 
   typename DoFHandler<dim>::active_cell_iterator
     cell = dof_handler.begin_active(),
     endc = dof_handler.end();
   for (; cell!=endc; ++cell)
     {
       fe_values.reinit (cell);
 
       cell_matrix = 0;
       cell_rhs = 0;
 
       for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                 fe_values.shape_grad (j, q_point) *
                                 fe_values.JxW (q_point));
 
       for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                          right_hand_side.value (fe_values.quadrature_point (q_point)) *
                          fe_values.JxW (q_point));
 
       cell->get_dof_indices (local_dof_indices);
 
       for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             cell_matrix(i,j));
 
       for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
     }
 
   hanging_node_constraints.condense (system_matrix);
   hanging_node_constraints.condense (system_rhs);

	std::map<unsigned int,double> boundary_values;
	//boundary_values[0] = 0;
	VectorTools::interpolate_boundary_values (dof_handler,
											0,
											ZeroFunction<dim>(),
											boundary_values);
	MatrixTools::apply_boundary_values (boundary_values,
										system_matrix,
										solution,
										system_rhs);
 }

template <int dim>
Point<dim> pointr(double r);

template<>
Point<2> pointr<2>(double r) {
	return Point<2>(r, 0);
}

template<> 
Point<3> pointr(double r) {
	return Point<3>(r, 0, 0);
}

template <int dim>
 void LaplaceProblem<dim>::solve () 
 {
   SolverControl           solver_control (5000, 1e-12);
   SolverCG<>              cg (solver_control);
 
   cg.solve (system_matrix, solution, system_rhs,
            PreconditionIdentity());
hanging_node_constraints.distribute (solution);

	std::cout << VectorTools::point_value (dof_handler, solution, pointr<dim>(0)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, pointr<dim>(5)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, pointr<dim>(10)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, pointr<dim>(15)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, pointr<dim>(20)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, pointr<dim>(40)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, pointr<dim>(0.1)) << std::endl;
	/*std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(40, 0)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0, 40)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0, 20)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(1,0)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0,1)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(2,0)) << std::endl;
	double a = VectorTools::point_value (dof_handler, solution, Point<2>(1,0))-VectorTools::point_value (dof_handler, solution, Point<2>(0,0));
	double b = VectorTools::point_value (dof_handler, solution, Point<2>(2,0))-VectorTools::point_value (dof_handler, solution, Point<2>(1,0));
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "a/b = " << (a/b) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0, 0)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0, 5)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0, 10)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0, 15)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0, 20)) << std::endl;
	std::cout << VectorTools::point_value (dof_handler, solution, Point<2>(0, 25)) << std::endl;*/

 }

template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
   Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
 
   KellyErrorEstimator<dim>::estimate (dof_handler,
                                       QGauss<dim-1>(3),
                                       typename FunctionMap<dim>::type(),
                                       solution,
                                       estimated_error_per_cell);
   GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                    estimated_error_per_cell,
                                                    0.3, 0.03);
triangulation.execute_coarsening_and_refinement ();

}

template <int dim>
void LaplaceProblem<dim>::output_results () const
 {
   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (solution, "solution");
   data_out.build_patches ();
 
   std::ofstream output (dim == 2 ? "solution2d.gpl" : "solution3d.gpl" );
   data_out.write_gnuplot (output);
   
	std::ofstream outgrid (dim == 2 ? "grid2d.eps" : "grid3d.eps");
	GridOut grid_out;
	grid_out.write_eps (triangulation, outgrid);

 
   std::ofstream output2(dim == 2 ? "solution2d.vtk" : "solution3d.vtk" );
   data_out.write_vtk (output2);

 }

int refines = 1;

template <int dim>
void LaplaceProblem<dim>::run () 
{
	make_grid();
	setup_system();
	assemble_system ();
	solve ();
	
	for(int i = 0; i < refines; i++) {
	
		refine_grid();
		setup_system();
		assemble_system ();
		solve ();
	}
	output_results ();
 }

 int main (int argc, char** argv) 
 {
	if(argc < 2) {
		printf("usage: %s <refinements>\n", argv[0]);
		exit(-1);
	}
	refines = atoi(argv[1]);
	LaplaceProblem<2> laplace_problem;
	laplace_problem.run ();

	LaplaceProblem<3> laplace_problem3d;
	laplace_problem3d.run ();
   return 0;
 }