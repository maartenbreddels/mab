#include "poisson_fem.hpp"
#include <boost/python.hpp>
//#include "vectorize.hpp"

namespace gd {

using namespace boost::python;
using namespace std;
using namespace dealii;

void py_export_poisson_fem() {
	class_<PoissonFEM<2>, boost::noncopyable >("PoissonFEM2d", init<int>())
		.def("make_grid", &PoissonFEM<2>::make_grid)
		.def("setup_system", &PoissonFEM<2>::setup_system)
		.def("refine_grid", &PoissonFEM<2>::refine_grid)
		.def("assemble_system", &PoissonFEM<2>::assemble_system)
		.def("solve", &PoissonFEM<2>::solve)
		.def("output_results", &PoissonFEM<2>::output_results)
		.def("output_grid", &PoissonFEM<2>::output_grid)
		.def("eval", (&PoissonFEM<2>::eval))
		.def("n_active_cells", (&PoissonFEM<2>::n_active_cells))
		.def("global_refine", (&PoissonFEM<2>::global_refine))
	;
	class_<PoissonFEM<1>, boost::noncopyable >("PoissonFEM1d", init<int>())
		.def("make_grid", &PoissonFEM<1>::make_grid)
		.def("setup_system", &PoissonFEM<1>::setup_system)
		.def("refine_grid", &PoissonFEM<1>::refine_grid)
		.def("assemble_system", &PoissonFEM<1>::assemble_system)
		.def("solve", &PoissonFEM<1>::solve)
		.def("output_results", &PoissonFEM<1>::output_results)
		.def("output_grid", &PoissonFEM<1>::output_grid)
		.def("eval", (&PoissonFEM<1>::eval))
		.def("getpoints", (&PoissonFEM<1>::getpoints))
		.def("n_active_cells", (&PoissonFEM<1>::n_active_cells))
		.def("grad", (&PoissonFEM<1>::grad))
		.def("global_refine", (&PoissonFEM<1>::global_refine))
			
	;
}


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
	//return -exp(-pow(r/4-phi, 2));
	//return ((r > 10) && (r < 15.0)) * 1.; //3/(4*M_PI) * pow(1+pow(r,2), -5./2);
	//return -1; 
	//return -1/(1+x/20.);
	return -3/(4*M_PI) * pow(1+pow(r/20,2), -5./2);
}


template <int dim>
PoissonFEM<dim>::PoissonFEM (int order) :
                 fe (order),
                dof_handler (triangulation)
 {}

template <int dim>
void PoissonFEM<dim>::make_grid()
{
	const Point<dim> center;

	if(dim == 1) {
		GridGenerator::hyper_cube(triangulation, 0, 40) ;//, 10);
	} else {
		GridGenerator::hyper_ball(triangulation, center, 40.) ;//, 10);
		//const HyperShellBoundary<dim> boundary_description(center);
		//triangulation.set_boundary (0, boundary_description);
	}
	//const HyperShellBoundary<dim> boundary_description(center);
	//const HyperShellBoundary<2> boundary_description2(center);
	//triangulation.set_boundary (0, boundary_description);
	//triangulation.set_boundary (1, boundary_description2);
	
	triangulation.refine_global (2);
	//triangulation.set_boundary (0);
	
	std::cout << "Number of active cells: "
				<< triangulation.n_active_cells()
				<< std::endl;
	std::cout << "Total number of cells: "
				<< triangulation.n_cells()
				<< std::endl;
	
}

template <int dim>
void  PoissonFEM<dim>::setup_system() {
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
void PoissonFEM<dim>::assemble_system () 
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
	for (; cell!=endc; ++cell)  {
		fe_values.reinit (cell);
	
		cell_matrix = 0;
		cell_rhs = 0;
 
	   for (unsigned int i=0; i<dofs_per_cell; ++i) {
			for (unsigned int j=0; j<dofs_per_cell; ++j) {
				for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
					//cout << "q: " << fe_values.quadrature_point (q_point) << endl;
					double r = fe_values.quadrature_point (q_point)(0);
					cell_matrix(i,j) += ( fe_values.shape_grad (i, q_point) *
									fe_values.shape_grad (j, q_point) *
					                      fe_values.JxW (q_point)) ;/* r * r; //*/
					cout << "-> " <<fe_values.JxW (q_point) << " " << r << " " << fe_values.shape_grad (i, q_point) << endl;
				}
			}
		}
		
		for (unsigned int i=0; i<dofs_per_cell; ++i) {
			for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
				double r = fe_values.quadrature_point (q_point)(0);
				cell_rhs(i) += (fe_values.shape_value (i, q_point) *
							right_hand_side.value (fe_values.quadrature_point (q_point)) *
				                fe_values.JxW (q_point)) ;/* r * r; //*/
	
			
				cell->get_dof_indices (local_dof_indices);
				
			}
		}
	
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
			system_matrix.add (local_dof_indices[i],
								local_dof_indices[j],
								cell_matrix(i,j));
	
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	//cout << system_matrix << endl;
	system_matrix.print(cout);
	cout << "next" << endl;
	system_rhs.print(cout);
	hanging_node_constraints.condense (system_matrix);
	hanging_node_constraints.condense (system_rhs);
	
	std::map<unsigned int,double> boundary_values;
	//boundary_values[0] = 0;
	if(dim == 1) {
		VectorTools::interpolate_boundary_values (dof_handler,
												1,
												ConstantFunction<dim>(0),
												boundary_values);
	} else {
		VectorTools::interpolate_boundary_values (dof_handler,
				0,
				ConstantFunction<dim>(0),
				boundary_values);
	}
	MatrixTools::apply_boundary_values (boundary_values,
										system_matrix,
										solution,
										system_rhs);
	cout << "boundary:" << endl;
	system_matrix.print(cout);
	cout << "next" << endl;
	system_rhs.print(cout);
}

template <int dim>
Point<dim> pointr(double r);

template<>
Point<1> pointr<1>(double r) {
	return Point<1>(r);
}

template<>
Point<2> pointr<2>(double r) {
	return Point<2>(r, 0);
}

template<> 
Point<3> pointr(double r) {
	return Point<3>(r, 0, 0);
}

template <int dim>
 void PoissonFEM<dim>::solve () 
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
void PoissonFEM<dim>::getpoints(double_vector p) // for 1d only
{
	double* pp = p.data().begin();
	/*Triangulation<2>::active_line_iterator
		cell = triangulation.begin_active_line(),
		endc = triangulation.end_line();

	double nextx;
	const Point<dim> center;
	cout << "starting.. " << (cell==endc) << endl; 
	for (; cell!=endc; ++cell) {
		cout << "-> " << cell->vertex(0)(0) << endl;
		*pp++ = (cell->vertex(0))(0);
		nextx =  (cell->vertex(1))(0);
	}
	*pp++ = nextx; // only include last point of last cell
	*/
	auto vertices = triangulation.get_vertices();
	auto used = triangulation.get_used_vertices();
	for(int i = 0; i < vertices.size(); i++) {
		if(used[i]) {
			*pp++ = vertices[i](0);
		}
	}
}

template <int dim>
double PoissonFEM<dim>::eval(double r) {
	Functions::FEFieldFunction<dim> ff(dof_handler, solution);
	return ff.value(pointr<dim>(r));
	//return VectorTools::point_value (dof_handler, solution, pointr<dim>(r));
}

template <int dim>
double PoissonFEM<dim>::grad(double r) {
	Functions::FEFieldFunction<dim> ff(dof_handler, solution);
	auto tensor = ff.gradient(pointr<dim>(r));
	return tensor[0];
}

template <int dim>
void PoissonFEM<dim>::refine_grid ()
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
void PoissonFEM<dim>::output_results () const
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
template <int dim>
void PoissonFEM<dim>::output_grid(char* filename) const
{
	//DataOut<dim> data_out;
	//data_out.attach_dof_handler (dof_handler);
	//data_out.add_data_vector (solution, "solution");
	//data_out.build_patches ();
	
	std::ofstream outgrid(filename);
	GridOut grid_out;
	grid_out.write_eps(triangulation, outgrid);
	
}

template <int dim>
int PoissonFEM<dim>::n_active_cells() const
{
	return triangulation.n_active_cells();
}
int refines = 1;


template <int dim>
void PoissonFEM<dim>::global_refine(int refine)
{
	return triangulation.refine_global(refine);
}

/*
template <int dim>
void PoissonFEM<dim>::run () 
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
*/
 /*int main (int argc, char** argv) 
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
 }*/

}