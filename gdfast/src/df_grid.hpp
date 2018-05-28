#include "datatypes.hpp"

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <fe/fe_q.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include "galaxy.hpp"

#include <fstream>


namespace gd {
using namespace dealii;
using namespace std;
 
class DFGrid {
public:
	DFGrid(int n_I1, int n_I2, double r1, double r2, Galaxy* galaxy) : galaxy(galaxy), r1(r1), r2(r2), dof_handler (triangulation), fe(1) {
		const Point<2> p1(r1, 0);
		const Point<2> p2(r2, 1);
		std::vector<unsigned int> repetitions = {n_I1, n_I2};
		GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p1, p2);
		dof_handler.distribute_dofs (fe);
	}
	
	int n_dofs() {
		return dof_handler.n_dofs();
	}
	
	void print_dof_indices_per_face() {
		DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
		const unsigned int   dofs_per_cell = fe.dofs_per_cell;
		std::vector<unsigned int> local_dof_indices (dofs_per_cell);
		cout << "dof number per cell" << endl;
		for (; cell!=endc; cell++) {
			cell->get_dof_indices (local_dof_indices);
			cout << "cell: " << (cell->used()  ? "  used" : "unused") << " " << (cell->active()  ? "  active" : "inactive") << " " ;
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				cout << " " << local_dof_indices[i];
			cout << endl; 
		}
	}
	
	void print_vertices() {
		std::vector<Point<2>> vertices = triangulation.get_vertices();
		cout << "# of used vertices " << vertices.size();
		auto v = vertices.begin();
		while(v != vertices.end()) {
			double x = (*v)(0);
			double y = (*v)(1);
			cout << " " << x << "," << y << endl;
			v++;
		}
	}
	
	void print_dof_indices_per_face_on_level(int level) {
		DoFHandler<2>::raw_cell_iterator cell = dof_handler.begin_raw(level), endc = dof_handler.end(level);
		const unsigned int   dofs_per_cell = fe.dofs_per_cell;
		std::vector<unsigned int> local_dof_indices (dofs_per_cell);
		cout << "dof number per cell on level: "  << level << endl;
		//for (; cell!=endc; ++cell) {
		while(cell != endc) {
		//if(1) {
			cell->get_dof_indices (local_dof_indices);
			cout << "cell: ";
			for (unsigned int i=0; i<dofs_per_cell; ++i)
				cout << " " << local_dof_indices[i];
			cout << endl; 
			cell++;
		}
	}
	
	double shape_value(unsigned int i, double x, double y) {
		Point<2> p(x, y);
		fe.shape_value(i, p);
	} 
	
	void refine_global(int count) {
		triangulation.refine_global(count);
		dof_handler.distribute_dofs (fe);
	}
	
	void output_grid( char* filename) {
		std::ofstream out (filename);
		GridOut grid_out;
		grid_out.write_eps (triangulation, out);
	}

	double r1;
	double r2;
	Galaxy* galaxy;
	Triangulation<2>     triangulation;
	FE_Q<2>              fe;
	DoFHandler<2>        dof_handler;
	
};


};