#include "hyperquadtree.hpp"
//#include <boost/python.hpp>
#include <algorithm>
#include <iterator>
namespace gd {
//using namespace boost::python;

void py_export_schw_hyperquadtree() {
	//NTree<2, double> quadtree(1.,2.,4.,6.);
	/*class_<DFGrid, boost::noncopyable >("DFGrid", init<int, int, double, double, Galaxy*>())
		.def("output_grid", &DFGrid::output_grid) 
		.def("n_dofs", &DFGrid::n_dofs) 
		.def("refine_global", &DFGrid::refine_global)
		.def("shape_value", &DFGrid::shape_value)
		.def("print_dof_indices_per_face", &DFGrid::print_dof_indices_per_face)
		.def("print_dof_indices_per_face_on_level", &DFGrid::print_dof_indices_per_face_on_level)
		.def("print_vertices", &DFGrid::print_vertices)
			;*/
}

//template<size_t... Indices>
//tuple<> get_indices(

struct DoSomething
{
    template<typename T>
    void operator()(T& t) const
    {
        //t.do_sth();
		//printf("blaat");
    }
};

const int gs_scale = 300;

// this makes gcc happy
template<size_t N, typename... Tail>
struct action;

/*template<int A, int B>
struct static_min {
	enum { value = A < B : A ? B };
};*/
template<int N>
struct gridder {
	template<class Tree, class P>
	void grid(Tree *tree, ofstream& f, P p1, P p2, int *gridsizes) {
		cout << "Error: gridding not implemented for N = " << N << endl;
	}
};

template<>
struct gridder<2>
{
	template<class Tree, class P>
	void grid(Tree *tree, ofstream& f, P p1, P p2, int *gridsizes) {
		double x1 = get<0>(p1);
		double x2 = get<0>(p2);
		double y1 = get<1>(p1);
		double y2 = get<1>(p2);
		cout << "x range" << x1 << " - " << x2 << endl;
		cout << "y range" << y1 << " - " << y2 << endl;
		f << "[";
		for(int xi = 0; xi < gridsizes[0]; xi++) {
			double x = x1 + xi * (x2-x1) / (gridsizes[0]-1);
			f << "[";
			for(int yi = 0; yi < gridsizes[1]; yi++) {
				double y = y1 + yi * (y2-y1) / (gridsizes[1]-1);
				//cout << "x,y" << x << "\t" << y << "\t" << tree->eval(x, y) << endl;
				//cout << "x,y" << x << "\t" << y << "\t" << tree->eval(x, y) << endl;
				f << tree->eval(x, y) << ",";
			}
			f << "],";
			//cout << endl;
		}
		f << "]";
		f.close();
	}
};
template<>
struct gridder<3>
{
	template<class Tree, class P>
	void grid(Tree *tree, ofstream& f, P p1, P p2, int *gridsizes) {
		double x1 = get<0>(p1);
		double x2 = get<0>(p2);
		double y1 = get<1>(p1);
		double y2 = get<1>(p2);
		double z1 = get<2>(p1);
		double z2 = get<2>(p2);
		cout << "x range" << x1 << " - " << x2 << endl;
		cout << "y range" << y1 << " - " << y2 << endl;
		cout << "z range" << z1 << " - " << z2 << endl;
		f << "[";
		for(int xi = 0; xi < gridsizes[0]; xi++) {
			f << "[";
			double x = x1 + xi * (x2-x1) / (gridsizes[0]-1);
			for(int yi = 0; yi < gridsizes[1]; yi++) {
				f << "[";
				double y = y1 + yi * (y2-y1) / (gridsizes[1]-1);
				for(int zi = 0; zi < gridsizes[2]; zi++) {
					double z = z1 + zi * (z2-z1) / (gridsizes[2]-1);
					//cout << "x,y" << x << "\t" << y << "\t" << tree->eval(x, y) << endl;
					//cout << "x,y" << x << "\t" << y << "\t" << tree->eval(x, y) << endl;
					f << tree->eval(x, y, z) << ",";
				}
				f << "],";
			}
			//cout << endl;
			f << "],";
		}
		f << "]";
		f.close();
	}
};

template<size_t N, typename T, typename... Tail>
struct action<N, T, Tail...>  {
	typedef NTree<N, T, T, Tail...> HTree;
	typedef typename HTree::Point Point;
	typedef action<N-1, Tail...> sub_action_type;
	int dim;
	sub_action_type sub_action;
	HTree* tree;
	Point p1, p2;

	action(int dim) : dim(dim), sub_action(dim), tree(NULL) {}
	void init(int subdivides, string name) {
		//QuadTree::Point p1(0.0,0.0);		
		if(N == dim) {
			tree = new HTree(p1, p2);
			//cout << "output: " << name << endl;
			while(subdivides--) tree->rootNode.split();
			//tree->output(name);
			//tree.output_unknown(name);
			//tree->output_ps((name+".ps").c_str(), gs_scale, gs_scale);
		}
		else if(N >= 1) {
			sub_action.init(subdivides, name);
		}
	}
	void optimize(double absfraction=0.2, double relfraction=0.05, int Nabsmin=1, int Nrelmin=0, int maxlevel=8, int max_level_difference=0) {
		if(N == dim) {
			double error = 0;
			cout << "integral is: " << tree->integrate(error) << endl;
			cout << "error    is: " << error << endl;
			tree->optimize(absfraction, relfraction, Nabsmin, Nrelmin, maxlevel, max_level_difference);
		} else {
			sub_action.optimize(absfraction, relfraction, Nabsmin, Nrelmin, maxlevel, max_level_difference);
		}
	}
	void limit_level_difference(int max_level_difference) {
		if(N == dim) {
			//double error = 0;
			tree->limit_level_difference(max_level_difference);
		} else {
			sub_action.limit_level_difference(max_level_difference);
		}
	} 
	void read(string inputname) {
		if(N == dim) {
			ifstream f;
			string filename_state = inputname + ".state";
			f.open(filename_state.c_str());
			cout << "reading from: " << filename_state << endl; 
			if(!f.is_open())
				throw runtime_error((string("file doesn't exists: ") + filename_state));
			
			//Point p1, p2;
			f >> p1 >> p2;
			cout << p1 << ", " << p2 << endl;

			vector<bool> splits;
			copy(istream_iterator<bool>(f), istream_iterator<bool>(), back_inserter(splits));
			//copy(splits.begin(), splits.end(), ostream_iterator<bool>(cout, " "));
			cout << "size: " << splits.size() << endl;
			
			//HTree tree(p1, p2);
			tree = new HTree(p1, p2);
			auto b = splits.begin();
			tree->rootNode.split(&b);
			assert(b == splits.end());
			tree->readpoints(inputname+".known");
			tree->readpoints(inputname+".solved");
			/*tree.output_ps((outputname+".check.ps").c_str(), gs_scale, gs_scale);
			tree.readpoints(inputname+".known");
			tree.readpoints(inputname+".solved");
			double error = 0;
			cout << "integral is: " << tree.integrate(error) << endl;
			cout << "error    is: " << error << endl;
			tree.optimize(0.0, 0.0, 10, 0);
			tree.output(outputname);
			tree.output_ps((outputname+".ps").c_str(), gs_scale, gs_scale);*/

		} else if(N >= 1) {
			sub_action.read(inputname);
		}
	}
	void write(string filename) {
		if(N == dim) {
			tree->output(filename);
			tree->output_ps((filename+".eps").c_str(), gs_scale, gs_scale);
		}
		else if(N >= 1) {
			sub_action.write(filename);
		}
	}
	void grid(string name, int* gridsizes) {
		if(N == dim) {
			ofstream f;
			string filename = name+".grid";
			f.open(filename);
			cout << "output filename: " << filename << endl;
			gridder<N> g;
			g.grid(tree, f, p1, p2, gridsizes);
		}
		else if(N >= 1) {
			sub_action.grid(name, gridsizes);
		}
	}
	void set_point1(T value, int point_dim) {
		if(N == dim) {
			do_tuple(p1, point_dim, [&](T& tuple_ref, bool last) { tuple_ref = value;} ); 
		}
		else if(N >= 1) {
			sub_action.set_point1(value, point_dim);
		}
	}
	void set_point2(T value, int point_dim) {
		if(N == dim) {
			do_tuple(p2, point_dim, [&](T& tuple_ref, bool last) { tuple_ref = value;} ); 
		}
		else if(N >= 1) {
			sub_action.set_point2(value, point_dim);
		}
	}
	
}; 
template<typename T>
struct action<1, T> {
	action(int dim) {}
	
	void set_point1(T value, int point_dim) {
	}
	void set_point2(T value, int point_dim) {
	}
	void init(int subdivides, string name) {
	}
	void read(string inputname) {
	}
	void write(string inputname) {
	}
	void optimize(double absfraction=0.2, double relfraction=0.05, int Nabsmin=1, int Nrelmin=0, int maxlevel=8, int max_level_difference=0) {}
	void grid(string name, int* gridsizes) {}
	void limit_level_difference(int max_level_difference) {
	} 
};


#include <unistd.h>
#include <getopt.h>

int option_do_init = 0;
int option_do_help = 0;
option options[] = {
	{"init", no_argument, 0, 'i'},
	{"grid", no_argument, 0, 'g'},
	{"initialize", no_argument, 0, 'i'},
	{"optimize", no_argument, 0, 't'},
	{"help", no_argument, 0, 'h'},
	{"dimension", required_argument, 0, 'd'},
	{"limit-level-difference", required_argument, 0, 'l'},
	{"subdivides", required_argument, 0, 's'},
	{"output", required_argument, 0, 'o'},
	{"input", optional_argument, 0, 'n'},
	{0}
};

extern "C" int main(int argc, char* const * argv) {
	//tuple<double, double, double> t(1,2,3);
	//tuple<double, double> t2 = get_indices<1,2>(t);
	//tuple<double, double> t(1,2);
	//std::cout << *boost::fusion::begin(t) << '\n';

	//boost::fusion::for_each(t, DoSomething());
	//return 0;

	/*cout << "test: " << log(0) << endl;
	cout << "test: " << exp(log(0)) << endl;
	cout << "test: " << isinf(log(0)) << endl;
	cout << "test: " << (log(0) < 0) << endl;
	cout << "test: " << save_exp(log(0)) << endl;*/
	int indexptr;
	int c;
	int dim = 0;
	int subdivides = 0;
	int max_level_difference = 0;
	bool init = false;
	bool optimize = false;
	bool do_grid = false;
	//bool do_limit_level_difference = false;
	string outputname;
	string inputname;
	while((c = getopt_long(argc, argv, "gn:o:ihd:s:tl:", &options[0], &indexptr)) != -1) {
		switch(c) {
			case 0:
				printf("long option: %s", options[indexptr].name);
				if(optarg)
					printf("with arg: %s", optarg);
				printf("\n");
				break;
			case 'h':
				printf("usage: %s -d <dim> [options]\n", argv[0]);
				break;
			case 'i':
				printf("initialize\n");
				init = true;
				break;
			case 'g':
				printf("grid\n");
				do_grid = true;
				break;
			case 'd':
				dim = atoi(optarg);
				printf("dimension: %i\n", dim);
				break;
			case 'l':
				max_level_difference = atoi(optarg);
				printf("maximum limit difference: %i\n", max_level_difference);
				//do_limit_level_difference = true;
				break;
			case 'o':
				outputname = string(optarg);
				printf("output: %s\n", outputname.c_str());
				break;
			case 'n':
				inputname = string(optarg);
				printf("input: %s\n", inputname.c_str());
				break;
			case 't':
				optimize = true;
				break;
			case 's':
				if(optarg == NULL) {
					fprintf(stderr, "subdivides requires int argument\n");
					exit(-1);
				}
				subdivides = atoi(optarg);
				printf("subdivide: %d\n", subdivides);
				break;
			default:
				abort();
		}
	}

	typedef action<5, double, double, double, double, double> action_highest_type;
	action_highest_type action_highest(dim);
	if(init) {
		for (int i = optind; i < argc; i++) {
			int index = i - optind; // 0 based index
			if(index >= dim) {
				//printf ("do %s\n", argv[i]);
				double value = atof(argv[i]);
				action_highest.set_point2(value, index-dim);
			} else {
				//printf ("do %s\n", argv[i]);
				double value = atof(argv[i]);
				action_highest.set_point1(value, index);
			}
		}
		/*cout << "p1 " << action_init_highest.p1 << endl; 
		cout << "p2 " << action_init_highest.p2 << endl;
		cout << "p1 " << action_init_highest.sub_action.p1 << endl; 
		cout << "p2 " << action_init_highest.sub_action.p2 << endl;
		cout << "p1s " << tuple_size<typename action_init_highest_type::Point>::value << endl;*/ 
		//action_init_highest.subdivide(dim, outputname);
		action_highest.init(subdivides, outputname);
		action_highest.write(outputname);
	}
	if(optimize) {
		if((argc-optind) != 5) {
			fprintf(stderr, "please give optimization parameters\n");
			exit(2);
		}
		double absfraction=atof(argv[optind+0]);
		double relfraction=atof(argv[optind+1]);
		int Nabsmin=atoi(argv[optind+2]);
		int Nrelmin=atoi(argv[optind+3]);
		int maxlevel=atoi(argv[optind+4]);
		cout << "absfraction, relfraction, Nabsmin, Nrelmin" << absfraction << "," <<  relfraction << "," <<  Nabsmin << "," << Nrelmin << endl;
		action_highest.read(inputname);
		action_highest.optimize(absfraction, relfraction, Nabsmin, Nrelmin, maxlevel, max_level_difference);
		/*if(do_limit_level_difference) {
			action_highest.limit_level_difference(max_level_difference);
		}*/
		action_highest.write(outputname);
	}
	if(do_grid) {
		cout << "starting gridding..." << endl;
		action_highest.read(inputname);
		cout << "starting gridding..." << endl;
		if((argc-optind) != dim) {
			fprintf(stderr, "please give gridsize with proper dimension\n");
			exit(2);
		}
		cout << "starting gridding..." << endl;
		int * gridsizes = new int[argc-optind];
		for (int i = optind; i < argc; i++) {
			int index = i - optind; // 0 based index
			gridsizes[index] = atoi(argv[i]);
			printf ("gridsizes[%d] = %d\n", index, gridsizes[index]);
		}
		cout << "starting gridding..." << endl;
		action_highest.grid(inputname, gridsizes);
		delete gridsizes;
		//action_init_highest.do_grid(dim, subdivides, inputname, outputname);
	}


	return 0;

	int count = atoi(argv[1]);
	double absfraction = atof(argv[2]);
	int Nabsmin = atoi(argv[3]);
	double relfraction = atof(argv[4]);
	int Nrelmin = atoi(argv[5]);
	cout << "Nabsmin: " << Nabsmin << endl;
	double error = 0;
	//int optimize = atoi(argv[6]);
	NLinear<1, 2, double> linear(3., 2.);
	NLinear<1, 2, double> linear2(4., 6.);
	NLinear<2, 2, double> bilinear(linear, linear2);
	cout << "integral " << linear.integrate() << ", " << linear.integrate(0., 1., 0., 1.) << endl;
	cout << "integral " << linear2.integrate() << ", " << linear2.integrate(0., 1., 0., 1.)  << endl;
	cout << "integral " << linear2.integrate(0., 1.0) << endl;
	cout << "integral " << bilinear.integrate() << endl;
	cout << "integral " << bilinear.integrate(0., 1., 0., 1.) << endl;
	cout << "integral " << bilinear.integrate(0., 0.5, 0., 0.5) << endl;
	cout << "integral " << bilinear.integrate2(0., 0.5, 0., 1., 0., 0.5, 0., 1.0) << endl;
	cout << "integral " << bilinear.integrate2(0.5, 3., 4., 10., 1., 5., 3., 8.0) << endl;
	//cout << "integral " << bilinear.integrate2(1., 5., 3., 8.0, 0.5, 3., 4., 10.) << endl;
	typedef NTree<2, double, double, double> QuadTree;
	typedef NTree<3, double, double, double, double> OcTree;
	QuadTree::Point p1(0.0,0.0);
	QuadTree::Point p2(14.0, 14.);
	set<QuadTree::Point> s;
	s.insert(p1);
	s.insert(p1);
	s.insert(p2);
	QuadTree quadtree(p1, p2);
	OcTree::Point v1(0.0,0.0,0.0);
	OcTree::Point v2(14.0,14.0,14.0);
	OcTree octree(v1, v2);
	//quadtree.add_point(0., 1.);
	//quadtree.rootNode.printinfo();
	/*quadtree.rootNode.walkrange([](double x1, double x2, double y1, double y2) {
		cout << "quad x:[" << x1 << "," << x2 << "]" << " y:[" << y1 << "," << y2 << "]" << endl;
	});
	quadtree.rootNode.walk([](double x, double y, double value) {
		cout << "point x:[" << x << "]" << " y:[" << y << "] value:"  << value << endl;
	});
	quadtree.rootNode.eval([](double x, double y) {
		return x*y;
	});
	quadtree.rootNode.walk([](double x, double y, double value) {
		cout << "point x:[" << x << "]" << " y:[" << y << "] value:"  << value << endl;
	});*/
	/*for_each(quadtree.points.begin(), quadtree.points.end(), [](QuadTree::Point p) {
		cout << "p x:[" << get<0>(p) << "]" << " y:[" << get<1>(p) << "]" << endl;
	});*/
	cout << "split" << endl;
	fflush(stdout);
	cout << "octree: # vertices: " << octree.points.size() << endl;
	octree.rootNode.split();
	//octree.rootNode.split();
	cout << "octree: # vertices: " << octree.points.size() << endl;
	cout << "quadtree: # vertices: " << quadtree.points.size() << endl;
	quadtree.rootNode.split();
	cout << "quadtree: # vertices: " << quadtree.points.size() << endl;
	
	//octree.rootNode.split();
	//quadtree.rootNode.split();
	octree.eval([](OcTree::Point p) {
		double x = get<0>(p);
		double y = get<1>(p);
		double z = get<2>(p);
		return exp(-0.5*( pow(x-2,2) + pow(y-2, 2) + pow(z-9,2)) ) * 1./(2*M_PI); 
	});

	quadtree.eval([](QuadTree::Point p) {
		double x = get<0>(p);
		double y = get<1>(p);
		return 
			exp(-0.5*( pow(x-2,2) + pow(y-2, 2) ) ) * 1./(2*M_PI);// + 
			//exp(-0.5*( pow(x-8,2) + pow(y-8, 2) ) ) * 1./(2*M_PI); 
	});
	cout << "integrated: " << quadtree.integrate(error) << endl;
	//cout << "integrated(oc): " << octree.integrate() << endl;
	for(int i = 0; i < count; i++) {
		//quadtree.rootNode.split();
		if(optimize) {
			quadtree.optimize(absfraction, relfraction, Nabsmin, Nrelmin);
			octree.optimize(absfraction, relfraction, Nabsmin, Nrelmin);
		}
		else {
			quadtree.rootNode.split();
			octree.rootNode.split();
		}
		quadtree.eval([](QuadTree::Point p) {
			double x = get<0>(p);
			double y = get<1>(p);
			return 
				exp(-0.5*( pow(x-2,2) + pow(y-2, 2) ) ) * 1./(2*M_PI);// + 
				//exp(-0.5*( pow(x-8,2) + pow(y-8, 2) ) ) * 1./(2*M_PI); 
		});
		octree.eval([](OcTree::Point p) {
			double x = get<0>(p);
			double y = get<1>(p);
			double z = get<2>(p);
			return exp(-0.5*( pow(x-2,2) + pow(y-2, 2) + pow(z-9,2)) ) * 1./(2*M_PI); 
		});
		error = 0.;
		cout << "# of points: " << quadtree.points.size() << endl;
		cout << "integrated: " << quadtree.integrate(error);
		cout << " error = " << error << endl;
		error = 0.;
		cout << "# of points(oc): " << octree.points.size() << endl;
		cout << "integrated(oc): " << octree.integrate(error);
		cout << " error = " << error << endl;
	}
/*
	quadtree.rootNode.split();
	quadtree.rootNode.left.subleft->split();
	quadtree.rootNode.left.subleft->right.subright->split();
	quadtree.rootNode.left.subleft->right.subright->left.subleft->split();
	quadtree.rootNode.left.subleft->right.subright->left.subleft->right.subright->split();
	quadtree.rootNode.split();*/
/*quadtree.rootNode.walk([](double x, double y, double value) {
		cout << "point x:[" << x << "]" << " y:[" << y << "] value:"  << value << endl;
	});*/
	/*for_each(quadtree.points.begin(), quadtree.points.end(), [](QuadTree::Point p) {
		cout << "p x:[" << get<0>(p) << "]" << " y:[" << get<1>(p) << "]" << endl;
	});*/
	quadtree.output_ps("quadtree.ps", 30., 30., 1.);
	quadtree.output_gpl("quadtree.gpl", 30., 30., 1.);
	octree.output_ps<0,1>("octree_xy.ps", 30., 30., 1.);
	octree.output_ps<0,2>("octree_xz.ps", 30., 30., 1.);
	quadtree.output("quadtree.ntree");
	/*cout << "size: " << quadtree.points.size() << endl;
	quadtree.rootNode.split();
	cout << "size: " << quadtree.points.size() << endl;
	quadtree.rootNode.split();
	cout << "size: " << quadtree.points.size() << endl;*/
	//quadtree.rootNode.left.split()
}

};