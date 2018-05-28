#include "aperture.hpp"
#include "integration.hpp"
#include <boost/python.hpp>
#include <stdio.h>
#include <algorithm>
#include "vectorize.hpp"


namespace gd {
using namespace boost::python;

using namespace boost::numeric::ublas;
using namespace boost::numeric;
using namespace boost;
using namespace std;

template<class ApertureStorage>
void add_aperture2d_moments(char* gridx, char* gridy, char* aperture, char* aperturestorage) {
		/*class_<GridX, boost::noncopyable >("Aperture1d+")
			.def("findindex", (&GridX::findindex))
			.def("inrange", (&GridX::inrange))
			.def("length", (&GridX::length))
			;*/
};

template<class Aperture1dType, class cc>
void add_aperture1d(const char* name) {
	class_<Aperture1dType, boost::noncopyable >(name, cc())
		.def("findindex", (&Aperture1dType::findindex))
		.def("inrange", (&Aperture1dType::inrange))
		.def("length", (&Aperture1dType::length))
		;
};

template<class Aperture2dType, class cc>
void add_aperture2d(const char* name) {
	class_<Aperture2dType, boost::noncopyable >(name, cc())
		.def("findindex", (&Aperture2dType::findindex))
		.def("inrange", (&Aperture2dType::inrange))
		.def("length", (&Aperture2dType::length))
		;
};

template<class ApertureStorage, class cc>
void add_aperture_storage(const char* name) {
	class_<ApertureStorage, boost::noncopyable >(name, cc())
		.def("bin", (&ApertureStorage::bin))
		;
};

void py_export_aperture() {
	class_<Binner>("Binner", init<double, double>())
			.def("bin", (&Binner::bin))
			.def("bin2d", (&Binner::bin2d))
			.def("bin3d", (&Binner::bin3d))
		;
	class_<ApertureGridMoments>("ApertureGridMoments", init<int_matrix, double, double, double, double, double_vector>())
			.def("bin", (&ApertureGridMoments::bin))
			.def("masses", (&ApertureGridMoments::masses))
			;
	class_<ApertureVelocityGrid>("ApertureVelocityGrid", init<int_matrix, double, double, double, double, double_vector, double>())
			.def("bin", (&ApertureVelocityGrid::bin))
			;
	class_<ApertureLikelihood1d>("ApertureLikelihood1d", init<double_matrix, double_vector, double_vector, double_vector, double, double>())
			.def("add", (&ApertureLikelihood1d::add))
			;
	class_<ApertureVelocityGridLinear>("ApertureVelocityGridLinear", init<int_vector, double, double, double>())
			.def("bin", (&ApertureVelocityGridLinear::bin))
			;
	class_<Aperture1d, boost::noncopyable >("Aperture1d")
		.def("findindex", vectorize(&Aperture1d::findindex))
		.def("inrange", (&Aperture1d::inrange))
		.def("length", (&Aperture1d::length))
		.def("uniform_transform", vectorize(&Aperture1d::uniform_transform))
		;
	class_<Grid1dRegular<>, bases<Aperture1d> >("Grid1dRegular", init<double, double, int>((boost::python::arg("x1"), boost::python::arg("x2"), boost::python::arg("length"))))
					;
	class_<Grid1dRegularLog<>, bases<Aperture1d> >("Grid1dRegularLog", init<double, double, int>((boost::python::arg("u1"), boost::python::arg("u2"), boost::python::arg("length"))))
					;
	class_<Grid1dIrregular<>, bases<Aperture1d> >("Grid1dIrregular", init<double_vector>())
			;
	/*typedef AperturePolar2d<Grid1dIrregular<NillClass>, Grid1dRegular<NillClass>> AperturePolar2d_IrrReg;
	class_<AperturePolar2d_IrrReg, boost::noncopyable>("AperturePolar2d_IrrReg", init<Grid1dIrregular<NillClass>*, Grid1dRegular<NillClass>*>())
			.def("findindex", (&AperturePolar2d_IrrReg::findindex))
			.def("length", (&AperturePolar2d_IrrReg::length))
			.def("inrange", (&AperturePolar2d_IrrReg::inrange))
	;*/
	
	class_<Aperture2d, boost::noncopyable>("Aperture2d")
			.def("findindex", (&Aperture2d::findindex))
			.def("length", (&Aperture2d::length))
			.def("inrange", (&Aperture2d::inrange))
			;
	class_<Aperture3d, boost::noncopyable>("Aperture3d")
			.def("findindex", (&Aperture3d::findindex))
			.def("length", (&Aperture3d::length))
			.def("inrange", (&Aperture3d::inrange))
			;
	class_<AperturePolar2d<>, bases<Aperture2d>, boost::noncopyable>("AperturePolar2d", init<Aperture1d*, Aperture1d*, double, double>())
			;
	class_<AperturePolar2dNoAngle<>, bases<Aperture2d>, boost::noncopyable>("AperturePolar2dNoAngle", init<Aperture1d*>())
		.def("findindex", (&AperturePolar2dNoAngle<>::findindex))
		.def("length", (&AperturePolar2dNoAngle<>::length))
		.def("inrange", (&AperturePolar2dNoAngle<>::inrange))
		;
	class_<ApertureGrid2d<>, bases<Aperture2d>, boost::noncopyable>("ApertureGrid2d", init<Aperture1d*,Aperture1d*>((boost::python::arg("gridx"), boost::python::arg("gridy"))))
		.def("findindex", (&ApertureGrid2d<>::findindex))
		.def("length", (&ApertureGrid2d<>::length))
		.def("lengthx", (&ApertureGrid2d<>::lengthx))
		.def("lengthy", (&ApertureGrid2d<>::lengthy))
		.def("inrange", (&ApertureGrid2d<>::inrange))
		;
	class_<ApertureSpherical3dNoAngle<>, bases<Aperture3d>, boost::noncopyable>("ApertureSpherical3dNoAngle", init<Aperture1d*>())
			;
	
	/*typedef ApertureStorageMoments<AperturePolar2d_IrrReg> ApertureStorageMomentsPolar_IrrReg;
	class_<ApertureStorageMomentsPolar_IrrReg, boost::noncopyable>("ApertureStorageMomentsPolar_IrrReg", init<AperturePolar2d_IrrReg*, int_vector, double_vector>())
			.def("bin", (&ApertureStorageMomentsPolar_IrrReg::bin))
	;*/
	 
	class_<ApertureStorageMoments<>, boost::noncopyable>("ApertureStorageMoments", init<Aperture2d*, int_vector, double_vector>())
			.def("bin", (&ApertureStorageMoments<>::bin))
			;
	 
	class_<ApertureStorageLosvd<>, boost::noncopyable>("ApertureStorageLosvd", init<Aperture2d*, double_vector, double_vector, double>())
			.def("bin", (&ApertureStorageLosvd<>::bin))
			.def("bin_rotate", (&ApertureStorageLosvd<>::bin_rotate))
			;
	class_<ApertureStorageLosvd2d<>, boost::noncopyable>("ApertureStorageLosvd2d", init<Aperture2d*, double_vector, double_vector, double>())
			//.def("bin", (&ApertureStorageLosvd2d<>::bin))
			.def("bin_rotate", (&ApertureStorageLosvd2d<>::bin_rotate))
			;
	class_<ApertureStorageMomentsIntrinsic<>, boost::noncopyable>("ApertureStorageMomentsIntrinsic", init<Aperture3d*, int_vector, double_vector>())
			.def("bin", (&ApertureStorageMomentsIntrinsic<>::bin))
			;
	 
	
	class_<ApertureGridMoments2d<Grid1dIrregular<>, Grid1dRegular<>> >("ApertureGridMoments2dPolarIrrigularR", init<Grid1dIrregular<>*, Grid1dRegular<>*, int_matrix, double_vector>())
			.def("bin", (&ApertureGridMoments2d<Grid1dIrregular<>, Grid1dRegular<>>::bin))
			;

	typedef Grid1dIrregular<NillClass> Grid1dIrregular_s; 
	typedef Grid1dRegular<NillClass> Grid1dRegular_s; 
	add_aperture1d<Grid1dIrregular_s, init<double_vector>>("Grid1dIrregular_s");
	add_aperture1d<Grid1dRegular_s, init<double, double, int>>("Grid1dRegular_s");

	typedef AperturePolar2d< Grid1dIrregular_s, Grid1dRegular_s, NillClass> AperturePolar2d_IrrReg_s;
	add_aperture2d<AperturePolar2d_IrrReg_s, init<Grid1dIrregular_s*, Grid1dRegular_s*, double, double>>("AperturePolar2d_IrrReg_s");

	typedef ApertureStorageMoments<AperturePolar2d_IrrReg_s> ApertureStorageMoments_Polar2d_IrrReg_s;
	add_aperture_storage<ApertureStorageMoments_Polar2d_IrrReg_s, init<AperturePolar2d_IrrReg_s*, int_vector, double_vector>>("ApertureStorageMoments_Polar2d_IrrReg_s");
	//typedef AperturePolar2d< Grid1dIrregular<NillClass>, Grid1dRegular<NillClass>, NillClass> aaa;
	//add_aperture2d_moments< ApertureStorageMoments<aaa> >("Grid1dIrregular_s", "Grid1dRegular_s", "AperturePolar2d_s_IrrReg", "ApertureStorageMoments_s_Polar2dIrrReg");
}

ApertureLikelihood1d::ApertureLikelihood1d(double_matrix likelihood_matrix, double_vector rs, double_vector vs, double_vector sigma_vs, double extra_sigma_r, double extra_sigma_v) : likelihood_matrix(likelihood_matrix), rs(rs), vs(vs), sigma_vs(sigma_vs), extra_sigma_r(extra_sigma_r), extra_sigma_v(extra_sigma_v) {
	assert(likelihood_matrix.size2() == rs.size());
	printf("sizes: %u %u\n", (int)likelihood_matrix.size1(), (int)likelihood_matrix.size2()); 
} 

void ApertureLikelihood1d::add(double_vector orbitrs, double_vector orbitvelocities, int orbitnr) {
	int stars = rs.size();
	int orbit_samples = orbitrs.size();
	assert((int)likelihood_matrix.size1() > orbitnr);
	
	
	double* starr = rs.data().begin();
	double* starv = vs.data().begin();
	double* starsigma_v = sigma_vs.data().begin();
	
	double* likelihoodp = likelihood_matrix.data().begin();
	//int likelihoodp_stride = likelihood_matrix.size2();
	int stride0 = likelihood_matrix.size2();
	likelihoodp += stride0*orbitnr; 
	
	//double* momentgridp = moment_grid.data().begin();
	//int stride0 = moment_grid.strides()[0]/moment_grid.itemsize();
	//int stride1 = moment_grid.strides()[1]/moment_grid.itemsize();
	//int stride2 = moment_grid.strides()[2]/moment_grid.itemsize();
	
	for(int i = 0; i < stars; i++) {
		double v = *starv++;
		double r = *starr++;
		double sigma_r = extra_sigma_r;
		double sigma_v = sqrt(pow(extra_sigma_v, 2) + pow(*starsigma_v++, 2)); 
		
		double* orbitr_ = orbitrs.data().begin();
		double* orbitv_ = orbitvelocities.data().begin();
		
		for(int j = 0; j < orbit_samples; j++) {
			double orbitr = *orbitr_++;
			double orbitv = *orbitv_++;
			
			double L = 1./(2.*M_PI*sigma_r*sigma_v) * exp(-0.5 *( pow((r-orbitr)/sigma_r, 2) + pow((v-orbitv)/sigma_v, 2) ) );
			*likelihoodp += L;
		}
		likelihoodp++;
	}
}


template<class Aperture>
ApertureStorageMoments<Aperture>::ApertureStorageMoments(Aperture *aperture, int_vector aperture_index_to_constraint, double_vector moment_grid) : aperture(aperture), aperture_index_to_constraint(aperture_index_to_constraint), moment_grid(moment_grid) {
}

template<class Aperture>
void ApertureStorageMoments<Aperture>::bin(double_vector x, double_vector y, double_vector velocity, int orbitnr) {
	int size = x.size();
	// momentgrid is in this order [orbitnr, moment, constaintnr]
	int* aperture_index_to_constraintp = aperture_index_to_constraint.data().begin();
	//int grid_stride = constraintgrid.size2();
	double* momentgridp = moment_grid.data().begin();
	unsigned int stride0 = moment_grid.strides()[0]/moment_grid.itemsize();
	unsigned int stride1 = moment_grid.strides()[1]/moment_grid.itemsize();
	//unsigned int stride2 = moment_grid.strides()[2]/moment_grid.itemsize();
	//printf("%d %d %d\n", stride0, stride1, stride2);
	int moments = moment_grid.dims()[1];
	double* xp = x.data().begin();
	double* yp = y.data().begin();
	double* velocityp = velocity.data().begin(); 
	//int gridwidth = constraintgrid.size1();
	//int gridheight = constraintgrid.size2();

	for(int i = 0; i < size; i++) {
		double xi = *xp++;
		double yi = *yp++;
		double v = *velocityp++;
		//printf("%f %f %f %d %d\n", xi, yi, v, i, size);
		// if position falls in grid
		if( aperture->inrange(xi, yi) ) {
			//int xgrid;// = gridx->findindex(xi);
			//int ygrid;// = gridy->findindex(yi);
			int aperture_index = aperture->findindex(xi, yi);
			//int constraintnr = aperture_index_to_constraintp[xgrid*grid_stride+ygrid];
			int constraintindex = aperture_index_to_constraintp[aperture_index];
			//int constraintindex = constraintnr -1;
			//printf("%d %d %d\n", xgrid, ygrid, constraintnr);
			// if position falls in an aperture we care about
			if(constraintindex >= 0) { // 'no value' is -1
				int offset = orbitnr*stride0+constraintindex;
				// store moments 0 always, and the higher moments as requested
				double moment = 1; // 0th moment is always 1
				momentgridp[offset+0*stride1] += moment;
				for(int j = 1; j < moments; j++) {
					moment *= v;
					momentgridp[offset+j*stride1] += moment; 
				}
			}
		}
	}
}

template<class Aperture>
ApertureStorageLosvd<Aperture>::ApertureStorageLosvd(Aperture *aperture, double_vector losvds_grid, double_vector masses_grid, double vmax) : aperture(aperture), losvds_grid(losvds_grid), masses_grid(masses_grid), vmax(vmax) {
}

template<class Aperture>
void ApertureStorageLosvd<Aperture>::bin(double_vector x, double_vector y, double_vector velocity, int orbitnr, int rotations) {
	int size = x.size();
	// losvds_grid is in this order [orbitnr, vel, r or aperture index]
	double* losvds_gridp = losvds_grid.data().begin();
	double* masses_gridp = masses_grid.data().begin();
	unsigned int stride0 = losvds_grid.strides()[0]/losvds_grid.itemsize();
	unsigned int stride1 = losvds_grid.strides()[1]/losvds_grid.itemsize();
	unsigned int masses_stride = masses_grid.strides()[0]/masses_grid.itemsize();
	int Nv = losvds_grid.dims()[1];
	int Naperture = losvds_grid.dims()[2];
	double* xp = x.data().begin();
	double* yp = y.data().begin();
	double* velocityp = velocity.data().begin(); 

	for(int i = 0; i < size; i++) {
		double xi = *xp++;
		double yi = *yp++;
		double v = *velocityp++;
		if(rotations == 0) {
			if( aperture->inrange(xi, yi) ) {
				int aperture_index = aperture->findindex(xi, yi);
				int offset = orbitnr*stride0+aperture_index;
				int vindex = (int)(((v+vmax)/(2*vmax)) * Nv);
				if((vindex >= 0) and (vindex < Nv)) { 
					losvds_gridp[offset+vindex*stride1] += 1;
				}
				masses_gridp[aperture_index+orbitnr*masses_stride] += 1;
				assert(aperture_index < Naperture);
				//assert(vindex < Nv);
			}
		}/* else {
			double xr = xi;
			double yi = yi;
			double vr = v;
			for(int j = 0; j < rotations; j++) {
				double costheta = distribution(rng_engine);
				double phi = distribution(rng_engine) * 2 * M_PI;
				double eta = distribution(rng_engine) * 2 * M_PI;
				double 
				if( aperture->inrange(xi, yi) ) {
					int aperture_index = aperture->findindex(xi, yi);
					int offset = orbitnr*stride0+aperture_index;
					int vindex = (int)(((v+vmax)/(2*vmax)) * Nv);
					if((vindex >= 0) and (vindex < Nv)) { 
						losvds_gridp[offset+vindex*stride1] += 1;
					}
					assert(aperture_index < Naperture);
					//assert(vindex < Nv);
				}
			}*/
		//}
	}
}


template<class Aperture>
void ApertureStorageLosvd<Aperture>::bin_rotate(double_vector r, double_vector vr, double_vector vt, int orbitnr, int rotations) {
	int size = r.size();
	// losvds_grid is in this order [orbitnr, vel, r or aperture index]
	double* losvds_gridp = losvds_grid.data().begin();
	double* masses_gridp = masses_grid.data().begin();
	unsigned int stride0 = losvds_grid.strides()[0]/losvds_grid.itemsize();
	unsigned int stride1 = losvds_grid.strides()[1]/losvds_grid.itemsize();
	unsigned int masses_stride = masses_grid.strides()[0]/masses_grid.itemsize();
	int Nv = losvds_grid.dims()[1];
	int Naperture = losvds_grid.dims()[2];
	double* rp = r.data().begin();
	//double* yp = y.data().begin();
	double* vrp = vr.data().begin(); 
	double* vtp = vt.data().begin(); 

	for(int i = 0; i < size; i++) {
		double ri = *rp++;
		//double yi = *yp++;
		double vr = *vrp++;
		double vt = *vtp++;
		
		// do random projections
		for(int j = 0; j < rotations; j++) {
			double costheta = distribution(rng_engine) * 2. - 1.;
			//double phi = distribution(rng_engine) * 2 * M_PI;
			double theta = acos(costheta);
			//sintheta = numpy.sqrt(1-costheta**2)
			double sintheta = sin(theta);
			double R = ri * sqrt(sintheta*sintheta);
			if( aperture->inrange(R, 0.) ) {
			//double Rmax = 1.5 * 60 * 60;
			//if(R < Rmax) {
				int aperture_index = aperture->findindex(R, 0.);
				//int aperture_index = (int)(R * Naperture / Rmax);
				int offset = orbitnr*stride0+aperture_index;
				for(int k = 0; k < 1; k++) {
					double eta = distribution(rng_engine) * 2 * M_PI;
					//double eta = random() * 2 * M_PI / RAND_MAX;
					double vlos = vr * costheta - vt * cos(eta) * sintheta;
					int vindex = (int)(((vlos+vmax)/(2*vmax)) * Nv);
					if((vindex >= 0) and (vindex < Nv)) { 
						losvds_gridp[offset+vindex*stride1] += 1;
					}
				}
				masses_gridp[aperture_index+orbitnr*masses_stride] += 1;
				assert(aperture_index < Naperture);
				//assert(vindex < Nv);
			}
		}
	}
}



template<class Aperture>
ApertureStorageLosvd2d<Aperture>::ApertureStorageLosvd2d(Aperture *aperture, double_vector losvds_grid, double_vector masses_grid, double vmax) : aperture(aperture), losvds_grid(losvds_grid), masses_grid(masses_grid), vmax(vmax) {
}



template<class Aperture>
void ApertureStorageLosvd2d<Aperture>::bin_rotate(double_vector r, double_vector vr, double_vector vt, int orbitnr, int rotations) {
	int size = r.size();
	// losvds_grid is in this order [orbitnr, vel, r or aperture index]
	double* losvds_gridp = losvds_grid.data().begin();
	double* masses_gridp = masses_grid.data().begin();
	unsigned int stride0 = losvds_grid.strides()[0]/losvds_grid.itemsize();
	unsigned int stride1 = losvds_grid.strides()[1]/losvds_grid.itemsize();
	unsigned int stride2 = losvds_grid.strides()[2]/losvds_grid.itemsize();
	unsigned int masses_stride0 = masses_grid.strides()[0]/masses_grid.itemsize();
	unsigned int masses_stride1 = masses_grid.strides()[1]/masses_grid.itemsize();
	int NLz			= losvds_grid.dims()[1];
	int Nv			= losvds_grid.dims()[2];
	int Naperture	= losvds_grid.dims()[3];
	double* rp = r.data().begin();
	//double* yp = y.data().begin();
	double* vrp = vr.data().begin(); 
	double* vtp = vt.data().begin(); 

	for(int i = 0; i < size; i++) {
		double ri = *rp++;
		//double yi = *yp++;
		double vr = *vrp++;
		double vt = *vtp++;
		
		for(int Lz = 0; Lz < NLz; Lz++) {
		// do random projections
			for(int j = 0; j < rotations; j++) {
				double costheta = 0;
				double theta = 0;
				double eta = 0;
				bool ok = false;
				while(not ok) { // ugly rejection sampling
					costheta = distribution(rng_engine) * 2. - 1.;
					theta = acos(costheta);
					eta = distribution(rng_engine) * 1 * M_PI;
					// cos alpha = Lz/L
					double cosalpha = sin(theta) * sin(eta);
					if( (cosalpha >= Lz*1./NLz ) && (cosalpha < (Lz+1)*1./NLz) ) 
						ok = true;
				}
				double phi = distribution(rng_engine) * 2 * M_PI;
				//sintheta = numpy.sqrt(1-costheta**2)
				double sintheta = sin(theta);
				double cosphi = cos(phi);
				double sinphi = sin(phi);
				//double R = ri * sqrt(sintheta*sintheta);
				double x = ri * sintheta*cos(phi);
				double y = ri * sintheta*sin(phi);
				double z = ri * costheta;
				double vx = vr * sintheta* cosphi - vt * sin(eta) * sinphi;
				//double vz = vr * costheta - vt * cos(eta) * sintheta;
				double xp = y;
				double yp = z;
				double vlos = -vx;
				
				//double vzp = vx;
				if( aperture->inrange(xp, yp) ) {
				//double Rmax = 1.5 * 60 * 60;
				//if(R < Rmax) {
					int aperture_index = aperture->findindex(xp, yp);
					//int aperture_index = (int)(R * Naperture / Rmax);
					int offset = Lz*stride1 + orbitnr*stride0+aperture_index;
					for(int k = 0; k < 1; k++) {
						//double eta = distribution(rng_engine) * 1 * M_PI;
						//double eta = random() * 2 * M_PI / RAND_MAX;
						//double vlos = vr * costheta - vt * cos(eta) * sintheta;
						int vindex = (int)(((vlos+vmax)/(2*vmax)) * Nv);
						if((vindex >= 0) and (vindex < Nv)) { 
							losvds_gridp[offset+vindex*stride2] += 1;
						}
					}
					masses_gridp[aperture_index+orbitnr*masses_stride0 + Lz*masses_stride1] += 1;
					assert(aperture_index < Naperture);
					//assert(vindex < Nv);
				}
			}
		}
	}
}

template<class Aperture>
ApertureStorageMomentsIntrinsic<Aperture>::ApertureStorageMomentsIntrinsic(Aperture *aperture, int_vector aperture_index_to_constraint, double_vector moment_grid) : aperture(aperture), aperture_index_to_constraint(aperture_index_to_constraint), moment_grid(moment_grid) {
	assert( moment_grid.dims()[1] == 13); // moment grid's 2nd dimension should be 13 to containt moments 0 to 4
}

template<class Aperture>
void ApertureStorageMomentsIntrinsic<Aperture>::bin(double_vector x, double_vector y, double_vector z, double_vector v1, double_vector v2, double_vector v3, int orbitnr) {
	int size = x.size();
	// momentgrid is in this order [orbitnr, moment, constaintnr]
	int* aperture_index_to_constraintp = aperture_index_to_constraint.data().begin();
	//int grid_stride = constraintgrid.size2();
	double* momentgridp = moment_grid.data().begin();
	unsigned int stride0 = moment_grid.strides()[0]/moment_grid.itemsize();
	unsigned int stride1 = moment_grid.strides()[1]/moment_grid.itemsize();
	//unsigned int stride2 = moment_grid.strides()[2]/moment_grid.itemsize();
	//printf("%d %d %d\n", stride0, stride1, stride2);
	//int moments = moment_grid.dims()[1];
	double* xp = x.data().begin();
	double* yp = y.data().begin();
	double* zp = z.data().begin();
	double* v1p = v1.data().begin(); 
	double* v2p = v2.data().begin(); 
	double* v3p = v3.data().begin(); 
	//int gridwidth = constraintgrid.size1();
	//int gridheight = constraintgrid.size2();

	for(int i = 0; i < size; i++) {
		double xi = *xp++;
		double yi = *yp++;
		double zi = *zp++;
		double v1value = *v1p++;
		double v2value = *v2p++;
		double v3value = *v3p++;
		double v1sq = v1value*v1value;
		double v2sq = v2value*v2value;
		double v3sq = v3value*v3value;
		//printf("%f %f %f %d %d\n", xi, yi, v, i, size);
		// if position falls in grid
		if( aperture->inrange(xi, yi, zi) ) {
			//int xgrid;// = gridx->findindex(xi);
			//int ygrid;// = gridy->findindex(yi);
			int aperture_index = aperture->findindex(xi, yi, zi);
			//int constraintnr = aperture_index_to_constraintp[xgrid*grid_stride+ygrid];
			int constraintindex = aperture_index_to_constraintp[aperture_index];
			//int constraintindex = constraintnr -1;
			//printf("%d %d %d\n", xgrid, ygrid, constraintnr);
			// if position falls in an aperture we care about
			if(constraintindex >= 0) { // 'no value' is -1
				int offset = orbitnr*stride0+constraintindex;
				// store moments 0 always, and the higher moments as requested
				//double moment = 1; // 0th moment is always 1
				
				momentgridp[offset+0*stride1] += 1; // 0th moment
				momentgridp[offset+1*stride1] += v1value; // 1st moment
				momentgridp[offset+2*stride1] += v2value;
				momentgridp[offset+3*stride1] += v3value;
				momentgridp[offset+4*stride1] += v1sq; // 2nd moment
				momentgridp[offset+5*stride1] += v2sq;
				momentgridp[offset+6*stride1] += v3sq;
				momentgridp[offset+7*stride1] += v1sq*v1value; // 3rd moment
				momentgridp[offset+8*stride1] += v2sq*v2value;
				momentgridp[offset+9*stride1] += v3sq*v3value;
				momentgridp[offset+10*stride1] += v1sq*v1sq; // 4th moment
				momentgridp[offset+11*stride1] += v2sq*v2sq;
				momentgridp[offset+12*stride1] += v3sq*v3sq;
			}
		}
	}
}


template<class GridX, class GridY>
ApertureGridMoments2d<GridX, GridY>::ApertureGridMoments2d(GridX* gridx, GridY* gridy, int_matrix constraintgrid, double_vector moment_grid) : gridx(gridx), gridy(gridy), constraintgrid(constraintgrid), moment_grid(moment_grid) {
	
}

template<class GridX, class GridY>
void ApertureGridMoments2d<GridX, GridY>::bin(double_vector x, double_vector y, double_vector velocity, int orbitnr) {
	int size = x.size();
	// momentgrid is in this order [orbitnr, moment, constaintnr]
	int* gridp = constraintgrid.data().begin();
	int grid_stride = constraintgrid.size2();
	double* momentgridp = moment_grid.data().begin();
	int stride0 = moment_grid.strides()[0]/moment_grid.itemsize();
	int stride1 = moment_grid.strides()[1]/moment_grid.itemsize();
	//int stride2 = moment_grid.strides()[2]/moment_grid.itemsize();
	//printf("%d %d %d\n", stride0, stride1, stride2);
	int moments = moment_grid.dims()[1];
	double* xp = x.data().begin();
	double* yp = y.data().begin();
	double* velocityp = velocity.data().begin(); 
	//int gridwidth = constraintgrid.size1();
	//int gridheight = constraintgrid.size2();
	
	for(int i = 0; i < size; i++) {
		double xi = *xp++;
		double yi = *yp++;
		double v = *velocityp++;
		//printf("%f %f %f %d %d\n", xi, yi, v, i, size);
		// if position falls in grid
		if( gridx->inrange(xi) && gridy->inrange(yi) ) {
			int xgrid = gridx->findindex(xi);
			int ygrid = gridy->findindex(yi);
			// constraint # start from 1, but indices from 0
			int constraintnr = gridp[xgrid*grid_stride+ygrid];
			int constraintindex = constraintnr -1;
			//printf("%d %d %d\n", xgrid, ygrid, constraintnr);
			// if position falls in an aperture we care about
			if(constraintindex >= 0) { // 'no value' is 0 for constraint, so -1 for index
				int offset = orbitnr*stride0+constraintindex;
				// store moments 0 always, and the higher moments as requested
				double moment = 1; // 0th moment is always 1
				momentgridp[offset+0*stride1] += moment;
				for(int j = 1; j < moments; j++) {
					moment *= v;
					momentgridp[offset+j*stride1] += moment; 
				}
			}
		}
	}
}


ApertureGridMoments::ApertureGridMoments(int_matrix grid, double x0, double y0, double width, double height, double_vector moment_grid) : grid(grid), x0(x0), y0(y0), width(width), height(height), moment_grid(moment_grid) {
	maxbin = *std::max_element(grid.data().begin(), grid.data().end());
	//printf("size: %d %d, max: %d", moment_grid.size1(), moment_grid.size2(), maxbin);
	assert(moment_grid.dims()[2] == maxbin);
	assert(width == height);
	assert(grid.size1() == grid.size2());
	
}





double_vector ApertureGridMoments::masses(Galaxy* galaxy) {
	double_vector masses(maxbin);
	Profile *profile = galaxy->getStellarProfile();
	double kpc_per_arcsec = 1./galaxy->kpc_to_arcsec(1.);
	for(unsigned int xi = 0; xi < grid.size1(); xi++) {
		for(unsigned int yi = 0; yi < grid.size2(); yi++) {
			int constraint = grid(xi, yi);
			if(constraint > 0) { 
				double x1 = x0 + xi * width/grid.size1();
				double y1 = y0 + yi * height/grid.size2();
				double x2 = x0 + (xi+1) * width/grid.size1();
				double y2 = y0 + (yi+1) * height/grid.size2();
				// double integral 
				auto dmassx = [&](double x){
					auto ddmassxy = [&](double y) {
						return profile->densityR(sqrt(x*x+y*y)*kpc_per_arcsec);
					};
					IntegratorGSL<> integratorGSL2(ddmassxy);
					return integratorGSL2.integrate(y1, y2);
				};
				IntegratorGSL<> integratorGSL(dmassx);
				double mass = integratorGSL.integrate(x1, x2);
				masses(constraint-1) += mass * ( kpc_per_arcsec * kpc_per_arcsec);
			}
			
		}
	}
	return masses;
}

void ApertureGridMoments::bin(double_vector x, double_vector y, double_vector velocity, int orbitnr) {
	int size = x.size();
	// momentgrid is in this order [orbitnr, moment, constaintnr]
	int* gridp = grid.data().begin();
	int grid_stride = grid.size2();
	double* momentgridp = moment_grid.data().begin();
	int stride0 = moment_grid.strides()[0]/moment_grid.itemsize();
	int stride1 = moment_grid.strides()[1]/moment_grid.itemsize();
	//int stride2 = moment_grid.strides()[2]/moment_grid.itemsize();
	//printf("%d %d %d\n", stride0, stride1, stride2);
	int moments = moment_grid.dims()[1];
	double* xp = x.data().begin();
	double* yp = y.data().begin();
	double* velocityp = velocity.data().begin(); 
	int gridwidth = grid.size1();
	int gridheight = grid.size2();
	
	for(int i = 0; i < size; i++) {
		double xi = *xp++;
		double yi = *yp++;
		double v = *velocityp++;
		int xgrid = (int)((xi - x0)/width*gridwidth);
		int ygrid = (int)((yi - y0)/height*gridheight);
		//printf("%f %f %f %d %d\n", xi, yi, v, i, size);
		// if position falls in grid
		if( (xgrid >= 0) && (ygrid >= 0) && (xgrid < gridwidth) && (ygrid < gridheight) ) {
			// constraint # start from 1, but indices from 0
			int constraintnr = gridp[xgrid*grid_stride+ygrid];
			int constraintindex = constraintnr -1;
			//printf("%d %d %d\n", xgrid, ygrid, constraintnr);
			// if position falls in an aperture we care about
			if(constraintindex >= 0) { // 'no value' is 0 for constraint, so -1 for index
				int offset = orbitnr*stride0+constraintindex;
				// store moments 0 always, and the higher moments as requested
				double moment = 1; // 0th moment is always 1
				momentgridp[offset+0*stride1] += moment;
				for(int j = 1; j < moments; j++) {
					moment *= v;
					momentgridp[offset+j*stride1] += moment; 
				}
			}
		}
	}
}


ApertureVelocityGrid::ApertureVelocityGrid(int_matrix grid, double x0, double y0, double width, double height, double_vector velocity_grid, double vmax) : grid(grid), x0(x0), y0(y0), width(width), height(height), velocity_grid(velocity_grid), vmax(vmax) {
	int maxbin = *std::max_element(grid.data().begin(), grid.data().end());
	//printf("size: %d %d, max: %d", moment_grid.size1(), moment_grid.size2(), maxbin);
	assert(velocity_grid.dims()[2] == maxbin);
	assert(width == height);
	assert(grid.size1() == grid.size2());
	
}

void ApertureVelocityGrid::bin(double_vector x, double_vector y, double_vector velocity, int orbitnr) {
	int size = x.size();
	// momentgrid is in this order [orbitnr, moment, constaintnr]
	int* gridp = grid.data().begin();
	int grid_stride = grid.size2();
	double* velocitygridp = velocity_grid.data().begin();
	int stride0 = velocity_grid.strides()[0]/velocity_grid.itemsize();
	int stride1 = velocity_grid.strides()[1]/velocity_grid.itemsize();
	//int stride2 = velocity_grid.strides()[2]/velocity_grid.itemsize();
	//printf("%d %d %d\n", stride0, stride1, stride2);
	int vbins = velocity_grid.dims()[1];
	double* xp = x.data().begin();
	double* yp = y.data().begin();
	double* velocityp = velocity.data().begin(); 
	int gridwidth = grid.size1();
	int gridheight = grid.size2();
	
	for(int i = 0; i < size; i++) {
		double xi = *xp++;
		double yi = *yp++;
		double v = *velocityp++;
		int xgrid = (int)((xi - x0)/width*gridwidth);
		int ygrid = (int)((yi - y0)/height*gridheight);
		//printf("%f %f %f %d %d\n", xi, yi, v, i, size);
		// if position falls in grid
		if( (xgrid >= 0) && (ygrid >= 0) && (xgrid < gridwidth) && (ygrid < gridheight) ) {
			// constraint # start from 1, but indices from 0
			int constraintnr = gridp[xgrid*grid_stride+ygrid];
			int constraintindex = constraintnr -1;
			//printf("%d %d %d\n", xgrid, ygrid, constraintnr);
			// if position falls in an aperture we care about
			if(constraintindex >= 0) { // 'no value' is 0 for constraint, so -1 for index
				int offset = orbitnr*stride0+constraintindex;
				// v
				int vindex = (int)((v+vmax)/(2*vmax) * vbins);
				// if velocity falls in range [-vmax, vmax], add '1'
				if((vindex >= 0) && (vindex < vbins))
					velocitygridp[offset+vindex*stride1] += 1; 
			}
		}
	}
}



ApertureVelocityGridLinear::ApertureVelocityGridLinear(int_vector velocity_grid, double Rmin, double Rmax, double vmax) : velocity_grid(velocity_grid), Rmin(Rmin), Rmax(Rmax), vmax(vmax) {
	//int maxbin = *std::max_element(grid.data().begin(), grid.data().end());
	//printf("size: %d %d, max: %d", moment_grid.size1(), moment_grid.size2(), maxbin);
	/*assert(velocity_grid.dims()[2] == maxbin);
	assert(width == height);
	assert(grid.size1() == grid.size2());*/
	
}

void ApertureVelocityGridLinear::bin(double_vector Rs, double_vector los_velocities, int orbitnr) {
	int size = Rs.size();
	// velocity_grid is in this order [orbitnr, vindex, Rindex]
	int* velocitygridp = velocity_grid.data().begin();
	int stride0 = velocity_grid.strides()[0]/velocity_grid.itemsize();
	int stride1 = velocity_grid.strides()[1]/velocity_grid.itemsize();
	//int stride2 = velocity_grid.strides()[2]/velocity_grid.itemsize();
	//printf("%d %d %d\n", stride0, stride1, stride2);
	int vbins = velocity_grid.dims()[1];
	int Rbins = velocity_grid.dims()[2];
	
	double* Rp = Rs.data().begin();
	double* velocityp = los_velocities.data().begin(); 
	
	//double* yp = y.data().begin();
	//int gridwidth = grid.size1();
	//int gridheight = grid.size2();
	
	for(int i = 0; i < size; i++) {
		double R = *Rp++;
		double v = *velocityp++;
		int Rindex = (int)((R - Rmin)/(Rmax-Rmin)*Rbins);
		// if position falls in grid
		if( (Rindex >= 0) && (Rindex < Rbins) ) {
			int offset = orbitnr*stride0+Rindex;
			int vindex = (int)((v+vmax)/(2*vmax) * vbins);
			//printf("%f %f %d %d ", R, v, Rindex, Rbins);
			//printf("%d %d %d\n", offset, vindex, vbins);
			if((vindex >= 0) && (vindex < vbins))
				velocitygridp[offset+vindex*stride1] += 1; 
		}
	}
}


}
