#pragma once
#include "datatypes.hpp"
#include "galaxy.hpp"
#include <math.h>
#include <random>


namespace gd {
	
void py_export_aperture();
/*! \addtogroup apertures
*  @{
*/
/**

more about apertures...
*/

/**@} */


//! Aperture that stores likelihoods
/** 

This aperture stores the likelihoods...
*/
class ApertureLikelihood1d {
public:
	ApertureLikelihood1d(double_matrix likelihood_matrix, double_vector rs, double_vector vs, double_vector sigma_vs, double extra_sigma_r, double extra_sigma_v);
	void add(double_vector r, double_vector velocity, int orbitnr);
	double_matrix likelihood_matrix;
	double_vector rs, vs, sigma_vs;
	double extra_sigma_r;
	double extra_sigma_v;
	//double_vector masses(Galaxy* profile);
	//int_matrix grid;
	//double_vector moment_grid;
	//double x0, y0, width, height;
	//int maxbin;
};



class NillClass {};
class Aperture1d {
public:
	virtual int findindex(double) { return 0; }
	virtual bool inrange(double) { return false; }
	virtual int length() {return 0;}
	virtual	double uniform_transform(double) {
			return 0;
		}
};

template<class Base=Aperture1d>
class Grid1dIrregular : public Base {
	public:
		Grid1dIrregular(double_vector rborders) : rborders(rborders), index(0) {
			rbordersp = rborders.data().begin();
			border_length = rborders.size();
		}
		int findindex(double r) {
			if( (r >= rbordersp[index]) && (r < rbordersp[index+1]))
				return index;
			if (r < rbordersp[0])
				return -1;
			if (r >= rbordersp[border_length-1])
				return -1;
	
			if(r < rbordersp[index]) // we should go to the left 
				while(r < rbordersp[index]) // if r >= rbordersp[index], we found it
					index--;
			else
				while(r >= rbordersp[index+1]) // if r < rbordersp[index+1], we found it
					index++;
			return index;
		}
		bool inrange(double r) {
			return (r >= rbordersp[0]) && (r < rbordersp[border_length-1]); 
		}
		int length() { return border_length-1; } /* n bins, need n+1 borders */
	private:
		double_vector rborders;
		int border_length;
		int index; // remember the old position (this is slow for random r, but slow for slowly changing r), index refers to the 'left' edge
		double* rbordersp;
};

template<class Base=Aperture1d>
class Grid1dRegular : public Base {
	public:
		Grid1dRegular(double x1, double x2, int length) : x1(x1), x2(x2), _length(length) {}
		int findindex(double x) {
			
			return x < x1 ? 0 :
					(x >= x2 ? _length-1 : (int)((x-x1)/(x2-x1)*(_length)));
		}
		bool inrange(double r) {
			return (r >= x1) && (r < x2); 
		}
		int length() { return _length; }
		double x1, x2;
		int _length;
};
template<class Base=Aperture1d>
class Grid1dRegularLog : public Base {
	public:
		Grid1dRegularLog(double u1, double u2, int length) : u1(u1), u2(u2), _length(length) {}
		int findindex(double x) {
			double u = log10(x);
			return u < u1 ? 0 :
					(u >= u2 ? -1 : (int)((u-u1)/(u2-u1)*(_length)));
		}
		bool inrange(double x) {
			double u = log10(x);
			return (u >= u1) && (u < u2); 
		}
		double uniform_transform(double uniform) {
			return pow(10, uniform * (u2-u1) + u1);
		}
		int length() { return _length; }
		double u1, u2;
		int _length;
};


class Aperture2d {
public:
	virtual int findindex(double, double) { return 0; }
	virtual bool inrange(double, double) { return false; }
	virtual int length() {return 0;}
};

class Aperture3d {
public:
	virtual int findindex(double, double, double) { return 0; }
	virtual bool inrange(double, double, double) { return false; }
	virtual int length() {return 0;}
};

template<class GridX=Aperture1d, class GridY=Aperture1d, class Base=Aperture2d>
class AperturePolar2d : public Base {
public:
	typedef GridX GridXType;
	typedef GridY GridYType;
	AperturePolar2d(GridX* gridx, GridY* gridy, double ellipticity, double position_angle) : gridx(gridx), gridy(gridy), ellipticity(ellipticity), position_angle(position_angle) {}
	void findindex2d(double x, double y, int& xi, int &yi) {
			double xe = (x * cos(position_angle) + y * sin(position_angle));
			double ye = (-x * sin(position_angle) + y * cos(position_angle))/(1 - ellipticity);
			double re = sqrt(xe*xe + ye*ye);
			//double r = sqrt(x*x+y*y);
			double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI); 
			xi = gridx->findindex(re);
			yi = gridy->findindex(phi);
	} 
	int findindex(double x, double y) {
		double xe = (x * cos(position_angle) + y * sin(position_angle));
		double ye = (-x * sin(position_angle) + y * cos(position_angle))/(1 - ellipticity);
		double re = sqrt(xe*xe + ye*ye);
		//double r = sqrt(x*x+y*y);
		double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI); 
		return gridx->findindex(re) * gridy->length() + gridy->findindex(phi);;
	} 
	bool inrange(double x, double y) {
		double xe = (x * cos(position_angle) + y * sin(position_angle));
		double ye = (-x * sin(position_angle) + y * cos(position_angle))/(1 - ellipticity);
		double re = sqrt(xe*xe + ye*ye);
		//double r = sqrt(x*x+y*y);
		double phi = fmod(atan2(y,x)+2*M_PI, 2*M_PI);
		return gridx->inrange(re) && gridy->inrange(phi);
	}
	int length() {
		return gridx->length() * gridy->length();
	}
	GridX* gridx;
	GridY* gridy;
	double ellipticity;
	double position_angle;
};

template<class GridX=Aperture1d, class GridY=Aperture1d, class Base=Aperture2d>
class ApertureGrid2d : public Base {
public:
	typedef GridX GridXType;
	typedef GridY GridYType;
	ApertureGrid2d(GridX* gridx, GridY* gridy) : gridx(gridx), gridy(gridy) {}
	void findindex2d(double x, double y, int& xi, int &yi) {
		xi = gridx->findindex(x);
		yi = gridy->findindex(y);
	} 
	int findindex(double x, double y) {
		return gridx->findindex(x) * gridy->length() + gridy->findindex(y);;
	} 
	bool inrange(double x, double y) {
		return gridx->inrange(x) && gridy->inrange(y);
	}
	int length() {
		return gridx->length() * gridy->length();
	}
	int lengthx() {
		return gridx->length();
	}
	int lengthy() {
		return gridy->length();
	}
	GridX* gridx;
	GridY* gridy;
};

template<class GridR=Aperture1d, class Base=Aperture2d>
class AperturePolar2dNoAngle : public Base {
public:
	typedef GridR GridRType;
	AperturePolar2dNoAngle(GridR* gridr) : gridr(gridr) {}
	void findindex2d(double x, double y, int& xi, int &yi) {
		double r = sqrt(x*x+y*y);
		xi = gridr->findindex(r);
		yi = 0;
	} 
	int findindex(double x, double y) {
		double r = sqrt(x*x+y*y);
		return gridr->findindex(r);
	} 
	bool inrange(double x, double y) {
		double r = sqrt(x*x+y*y);
		return gridr->inrange(r);
	}
	int length() {
		return gridr->length();
	}
	GridR* gridr;
};

template<class GridR=Aperture1d, class Base=Aperture3d>
class ApertureSpherical3dNoAngle : public Base {
public:
	typedef GridR GridRType;
	ApertureSpherical3dNoAngle(GridR* gridr) : gridr(gridr) {}
	int findindex(double x, double y, double z) {
		double r = sqrt(x*x+y*y+z*z);
		return gridr->findindex(r);
	} 
	bool inrange(double x, double y, double z) {
		double r = sqrt(x*x+y*y+z*z);
		return gridr->inrange(r);
	}
	int length() {
		return gridr->length();
	}
	GridR* gridr;
};

template<class Aperture=Aperture2d>
class ApertureStorageMoments {
	public:
		typedef Aperture ApertureType;
		ApertureStorageMoments(Aperture *aperture, int_vector aperture_index_to_constraint, double_vector moment_grid);
		void bin(double_vector x, double_vector y, double_vector los_velocities, int orbitnr);
		Aperture *aperture;
		int_vector aperture_index_to_constraint;
		double_vector moment_grid;
};


template<class Aperture=Aperture2d>
class ApertureStorageLosvd {
	public:
		typedef Aperture ApertureType;
		ApertureStorageLosvd(Aperture *aperture, double_vector losvds_grid, double_vector masses_grid, double vmax);
		void bin(double_vector x, double_vector y, double_vector los_velocities, int orbitnr, int rotations=0);
		void bin_rotate(double_vector x, double_vector y, double_vector los_velocities, int orbitnr, int rotations=1);
		void seed(int s);
		Aperture *aperture;
		double_vector losvds_grid;
		double_vector masses_grid;
		double vmax;
		std::uniform_real_distribution<double> distribution;
		std::mt19937 rng_engine;
		//int rotations;
};


template<class Aperture=Aperture2d>
class ApertureStorageLosvd2d {
	public:
		typedef Aperture ApertureType;
		ApertureStorageLosvd2d(Aperture *aperture, double_vector losvds_grid, double_vector masses_grid, double vmax);
		//void bin(double_vector x, double_vector y, double_vector los_velocities, int orbitnr, int rotations=0);
		void bin_rotate(double_vector x, double_vector y, double_vector los_velocities, int orbitnr, int rotations=1);
		void seed(int s);
		Aperture *aperture;
		double_vector losvds_grid;
		double_vector masses_grid;
		double vmax;
		std::uniform_real_distribution<double> distribution;
		std::mt19937 rng_engine;
		//int rotations;
};

template<class GridX, class GridY>
class ApertureGridMoments2d {
	public:
		ApertureGridMoments2d(GridX* gridx, GridY* gridy, int_matrix constraintgrid, double_vector moment_grid);
		void bin(double_vector x, double_vector y, double_vector los_velocities, int orbitnr);
		GridX* gridx;
		GridY* gridy;
		int_matrix constraintgrid;
		double_vector moment_grid;
};


class ApertureVelocityGridLinear {
	public:
		ApertureVelocityGridLinear(int_vector velocity_grid, double Rmin, double Rmax, double vmax);
		void bin(double_vector Rs, double_vector los_velocities, int orbitnr);
		int_vector velocity_grid;
		double Rmin, Rmax, vmax; 
};


class ApertureGridMoments {
public:
	ApertureGridMoments(int_matrix grid, double x0, double y0, double width, double height, double_vector moment_grid);
	void bin(double_vector x, double_vector y, double_vector velocity, int orbitnr);
	double_vector masses(Galaxy* profile);
	int_matrix grid;
	double x0, y0, width, height;
	double_vector moment_grid;
	int maxbin;
};

class ApertureVelocityGrid {
public:
	ApertureVelocityGrid(int_matrix grid, double x0, double y0, double width, double height, double_vector velocity_grid, double vmax);
	void bin(double_vector x, double_vector y, double_vector velocity, int orbitnr);
	int_matrix grid;
	double x0, y0, width, height; 
	double_vector velocity_grid;
	double vmax;
};

class Binner {
public:
	Binner(double xmin, double xmax) : xmin(xmin), xmax(xmax) {
	}
	void bin(double_vector x, double_vector data, double_vector binned) {
		double* xp = x.data().begin();
		double* datap = data.data().begin();
		double* binnedp = binned.data().begin();
		int size = x.size();
		int bincount = binned.size();
		for(int i=0; i < size; i++) {
			double xi = xp[i];
			if((xi >= xmin) && (xi < xmax)) {
				int index = (int)((xi-xmin)/(xmax-xmin)*bincount);
				binnedp[index] += datap[i];
			}
		}
		
	}
	void bin2d(double_vector x, double_vector data1, double_vector data2, double_matrix momentgrid, int_vector pow1, int_vector pow2) {
		double* xp = x.data().begin();
		double* data1p = data1.data().begin();
		double* data2p = data2.data().begin();
		int* pow1p = pow1.data().begin();
		int* pow2p = pow2.data().begin();
		double* momentgridp = momentgrid.data().begin();
		int stride = momentgrid.size2();
		int moments = momentgrid.size1();
		assert(moments == (int)pow1.size());
		assert(moments == (int)pow2.size());
		assert(x.size() == data1.size());
		assert(x.size() == data2.size());
		int size = x.size();
		int bincount = momentgrid.size2();
		for(int i=0; i < size; i++) {
			double xi = *xp++;
			double data1i = *data1p++;
			double data2i = *data2p++;
			if((xi >= xmin) && (xi < xmax)) {
				int index = (int)((xi-xmin)/(xmax-xmin)*bincount);
				if( (index >= 0) && (index < bincount)) {
					for(int m = 0; m < moments; m++)
						momentgridp[index+m*stride] += pow(data1i, pow1p[m]) * pow(data2i, pow2p[m]);
				}
			}
		}
		
	}	
	void bin3d(double_vector x, double_vector x1, double_vector x2, double_vector x3, double_matrix momentgrid) {
		double* xp = x.data().begin();
		double* x1p = x1.data().begin();
		double* x2p = x2.data().begin();
		double* x3p = x3.data().begin();
		double* momentgridp = momentgrid.data().begin();
		int stride = momentgrid.size2();
		int moments = momentgrid.size1();
		assert(moments == (1+3+3));
		assert(x.size() == x1.size());
		assert(x.size() == x2.size());
		assert(x.size() == x3.size());
		int size = x.size();
		int bincount = momentgrid.size2();
		for(int i=0; i < size; i++) {
			double xi = *xp++;
			double x1i = *x1p++;
			double x2i = *x2p++;
			double x3i = *x3p++;
			if((xi >= xmin) && (xi < xmax)) {
				int index = (int)((xi-xmin)/(xmax-xmin)*bincount);
				if( (index >= 0) && (index < bincount)) {
					momentgridp[index+0*stride] += 1;
					momentgridp[index+1*stride] += x1i;
					momentgridp[index+2*stride] += x2i;
					momentgridp[index+3*stride] += x3i;
					momentgridp[index+4*stride] += x1i*x1i;
					momentgridp[index+5*stride] += x2i*x2i;
					momentgridp[index+6*stride] += x3i*x3i;
				}
			}
		}
		
	}
private:
	double xmin, xmax;
};

template<class Aperture=Aperture3d>
class ApertureStorageMomentsIntrinsic {
	public:
		typedef Aperture ApertureType;
		ApertureStorageMomentsIntrinsic(Aperture *aperture, int_vector aperture_index_to_constraint, double_vector moment_grid);
		void bin(double_vector x, double_vector y, double_vector z, double_vector v1, double_vector v2, double_vector v3, int orbitnr);
		Aperture *aperture;
		int_vector aperture_index_to_constraint;
		double_vector moment_grid;
};

	
}


