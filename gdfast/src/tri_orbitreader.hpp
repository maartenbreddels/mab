#pragma once
#include <qfile.h>
#include <qdatastream.h>
#include <vector>
#include <tr1/tuple>
#include <map>
#include <pyublas/numpy.hpp>
#include <inttypes.h>

using namespace std;

namespace gd {

typedef pyublas::numpy_vector<double> double_vector;
typedef pyublas::numpy_matrix<double> double_matrix;
typedef pyublas::numpy_vector<float> float_vector;
typedef pyublas::numpy_matrix<float> float_matrix;

class VelocityHistogram {
public:
	VelocityHistogram(int n) : n(n) {
		amplitudes = new double[n];
		fill_n(amplitudes, n, 0);
	}
	double* amplitudes;
	int n;
};

class Orbit {
public:
	Orbit(int index, int index1, int index2, int index3, int constraints, int noHistogramValues, int noMassConstraints) /*: histograms(constraints, VelocityHistogram(noHistogramValues))*/ {
		//histograms = new VelocityHistogram[constraints](constraints);
		massconstraints = new double[noMassConstraints];
	}
	vector<VelocityHistogram*> histograms;
	double maxHistogramValue;
	double* light;
	int noMassConstraints;
	double* massconstraints;
};

class Orblib {
public:
	Orblib(char* filename);
	void read();
	void readInfoOnly();
	void fillorbit(double_vector);
	void fillmassconstraints(double_vector constraintmatrix);
	char* filename;
	Q_INT32 noOrbits, noI1, noI2, noI3, noDithering;
	Q_INT32 noOrbitsPerIndex;
	Q_INT32 smom1; // what is this var?
	Q_INT32 noPhi, noTheta, noRadii;
	double *quad_lradii, *quad_lphi, *quad_ltheta;
	Q_INT32 noConstraints, noMaxVelocityHistograms;
	vector<Orbit*> orbits;
	double maxHistogramValue;
	double dvhist;
	Q_INT32 noMassConstraints;
	
	//double* massconstraints;
	double* velocityHistogram;
	Q_UINT32* orbittypes;
	typedef tr1::tuple<int, int, int> orbitidentifier;
	std::map< orbitidentifier, int > orbitToIndex;
protected:
	void readHeader(QDataStream& stream);
	void readOrbit(QDataStream& stream, int orbitIndex);
	
};

}