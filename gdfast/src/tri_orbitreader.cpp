#include "tri_orbitreader.hpp"

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <tr1/functional>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/ref.hpp>
#include <boost/python.hpp>



namespace gd {
//using namespace boost::python;

void py_export_tri_orbitreader() {
	boost::python::class_<Orblib >("Orblib", boost::python::init<char*>())
			.def("read", (&Orblib::read))
			.def("readInfoOnly", (&Orblib::readInfoOnly))
			.def("fillorbit", (&Orblib::fillorbit))
			.def("fillmassconstraints", &Orblib::fillmassconstraints)
			.def_readonly("noI1", &Orblib::noI1)
			.def_readonly("noI2", &Orblib::noI2)
			.def_readonly("noI3", &Orblib::noI3)
			.def_readonly("noPhi", &Orblib::noPhi)
			.def_readonly("noTheta", &Orblib::noTheta)
			.def_readonly("noRadii", &Orblib::noRadii)
			.def_readonly("noMassConstraints", &Orblib::noMassConstraints)
			.def_readonly("noMaxVelocityHistograms", &Orblib::noMaxVelocityHistograms)
			.def_readonly("noConstraints", &Orblib::noConstraints)
			.def_readonly("dvhist", &Orblib::dvhist)
			;

}
using namespace std;
//using namespace std::tr1;
//using namespace std::tr1::placeholders;
using namespace boost;
using namespace boost::lambda;
//using namespace boost::numeric::ublas;


struct QDataStreamReader {
	QDataStream& stream;
	QDataStreamReader(QDataStream& stream) : stream(stream) {}
	template<typename T>
	void operator()(T& t) {
		stream >> t;
	}
};


Orblib::Orblib(char* filename) : filename(filename), quad_lradii(0), quad_lphi(0), quad_ltheta(0), velocityHistogram(0), orbittypes(0) {
}

void Orblib::fillmassconstraints(double_vector constraintmatrix) {
	for(int orbitIndex = 0; orbitIndex < noOrbits; orbitIndex++) {
		Orbit* orbit = orbits[orbitIndex];
		for(int r = 0; r < noRadii; r++)
			for(int theta = 0; theta < noTheta; theta++)
				for(int phi = 0; phi < noPhi; phi++)
					for(int i = 0; i < 16; i++) {
						//constraintmatrix.sub(orbitIndex, r, theta, phi, i) = orbit->massconstraints[r+noRadii*(theta + noTheta*(phi + noPhi*i))];
						/*if(i>0)
							constraintmatrix.sub(orbitIndex, r, theta, phi, i) = orbit->massconstraints[i+16*(phi + noPhi*(theta + noTheta*r))] * orbit->massconstraints[16*(phi + noPhi*(theta + noTheta*r))];
						else
							constraintmatrix.sub(orbitIndex, r, theta, phi, i) = orbit->massconstraints[i+16*(phi + noPhi*(theta + noTheta*r))];
						*/
					}
	} 
}
void Orblib::fillorbit(double_vector orbitmatrix) {
	for(int orbitIndex = 0; orbitIndex < noOrbits; orbitIndex++) {
		Orbit* orbit = orbits[orbitIndex];
		for(int constraint = 0; constraint < noConstraints; constraint++) {
			VelocityHistogram* histogram = orbit->histograms[constraint];
			for(int index = 0; index < noMaxVelocityHistograms; index++) {
				double value = histogram->amplitudes[index];
				orbitmatrix.sub(orbitIndex, constraint, index) = value;
			} 
		}
	}
}

void Orblib::read() {
	QFile file(filename);
	file.open(IO_ReadOnly);
	QDataStream stream(&file);
	stream.setByteOrder(QDataStream::LittleEndian);
	readHeader(stream);
	
	maxHistogramValue = 0;
	velocityHistogram = new double[noMaxVelocityHistograms*2+1];
	orbittypes = new Q_UINT32[noOrbitsPerIndex];
	
	for(int orbitIndex = 0; orbitIndex < noOrbits; orbitIndex++) {
		readOrbit(stream, orbitIndex+1);
	}
	//readOrbit(stream, 1);
	/*for(int i = 0; i < 10; i++) {
		Q_UINT32 t;
		stream >> t;
		cout << "t = " << t << endl;
			//cout << file.getch() << endl;
}*/
	file.close();
}


void Orblib::readInfoOnly() {
	QFile file(filename);
	file.open(IO_ReadOnly);
	QDataStream stream(&file);
	stream.setByteOrder(QDataStream::LittleEndian);
	readHeader(stream);
	
	file.close();
}

void Orblib::readHeader(QDataStream& stream) {
	Q_UINT32 size1, size2;
	
	// part
	stream >> size1;
	stream >> noOrbits >> noI1 >> noI2 >> noI3 >> noDithering;
	cout << "number of orbits: " << noOrbits << endl;
	cout << "number of I1-I3: " << noI1 << ", " << noI2 << ", " << noI3 << endl;
	noOrbitsPerIndex = pow(noDithering, 3);
	cout << "number of dithering: " << noDithering << "^3 = " << noOrbitsPerIndex << endl;
	stream >> size2;
	ASSERT(size1 == size2);
	
	// part
	stream >> size1 >> smom1 >> noPhi >> noTheta >> noRadii >> size2;
	cout << "number of grids in the r, phi, theta direction: " << noRadii << ", " << noPhi << ", " << noTheta << "(smom1 = " << smom1 << ")" << endl;
	ASSERT(size1 == size2);
	noMassConstraints = smom1*noPhi*noTheta*noRadii;
		
	quad_lradii = new double[noRadii + 1];
	quad_lphi = new double[noPhi + 1];
	quad_ltheta = new double[noTheta + 1];
	
	// part
	stream >> size1;
	//cout << "size1 = " << size1 << endl;
	QDataStreamReader sr(stream);
	//for_each(quad_lradii, quad_lradii+noRadii + 1, cout << _1 << '\n');
	for_each(quad_lradii, quad_lradii+noRadii + 1, sr); //(*s2) >> _1);
	//for_each(quad_lradii, quad_lradii+noRadii + 1, cout<< constant("r = ")  << _1 << '\n');
	stream >> size2;
	//cout << "size2 = " << size2 << endl;
	ASSERT(size1 == size2);
		
	// part
	stream >> size1;
	//cout << "size1 = " << size1 << endl;
	for_each(quad_ltheta, quad_ltheta+noTheta+ 1, sr); //(*s2) >> _1);
	for_each(quad_ltheta, quad_ltheta+noTheta + 1, cout << constant("theta = ") << boost::lambda::_1 << '\n');
	stream >> size2;
	//cout << "size2 = " << size2 << endl;
	ASSERT(size1 == size2);
		
	// part
	stream >> size1;
	//cout << "size1 = " << size1 << endl;
	for_each(quad_lphi, quad_lphi+noPhi+ 1, sr); //(*s2) >> _1);
	for_each(quad_lphi, quad_lphi+noPhi + 1, cout << constant("phi = ") << boost::lambda::_1 << '\n');
	stream >> size2;
	//cout << "size2 = " << size2 << endl;
	ASSERT(size1 == size2);
		
	//noConstraints, noMaxVelocityHistograms;
	/*
	t1 = hist_basic(1,3)/2.0_sp ! corrected by Remco 20/JAN/2003
	write (unit=handle) h_nconstr, t1, hist_basic(1,1)/hist_basic(1,3)
	*/
	stream >> size1 >> noConstraints >> noMaxVelocityHistograms >> dvhist >> size2;
	noMaxVelocityHistograms = noMaxVelocityHistograms * 2 + 1; // TODO no idea why this is stored like this 
	cout << "number of constraints: " <<  noConstraints << endl;
	cout << "noMaxVelocityHistograms, dvhist: " <<  noMaxVelocityHistograms << ", " << dvhist << endl;
	ASSERT(size1 == size2);
}
void Orblib::readOrbit(QDataStream& stream, int orbitIndex){
	QDataStreamReader sr(stream);
	Q_INT32 size1, size2;
	Q_INT32 readOrbitIndex, orbitNumbers[4];
	stream >> size1 >> readOrbitIndex >> orbitNumbers[0] >> orbitNumbers[1] >> orbitNumbers[2] >> orbitNumbers[3] >> size2;
	ASSERT(readOrbitIndex == orbitIndex);
	//cout << "orbitIndex: " <<  orbitIndex << "(read: " << readOrbitIndex << ")" << endl;
	
	orbitidentifier orbitId(orbitNumbers[0], orbitNumbers[1], orbitNumbers[2]);
	orbitToIndex[orbitId] = orbitIndex;
	

	
	//cout << "orbitNumbers: " <<  orbitNumbers[0] << ", " << orbitNumbers[1]  << ", " << orbitNumbers[2]  << ", " << orbitNumbers[3] << endl;
	ASSERT(size1 == size2);
	ASSERT(orbitNumbers[3] == 0);
	Orbit* orbit = new Orbit(orbitIndex, orbitNumbers[0], orbitNumbers[1], orbitNumbers[2], noConstraints, noMaxVelocityHistograms, noMassConstraints);
	orbits.push_back(orbit);
	
	stream >> size1;
	//cout << "size1(a) = " << size1 << endl;
	for_each(orbittypes, orbittypes+noOrbitsPerIndex, sr); //(*s2) >> _1);
	//for_each(orbittypes, orbittypes+noOrbitsPerIndex, cout << constant("orbittype = ") << _1 << '\n');
	stream >> size2;
	//cout << "size2(a) = " << size2 << endl;
	ASSERT(size1 == size2);
		
		//noPhi, noTheta, noRadii
		//double* quad_light = new double[smom1,noPhi];
		//matrix<double> m(smom1,noPhi,noTheta);
		
	stream >> size1;
	//cout << "size1(b) = " << size1 << endl;
	for_each(orbit->massconstraints, orbit->massconstraints+noMassConstraints, sr); //(*s2) >> _1);
	orbit->light = new double[noRadii];
	for(int radiusIndex = 0; radiusIndex < noRadii; radiusIndex++) {
		orbit->light[radiusIndex] = 0;
		for(int thetaIndex = 0; thetaIndex < noTheta; thetaIndex++) {
			for(int phiIndex = 0; phiIndex < noPhi; phiIndex++) {
				orbit->light[radiusIndex] +=  orbit->massconstraints[ ((radiusIndex*noTheta+ thetaIndex)*noPhi + phiIndex)*16 + 0]; 
			}
		}
	}
		//for_each(massconstraints, massconstraints+noMassConstraints, cout << constant("mass = ") << _1 << '\n');
	stream >> size2;
	
	//cout << "size2(b) = " << size2 << endl;
	ASSERT(size1 == size2);
		
	maxHistogramValue = 0;
	for(int constraint = 0; constraint < noConstraints; constraint++) {
		orbit->histograms.push_back(new VelocityHistogram(noMaxVelocityHistograms));
		Q_INT32 ivmin, ivmax;
		stream >> size1 >> ivmin >> ivmax >> size2;
		//cout << "size(1,1) " << size1 << "/" << size2;
		//cout << "ivmin/ivmax: " <<  ivmin << "/" << ivmax << endl;
			//cout << "noMaxVelocityHistograms, dvhist: " <<  noMaxVelocityHistograms << ", " << dvhist << endl;
		ASSERT(size1 == size2);
		
		if(ivmin <= ivmax) {
			int noVelocityHistograms = (ivmax-ivmin+1);
			stream >> size1;
			//cout << noVelocityHistograms << "/";
			//cout << "size1 = " << size1 << endl;
			for_each(velocityHistogram, velocityHistogram+noVelocityHistograms, sr);
			//for_each(velocityHistogram, velocityHistogram+noVelocityHistograms, cout << constant("value = ") << _1 << '\n'); 
			stream >> size2;
			//cout << "size2 = " << size2 << endl;
			ASSERT(size1 == size2);
			for(int i = 0; i < noVelocityHistograms; i++) {
				int index = noMaxVelocityHistograms/2+ivmin+i;
				ASSERT(index >= 0 && index < noMaxVelocityHistograms);
				if(!(index >= 0 && index < noMaxVelocityHistograms)) {
					size2 = 1;
				}
				//cout << index << " ";
				orbit->histograms[constraint]->amplitudes[index] = velocityHistogram[i];
				maxHistogramValue = velocityHistogram[i] > maxHistogramValue ? velocityHistogram[i] : maxHistogramValue;
				orbit->maxHistogramValue = maxHistogramValue;
			}
		}
	}
		
}
}