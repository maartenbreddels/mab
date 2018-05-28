#pragma once
#include <boost/python.hpp>
#include "datatypes.hpp"
#include "vectorize.hpp"

namespace gd {

using namespace boost::python;

template<class ProfileType>
void fill(ProfileType *profile, double_vector r, double_vector target) {
	double* rp = r.data().begin();
	double* tp = target.data().begin();
	int size = r.size();
	for(int i = 0; i < size; i++) {
		tp[i] = profile->dphidr(rp[i]);
	}
	
}

template<class ProfileType, class Ctor>
	void py_export_profile_profile(string name) {
		class_<ProfileType, bases<Profile, Density> >(name.c_str(), Ctor())
			.def("potentialr", vectorize(&ProfileType::potentialr))
			.def("densityr", vectorize(&ProfileType::densityr))
			.def("densityR", vectorize(&ProfileType::densityR))
			.def("dphidr", vectorize(&ProfileType::dphidr))
			.def("fill", (&fill<ProfileType>))
		//.def("dphidr2", (&Plummer::dphidr2))
			;
	}

template<class ProfileType, class Ctor>
	void py_export_profile_profile_kw(string name, Ctor ctor) {
		class_<ProfileType, bases<Profile, Density> >(name.c_str(), ctor)
			.def("potentialr", vectorize(&ProfileType::potentialr))
			.def("densityr", vectorize(&ProfileType::densityr))
			.def("densityR", vectorize(&ProfileType::densityR))
			.def("dphidr", vectorize(&ProfileType::dphidr))
			.def("enclosed_mass", vectorize(&ProfileType::enclosed_mass))
			.def("fill", (&fill<ProfileType>))
		//.def("dphidr2", (&Plummer::dphidr2))
			;
	}

}