#include "df.hpp"
#include "integration.hpp"
#include <boost/python.hpp>

namespace gd {
using namespace boost::python;

void py_export_df() {
	class_<DF >("DF", init<ProfileModel*, double_vector, double_vector>())
		.def("moments", &DF::moments)
		.def("momentsL", &DF::momentsL)
		.def("momentsL2", &DF::momentsL2)
		.def("momentsE", &DF::momentsE)
		//.def("moments_projected", &DF::moments_projected)
	;

}
/*
double DF::moments_projected(double_vector mass3d, double_vector density, double_vector v, double E, double L, double rp, double ra, double N) {
	int nR = Rs.size();
	int nr = rs.size();
	double* Rsp = Rs.data().begin();
	double* rbordersp = rborders.data().begin();
	double* rsp = rs.data().begin();
	double* densityp = density.data().begin();
	double* mass3dp = mass3d.data().begin();*/
	//double* vp = v.data().begin();
	/*
	auto density_integrand_norm = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	*/
	//IntegratorGSL<> density_integrator_norm(density_integrand_norm);
/*double R = 0;
	int rindex = 0;
	double rrho = 0;
	auto density_integrand = [&](double r) -> double {
		double rhor = mass3dp[rindex]/(4*M_PI*rrho*rrho*rrho);//;;1/N/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
		//double rhor = 1/N/sqrt(r*r*2*(E-profile_model->potentialr(r))-L*L);
		//return rhor /r/sqrt(r*r-R*R)*R;
		//printf("r = %10f : %10f : %10f : %10f\n", r, rhor, sqrt(r*r-R*R)/r, rhor*sqrt(r*r-R*R)*r);
		//return rhor /sqrt(r*r-R*R)/r*R;
		//return rhor /sqrt(r*r/(R*R)-1)*r;
		return rhor/sqrt(r*r-R*R)*r*R;
		//return rhor/sqrt(r*r-R*R)/r*R;
	};
	IntegratorGSL<> density_integrator(density_integrand, 1000000);
	//double N = density_integrator.integrate(rp, ra);
	
	/*auto varvr_integrand = [&](double r) -> double {
		return 1*sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> varvr_integrator(varvr_integrand);
	*/
	
/*for(int i = 0; i < nR; i++) {
		//double r1 = rbordersp[i];
		//double r2 = rbordersp[i+1];
		R = Rsp[i]; //(r2+r1)/2;
		//double Er1 = 2*(E-profile_model->potentialr(r1))-L*L/(r1*r1);
		//double Er2 = 2*(E-profile_model->potentialr(r2))-L*L/(r2*r2);
		if(R < ra) { // orbits only contribute to R < ra
			/*double R1 = R;//min(R, rp); // integral is 0 for r<rp
			if(rp > R1)
				R1 = rp;//*1.0001;
			double R2 = ra;//*0.999;
			*/
/*for(int j = 0; j < nr; j++) {
				double r1 = rbordersp[j];
				double r2 = rbordersp[j+1];
				rrho = rsp[j]; //(r1+r2)/2;
				rindex = j;
				if((mass3dp[rindex]>0) & (r2 > R)) { // if r2 < R, no mass
					if(r1<R)
						r1 = R; 
					densityp[i] += density_integrator.integrate(r1,r2);//(R1+R2)/2, R2);
				}
			}
			/**/
			/*
			printf("R=%f r1=%f r2=%f R1=%f R2=%f rp=%f ra=%f\n", R, r1, r2, R1, R2, rp, ra);
			double s = 1.001;
			printf("%f %f %f %f\n", density_integrand(R1), density_integrand(R1*s), density_integrand(R2/s), density_integrand(R2));
			printf("%f %f %f %f\n", density_integrand((R1+R2)/2), density_integrand(R1*s), density_integrand(R2/s), density_integrand(R2));
			/**/
			
			/*
			 * try{
				//densityp[i] = density_integrator.integrate_to_sigularity(R1, R2);//(R1+R2)/2, R2);
				densityp[i] = density_integrator.integrate_to_sigularity((R1+R2)/2, R2);//(R1+R2)/2, R2);
				densityp[i] += -density_integrator.integrate((R1+R2)/2, R1);
			} catch(std::range_error& e) {
			}
			/**/
			//varvrp[i] = varvr_integrator.integrate(r1, r2)/N;
//} else {
			/*if((r1<=rp) & (r2>=ra)) { // peri and api are inside the bin
				r1 = rp;
				r2 = ra;
				densityp[i] = 1; //density_integrator.integrate(r1, r2)/N;
				varvrp[i] = varvr_integrator.integrate(r1, r2)/N;
			}*/
//	}
//	}
//	double total = 0;
	/*for(int i = 0; i < n; i++) {
		total += densityp[i];
	}
	for(int i = 0; i < n; i++) {
		if(total > 0)
			densityp[i] /= total;
	}*/
//return 0;
//}

double DF::velocities(double_vector density, double_vector vr, double_vector vt, double E, double L, double rp, double ra) {
	int n = density.size();
	double* rbordersp = rborders.data().begin();
	double* densityp = density.data().begin();
	double* vrp = vr.data().begin();
	double* vtp = vt.data().begin();
	/*
	auto density_integrand_norm = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	*/
	//IntegratorGSL<> density_integrator_norm(density_integrand_norm);

	auto density_integrand = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> density_integrator(density_integrand);
	double N = density_integrator.integrate(rp, ra);
	
	auto vr_integrand = [&](double r) -> double {
		return sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> vr_integrator(vr_integrand);
	
	
	for(int i = 0; i < n; i++) {
		double r1 = rbordersp[i];
		double r2 = rbordersp[i+1];
		double Er1 = 2*(E-profile_model->potentialr(r1))-L*L/(r1*r1);
		double Er2 = 2*(E-profile_model->potentialr(r2))-L*L/(r2*r2);
		if((Er1>0) | (Er2>0)) {
			// only integrator over the allowable region
			r1 = Er1 > 0 ? r1 : rp;
			r2 = Er2 > 0 ? r2 : ra;
			densityp[i] = density_integrator.integrate(r1, r2)/N;
			vrp[i] = vr_integrator.integrate(r1, r2);
			vtp[i] = (L*log(r2)-L*log(r1)); // integrated analytically
		} else {
			if((r1<=rp) & (r2>=ra)) { // peri and api are inside the bin
				r1 = rp;
				r2 = ra;
				densityp[i] = 1; //density_integrator.integrate(r1, r2)/N;
				vrp[i] = vr_integrator.integrate(r1, r2);
				vtp[i] = (L*log(r2)-L*log(r1)); // integrated analytically
			}
		}
	}
	return N;
}

double DF::moments(double_vector density, double_vector varvr, double_vector varvt, double E, double L, double rp, double ra) {
	int n = density.size();
	double* rbordersp = rborders.data().begin();
	double* densityp = density.data().begin();
	double* varvrp = varvr.data().begin();
	double* varvtp = varvt.data().begin();
	/*
	auto density_integrand_norm = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	*/
	//IntegratorGSL<> density_integrator_norm(density_integrand_norm);

	auto density_integrand = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> density_integrator(density_integrand);
	double N = density_integrator.integrate(rp, ra);
	
	auto varvr_integrand = [&](double r) -> double {
		return 1*sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> varvr_integrator(varvr_integrand);
	
	auto varvt_integrand = [&](double r) -> double {
		return L*L/(r*r)/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r))/2;
	};
	IntegratorGSL<> varvt_integrator(varvt_integrand);
	int key = GSL_INTEG_GAUSS41;
	for(int i = 0; i < n; i++) {
		double r1 = rbordersp[i];
		double r2 = rbordersp[i+1];
		double Er1 = 2*(E-profile_model->potentialr(r1))-L*L/(r1*r1);
		double Er2 = 2*(E-profile_model->potentialr(r2))-L*L/(r2*r2);
		if((Er1>0) | (Er2>0)) {
			// only integrator over the allowable region
			if( (Er1 > 0) & (Er2 > 0) ) {
				densityp[i] = density_integrator.integrate_no_singularities(r1, r2, key)/N;
				varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				varvtp[i] = varvt_integrator.integrate_no_singularities(r1, r2, key)/N;
			} else {
				r1 = Er1 > 0 ? r1 : rp;
				r2 = Er2 > 0 ? r2 : ra;
				densityp[i] = density_integrator.integrate(r1, r2)/N;
				varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				varvtp[i] = varvt_integrator.integrate(r1, r2)/N;
			}
		} else {
			if((r1<=rp) & (r2>=ra)) { // peri and api are inside the bin
				r1 = rp;
				r2 = ra;
				densityp[i] = 1; //density_integrator.integrate(r1, r2)/N;
				varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				varvtp[i] = varvt_integrator.integrate(r1, r2)/N;
			}
		}
	}
	return N;
	
}

double DF::momentsL(double_vector density, double_vector varvr, double_vector varvt, double E, double L1, double L2, double rp1, double ra1, double rp2, double ra2, double L, double rp, double ra) {
	int n = density.size();
	double* rbordersp = rborders.data().begin();
	double* densityp = density.data().begin();
	double* varvrp = varvr.data().begin();
	double* varvtp = varvt.data().begin();
	/*
	auto density_integrand_norm = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	*/
	//IntegratorGSL<> density_integrator_norm(density_integrand_norm);

	auto density_integrand_L = [&](double r) -> double {
		// compute the max L that is possible
		double Lclip = sqrt(2*(E-profile_model->potentialr(r)))*r * 0.9999;
		// clip the L's to this value
		double L1p = min(L1, Lclip);
		double L2p = min(L2, Lclip);
		double Er1 = 2*(E-profile_model->potentialr(r))-L1p*L1p/(r*r);
		double Er2 = 2*(E-profile_model->potentialr(r))-L2p*L2p/(r*r);
		double v =  r* ( atan(L2p/sqrt(Er2)/r) - atan(L1p/sqrt(Er1)/r) );
		//printf("%f %f %f %f : L %f %f %f %f %f: E %f %f\n", Er1, Er2, r, v, L1, L2, L1p, L2p, Lclip, Er1, Er2);
		return v;
	};
	auto density_integrand = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> density_integrator(density_integrand_L);
	//printf("N rs: %f %f %f %f\n", rp1, rp2, ra1, ra2); 
	double N = 1;//density_integrator.integrate(min(rp1, rp2), max(ra1,ra2));
	
	auto varvr_integrand = [&](double r) -> double {
		return 1*sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> varvr_integrator(varvr_integrand);
	
	auto varvt_integrand = [&](double r) -> double {
		return L*L/(r*r)/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r))/2;
	};
	IntegratorGSL<> varvt_integrator(varvt_integrand);
	int key = GSL_INTEG_GAUSS41;
	for(int i = 0; i < n; i++) {
		double r1 = rbordersp[i];
		double r2 = rbordersp[i+1];
		//double Er1 = 2*(E-profile_model->potentialr(r1))-L*L/(r1*r1);
		//double Er2 = 2*(E-profile_model->potentialr(r2))-L*L/(r2*r2);
		varvrp[i] = i;
		varvtp[i] = i;
		if( (r1 > min(rp1, rp2)) | (r2 < max(ra1, ra2)) ) {// (Er1>0) | (Er2>0)) {
			// only integrator over the allowable region
			//double r1d = max(min(rp1, rp2), r1);
			//double r2d = min(max(ra2, ra2), r2);
			//printf("rs: %f %f\n", r1d, r2d);
			//if((r1d < r2d)) { // (Er1 > 0) & (Er2 > 0) ) {
				densityp[i] = density_integrator.integrate(r1, r2)/N;
				//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				//varvtp[i] = varvt_integrator.integrate_no_singularities(r1, r2, key)/N;
			//}
				/* else {
				r1 = Er1 > 0 ? r1 : rp;
				r2 = Er2 > 0 ? r2 : ra;
				densityp[i] = density_integrator.integrate(r1d, r2d)/N;
				//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				//varvtp[i] = varvt_integrator.integrate(r1, r2)/N;
			}*/
		} else {
			if((r1<=rp) & (r2>=ra)) { // peri and api are inside the bin
				r1 = rp;
				r2 = ra;
				//densityp[i] = 1; //density_integrator.integrate(r1, r2)/N;
				//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				//varvtp[i] = varvt_integrator.integrate(r1, r2)/N;
			}
		}
	}
	return N;
	
}	

double DF::momentsE(double_vector density, double_vector varvr, double_vector varvt, double E1, double E2, double L, double rp1, double ra1, double rp2, double ra2, double E, double rp, double ra) {
	int n = density.size();
	double* rbordersp = rborders.data().begin();
	double* densityp = density.data().begin();
	double* varvrp = varvr.data().begin();
	double* varvtp = varvt.data().begin();
	/*
	auto density_integrand_norm = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	*/
	//IntegratorGSL<> density_integrator_norm(density_integrand_norm);

	auto density_integrand_E = [&](double r) -> double {
		// compute the max E that is possible
		//double Eclip = profile_model->potentialr(r) + L*L/(2*r*r) * 1.001; //9999;
		// clip the E's to this value
		//double E1p = max(E1, Eclip);
		//double E2p = max(E2, Eclip);
		double Er1 = 2*(E1-profile_model->potentialr(r))-L*L/(r*r);
		double Er2 = 2*(E2-profile_model->potentialr(r))-L*L/(r*r);
		double v = 0;
		if(Er1 > 0)
			v += sqrt(Er1);
		if(Er2 > 0)
			v -= sqrt(Er2);
		v *= 2;
		//double v =  r* ( atan(L2p/sqrt(Er2)/r) - atan(L1p/sqrt(Er1)/r) );
		//double v = 2 * (sqrt(Er1) - sqrt(Er2));
		//printf("%f %f %f %f : L %f %f\n", Er1, Er2, r, v, E1, E2);
		return -v;
	};
	auto density_integrand = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> density_integrator(density_integrand_E);
	//printf("N rs: %f %f %f %f\n", rp1, rp2, ra1, ra2); 
	double N = 1;//density_integrator.integrate(min(rp1, rp2), max(ra1,ra2));
	
	auto varvr_integrand = [&](double r) -> double {
		return 1*sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	IntegratorGSL<> varvr_integrator(varvr_integrand);
	
	auto varvt_integrand = [&](double r) -> double {
		return L*L/(r*r)/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r))/2;
	};
	IntegratorGSL<> varvt_integrator(varvt_integrand);
	int key = GSL_INTEG_GAUSS41;
	for(int i = 0; i < n; i++) {
		double r1 = rbordersp[i];
		double r2 = rbordersp[i+1];
		//double Er1 = 2*(E-profile_model->potentialr(r1))-L*L/(r1*r1);
		//double Er2 = 2*(E-profile_model->potentialr(r2))-L*L/(r2*r2);
		varvrp[i] = i;
		varvtp[i] = i;
		if( (r1 > min(rp1, rp2)) & (r2 < max(ra1, ra2)) ) {// (Er1>0) | (Er2>0)) {
			// only integrator over the allowable region
			//double r1d = max(min(rp1, rp2), r1);
			//double r2d = min(max(ra2, ra2), r2);
			//printf("rs: %f %f : %f %f : %f %f \n", r1, r2, rp1, ra1, rp2, ra2);
			//if((r1d < r2d)) { // (Er1 > 0) & (Er2 > 0) ) {
			// TODO: clip values r1 and r2 such that there is al least one positive energyy
				{	densityp[i] = density_integrator.integrate(r1, r2)/N;
				//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				//varvtp[i] = varvt_integrator.integrate_no_singularities(r1, r2, key)/N;
			}/* else {
				r1 = Er1 > 0 ? r1 : rp;
				r2 = Er2 > 0 ? r2 : ra;
				densityp[i] = density_integrator.integrate(r1d, r2d)/N;
				//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				//varvtp[i] = varvt_integrator.integrate(r1, r2)/N;
			}*/
		} else {
			if((r1<=rp) & (r2>=ra)) { // peri and api are inside the bin
				r1 = rp;
				r2 = ra;
				//densityp[i] = 1; //density_integrator.integrate(r1, r2)/N;
				//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				//varvtp[i] = varvt_integrator.integrate(r1, r2)/N;
			}
		}
	}
	return N;
	
}

double DF::momentsL2(double_vector density, double_vector varvr, double_vector varvt, double_vector v22, double_vector v40, double_vector v04, double E, double L1, double L2, double rp, double ra, bool debug) {
	int n = density.size();
	double* rbordersp = rborders.data().begin();
	double* densityp = density.data().begin();
	double* varvrp = varvr.data().begin();
	double* varvtp = varvt.data().begin();
	double* v40p = v40.data().begin();
	double* v04p = v04.data().begin();
	double* v22p = v22.data().begin();
/*
	auto density_integrand_norm = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};
	*/
	//IntegratorGSL<> density_integrator_norm(density_integrand_norm);
	double N = 1;
	double epsabs = 0;//1e-6;
	double epsrel = 1e-6;
	auto density_integrand_L = [&](double r) -> double {
		// compute the max L that is possible
		double Lmax = sqrt(2*(E-profile_model->potentialr(r)))*r;// * 0.9999;
		// clip the L's to this value
		double v1 = 0;
		double v2 = 0;
		if(L1 >= Lmax)
			return 0;
		if(L1 < Lmax) {
			double Er = 2*(E-profile_model->potentialr(r))-L1*L1/(r*r);
			v1 = atan(L1/(sqrt(Er)*r));
		} else {
			v1 = M_PI/2; //atan(1);
		}
		if(L2 < Lmax) {
			double Er = 2*(E-profile_model->potentialr(r))-L2*L2/(r*r);
			v2 = atan(L2/(sqrt(Er)*r));
		} else {
			v2 = M_PI/2;//atan(1);
		}
			
		//double L1p = min(L1, Lclip);
		//double L2p = min(L2, Lclip);
		//double Er1 = 2*(E-profile_model->potentialr(r))-L1p*L1p/(r*r);
		//double Er2 = 2*(E-profile_model->potentialr(r))-L2p*L2p/(r*r);
		//double v =  r* ( atan(L2p/(sqrt(Er2)*r)) - atan(L1p/(sqrt(Er1)*r)) );
		double v =  r* ( v2 - v1 );
		if((v < 0) & (fabs((v2 - v1)/(v2+v1)) < 1e-8))
			v = 0;
		if(v < 0) {
			printf("L %f %f %f %f %f]\n", L1, L2, Lmax, v1, v2);
			assert(v >= 0);
		}
		if(debug)
			printf("L1=%f L2=%f Lmax=%f E=%f v1=%f v2=%f dv=%f r=%f\n", L1, L2, Lmax, E, v1, v2, v2-v1, r);
		
		//printf("[%f %f %f %f : L %f %f %f %f %f: E %f %f]\n", Er1, Er2, r, v, L1, L2, L1p, L2p, Lclip, Er1, Er2);
		return v;
	};
	/*auto density_integrand = [&](double r) -> double {
		return 1/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r));
	};*/
	IntegratorGSL<> density_integrator(density_integrand_L, 10000, epsrel, epsabs);
	//printf("N rs: %f %f %f %f\n", rp1, rp2, ra1, ra2); 
	//double N = density_integrator.integrate(min(rp1, rp2), max(ra1,ra2));
	try {
		N = density_integrator.integrate_to_sigularity(rp, ra);
	} catch(std::exception& e) {
		printf("could not compute normalization rp=%f ra=%f\n", rp, ra);
		throw;
	}
	
	if(debug)
		printf("N = %f\n", N);
	
	auto varvr_integrand = [&](double r) -> double {
		double Lmax = sqrt(2*(E-profile_model->potentialr(r)))*r;
		double v1 = 0;
		double v2 = 0;
		double c = 2*(E-profile_model->potentialr(r)) * r * r;
		if(L1 >= Lmax)
			return 0;
		if(L1 < Lmax) {
			//double Er = 2*(E-profile_model->potentialr(r))-L1*L1/(r*r);
			v1 = (L1 * sqrt(c-L1*L1) + c * atan(L1/sqrt(c-L1*L1)))/r/2;
		} else {
			v1 = M_PI/4 * c/r;
		}
		if(L2 < Lmax) {
			v2 = (L2 * sqrt(c-L2*L2) + c * atan(L2/sqrt(c-L2*L2)))/r/2;
		} else {
			v2 = M_PI/4 * c/r;
		}
		
		double v =  ( v2 - v1 );
		if((v < 0) & (fabs((v2 - v1)/(v2+v1)) < 1e-8))
			v = 0;
		if(v < 0) {
			printf("L %f %f %f %f %f]\n", L1, L2, Lmax, v1, v2);
			assert(v >= 0);
		}
		//printf("[%f %f %f %f : L %f %f %f %f %f: E %f %f]\n", Er1, Er2, r, v, L1, L2, L1p, L2p, Lclip, Er1, Er2);
		return v;
	};
	IntegratorGSL<> varvr_integrator(varvr_integrand, 10000, epsrel, epsabs);
	
	
	auto v40_integrand = [&](double r) -> double {
		double Lmax = sqrt(2*(E-profile_model->potentialr(r)))*r;
		double v1 = 0;
		double v2 = 0;
		double c = 2*(E-profile_model->potentialr(r)) * r * r;
		assert(L1 < L2);
		if(L1 >= Lmax)
			return 0;
		if(L1 < Lmax) {
			//double Er = 2*(E-profile_model->potentialr(r))-L1*L1/(r*r);
			v1 = (L1 * (5*c - 2*L1*L1) * sqrt(c-L1*L1) + 3 * c*c * atan(L1/sqrt(c-L1*L1)))/(8) /r/ r/r;
		} else {
			v1 = M_PI/16 * 3 * c * c /r/r/r;
		}
		if(L2 < Lmax) {
			v2 = (L2 * (5*c - 2*L2*L2) * sqrt(c-L2*L2) + 3 * c*c * atan(L2/sqrt(c-L2*L2)))/(8) /r/r/r;
			//v2 = (L2 * sqrt(c-L2*L2) + c * atan(L2/sqrt(c-L2*L2)))/r/2;
		} else {
			v2 = M_PI/16 * 3 * c * c/r/r/r;
		}
		
		double v =  ( v2 - v1 );
		if((v < 0) & (fabs((v2 - v1)/(v2+v1)) < 1e-8))
		   v = 0;
		if(v < 0) {
			printf("L %e %e %e %e %e %e]\n", L1, L2, Lmax, v1, v2, v);
			assert(v >= 0);
		}
		//printf("[%f %f %f %f : L %f %f %f %f %f: E %f %f]\n", Er1, Er2, r, v, L1, L2, L1p, L2p, Lclip, Er1, Er2);
		return v;
	};
	IntegratorGSL<> v40_integrator(v40_integrand, 10000, epsrel, epsabs);	
	
	
	auto v04_integrand = [&](double r) -> double {
		double Lmax = sqrt(2*(E-profile_model->potentialr(r)))*r;
		double v1 = 0;
		double v2 = 0;
		double c = 2*(E-profile_model->potentialr(r)) * r * r;
		if(L1 >= Lmax)
			return 0;
		if(L1 < Lmax) {
			//double Er = 2*(E-profile_model->potentialr(r))-L1*L1/(r*r);
			v1 = (-L1 * (3*c + 2*L1*L1) * sqrt(c-L1*L1) + 3 * c*c * atan(L1/sqrt(c-L1*L1)))/(8) /r/r/r;
		} else {
			v1 = M_PI/16 * 3 * c * c /r/r/r;
		}
		if(L2 < Lmax) {
			v2 = (-L2 * (3*c + 2*L2*L2) * sqrt(c-L2*L2) + 3 * c*c * atan(L2/sqrt(c-L2*L2)))/(8) /r/r/r;
			//v2 = (L2 * sqrt(c-L2*L2) + c * atan(L2/sqrt(c-L2*L2)))/r/2;
		} else {
			v2 = M_PI/16 * 3 * c * c/r/r/r;
		}
		
		double v =  ( v2 - v1 );
		if((v < 0) & (fabs((v2 - v1)/(v2+v1)) < 1e-8))
			v = 0;
		if(v < 0) {
			printf("L %f %f %f %f %f]\n", L1, L2, Lmax, v1, v2);
			assert(v >= 0);
		}
		if(debug)
			printf("L1=%f L2=%f Lmax=%f E=%f v1=%f v2=%f dv=%f c=%f r=%f\n", L1, L2, Lmax, E, v1, v2, v2-v1, c, r);
		//printf("[%f %f %f %f : L %f %f %f %f %f: E %f %f]\n", Er1, Er2, r, v, L1, L2, L1p, L2p, Lclip, Er1, Er2);
		return v;
	};
	IntegratorGSL<> v04_integrator(v04_integrand, 10000, epsrel, epsabs);
	
	auto v22_integrand = [&](double r) -> double {
		double Lmax = sqrt(2*(E-profile_model->potentialr(r)))*r;
		double v1 = 0;
		double v2 = 0;
		double c = 2*(E-profile_model->potentialr(r)) * r * r;
		const int p = 3;
		if(L1 >= Lmax)
			return 0;
		if(L1 < Lmax) {
			//double Er = 2*(E-profile_model->potentialr(r))-L1*L1/(r*r);
			v1 = (L1 * (-c +2*L1*L1) * sqrt(c-L1*L1) + c*c * atan(L1/sqrt(c-L1*L1)))/(8) /pow(r, p);
		} else {
			v1 = M_PI/16 * c*c /pow(r, p);
		}
		if(L2 < Lmax) {
			v2 = (L2 * (-c +2*L2*L2) * sqrt(c-L2*L2) + c*c * atan(L2/sqrt(c-L2*L2)))/(8) /pow(r, p);
			//v2 = (L2 * sqrt(c-L2*L2) + c * atan(L2/sqrt(c-L2*L2)))/r/2;
		} else {
			v2 = M_PI/16 * c*c /pow(r, p);
		}
		
		double v =  ( v2 - v1 );
		if((v < 0) & (fabs((v2 - v1)/(v2+v1)) < 1e-8))
			v = 0;
		if(v < 0) {
			printf("L %f %f %f %f %f]\n", L1, L2, Lmax, v1, v2);
			assert(v >= 0);
		}
		//printf("[%f %f %f %f : L %f %f %f %f %f: E %f %f]\n", Er1, Er2, r, v, L1, L2, L1p, L2p, Lclip, Er1, Er2);
		return v;
	};
	IntegratorGSL<> v22_integrator(v22_integrand, 10000, epsrel, epsabs);
	
	/*auto varvt_integrand = [&](double r) -> double {
		return L*L/(r*r)/sqrt(2*(E-profile_model->potentialr(r))-L*L/(r*r))/2;
	};*/
	
	auto varvt_integrand = [&](double r) -> double {
		double Lmax = sqrt(2*(E-profile_model->potentialr(r)))*r;
		double v1 = 0;
		double v2 = 0;
		double c = 2*(E-profile_model->potentialr(r)) * r * r;
		if(L1 < Lmax) {
			v1 = (-L1 * sqrt(c-L1*L1) + c * atan(L1/sqrt(c-L1*L1)))/r/2;
		} else {
			v1 = M_PI/4 * c/r;
		}
		if(L2 < Lmax) {
			v2 = (-L2 * sqrt(c-L2*L2) + c * atan(L2/sqrt(c-L2*L2)))/r/2;
		} else {
			v2 = M_PI/4 * c/r;
		}
		
		double v =  ( v2 - v1 );
		if((v < 0) & (fabs((v2 - v1)/(v2+v1)) < 1e-8))
			v = 0;
		if(v < 0) {
			printf("L %f %f %f %f %f]\n", L1, L2, Lmax, v1, v2);
			assert(v >= 0);
		}
		//printf("[%f %f %f %f : L %f %f %f %f %f: E %f %f]\n", Er1, Er2, r, v, L1, L2, L1p, L2p, Lclip, Er1, Er2);
		return v;
	};
	IntegratorGSL<> varvt_integrator(varvt_integrand, 10000, epsrel, epsabs);
	
	
	int key = GSL_INTEG_GAUSS41;
	//double r_apo = rbordersp[n]; / r borders array should have length n+1 
	//double r_apo = rbordersp[n];
	for(int i = 0; i < n; i++) {
		double r1 = rbordersp[i];
		double r2 = rbordersp[i+1];
		double Er1 = 2*(E-profile_model->potentialr(r1))-L1*L1/(r1*r1);
		double Er2 = 2*(E-profile_model->potentialr(r2))-L1*L1/(r2*r2);
		varvrp[i] = 0;
		varvtp[i] = 0;
		if(debug)
			printf("r1=%f r2=%f\n", r1, r2);
		if((Er1 > 0) | (Er2 > 0)) { //  (r1 > min(rp1, rp2)) | (r2 < max(ra1, ra2)) ) {// (Er1>0) | (Er2>0)) {
			if(Er1 < 0) {
				if(debug)
					printf("looking for root r1 [%f, %f]\n", r1, r2);
				auto ekinrad = [&](double r) -> double { // kinetic energy in radial direction
					return (-L1*L1/(2*r*r) - (this->profile_model->potentialr(r) - E));
				};
				RootFinderGSL<> rootfinder(ekinrad);
				r1 = rootfinder.findRoot(r1, r2);
				double scale = (1+1e-5);
				while(ekinrad(r1) < 0) {
					r1 *= scale;
				}
				if(debug)
					printf("found root r1 [%f, %f]\n", r1, r2);
			}
			if(Er2 < 0) {
				if(debug)
					printf("looking for root r2 [%f, %f]\n", r1, r2);
				auto ekinrad = [&](double r) -> double { // kinetic energy in radial direction
					return (-L1*L1/(2*r*r) - (this->profile_model->potentialr(r) - E));
				};
				RootFinderGSL<> rootfinder(ekinrad);
				double scale = (1+1e-5);
				r2 = rootfinder.findRoot(r1, r2);
				while(ekinrad(r2) < 0) {
					r2 /= scale;
				}
				if(debug)
					printf("found root r2 [%f, %f]\n", r1, r2);
			}
			if(r1<r2) {
				//printf("v00\n"); fflush(stdout);
				double rc = (r1+r2)/2;
				try {
					//double rho = density_integrator.integrate(r1, r2)/N;
					densityp[i] = -density_integrator.integrate_to_sigularity(rc,r1)/N*2;
					densityp[i] += density_integrator.integrate_to_sigularity(rc,r2)/N*2;
					//densityp[i] = rho;
				} catch(std::exception& e) {
					printf("could not compute density for r[%f,%f]\n", r1, r2);
					throw;
				}
				//printf("rs: %f %f, rho=%f\n", r1, r2, rho);
				//printf("v20\n"); fflush(stdout);
				try {
					varvrp[i] = -varvr_integrator.integrate_to_sigularity(rc,r1)/N*2;
					varvrp[i] += varvr_integrator.integrate_to_sigularity(rc,r2)/N*2;
					//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				} catch(std::exception& e) {
					printf("could not compute v20 for r[%f,%f]\n", r1, r2);
					throw;
				}
				//printf("v02\n"); fflush(stdout);
				try {
					varvtp[i] = -varvt_integrator.integrate_to_sigularity(rc,r1)/N*2;
					varvtp[i] += varvt_integrator.integrate_to_sigularity(rc,r2)/N*2;
					//varvtp[i] = varvt_integrator.integrate_no_singularities(r1, r2, key)/N;
				} catch(std::exception& e) {
					printf("could not compute v02 for r[%f,%f]\n", r1, r2);
					throw;
				}
				//printf("v40\n"); fflush(stdout);
				try {
					v40p[i] = -v40_integrator.integrate_to_sigularity(rc,r1)/N*2;
					v40p[i] += v40_integrator.integrate_to_sigularity(rc,r2)/N*2;
					//v40p[i] = v40_integrator.integrate_no_singularities(r1, r2, key)/N;
				} catch(std::exception& e) {
					printf("could not compute v40 for r[%f,%f]\n", r1, r2);
					throw;
				}
				//printf("v04\n"); fflush(stdout);
				try {
					v04p[i] = -v04_integrator.integrate_to_sigularity(rc,r1)/N*2;
					v04p[i] += v04_integrator.integrate_to_sigularity(rc,r2)/N*2;
					//v04p[i] = v04_integrator.integrate_no_singularities(r1, r2, key)/N;
				} catch(std::exception& e) {
					printf("could not compute v04 for r[%f,%f] rp=%f ra=%f\n", r1, r2, rp, ra);
					throw;
				}
				//printf("v22\n"); fflush(stdout);
				try {
					v22p[i] = -v22_integrator.integrate_to_sigularity(rc,r1)/N*2;
					v22p[i] += v22_integrator.integrate_to_sigularity(rc,r2)/N*2;
					//v22p[i] = v22_integrator.integrate_no_singularities(r1, r2, key)/N;
				} catch(std::exception& e) {
					printf("could not compute v22 for r[%f,%f]\n", r1, r2);
					throw;
				}
			}
			//printf("done\n"); fflush(stdout);
//}
				/* else {
				r1 = Er1 > 0 ? r1 : rp;
				r2 = Er2 > 0 ? r2 : ra;
				densityp[i] = density_integrator.integrate(r1d, r2d)/N;
				//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				//varvtp[i] = varvt_integrator.integrate(r1, r2)/N;
			}*/
		} else {
			//if((r1<=rp) & (r2>=ra)) { // peri and api are inside the bin
			//	r1 = rp;
			//	r2 = ra;
				//densityp[i] = 1; //density_integrator.integrate(r1, r2)/N;
				//varvrp[i] = varvr_integrator.integrate_no_singularities(r1, r2, key)/N;
				//varvtp[i] = varvt_integrator.integrate(r1, r2)/N;
			//}
		}
	}
	return N;
	
}	

	
};

