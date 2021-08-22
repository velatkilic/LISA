
#include "ParticleDist.h"

//! @brief Calculates number of droplets [m^-3 mm^-1] for a given diameter and rain rate
//! @param d : Rain droplet diameter [mm] 
//! @param Rr : Rain rate [mm/hr]
//! @return Number of droplets per m^3 volume per mm diameter
double MarshallPalmerRain::N_model(double d, double Rr) {
	double out = 8000.0 * exp(-4.1 * pow(Rr, -0.21) * d);
	return out;
}

//! @brief Calculates particle density integrated from dst to inifinity for the Monte Carlo step
//! @param Rr : Rain rate [mm/hr]
//! @param dst : Particle diameter [mm] from where integration starts
//! @return Particle density (number of particles per m^3 of volume)
double MarshallPalmerRain::N_total(double Rr, double dst) {
	double lam = 4.1 * pow(Rr, -0.21);
	double out = 8000.0 * exp(-lam * dst) / lam;
	return out;
}

//! @brief Sample N particle diameters from Marshall-Palmer distribution whose diameter are above dst for the Monte Carlo step
//! @param Rr : Rain rate [mm/hr]
//! @param dst : Only particles with diameters above dst [mm] are considered
//! @param N : Number of particle diameters to sample
//! @return Particle diameters in units of millimeters (vector of size N) drawn from an exponential distribution (Marshall-Palmer model)
std::vector<double> MarshallPalmerRain::N_sample(double Rr, double dst, int N) {
	double lam = 4.1 * pow(Rr, -0.21);
	double r(0);
	std::vector<double> diameters(N);
	for (int i = 0; i < N; i++) {
		r = dis(gen); // generate a random number [0, 1]
		diameters[i] = (-log(1.0 - r) / lam) + dst;
	}
	return diameters;
}