#ifndef PARTICLEDIST_H
#define PARTICLEDIST_H

#include<vector>
#include<math.h>
#include<random>

class ParticleDist {
	public:
		ParticleDist() {
			std::random_device rd;
			gen = std::mt19937(rd());
			dis = std::uniform_real_distribution<double>(0.0, 1.0);
		}

		//! @brief Calculates number of droplets [m^-3 mm^-1] for a given diameter and "rain rate" (or any other intensity correlated number)
		//! @param d : Rain droplet diameter [mm] 
		//! @param Rr : Rain rate [mm/hr] or any other weather intensity correlated number
		//! @return Number of droplets per m^3 volume per mm diameter
		virtual double N_model(double d, double Rr) = 0;

		//! @brief Calculates particle density integrated from dst to inifinity for the Monte Carlo step
		//! @param Rr : Rain rate [mm/hr] or any other weather intensity correlated number
		//! @param dst : Particle diameter [mm] from where integration starts
		//! @return Particle density (number of particles per m^3 of volume)
		virtual double N_total(double Rr, double dst) = 0;

		//! @brief Sample N particle diameters from N_model whose diameter are above dst for the Monte Carlo step
		//! @param Rr : Rain rate [mm/hr] or any other weather intensity correlated number
		//! @param dst : Only particles with diameters above dst [mm] are considered
		//! @param N : Number of particle diameters to sample
		//! @return Particle diameters in units of millimeters (vector of size N) drawn from N_model
		virtual std::vector<double> N_sample(double Rr, double dst, int N) = 0; // 
	protected:

		//! @brief Random number generator
		std::mt19937 gen;

		//! @brief Uniform real distribution (initialized to [0, 1])
		std::uniform_real_distribution<double> dis;
};

class MarshallPalmerRain : public ParticleDist {
	public:
		double N_model(double d, double Rr) override;
		double N_total(double Rr, double dst) override;
		std::vector<double> N_sample(double Rr, double dst, int N) override;
};

#endif // PARTICLEDIST_H
