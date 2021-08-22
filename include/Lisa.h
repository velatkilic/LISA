#ifndef LISA_H
#define LISA_H

#include "MiniLisa.h"

class Lisa : public MiniLisa {
	public:

		//! @brief Lisa default constructor
		Lisa() : MiniLisa(), dst(0.05) { 
			ref_r = pow(std::abs((nref - 1.0) / (nref + 1.0)), 2.0); 
		}

		//! @brief Construct Lisa object with given parameters assuming minimum particle diameter of 50 um for the Monte Carlo part
		//! @param lidar : Lidar object (i.e VLS_128 or user defined lidars)
		//! @param inclusion : Material object for the inclusion (i.e water for atmospheric phenomena such as snow, rain and fog)
		//! @param pdist : ParticleDist object for particle distibution (i.e Marshall-Palmer distribution for rain)
		Lisa(Lidar& lidar, Material& inclusion, ParticleDist& pdist) : MiniLisa(lidar, inclusion, pdist), dst(0.05) { 
			ref_r = pow(std::abs((nref - 1.0) / (nref + 1.0)), 2.0); 
		}

		std::vector<std::vector<double>> augment(std::vector<std::vector<double>>& pc) override;
		std::vector<std::vector<double>> augment(std::vector<std::vector<double>>& pc, double Rr) override;

		
	private:
		//! @brief Only particles with diameters larger than dst (in millimeters, default is 0.05 mm)
		//! are tracked during Monte Carlo to keep computations simple.
		double dst;

		//! @brief Fresnel reflectivity at normal incidence.
		double ref_r;
		std::vector<double> augment(std::vector<double>& pc, double alpha, double Rr);
};

#endif // LISA_H
