#ifndef MINILISA_H
#define MINILISA_H

#include "Augmenter.h"
#include "Lidar.h"
#include "Material.h"
#include "Utils.h"
#include "ParticleDist.h"
#include <random>
#include <future>

class MiniLisa : public Augmenter {
public:
	//! @brief MiniLisa default constructor
	MiniLisa():MiniLisa(Lidar(), Water(), MarshallPalmerRain()) {}


	//! @brief Construct Lisa object with given parameters
	//! @param lidar : Lidar object (i.e VLS_128 or user defined lidars)
	//! @param inclusion : Material object for the inclusion (i.e water for atmospheric phenomena such as snow, rain and fog)
	//! @param pdist : ParticleDist object for particle distibution (i.e Marshall-Palmer distribution for rain)
	MiniLisa(Lidar& lidar, Material& inclusion, ParticleDist& pdist);

	std::vector<std::vector<double>> augment(std::vector<std::vector<double>>& pc) override;
	virtual std::vector<std::vector<double>> augment(std::vector<std::vector<double>>& pc, double Rr);

	double calc_alpha(double Rr);
protected:
	//! @brief Minimum detectable power in arbitrary units. Inferred from maximum lidar range
	double Pmin;

	//! @brief Refractive index of the inclusion/scatterer (i.e water droplets)
	std::complex<double> nref;

	//! @brief Vector of particle diameters (in millimeters) for calculating average extinction coefficient
	std::vector<double> d;

	//! @brief Extinction efficiency as calculated from Mie theory at diameters d given inclusion refractive index and lidar laser wavelength
	std::vector<double> qext;

	//! @brief ParticleDist object (i.e Marshall-Palmer Rain)
	ParticleDist& pdist;

	//! @brief Lidar object (i.e VLS_128)
	Lidar& lidar;

	//! @brief Material object (i.e water)
	Material& inclusion;

	//! @brief Random number generator
	std::mt19937 gen;

	void _calc_qexts();

private:
	std::vector<double> augment(std::vector<double>& pc, double alpha);
};
#endif // MINILISA_H