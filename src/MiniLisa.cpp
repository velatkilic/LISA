#include "MiniLisa.h"

MiniLisa::MiniLisa(Lidar& lidar, Material& inclusion, ParticleDist& pdist) :
	pdist(pdist), lidar(lidar), inclusion(inclusion) {

	Pmin = 0.9 / pow(lidar.get_max_range(), 2.0);
	nref = inclusion.ref_ind(lidar.get_wavelength());
	d    = logspace(-6.0, 1.5, 2000);
	qext.resize(d.size());

	std::random_device rd;
	gen = std::mt19937(rd());
	_calc_qexts();
}

//! @brief Augment a point (vector of doubles) given extinction coefficient alpha. 
//! @param pc : Point (vector of doubles) of size 4: (x, y, z, reflectance)
//! @param alpha : Extinction coefficient in units of m^-1
//! @param Rr : Rain rate in units of mm/hr
//! @return Augmented point (vector of doubles) of size 5: (x_new, y_new, z_new, reflectance_new).
std::vector<double> MiniLisa::augment(std::vector<double>& pc, double alpha) {
	double& x = pc[0];
	double& y = pc[1];
	double& z = pc[2];
	double& ref = pc[3];

	double range = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0)); // calculate range
	std::vector<double> pcnew(4, 0.0);

	// If range smaller than lidar min range, return all zeros
	if (range < lidar.get_min_range()) {
		return pcnew;
	}

	// Calculate SNR
	double P0 = ref * exp(-2.0 * alpha * range) / pow(range, 2.0); // back reflected power in arbitrary units
	double SNR = P0 / Pmin; // signal to noise ratio

	// If back reflected power is below the noise floor, then return all zeros
	if (SNR < 1.0) {
		return pcnew;
	}

	// Calculate angles
	double phi = atan2(y, x);
	double the = acos(z / range);

	// Add SNR dependent range uncertainty and return "fuzzy" point cloud
	double sig = lidar.get_ran_uncer() / sqrt(2.0 * SNR); // range uncertainty
	std::normal_distribution<double> dist{ 0.0, sig }; // with zero mean and sig standard deviation
	double range_new = range + dist(gen); // add range noise

	// Calculate new point cloud based on new range
	pcnew[0] = range_new * sin(the) * cos(phi);
	pcnew[1] = range_new * sin(the) * sin(phi);
	pcnew[2] = range_new * cos(the);
	pcnew[3] = ref * exp(-2.0 * alpha * range);

	return pcnew;
}

//! @brief Augment a 2D point cloud (vector of vectors) given a rain rate.
//! @param pc : Point cloud
//! @param Rr : Rain rate in mm/hr
//! @return Augmented point cloud
std::vector<std::vector<double>> MiniLisa::augment(std::vector<std::vector<double>>& pc, double Rr) {
	// calculate extinction coefficient from the randomly sampled rainrate
	double alpha = calc_alpha(Rr);

	// Augment the lidar scan with the given extinction coefficient
	std::vector<std::vector<double>> pcnew(pc.size());
	for (size_t i = 0; i < pc.size(); i++) {
		pcnew[i] = augment(pc[i], alpha);
	}
	return pcnew;
}

//! @brief Augment a 2D point cloud (vector of vectors). Rain rate is selected from an exponential distribution.
//! @param pc : Point cloud
//! @return Augmented point cloud
std::vector<std::vector<double>> MiniLisa::augment(std::vector<std::vector<double>>& pc) {
	// Sample rain rate from an exponential distribution
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	double r = dis(gen);
	double Rr = -log(1. - r) / 0.05;

	// augment with the random rainrate
	std::vector<std::vector<double>> pcnew = augment(pc, Rr);

	return pcnew;
}

//! @brief Calculates extinction coefficients for given MiniLisa parameters. This methods might take some time but 
//! it is multi-threaded and runs only once when a MiniLisa object is initialized.
void MiniLisa::_calc_qexts() {
	std::cout << std::endl << "Calculating extinction efficiencies using Mie theory to init augmenter class ... " << std::endl;
	// calculate size parameter
	std::vector<double> x(d.size());
	for (size_t i = 0; i < d.size(); i++) {
		x[i] = 3.14159265359 * d[i] * 1.0e3 / lidar.get_wavelength(); // unit conversion needed for d[mm] to compare against wavelength[um]
	}

	// calculate extinction efficiency asynch
	std::vector<std::future<double>> future(d.size()); 
	for (size_t i = 0; i < d.size(); i++) {
		future[i] = std::async(calc_qext, x[i], nref);
	}

	// extract results
	for (size_t i = 0; i < d.size(); i++) {
		qext[i] = (future[i]).get();
	}
}

//! @brief Calculates extinction coefficient for a given rain rate
//! @param Rr : Rain rate in mm/hr
//! @return Extinction coefficient in m^-1
double MiniLisa::calc_alpha(double Rr) {
	std::vector<double> f(d.size());
	for (size_t i = 0; i < d.size(); i++) {
		f[i] = (pdist.N_model(d[i], Rr)) * pow(d[i], 2.0) * qext[i];
	}
	double alpha = 1.0e-6 * trapz(f, d) * 3.14159265359 / 4; // m^-1
	return alpha;
}
