#include "Lisa.h"

//! @brief Augment a 2D point cloud (vector of vectors). Rain rate is selected from an exponential distribution.
//! @param pc : Point cloud
//! @return Augmented point cloud
std::vector<std::vector<double>> Lisa::augment(std::vector<std::vector<double>>& pc)
{
	// Sample rain rate from an exponential distribution
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	double r = dis(gen);
	double Rr = -log(1. - r) / 0.05;

	// augment with the random rainrate
	std::vector<std::vector<double>> pcnew = augment(pc, Rr);

	return pcnew;
}

//! @brief Augment a 2D point cloud (vector of vectors) given a rain rate.
//! @param pc : Point cloud
//! @param Rr : Rain rate in mm/hr
//! @return Augmented point cloud
std::vector<std::vector<double>> Lisa::augment(std::vector<std::vector<double>>& pc, double Rr)
{
	// calculate extinction coefficient from the randomly sampled rainrate
	double alpha = calc_alpha(Rr);

	// Augment the lidar scan with the given extinction coefficient
	std::vector<std::vector<double>> pcnew(pc.size());
	for (size_t i = 0; i < pc.size(); i++) {
		pcnew[i] = augment(pc[i], alpha, Rr);
	}
	return pcnew;
}

//! @brief Augment a point (vector of doubles) given extinction coefficient alpha and a rain rate Rr. Although alpha can be calculated from rain rate Rr,
//! in many cases, the entire point cloud is augmented with a single rain rate and therefore re-calculating it each time would have been inefficient.
//! @param pc : Point (vector of doubles) of size 4: (x, y, z, reflectance)
//! @param alpha : Extinction coefficient in units of m^-1
//! @param Rr : Rain rate in units of mm/hr
//! @return Augmented point (vector of doubles) of size 5: (x_new, y_new, z_new, reflectance_new, label).
//! Label can take 3 values: 0 (lost point), 1 (randomly scattered point), 2 (original point with average effects)
std::vector<double> Lisa::augment(std::vector<double>& pc, double alpha, double Rr)
{
	double& x = pc[0];
	double& y = pc[1];
	double& z = pc[2];
	double& ref = pc[3];

	// calculate range and angles
	double range = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0)); 
	double phi   = atan2(y, x);
	double the   = acos(z / range);

	std::vector<double> pcnew(5, 0.0);

	// Calculate backreflected power from the target (in arbitrary units)
	double P0 = ref * exp(-2.0 * alpha * range) / pow(range, 2.0);

	// Beam calculations
	double Db   = range * tan(lidar.get_b_div());                   // beam diameter
	double bvol = 3.14159265359 * range * pow(Db / 2.0, 2.0) / 3.0; // conic beam volume
	double Nt   = bvol * pdist.N_total(Rr, dst);                    // total number of particles in the beam path

	// probabilistic rounding of the number of particles
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	double r   = dis(gen);
	int Nint   = (int) std::floor(Nt);
	double dif = Nt - (double) Nint;
	if (dif > r) ++Nint;

	// Monte Carlo: draw range from a quadratic PDF
	std::vector<double> ran_r;
	ran_r.reserve(Nt);
	double temp;
	for (size_t i = 0; i < Nt; i++) {
		temp = range * pow(dis(gen), 0.3333);
		if (temp > lidar.get_min_range()) ran_r.push_back(temp);
	}
	Nt = ran_r.size();

	// Monte Carlo: calculate back reflected power from randomly drawn particles
	double Pr(0), Prmax(0); // back reflected power for randomly drawn particles
	double ran_rmax(0);     // range of the particle with the largest return
	double Db_r(0);         // temp variable to hold beam diameter at a randomly drawn range

	if (Nt > 0) {
		std::vector<double> Dr = pdist.N_sample(Rr, dst, Nt); // diameter for randomly drawn particles
		for (size_t i = 0; i < Nt; i++) {
			Db_r = ran_r[i] * tan(lidar.get_b_div()); // beam diameter
			Pr   = ref_r * exp(-2.0 * alpha * ran_r[i]) * std::fmin(pow(Dr[i]/Db_r,2.0),1.0) / pow(ran_r[i], 2.0);

			// keep particle with the largest power
			if (Pr > Prmax) {
				Prmax    = Pr;
				ran_rmax = ran_r[i];
			}
		}
	}

	// Compare max random power to power from target
	double range_new(0), label(0), P_new(0), ref_new(0);

	if (range < lidar.get_min_range()) {
		// If range smaller than lidar min range and Prmax is larger than Pmin,
		// then randomly scatter a point else return all zeros
		if (Prmax > Pmin) {
			range_new = ran_rmax;
			P_new     = Prmax;
			ref_new   = ref_r * exp(-2.0 * alpha * range_new);
			label     = 1; // label for randomly scattered point
		} else {
			return pcnew;
		}
	} else {
		// If range is not smaller than min range
		if (P0 > Prmax) {
			// if object reflects more power than particles and SNR is larger than 1,
			// then keep the original points else return all zeros
			if (P0 > Pmin) {
				range_new = range;
				P_new     = P0;
				ref_new   = ref * exp(-2.0 * alpha * range_new);
				label     = 2; // label for original points
			} else {
				return pcnew;
			}
			
		} else {
			// if a random droplet backreflects more power than the object and that power is larger than Pmin,
			// then randomly scatter points else return all zeros
			if (Prmax > Pmin) {
				range_new = ran_rmax;
				P_new     = Prmax;
				ref_new   = ref_r * exp(-2.0 * alpha * range_new);
				label     = 1; // label for randomly scattered point
			} else {
				return pcnew;
			}
			
		}
	}

	// Add SNR dependent range uncertainty and return "fuzzy" point cloud
	double SNR = P_new / Pmin; // signal to noise ratio
	double sig = lidar.get_ran_uncer() / sqrt(2.0 * SNR); // range uncertainty
	std::normal_distribution<double> dist{ 0.0, sig }; // with zero mean and sig standard deviation
	range_new = range_new + dist(gen); // add range noise

	// Calculate new point cloud based on new range
	pcnew[0] = range_new * sin(the) * cos(phi);
	pcnew[1] = range_new * sin(the) * sin(phi);
	pcnew[2] = range_new * cos(the);
	pcnew[3] = ref_new;
	pcnew[4] = label;

	return pcnew;
}
