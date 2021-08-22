#ifndef LIDAR_H
#define LIDAR_H

#include <iostream>

class Lidar {
public:
	//! @brief Construct Lidar object with default values:
	//! range uncertainty = +/- 0.06 m
	//! maximum range = 200 m
	//! minimum range = 1.5 m
	//! wavelength = 0.903 um
	//! beam divergence = 0.5 mrad
	Lidar() : ran_uncer(0.06), max_range(200.0), min_range(1.5),
		wavelength(0.903), b_div(0.5e-3) {}

	//! @brief Constructs Lidar object with given parameters:
	//! @param ran_uncer : Range uncertainty in meters
	//! @param max_range : Maximum range in meters
	//! @param min_range : Minimum range in meters
	//! @param wavelength : Lidar wavelength in micro meters
	//! @param b_div : Beam divergence in radians
	Lidar(double ran_uncer, double max_range, double min_range, double wavelength, double b_div):
		ran_uncer(ran_uncer),
		max_range(max_range),
		min_range(min_range),
		wavelength(wavelength),
		b_div(b_div)
	{}

	//! @brief Constructs Lidar object with given parameters assuming beam divergence is 0.5 milliradians:
	//! @param ran_uncer : Range uncertainty in meters
	//! @param max_range : Maximum range in meters
	//! @param min_range : Minimum range in meters
	//! @param wavelength : Lidar wavelength in micro meters
	Lidar(double ran_uncer, double max_range, double min_range, double wavelength):
		ran_uncer(ran_uncer),
		max_range(max_range),
		min_range(min_range),
		wavelength(wavelength),
		b_div(0.5e-3)
	{}
	
	//! @brief Get lidar range uncertainty
	//! @return Range uncertainty in meters
	double get_ran_uncer() { return ran_uncer; }

	//! @brief Get maximum range
	//! @return Maximum lidar range in meters
	double get_max_range() { return max_range; }

	//! @brief Get minimum range
	//! @return Minimum range in meters
	double get_min_range() { return min_range; }

	//! @brief Get lidar laser wavelength
	//! @return Lidar laser wavelength in micrometers
	double get_wavelength() { return wavelength; }

	//! @brief Get lidar laser beam divergence
	//! @return Beam divergence in radians
	double get_b_div() { return b_div; }

	//! @brief Set lidar range uncertainty
	//! @param x : Range uncertainty in meters
	void set_ran_uncer(double x) { ran_uncer = x; }

	//! @brief Set maximum range
	//! @param x : Maximum lidar range in meters
	void set_max_range(double x) { max_range = x; }

	//! @brief Set minimum range
	//! @param x : Minimum range in meters
	void set_min_range(double x) { min_range = x; }

	//! @brief Set lidar laser wavelength
	//! @param x : Lidar laser wavelength in micrometers
	void set_wavelength(double x) { wavelength = x; }

	//! @brief Set lidar laser beam divergence
	//! @param x : Beam divergence in radians
	void set_b_div(double x) { b_div = x; }

	//! @brief Print Lidar object properties
	void print() {
		std::cout << "Range uncertainty [m] " << ran_uncer << std::endl
			<< "Max range [m] " << max_range << std::endl
			<< "Min range [m] " << min_range << std::endl
			<< "Laser wavelength [um] " << wavelength << std::endl
			<< "Beam divergence [rad] " << b_div << std::endl;
	}

	friend std::ostream& operator<<(std::ostream& os, const Lidar&);

private:
	//! @brief Range uncertainty in meters
	double ran_uncer;

	//! @brief Maximum range in meters
	double max_range;

	//! @brief Minimum range in meters
	double min_range;

	//! @brief Laser wavelength in micrometers
	double wavelength;

	//! @brief Beam divergence in radians
	double b_div;
};

/*
Define some lidars

Following data is taken from:
Carballo et al., "LIBRE: The Multiple 3D LiDAR Dataset"
https://arxiv.org/pdf/2003.06129.pdf
*/

// Veldoyne
static Lidar VLS_128  = Lidar(0.06, 245.0, 3.0, 0.903);
static Lidar HDL_64S2 = Lidar(0.04, 120.0, 3.0, 0.903);
static Lidar HDL_32E  = Lidar(0.04, 100.0, 2.0, 0.903);
static Lidar VLP_32C  = Lidar(0.06, 200.0, 1.0, 0.903);
static Lidar VLP_16   = Lidar(0.06, 100.0, 1.0, 0.903);

// Hesai
static Lidar PANDAR_64  = Lidar(0.04, 200.0, 0.3, 0.905);
static Lidar PANDAR_40P = Lidar(0.04, 200.0, 0.3, 0.905);

// Ouster
static Lidar OS1_64 = Lidar(0.06, 120.0, 0.8, 0.850);
static Lidar OS1_16 = Lidar(0.06, 120.0, 0.8, 0.850);

// RoboSense
static Lidar RS_LIDAR32 = Lidar(0.06, 200.0, 0.4, 0.905);

#endif // LIDAR_H