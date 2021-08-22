#ifndef MATERIAL_H
#define MATERIAL_H

#include <math.h>
#include <iostream>
#include <complex>

class Material {
public:
	//! @brief Calculates refractive index at a given wavelength.
	//! Dispersion curve for materials (Sellmeier, Lorentz-Drude, interpolated tables, etc) is implemented
	//! @param wavelength : Wavelength in micrometers
	//! @return Refractive index at given wavelength
	virtual std::complex<double> ref_ind(double wavelength) = 0;
};


class Water : public Material {
public:

	//! @brief Refractive index of water from a Sellmeier equation valid between [0.182, 1.129] um wavelength.
	//! Source: M. Daimon and A. Masumura. Measurement of the refractive index of distilled water
	//! from the near-infrared region to the ultraviolet region, Appl. Opt. 46, 3811-3820 (2007).
	//! https://refractiveindex.info/?shelf=main&book=H2O&page=Daimon-19.0C
	//! @param wavelength : Wavelength in micrometers
	//! @return Refractive index of water at given wavelength
	std::complex<double> ref_ind(double wavelength) override {
		// Make sure wavelength range is correct
		double lambda2 = wavelength * wavelength;
		if (wavelength <= 0.182) {
			lambda2 = 0.182 * 0.182;
			std::cout << "Wavelength needs to be larger than 0.182 um. Picking 0.182 um instead." << std::endl;
		}
		else if (wavelength >= 1.129) {
			std::cout << "Wavelength needs to be smaller than 1.129 um. Picking 1.129 um instead." << std::endl;
		}

		// Calculate permittivity from Sellmeier equation
		std::complex<double> epsilon = 1.0 + (5.672526103e-1 * lambda2) / (lambda2 - 5.085550461e-3)
				+ (1.736581125e-1 * lambda2) / (lambda2 - 1.814938654e-2)
				+ (2.121531502e-2 * lambda2) / (lambda2 - 2.617260739e-2)
				+ (1.138493213e-1 * lambda2) / (lambda2 - 1.073888649e1);
		return sqrt(epsilon);
	}
};

#endif // MATERIAL_H

