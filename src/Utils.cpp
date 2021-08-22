#include "Utils.h"

//! @brief Numerical integration using the trapezoidal rule: integrate f over x
//! @param f : Integrand (vector of doubles)
//! @param x : Independent variable (vector of doubles)
//! @return Integral of f over x
double trapz(std::vector<double> f, std::vector<double> x) {

	if (f.size() != x.size()) throw "f and x should have the same size";

	double out(0.0), delta(0.0);
	size_t size = f.size();
	for (size_t i = 0; i < size-1; i++) {
		delta = x[i + 1] - x[i];
		out   += (f[i + 1] + f[i]) * delta / 2.0;
	}

	return out;
}

//! @brief Calculate a vector of numbers evenly spaced in a logarithmic scale
//! @param start : exponent starting value (i.e start=0 corresponds to 10^0 = 1)
//! @param end : exponent ending value (i.e end=2 corresponds to 10^2 = 100)
//! @param N : Number of points in range [10^start, 10^end]
//! @return  A vector of doubles evenly spaced in a logarithmic scale
std::vector<double> logspace(double start, double end, int N) {
	double delta = (end - start)/((double) N - 1);
	std::vector<double> out(N);
	for (int i = 0; i < N; i++) {
		out[i] = pow(10.0, i*delta + start);
	}
	return out;
}

//! @brief Calculate extinction efficiency using Mie theory assuming the material is in air. Method is modified from
//! Bohren and Huffman, "Absorption and Scattering of Light by Small Particles"
//! @param x : "Size parameter" defined as x = k*a where k is the wave-vector and a is the particle diameter
//! @param nref : Refractive index of the inclusion
//! @return Extinction efficiency
double calc_qext(double x, std::complex<double> nref) {
	// modified from Bohren's implementation
	using namespace std::complex_literals;

	// optimal number of iterations from Wiscombe
	int nstop = (int)std::round(x + 4.0 * pow(x, 0.3333) + 2.0);
	std::complex<double> psi0, psi1, psi, chi0, chi1, chi, xi1, xi;

	// downward recurrence
	int nmx = (int)std::round(fmax(nstop, std::abs(nref) * x) + 15.0);
	std::vector<std::complex<double>> D(nmx);
	D[nmx - 1] = 0.0 + 1i * 0.0;

	std::complex<double> rho = nref * x;
	std::complex<double> temp = 0.0 + 1i * 0.0;
	for (int n = nmx - 1; n > 0; n--) {
		temp = (n + 1.0) / rho;
		D[n - 1] = temp - 1.0 / (D[n] + temp);
	}

	// upward recurrence
	psi0 = cos(x);   psi1 = sin(x);
	chi0 = -sin(x);  chi1 = cos(x);
	xi1 = psi1 - 1i * chi1;

	double qext(0), en(0), fn(0), tmp1(0), tau(0), p0(0), p1(1), p(0);
	std::complex<double> an, bn, tmp2;
	std::complex<double> s1_1 = 0.0 + 1i * 0.0;

	for (int n = 0; n < nstop; n++) {
		en = n + 1.0;
		fn = (2.0 * en + 1.0) / (en * (en + 1.0));

		tmp1 = (2.0 * en - 1.0);
		psi = tmp1 * psi1 / x - psi0;
		chi = tmp1 * chi1 / x - chi0;
		xi = psi - 1i * chi;

		tmp2 = (D[n] / nref + en / x);
		an = tmp2 * psi - psi1;
		an = an / (tmp2 * xi - xi1);

		tmp2 = (nref * D[n] + en / x);
		bn = tmp2 * psi - psi1;
		bn = bn / (tmp2 * xi - xi1);

		// Calculate s1_1 (we only need s1_1 for qext)
		tau = en * p1 - (en + 1.0) * p0;
		s1_1 += fn * (an * p1 + bn * tau);

		// prev becomes next
		psi0 = psi1;
		psi1 = psi;
		chi0 = chi1;
		chi1 = chi;
		xi1 = psi1 - 1i * chi1;
		p = p1; // save in temporary variable
		p1 = ((2.0 * en + 1.0) * p1 - (en + 1.0) * p0) / en;
		p0 = p;
	}
	qext = 4 * std::real(s1_1) / (x * x);
	return qext;
}