#include <iostream>
#include <Lidar.h>
#include <Material.h>
#include <Utils.h>
#include <Lisa.h>
#include <MiniLisa.h>

int main() {
	// test predefine lidars
	std::cout << "Test lidar" << std::endl;
	std::cout << VLS_128 << std::endl;
	
	// test water
	std::cout << "Test water dispersion" << std::endl;
	Water water;
	std::cout << "n = " << water.ref_ind(0.5876) << " which should be 1.3335" << std::endl;

	// test Mie
	double xm = 2.0 * 3.14159265359 * .525 / .6328;
	double qext = calc_qext(xm, 1.55);
	std::cout << std::endl << "QEXT = " << qext << " which should be 3.10543" << std::endl;

	xm = 1.0;
	std::complex<double> m(1.5, 1);
	qext = calc_qext(xm, m);
	std::cout << std::endl << "QEXT = " << qext << " which should be 2.336" << std::endl;

	// test trapezoidal rule
	std::vector<double> f(101), x(101);
	for (size_t i = 0; i < 101; i++) {
		x[i] = i * 0.01;
		f[i] = 2.0*x[i];
	}
	double tst = trapz(f, x);
	std::cout << std::endl << "Integral of 2x in range [0, 1] is " << tst << " which should be 1.0 " << std::endl << std::endl;
	
	
	// test logspace
	std::vector<double> tst2 = logspace(0.0, 1.0, 10);
	std::cout << "logspace(0,1,10) = ";
	for (int i = 0; i < 10; i++) {
		std::cout << tst2[i] << " ";
	}
	std::cout << std::endl << "Output should be: 1,  1.29154967,  1.66810054,  2.15443469,  2.7825594 , 3.59381366, 4.64158883, 5.9948425, 7.74263683, 10." << std::endl;

	// Test MP rain model
	MarshallPalmerRain rain;
	double Nd = rain.N_model(0.01, 30);
	std::cout << std::endl << "N_model(0.01, 30) = " << Nd << " which should be 7841.0256295643185" << std::endl;
	double Nt = rain.N_total(30, 0.050e-3);
	std::cout << std::endl << "N_total(30, 0.050e-3) = " << Nt << " which should be 3985.2723697019424" << std::endl;
	std::vector<double> dsam = rain.N_sample(30, 0.050e-3, 1000);
	double avg(0);
	for (size_t i = 0; i < 1000; i++) {
		avg += dsam[i];
	}
	avg /= 1000;
	std::cout << std::endl << "Average of N_sample(30, 0.050e-3, 1000) is " << avg << " which should be 0.50 +/- 0.02" << std::endl;

	// Test alpha
	MiniLisa augmenter(VLS_128, water, rain);
	double alpha = augmenter.calc_alpha(10.0); // m^-1
	alpha = alpha * 10000.0 / log(10.0); // convert to dB/km
	double alphat = 1.45 * pow(10.0,.64);
	std::cout << std::endl << "Alpha for 10 mm/hr rain is " << alpha << " which should be rougly " << alphat << std::endl;

	alpha = augmenter.calc_alpha(30.0); // m^-1
	alpha = alpha * 10000.0 / log(10.0); // convert to dB/km
	alphat = 1.45 * pow(30.0, .64);
	std::cout << std::endl << "Alpha for 30 mm/hr rain is " << alpha << " which should be rougly " << alphat << std::endl;

	// MiniLisa tests

	// test rmin
	std::vector<std::vector<double>> pc0(1);
	pc0[0] = std::vector<double>({ 1,1,1,0.5 });
	std::vector<std::vector<double>> tst_pc0 = augmenter.augment(pc0, 30);

	std::cout << std::endl << "When pc range is smaller than rmin: ";
	for (size_t i = 0; i < 4; i++) { std::cout << tst_pc0[0][i] << " "; } 
	std::cout << std::endl << "which should be all zeros";

	// test SNR calculation
	std::vector<std::vector<double>> pc1(1);
	pc1[0] = std::vector<double>({ 100,0,0,0.1 });
	std::vector<std::vector<double>> tst_pc1 = augmenter.augment(pc1, 30);
	std::cout << std::endl << "When SNR is smaller than 1, return all zeros: ";
	for (size_t i = 0; i < 4; i++) { std::cout << tst_pc1[0][i] << " "; }
	std::cout << std::endl;

	// test range uncertainty
	std::vector<std::vector<double>> pc2(1);
	pc2[0] = std::vector<double>({ 100,0,0,0.5 });
	std::vector<std::vector<double>> tst_pc2 = augmenter.augment(pc2, 30);
	double range = sqrt(pow(tst_pc2[0][0],2.0) + pow(tst_pc2[0][1],2.0) + pow(tst_pc2[0][2],2.0));
	std::cout << std::endl << "New range is " << range << " which should be 100 +/- ran uncer(sigma about 0.03 in this case)";

	// Lisa tests
	Lisa lisa(VLS_128, water, rain);

	// test rmin
	std::vector<std::vector<double>> tst_lisa0 = lisa.augment(pc0, 30);

	for (size_t i = 0; i < 5; i++) { std::cout << tst_lisa0[0][i] << " "; }
	std::cout << std::endl;

	// 
	std::vector<std::vector<double>> tst_lisa1 = lisa.augment(pc1, 30);
	for (size_t i = 0; i < 5; i++) { std::cout << tst_lisa1[0][i] << " "; }
	std::cout << std::endl;

	// test range uncertainty
	std::vector<std::vector<double>> tst_lisa2 = lisa.augment(pc2, 30);
	range = sqrt(pow(tst_lisa2[0][0], 2.0) + pow(tst_lisa2[0][1], 2.0) + pow(tst_lisa2[0][2], 2.0));
	std::cout << std::endl << "New range is " << range << std::endl;
	
	return 0;
}