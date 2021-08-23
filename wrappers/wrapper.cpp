#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include <Lisa.h>
#include <MiniLisa.h>
#include <ParticleDist.h>
#include <Material.h>
#include <Lidar.h>
#include <Utils.h>

namespace py = pybind11;

// Trampoline classes 
// see: https://pybind11.readthedocs.io/en/stable/advanced/classes.html
class PyParticleDist : public ParticleDist {
public:
	
    double N_model(double d, double Rr) override {
        PYBIND11_OVERRIDE_PURE(
            double, 		/* Return type */
            ParticleDist,   /* Parent class */
            N_model,        /* Name of function in C++ (must match Python name) */
            d, Rr           /* Argument(s) */
        );
    }
	
	double N_total(double Rr, double dst) override {
        PYBIND11_OVERRIDE_PURE(
            double, 		 /* Return type */
            ParticleDist,    /* Parent class */
            N_total,         /* Name of function in C++ (must match Python name) */
            Rr, dst          /* Argument(s) */
        );
    }
	
	std::vector<double> N_sample(double Rr, double dst, int N) override {
		PYBIND11_OVERRIDE_PURE(
            std::vector<double>,  /* Return type */
            ParticleDist,         /* Parent class */
            N_sample,             /* Name of function in C++ (must match Python name) */
            Rr, dst, N            /* Argument(s) */
        );
	}
};

class PyMaterial : public Material {
public:
    std::complex<double> ref_ind(double wavelength) override {
        PYBIND11_OVERRIDE_PURE(
            std::complex<double>, /* Return type */
            Material,             /* Parent class */
            ref_ind,              /* Name of function in C++ (must match Python name) */
            wavelength            /* Argument(s) */
        );
    }
};

PYBIND11_MODULE(pylisa, m) {
		
	// wrap mini lisa
	py::class_<MiniLisa>(m, "MiniLisa")
	.def(py::init<>())
	.def(py::init<Lidar &, Material &, ParticleDist &>())
	.def("augment",py::overload_cast<std::vector<std::vector<double>> &>(&MiniLisa::augment))
	.def("augment",py::overload_cast<std::vector<std::vector<double>> &, double>(&MiniLisa::augment))
	.def("calc_alpha",&MiniLisa::calc_alpha);
	
	// wrap lisa
	py::class_<Lisa, MiniLisa>(m, "Lisa")
	.def(py::init<>())
	.def(py::init<Lidar &, Material &, ParticleDist &>())
	.def("augment",py::overload_cast<std::vector<std::vector<double>> &>(&Lisa::augment))
	.def("augment",py::overload_cast<std::vector<std::vector<double>> &, double>(&Lisa::augment));

	
	// wrap particle dist
	py::class_<ParticleDist, PyParticleDist>(m, "ParticleDist")
	.def(py::init<>())
	.def("N_model",&ParticleDist::N_model)
	.def("N_total",&ParticleDist::N_total)
	.def("N_sample", &ParticleDist::N_sample);
	
	py::class_<MarshallPalmerRain, ParticleDist>(m, "MarshallPalmerRain")
	.def(py::init<>())
	.def("N_model",&MarshallPalmerRain::N_model)
	.def("N_total",&MarshallPalmerRain::N_total)
	.def("N_sample", &MarshallPalmerRain::N_sample);
	
	// wrap material
	py::class_<Material, PyMaterial>(m, "Material")
    .def(py::init<>())
	.def("ref_ind",&Material::ref_ind);
	
	py::class_<Water, Material>(m, "Water")
    .def(py::init<>())
    .def("ref_ind", &Water::ref_ind);
	
	// wrap lidar
	py::class_<Lidar>(m, "Lidar")
        .def(py::init<>())
        .def("get_ran_uncer", &Lidar::get_ran_uncer)
		.def("get_max_range", &Lidar::get_max_range)
		.def("get_min_range", &Lidar::get_min_range)
		.def("get_wavelength", &Lidar::get_wavelength)
		.def("get_b_div", &Lidar::get_b_div)
		.def("set_ran_uncer", &Lidar::set_ran_uncer)
		.def("set_max_range", &Lidar::set_max_range)
		.def("set_min_range", &Lidar::set_min_range)
		.def("set_wavelength", &Lidar::set_wavelength)
		.def("set_b_div", &Lidar::set_b_div)
		;
	
	// wrap utils
    m.def("trapz", &trapz,
	"Numerical integration using the trapezoidal rule");
	
	m.def("logspace", &logspace,
	"Calculate a vector of numbers evenly spaced in a logarithmic scale");

	m.def("calc_qext", &calc_qext,
	"Calculate extinction efficiency using Mie theory assuming the material is in air");
}