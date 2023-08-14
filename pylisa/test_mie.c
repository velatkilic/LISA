#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mie.h"

// test functions should return void and take no arguments
typedef void (*test_func_t)();

// printf decorator for testing functions
void test_func(test_func_t func, char *msg) {
    printf("Testing: %s", msg);
    func();
    printf("\nPASSED\n--------------------\n");
}

// assert equal within a threshold
void assert_almost_equal(double a, double b, double eps) {
    assert(fabs(a - b) < eps);
}

// test extinction efficiency calculation
void test_calc_qext() {
    float xm = 2.0 * 3.14159265359 * .525 / .6328;
    complex ref_rel = Cform(1.55, 0.);
    int n_ang = 2;
    complex S1[10]; complex S2[10];
    float q_ext, q_sca, q_back, g_aniso = 0.;
    BHMie(&xm, &ref_rel, &n_ang, S1, S2, &q_ext, &q_sca, &q_back, &g_aniso);
    printf("\nqext = %f\n", q_ext);
    assert_almost_equal(q_ext, 3.10543, 1e-5);
}

// test extinction efficiency calculation
void test_calc_BHMie_for_large_size_param() {
    float xm = 20000.;
    complex ref_rel = Cform(1.55, 0.);
    int n_ang = 2;
    complex S1[10]; complex S2[10];
    float q_ext, q_sca, q_back, g_aniso = 0.;
    BHMie(&xm, &ref_rel, &n_ang, S1, S2, &q_ext, &q_sca, &q_back, &g_aniso);
    // original code causes a segfault here, solution was to dynamically allocate cd
}

int main() {
    printf("\n--------------------\n");
    test_func(test_calc_qext, "extinction coefficient for a=.525, lam=.6328, nref=1.55");
    test_func(test_calc_BHMie_for_large_size_param, "Large size param does not cause segfault");
    return 0;
}