#ifndef __NAUNET_DATA_H__
#define __NAUNET_DATA_H__

// Struct for holding the nessesary additional variables for the problem.
struct NaunetData {
    // clang-format off
    double nH;
    double Tgas;
    double zeta_cr = 1.300e-17;
    double zeta_xr;
    double Tdust = 1.500e+01;
    double Av = 1.000e+00;
    double G0 = 1.000e+00;
    double omega = 5.000e-01;
    double rG = 1.000e-05;
    double barr = 1.500e-08;
    double sites = 1.000e+15;
    double hop = 3.000e-01;
    double nMono = 2.000e+00;
    double duty = 3.160e-19;
    double Tcr = 7.000e+01;
    double branch = 1.000e-02;
    double opt_frz = 1.000e+00;
    double opt_thd = 1.000e+00;
    double opt_uvd = 1.000e+00;
    double opt_crd = 1.000e+00;
    double opt_rcd = 1.000e+00;
    
    // clang-format on
};

#endif