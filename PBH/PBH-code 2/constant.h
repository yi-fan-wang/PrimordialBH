//
//  constant.h
//  SGWB
//
//  Created by Wang Yifan on 9/21/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#ifndef constant_h
#define constant_h


#endif /* constant_h */

/* constant */
double const pi = 3.1415926536;
double const G = 6.6726E-11;                  //SI unit
double const c = 299792458;                 //m*s^-1

/* cosmological parameters */
double const h_0=0.678;       //Hubble constant in unit: 100 km/s/Mpc
double const Omega_m =0.3070;  //matter content
double const Omega_r = 4.165e-5/h_0/h_0;    //radiation content (http://www.astro.ucla.edu/%7Ewright/CC.python)
double const Omega_b = 0.045;
double const Omega_dm = 0.27;
double const Omega_lambda=1-Omega_m-Omega_r;    //dark energy
double const universe_age_yr = 13.813e9;    //universe age = 1.381784e+10, according to current program
                                                //universe age = 13.809 Gyr, http://www.astro.ucla.edu/%7Ewright/CosmoCalc.html
                                                //universe age = 13.813 Gyr, planck 2015
double const solar_mass = 1.989E30;         //SI unit: kg

/* transformation coefficient */
double const yr2s = 3.154e7; //1yr = 3.154e7 s  (http://www.wolframalpha.com/input/?i=1+year+%3D+s)
double const pc2m = 3.086e16; //1parsec = 3.086e16 m (http://www.wolframalpha.com/input/?i=1+parsec+in+m)

/* eps */
//double const eps_z2t_yr = 0.1 * 1e6;

double H_0_SI = 100*h_0 * 1000 / (1e6 * pc2m);
//Hubble parameter H(z), unit: s^-1
double H_z_SI(double z)
{
     //H(z) = H(0) * sqrt{ Omega_M * (1+z)^3 + Omega_lambda }
    return 100*h_0 * 1000 / (1e6 * pc2m) * sqrt(Omega_lambda+Omega_m*(1+z)*(1+z)*(1+z)+Omega_r*(1+z)*(1+z)*(1+z)*(1+z));
}

//Hubble parameter H(z), unit: yr^-1
double H_z_peryr(double z)
{
    /*
     H(z) = H(0) * sqrt{ Omega_R * (1+z)^4 + Omega_M * (1+z)^3 + Omega_lambda }
     
     1 parsec = 3.26 light year = 3.26 * speed of light m/s * 1year
     H_0 = 100h_0 km * s^-1 * Mpc^-1 = 100h_0 km * s^-1 * (3.26 * speed of light * 1e6 yr)^-1
         = 100h_0 *1000 / speed of light / 3.26 / 1e6
     */
    return 100*h_0*1000 / c / 3.26  / 1e6 * sqrt(Omega_lambda+Omega_m*(1+z)*(1+z)*(1+z)+Omega_r*(1+z)*(1+z)*(1+z)*(1+z));
}
