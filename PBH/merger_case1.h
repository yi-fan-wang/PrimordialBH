//
//  merger_case1.h
//  PBH
//
//  Created by Wang Yifan on 20/03/2017.
//  Copyright © 2017 Wang Yifan. All rights reserved.
//

#ifndef merger_case1_h
#define merger_case1_h

#include <stdio.h>

#endif /* merger_case1_h */

double PBH_Merger_rate_case1(double z,double fraction,double alpha, double beta, double M1,double M2,double M3,double Mbar) // unit: SI
{
    // This is the first try of calculating PBH merger rate in the existence of tidal force enhancement.
    // We here only consider the 1<M</xi/f and M^(10/7)</xi/f case
    /* alpha and beta are numerical coefficients */
    double xi =2*Mbar/(M1+M2);
    double eta = 2*M3/(M1+M2);
    double N = alpha *xi;
    double M = beta * eta;
    double c5=c*c*c*c*c,fraction4=fraction*fraction*fraction*fraction;
    
    double cosmic_time_SI = z2t_SI(0, z);              //looking-forward time, unit:s
    double dp_dt_1,dp_dt_2;                             //dp_dt_1 is for b<a, dp_dt_2 is for b>a
    double M_pbh_SI = Mbar * solar_mass;
    
    double z_eq=2.4e4*Omega_m*h_0*h_0-1;           //z_eq: matter radiation equality epoch
                                                //http://www.tapir.caltech.edu/~chirata/ph217/lec06.pdf (eq.9)
    
    /* value of certain parameters */
    /* !!! Note here you should use SI mass !!! */
    double Q=3.0/170*pow(G*solar_mass, -3.0)*2./M1/M2/(M1+M2)*c5;       //multiply c5 to correct the dimension
    double xbar=1.0/(1+z_eq)/pow(fraction, 1.0/3)*pow(8.0*pi*G*M_pbh_SI/3/H_0_SI/H_0_SI/Omega_dm, 1.0/3);
    
    //physical mean separation, unit: m
    double xbar4 = xbar*xbar*xbar*xbar;
    double Tprime = xbar4*Q*pow(N,4.0)/fraction4;
    
    //some time point
    double tc = Tprime*pow(M, 7.)*pow(fraction/xi,37./3);
    double td = Tprime * pow(M, -3.)*pow(fraction/xi,16./3);
    double te = pow(alpha, 4.)*xbar4*Q*(fraction/xi,4./3);
    
    printf("%e\n",Q);
    double n_BH = 3.0*H_0_SI*H_0_SI/8/pi/G*fraction*Omega_dm/M_pbh_SI;        //number density in comoving volume
    
    double dp_dt_1_factor = 3./58*M*pow(cosmic_time_SI/Tprime,3./8)/cosmic_time_SI;
    double dp_dt_2_factor = 3./2*M*pow(cosmic_time_SI/Tprime, 3./8)/cosmic_time_SI;
    if(cosmic_time_SI<tc)
        dp_dt_1 = dp_dt_1_factor * (pow(pow(M, 32./37)*pow(cosmic_time_SI/Tprime, 6./37), -29./16) - pow(M*M, -29./16));
    else if(cosmic_time_SI<te)
        dp_dt_1 = dp_dt_1_factor * (pow(pow(fraction/xi, -32./21)*pow(cosmic_time_SI/Tprime, 2./7), -29./16) - pow(M*M, -29./16));
    else
        dp_dt_1 = 0;
    
    if(cosmic_time_SI<td)
        dp_dt_2 = dp_dt_2_factor * (pow(1./M/M, -1./16) - pow(1., -1./16));
    else if(cosmic_time_SI<te)
        dp_dt_2 = dp_dt_2_factor * (pow(pow(fraction/xi, -28./9)*pow(cosmic_time_SI/Tprime, 2./3), -1./16) - pow(1., -1./16));
    else
        dp_dt_2 = 0;
    
    
    return n_BH*(dp_dt_1+dp_dt_2);
}
