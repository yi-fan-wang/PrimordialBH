//
//  Unequal_pbh_merger_rate.h
//  PBH
//
//  Created by Wang Yifan on 19/01/2017.
//  Copyright © 2017 Wang Yifan. All rights reserved.
//

#ifndef Unequal_pbh_merger_rate_h
#define Unequal_pbh_merger_rate_h


#endif /* Unequal_pbh_merger_rate_h */

double unequal_pbh_merger_rate(double mbar_solar,double fraction,double m1_solar,double m2_solar,double m3_solar,double z)
{
    //currently, a problem is to calculate the beta, beta initially is 0.8, and now is 50, we observe the difference
    double alpha=0.4,beta=25,zeta=2*mbar_solar/(m1_solar+m2_solar),eta=2*m3_solar/(m1_solar+m2_solar);
    double paraN=alpha*zeta,paraM=eta*beta;
    double alpha4=alpha*alpha*alpha*alpha;
    double paraM7=paraM*paraM*paraM*paraM*paraM*paraM*paraM;
    double paraN4=paraN*paraN*paraN*paraN;
    
    double c5=c*c*c*c*c,fraction4=fraction*fraction*fraction*fraction;
    double Q,xbar,xbar4,t_c,T;  //the parameter in Sasaki's paper
    double z_eq,n_BH;           //z_eq: matter radiation equality epoch, n_BH: comoving number density of BH
    double cosmic_time_SI;         // the looking forward time
    double n_dp_dt;
    double mbar_SI = mbar_solar * solar_mass;
    
    z_eq=2.4e4*Omega_m*h_0*h_0-1;           //http://www.tapir.caltech.edu/~chirata/ph217/lec06.pdf (eq.9)
    
    /* value of certain parameters */
    Q=3.0/170*pow(G*mbar_SI, -3.0)*c5;                                  //multiply c5 to correct the dimension
    xbar=1.0/(1+z_eq)/pow(fraction, 1.0/3)*pow(8.0*pi*G*mbar_SI/3/H_0_SI/H_0_SI/Omega_dm, 1.0/3);
    //physical mean separation, unit: m
    xbar4 = xbar*xbar*xbar*xbar;
    t_c = Q*xbar4*pow(fraction, 25.0/3)*alpha4*pow(zeta, -25.0/3)*paraM7;                              //the t_c in sasaki's paper
    T = xbar4*Q/fraction4;                                            //the capital T in sasaki's paper
    
    n_BH = 3.0*H_0_SI*H_0_SI/8/pi/G*fraction*Omega_dm/mbar_SI;        //number density in comoving volume
    
    cosmic_time_SI = z2t_SI(0, z);              //looking-forward time, unit:s
    //printf("%e,%e,%e\n",cosmic_time_SI,t_c,pow(paraM,-58.0/37));
    
    if(cosmic_time_SI<t_c)
    n_dp_dt = n_BH*3.0/58 * pow(cosmic_time_SI/T,3.0/8) / cosmic_time_SI * paraM*pow(paraN,-3.0/2) * (pow(paraM,-58.0/37)* pow(cosmic_time_SI/T/paraN4,-87.0/296)-pow(paraM,-29.0/8) );
    else
    n_dp_dt = n_BH*3.0/58 * pow(cosmic_time_SI/T,3.0/8) / cosmic_time_SI * paraM*pow(paraN,-3.0/2) * (pow(cosmic_time_SI/paraN4/T,-29.0/56)*pow(fraction/zeta,58.0/21)-pow(paraM,-29.0/8));
    
    return n_dp_dt;
}

double gauss_massdist(double mbar,double sigma,double m)
{
    //return 1.0/sqrt(2*pi)/sigma*exp(-(m-mbar)*(m-mbar)/2/sigma/sigma);
    return 1./1000;
}
