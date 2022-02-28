//
//  Merger_rate.h
//  PBH
//
//  Created by Wang Yifan on 9/20/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#ifndef Merger_rate_h
#define Merger_rate_h


#endif /* Merger_rate_h */


double PBH_Merger_rate(double z,double fraction,double M_phb_solarunit) // unit: SI
{
    double c5=c*c*c*c*c,fraction4=fraction*fraction*fraction*fraction;
    double Q,xbar,xbar4,t_c,T;  //the parameter in Sasaki's paper
    double z_eq,n_BH;           //z_eq: matter radiation equality epoch, n_BH: comoving number density of BH
    double cosmic_time_SI;         // the looking forward time
    double n_dp_dt;
    double M_pbh_SI = M_phb_solarunit * solar_mass;
    
    z_eq=2.4e4*Omega_m*h_0*h_0-1;           //http://www.tapir.caltech.edu/~chirata/ph217/lec06.pdf (eq.9)
    
    /* value of certain parameters */
    Q=3.0/170*pow(G*M_pbh_SI, -3.0)*c5;                                  //multiply c5 to correct the dimension
    xbar=1.0/(1+z_eq)/pow(fraction, 1.0/3)*pow(8.0*pi*G*M_pbh_SI/3/H_0_SI/H_0_SI/Omega_dm, 1.0/3);
                                                                    //physical mean separation, unit: m
    xbar4 = xbar*xbar*xbar*xbar;
    t_c = Q*xbar4*pow(fraction, 25.0/3);                              //the t_c in sasaki's paper
    T = xbar4*Q/fraction4;                                            //the capital T in sasaki's paper
    
    n_BH = 3.0*H_0_SI*H_0_SI/8/pi/G*fraction*Omega_dm/M_pbh_SI;        //number density in comoving volume
    
    cosmic_time_SI = z2t_SI(0, z);              //looking-forward time, unit:s
    
    if(cosmic_time_SI<t_c)
        n_dp_dt = n_BH*3.0/58*( -pow(cosmic_time_SI/T,3.0/8)+pow(cosmic_time_SI/T,3.0/37) ) / cosmic_time_SI;
    else
        n_dp_dt = n_BH*3.0/58* pow(cosmic_time_SI/T,3.0/8) * (-1.+pow(cosmic_time_SI/t_c, -29.0/56)*pow(fraction, -29.0/8))/cosmic_time_SI;

    return n_dp_dt;
}



