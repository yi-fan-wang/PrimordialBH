//
//  astro_Omega_gw.h
//  PBH
//
//  Created by Wang Yifan on 10/9/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#ifndef astro_Omega_gw_h
#define astro_Omega_gw_h


#endif /* astro_Omega_gw_h */

double Omega_gw_bbh_deltamass_integrand(double f_obs,double z,double delay_time_Myr,double m1_solarunit,double m2_solarunit,double xi)
// Rv(z)/(1+z)/E * dE/df(source frame)
// variable : f_obs and z
{
    double f_source = f_obs * (1+z);
    double E = sqrt(Omega_m * (1+z) * (1+z) * (1+z) + Omega_lambda);
    //double Rv_SI = Rv(z, delay_time_Myr) / (1e6 * pc2m)/ (1e6 * pc2m)/ (1e6 * pc2m)/yr2s  ;
    return Rv(z, delay_time_Myr)/E * dE_over_dnu(f_source, m1_solarunit, m2_solarunit,xi);
}

double Omega_gw_bbh_deltamass(double f, double local_merger_rate_ligounit, double delay_time_Myr, double m1_solarunit, double m2_solarunit,double xi)
//local_merger_rate in unit of : Gpc^-3 yr^-1
{
    double local_merger_rate_SI = local_merger_rate_ligounit / (1e9 * pc2m) / (1e9 * pc2m) / (1e9 * pc2m) / yr2s;
    double rhoc = 3*H_0_SI*H_0_SI * c * c / 8 / pi / G ;
    double factor = f / rhoc * local_merger_rate_SI/H_0_SI;
    // unit : same with SFR, M_solar Mpc^-3 yr^-1
    //double Rv0_SI = Rv(0, delay_time_Myr)*solar_mass/( 1e6 * pc2m)/( 1e6 * pc2m)/( 1e6 * pc2m)/yr2s;
    
    double z_max = 10;
    double sup,inf=0;
    double f_max = f_cut(m1_solarunit, m2_solarunit,xi);
    //printf("fcut = %f\n",f_max);
    if (f>f_max)
        return 0;
    if (f<f_max/(1.0+z_max))
        sup=z_max;
    else
        sup=f_max/f-1.0;
    
    int i,n;
    double h,e,s;
    double t1,s1,s2,x;
    double eps = t1/1e5;//3 order of magnitude smaller than the estimated order(i.e estimated value:0.01)
    
    n=1;
    h=sup-inf;
    t1=h*(Omega_gw_bbh_deltamass_integrand(f, inf, delay_time_Myr, m1_solarunit, m2_solarunit,xi)+Omega_gw_bbh_deltamass_integrand(f, sup, delay_time_Myr, m1_solarunit, m2_solarunit,xi))/2.0;
    //printf("f=%f,t1=%e\n",f,t1);
    eps = t1/1e4;
    s1=t1;
    e=eps+1.0;
    while(e>=eps&&n<4096){
        s=0.0;
        for(i=0;i<=n-1;i++){
            x=inf+(i+0.5)*h;
            s=s+Omega_gw_bbh_deltamass_integrand(f, x, delay_time_Myr, m1_solarunit, m2_solarunit,xi);
        }
        s2=(t1+2*h*s)/3.0;
        e=fabs(s2-s1);
        t1=(t1+h*s)/2.0;
        s1=s2;
        n=n+n;
        h=h/2.0;
    }
    // printf("n=%d,s2=%e\n",n,s2);
    return s2/ Rv(0, delay_time_Myr) * factor;
    
    // Omega_gw(f) = f * (localmergerrate) / rho_c / H_0 * \int Rv(z)/Rv(0) * 1/E * dE/df(sourceframe) dz
    // Rv(z) = \int sfr(z)*metal_fraction(z)/(1+z_f)/(1+z_f)/H(z)
    
    
}
