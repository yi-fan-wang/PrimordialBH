//
//  PBH_Omega_gw.h
//  SGWB
//
//  Created by Wang Yifan on 9/27/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#ifndef PBH_Omega_gw_h
#define PBH_Omega_gw_h


#endif /* PBH_Omega_gw_h */

double PBH_Omega_gw_integrand(double f_obs,double z,double m1_pbh_solarunit,double m2_pbh_solarunit,double xi,double fraction)
// Rv(z)/(1+z)/E * dE/df(source frame)
// variable : f_obs and z
{
    double m_pbh_solarunit = m1_pbh_solarunit;             //这一行之前算错了,for pbh, the two component mass should be identical
    double f_source = f_obs * (1+z);
    double E = sqrt(Omega_r*(1+z)*(1+z)*(1+z)*(1+z) + Omega_m * (1+z) * (1+z) * (1+z) + Omega_lambda);
    return PBH_Merger_rate(z, fraction, m_pbh_solarunit)/(1+z)/E * dE_over_dnu(f_source, m1_pbh_solarunit, m2_pbh_solarunit,xi);
}

double PBH_Omega_gw(double freq, double fraction, double m1_pbh_solarunit,double m2_pbh_solarunit,double xi)
//local_merger_rate in unit of : Gpc^-3 yr^-1
{
    double rhoc = 3*H_0_SI*H_0_SI * c * c / 8 / pi / G ;
    double factor = freq / rhoc /H_0_SI;
    //printf("factor=%e\n",factor);
    
    double z_max = 1e6; // This is the z_max
    double sup,inf=0;
    double f_max = f_cut(m1_pbh_solarunit, m2_pbh_solarunit,xi); // Two BH have the same mass
    
    if (freq>f_max)
        return 0;
    if (freq<f_max/(1.0+z_max))
        sup=z_max;
    else
        sup=f_max/freq-1.0;
    
    int i,n;
    double h,e,s;
    double t1,s1,s2,x;
    double eps;
    
    // From here, I divide the interval into two part, z:0~1000 and z>1000
    double sup1,sup2;
    double inf1,inf2;
    double node,sum;
    double z_middle=1000;
    double n1=10000;    //number of interval z: 0~z_middle
    double n2=10000;   //number of interval z: >z_middle
    
    if(sup>z_middle){
        sup1 = z_middle;
        inf1 = 0;
        sum = 0;
        h = (sup1-inf1)/n1;
        for(i=0;i<=n1;i++)
        {
            node = inf1+i*h;
            sum += PBH_Omega_gw_integrand(freq, node, m1_pbh_solarunit, m2_pbh_solarunit,xi, fraction);
        }
        sum = sum-0.5*PBH_Omega_gw_integrand(freq, inf1, m1_pbh_solarunit, m2_pbh_solarunit,xi, fraction)-0.5*PBH_Omega_gw_integrand(freq, sup1, m1_pbh_solarunit, m2_pbh_solarunit,xi, fraction);
        sum = sum*h;
  // the z>1000 part, use the adaptive simpson method
    
        n=1;
        inf=n1;
        h=sup-inf;
        t1=h*(PBH_Omega_gw_integrand(freq, inf, m1_pbh_solarunit,m2_pbh_solarunit,xi, fraction)+PBH_Omega_gw_integrand(freq, sup, m1_pbh_solarunit, m2_pbh_solarunit,xi, fraction))/2.0;

        eps = t1/1e6;
        s1=t1;
        e=eps+1.0;
        while(e>eps&&n<n2){
            s=0.0;
            for(i=0;i<=n-1;i++){
                x=inf+(i+0.5)*h;
                s=s+PBH_Omega_gw_integrand(freq, x, m1_pbh_solarunit, m2_pbh_solarunit,xi, fraction);
            }
            s2=(t1+2*h*s)/3.0;
            e=fabs(s2-s1);
            t1=(t1+h*s)/2.0;
            s1=s2;
            n=n+n;
            h=h/2.0;
            printf("n=%d\n",n);
        }
        
        return (s2+sum) * factor;
    }
    else
    {   sup2 = sup;
        inf2 = 0;
        sum = 0;
        h = (sup2-inf2)/n1;
        for(i=0;i<=n1;i++)
        {
            node = inf2+i*h;
            sum += PBH_Omega_gw_integrand(freq, node, m1_pbh_solarunit, m2_pbh_solarunit,xi, fraction);
        }
        sum = sum-0.5*PBH_Omega_gw_integrand(freq, inf2, m1_pbh_solarunit, m2_pbh_solarunit,xi, fraction)-0.5*PBH_Omega_gw_integrand(freq, sup2, m1_pbh_solarunit,m2_pbh_solarunit, xi, fraction);
        sum = sum*h;
    
        return sum * factor;
    }
    
    // Omega_gw(f) = f * (localmergerrate) / rho_c / H_0 * \int Rv(z)/Rv(0) * 1/E * dE/df(sourceframe) dz
    
    
}

