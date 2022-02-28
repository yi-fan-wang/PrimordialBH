//
//  Unequal_PBH_Omega.h
//  PBH
//
//  Created by Wang Yifan on 19/01/2017.
//  Copyright © 2017 Wang Yifan. All rights reserved.
//

#ifndef Unequal_PBH_Omega_h
#define Unequal_PBH_Omega_h


#endif /* Unequal_PBH_Omega_h */

double unequal_pbh_omega_integrand(double mbar_solar,double sigma, double fraction,double f_obs,double m1_solar,double m2_solar,double m3_solar,double z)
// R_PBH(z)/(1+z)/E * dE/df(source frame)*F(M1)*F(M2)*F(M3)
// variable : f_obs and z
{
    double xi=0;
    double f_source = f_obs * (1+z);
    double E = sqrt(Omega_r*(1+z)*(1+z)*(1+z)*(1+z) + Omega_m * (1+z) * (1+z) * (1+z) + Omega_lambda);
    //printf("  mergerrate=%e\n",unequal_pbh_merger_rate(mbar_solar, fraction, m1_solar, m2_solar, m3_solar, z));
    return unequal_pbh_merger_rate(mbar_solar, fraction, m1_solar, m2_solar, m3_solar, z)/(1+z)/E * dE_over_dnu(f_source, m1_solar, m2_solar,xi) * gauss_massdist(mbar_solar, sigma, m1_solar) * gauss_massdist(mbar_solar, sigma, m2_solar) * gauss_massdist(mbar_solar, sigma, m3_solar);
}

double unequal_pbh_omega(double mbar_solar,double sigma,double fraction,double f_obs)
{
    double rhoc = 3*H_0_SI*H_0_SI * c * c / 8 / pi / G ;
    double factor = f_obs / rhoc /H_0_SI;
    
    double xi=0;
    
    int n_m=20;
    double inf_m=mbar_solar-3*sigma;
    double sup_m=mbar_solar+3*sigma;// The upper limit and lower limit for the mass integration
    double h_m = (sup_m-inf_m)/n_m;
    
    int n_z=300;     // number of the discretization interval
    double inf_z=0;
    double sup_z;
    double z_max=30;   // The upper limit and lower limit for the redshift integration, z_max = 30, because we are interested in the region f=[10Hz,300Hz]. The source locating larger than z=30 will not contribute to this region
    double h_z;
    
    double f_max;
    
    int i,j,k,l;
    double m3,m2,m1,z;
    double sum1,sum2,sum3,sum4;


    
    sum4=0;
    for(i=0;i<=n_m;i++) // Only m3 is variable
    {
        printf("i=%d\n",i);
        
        m3=inf_m + h_m * i;
        
        sum3=0;
        for(j=0;j<=n_m;j++) // (m3,m2) are variables
        {
            m2 = inf_m + h_m * j;
            
            sum2=0;
            for(k=0;k<=n_m;k++) // (m3,m2,m1) are variables
            {
                m1 = inf_m + h_m * k;
                // The integration of z
                f_max = f_cut(m1, m2,xi); // The cut-off frequency of the binary with mass m1 and m2, no spin
                
                if (f_obs>f_max)
                    sum1 = 0;
                else
                {
                if (f_obs<f_max/(1.0+z_max))
                    sup_z=z_max;
                else
                    sup_z=f_max/f_obs-1.0;
                h_z = (sup_z-inf_z)/n_z;
                sum1 = 0;
                
                for(l=0;l<=n_z;l++)
                {
                    z = inf_z +h_z * l;
                    
                    if(l==0||l==n_z)
                        sum1 += 1./2*unequal_pbh_omega_integrand(mbar_solar, sigma, fraction, f_obs, m1, m2, m3, z);
                    else
                        sum1 += unequal_pbh_omega_integrand(mbar_solar, sigma, fraction, f_obs, m1, m2, m3, z);
                    //printf("sup_z=%f,z=%f,sum1=%e\n",sup_z,z,sum1);
                }
                sum1 = sum1 * h_z;
                }
                
                if(k==0||k==n_m)
                sum2 += 1./2*sum1;
                else
                sum2 += sum1;
            }
            sum2 = sum2*h_m;
            
            if(j==0||j==n_m)
                sum3 += 1./2*sum2;
            else
                sum3 += sum2;
        }
        sum3 = sum3*h_m;
        
        
        if(i==0||i==n_m)
            sum4 += 1./2*sum3;
        else
            sum4 += sum3;
    }
    sum4 = sum4*h_m;
    
    return factor*sum4;
}
