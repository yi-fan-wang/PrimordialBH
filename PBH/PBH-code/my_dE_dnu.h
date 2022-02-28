//
//  my_dE_dnu.h
//  PBH
//
//  Created by Wang Yifan on 10/9/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#ifndef my_dE_dnu_h
#define my_dE_dnu_h


#endif /* my_dE_dnu_h */

double dE_over_dnu(double freq, double m1_solarunit,double m2_solarunit,double chi)
{
    double M_c_solarunit = pow(m1_solarunit*m2_solarunit, 3.0/5)/pow(m1_solarunit + m2_solarunit, 1.0/5);
    double M_c_SI = M_c_solarunit * solar_mass;
    double M_total = m1_solarunit + m2_solarunit;
    double sratio = m1_solarunit * m2_solarunit / M_total / M_total ;
    
    double dE_over_dnu;
    double factor = pow(G*pi,2.0/3) * pow(M_c_SI, 5.0/3)/3;
    
    double mu_k0_f1 = 1-4.455*pow(1-chi, 0.217)+3.521*pow(1-chi, 0.26);
    double mu_k0_f2 = (1-0.63*pow(1-chi, 0.3))/2;
    double mu_k0_sigma = (1-0.63*pow(1-chi, 0.3))*pow(1-chi, 0.45)/4;
    double mu_k0_f3 = 0.3236+0.04894*chi+0.01346*chi*chi;
    
    double y_10_f1=0.6437,y_11_f1=0.827,y_12_f1=-0.2726,y_20_f1=-0.05822,y_21_f1=-3.935,y_30_f1=-7.092;
    double y_10_f2=0.1469,y_11_f2=-0.1228,y_12_f2=-0.02609,y_20_f2=-0.0249,y_21_f2=0.1701,y_30_f2=2.325;
    double y_10_sigma=-0.4098,y_11_sigma=-0.03523,y_12_sigma=0.1008,y_20_sigma=1.829,y_21_sigma=-0.02017,y_30_sigma=-2.87;
    double y_10_f3=-0.1331,y_11_f3=-0.08172,y_12_f3=0.1451,y_20_f3=-0.2714,y_21_f3=0.1279,y_30_f3=4.922;
    
    double unit_correction = c*c*c/G/solar_mass;
    
    double v1=(mu_k0_f1+y_10_f1*sratio+y_11_f1*sratio*chi+y_12_f1*sratio*chi*chi+y_20_f1*sratio*sratio+y_21_f1*sratio*sratio*chi+y_30_f1*sratio*sratio*sratio)/pi/M_total*unit_correction;
    
    double v2=(mu_k0_f2+y_10_f2*sratio+y_11_f2*sratio*chi+y_12_f2*sratio*chi*chi+y_20_f2*sratio*sratio+y_21_f2*sratio*sratio*chi+y_30_f2*sratio*sratio*sratio)/pi/M_total*unit_correction;
    
    double sigma=(mu_k0_sigma+y_10_sigma*sratio+y_11_sigma*sratio*chi+y_12_sigma*sratio*chi*chi+y_20_sigma*sratio*sratio+y_21_sigma*sratio*sratio*chi+y_30_sigma*sratio*sratio*sratio)/pi/M_total*unit_correction;
    
    double fgmax=(mu_k0_f3+y_10_f3*sratio+y_11_f3*sratio*chi+y_12_f3*sratio*chi*chi+y_20_f3*sratio*sratio+y_21_f3*sratio*sratio*chi+y_30_f3*sratio*sratio*sratio)/pi/M_total*unit_correction;
    
    double v = pow(pi*M_total*freq/unit_correction,1./3);
    double v_1 = pow(pi*M_total*v1/unit_correction,1./3);
    double v_2 = pow(pi*M_total*v2/unit_correction,1./3);
    
    double mod1 = pow((1. + (-323./224. + 451.*sratio/168)*v*v + chi*(27./8. - 11.*sratio/6)*v*v*v),2.);
    double mod2 = pow((1. + (1.4547*chi - 1.8897)*v + (-1.8153*chi + 1.6557)*v*v),2.);
    
    /* mod1_v1,mod1_v2,mod2_v1,mod_v2 should be constant  */
    double mod1_v1 = pow((1. + (-323./224. + 451.*sratio/168.)*v_1*v_1 + chi*(27./8. - 11.*sratio/6.)*v_1*v_1*v_1),2.);
    double mod2_v1 = pow((1. + (1.4547*chi - 1.8897)*v_1 + (-1.8153*chi + 1.6557)*v_1*v_1),2.);
    double mod1_v2 = pow((1. + (-323./224. + 451.*sratio/168.)*v_2*v_2 + chi*(27./8. - 11.*sratio/6.)*v_2*v_2*v_2),2.);
    double mod2_v2 = pow((1. + (1.4547*chi - 1.8897)*v_2 + (-1.8153*chi + 1.6557)*v_2*v_2),2.);
    
    double transition1 = 1.0/v1 * mod1_v1/mod2_v1;
    double transition2 = 1.0/v1 * pow(v2,-4.0/3) * mod2_v2*mod1_v1/mod2_v1;
    
    double fcbrt = cbrt(freq);
    double f23 = fcbrt*fcbrt;
    
    if(freq<v1)
        dE_over_dnu = factor * 1./fcbrt * mod1;
    else if (freq<v2)
        dE_over_dnu = factor * transition1 * f23 * mod2;
    else if (freq<fgmax)
    {
        dE_over_dnu = factor * transition2 * freq * freq / pow((1+4*(freq-v2)*(freq-v2)/sigma/sigma),2.);
        //printf("%e\n",v2);
    }
    else
        dE_over_dnu = 0.0;
    return (dE_over_dnu);
}

double f_cut(double m1_solarunit,double m2_solarunit,double chi)
{
    double M_total = m1_solarunit + m2_solarunit;
    double sratio = m1_solarunit * m2_solarunit / M_total / M_total ;
    
    double mu_k0_f3 = 0.3236+0.04894*chi+0.01346*chi*chi;
    
    double y_10_f3=-0.1331,y_11_f3=-0.08172,y_12_f3=0.1451,y_20_f3=-0.2714,y_21_f3=0.1279,y_30_f3=4.922;
    
    double unit_correction = c*c*c/G/solar_mass;
    double f_cut=(mu_k0_f3+y_10_f3*sratio+y_11_f3*sratio*chi+y_12_f3*sratio*chi*chi+y_20_f3*sratio*sratio+y_21_f3*sratio*sratio*chi+y_30_f3*sratio*sratio*sratio)/pi/M_total*unit_correction;
    return f_cut;
    
}
