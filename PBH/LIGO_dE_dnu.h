//
//  LIGO_dE_dnu.h
//  SGWB
//
//  Created by Wang Yifan on 9/27/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#ifndef LIGO_dE_dnu_h
#define LIGO_dE_dnu_h


#endif /* LIGO_dE_dnu_h */

double dE_over_dnu(double freq, double m1_solarunit,double m2_solarunit,double xi)
{
    double LAL_MTSUN_SI = 4.9254923218988636432342917247829673e-6; /**< Geometrized solar mass, s. = LAL_MSUN_SI / LAL_MPL_SI * LAL_TPL_SI */
    double M_tot = m1_solarunit + m2_solarunit;
    double M_c_SI = pow(m1_solarunit * m2_solarunit,3./5)/pow(m1_solarunit + m2_solarunit,1./5)*solar_mass;
    double etha = (m1_solarunit*m2_solarunit)/M_tot/M_tot;
    
    double dE_over_dnu;
    double factor = pow(G*pi,2.0/3) * pow(M_c_SI, 5.0/3)/3;
    
    double v1 = 404.*((0.66389*etha*etha-0.10321*etha+0.10979)/0.125481)*(20./M_tot);
    double v2 = 807.*((1.3278*etha*etha-0.20642*etha+0.21957)/0.250953)*(20./M_tot) ;
    double fgmax = 1153.*((1.7086*etha*etha-0.26592*etha+0.28236)/0.322668)*(20./M_tot) ;
    double sigma = 237.*((1.1383*etha*etha-0.177*etha+0.046834)/0.0737278)*(20./M_tot);
    
    double v = pow(pi*M_tot*freq*LAL_MTSUN_SI,1./3);
    double v_1 = pow(pi*M_tot*v1*LAL_MTSUN_SI,1./3);
    double v_2 = pow(pi*M_tot*v2*LAL_MTSUN_SI,1./3);
    
    double mod1 = pow((1 + (-323./224. + 451.*etha/168)*v*v + xi*(27./8 - 11.*etha/6)*v*v*v),2.);
    double mod2 = pow((1 + (1.4547*xi - 1.8897)*v + (-1.8153*xi + 1.6557)*v*v),2.);
    
    /* mod1_v1,mod1_v2,mod2_v1,mod_v2 should be constant  */
    double mod1_v1 = pow((1 + (-323./224. + 451.*etha/168)*v_1*v_1 + xi*(27./8 - 11.*etha/6)*v_1*v_1*v_1),2.);
    double mod2_v1 = pow((1 + (1.4547*xi - 1.8897)*v_1 + (-1.8153*xi + 1.6557)*v_1*v_1),2.);
    double mod1_v2 = pow((1 + (-323./224. + 451.*etha/168)*v_2*v_2 + xi*(27./8 - 11.*etha/6)*v_2*v_2*v_2),2.);
    double mod2_v2 = pow((1 + (1.4547*xi - 1.8897)*v_2 + (-1.8153*xi + 1.6557)*v_2*v_2),2.);

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
        dE_over_dnu = 0;
    return (dE_over_dnu);
}

double f_cut(double m1_solarunit,double m2_solarunit,double xi)
{
    //double LAL_MTSUN_SI = 4.9254923218988636432342917247829673e-6; /**< Geometrized solar mass, s. = LAL_MSUN_SI / LAL_MPL_SI * LAL_TPL_SI */
    double M_tot = m1_solarunit + m2_solarunit;
    //double M_c_SI = pow(m1_solarunit * m2_solarunit,3./5)/pow(m1_solarunit + m2_solarunit,1./5)*solar_mass;
    double etha = (m1_solarunit*m2_solarunit)/M_tot/M_tot;
    
    double v1 = 404.*((0.66389*etha*etha-0.10321*etha+0.10979)/0.125481)*(20./M_tot);
    double v2 = 807.*((1.3278*etha*etha-0.20642*etha+0.21957)/0.250953)*(20./M_tot) ;
    double fgmax = 1153.*((1.7086*etha*etha-0.26592*etha+0.28236)/0.322668)*(20./M_tot) ;
    double sigma = 237.*((1.1383*etha*etha-0.177*etha+0.046834)/0.0737278)*(20./M_tot);
    
    //printf("LIGO code, v1=%f,v2=%f,fmax=%f,sigma=%f\n",v1,v2,fgmax,sigma);
    return fgmax;

}

double my_spin_f_cut(double m1_solarunit,double m2_solarunit,double chi)
{
    double M_c_solarunit = pow(m1_solarunit*m2_solarunit, 3.0/5)/pow(m1_solarunit + m2_solarunit, 1.0/5);
    double M_total_solarunit = m1_solarunit + m2_solarunit;
    double sratio = m1_solarunit * m2_solarunit / M_total_solarunit / M_total_solarunit ;
    
    //    double v = pow(pi*M_total_solarunit*solar_mass*f,1.0/3);//v is dimensionless by assuming G=c=1
    //    double alpha2 = -323.0/224 + 451.0 * sratio/168;
    //    double alpha3 = (27.0/8 - 11.0 * sratio/6)/chi;         //chi is also dimensionless
    //    double epsilon1 = 1.4547*chi - 1.8897;
    //    double epsilon2 = -1.8153*chi + 1.6557;
    //    double addition1 = (1 + alpha2*v*v + alpha3*v*v*v) * (1 + alpha2*v*v + alpha3*v*v*v);
    //    double addition2 = (1 + epsilon1*v + epsilon2*v*v) * (1 + epsilon1*v + epsilon2*v*v);
    
    double mu_k0_f1 = 1-4.455*pow(1-chi, 0.217)+3.521*pow(1-chi, 0.26);
    double mu_k0_f2 = (1-0.63*pow(1-chi, 0.3))/2;
    double mu_k0_sigma = (1-0.63*pow(1-chi, 0.3))*pow(1-chi, 0.45)/4;
    double mu_k0_f3 = 0.3236+0.04894*chi+0.01346*chi*chi;
    
    double y_10_f1=0.6437,y_11_f1=0.827,y_12_f1=-0.2726,y_20_f1=-0.05822,y_21_f1=-3.935,y_30_f1=-7.092;
    double y_10_f2=0.1469,y_11_f2=-0.1228,y_12_f2=-0.02609,y_20_f2=-0.0249,y_21_f2=0.1701,y_30_f2=2.325;
    double y_10_sigma=-0.4098,y_11_sigma=-0.03523,y_12_sigma=0.1008,y_20_sigma=1.829,y_21_sigma=-0.02017,y_30_sigma=-2.87;
    double y_10_f3=-0.1331,y_11_f3=-0.08172,y_12_f3=0.1451,y_20_f3=-0.2714,y_21_f3=0.1279,y_30_f3=4.922;
    
    double unit_correction = c*c*c/G/solar_mass;
    
    double f_merg=(mu_k0_f1+y_10_f1*sratio+y_11_f1*sratio*chi+y_12_f1*sratio*chi*chi+y_20_f1*sratio*sratio+y_21_f1*sratio*sratio*chi+y_30_f1*sratio*sratio*sratio)/pi/M_total_solarunit*unit_correction;
    
    double f_ring=(mu_k0_f2+y_10_f2*sratio+y_11_f2*sratio*chi+y_12_f2*sratio*chi*chi+y_20_f2*sratio*sratio+y_21_f2*sratio*sratio*chi+y_30_f2*sratio*sratio*sratio)/pi/M_total_solarunit*unit_correction;
    
    double sigma=(mu_k0_sigma+y_10_sigma*sratio+y_11_sigma*sratio*chi+y_12_sigma*sratio*chi*chi+y_20_sigma*sratio*sratio+y_21_sigma*sratio*sratio*chi+y_30_sigma*sratio*sratio*sratio)/pi/M_total_solarunit*unit_correction;
    
    double f_cut=(mu_k0_f3+y_10_f3*sratio+y_11_f3*sratio*chi+y_12_f3*sratio*chi*chi+y_20_f3*sratio*sratio+y_21_f3*sratio*sratio*chi+y_30_f3*sratio*sratio*sratio)/pi/M_total_solarunit*unit_correction;
    
    printf("my   code, v1=%f,v2=%f,fmax=%f,sigma=%f\n",f_merg,f_ring,f_cut,sigma);
    return f_cut;
}

double nospin_fcut(double m1_solarunit,double m2_solarunit)
{
    double a_k1=2.9740e-1,b_k1=4.4810e-2,c_k1=9.5560e-2;
    double a_k2=5.9411e-1,b_k2=8.9794e-2,c_k2=1.9111e-1;
    double a_k3=5.0801e-1,b_k3=7.7515e-2,c_k3=2.2369e-2;
    double a_k4=8.4845e-1,b_k4=1.2848e-1,c_k4=2.7299e-1;
    double unit_correction = c*c*c/G/solar_mass;
    
    double M_c_solarunit = pow(m1_solarunit*m2_solarunit, 3.0/5)/pow(m1_solarunit + m2_solarunit, 1.0/5);
    double M_total_solarunit = m1_solarunit + m2_solarunit;
    double sratio = m1_solarunit * m2_solarunit / M_total_solarunit / M_total_solarunit ;
    
    double f_merg=(a_k1*sratio*sratio+b_k1*sratio+c_k1)/pi/M_total_solarunit*unit_correction;
    double f_ring=(a_k2*sratio*sratio+b_k2*sratio+c_k2)/pi/M_total_solarunit*unit_correction;
    double sigma=(a_k3*sratio*sratio+b_k3*sratio+c_k3)/pi/M_total_solarunit*unit_correction;
    double f_cut=(a_k4*sratio*sratio+b_k4*sratio+c_k4)/pi/M_total_solarunit*unit_correction;
    
    printf("nospin my, v1=%f,v2=%f,fmax=%f,sigma=%f\n",f_merg,f_ring,f_cut,sigma);
    return 0;
}
