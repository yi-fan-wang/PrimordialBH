
//Star formation rate R*(z), unit: M_solar*Mpc^-3*yr^-1

/*
double sfr(double z)    //SFR, HB06
{
    double a,b,c,d;
    a=0.0170;
    b=0.13;
    c=3.3;
    d=5.3;
    return (a+b*z)*h_0/(1+pow(z/c, d));
}*/

double sfr(double z)    // SFR, Vangioni
{
    double nu,a,b,z_m;
    nu = 0.146;
    a = 2.80;
    b = 2.46;
    z_m = 1.72;
    return nu * a * exp(b*(z-z_m)) / (a - b + b * exp(a * (z - z_m)) );
}

double metal(double z)  // the metallicity < Z_solar/2 fraction
{
    double z_solar = 0.02,z_sigma=0.5;
    double y=0.019,R=0.27,rhob=2.77e11*Omega_b*h_0*h_0;
    double factor = y*(1-R)/rhob*1e6*pc2m/1000/yr2s/(h_0*100);
    
    double sup=20;
    double inf=z;
    
    int i,n;
    double h,e,s,eps;
    double t1,s1,s2,x;
    
    // integrand: sfr(z)/(1+z)/E(z)
    n=1;
    h=sup-inf;
    t1=h*(sfr(sup)/(1+sup)/sqrt((1+sup)*(1+sup)*(1+sup)*Omega_m+Omega_lambda)+sfr(inf)/(1+inf)/sqrt((1+inf)*(1+inf)*(1+inf)*Omega_m+Omega_lambda))/2.0;
    eps = t1/1e5;
    s1=t1;
    e=eps+1.0;
    while(e>=eps){
        s=0.0;
        for(i=0;i<=n-1;i++){
            x=inf+(i+0.5)*h;
            s=s+sfr(x)/(1+x)/sqrt((1+x)*(1+x)*(1+x)*Omega_m+Omega_lambda);
        }
        s2=(t1+2*h*s)/3.0;
        e=fabs(s2-s1);
        t1=(t1+h*s)/2.0;
        s1=s2;
        n=n+n;
        h=h/2.0;
    }
    return 0.5+0.5*gsl_sf_erf(log10(z_solar/2/(s2*factor*3))/sqrt(2)/z_sigma);//return (1+erf)/2
}

//the integrand of Rv, i.e. :
//Rv = \int R*(z_f) / [ (1+z_f) * H(z_f) * t_delay ]
//with the assumption of p(t_delay)= 1 / t_delay
double Rv_integrand(double z_f,double z)
{
    //z_f = z at formation epoch,   z = variable = the epoch emitting GW
    
    /*  Note here the unit of star formation rate could be not consistent with the delay time, but doesn't make a difference finally
     SFR unit: M_solar Mpc^-3 yr^-1
     delay time unit: yr (Sometimes could be Myr, but doesn't matter, cancels with H(z) unit
        But delay time unit must be consistent with H(z) unit)
     */
    if(metal(z_f)>1)
        printf("ERROR!");
    return sfr(z_f)*0.0000000000001/(1+z_f)/(1+z_f)/H_z_peryr(z_f)/(z2t_SI(0, z_f)/yr2s-z2t_SI(0, z)/yr2s);

    
    //It doesn't matter what unit we chose for H(z), because we'll calculte the dimensionless Rv(z)/Rv(0);
}


double Rv(double z,double delay_time_Myr)
// unit : same with SFR, M_solar Mpc^-3 yr^-1
//  discretize the interval by hand, because the function is too sharp!
//  trapzoid segmentation
{   double delay_time_yr = delay_time_Myr * 1e6;    // Delay time!!!    unit: Myr
    double inf,sup;
    double h1=0,h2=0,h3=0;//divide the interval to 3 parts
    double fx;
    
   // double norm_factor = log(universe_age_yr)-log(delay_time_yr);
    /*  normalization factor of delay time. 
        Don't understand why the maximus delay time is set to be Hubble time. Shouldn't it be Hubble time minus the merger time(lookback)?
     */
    
    double norm_factor = (log(universe_age_yr-z2t_SI(0,z)/yr2s)-log(delay_time_yr));
    
    int i;
    double sum1=0,sum2=0,sum3=0;
    
    /* integral variable : the epoch of formation
        lower limit: the epoch of [t(z) + minimum delay time ](looking back time)
        upper limit: large enough, z=8 maybe enough
     */
    //sup=6;
    inf=t_yr2z(0,z2t_yr(0, z)+delay_time_yr,10);
    
    /* We will integrate Rv(z) again, but only from z=0~6, so in this function, the input z will not exceed 6, so it's appropriate to divide the interval into three part by hand, which are [inf, inf + 0.1] [ inf + 0.1, inf + 1.1]
        [int + 1.1, 8]
     */
    
    h1=0.1/100;//first interval: z = [inf , inf + 0.1]
    for(i=0;i<100;i++)
    {   fx=Rv_integrand(inf+i*h1, z);
        sum1=sum1+fx;
    }
    
    h2=1.0/100;//second interval: z = [inf + 0.1 , inf + 1.1]
    for(i=0;i<100;i++)
    {   fx=Rv_integrand(inf+0.1+i*h2, z);
        sum2=sum2+fx;
    }
    
    h3=(12.0-(inf+1.1))/100;//second interval: z = [inf + 1.1 , 12]
    for(i=0;i<100;i++)
    {   fx=Rv_integrand(inf+1.1+i*h3, z);
        sum3=sum3+fx;
    }
    
    sum1=sum1-0.5*Rv_integrand(inf, z)-0.5*Rv_integrand(inf+0.1, z);
    sum2=sum2-0.5*Rv_integrand(inf+0.1, z)-0.5*Rv_integrand(inf+1.1, z);
    sum3=sum3-0.5*Rv_integrand(inf+1.1, z)-0.5*Rv_integrand(12.0, z);
    
    return (sum1*h1+sum2*h2+sum3*h3)/norm_factor;
    //return (log(universe_age_myr)-log(delay_time))/(log(universe_age_myr-z2t(0,z))-log(delay_time));
}

/*      The following is for binary neutron star
double Iv(double z,double delay_time_Myr)//unit: same with SFR: M_solar Mpc^-3 yr^-1
{
    double E = sqrt(Omega_m * (1+z) * (1+z) * (1+z) + Omega_lambda);
    return Rv(z, delay_time_Myr)/E/cbrt(1+z);
}

//the integral of Iv(z)
//adaptive simpson method
// two variable, one is observed frequency, another is f_max, depends on M_c
double Integral_Iv(double f,double f_max,double delay_time_Myr)
{
    if(f>f_max)
        return 0;//However, the input f has already been made smaller than f_max
    
    double z_max = 6;
    
    double sup;
    double inf=0;
    
    int i,n;
    double h,e,s;
    double t1,s1,s2,x;
    double eps = 1e-5;//3 order of magnitude smaller than the estimated order(i.e estimated value:0.01)
    
    if (f<f_max/(1.0+z_max))
        sup=z_max;
    else
        sup=f_max/f-1.0;
    
    n=1;
    h=sup-inf;
    t1=h*(Iv(inf,delay_time_Myr)+Iv(sup,delay_time_Myr))/2.0;
    s1=t1;
    e=eps+1.0;
    while(e>=eps){
        s=0.0;
        for(i=0;i<=n-1;i++){
            x=inf+(i+0.5)*h;
            s=s+Iv(x,delay_time_Myr);
        }
        s2=(t1+2*h*s)/3.0;
        e=fabs(s2-s1);
        t1=(t1+h*s)/2.0;
        s1=s2;
        n=n+n;
        h=h/2.0;
    }
    //printf("s2=%e\n",s2);
    return(s2);
}   */


