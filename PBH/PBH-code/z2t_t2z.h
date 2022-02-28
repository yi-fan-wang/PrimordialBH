//
//  z2t_t2z.h
//  SGWB
//
//  Created by Wang Yifan on 9/22/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#ifndef z2t_t2z_h
#define z2t_t2z_h


#endif /* z2t_t2z_h */

// (z_0, z_f) --> look back time
// convert z to a(scale factor) and do integral
// a trick for doing integral is narrow the integral interval by utilising the variable transformation


double z2t_integrand(double a)
{
    return 1/H_0_SI/sqrt(Omega_m/a+Omega_r/a/a+Omega_lambda*a*a);
}

// Look forward time , SI
double z2t_SI(double z_0,double z_f)
{
    double a_sup = 1./(1.+z_f); //integrate from 0 to a_sup
    double n = 1000;            //number of points in integrals
    double i,az,sum=0;
    for(i=0;i<=n;i++)
    {
        az = i*a_sup/n;         // divide a_sup into n segments
        sum += z2t_integrand(az);
    }
    sum = sum - 0.5*z2t_integrand(0) - 0.5*z2t_integrand(a_sup);    //trapzoid
    return sum*a_sup/n;         // sum should be multiplied by each segment's length
}

/*
// (z_0, z_f) --> t
// t=\int_{z_0}^{z_f} 1/(1+z)/H(z) dz
// adapative simpson method, sum = h/3 * [ f(a) + 4 *odd + 2 *even + f(b)]
double z2t_yr(double z_0,double z_f)// unit: Myr (consistent with the unit of H(z) )
{
    double sup=z_f;
    double inf=z_0;
    int i,n;
    double h,e,s;
    double t1,s1,s2,x;
    double eps;
    
    n=1;
    h=sup-inf;
    t1=h*(1.0/(1+inf)/H_z_peryr(inf)+1.0/(1+sup)/H_z_peryr(sup))/2.0;
    eps = 0.01 * 1e6;    //uncertainty: 0.01 Myr
    s1=t1;
    e=eps+1.0;
    while(e>=eps){
        s=0.0; // initialize s
        // the total number of nodes is 2n
        // the distance between nodes is h/2
        for(i=0;i<=n-1;i++){
            x=inf+(i+0.5)*h;
            s=s+1.0/(1+x)/H_z_peryr(x); // sum of the new nodes' value
            // all the new nodes are odd nodes
        }
        s2=(t1+2*h*s)/3.0;
        e=fabs(s2-s1);
        t1=(t1+h*s)/2.0;
        //t1 is an auxilary value, equals to: h/2 * [ f(a)/2 + f(1) + ... + f(2n-1) + f(b)/2 ]
        s1=s2;
        n=n+n;
        h=h/2.0;
        //printf("n=%d\n",n);
    }
    return(s2);
}*/

/*
// (t, z_0) --> z_f (suitable for z_f < 600 )
double t_yr2z(double z_0,double t,double up) // t in unit of yr because we use H_z_peryr(z)
{
    // use the value of upper limit and lower limit to squeeze the true value middle
    // 'middle' is like a tentative value, waiting to be test
    
    //double x0=10.0; // variable
    //double x1=z_0;  // lower limit
    //double x2=10.0; // upper limit
    
    
    //double middle = 10.0; //initialize middle = upper limit, becuase if no loop is implemented, we return middle after all
    
    double low = z_0; //lower limit
    //double up = 60000.0; // upper limit
    
    double middle = (low+up)/2; // first test value
    
    double eps = 1e-3;
    
    while (fabs(z2t_yr(z_0, middle)-t) > eps)
    {
        if (z2t_yr(z_0, middle) > t)
            up = middle;
        else
            low = middle;
    
        middle = (low+up)/2; // change the test value to a new, more accurate position
    }
    
    return middle;
}*/
