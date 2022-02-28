//
//  constraint.h
//  PBH
//
//  Created by Wang Yifan on 10/10/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#ifndef constraint_h
#define constraint_h


#endif /* constraint_h */

void constraint(FILE *fp)
{
    double frac1,frac2,frac3;
    double o1=3.25e-8,o2=5.17e-9,o5=6.55e-10;   // The minimum of the sensitivity curve
    double f1 = 28., f2 = 28., f3 = 27.;        // The corresponding freq of the local minimum of the sens curve
    double index,mass;
    double low,up,middle;
    double eps = 1e-3;
    
    for(index=2;index<=2.7;index=index+0.1)
    {
        mass=pow(10., index);
        printf("%f\n",index);
        low = 0.; //lower limit
        up = 1.; // upper limit
        middle = (low+up)/2.; // first test value
        
        while (fabs((o1-PBH_Omega_gw(28., middle, mass, mass,0))/o1) > eps && middle<=0.9999)
        {
            if (PBH_Omega_gw(28., middle, mass, mass,0) > o1)
                up = middle;
            else
                low = middle;
            
            middle = (low+up)/2; // change the test value to a new, more accurate position
            //printf("o1 middle=%f\n",middle);
        }
        
        frac1=middle;
        
        low = 0; //lower limit
        up = frac1; // upper limit
        middle = (low+up)/2; // first test value
        
        while (fabs((o2-PBH_Omega_gw(28, middle, mass, mass,0))/o2) > eps)
        {
            if (PBH_Omega_gw(28, middle, mass, mass,0) > o2)
                up = middle;
            else
                low = middle;
            
            middle = (low+up)/2; // change the test value to a new, more accurate position
             //printf("o2 middle=%f\n",middle);
        }
        
        frac2=middle;
        
        low = 0; //lower limit
        up = frac2; // upper limit
        middle = (low+up)/2; // first test value
        
        while (fabs((o5-PBH_Omega_gw(27, middle, mass, mass,0))/o5) > eps)
        {
            if (PBH_Omega_gw(27, middle, mass, mass,0) > o5)
                up = middle;
            else
                low = middle;
            
            middle = (low+up)/2; // change the test value to a new, more accurate position
            //printf("o5 middle=%f\n",middle);
        }
        
        frac3=middle;
        
        
        fprintf(fp,"%f %e %e %e\n",mass,frac1,frac2,frac3);
    }
}



//O1: f=28      Omega=3.25e-8
//O2: f=28      Omega=5.17e-9
//O5: f=27      Omega=6.55e-10
