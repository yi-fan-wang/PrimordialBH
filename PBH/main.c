//
//  main.c
//  SGWB
//
//  Created by Wang Yifan on 7/3/16.
//  Copyright © 2016 Wang Yifan. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_sf_erf.h>

#include "constant.h"
#include "z2t_t2z.h"
//#include "LIGO_dE_dnu.h"
#include "my_dE_dnu.h"
#include "PBH_Merger_rate.h"
//#include "merger_case1.h"
//#include "lensing_merger_rate.h"
#include "PBH_Omega_gw.h"
//#include "Unequal_pbh_merger_rate.h"
//#include "Unequal_PBH_Omega.h"

#include "constraint.h"
//#include "astro_merger_rate.h"
//#include "astro_Omega_gw.h"

void main(){
    //double m1=36.2,m2=29.1;
    //double m_pbh = 100*pow(2.,1./5);
    //double xi = -0.06;      //PBH effective spin
    //double frac1 = 2.46e-3,frac2 = 1.01e-3, frac3=4.8e-4; //28.1 mass PBH
    
    //double f;
    //double m_pbh = 10 * pow(2.,1.5);
    //double m1 = 10,m2=10;
    //double xi = 0;
    //double frac = 2.560887e-02;
    //double f;
    //double m1=14.2,m2=7.5;  //PBH mass=32.3
    //double m_pbh = 8.9*pow(2., 1./5);
    //double frac1 = 6.35e-3, frac2=2.66e-3,frac3=7.5e-4;  // 8.9 mass PBH
    //double xi = 0.21;
    
    double m1=0.01,m2=0.01;
    double frac=0.0348173;
    double f;
    //double m1 = 30, m2 =30, m3 =30, mbar = 30,z,f;
    //double fraction = 1e-3;
    //double alpha = 1, beta =1.1;
    //double unit_con = yr2s*(pc2m*1e9)*(pc2m*1e9)*(pc2m*1e9);
    printf("%e\n",PBH_Merger_rate(1, 0.1, 1));
    //char name[100] = "";
    //sprintf(name,"alpha-%.1f-beta-%.1f-f-%.3f-m-%.1f.txt",alpha,beta,fraction,mbar);
    //FILE *fp;
    //fp=fopen(name,"w");
    //fp = fopen("Wangsai.txt","w");
    //fp = fopen("O2-results/constraints3.txt", "w");
    //fp=fopen("O2-results/0p1Msun.txt","w");
    //for(f=-4;f<1;f=f+0.1)
    //    fprintf(fp,"%f %e\n",pow(10.,f),PBH_Merger_rate(0, pow(10.,f), 1)*yr2s*(pc2m*1e9)*(pc2m*1e9)*(pc2m*1e9));
    //for(f=0;f<6;f=f+0.1)
    //    fprintf(fp, "%f %e\n",pow(10.,f),PBH_Omega_gw(pow(10.,f), frac, m1, m2, 0.));
    /*
    double upper = 3e-8;
    frac=1.;
    double m1ind;
    for(m1ind=-3;m1ind<-1;m1ind=m1ind+0.1)
    {
        m1 = pow(10., m1ind);
        printf("%f\n",m1);
        if (m1<100)
        {
            while(PBH_Omega_gw(25., frac, m1, m1, 0.)>upper)
            frac =frac-0.001;
        }
        else
        {
            frac=1.;
            while(PBH_Omega_gw(25., frac, m1, m1, 0.)>upper)
            frac =frac-0.001;
        }
        fprintf(fp,"%f %f\n",m1,frac);
    }*/



}



/*
Chirp mass = 28.1, mass = 32.3, Local merger rate = 3.4(0.6~12) ==> fraction
                                                                    2.46e-3         12.268764
                                                                    1.01e-3         3.4377
                                                                    4.8e-4           0.59
Chipr mass = 8.9 , component mass = 10.2, Local merger rate = 37(6~129)  ==> fraction             local mergerate
                                                                                6.35e-3             128.995
                                                                                2.66e-3             37.09
                                                                                7.5e-4             5.96
 Fig 4: 0.1M_odot f=0.06
        1,10,20, f= Huang paper

*/
