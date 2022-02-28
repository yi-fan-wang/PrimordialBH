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

#include "constant.h"
#include "z2t_t2z.h"
#include "my_dE_dnu.h"
#include "PBH_Merger_rate.h"

#include "PBH_Omega_gw.h"

void main(){

    double m_pbh = 10 * pow(2.,1.5);
    double m1 = 10,m2=10;
    double xi = 0;
    double frac = 2.560887e-02;
    double f;
    FILE *fp;
    //fp=fopen(name,"w");
    fp = fopen("Sasaki-localrate-1Msun.txt","w");
    
    for(f=-4;f<1;f=f+0.1)
        fprintf(fp,"%f %e\n",pow(10.,f),PBH_Merger_rate(0, pow(10.,f), 1)*yr2s*(pc2m*1e9)*(pc2m*1e9)*(pc2m*1e9));
}
