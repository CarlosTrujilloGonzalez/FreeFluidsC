/*
 * FFtests.c
 *
 *  Created on: 17/04/2015
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2018  Carlos Trujillo Gonzalez

 *This program is free software; you can redistribute it and/or
 *modify it under the terms of the GNU General Public License version 3
 *as published by the Free Software Foundation

 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.

 *You should have received a copy of the GNU General Public License
 *along with this program; if not, write to the Free Software
 *Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include <math.h>
#include <stdio.h>
#include "FFeosPure.h"
#include "FFeosMix.h"
#include "FFtools.h"
#include "FFactivity.h"

//Some test with a single substance
void TestEosPure(){
    //first we define some variables for the tests
    enum FF_EOS eos=FF_IdealGas;
    enum FF_EosType eosType;
    FF_ThermoProperties thR;
    double refT=298.15;
    double refP=1.01325e5;
    double liqFraction=10;
    double Arr,Z;
    FF_CubicParam paramC;
    double Tb,Vp;
    char option;
    char state;
    char var;//the type of calculation from pressure: H=P,H, U=P,U, S=P,S
    double answerL[3],answerG[3];
    double answer2[2],answer3[3],answer4[4],answer6[6];
    //double answer2[2],answer3[3],aFCl[numSubs],aFCg[numSubs],aACl[numSubs],aACg[numSubs];

    //Some data for single substance test
    //***********************************

    //First some cubic eos data. Uncomment the one you want to use
    FF_CubicEOSdata data1C={0,FF_PRSV1,18,647.1,22060000,0.23,0.344,0,0,-0.0655854,0,0,0};//EOS constants for water PRSV1
    //FF_CubicEOSdata data1C={0,FF_PRMELHEM,18,647.1,22060000,0.23,0.344,0,0,0.88917,-0.00155,0,0};//EOS constants for water PRMELHEM
    //FF_CubicEOSdata data1C={0,FF_PRMC,18,647.1,22060000,0.23,0.344,0,0,0.919393,-0.355014,0.38726,0};//EOS constants for water PRMC
    //FF_CubicEOSdata data1C={0,FF_PRALMEIDA,18,647.1,22060000,0.23,0.344,0,0,0.81473,0.02707,0.96611,0};//EOS constants for water ALMEIDA
    //FF_CubicEOSdata data1C={0,FF_SRKTWU91,18,647.1,22060000,0.23,0.344,0,0,0.413297,0.874988,2.19435,0};//EOS constants for water SRKTWU91
    //FF_CubicEOSdata data1C={0,FF_PRSV1,58.12,425.1,3797000,0.0,0.199,0,0,-0.020735,0,0,0};//EOS constants for butane PRSV1
    //FF_CubicEOSdata data1C={0,FF_PRMC,58.12,425.1,3797000,0.0,0.199,0,0,0.702394,-0.418049,1.287094,0};//EOS constants for butane PRMC
    //FF_CubicEOSdata data1C={0,FF_PRALMEIDA,58.12,425.1,3797000,0.0,0.199,0,0,0.58249,0.06849,0.98301,0};//EOS constants for butane PRALMEIDA
    //FF_CubicEOSdata data1C={0,FF_PRSV1,106.165,616.89,3534600,0,0.325,0.0161,0,0,0};//m-xylene PRSV1
    //FF_CubicEOSdata data1C={0,FF_SRKSOF,106.165,616.89,3534600,0,0.325,-2.026244,2.300990,0,0};//m-xylene SRKSOF

    //Now some data for PCSAFT eos. The same
    FF_SaftEOSdata data1S={0,FF_PCSAFT,18,647.1,22060000,0.23,0,2.6273,1.5,180.3,0.0942,1804.22,0,0,2,2,0};//EOS constants for water PCSAFT 4C
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,18,647.1,22060000,0.23,0,2.6273,1.5,180.3,0.1799,1804.22,0,0,2,2,0};//EOS constants for water PCSAFT 4C
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,18,647.1,22060000,0.23,0,2.229,2.1945,141.6,0.2039,1804.17,0,0,2,2,0};//EOS constants for water PCSAFT 4C
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_JC,18,647.1,22060000,0.23,0,2.9657,1.0405,175.15,0.08924,2706.7,1.85,0.66245,1,1,0};//EOS constants for water 2B, polar Al Saifi JC
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_GV,18,2.7920,1.2255,203.00,0.07172,2406.7,1.85,0,1,1,0,647.1,22060000,0.23,0.344,-0.0655854,0};//EOS constants for water 2B, polar Al Saifi GV
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,32.05,2.6510,2.7920,186.60,0.14600,2090.20,0,0,1,1,0,512.5,80.8e5,0,0.556,-0.168,0.179};//Methanol 2B, non polar
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_JC,32.05,3.1369,1.7266,168.84,0.06311,2585.9,1.7,0.3513,1,1,0,512.5,80.8e5,0,0.556,-0.168,0.179};//Methanol 2B, polar JC
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_JC,46.068,516.2,6383475,0.248,0.635,3.2774,2.2049,187.24,0.03363,2652.7,1.7,0.29466,1,1,0};//PCSAFT constants for ethanol 2B, polar Al Saifi JC
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,58.08,508.1,4760000,0.232,0.304,3.2557,2.77409,253.406,0,0,0,0,0,0,0};//PCSAFT constants for acetone non associating
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,58.08,508.1,4760000,0.232,0.304,3.0848,3.0925,168.32,0.964,1321.2,2.88,0,1,1,0};//PCSAFT constants for acetone 2B
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_GV,58.08,508.1,4760000,0.232,0.304,3.2742,2.7447,232.99,0,0,2.88,0,0,0,0};//PCSAFT constants for acetone non associating, polar Gross-Vrabeck GV
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_GV,58.08,508.1,4760000,0.232,0.304,3.3659,2.5060,244.18,0,0,2.88,0,0,0,0};//PCSAFT constants for acetone non associating, polar Sadowski GV
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_GV,58.08,508.1,4760000,0.232,0.304,3.228,2.7855,210.14,0.0,0.0,2.88,0,0,0,0};//PCSAFT constants for acetone non associating, polar De Villiers GV
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_JC,58.08,508.1,4760000,0.232,0.304,3.6028,2.1873,245.49,0.0,0.0,2.72,0.2969,0,0,0};//PCSAFT constants for acetone non associating, polar De Villiers JC
    //FF_SaftEOSdata data1S={0,FF_DPCSAFT_GV,72.11,535.6,4154325,0.249,0.329,3.7334,2.3870,259.87,0,0,3.3,0,0,0,0};//PCSAFT constants for 2-butanone non associating, polar De Villiers
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,114.23,568.7,2490000,0.256,0.4,3.8373,3.8176,242.78,0.0,0.0,0,0,0,0,0};//Octane
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,226.448,722.0,14.1e5,0,0.742,3.9552,6.6485,254.70,0,0,0,0,0,0,0};//Hexadecano
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,58.12,425.1,3796000,0.0,0,3.7024,2.3421,222.38,0.0,0.0,0,0,0,0,0};//EOS constants for BUTANE PCSAFT
    //FF_SaftEOSdata data1S={0,FF_PCSAFT,106.165,616.89,3534600,0.0,0,3.756,3.186,284,0.0,0.0,0,0,0,0,0};//EOS constants for m-xylene PCSAFT

    //struct SWEOSdata sample={id,eos,MW,Tc,Pc,Zc,w,tRef,rhoRef,{n},{t},{a},{b},{e},{g},{d},{c},7,9,0};//Sample
    FF_SWEOSdata data1w={0,FF_SW,18.01528,647.096,22060000,0.23,0,647.096,17873.71609,{0.82728408749586,-0.18602220416584e1,-0.11199009613744e1,0.15635753976056,0.87375844859025,-0.36674403715731,0.53987893432436e-1,0.10957690214499e1,
              0.53213037828563e-1,0.13050533930825e-1,-0.41079520434476,0.14637443344120,-0.55726838623719e-1,-0.112017741438e-1,-0.66062758068099e-2,0.46918522004538e-2},
    {0.5,1.25,1.875,0.125,1.5,1.0,0.75,1.5,0.625,2.625,5.0,4.0,4.5,3.0,4.0,6.0},{0},{0},{0},{0},{1,1,1,2,2,3,4,1,5,5,1,2,4,4,1,1},{0,0,0,0,0,0,0,1,1,1,2,2,2,3,5,5},7,9,0};//Water
    FF_SWEOSdata data1x={0,FF_SW,18,647.096,22060000,0.23,0,647.096,17873.71609,{3.46821920e-1,5.03423025e-1,-3.51059570e-1,5.07004866e-2,1.99939129e-4,-5.69888763e-1,
              -1.96198912e-1,-2.02509554,-1.09353609,7.25785202e-2,2.16072642e-1,-1.01542630e-1,7.46926106e-1,2.18830463e-3},{1.5,0.25,1.25,0.25,0.875,1.375,0,2.375,2,2.125,3.5,6.5,4.75,12.5},
              {0},{0},{0},{0},{1,1,1,3,7,2,1,1,2,5,1,1,4,2},{0,0,0,0,0,0,1,1,1,1,2,2,2,3},6,8};//water
    //printf("%f %f %d %d\n",data1w.n[15],data1w.t[15],data1w.d[15],data1w.c[15]);

    FF_SWEOSdata data1Y={106.165,616.89,3534600,0,616.89,2665,{1.2791017E-05,0.041063111,1.505996,-2.3095875,-0.46969,0.171031,-1.001728,-0.3945766,0.6970578,-0.3002876,
            -0.024311,0.815488,-0.330647,-0.123393,-0.54661},{1,0.91,0.231,0.772,1.205,0.323,2.7,3.11,0.768,4.1,0.818,2,2.9,3.83,0.5 },{1.0244,1.3788,0.9806,6.3563},
    {0.713,0.9169,0.6897,0.7245},{1.66,1.9354,1.0323,78},{1.1013,0.6515,0.4975,1.26},{8,4,1,1,2,3,1,3,2,2,7,1,1,3,3},
    {0,0,0,0,0,0,2,2,1,2,1},6,5,4};//m-Xylene

    //We need Cp0 correlation for some thermodynamic properties
    FF_Correlation cp={0,2,{33.7633590698,-0.00594595819712,2.2357540729e-05,-9.9620089955e-09,1.09748701013e-12},50,1500};//Cp0 parameters for water
    //FF_Correlation cp={0,3,7.97183,6.27078,2572.63,2.0501,1156.72,50,1500};//Cp0 parameters for water DIPPR 107
    //FF_Correlation cp={0,2,24.258,0.23351,0.0001279,-2.4369e-07,8.5524e-11,50,1500};//Cp0 parameters for butane
    //FF_Correlation cp={0,3,18.6383,57.4178,1792.73,38.6599,814.151,50,1500};//Cp0 parameters for butane DIPPR 107
    //FF_Correlation cp={0,2,-33.28856,0.635358,-0.0003678323,6.869453E-08,7.498678E-12,50,1500};//Cp0 parameters for m-xylene

    //And now some calculations for single component
    //**********************************************

    //First do some direct cubic eos calculations
    printf("Cubic EOS calculations for the uncommented substance\n");
    printf("----------------------------------------------------\n");
    thR.MW=data1C.MW;
    thR.T=372.15;
    thR.P=101325;
    FF_FixedParamCubic(&data1C,&paramC);//We find the fixed parameter of the EOS
    FF_ThetaDerivCubic(&thR.T,&data1C,&paramC);//And now the temperature dependent ones
    printf("Some cubic parameters. b:%f  Theta:%f\n\n",paramC.b,paramC.Theta);//here we have some of the results

    //Now we will obtain the volume from T and P and, with it, we will obtain some other properties using direct call to cubic eos functions
    option='b';//we want the liquid and the gas answer, if possible
    FF_VfromTPcubic(&thR.T,&thR.P,&paramC,&option,answerL,answerG,&state);//We find the V,Arr and Z
    printf("Direct cubic EOS calc of V from T,P. T:%f, P:%f, Gas vol.:%f, Gas Arr:%f, Gas Z:%f\n\n",thR.T,thR.P,answerG[0],answerG[1],answerG[2]);

    //Now we will do the same but with call to generic functions, that will call the specific one
    eosType=FF_CubicType;//we need to indicate the type of eos in order to interpret the void pointer passed
    option='s';//We ask for the stable phase
    FF_VfromTPeos(&eosType,&thR.T,&thR.P,&data1C,&option,answerL,answerG,&state);
    if (state=='G')  thR.V=answerG[0];//we copy the gas volume
    else if (state=='L') thR.V=answerL[0];//we copy the liquid volume
    printf("Generic EOS calc of V from T,P determining the stable phase. V:%f, state:%c\n",thR.V,state);//'L' will indicate liquid, 'G' gas

    FF_ArrZfromTVeos(&eosType,&thR.T,&thR.V,&data1C,&Arr,&Z);
    printf("Arr and Z calculation from the obtained volume, T:%f V:%f Arr:%f Z:%f\n\n",thR.T,thR.V,Arr,Z);

    //Some residual, ideal, and consolidated, properties obtained from T and V
    FF_ExtResidualThermoEOS(&eosType,&data1C,&thR);
    printf("Residual properties from T,V. T:%f, V:%f, P:%f H:%f\n\n",thR.T,thR.V,thR.P,thR.H);
    FF_IdealThermoEOS(&cp.form,cp.coef,&refT,&refP,&thR);
    printf("Ideal properties from T,V. P:%f H:%f\n\n",thR.P,thR.H);
    FF_ThermoEOS(&eosType,&data1C,&cp.form,cp.coef,&refT,&refP,&thR);
    printf("Some consolidated properties from T,V. T:%f, V:%f, P:%f, H:%f, S:%f, SS:%f\n\n",thR.T,thR.V,thR.P,thR.H,thR.S,thR.SS);

    //vapor pressure and boiling point at given T and P
    FF_VpEOS(&eosType,&thR.T,&data1C,&Vp);
    printf("Vapor pressure:%f \n",Vp);
    FF_TbEOS(&eosType,&thR.P,&data1C,&Tb);
    printf("Boiling temperature:%f \n\n",Tb);

    //Consolidated properties from P and H
    var='H';
    thR.H=thR.H+2000;
    FF_ThermoEOSfromPX(&eosType,&data1C,&cp.form,cp.coef,&refT,&refP,&var,&thR,&liqFraction);
    printf("Properties from P,H. P:%f, H:%f, T:%f, Liq.fract.:%f\n\n",thR.P,thR.H,thR.T,liqFraction);



    //Now with PCSAFT calculations using generic functions
    printf("Calculations with PCSAFT\n");
    printf("------------------------\n");
    eosType=FF_SAFTtype;
    thR.MW=data1S.MW;
    FF_VfromTPeos(&eosType,&thR.T,&thR.P,&data1S,&option,answerL,answerG,&state);
    if (state=='G'){
        thR.V=answerG[0];//we copy the gas volume
        Arr=answerG[1];
        Z=answerG[2];
    }
    else if (state=='L'){
        thR.V=answerL[0];//we copy the liquid volume
        Arr=answerL[1];
        Z=answerL[2];
    }
    printf("V calculation from T,P. State:%c, T:%f, V:%f, Arr:%f, Z:%f\n\n",state,thR.T,thR.V,Arr,Z);//'L' will indicate liquid, 'G' gas

    FF_ExtResidualThermoEOS(&eosType,&data1S,&thR);
    printf("Residual properties from T,V. T:%f, V:%f, P:%f H:%f\n\n",thR.T,thR.V,thR.P,thR.H);
    FF_IdealThermoEOS(&cp.form,cp.coef,&refT,&refP,&thR);
    printf("Ideal properties. P:%f H:%f\n\n",thR.P,thR.H);
    FF_ThermoEOS(&eosType,&data1S,&cp.form,cp.coef,&refT,&refP,&thR);
    printf("Consolidated thermo prop from T,V. T:%f, V:%f, P:%f, H:%f, S:%f, SS:%f\n\n",thR.T,thR.V,thR.P,thR.H,thR.S,thR.SS);

    FF_VpEOS(&eosType,&thR.T,&data1S,&Vp);
    printf("T:%f, Vapor pressure:%f \n\n",thR.T,Vp);
    FF_TbEOS(&eosType,&thR.P,&data1S,&Tb);
    printf("P:%f, Boiling temperature:%f \n\n",thR.P,Tb);

    var='H';
    thR.H=thR.H+2000;
    //FF_ThermoEOSfromPX(const enum FF_EosType *eosType,const void *data,const int *equation,const double coef[],double *refT,double *refP,char *var, FF_ThermoProperties *th,double *liqFraction)
    FF_ThermoEOSfromPX(&eosType,&data1S,&cp.form,cp.coef,&refT,&refP,&var,&thR,&liqFraction);
    printf("Properties from P,H. P:%f, H:%f, T:%f, Liq.fract.:%f\n\n",thR.P,thR.H,thR.T,liqFraction);


}


//Some tests with mixtures
void TestEosMix()
{
    int numSubs=2;
    enum FF_EosType eosType;
    enum FF_MixingRule rule;
    FF_ThermoProperties thR;
    char option;
    char state;
    double answerL[3],answerG[3];
    double refT=298.15;
    double refP=1.01325e5;
    double bT,bP;//boiling temperature, pressure
    double dT,dP;//dew temperature, pressure
    double substPhiL[2],substPhiG[2];
    double bTguess, dTguess,bPguess, dPguess;

    printf("n-Butane-n-Pentane at 390K and 11 bar using PR78 EOS. Kij=0\n");
    printf("-----------------------------------------------------------\n");
    //first the needed declarations
    eosType=FF_CubicType;
    rule=FF_VdW;
    thR.T=390;
    thR.P=11e5;
    option='s';
    bTguess=300;
    dTguess=400;
    FF_CubicEOSdata dataC[2]={{0,FF_PR78,58.12,425.1,3797000,0.0,0.199,0,-1,0,0,0,0},{0,FF_PR78,72.15,469.7,3370000,0.263,0.251,0,-1,0,0,0,0}};//n-Butane-n-Pentane
    double pintParamC[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};//No interaction parameters used
    FF_Correlation cpC[2]={{0,2,24.2584095001,0.233508005738,0.000127899998915,-2.43686088197e-07,8.55236576003e-11,0,0},
                        {0,2,33.7804489136,0.248501002789,0.000253480189713,-3.83800198733e-07,1.29766600243e-10,0,0}};//n-Butane-n-Pentane
    double xC[2]={0.3563,1-0.3563};
    double yC[2]={0.3563,1-0.3563};
    FF_MixVfromTPeos(&eosType,&rule,&thR.T,&thR.P,&numSubs,dataC,pintParamC,xC,&option,answerL,answerG,&state);
    printf("T:%f, P:%f, Vliquid:%f, Vgas:%f, Zliquid:%f, Zgas:%f, Phase:%c\n",thR.T,thR.P,answerL[0],answerG[0],answerL[2],answerG[2],state);
    if (state=='G') thR.V=answerG[0];
    else if (state=='L') thR.V=answerL[0];
    FF_BubbleT(&eosType,&rule,&thR.P,&numSubs,dataC,pintParamC,xC,&bTguess,&bT,yC,substPhiL,substPhiG);
    FF_DewT(&eosType,&rule,&thR.P,&numSubs,dataC,pintParamC,yC,&bTguess,&dT,xC,substPhiL,substPhiG);
    printf("Bubble T:%f, Dew T.:%f\n",bT,dT);
    FF_MixThermoEOS(&eosType,&rule,&numSubs,dataC,pintParamC,xC,cpC,&refT,&refP,&thR);
    printf("Cp:%f, Cv:%f, Sound.S:%f\n\n",thR.Cp,thR.Cv,thR.SS);

    printf("Acetone-Hexane(x1=0.255)at 318.15K,PR78 EOS, Kij=-0.0218+3.142e-4*T+10.155/T\n");
    printf("----------------------------------------------------------------------------\n");
    //first the needed declarations
    eosType=FF_CubicType;
    rule=FF_VdW;
    thR.T=318.15;
    option='s';
    double xC3[2]={0.255,1-0.255};
    double yC3[2]={0.255,1-0.255};
    bPguess=1e5;
    FF_CubicEOSdata dataC3[2]={{0,FF_PR78,58.08,508.1,4760000,0.232,0.304,0,0,0,0,0},{0,FF_PR78,86.1766,507.6,3025000,0.2638,0.299,0,0,0,0,0}};
    double pintParamC3[24]={0,0,0,0,0,0,-0.0218,3.142e-4,10.155,0,0,0,-0.0218,3.142e-4,10.155,0,0,0,0,0,0,0,0,0};//No interaction parameters used
    FF_BubbleP(&eosType,&rule,&thR.T,&numSubs,dataC3,pintParamC3,xC3,&bPguess,&bP,yC3,substPhiL,substPhiG);
    printf("Bubble P:%f, acetone fraction in gas:%f\n",bP,yC3[0]);
    printf("Results according to DDBST GmbH: P: 78687 Pa, y1:0.516\n\n");

    printf("Acetone-Methanol(x1=0.1146) at 288.95K using PR78 EOS, Kij=0\n");
    printf("------------------------------------------------------------\n");
    //first the needed declarations
    eosType=FF_CubicType;
    rule=FF_VdW;
    thR.T=288.95;
    option='s';
    double xC2[2]={0.1146,1-0.1146};
    double yC2[2]={0.1146,1-0.1146};
    bPguess=1e5;
    FF_CubicEOSdata dataC2[2]={{0,FF_PR78,58.08,508.1,4760000,0.232,0.304,0,0,0,0,0},{0,FF_PR78,32.05,512.5,80.8e5,0.224,0.556,0,0,0,0,0}};
    double pintParamC2[24]={0,0,0,0,0,0,0,0,00,0,0,0,0,0,0,0,0,0,0,0,0,0,0};//No interaction parameters used
    FF_BubbleP(&eosType,&rule,&thR.T,&numSubs,dataC2,pintParamC2,xC2,&bPguess,&bP,yC2,substPhiL,substPhiG);
    printf("Bubble P:%f, acetone fraction in gas:%f\n",bP,yC2[0]);
    printf("Results according to DDBST GmbH: P: 13332 Pa, y1:0.2368\n\n");


    printf("Acetone-Methanol(x1=0.1146) at 288.95K using associating PCSAFT EOS, Kij=0\n");
    printf("--------------------------------------------------------------------------\n");
    //first the needed declarations
    eosType=FF_SAFTtype;
    rule=FF_NoMixRul;
    thR.T=288.95;
    option='s';
    double xS[2]={0.1146,1-0.1146};
    double yS[2]={0.1146,1-0.1146};
    bPguess=1e5;
    FF_SaftEOSdata dataS[2]={{0,FF_PCSAFT,58.08,508.1,4760000,0.232,0.304,3.0848,3.0925,168.32,0.964,1321.2,2.88,0,1,1,0},
                    {0,FF_PCSAFT,32.05,512.5,80.8e5,0.224,0.556,2.6510,2.7920,186.60,0.14600,2090.20,0,0,1,1,0}};
    double pintParamS[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};//No interaction parameters used
    FF_BubbleP(&eosType,&rule,&thR.T,&numSubs,dataS,pintParamS,xS,&bPguess,&bP,yS,substPhiL,substPhiG);
    printf("Bubble P:%f, acetone fraction in gas:%f\n",bP,yS[0]);
    printf("Results according to DDBST GmbH: P: 13332 Pa, y1:0.2368\n\n");
}



//Allows to test the correlation parameters optimization using water data as input. Takes 20 seconds
void TestCorrelationOptimization()
{
    unsigned nVar=5;
    FF_CorrelationData corrData;
//{enum FF_CorrEquation eq ;double Tc,Pc,rhoC;unsigned nPoints;double x[40],y[40];}FF_CorrelationData
    corrData.eq=FF_DIPPR101Vp;//vapor pressure of water
    corrData.nPoints=7;
    corrData.x[0]=273.15;
    corrData.y[0]=6.082e2;
    corrData.x[1]=373.15;
    corrData.y[1]=1.014e5;
    corrData.x[2]=473.15;
    corrData.y[2]=1.555e6;
    corrData.x[3]=573.15;
    corrData.y[3]=8.605e6;
    corrData.x[4]=423.15;
    corrData.y[4]=4.763e5;
    corrData.x[5]=523.15;
    corrData.y[5]=3.98e6;
    corrData.x[6]=623.15;
    corrData.y[6]=1.654e7;

    double var[nVar],lb[nVar],ub[nVar],grad[nVar];
    char enforceLimits[nVar];
    int i;
    for (i=0;i<nVar;i++) enforceLimits[i]='n';
    double error;
    printf("Coefficients optimization for water vapor pressure using DIPPR 101 equation\n");
    printf("The calculation takes 20 seconds\n");
    printf("---------------------------------------------------------------------------\n");
    FF_OptCorrelation(nVar,lb,ub,enforceLimits,&corrData,var,&error);
    printf("Error:%f %%.\n\n",error*100);
}


//Allows to test the EOS parameters optimization using water data as input
void TestEOSoptimization()
{
    //enum OptModel optModel=cubicParam;
    enum FF_OptModel optModel=SAFTParam;
    //unsigned nVar=3;
    unsigned nVar=5;
    unsigned optTime=30;//Modify time as necessary
    FF_EOSPvRhoData EOSData;

    //EOSData.eos=FF_PRFIT4;//Data for water. Use the EOS needed. Most cubics oiptimize only vapor pressure
    //EOSData.eos=FF_PRMC;
    EOSData.eos=FF_PCSAFT2B;
    //EOSData.eos=FF_PCSAFT4C;

    EOSData.zcFilter=0.07;//This allows to reject solutions that give Zc error > than the fraction indicated
    EOSData.ldensFilter=0.04;//the same for liquid density error

    EOSData.MW=18;
    EOSData.Tc=647;
    EOSData.Pc=220e5;
    EOSData.Zc=0.23;
    EOSData.w=0.344;
    EOSData.nPoints=8;
    EOSData.points[0][0]=273.15;
    EOSData.points[0][1]=6.082e2;
    EOSData.points[0][2]=998.6;
    EOSData.points[1][0]=373.15;
    EOSData.points[1][1]=1.014e5;
    EOSData.points[1][2]=957.3;
    EOSData.points[2][0]=473.15;
    EOSData.points[2][1]=1.555e6;
    EOSData.points[2][2]=863.4;
    EOSData.points[3][0]=573.15;
    EOSData.points[3][1]=8.605e6;
    EOSData.points[3][2]=695;
    EOSData.points[4][0]=423.15;
    EOSData.points[4][1]=4.763e5;
    EOSData.points[4][2]=915.9;
    EOSData.points[5][0]=523.15;
    EOSData.points[5][1]=3.98e6;
    EOSData.points[5][2]=793.9;
    EOSData.points[6][0]=323.15;
    EOSData.points[6][1]=1.235e4;
    EOSData.points[6][2]=987.2;
    EOSData.points[7][0]=623.15;
    EOSData.points[7][1]=1.654e7;
    EOSData.points[7][2]=5.479e+02;

    double var[nVar],lb[nVar],ub[nVar],grad[nVar];
    char enforceLimits[nVar];
    double error;
    int i;
    for (i=0;i<nVar;i++)
    {
        lb[i]=-1e10;
        ub[i]=1e10;
        var[i]=1;//(lb[i]+ub[i])/2;
        enforceLimits[i]='n';
    }
    printf("Coefficients optimization for water using different EOS. Select in the function\n");
    printf("Meaning of the coefficients depends on the selected EOS. It takes 30 seconds\n");
    printf("-------------------------------------------------------------------------------\n");
    FF_OptEOSparam(optTime,nVar,lb,ub,enforceLimits,&EOSData,var,&error);
    printf("Errors total:%f%%, Vp:%f%%, Liq.dens.:%f%%\n",error*100,EOSData.vpError*100,EOSData.ldensError*100);
}

//Allows to test Unifac calculation using Unifac standard and Unifac Dortmund. The model can be
TestUnifac(){
    int i;
    FF_UnifacData *uni =(FF_UnifacData*) calloc(1,sizeof(FF_UnifacData));
    uni->model=FF_UNIFACStd;
    int numData=9;
    int data1[9][3]={{0,1,2},{0,2,4},{1,1,1},{1,2,1},{1,14,1},{2,1,1},{2,2,4},{2,3,1},{3,9,6}};//hexane,ethanol,methylcyclopentane,benzene Unifac Std subgroups
    FF_UNIFACParams(numData, data1, uni);
    double T=334.85;//Temperature
    double x1[4]={0.162,0.068,0.656,0.114};//composition
    double lnGammaC[uni->numSubs],lnGammaR[uni->numSubs],gE;
    FF_ActivityUNIFAC(uni,&T,x1,lnGammaC,lnGammaR,&gE);
    printf("\nSystem hexane/ethanol/metylcylopentane/benzene(0.162,0.068,0.656,0.114) at 334.85 K using UNIFAC standard:\n");
       printf("------------------------------------------------------------------------------\n");
    for (i=0;i<uni->numSubs;i++) printf("Substance:%i Activity:%f\n",i, exp(lnGammaC[i]+lnGammaR[i]));
    printf("\n");

    uni->model=FF_UNIFACDort;
    numData=4;
    int data2[4][3]={{0,1,1},{0,18,1},{1,1,2},{1,2,4}};//acetone-hexane at 0.1/0.9 and 313ºK should give activity=3.6/1.0
    FF_UNIFACParams(numData, data2, uni);
    T=313;
    double x2[2]={0.1,0.9};
    FF_ActivityUNIFAC(uni,&T,x2,lnGammaC,lnGammaR,&gE);
    printf("\nSystem acetone/hexane(0.1,0.9) at 313 K using UNIFAC Dortmund:\n");
    printf("----------------------------------------------------------------\n");
    for (i=0;i<uni->numSubs;i++) printf("Substance:%i Activity:%f\n",i, exp(lnGammaC[i]+lnGammaR[i]));
    printf("\n");

    //Activity tests in polymers
    uni->model=FF_UNIFACZM;
    numData=4;
    //numData=6;
    int data3[4][3]={{0,9,6},{1,2,1*557},{1,3,1*557},{1,21,1*557}};//benzene-PVAc MW=48000 (557 monomer units)
    //int data3[6][3]={{0,14,1},{0,2,2},{0,1,1},{1,2,1*1974},{1,3,1*1974},{1,21,1*1974}};//propanol-PVAc MW=170000 (1974 monomer units)
    FF_UNIFACParams(numData, data3, uni);
    T=303.15;
    double MW[2]={78.114,86.1*557}, q[2]={0.24,1-0.24},massFract[2];//mass fractionn of benzene and PVAc
    //double MW[2]={60.096,86.1*1974}, q[2]={0.018,1-0.018},massFract[2];//mass fraction propanol and PVAc
    FF_FractionsCalculation(2, MW, q, 1, massFract, x2);
    FF_ActivityUNIFAC(uni,&T,x2,lnGammaC,lnGammaR,&gE);
    printf("\nSystem benzene/PVAc(0.24,0.76)in mass at 303 K using UNIFACZM Experimental=0.72\n");
    printf("-------------------------------------------------------------------------------\n");
    for (i=0;i<uni->numSubs;i++) printf("Substance:%i Activity:%f\n",i, exp(lnGammaC[i]+lnGammaR[i]));
    printf("\n");
    uni->FV[0]=MW[0]/(867.9*1000)-4.84e-5*1.2;//real molar volume less Van der Waals volume multiply by 1.2
    uni->FV[1]=MW[1]/(1185.6*1000)-4.79e-5*557*1.2;
    //uni->FV[0]=MW[0]/(804.6*1000)-4.217e-5*1.2;//real molar volume less Van der Waals volume multiply by 1.2
    //uni->FV[1]=MW[1]/(1185.6*1000)-4.79e-5*1974*1.2;
    uni->model=FF_EntropicFV;
    FF_ActivityUNIFAC(uni,&T,x2,lnGammaC,lnGammaR,&gE);
    printf("\nSystem benzene/PVAc(0.24,0.76)in mass at 303 K using EntropicFV 1.2\n");
    printf("-------------------------------------------------------------------\n");
    for (i=0;i<uni->numSubs;i++) printf("Substance:%i Activity:%f\n",i, exp(lnGammaC[i]+lnGammaR[i]));
    printf("\n");

    if (uni==NULL) return;
    free(uni);
}

//Load data for the system acetone-water
void loadAcetoneWater( FF_BaseProp baseProp[], FF_CubicEOSdata dataC[], FF_Correlation vpCorr[],double pintParamWilson[2][2][6],
                        double pintParamNRTL[2][2][6],double pintParamUNIQUAC[2][2][6]){
    //Acetone-Water
    baseProp[0].Tc=508.1;
    baseProp[0].Pc=4700000;
    baseProp[0].Zra=0.2477;
    baseProp[0].r=2.5734;
    baseProp[0].q=2.3359;
    baseProp[0].qRes=baseProp[0].q;
    baseProp[1].Tc=647.096;
    baseProp[1].Pc=2.2064E+07;
    baseProp[1].Zra=0.2372;
    baseProp[1].r=0.92;
    baseProp[1].q=1.4;
    baseProp[1].qRes=1;
    dataC[0].eos=FF_PR78;
    dataC[0].Tc=508.1;
    dataC[0].Pc=4700000;
    dataC[0].Zc=0.232;
    dataC[0].w=0.304;
    dataC[0].c=0.00;
    dataC[1].eos=FF_PR78;
    dataC[1].Tc=647.096;
    dataC[1].Pc=2.2064E+07;
    dataC[1].Zc=0.2294;
    dataC[1].w=0.344;
    dataC[1].c=0.00;
    vpCorr[0].form=21;
    vpCorr[0].coef[0]=6.696925E+01;
    vpCorr[0].coef[1]=-5.784592E+03;
    vpCorr[0].coef[2]=-7.858812E+00;
    vpCorr[0].coef[3]=7.104956E-06;
    vpCorr[0].coef[4]=2.0;
    vpCorr[1].form=21;
    vpCorr[1].coef[0]=6.702455E+01;
    vpCorr[1].coef[1]=-7.276391E+03;
    vpCorr[1].coef[2]=-7.342973E+00;
    vpCorr[1].coef[3]=4.161914E-06;
    vpCorr[1].coef[4]=2.0;
    pintParamWilson[0][1][0]=393.27;//lamda12-lambda11 cal/mol PychemQt. Usar Pol1C
    pintParamWilson[1][0][0]=1430;//PychemQt
    //pintParamWilson[0][1][0]=469.02;//Bruin 1970
    //pintParamWilson[1][0][0]=1489.07;//Bruin 1970

    /*pintParamNRTL[0][1][0]=377.577;//PychemQt. Usar Pol1K
    pintParamNRTL[0][1][4]=0.5856;
    pintParamNRTL[1][0][0]=653.885;
    pintParamNRTL[1][0][4]=0.5856;*/

    pintParamNRTL[0][1][0]=3489.3;//Tochigi et al. J/mol Pol2J
    pintParamNRTL[0][1][1]=3477.2e-3;
    pintParamNRTL[0][1][2]=-1582.2e-5;
    pintParamNRTL[0][1][4]=0.5466;

    pintParamNRTL[1][0][0]=-13165.1;
    pintParamNRTL[1][0][1]=73040.1e-3;
    pintParamNRTL[1][0][2]=-5699.4e-5;
    pintParamNRTL[1][0][4]=0.5466;

    pintParamUNIQUAC[0][1][0]=601.61;//PychemQt. Pol1C
    pintParamUNIQUAC[1][0][0]=-52.302;//PychemQt. Pol1C
}

//Load data for the system water-ethanol
void loadWaterEthanol(double Tc[2],double Pc[2],double Zra[2],double r[2],double q[2],double qRes[2],double pintParamWilson[2][2][6],
                        double pintParamNRTL[2][2][6],double pintParamUNIQUAC[2][2][6]){
    //water-ethanol
    Tc[0]=647.096;
    Pc[0]=2.2064E+07;
    Zra[0]=0.2372;
    r[0]=0.92;
    q[0]=1.4;
    qRes[0]=1;

    Tc[1]=516.2;
    Pc[1]=6383475;
    Zra[1]=0.2502;
    r[1]=2.1054;
    q[1]=1.972;
    qRes[1]=0.92;

    //pintParamWilson[0][1][0]=58.68/18.67;
    pintParamWilson[0][1][0]=975.49;//lamda12-lambda11 cal/mol PychemQt. Usar PolC
    pintParamWilson[1][0][0]=276.76;//lamda12-lambda11 cal/mol PychemQt

    pintParamNRTL[0][1][0]=670.441;//PychemQt. Usar PolK
    pintParamNRTL[0][1][4]=0.3031;
    pintParamNRTL[1][0][0]=-55.1681;
    pintParamNRTL[1][0][4]=0.3031;

    pintParamUNIQUAC[0][1][0]=232.01;//PychemQt. PolC
    pintParamUNIQUAC[1][0][0]=50.88;//PychemQt. PolC

    //pintParamUNIQUAC[0][1][0]=164.24;//Nagata 1989. PolK
    //pintParamUNIQUAC[1][0][0]=49.2;//PNagata 1989. PolK

    pintParamUNIQUAC[0][1][0]=-1881.18;//Voutsa 2011 Pol2K
    pintParamUNIQUAC[0][1][1]=9.2619;
    pintParamUNIQUAC[0][1][2]=-0.00938;
    pintParamUNIQUAC[1][0][0]=1030.38;
    pintParamUNIQUAC[1][0][1]=-4.4853;
    pintParamUNIQUAC[1][0][2]=0.004108;

    //pintParamUNIQUAC[0][1][0]=378.1;//en cal/mol Prausnitz
    //pintParamUNIQUAC[1][0][0]=258.4;//en cal/mol Prausnitz
}



/*
double f1(double x){
    return 3*pow(x,5)-4*pow(x,2)+2*x+1;//-1,2
}

double d1(double x){
    return 15*pow(x,4)-8*x+2;
}

double f2(double x){
    return cos(x)-x*x*x;//0,1
}

double d2(double x){
    return -sin(x)-3*x*x;
}

double f3(double x){
    return x*x-2;//0,2
}

double d3(double x){
    return 2*x;
}

int main(void){
    double x[3]={4.0,8.0};
    double y[3]={0,0};
    int n=2;
    FFnormalize(n,x);
    printf("x0:%f x1:%f\n",x[0],x[1]);
    double x2,y2,ftol=0.00000001;
    int niter=200;
    //FFsolverBisection(f2,&x2,&y2,&n,-5.0, 5.0,&ftol, &niter);
    //printf("%f %f %i\n",x2,y2,n);
    //FFsolverRegula(f2,&x2,&y2,&n,-5, 5.0,&ftol, &niter);
    //printf("%f %f %i\n",x2,y2,n);
    FFsolverRegulaI(f3,&x2,&y2,&n,-5, 5.0,&ftol, &niter);
    printf("%f %f %i\n",x2,y2,n);
    FFsolverRegulaA(f3,&x2,&y2,&n,-5, 5.0,&ftol, &niter);
    printf("%f %f %i\n",x2,y2,n);
    //FFsolverNewtonND(f2,&x2,&y2,&n,1,&ftol, &niter);
    //printf("%f %f %i\n",x2,y2,n);
    /*FFsolverNewton(f1,d3,&x2,&y2,&n,1,&ftol, &niter);
    printf("%f %f %i\n",x2,y2,n);
    FFsolverIntervalNewtonND(f3,x,y,&n,-5.0, 5.0,&ftol, &niter);
    printf("%f %f %f %f %f %f %i\n",x[0],y[0],x[1],y[1],x[2],y[2],n);
    FFsolverIntervalNewtonND2(f2,x,y,&n,-1.0, 0.0,&ftol, &niter);
    printf("%f %f %f %f %f %f %i\n",x[0],y[0],x[1],y[1],x[2],y[2],n);
    return 0;
}
*/

//Activity coeff test
void FF_ActivityTest(){
    int numSubs=2;
    double c[numSubs];
    double T;
    double P;
    FF_BaseProp baseProp[numSubs];
    FF_CubicEOSdata dataC[numSubs];
    double Tc[numSubs],Pc[numSubs],Zra[numSubs];
    double r[numSubs],q[numSubs],qRes[numSubs];
    double pintParamWilson[numSubs][numSubs][6];
    double pintParamNRTL[numSubs][numSubs][6];
    double pintParamUNIQUAC[numSubs][numSubs][6];
    double gamma[numSubs];
    double gE;
    double phi[numSubs];
    int i,j,k;
    enum FF_IntParamForm form;
    enum FF_ActModel model=FF_NRTL;
    bool useVp;
    enum FF_EOS eos[numSubs];
    eos[0]=FF_PR78;
    eos[1]=FF_PR78;
     FF_Correlation vpCorr[numSubs];
    for (i=0;i<numSubs;i++) for(j=0;j<numSubs;j++) for(k=0;k<6;k++) pintParamWilson[i][j][k]=0;
    //for(i=0;i<numSubs*numSubs*6;i++) pintParamNRTL[i]=0;
    for (i=0;i<numSubs;i++) for(j=0;j<numSubs;j++) for(k=0;k<6;k++) pintParamNRTL[i][j][k]=0;
    for (i=0;i<numSubs;i++) for(j=0;j<numSubs;j++) for(k=0;k<6;k++) pintParamUNIQUAC[i][j][k]=0;

    /*
   //water-PG

    pintParamNRTL[0][1][0]=-476.112;//PychemQt Pol1K
    pintParamNRTL[0][1][4]=0.3273;

    pintParamNRTL[1][0][0]=618.458;
    pintParamNRTL[1][0][4]=0.3273;

    pintParamNRTL[0][1][0]=2515.3;//C.Trujillo 2011 Pol1K
    pintParamNRTL[0][1][1]=-0.038797;
    pintParamNRTL[0][1][4]=0.3;

    pintParamNRTL[1][0][0]=-406.84;
    pintParamNRTL[1][0][1]=-0.02265;
    pintParamNRTL[1][0][4]=0.3;
    */

    loadAcetoneWater(baseProp,dataC,vpCorr,pintParamWilson,pintParamNRTL,pintParamUNIQUAC);
    //loadWaterEthanol(Tc,Pc,Zra,r,q,qRes,pintParamWilson,pintParamNRTL,pintParamUNIQUAC);

    c[0]=0.2;
    c[1]=1-c[0];
    T=273.15+150;
    P=9.86e5;

    form=FF_Pol2J;
    useVp=false;
    FF_PhiFromActivity(&numSubs,&model,baseProp,pintParamNRTL,&form,&useVp,eos,dataC,&T,&P,c,gamma,phi,&gE);
    printf("gamma: %f %f phi: %f %f\n",gamma[0],gamma[1],phi[0],phi[1]);
    useVp=true;
    FF_PhiFromActivity(&numSubs,&model,baseProp,pintParamNRTL,&form,&useVp,eos,vpCorr,&T,&P,c,gamma,phi,&gE);
    printf("gamma: %f %f phi: %f %f\n",gamma[0],gamma[1],phi[0],phi[1]);

    form=FF_Pol1C;
    FF_ActivityWilson(&numSubs,baseProp,pintParamWilson,&form,&T,c,gamma,&gE);
    printf("%f %f %f\n",gamma[0],gamma[1],gE);

    form=FF_Pol2J;
    FF_ActivityNRTL(&numSubs,pintParamNRTL,&form,&T,c,gamma,&gE);
    printf("%f %f %f\n",gamma[0],gamma[1],gE);

    form=FF_Pol2K;
    FF_ActivityUNIQUAC(&numSubs,baseProp,pintParamUNIQUAC,&form,&T,c,gamma,&gE);
    printf("%f %f %f\n",gamma[0],gamma[1],gE);
    return 0;
}
