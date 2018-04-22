/*
 * FFphysprop.c
 *
 *  Created on: 26/12/2015
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

// contains mainly calculations for pure substances physical properties, by non-eos methods
//Single substance, physical properties correlations
//==================================================


#include <math.h>
#include <stdio.h>
#include "FFbasic.h"
#include "FFphysprop.h"

//Calculates the result of the given equation
void CALLCONV FF_CorrelationResult(const int *eq,const double coef[],const int *nPoints,double x[],double y[]){//x contains the in variable, and y the out variable
    int i,j;
    double Tm;
    switch (*eq)
    {
    case FF_DIPPR100://DIPPR-100. Polynomial a+b*T+c*T^2+d*T^3+e*T^4. Liquid density,heat capacity,thermal conductivity. Vapor viscosity. Ideal gas heat capacity.Dielectric constant
    case FF_DIPPR100Ld:
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]*x[i]+coef[2]*pow(x[i],2)+coef[3]*pow(x[i],3)+coef[4]*pow(x[i],4);
        break;
    case FF_Polynomial://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5.Ideal gas heat capacity
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]*x[i]+coef[2]*pow(x[i],2)+coef[3]*pow(x[i],3)+coef[4]*pow(x[i],4)+coef[5]*pow(x[i],2);
        //printf("%f\n",y[0]);
        break;
    case FF_Polynomial2://Polynomial a+b*T^0.25+c*T^0.5+d*T+e*T^2+f*T^3.
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]*pow(x[i],0.125)+coef[2]*pow(x[i],0.25)+coef[3]*pow(x[i],0.5)+coef[4]*x[i];
        //printf("%f\n",y[0]);
        break;
	case FF_expDIPPR100://exp(a + b*T + c*T^2 + d*T^3 + e*T^4). Used for liquid viscosity
		for (i=0;i<*nPoints;i++) y[i]=exp(coef[0]+coef[1]*x[i]+coef[2]*pow(x[i],2)+coef[3]*pow(x[i],3)+coef[4]*pow(x[i],4));
        break;
    case FF_DIPPR101://DIPPR-101. exp{a+b/T+c*ln(T)+d*T^e} Liquid viscosity. Vapor pressure.
    case FF_DIPPR101Lv:
    case FF_DIPPR101Vp:
        for (i=0;i<*nPoints;i++) y[i]=exp(coef[0]+coef[1]/x[i]+coef[2]*log(x[i])+coef[3]*pow(x[i],coef[4]));
          //printf("y[i]:%f\n",y[i]);
        break;
    case FF_logDIPPR101:
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]/x[i]+coef[2]*log(x[i])+coef[3]*pow(x[i],coef[4]);
        break;
    case FF_DIPPR102://DIPPR-102. a*T^b/(1+c/T+d/T^2). Vapor viscosity. solid heat capacity
        for (i=0;i<*nPoints;i++) y[i]=coef[0]*pow(x[i],coef[1])/(1+coef[2]/x[i]+coef[3]/pow(x[i],2));
        break;
    case FF_DIPPR103://DIPPR-103. a + b*exp(-c/T^d)
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]*exp(-coef[2]/pow(x[i],coef[3]));
        break;
    case FF_DIPPR104://DIPPR-104 Polynomial. a+b/T+c/T^3+d/T^8+e/T^9. Second virial coefficient
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]/x[i]+coef[2]/pow(x[i],3)+coef[3]/pow(x[i],8)+coef[4]/pow(x[i],9);
        break;
    case FF_DIPPR105://DIPPR-105. a/b^{1+(1-T/c)^d}. Liquid density
        for (i=0;i<*nPoints;i++) y[i]=coef[0]/pow(coef[1],1+pow(1-x[i]/coef[2],coef[3]));
        break;
    case FF_DIPPR106://DIPPR-106. a*(1-Tr)^(b+c*Tr+d*Tr^2+e*Tr^3). Surface tension. Heat of vaporization. Liquid density. Independent variable is always Tr
    case FF_DIPPR106Hv:
    case FF_DIPPR106Ld:
    case FF_DIPPR106SurfT:
        for (i=0;i<*nPoints;i++){
            Tm=x[i]/coef[5];
            y[i]=coef[0]*pow((1-Tm),(coef[1]+coef[2]*Tm+coef[3]*pow(Tm,2)+coef[4]*pow(Tm,3)));
        }
        break;
    case FF_DIPPR107://DIPPR-107. Aly & Lee. a+b*{(c/T)/sinh(c/T)}^2+d*{(e/T)/cosh(e/T)}^2. Ideal gas heat capacity
    case FF_DIPPR107Cp:
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]*pow(coef[2]/x[i]/sinh(coef[2]/x[i]),2)+coef[3]*pow(coef[4]/x[i]/cosh(coef[4]/x[i]),2);
        break;
    case FF_DIPPR115://DIPPR115. exp{a+b/T+c*ln(T)+d*T^2+e/T^2}. Vapor pressure
        for (i=0;i<*nPoints;i++) y[i]=exp(coef[0]+coef[1]/x[i]+coef[2]*log(x[i])+coef[3]*pow(x[i],2)+coef[4]*pow(x[i],-2));
        break;
    case FF_DIPPR116://DIPPR116. a+b*Tr^0.35+c*Tr^(2/3)+ d*Tr + e*Tr^(4/3) with a=Rho crit and Tr=T/f with f=Tc. Liquid density
    case FF_DIPPR116Ld:
    case FF_PPDS10:
        for (i=0;i<*nPoints;i++){
            Tm=1-x[i]/coef[5];
            y[i]=coef[0]+coef[1]*pow(Tm,0.35)+coef[2]*pow(Tm,0.666667)+coef[3]*Tm+coef[4]*pow(Tm,1.333333);
        }
        break;
    case FF_PPDS9://e*exp[a*((c-T)/(T-d))^1/3+b*((c-T)/(T-d))^4/3]
        for (i=0;i<*nPoints;i++){

            y[i]=coef[4]*exp(coef[0]*pow((coef[2]-x[i])/(x[i]-coef[3]),0.33333)+coef[1]*pow((coef[2]-x[i])/(x[i]-coef[3]),1.33333));
        }
        break;
    case FF_Wilhoit://Wilhoit equation Cp0 J/mol·K (8 coefficients)
        for (i=0;i<*nPoints;i++){
            //if (x[i]>coef[7]) Tm=(x[i]-coef[7])/(x[i]+coef[6]);
            //else Tm=0;
            Tm=(x[i]-coef[7])/(x[i]+coef[6]);
            y[i]=R*(coef[0] + coef[1]/pow(x[i],2)*exp(-coef[2]/x[i]) + coef[3]*pow(Tm,2) + (coef[4] - coef[5] /pow((x[i] -coef[7]),2))*pow(Tm,8));
        }
        break;
    case FF_Cooper://Cooper (11 coefficients used in IAPWS95 and CO2),with potential term is used in short fundamental equations with 11 coefficients also(lacks last exp terms)
        for (i=0;i<*nPoints;i++){
            y[i]=coef[0]+coef[1]*pow(x[i],coef[2]);
            for (j=3;j<13;j=j+2){
                if (coef[j]>0) y[i]=y[i]+coef[j]*pow((coef[j+1]/x[i]),2)*exp(coef[j+1]/x[i])/pow((exp(coef[j+1]/x[i])-1),2);
            }
            y[i]=y[i]*R;
        }
        break;
    case FF_Jaechske://Jaeschke and Schley equation (9 coefficients). Used by GERG2004
        for (i=0;i<*nPoints;i++){
            y[i]=1+coef[0];
            for (j=1;j<9;j=j+4) if (coef[j]>0) y[i]=y[i]+coef[j]*pow(coef[j+1]/x[i]/sinh(coef[j+1]/x[i]),2)+coef[j+2]*pow(coef[j+3]/x[i]/cosh(coef[j+3]/x[i]),2);
            y[i]=y[i]*R;
        }
        break;
    case FF_ChemSep16://ChemSep nº16. a+exp(b/T+c+d*T+e*T^2). Ideal gas heat capacity. Liquid heat capacity,thermal conductivity,surface tension
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+exp(coef[1]/x[i]+coef[2]+coef[3]*x[i]+coef[4]*pow(x[i],2));
        break;
    case FF_Antoine1://Antoine. 10^{a-b/(T+c)}. Vapor pressure
        for (i=0;i<*nPoints;i++) y[i]=pow(10,(coef[0]-coef[1]/(x[i]+coef[2])));
        break;
    case FF_Antoine2://Antoine. exp{a-b/(T+c)}. Vapor pressure
        for (i=0;i<*nPoints;i++) y[i]=exp(coef[0]-coef[1]/(x[i]+coef[2]));
        break;
    case FF_Wagner25://Modified Wagner equation. a*exp{(b*Tm+c*Tm^1.5+d*Tm^2.5+e*Tm^5)/(1-Tm)} with a=Pc, and Tm=1-T/f; f=Tc. Vapor pressure
        for (i=0;i<*nPoints;i++){
            Tm=1-x[i]/coef[5];
            y[i]=coef[0]*exp((coef[1]*Tm+coef[2]*pow(Tm,1.5)+coef[3]*pow(Tm,2.5)+coef[4]*pow(Tm,5))/(1-Tm));
        }
    case FF_Wagner36://Original Wagner equation. a*exp{(b*Tm+c*Tm^1.5+d*Tm^3+e*Tm^6)/(1-Tm)} with a=Pc, and Tm=1-T/f; f=Tc. Vapor pressure
        for (i=0;i<*nPoints;i++){
            Tm=1-x[i]/coef[5];
            y[i]=coef[0]*exp((coef[1]*Tm+coef[2]*pow(Tm,1.5)+coef[3]*pow(Tm,3)+coef[4]*pow(Tm,6))/(1-Tm));
        }
        break;
    case FF_PCWIN://PCWIN. a+b*(1-T/e)+c*ln(1-T/e)+d*(1-T/e)^3 with e=Tc, liquid density
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]*(1-x[i]/coef[4])+coef[2]*log(1-x[i]/coef[4])+coef[3]*pow((1-x[i]/coef[4]),3);
        break;
    case FF_Rackett://Rackett. a/b^{(1-T/c)^d}. Liquid density
        for (i=0;i<*nPoints;i++) y[i]=coef[0]/pow(coef[1],pow(1-x[i]/coef[2],coef[3]));
        break;
    case FF_ExtAndrade1://Extended Andrade 1. exp(a+b/T+c*T+d*T^2). Surface tension. Liquid viscosity
        for (i=0;i<*nPoints;i++) y[i]=exp(coef[0]+coef[1]/x[i]+coef[2]*x[i]+coef[3]*pow(x[i],2)+coef[4]*pow(x[i],3));
        break;
    case FF_ExtAndrade2://Extended Andrade 2. 10^(a+b/T+c*T+d*T^2). Surface tension
        for (i=0;i<*nPoints;i++) y[i]=pow(10,(coef[0]+coef[1]/x[i]+coef[2]*x[i]+coef[3]*pow(x[i],2)+coef[4]*pow(x[i],3)));
        break;
    case FF_WagnerGd://Original Wagner equation. a*exp{(b*Tm+c*Tm^1.5+d*Tm^3+e*Tm^6)/(1-Tm)} with a=Pc, and Tm=1-T/f; f=Tc. Vapor pressure
        for (i=0;i<*nPoints;i++){
            Tm=1-x[i]/coef[5];
            y[i]=coef[0]*exp(coef[1]*pow(Tm,0.4)+coef[2]*Tm+coef[3]*pow(Tm,2.1)+coef[4]*pow(Tm,5.6));
        }
        break;
    }
}

//Calculates physical property, with input and output in SI units(kgr, not moles), using the given correlation, that may or not be in SI units
EXP_IMP void CALLCONV FF_PhysPropCorr(const int *cor,const double coef[],const double *MW,const int *nPoints,double x[],double y[]){
    //printf("%i %f %f %f %f %f %f %f\n",*cor,*MW,coef[0],coef[1],coef[2],coef[3],coef[4],coef[5]);
    int i;
    int eq;
    //First we convert the input variable if necessary
    if (*cor==22) for (i=0;i<*nPoints;i++) x[i]=x[i]-273.15;
    //Second we chose the calculation equation to use
    switch (*cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
    case 2://DIPPR 100 Cp0 in J/mol·K
    case 6://DIPPR 100 Cp0 in J/kgr·K
    case 15://DIPPR 100 Liquid Cp in J/mol·K
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 17://DIPPR 100 Liquid Cp in J/kgr·K
    case 40://DIPPR 100 Liquid density in kgr/m3
    case 50://DIPPR 100 Liquid thermal conductivity W/(m·K)
    case 60://DIPPR 100 Liquid surface tension N/m
    case 63://DIPPR 100 Liquid surface tension dyna/cm
    case 70://DIPPR 100 Solid density in Kmol/m3
    case 80://DIPPR 100 Solid Cp in J/Kmol·K
    case 111://DIPPR100 Gas viscosity in Pa·s
    case 121://DIPPR100 Gas thermal conductivity in W/(m·K)
        eq=FF_DIPPR100;
        break;
    case 10://Polynomial Cp0 en cal/(mol·K)
        eq=FF_Polynomial;
        //printf("Equation:%i\n",eq);
        break;
    case 20://DIPPR101 Vp Pa
    case 21://DIPPR101 Vp KPa
    case 30://DIPPR101 Liquid viscosity Pa·s
        eq=FF_DIPPR101;
        break;
    case 81://DIPPR 102 Solid Cp in J/Kmol·K
    case 110://DIPPR 102 Gas viscosity in Pa·s
    case 120://DIPPR 102 Gas thermal conductivity in W/(m·K)
        eq=FF_DIPPR102;
        break;
    case 41://DIPPR105 Liquid density in mol/dm3
    case 42://DIPPR105 Liquid density in kgr/m3
        eq=FF_DIPPR105;
        break;
    case 45://DIPPR106 Liquid density in mol/dm3
    case 47://DIPPR106 Liquid density in kg/m3
    case 61://DIPPR 106 Liquid surface tension N/m
    case 90://DIPPR 106 HvSat J/Kmol
    case 91://DIPPR 106 HvSat J/kgr
        eq=FF_DIPPR106;
        break;
    case 3://DIPPR107 Cp in cal/mol·K
    case 4://DIPPR107 Cp in J/Kmol·K
        eq=FF_DIPPR107;
        break;
    case 5://Wilhoit Cp0 J/mol·K
        eq=FF_Wilhoit;
        break;
    case 7://Cooper Cp0 J/mol·K
        eq=FF_Cooper;
        break;
    case 8://Jaechske Cp0 J/mol·K
        eq=FF_Jaechske;
        break;
    case 9://ChemSep nº16 Ideal gas heat capacity in J/Kmol·K
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
    case 51://ChemSep nº16 Liquid thermal conductivity W/(m·K)
    case 62://chemSep nº16 Liquid surface tension N/m
        eq=FF_ChemSep16;
        break;
    case 22: //Antoine base 10 in C and mmHg
    case 36: //Antoine base 10 in K and cP
        eq=FF_Antoine1;
        break;
    case 23://Antoine base e in K and Pa
        eq=FF_Antoine2;
        break;
    case 24://Wagner 25 in K and Pa
        eq=FF_Wagner25;
        break;
    case 25://Wagner 36 in K and Pa
        eq=FF_Wagner36;
        break;
    case 31://Extended Andrade 1 in cP
    case 33://Extended Andrade 1 in cP/MW
        eq=FF_ExtAndrade1;
        break;
    case 32://Extended Andrade 2 in cP
        eq=FF_ExtAndrade2;
        break;
    case 34://Cheric viscosity in cP
        eq=FF_ChericVisc;
        break;
    case 37://PPDS9 in Pa·s
        eq=FF_PPDS9;
        break;
    case 43://PCWIN liquid density in kgr/m3
        eq=FF_PCWIN;
        break;
    case 48://Rackett liquid density in kgr/dm3
        eq=FF_Rackett;
        break;
    case 44://DIPPR116 liquid density in mol/dm3
    case 46://DIPPR116 liquid density in kg/m3
        eq=FF_DIPPR116;
        break;
    case 101:
        eq=FF_WagnerGd;
        break;
    }

    FF_CorrelationResult(&eq,coef,nPoints,x,y);//We make the calculations

    switch (*cor){//Now we convert to SI units
    case 1://DIPPR 100 Cp in KJ/kgr·K
    case 21://DIPPR101 Vp KPa
    case 48://Rackett liquid density in kgr/dm3
        for (i=0;i<*nPoints;i++) y[i]=y[i]*1e3;
        break;
    case 2://DIPPR 100 Cp in J/mol·K
    case 5://Wilhoit Cp0 J/mol·K
    case 7://Cooper Cp0 J/mol·K
    case 8://Jaechske Cp0 J/mol·K
    case 15://DIPPR 100 Liquid Cp in J/mol·K
        for (i=0;i<*nPoints;i++) y[i]=y[i]/ *MW*1e3;
        break;
    case 3://DIPPR107 Cp in cal/mol·K
    case 10://Polynomial Cp0 en cal/(mol·K)
        for (i=0;i<*nPoints;i++) y[i]=y[i]/ *MW*1e3*4.1868;
        break;
    case 4://DIPPR107 Cp in J/Kmol·K
    case 9://ChemSep nº16 Ideal gas heat capacity in J/Kmol·K
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
    case 80://DIPPR 100 Solid Cp in J/Kmol·K
    case 81://DIPPR 102 Solid Cp in J/Kmol·K
    case 90://DIPPR 106 HvSat J/Kmol
        for (i=0;i<*nPoints;i++) y[i]=y[i]/ *MW;
        break;
    case 41://DIPPR105 Liquid density in mol/dm3
    case 44://DIPPR116 liquid density in mol/dm3
    case 45://DIPPR106 liquid density in mol/dm3
    case 70://DIPPR 100 Solid density in Kmol/m3
        for (i=0;i<*nPoints;i++) y[i]=y[i]* *MW;
        break;
    case 22://Antoine base 10 in C and mmHg
        for (i=0;i<*nPoints;i++) y[i]=y[i]*133.32239;
        break;
    case 31://Extended Andrade 1 in cP
    case 32://Extended Andrade 2 in cP
    case 34://Cheric viscosity in cP
    case 36://Antoine1 viscosity in cP
    case 63://DIPPR 100 Liquid surface tension dyna/cm
        for (i=0;i<*nPoints;i++) y[i]=y[i]*1e-3;
        break;
    case 33://Extended Andrade 1 in cP/MW
        for (i=0;i<*nPoints;i++) y[i]=y[i]* *MW*1e-3;
    }
    //Last we reconvert the input variable if necessary
    if (*cor==22) for (i=0;i<*nPoints;i++) x[i]=x[i]+273.15;
}

//Vapor pressure related calculations
//-----------------------------------
//Acentric factor calculation by the definition equation
void CALLCONV FF_WfromDefinition(const double *Pc,const double *Vp,double *w){
    *w=-log(*Vp/ *Pc)/log(10)-1;
}

//Acentric factor calculation from one vapor pressure
void CALLCONV FF_WfromOneVap(const double *Tc,const double *Pc,const double *T,const double *Vp,double *w){
    double Tr,tau,f0,f1,f2;
    Tr=*T/ *Tc;
    tau=1-Tr;
    f0=(-5.97616*tau+1.29874*pow(tau,1.5)-0.60394*pow(tau,2.5)-1.06841*pow(tau,5))/Tr;
    f1=(-5.03365*tau+1.11505*pow(tau,1.5)-5.41217*pow(tau,2.5)-7.46628*pow(tau,5))/Tr;
    *w=(f0-log(*Pc/ *Vp))/f1;
}

//Vapor pressure using Ambrose-Walton equation. Needs only Tc,Pc and w
void CALLCONV FF_VpAmbroseWalton(const  FF_BaseProp *baseProp,const int *nPoints,const double T[],double Vp[]){
    int i;
    double Tr[*nPoints],tau[*nPoints],f0[*nPoints],f1[*nPoints],f2[*nPoints];
    for (i=0;i<*nPoints;i++){
        Tr[i]=T[i]/ baseProp->Tc;
        tau[i]=1-Tr[i];
        f0[i]=(-5.97616*tau[i]+1.29874*pow(tau[i],1.5)-0.60394*pow(tau[i],2.5)-1.06841*pow(tau[i],5))/Tr[i];
        f1[i]=(-5.03365*tau[i]+1.11505*pow(tau[i],1.5)-5.41217*pow(tau[i],2.5)-7.46628*pow(tau[i],5))/Tr[i];
        f2[i]=(-0.64771*tau[i]+2.41539*pow(tau[i],1.5)-4.26979*pow(tau[i],2.5)+3.25259*pow(tau[i],5))/Tr[i];
        Vp[i]=baseProp->Pc*exp(f0[i]+f1[i]* baseProp->w+f2[i]*pow(baseProp->w,2));
    }
}

//Vapor pressure using Riedel-Vetere equation. Needs only Tc,Pc and one boiling point
void CALLCONV FF_VpRiedelVetere(const double *Tc,const double *Pc,const double *Tref,const double *VpRef,const int *type,
                               const int *nPoints,const double T[],double Vp[]){
    double Tr[*nPoints],TrefR,PrefR,h,K,psi,alpha,Q;
    int i;
    TrefR=*Tref/ *Tc;
    PrefR=*VpRef/ *Pc;
    //first we proceed to calculate h and K
    if(*VpRef>0){
        //first we proceed to calculate h and K
        h=-TrefR*log(PrefR)/(1-TrefR);
        if(*type==1)K=0.373-0.03*h;
        else if(*type==2)K=-0.12+0.025*h;
        else K=0.0838;
        //Now psi,alpha
        psi=-35+36/TrefR+42*log(TrefR)-pow(TrefR,6);
        alpha=(3.758*K*psi-log(PrefR)) /(K*psi-log(TrefR));
        Q=K*(3.758-alpha);
        //printf("K,psi,alpha,Q: %f %f %f %f\n",K,psi,alpha,Q);
    }
    //And now Vp
    for (i=0;i<*nPoints;i++){
        Tr[i]=T[i]/ *Tc;
        Vp[i]=*Pc*exp(-35*Q+36*Q/Tr[i]+(42*Q+alpha)*log(Tr[i])-Q*pow(Tr[i],6));
    }
}

//Liquid density calculations
//---------------------------
//Rackett equation for saturated liquid density. If you supply ref rho, Vc is not used. It is better to use w than Zra, and Zra than Zc
void CALLCONV FF_LiqDensSatRackett(const  FF_BaseProp *baseProp,const double *Tref,const double *rhoRef,
                                   const int *nPoints,const double T[],double rho[]){
    double Vc,Zc;//We define local variables in order to use them in the calculation
    int i;
    if (baseProp->Vc==0) Vc=baseProp->Zc*R*baseProp->Tc/baseProp->Pc;
    else Vc=baseProp->Vc;
    if (baseProp->w>0) Zc=0.29056-0.08775*baseProp->w;
    else if (baseProp->Zra>0) Zc=baseProp->Zra;
    else Zc=baseProp->Zc;
    //printf("%f %f %f %f %f %f\n",*MW,*Tc,*Pc,*Zc,*Vc,*rhoRef);
    for(i=0;i<*nPoints;i++){
        if (*rhoRef==0) rho[i]= baseProp->MW/(Vc* pow(Zc,pow((1- T[i]/ baseProp->Tc),0.28571)))/1000;//Rackett equation with no reference data in kgr/m3
        else rho[i]= *rhoRef * pow(Zc,(pow((1- *Tref / baseProp->Tc),0.28571)-pow((1- T[i]/ baseProp->Tc),0.28571)));//Rackett equation with reference data
    }
}

//Chueh-Prausnitz pressure correction for liquid density
void CALLCONV FF_LiqDensChuehPrausnitz(const  FF_BaseProp *baseProp, const int *nPoints,const double T[],const double P[],
                                       const double Vp[],double rhoin[], double rho[]){
    int i;
    double Tr,N;
    for(i=0;i<*nPoints;i++){
        Tr=T[i]/baseProp->Tc;
        N=(1-0.89*baseProp->w)*exp(6.9547-76.2853*Tr+191.306*Tr*Tr-203.5472*pow(Tr,3)+82.763*pow(Tr,4));
        rho[i]=pow((1+9*baseProp->Zc*N*(P[i]-Vp[i])/baseProp->Pc),0.11111)*rhoin[i];
    }
}
