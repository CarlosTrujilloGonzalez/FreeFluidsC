/*
 * FFphysprop.c
 *
 *  Created on: 26/12/2015
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2019  Carlos Trujillo Gonzalez

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
#include "FFeosPure.h"

//Calculates the result of the given equation
void CALLCONV FF_CorrelationResult(const int *eq,const double coef[],const int *nPoints,double x[],double y[]){//x contains the in variable, and y the out variable
    int i,j;
    double Tm;
    switch (*eq)
    {
    case FF_DIPPR100://DIPPR-100. Polynomial a+b*T+c*T^2+d*T^3+e*T^4. Liquid density,heat capacity,thermal conductivity. Vapor viscosity. Ideal gas heat capacity.Dielectric constant
    case FF_DIPPR100Ld:
        //for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]*x[i]+coef[2]*pow(x[i],2)+coef[3]*pow(x[i],3)+coef[4]*pow(x[i],4);
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+x[i]*(coef[1]+x[i]*(coef[2]+x[i]*(coef[3]+x[i]*coef[4])));
        break;
    case FF_Polynomial://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5.Ideal gas heat capacity
        for (i=0;i<*nPoints;i++) y[i]=coef[0]+coef[1]*x[i]+coef[2]*pow(x[i],2)+coef[3]*pow(x[i],3)+coef[4]*pow(x[i],4)+coef[5]*pow(x[i],5);
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
    case FF_PPDS15://a/Tm+b+c*Tm+d*Tm^2+e*Tm^3
        for (i=0;i<*nPoints;i++){
            Tm=1-x[i]/coef[5];
            y[i]=coef[0]/Tm+coef[1]+coef[2]*Tm+coef[3]*Tm*Tm+coef[4]*Tm*Tm*Tm;
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
        break;
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
    case FF_Tait://Tait equation with P=1e5, V=(a+b*T+c*T^2)*(1-0.0894*ln(1+P/(d*exp(-e*T))). Liquid volume
        for (i=0;i<*nPoints;i++) y[i]=(coef[0]+coef[1]*x[i]+coef[2]*pow(x[i],2))*(1-0.0894*log(1+1e5/(coef[3]*exp(-coef[4]*x[i]))));
        break;
    case FF_ExtWagner://Wagner equation with 6 terms
        for (i=0;i<*nPoints;i++){
            Tm=1-x[i]/coef[1];
            y[i]=coef[0]*exp(coef[1]*(coef[2]*pow(Tm,coef[3])+coef[4]*pow(Tm,coef[5])+coef[6]*pow(Tm,coef[7])+coef[8]*pow(Tm,coef[9])+
                    coef[10]*pow(Tm,coef[11])+coef[12]*pow(Tm,coef[13]))/x[i]);
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
    if ((*cor==22)||(*cor==63)||(*cor==240)) for (i=0;i<*nPoints;i++) x[i]=x[i]-273.15;
    //Second we chose the calculation equation to use
    switch (*cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
    case 2://DIPPR 100 Cp0 in J/mol·K
    case 6://DIPPR 100 Cp0 in J/kgr·K
    case 15://DIPPR 100 Liquid Cp in J/mol·K
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 17://DIPPR 100 Liquid Cp in J/kgr·K
    case 40://DIPPR 100 Liquid density in kgr/m3
    case 49://DIPPR 100 Liquid density in cm3/mol
    case 50://DIPPR 100 Liquid thermal conductivity W/(m·K)
    case 60://DIPPR 100 Liquid surface tension N/m
    case 63://DIPPR 100 Liquid surface tension dyna/cm
    case 70://DIPPR 100 Solid density in Kmol/m3
    case 71://DIPPR 100 Solid density in kgr/m3
    case 80://DIPPR 100 Solid Cp in J/(Kmol·K)
    case 82://DIPPR 100 Solid Cp in J/(kgr·K)
    case 111://DIPPR100 Gas viscosity in Pa·s
    case 112://DIPPR100 viscosity correction factor for NIST ECS calculation. Input reduced rho, output factor
    case 121://DIPPR100 Gas thermal conductivity in W/(m·K)
    case 140://DIPPR100 T in K from liquid enthalpy in J/kgr
    case 150://DIPPR100 Isothermal compressibility adimensional
        eq=FF_DIPPR100;
        break;
    case 10://Polynomial Cp0 en cal/(mol·K)
        eq=FF_Polynomial;
        //printf("Equation:%i\n",eq);
        break;
    case 130://Polynomial2 boil temp in K
        eq=FF_Polynomial2;
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
    case 200://DIPPR107 Cp in J/kg·K
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
    case 26://Wagner 36 in K and Pa
        eq=FF_ExtWagner;
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
    case 19://PPDS15 Cp in J/kg·K
    case 151://PPDS15 liquid isothermal compressibility adimensional
        eq=FF_PPDS15;
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
    case 240://Tait equation for polymer density in m3/kg at 1 bar
        eq=FF_Tait;
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
        break;
    case 49://DIPPR 100 Liquid density in cm3/mol
        for (i=0;i<*nPoints;i++) y[i]=*MW*1e3/y[i];
        break;
    case 240://Tait equation for polymer density in m3/kg at 1 bar
        for (i=0;i<*nPoints;i++) y[i]=1/y[i];
    }
    //Last we reconvert the input variable if necessary
    if ((*cor==22)||(*cor==63)||(*cor==240)) for (i=0;i<*nPoints;i++) x[i]=x[i]+273.15;
}

//Calculates specific enthalpy and entropy from a Cp correlation, with reference T=0 K
void FF_SpecificEnthalpyEntropyCorr(int *cor,const double coef[],double *MW,const int *nPoints,double x[],double H[],double S[]){
    int i;
    switch (*cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
    case 2://DIPPR 100 Cp0in J/mol·K
    case 6://DIPPR 100 Cp0 in J/kgr·K
    case 15://DIPPR 100 Liquid Cp in J/mol·K
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 17://DIPPR 100 Liquid Cp in J/kgr·K
        for (i=0;i<*nPoints;i++){
            H[i]=x[i]*(coef[0]+x[i]*(coef[1]/2+x[i]*(coef[2]/3+x[i]*(coef[3]/4+x[i]*coef[4]/5))));//This is the integration from Cp0
            S[i]=coef[0]*log(x[i])+x[i]*(coef[1]+x[i]*(coef[2]/2+x[i]*(coef[3]/3+x[i]*coef[4]/4)));
        }
        break;
    case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
        for (i=0;i<*nPoints;i++){
            H[i]=x[i]*(coef[0]+x[i]*(coef[1]/2+x[i]*(coef[2]/3+x[i]*(coef[3]/4+x[i]*(coef[4]/5+x[i]*coef[5]/6)))));//This is the integration from Cp0
            S[i]=coef[0]*log(x[i])+x[i]*(coef[1]+x[i]*(coef[2]/2+x[i]*(coef[3]/3+x[i]*(coef[4]/4+x[i]*coef[5]/5))));
        }

        break;
    case 3://DIPPR 107 correlation in calories/mol·K
    case 4://DIPPR 107 correlation in J/Kmol·K
    case 200://DIPPR107 Cp in J/kg·K
        for (i=0;i<*nPoints;i++){
            H[i]=coef[0]*x[i]+coef[1]*coef[2]*(1/tanh(coef[2]/x[i]))-coef[3]*coef[4]*tanh(coef[4]/x[i]);
            S[i]=coef[0]*log(x[i])+coef[1]*(coef[2]/x[i]/tanh(coef[2]/x[i])-log(sinh(coef[2]/x[i])))-coef[3]*(coef[4]/x[i]*tanh(coef[4]/x[i])-log(cosh(coef[4]/x[i])));
        }
        break;
    /*
    case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
        int j=20;
        double interval;
        double T,Cp1,Cp2;
        for (i=0;i<*nPoints;i++){
            interval=x[i]/i;
            Cp1=exp(coef[1]/ *refT+coef[2]+coef[3]* *refT+coef[4]*pow(*refT,2));
            for (j=1;j<21;j++){
                T=*refT+j*interval;
                Cp2=exp(coef[1]/T+coef[2]+coef[3]*T+coef[4]*pow(T,2));
                th0->H=th0->H+(Cp1+Cp2)/2*interval;
                th0->S=th0->S+(Cp1+Cp2)/(T+T-interval)*interval;
                Cp1=Cp2;
                }
        }
        break;
    */
    case 5:{//Wilhoit equation J/mol·K (8 coefficients)
        int j;
        double y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
        a7_a6=coef[7]/coef[6];
        a7_a6_2=a7_a6*a7_a6;
        a7_a6_4=a7_a6_2*a7_a6_2;
        x1=(coef[4]*coef[7]*coef[7] - coef[5])/(coef[6]*coef[6]);
        for (i=0;i<*nPoints;i++){
            if (x[i]<=coef[7]) y=0;
            else y=(x[i]-coef[7])/(x[i]+coef[6]);
            y2=y*y;
            y4=y2*y2;
            if (x[i]<=coef[7]) h=0;
            else h=(coef[6]+coef[7])*((2*coef[3]+8*coef[4])*log(1-y)+ (coef[3]*(1+1/(1-y))+coef[4]*(7+1/(1-y)))*y+
                    coef[4]*(3*y2+5*y*y2/3+y4+0.6*y4*y+y4*y2/3)+ (coef[4]-coef[5]/pow((coef[6]+coef[7]),2))*y4*y2*y/7);
            H[i]= R*x[i]*(coef[0]+coef[1]*exp(-coef[2]/x[i])/(coef[2]*x[i]))+R*h;

            if (x[i]<=coef[7]) s=0;
            else{
                z = x[i]*(coef[7] + coef[6])/((x[i] + coef[6])*coef[7]);
                w=0;
                for (j=1;j<8;j++) w=w+(x1*pow(-a7_a6,6-j) - coef[4])*pow(y,j)/j;
                s=(coef[3] + ((coef[4]*coef[7]*coef[7]-coef[5])*a7_a6_4/(coef[6]*coef[6])))*a7_a6_2*log(z)+
                    (coef[3] + coef[4])*log((x[i] + coef[6])/(coef[6] + coef[7]))-
                    (coef[3]*(coef[6] + coef[7])/coef[6] + coef[5]*y4*y2/(7.*coef[7]*(coef[6] + coef[7])))*y+w;
            }
            S[i]=R*(coef[0]*log(x[i])+coef[1]*(1+coef[2]/x[i])*exp(-coef[2]/x[i])/(coef[2]*coef[2])+s);
        }
    }
        break;
    case 7:{//Cooper (11 coefficients used in IAPWS95 and CO2) plus potential term  (used in short fundamental equations with 11 coefficients also,lacks last exp terms)
        int j;
        for (i=0;i<*nPoints;i++){
            /*
            th0->Cp=coef[0]+coef[1]*pow(th0->T,coef[2]);
            for (i=3;i<13;i=i+2){
                if (coef[i]>0) th0->Cp=th0->Cp+coef[i]*pow((coef[i+1]/th0->T),2)*exp(coef[i+1]/th0->T)/pow((exp(coef[i+1]/th0->T)-1),2);
            }
            th0->Cp=th0->Cp*R;*/
            H[i]=coef[0]*x[i]+coef[1]/(coef[2]+1)*pow(x[i],(coef[2]+1));
            for (j=3;j<13;j=j+2) if (coef[j]>0) H[i]=H[i]+coef[j]*coef[j+1] / (exp(coef[j+1]/x[i])-1);
            H[i]=H[i]*R;
            if (coef[1]>0) S[i]=coef[0]*log(x[i])+coef[1]/coef[2]*pow(x[i],coef[2]);
            else S[i]=coef[0]*log(x[i]);
            //printf("T:%f S0(0):%f\n",x[i],S[i]*R/th0->MW);
            for (j=3;j<13;j=j+2){
                if (coef[j]>0) S[i]=S[i]+coef[j]*coef[j+1]*(exp(coef[j+1]/x[i])/(x[i]*(exp(coef[j+1]/x[i])-1))-
                    log(exp(coef[j+1]/x[i])-1)/coef[j+1]);
                //printf("S0(%i):%f\n",i,S[i]*R/th0->MW);
            }
            S[i]=S[i]*R;
        }
    }
        break;
    case 8:{//Jaeschke and Schley equation (9 coefficients). Used by GERG2004
        int j;
        for (i=0;i<*nPoints;i++){
            H[i]=(1+coef[0])*x[i];
            for (j=1;j<9;j=j+4) if (coef[j]>0) H[i]=H[i]+2*coef[j]*coef[j+1] / (exp(2*coef[j+1]/x[i])-1)+
                    +2*coef[j+2]*coef[j+3] / (exp(2*coef[j+3]/x[i])+1);
            H[i]=H[i]*R;
            S[i]=(1+coef[0])*log(x[i]);
            for (j=1;j<9;j=j+4) if (coef[j]>0) S[i]=S[i]+coef[j]*(coef[j+1]/x[i]/tanh(coef[j+1]/x[i])-log(sinh(coef[j+1]/x[i])))-coef[j+2]*(coef[j+3]/x[i]*tanh(coef[j+3]/x[i])-log(cosh(coef[j+3]/x[i])));
            S[i]=S[i]*R;
        }
    }
        break;
    case 19:{//PPDS15 equation a/Tm+b+c*Tm+d*Tm^2+e*Tm^3 Cpl in J/(kgr·K)
        double Tm,Tm2,Tc2;
        Tc2=coef[5]*coef[5];
        for (i=0;i<*nPoints;i++){
            Tm=1-x[i]/coef[5];
            Tm2=Tm*Tm;
            H[i]=-coef[5]*(coef[0]*log(Tm)+coef[1]*Tm+0.5*coef[2]*Tm2+coef[3]*Tm2*Tm/3+0.25*coef[4]*Tm2*Tm2);
            S[i]=-(6*Tc2*x[i]*(coef[2]+2*coef[3]+3*coef[4])-3*coef[5]*x[i]*x[i]*(3*coef[4]+coef[3])+2*x[i]*x[i]*x[i]*coef[4])/(6*Tc2*coef[5])+
                    (coef[0]+coef[1]+coef[2]+coef[3]+coef[4])*log(x[i])-coef[0]*log(abs(x[i]-coef[5]));
        }
      }
        break;

    }
    //now is necessary to pass to J/kg·K
    switch (*cor){
    case 1://DIPPR 100 Cp0 in KJ/kgr·K
        for (i=0;i<*nPoints;i++){
            H[i]=H[i]*1e3;
            S[i]=S[i]*1e3;
        }
        break;
    case 2://DIPPR 100 Cp0in J/mol·K
    case 5://Wilhoit equation J/mol·K (8 coefficients)
    case 7://Cooper J/mol·K (11 coefficients)
    case 8://Jaeschke and Schley equation J/mol·K (9 coefficients)
    case 15://DIPPR 100 Liquid Cp in J/mol·K
        for (i=0;i<*nPoints;i++){
            H[i]=H[i]*1e3/ *MW;
            S[i]=S[i]*1e3/ *MW;
        }
        break;
    case 4://DIPPR 107 correlation in J/Kmol·K
    case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
    case 16://DIPPR 100 Liquid Cp in J/Kmol·K
    case 18://ChemSep nº16 Liquid Cp in J/Kmol·K
        for (i=0;i<*nPoints;i++){
            H[i]=H[i]/ *MW;
            S[i]=S[i]/ *MW;
        }
        break;
    case 3://DIPPR 107 correlation in calories/mol·K
    case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 Cp0 in cal/(mol·K)
        for (i=0;i<*nPoints;i++){
            H[i]=H[i]*4.1868*1e3/ *MW;
            S[i]=S[i]*4.1868*1e3/ *MW;
        }
        break;
    }
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
void CALLCONV FF_VpAmbroseWalton(const FF_BaseProp *baseProp,const double *T,double *Vp){
    int i;
    double Tr,tau,f0,f1,f2;
    *Vp=0;
    if((baseProp->Tc>0)&&(baseProp->Pc>0)&&(baseProp->w>0)){
        if((*T>0)&&(*T<baseProp->Tc)){
            Tr=*T/ baseProp->Tc;
            tau=1-Tr;
            f0=(-5.97616*tau+1.29874*pow(tau,1.5)-0.60394*pow(tau,2.5)-1.06841*pow(tau,5))/Tr;
            f1=(-5.03365*tau+1.11505*pow(tau,1.5)-5.41217*pow(tau,2.5)-7.46628*pow(tau,5))/Tr;
            f2=(-0.64771*tau+2.41539*pow(tau,1.5)-4.26979*pow(tau,2.5)+3.25259*pow(tau,5))/Tr;
            *Vp=baseProp->Pc*exp(f0+f1* baseProp->w+f2*pow(baseProp->w,2));
        }
    }
}

//Vapor pressure using Riedel-Vetere equation. Needs only Tc,Pc and one boiling point
void CALLCONV FF_VpRiedelVetere(const  FF_BaseProp *baseProp,const double *Tref,const double *VpRef,const double *T,double *Vp){
    double Tr,TrefR,PrefR,h,K,psi,alpha,Q;
    int i;
    *Vp=0;
    if((baseProp->Tc>0)&&(baseProp->Pc>0)&&(*Tref>0)&&(*VpRef>0)){
        TrefR=*Tref/ baseProp->Tc;
        PrefR=*VpRef/ baseProp->Pc;
        //first we proceed to calculate h and K
        //first we proceed to calculate h and K
        h=-TrefR*log(PrefR)/(1-TrefR);
        if(baseProp->type==FF_Alcohol)K=0.373-0.03*h;
        else if(baseProp->type==FF_Acid)K=-0.12+0.025*h;
        else K=0.0838;
        //Now psi,alpha
        psi=-35+36/TrefR+42*log(TrefR)-pow(TrefR,6);
        alpha=(3.758*K*psi-log(PrefR)) /(K*psi-log(TrefR));
        Q=K*(3.758-alpha);
        //printf("K,psi,alpha,Q: %f %f %f %f\n",K,psi,alpha,Q);
        //And now Vp
        if((*T>0)&&(*T<baseProp->Tc)){
            Tr=*T/ baseProp->Tc;
            *Vp=baseProp->Pc*exp(-35*Q+36*Q/Tr+(42*Q+alpha)*log(Tr)-Q*pow(Tr,6));
        }
    }
}

//Vapor pressure calculation by correlations
void CALLCONV FF_Vp(double *T,FF_SubstanceData *data,double *Vp){
    int nPoints=1;
    *Vp=0;
    if(data->vpCorr.form>0) FF_PhysPropCorr(&data->vpCorr.form,data->vpCorr.coef,&data->baseProp.MW,&nPoints,T,Vp);
    else if((data->baseProp.Tc>0)&&(data->baseProp.Pc>0)){
        if((data->vp.x>0)&&(data->vp.y>0)) FF_VpRiedelVetere(&data->baseProp,&data->vp.x,&data->vp.y,T,Vp);
        else if(data->baseProp.w>0) FF_VpAmbroseWalton(&data->baseProp,T,Vp);
    }
}

//Density equations
//-----------------
//Rackett equation for saturated liquid density. If you supply ref rho, Vc is not used. It is better to use w than Zra, and Zra than Zc
void CALLCONV FF_LiqDensSatRackett(const  FF_BaseProp *baseProp,const double *Tref,const double *rhoRef,const double *T,double *rho){
    double Vc,Zc;//We define local variables in order to use them in the calculation
    int i;
    if (baseProp->Vc==0) Vc=baseProp->Zc*R*baseProp->Tc/baseProp->Pc;
    else Vc=baseProp->Vc;
    if (baseProp->w>0) Zc=0.29056-0.08775*baseProp->w;
    else if (baseProp->Zra>0) Zc=baseProp->Zra;
    else Zc=baseProp->Zc;
    //printf("%f %f %f %f %f %f\n",*MW,*Tc,*Pc,*Zc,*Vc,*rhoRef);       
        if ((*Tref>0)&&(*rhoRef>0)) *rho= *rhoRef * pow(Zc,(pow((1- *Tref / baseProp->Tc),0.28571)-pow((1- *T/ baseProp->Tc),0.28571)));//Rackett equation with reference data
        else *rho= baseProp->MW/(Vc* pow(Zc,pow((1- *T/ baseProp->Tc),0.28571)))/1000;//Rackett equation with no reference data in kgr/m3
}

//Chueh-Prausnitz pressure correction for liquid density
void CALLCONV FF_LiqDensChuehPrausnitz(const  FF_BaseProp *baseProp, const double *T,const double *P,const double *Pin,double *rhoIn, double *rho){
    int i;
    double Tr,N;
    Tr=*T/baseProp->Tc;
    N=(1-0.89*baseProp->w)*exp(6.9547-76.2853*Tr+191.306*Tr*Tr-203.5472*pow(Tr,3)+82.763*pow(Tr,4));
    *rho=pow((1+9*baseProp->Zc*N*(*P- *Pin)/baseProp->Pc),0.11111)* *rhoIn;
}

//Tait equation for polymer density, with T and P dependence
void CALLCONV FF_LiqDensTait(const int *eq,const double coef[],const  FF_BaseProp *baseProp, const int *nPoints,const double T[],const double P[], double rho[]){
    int i;
    double Tu;
    for(i=0;i<*nPoints;i++){
        if (*eq==240) Tu=T[i]-273.15;
        else Tu=T[i];
        rho[i]=(coef[0]+coef[1]*Tu+coef[2]*pow(Tu,2))*(1-0.0894*log(1+P[i]/(coef[3]*exp(-coef[4]*Tu))));
        if (*eq==240) rho[i]=1/rho[i];//As Tait equation gives volume
    }
}

//Liquid density calculation
void CALLCONV FF_LiqDensTP(double *T,double *P,FF_SubstanceData *data,double *liqDens){
    int nPoints=1;
    double satLiqDens,Vp;
    *liqDens=0;
    if(data->lDensCorr.form>0){
        if(data->baseProp.type==FF_Polymer)FF_LiqDensTait(&data->lDensCorr.form,data->lDensCorr.coef,&data->baseProp,&nPoints,T,P,&satLiqDens);
        else FF_PhysPropCorr(&data->lDensCorr.form,data->lDensCorr.coef,&data->baseProp.MW,&nPoints,T,&satLiqDens);
    }
    else FF_LiqDensSatRackett(&data->baseProp,&data->lDens.x,&data->lDens.y,T,&satLiqDens);
    if((*P>20e5)&&(data->baseProp.type!=FF_Polymer)){
        FF_Vp(T,data,&Vp);
        if((*P>Vp)&&(Vp>0)) FF_LiqDensChuehPrausnitz(&data->baseProp,T,P,&Vp,&satLiqDens,liqDens);
        else *liqDens=satLiqDens;
    }
    else *liqDens=satLiqDens;
}

//Liquid Cp. Bondi method
void CALLCONV FF_LiqCpBondi(const  FF_SubstanceData *data,const double *T,double *Cp){
    int nPoints=1;
    double Cp0,Tr;
    if(data->cp0Corr.form>0){
        Tr=*T/data->baseProp.Tc;
        FF_PhysPropCorr(&data->cp0Corr.form,data->cp0Corr.coef,&data->baseProp.MW,&nPoints,T,&Cp0);
        *Cp=R*1000*(1.586+0.49/(1-Tr)+data->baseProp.w*(4.2775+6.3*pow((1-Tr),0.33333)/Tr+0.4355/(1-Tr)))/data->baseProp.MW+Cp0;
    }
}


//Transport properties of liquids
//-------------------------------

//Lucas liquid viscosity pressure correction.
void CALLCONV FF_LiqViscPcorLucas(double *T,double *P,double *Pref,FF_BaseProp *data,double *rpVisc,double *visc){
    double Tr,Tr2,Tr3,Tr4,A,C,D,dPr;
    Tr=*T/data->Tc;
    Tr2=Tr*Tr;
    Tr3=Tr2*Tr;
    Tr4=Tr3*Tr;
    dPr=(*P- *Pref)/data->Pc;
    A=0.9991-(4.674e-4/(1.0523*pow(Tr,- 0.03877)-1.0513));
    C=-0.07921 + 2.1616*Tr - 13.404*Tr2 + 44.1706*Tr3 - 84.8291*Tr4 + 96.1209*Tr3*Tr2 - 59.8127*Tr3*Tr3 + 15.6719*Tr4*Tr3;
    D=(0.3257/pow(1.0039 - pow(Tr,2.573),0.2906)) - 0.2086;
    *visc=*rpVisc*(1+D*pow((dPr/2.118),A)/(1+C*data->w*dPr));
    //printf("A:%f C:%f D:%f ratio:%f\n",A,C,D,*visc/ *rpVisc);
}

//Liquid viscosity from T,P
void CALLCONV FF_LiqViscTP(double *T,double *P,FF_SubstanceData *data,double *liqVisc){
    int nPoints=1;
    double satLiqVisc,Vp;
    *liqVisc=0;
    if(data->lViscCorr.form>0){
        FF_PhysPropCorr(&data->lViscCorr.form,data->lViscCorr.coef,&data->baseProp.MW,&nPoints,T,&satLiqVisc);
        if(*P>10e5){
            FF_Vp(T,data,&Vp);
            if((*P>Vp)&&(Vp>0)) FF_LiqViscPcorLucas(T,P,&Vp,&data->baseProp,&satLiqVisc,liqVisc);
            else *liqVisc=satLiqVisc;
        }
        else *liqVisc=satLiqVisc;
    }
}

//Grunberg-Nissan method for mixture liquid viscosity
//BIPs are missing
void CALLCONV FF_MixLiqViscGrunberg(FF_MixData *mix,double *T,double *P,double x[],double *visc){
    double satVisc,vp,sVisc[mix->numSubs];
    int i,j,nPoints=1;
    *visc=0;
    for(i=0;i<mix->numSubs;i++){
        if(mix->lViscCorr[i].form>0){
            FF_PhysPropCorr(mix->lViscCorr[i].form,mix->lViscCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&satVisc);
            if((mix->vpCorr[i].form>0)&&(*P>10e5)){
                FF_PhysPropCorr(mix->vpCorr[i].form,mix->vpCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&vp);
                FF_LiqViscPcorLucas(T,P,&vp,&mix->baseProp[i],&satVisc,&sVisc[i]);
            }
            else sVisc[i]=satVisc;
        }
        else return;
    }
    for(i=0;i<mix->numSubs;i++){
        *visc=*visc+x[i]*log(sVisc[i]);
        for(j=0;j<mix->numSubs;j++) *visc=*visc+0.5*x[i]*x[j];//Necessary to add BIPs
    }
}

//Teja-Rice method for mixture liquid viscosity
//BIP is missing
void CALLCONV FF_MixLiqViscTeja(FF_MixData *mix,double *T,double *P,double x[],double *visc){
    //printf("T:%f P:%f x[0]:%f x[1]:%f\n",*T,*P,x[0],x[1]);
    int i,j,first,second,nPoints=1;
    double Vcij[mix->numSubs][mix->numSubs],TVcij[mix->numSubs][mix->numSubs],Vcm=0,Tcm=0,MWm=0,wm=0,em,e1,e2,mu1,mu2,T1,T2,satVisc,vp;
    *visc=0;
    for(i=0;i<mix->numSubs;i++){
        for(j=0;j<mix->numSubs;j++){
            Vcij[i][j]=pow(mix->baseProp[i].Vc,0.33333)+pow(mix->baseProp[j].Vc,0.33333);
            Vcij[i][j]=Vcij[i][j]*Vcij[i][j]*Vcij[i][j]/8;
            //TVcij[i][j]=psi[i][j]*pow((mix->baseProp[i].Tc*mix->baseProp[i].Vc*mix->baseProp[j].Tc*mix->baseProp[j].Vc),0.5);//To use with BIPs
            TVcij[i][j]=(mix->baseProp[i].Tc*mix->baseProp[i].Vc+mix->baseProp[j].Tc*mix->baseProp[j].Vc)/2;//Better equation if no BIPs (psi) are used
        }
    }
    for(i=0;i<mix->numSubs;i++){
        MWm=MWm+x[i]*mix->baseProp[i].MW;
        wm=wm+x[i]*mix->baseProp[i].w;
        for(j=0;j<mix->numSubs;j++){
            Vcm=Vcm+x[i]*x[j]*Vcij[i][j];
            Tcm=Tcm+x[i]*x[j]*TVcij[i][j];
        }
    }
    Tcm=Tcm/Vcm;

    //Finding the two main substances
    if(x[0]>x[1]){
        first=0;
        second=1;
    }
    else{
        first=1;
        second=0;
    }
    for(i=2;i<mix->numSubs;i++){
        if(x[i]>x[first]){
            second=first;
            first=i;
        }
        else if(x[i]>second) second=i;
    }

    //final calculation
    em=pow(Vcm,0.66667)/pow((Tcm*MWm),0.5);
    e1=pow(mix->baseProp[first].Vc,0.66667)/pow((mix->baseProp[first].Tc*mix->baseProp[first].MW),0.5);
    e2=pow(mix->baseProp[second].Vc,0.66667)/pow((mix->baseProp[second].Tc*mix->baseProp[second].MW),0.5);
    T1=*T*mix->baseProp[first].Tc/Tcm;
    T2=*T*mix->baseProp[second].Tc/Tcm;
    //printf("Vcm:%f Tcm:%f T1:%f T2:%f em:%f e1:%f e2:%f\n",Vcm,Tcm,T1,T2,em,e1,e2);

    if(mix->lViscCorr[first].form>0){//Calculation of main componend viscosity at equivalent T
        FF_PhysPropCorr(&mix->lViscCorr[first].form,mix->lViscCorr[first].coef,&mix->baseProp[first].MW,&nPoints,&T1,&satVisc);
        if((mix->vpCorr[first].form>0)&&(*P>10e5)){
            FF_PhysPropCorr(&mix->vpCorr[first].form,mix->vpCorr[first].coef,&mix->baseProp[first].MW,&nPoints,&T1,&vp);
            FF_LiqViscPcorLucas(&T1,P,&vp,&mix->baseProp[first],&satVisc,&mu1);
        }
        else mu1=satVisc;//If we do not have the Vp we use the saturated pressure viscosity
    }
    else return;
    if(mix->lViscCorr[second].form>0){//Calculation of second main componend viscosity at equivalent T
        FF_PhysPropCorr(&mix->lViscCorr[second].form,mix->lViscCorr[second].coef,&mix->baseProp[second].MW,&nPoints,&T2,&satVisc);
        if((mix->vpCorr[second].form>0)&&(*P>10e5)){
            FF_PhysPropCorr(&mix->vpCorr[second].form,mix->vpCorr[second].coef,&mix->baseProp[second].MW,&nPoints,&T2,&vp);
            FF_LiqViscPcorLucas(&T2,P,&vp,&mix->baseProp[second],&satVisc,&mu2);
        }
        else mu2=satVisc;
    }
    else return;
    //printf("mu1:%f mu2:%f\n",mu1,mu2);
    *visc=exp(log(e1*mu1)+(log(e2*mu2)-log(e1*mu1))*(wm-mix->baseProp[first].w)/(mix->baseProp[second].w-mix->baseProp[first].w))/em;
}

//Thermal conductivity of liquids. Latini method
void CALLCONV FF_LiquidThCondLatini(double *T,FF_BaseProp *data,double *thCond){
    double Tr,A,alpha,beta,gamma;
    Tr=*T/data->Tc;
    switch(data->type){
    case FF_Alkane:
        A=0.00350;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    case FF_Alkene:
    case FF_Alkyne:
        A=0.0361;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_Cycloalkane:
        A=0.031;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_Aromatic:
        A=0.0346;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_Alcohol:
    case FF_Polyol:
        A=0.00339;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    case FF_Ether:
        A=0.0385;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_Ketone:
        A=0.00383;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    case FF_Acid:
        A=0.00319;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    case FF_Ester:
        A=0.0415;
        alpha=1.2;
        beta=1.0;
        gamma=0.167;
        break;
    case FF_HaloAlkane:
    case FF_HaloAlkene:
        A=0.562;
        alpha=0.0;
        beta=0.5;
        gamma=-0.167;
        break;
    default:
        A=0.00383;
        alpha=1.2;
        beta=0.5;
        gamma=0.167;
        break;
    }
    *thCond=A*pow(data->Tb,alpha)*pow((1-Tr),0.38)/(pow(data->MW,beta)*pow(data->Tc,gamma)*pow(Tr,(1/6)));
}

//Liquid thermal conductivity from T
void CALLCONV FF_LiqThCondT(double *T,FF_SubstanceData *data,double *liqThCond){
    int nPoints=1;
    *liqThCond=0;
    if(data->lThCCorr.form>0){
        FF_PhysPropCorr(&data->lThCCorr.form,data->lThCCorr.coef,&data->baseProp.MW,&nPoints,T,liqThCond);
    }
    else FF_LiquidThCondLatini(T,&data->baseProp,liqThCond);
}

//Li method for mixture liquid thermal conductivity
void CALLCONV FF_MixLiqThCondLi(FF_MixData *mix,double *T,double *P,double x[],double *thCond){
    double sThCond[mix->numSubs],sMolarVol[mix->numSubs],molarVol=0,Tref=0,rhoRef=0,volFract[mix->numSubs],lambda;
    int i,j,nPoints=1,overTc=0;
    *thCond=0;
    for(i=0;i<mix->numSubs;i++){
        if(*T>mix->baseProp[i].Tc) overTc=1;
    }
    for(i=0;i<mix->numSubs;i++){
        if(mix->lThCCorr[i].form>0){
            FF_PhysPropCorr(&mix->lThCCorr[i].form,mix->lThCCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&sThCond[i]);
        }
        else FF_LiquidThCondLatini(T,&mix->baseProp[i],&sThCond[i]);
        if(overTc==1) sMolarVol[i]=mix->baseProp[i].Vc;
        else if(mix->lDensCorr[i].form>0){
            FF_PhysPropCorr(&mix->lDensCorr[i].form,mix->lDensCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&sMolarVol[i]);
        }
        else FF_LiqDensSatRackett(&mix->baseProp[i],&Tref,&rhoRef,T,&sMolarVol[i]);
        sMolarVol[i]=mix->baseProp[i].MW*1e-3/sMolarVol[i];
        molarVol=molarVol+x[i]*sMolarVol[i];
    }


    for(i=0;i<mix->numSubs;i++){
        volFract[i]=x[i]*sMolarVol[i]/molarVol;
        //printf("sThCond[%i]:%f sMolarVol[%i]:%f\n",i,sThCond[i],i,sMolarVol[i]);
    }
    for(i=0;i<mix->numSubs;i++){
        for(j=0;j<mix->numSubs;j++){
            lambda=2/(1/sThCond[i]+1/sThCond[j]);
            *thCond=*thCond+volFract[i]*volFract[j]*lambda;
        }
    }
}

//SurfaceTension, MacLeod-Sugden method. Very sensible to Parachor value
void CALLCONV FF_SurfTensMcLeod(double *T,FF_SubstanceData *data,double *surfTens){
    int nPoints=1;
    double lDens,Vp,gDens;
    *surfTens=0;
    if(data->baseProp.Pa>0){
        if(data->lDensCorr.form>0){
            FF_PhysPropCorr(&data->lDensCorr.form,data->lDensCorr.coef,&data->baseProp.MW,&nPoints,T,&lDens);
        }
        else if((data->lDens.x>0)&&(data->lDens.y>0)) FF_LiqDensSatRackett(&data->baseProp,&data->lDens.x,&data->lDens.y,T,&lDens);
        else return;

        if(data->gDensCorr.form>0){
            FF_PhysPropCorr(&data->gDensCorr.form,data->gDensCorr.coef,&data->baseProp.MW,&nPoints,T,&gDens);
        }
        else if(data->vpCorr.form>0){
            FF_PhysPropCorr(&data->vpCorr.form,data->vpCorr.coef,&data->baseProp.MW,&nPoints,T,&Vp);
        }
        else FF_VpAmbroseWalton(&data->baseProp,T,&Vp);
        gDens=0.3*data->baseProp.MW*Vp/(R* *T)*1e-3;
        *surfTens=pow((data->baseProp.Pa*(lDens-gDens)*1000/data->baseProp.MW),4);
        //printf("lDens:%f gDens:%f\n",lDens,gDens);
    }
}

//SurfaceTension, Sastri-Rao method. Very good approximation is obtained
void CALLCONV FF_SurfTensSastri(double *T,FF_BaseProp *data,double *surfTens){
    //double Tm=1- *T/data->Tc;
    //*surfTens=kb*data->Tc*pow((Av/data->Vc),0.666667)*(4.35+4.14*data->w)*pow(Tm,1.26)*(1+0.19*pow(Tm,0.5)-0.25*Tm);//Not good for alcohols
    //Use Sastri-Rao:
    double k,x,y,z,m;
    switch(data->type){
    case FF_Alcohol:
    case FF_Polyol:
        k=2.28;
        x=0.25;
        y=0.175;
        z=0;
        m=0.8;
        break;
    case FF_Acid:
        k=0.125;
        x=0.5;
        y=-1.5;
        z=1.85;
        m=1.222;
        break;
    default:
        k=0.158;
        x=0.5;
        y=-1.5;
        z=1.85;
        m=1.222;
        break;
    }
    *surfTens=k*pow(data->Pc*1e-5,x)*pow(data->Tb,y)*pow(data->Tc,z)*pow(((1- *T/data->Tc)/(1-data->Tb/data->Tc)),m)*1e-3;
}

//Surface tension from T
void CALLCONV FF_SurfTensT(double *T,FF_SubstanceData *data,double *surfTens){
    int nPoints=1;
    *surfTens=0;
    if (data->lSurfTCorr.form>0){
        FF_PhysPropCorr(&data->lSurfTCorr.form,data->lSurfTCorr.coef,&data->baseProp.MW,&nPoints,T,surfTens);
    }
    else if (data->lDensCorr.form>0){
        FF_SurfTensMcLeod(T,data,surfTens);
    }
    else FF_SurfTensSastri(T,&data->baseProp,surfTens);
}

//Linear method for mixture surface tension. Not very exact
void CALLCONV FF_MixLiqSurfTensLinear(FF_MixData *mix,double *T,double *P,double x[],double *surfTens){
    int i,j,nPoints=1;
    double sSurfTens[mix->numSubs];
    *surfTens=0;
    for(i=0;i<mix->numSubs;i++){
        if(mix->lSurfTCorr[i].form>0){
            FF_PhysPropCorr(&mix->lSurfTCorr[i].form,mix->lSurfTCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&sSurfTens[i]);
        }
        else FF_SurfTensSastri(T,&mix->baseProp[i],&sSurfTens[i]);
        *surfTens=*surfTens+x[i]*sSurfTens[i];
    }
}

//Winterfeld method for mixture surface tension.
void CALLCONV FF_MixLiqSurfTensWinterfeld(FF_MixData *mix,double *T,double *P,double x[],double *surfTens){
    int i,j,nPoints=1;
    double sSurfTens[mix->numSubs],sLiqDens[mix->numSubs],A=0,B=0,aux;
    *surfTens=0;
    for(i=0;i<mix->numSubs;i++){
        if(mix->lDensCorr[i].form>0) FF_PhysPropCorr(&mix->lDensCorr[i].form,mix->lDensCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&sLiqDens[i]);
        else FF_LiqDensSatRackett(&mix->baseProp[i],&mix->lDens[i].x,&mix->lDens[i].y,T,&sLiqDens[i]);
        sLiqDens[i]=sLiqDens[i]*1000/mix->baseProp[i].MW;
        if(mix->lSurfTCorr[i].form>0) FF_PhysPropCorr(&mix->lSurfTCorr[i].form,mix->lSurfTCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&sSurfTens[i]);
        else FF_SurfTensSastri(T,&mix->baseProp[i],&sSurfTens[i]);
        //printf("sSurfTens[%i]:%f form:%i sLiqDens:%f\n",i,sSurfTens[i],mix->lSurfTCorr[i].form,sLiqDens[i]);
    }
    for(i=0;i<mix->numSubs;i++){
        aux=x[i]/sLiqDens[i];
        B=B+aux;
        for(j=0;j<mix->numSubs;j++){
            A=A+x[i]*x[j]*pow((sSurfTens[i]*sSurfTens[j]),0.5)/(sLiqDens[i]*sLiqDens[j]);
        }
    }
    //printf("A:%f B:%f\n",A,B);
    *surfTens=A/(B*B);
}


//McLeod-Sugden method for mixture surface tension.
void CALLCONV FF_MixLiqSurfTensMcLeod(FF_MixData *mix,double *rhoL,double *rhoG,double x[],double y[],double *surfTens){
    int i;
    double Pl=0,Pg=0;
    for(i=0;i<mix->numSubs;i++){
        if(!(mix->baseProp[i].Pa>0)) return;
        Pl=Pl+x[i]*mix->baseProp[i].Pa;
        Pg=Pg+y[i]*mix->baseProp[i].Pa;
    }
    *surfTens=(*rhoL*Pl- *rhoG*Pg);
    *surfTens=pow(*surfTens,3.8);
}

//Transport properties of gases
//-----------------------------

//Gas viscosity pressure prediction/correction. Chung method
void CALLCONV FF_GasViscTVcpChung(double *T,double *V,FF_BaseProp *data,double *ldVisc,double *visc){
    double omega,Ta,sigma,Fc,muR,muR4,k;
    double E[10],y,G1,G2,eta2,eta1;
    int i;
    //Chung method
    Ta=1.2593* *T/data->Tc;
    sigma= 0.809*pow(data->Vc,0.33333);
    omega= 1.16145*pow(Ta,- 0.14874)+0.52487*exp(- 0.7732*Ta)+2.16178*exp(- 2.43787*Ta);
    if ((data->mu<99)&&(data->mu>0)) muR=131.3*data->mu/pow(data->Vc*1e6*data->Tc,0.5);
    else muR=0;
    muR4=muR*muR*muR*muR;
    if(data->type==FF_Water) k=0.076;//Water
    else if(data->type==FF_Alcohol) k=0.0682+4.74/data->MW;//alcohol
    else if(data->type==FF_Polyol) k=0.0682+4.74*2/data->MW;//Polyol
    else k=0;
    Fc=1-0.2756*data->w+0.059035*muR4+k;
    if(*ldVisc<=0){
        *ldVisc=26.692*Fc*pow((data->MW* *T),0.5)*1e-11/(sigma*sigma*omega);
    }

    double coef[10][4]={{6.324,50.412,-51.680,1189.0},{1.210e-3,-1.154e-3,-6.257e-3,0.03728},{5.283,254.209,-168.48,3898},{6.623,38.096,-8.464,31.42},
                        {19.745,7.63,-14.354,31.53},{-1.9,-12.537,4.985,-18.15},{24.275,3.45,-11.291,69.35},{0.7972,1.117,0.01235,-4.117},
                        {-0.2382,0.0677,-0.8163,4.025},{0.06863,0.3479,0.5926,-0.727}};
    for(i=0;i<10;i++) E[i]=coef[i][0]+coef[i][1]*data->w+coef[i][2]*muR4+coef[i][3]*k;
    y=data->Vc/(6* *V);
    G1=(1-0.5*y)/(pow((1-y),3));
    G2=(E[0]*((1-exp(-E[3]*y))/y)+E[1]*G1*exp(E[4]*y)+E[2]*G1)/(E[0]*E[3]+E[1]+E[2]);
    eta2=E[6]*y*y*G2*exp(E[7]+E[8]/Ta+E[9]/(Ta*Ta));
    eta1=pow(Ta,0.5)*(Fc*(1/G2+E[5]*y))/omega+eta2;
    if(eta1>1) *visc=*ldVisc*eta1;
    else *visc=*ldVisc;
    //printf("Ta:%f omega:%f muR:%f Fc:%f ldVisc:%f y:%f G1:%f G2:%f eta2:%f eta1:%f\n",Ta,omega,muR,Fc,*ldVisc,y,G1,G2,eta2,eta1);
    //for(i=0;i<10;i++) printf("E[%i]:%f\n",i,E[i]);
}

//Gas viscosity pressure prediction/correction. Lucas method
void CALLCONV FF_GasViscTPcpLucas(const double *T,const double *P,const FF_BaseProp *data,double *lpVisc,double *visc){
    double Z1,mur,Q,xi,Tr,Fp0,Fq0,Z2,Pr,a,b,c,d,e,f,Y,Fp,Fq;
    Tr=*T/data->Tc;
    Pr=*P/data->Pc;
    if((data->mu<999)&&(data->mu>0)) mur=52.46*data->mu*data->mu*data->Pc*1e-5/(data->Tc*data->Tc);
    else mur=0;
    if(mur<0.022) Fp0=1;
    else if(mur<0.075) Fp0=1+30.55*pow((0.292-data->Zc),1.72);
    else Fp0=1+30.55*pow((0.292-data->Zc),1.72)*fabs(0.96+0.1*(Tr-0.7));
    if(data->MW<4.04){
        if(data->MW<3) Q=0.76;//Hydrogen
        else if(data->MW<4.02) Q=1.38;//He
        Fq0=1.22*pow(Q,0.15);
    }
    else Fq0=1;
    xi=0.176*pow(data->Tc/(data->MW*data->MW*data->MW*data->Pc*data->Pc*data->Pc*data->Pc*1e-20),0.1666667)*1e-7;
    if(*lpVisc<=0){//If low pressure density is not supplied it is calculated
        *lpVisc=Fp0*Fq0*(0.807*pow(Tr,0.618)-0.357*exp(-0.449*Tr)+0.34*exp(-4.058*Tr)+0.018)*1e-14/xi;
    }
    //printf("lpVisc:%f\n",*lpVisc);
    Z1=*lpVisc*xi;
    if(Tr<=1){//The original method has been modified acording to Chemsep book
        Tr=1;
        //alpha=3.262+14.98*pow(Pr,5.508);
        //beta=1.39+5.746*Pr;
        //Z2=0.6+0.76*pow(Pr,alpha)+(6.99*pow(Pr,beta)-0.6)*(1-Tr);
    }
    a=1.245e-3*exp(5.1726*pow(Tr,-0.3286))/Tr;
    b=a*(1.6553*Tr-1.2723);
    c=0.4489*exp(3.0578*pow(Tr,-37.7332))/Tr;
    d=1.7368*exp(2.231*pow(Tr,-7.6351))/Tr;
    e=1.3088;
    f=0.9425*exp(-0.1853*pow(Tr,0.4489));
    Z2=Z1*(1+(a*pow(Pr,e))/(b*pow(Pr,f)+1/(1+c*pow(Pr,d))));

    Y=Z2/Z1;
    Fp=(1+(Fp0-1)*pow(Y,-3))/Fp0;
    Fq=(1+(Fq0-1)*(1/Y-0.007*pow(log(Y),4)))/Fq0;
    *visc=Z2*Fp*Fq/xi;
    //printf("lpVisc:%f visc:%f\n",*lpVisc,*visc);
}

//Gas viscosity from T,P
void CALLCONV FF_GasViscTP(double *T,double *P,FF_SubstanceData *data,double *gasVisc){
    int nPoints=1;
    double lpGasVisc=0;
    *gasVisc=0;
    if(data->gViscCorr.form>0){
        FF_PhysPropCorr(&data->gViscCorr.form,data->gViscCorr.coef,&data->baseProp.MW,&nPoints,T,&lpGasVisc);
    }
    FF_GasViscTPcpLucas(T,P,&data->baseProp,&lpGasVisc,gasVisc);
}

//Viscosity of gas mixtures. Wilke method
void CALLCONV FF_MixGasViscTPcpWilke(FF_MixData *mix,double *T,double *P,double y[],double *gVisc){
    int i,j,nPoints=1;
    double sGvisc[mix->numSubs],lpVisc,phi[mix->numSubs][mix->numSubs],sigma;
    lpVisc=0;
    *gVisc=0;
    for(i=0;i<mix->numSubs;i++){
        if(mix->gViscCorr[i].form>0){
            FF_PhysPropCorr(&mix->gViscCorr[i].form,mix->gViscCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&lpVisc);
        }
        FF_GasViscTPcpLucas(T,P,&mix->baseProp[i],&lpVisc,&sGvisc[i]);//Individual viscosities P corrected
        //printf("sGasVisc[%i]:%f\n",i,sGvisc[i]);
    }
    for(i=0;i<mix->numSubs;i++){
        for(j=0;j<mix->numSubs;j++){
            phi[i][j]=pow(sGvisc[i]/sGvisc[j],0.5)*pow(mix->baseProp[j].MW/mix->baseProp[i].MW,0.25);
            phi[i][j]=pow((1+phi[i][j]),2)/pow(8*(1+mix->baseProp[i].MW/mix->baseProp[j].MW),0.5);
        }
    }
    for(i=0;i<mix->numSubs;i++){
        sigma=0;
        for(j=0;j<mix->numSubs;j++){
            sigma=sigma+y[j]*phi[i][j];
        }
        *gVisc=*gVisc+y[i]*sGvisc[i]/sigma;
    }
}

//Viscosity of gas mixtures. Lucas method
void CALLCONV FF_MixGasViscTPcpLucas(FF_MixData *mix,double *T,double *P,double y[],double *gVisc){
    int i,imax;
    double Tcm=0,Vcm=0,Zcm=0,Pcm=0,MWm=0,Fp0m=0,Fq0m=0,MWmax=0,MWmin=1e10,mur,Trs,Fp0s,Qs,Fq0s,Tr,Pr,xi,lpVisc;
    double Z1,Z2,a,b,c,d,e,f,Ym,Fpm,Fqm;
    *gVisc=0;
    for(i=0;i<mix->numSubs;i++){
        if(mix->baseProp[i].Tc==0) return;
        else Tcm=Tcm+y[i]*mix->baseProp[i].Tc;
        if((mix->baseProp[i].Vc==0)&&(mix->baseProp[i].Zc>0)) mix->baseProp[i].Vc=mix->baseProp[i].Zc*R*mix->baseProp[i].Tc/mix->baseProp[i].Pc;
        if(mix->baseProp[i].Vc<=0) return;
        else Vcm=Vcm+y[i]*mix->baseProp[i].Vc;
        if((mix->baseProp[i].Zc==0)&&(mix->baseProp[i].Vc>0)) mix->baseProp[i].Zc=mix->baseProp[i].Pc*mix->baseProp[i].Vc/(R*mix->baseProp[i].Tc);
        if(mix->baseProp[i].Zc<=0) return;
        else Zcm=Zcm+y[i]*mix->baseProp[i].Zc;
        Pcm=R*Zcm*Tcm/Vcm;
        MWm=MWm+y[i]*mix->baseProp[i].MW;
        if((mix->baseProp[i].mu<999)&&(mix->baseProp[i].mu>0)) mur=52.46*mix->baseProp[i].mu*mix->baseProp[i].mu*mix->baseProp[i].Pc*1e-5/
                (mix->baseProp[i].Tc*mix->baseProp[i].Tc);
        else mur=0;
        Trs=*T/mix->baseProp[i].Tc;
        if(mur<0.022) Fp0s=1;
        else if(mur<0.075) Fp0s=1+30.55*pow((0.292-mix->baseProp[i].Zc),1.72);
        else Fp0s=1+30.55*pow((0.292-mix->baseProp[i].Zc),1.72)*fabs(0.96+0.1*(Trs-0.7));
        Fp0m=Fp0m+y[i]*Fp0s;
        if(mix->baseProp[i].MW<4.04){
            if(mix->baseProp[i].MW<3) Qs=0.76;//Hydrogen
            else if(mix->baseProp[i].MW<4.02) Qs=1.38;//He
            Fq0s=1.22*pow(Qs,0.15);
        }
        else Fq0s=1;
        Fq0m=Fq0m+y[i]*Fq0s;
        if(mix->baseProp[i].MW>MWmax){
            MWmax=mix->baseProp[i].MW;
            imax=i;
        }
        if(mix->baseProp[i].MW<MWmin) MWmin=mix->baseProp[i].MW;
    }
    if(((MWmax/MWmin)>9)&&(y[imax]>0.05)&&(y[imax]<0.7)) Fq0m=Fq0m*(1-0.01*pow((MWmax/MWmin),0.87));
    Tr=*T/Tcm;
    Pr=*P/Pcm;
    xi=0.176*pow(Tcm/(MWm*MWm*MWm*Pcm*Pcm*Pcm*Pcm*1e-20),0.1666667)*1e-7;
    lpVisc=Fp0m*Fq0m*(0.807*pow(Tr,0.618)-0.357*exp(-0.449*Tr)+0.34*exp(-4.058*Tr)+0.018)*1e-14/xi;

    Z1=lpVisc*xi;
    if(Tr<=1) Tr=1;
    a=1.245e-3*exp(5.1726*pow(Tr,-0.3286))/Tr;
    b=a*(1.6553*Tr-1.2723);
    c=0.4489*exp(3.0578*pow(Tr,-37.7332))/Tr;
    d=1.7368*exp(2.231*pow(Tr,-7.6351))/Tr;
    e=1.3088;
    f=0.9425*exp(-0.1853*pow(Tr,0.4489));
    Z2=Z1*(1+(a*pow(Pr,e))/(b*pow(Pr,f)+1/(1+c*pow(Pr,d))));

    Ym=Z2/Z1;
    Fpm=(1+(Fp0m-1)*pow(Ym,-3))/Fp0m;
    Fqm=(1+(Fq0m-1)*(1/Ym-0.007*pow(log(Ym),4)))/Fq0m;
    *gVisc=Z2*Fpm*Fqm/xi;
    //printf("Tcm:%f\n",Tcm);
}

//Gas low pressure thermal conductivity prediction
void CALLCONV FF_GasLpThCondTCpChung(double *T,double *Cp0,FF_BaseProp *data,double *ldThCond){
    double ldVisc,visc,alpha,Tr,beta,Z,psi,P;
    //Chung method
    ldVisc=0;
    P=1e5;
    FF_GasViscTPcpLucas(T,&P,data,&ldVisc,&visc);
    alpha=(*Cp0*data->MW*0.001-R)/R-1.5;
    beta=0.7862 - 0.7109*data->w + 1.3168 *data->w*data->w;
    Tr=*T/data->Tc;
    Z=2.0 + 10.5*Tr*Tr;
    psi=1+alpha*((0.215 + 0.28288 *alpha - 1.061 *beta + 0.26665*Z)/(0.6366 + beta*Z + 1.061*alpha*beta ));
    *ldThCond=3.75*psi*ldVisc*R*1000/data->MW;
}

//Gas thermal conductivity pressure correction
void CALLCONV FF_GasThCondTVcorChung(double *T,double *V,FF_BaseProp *data,double *ldThCond,double *thCond){
    int i;
    double Tr,y,B[7],G1,muR,muR4,k,G2;
    double coef[7][4]={{2.4166,7.4824e-1,-9.1858e-1,1.2172e2},{-5.0924e-1,-1.5094,-4.9991e1,6.9983e1},{6.6107,5.6207,6.4760e1,2.7039e1},{1.4543e1,-8.9139,-5.6379,7.4344e1},
                        {7.9274e-1,8.2019e-1,-6.9369e-1,6.3173},{-5.8634,1.2801e1,9.5893,6.5529e1},{9.1089e1,1.2811e2,-5.4217e1,5.2381e2}};
    Tr=*T/data->Tc;
    y=data->Vc/(6* *V);
    G1=(1-0.5*y)/(pow((1-y),3));
    if ((data->mu<99)&&(data->mu>0)) muR=131.3*data->mu/pow(data->Vc*1e6*data->Tc,0.5);
    else muR=0;
    muR4=muR*muR*muR*muR;
    if(data->type==FF_Water) k=0.076;//Water
    else if(data->type==FF_Alcohol) k=0.0682+4.74/data->MW;//alcohol
    else if(data->type==FF_Polyol) k=0.0682+4.74*2/data->MW;//Polyol
    else k=0;
    for(i=0;i<7;i++) B[i]=coef[i][0]+coef[i][1]*data->w+coef[i][2]*muR4+coef[i][3]*k;
    G2=(B[0]*(1-exp(-B[3]*y))/y+B[1]*G1*exp(B[4]*y)+B[2]*G1)/(B[0]*B[3]+B[1]+B[2]);
    *thCond=*ldThCond*(1/G2+B[5]*y)+3.586e-3*pow(data->Tc*1000/data->MW,0.5)*B[6]*y*y*pow(Tr,0.5)*G2/pow(data->Vc*1e6,0.6667);
    if(*thCond<*ldThCond) *thCond=*ldThCond;

    //printf("Ta:%f omega:%f muR:%f Fc:%f ldVisc:%f y:%f G1:%f G2:%f eta2:%f eta1:%f\n",Ta,omega,muR,Fc,*ldVisc,y,G1,G2,eta2,eta1);
    //for(i=0;i<10;i++) printf("E[%i]:%f\n",i,E[i]);
}

//Gas thermal conductivity from T,V
void CALLCONV FF_GasThCondTV(double *T,double *V,FF_SubstanceData *data,double *gasThCond){
    int nPoints=1;
    double Cp0,ldGasThCond=0;
    *gasThCond=0;
    if(data->gThCCorr.form>0){
        FF_PhysPropCorr(&data->gThCCorr.form,data->gThCCorr.coef,&data->baseProp.MW,&nPoints,T,&ldGasThCond);
    }
    else{
        if(data->cp0Corr.form>0){
            FF_PhysPropCorr(&data->cp0Corr.form,data->cp0Corr.coef,&data->baseProp.MW,&nPoints,T,&Cp0);
            FF_GasLpThCondTCpChung(T,&Cp0,&data->baseProp,&ldGasThCond);
        }
        else return;
    }
    FF_GasThCondTVcorChung(T,V,&data->baseProp,&ldGasThCond,gasThCond);
}

//Thermal conductivity of low pressure gas mixtures. Mason and Saxena method
//Correction for high pressure in Chemsep book
void CALLCONV FF_MixLpGasThCondTpMason(FF_MixData *mix,double *T,double y[],double *gThCond){
    int i,j,nPoints=1;
    double sGthCond[mix->numSubs],Cp0,Gamma[mix->numSubs],ratio,Tr[mix->numSubs],A[mix->numSubs][mix->numSubs],sigma;
    *gThCond=0;
    for(i=0;i<mix->numSubs;i++){
        if(mix->gThCCorr[i].form>0){
            FF_PhysPropCorr(&mix->gThCCorr[i].form,mix->gThCCorr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&sGthCond[i]);
        }
        else{
            if(mix->gThCCorr[i].form>0){
                FF_PhysPropCorr(&mix->cp0Corr[i].form,mix->cp0Corr[i].coef,&mix->baseProp[i].MW,&nPoints,T,&Cp0);
                FF_GasLpThCondTCpChung(T,&Cp0,&mix->baseProp[i],&sGthCond[i]);
            }
            else return;
        }
        //printf("lpGasThCond[%i]:%f\n",i,sGthCond[i]);
        Tr[i]=*T/mix->baseProp[i].Tc;
        Gamma[i]=210*pow(mix->baseProp[i].Tc*mix->baseProp[i].MW*mix->baseProp[i].MW*mix->baseProp[i].MW/
                         (mix->baseProp[i].Pc*mix->baseProp[i].Pc*mix->baseProp[i].Pc*mix->baseProp[i].Pc*1e-20),0.1666667);
     }
    for(i=0;i<mix->numSubs;i++){
        sigma=0;
        for(j=0;j<mix->numSubs;j++){
            ratio=Gamma[j]*(exp(0.0464*Tr[i])-exp(-0.2412*Tr[i]))/(Gamma[i]*(exp(0.0464*Tr[j])-exp(-0.2412*Tr[j])));
            A[i][j]=pow((1+pow(ratio,0.5)*pow(mix->baseProp[i].MW/mix->baseProp[j].MW,0.25)),2)/pow(8*(1+mix->baseProp[i].MW/mix->baseProp[j].MW),0.5);
            sigma=sigma+y[j]*A[i][j];
        }
        *gThCond=*gThCond+y[i]*sGthCond[i]/sigma;
    }
}


//Temperature,volume dependent, phase independent viscosity.
void CALLCONV FF_ViscosityTDens(int subsRef, double *T, double *molarDens,double eta[3]){
    int i,j;
    double MW;
    double rhoReducing,Treducing;
    double Tred,lnTred;
    double epsilon;//epsilon/kb
    double sigma;
    double CSred;//reduced cross section
    double b[9]={-19.572881,219.73999,-1015.3226,2471.01251,-3375.1717,2491.6597,-787.26086,14.085455,-0.34664158};//second virial coefficients
    double f1,g1,g2;
    double etai=0;//ideal part
    double etaf=0;//Rainwater-Friend part
    double etab=0;//Batschinski part
    double delta,delta0,tau;
    if (subsRef==1){//propane
        MW=44.0956;
        epsilon=263.88;
        sigma=0.49748;//nm
        Treducing=369.825;
        rhoReducing=5000;//mol/m3
        double a[4]={0.25104574,-0.47271238,0.0,0.060836515,0.0};//coefficients for propane
        double e[4][3]={{35.9873030195,-180.512188564,87.7124888223}, {-105.773052525,205.319740877,-129.21093261},
                        {58.9491587759,-129.7400331,76.6280419971},{-9.59407868475,21.0726986598,-14.3971968187}};
        Tred=*T/epsilon;
        lnTred=log(Tred);
        CSred=exp(a[0]+lnTred*(a[1]+lnTred*(a[2]+lnTred*(a[3]+lnTred*a[4]))));
        f1=1616.88405374;
        g1=2.50053938863;
        //g2=2.1517543007382995;//From CoolProp
        g2=0.860516059264;
        delta=*molarDens/rhoReducing;
        tau=*T/Treducing;
        delta0=g1*(1+g2*pow(tau,0.5));
        etab=0;
        for(i=0;i<4;i++)for(j=0;j<3;j++)etab=etab+e[i][j]*pow(delta,i+2)/pow(tau,j);
        etab=etab+f1*(delta/(delta0-delta)-delta/delta0);

    }
    else if(subsRef==2){//R134A
        MW=102.031;
        epsilon=299.363;
        sigma=4.68932;//nm
        Treducing=374.21;
        rhoReducing=5017.053;//mol/m3
        double a[2]={0.355404,-0.464337,0.0257353};
        CSred=exp(a[0]+lnTred*(a[1]+lnTred*a[2]));
        f1=214.76332;
        g1=3.163695636;
        g2=0.0;
    }

    etai=0.021357*pow(*T*MW,0.5)/(sigma*sigma*CSred);//micro Pa·s.
    etaf=b[0];
    for(i=1;i<7;i++){
        etaf=etaf+b[i]*pow(Tred,(-0.25*i));
        //printf("i%i pow:%f\n",i,pow(Tred,(-0.5*i)));
    }

    etaf=etaf+b[7]*pow(Tred,-2.5)+b[8]*pow(Tred,-5.5);//this is the reduced second viscosity Virial coefficient
    etaf=etai*etaf*Av*pow(sigma,3)*1e-27;//pass sigma^3 from nm to m3
    etaf=etaf* *molarDens;//microPa·s

    eta[0]=etai*1e-6;
    eta[1]=etaf*1e-6;
    eta[2]=etab*1e-6;
}

//Extended correponding states calculation for saturation properties
void CALLCONV FF_CorrespondingStatesSat(FF_SubstanceData *subs, FF_SubstanceData *ref, double *T, double result[2]){
    int nPoints=1;
    double Vl,lDens,target,T0,V0,lDens0,act,actPlus,der,error,lVisc0,ldVisc=0.0,ldVisc0=0.0,ldViscCorr=0.0,ldVisc0Corr=0.0,Vp,Vp0,lThC0=0.0,ldThC0=0.0,ldThC=0.0;
    double Cp0=0.0,Cp=0.0;
    double f=1.0,h=1.0,Fn=1.0;
    double rhoMolar;
    int i;
    result[0]=0.0;
    result[1]=0.0;

    if ((subs->swData.eos>0)&&(ref->swData.eos>0)&&(subs->gViscCorr.form==112)){//Nist method with correction coef.
        int eosType=FF_SWtype;
        double tau,delta,Arr,Z,tau0,delta0,answerL[3],answerG[3],ArrDer[6];
        char option='l';
        char state;
        double psi,eta[3];
        FF_VpEOS(&eosType,T,&subs->swData,&Vp);
        FF_VfromTPsw(T,&Vp,&subs->swData,&option,answerL,answerG,&state);
        tau=subs->swData.tRef/ *T;
        Vl=answerL[0];
        delta=1/(Vl*subs->swData.rhoRef);
        psi=subs->gViscCorr.coef[0]+delta*(subs->gViscCorr.coef[1]+delta*(subs->gViscCorr.coef[2]+delta*subs->gViscCorr.coef[3]));
        FF_ArrZfromTVsw(T,&Vl,&subs->swData,&Arr,&Z);//Arr and Z to match
        tau0=tau;//initial T for the reference
        T0=ref->swData.tRef/tau0;
        FF_VpEOS(&eosType,&T0,&ref->swData,&Vp0);
        FF_VfromTPsw(&T0,&Vp0,&ref->swData,&option,answerL,answerG,&state);
        V0=answerL[0];
        delta0=1/(V0*ref->swData.rhoRef);
        FF_ArrDerSW(&tau0,&delta0,&ref->swData,ArrDer);
        //printf("T:%f delta:%f Arr:%f, T0:%f Arr0:%f\n",*T,delta,Arr,T0,ArrDer[0]);
        while (fabs(Arr-ArrDer[0])>0.002){
            tau0=tau0+(Arr-ArrDer[0])/(1.5*ArrDer[3]);
            T0=ref->swData.tRef/tau0;
            FF_VpEOS(&eosType,&T0,&ref->swData,&Vp0);
            FF_VfromTPsw(&T0,&Vp0,&ref->swData,&option,answerL,answerG,&state);
            V0=answerL[0];
            delta0=1/(V0*ref->swData.rhoRef);
            FF_ArrDerSW(&tau0,&delta0,&ref->swData,ArrDer);
            printf("Arr:%f, Z:%f T:%f T0:%f delta0:%f Arr0:%f Z0:%f\n",Arr,Z,*T,T0,delta0,ArrDer[0],1+delta0*ArrDer[1]);
        }
        rhoMolar=psi*delta0*ref->swData.rhoRef;
        FF_ViscosityTDens(subs->gViscCorr.coef[5],&T0,&rhoMolar,eta);
        //printf("T0:%f P:%f psi:%f rhoMolar:%f Ideal:%f Friend:%f Batschinski:%f\n",T0,Vp0,psi,rhoMolar,eta[0],eta[1],eta[2]);
        ldVisc0=eta[0];
        lVisc0=eta[0]+eta[1]+eta[2];
    }
    else if((subs->vpCorr.form>0)&&(subs->lDensCorr.form>0)){
        FF_PhysPropCorr(&subs->vpCorr.form,&subs->vpCorr.coef,&subs->baseProp.MW,&nPoints,T,&Vp);
        FF_PhysPropCorr(&subs->lDensCorr.form,&subs->lDensCorr.coef,&subs->baseProp.MW,&nPoints,T,&lDens);
        Vl=subs->baseProp.MW/lDens*1000;
        target=*T*Vp*lDens*1000/subs->baseProp.MW;
        //printf("T:%f Target:%f\n",*T,target);
        T0=*T*ref->baseProp.Tc/subs->baseProp.Tc;
        //printf("T:%f refTc:%f subsTc:%f\n",*T,ref->baseProp.Tc,subs->baseProp.Tc);
        FF_PhysPropCorr(&ref->vpCorr.form,&ref->vpCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&Vp0);
        FF_PhysPropCorr(&ref->lDensCorr.form,&ref->lDensCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&lDens0);
        act=T0*Vp0*lDens0*1000/ref->baseProp.MW;
        error=(target-act)/target;
        //printf("T0:%f Actual:%f error:%f\n",T0,act,error);
        i=0;
        while ((fabs(error)>0.01)&&(i<10)){
            T0=T0+0.01;
            FF_PhysPropCorr(&ref->vpCorr.form,&ref->vpCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&Vp0);
            FF_PhysPropCorr(&ref->lDensCorr.form,&ref->lDensCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&lDens0);
            actPlus=T0*Vp0*lDens0*1000/ref->baseProp.MW;
            der=(actPlus-act)/0.01;
            if(der>0) T0=T0-0.01+0.6*error*target/der;
            else break;
            FF_PhysPropCorr(&ref->vpCorr.form,&ref->vpCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&Vp0);
            FF_PhysPropCorr(&ref->lDensCorr.form,&ref->lDensCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&lDens0);
            act=T0*Vp0*lDens0*1000/ref->baseProp.MW;
            error=(target-act)/target;
            //printf("T0:%f Actual:%f error:%f Previous der:%f\n",T0,act,error,der);
            i++;
        }
        V0=ref->baseProp.MW/lDens0*1000;
        if(ref->baseProp.Vc>0) FF_GasViscTVcpChung(&T0,&V0,&ref->baseProp,&ldVisc0,&ldVisc0Corr);//calculation of the ideal part of the reference
        else if(ref->baseProp.Pc>0) FF_GasViscTPcpLucas(&T0,&Vp0,&ref->baseProp,&ldVisc0,&ldVisc0Corr);
        else if(ref->gViscCorr.form>0) FF_PhysPropCorr(&ref->gViscCorr.form,ref->gViscCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&ldVisc0);
        if (ref->lViscCorr.form>0){//total viscosity of the reference
            FF_PhysPropCorr(&ref->lViscCorr.form,ref->lViscCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&lVisc0);
        }
    }
    else return;

    f=*T/T0;
    h=Vl/V0;
    Fn=pow((f*subs->baseProp.MW/ref->baseProp.MW),0.5)*pow(h,-0.66667);
    if(subs->baseProp.Vc>0) FF_GasViscTVcpChung(T,&Vl,&subs->baseProp,&ldVisc,&ldViscCorr);//Calculation of ideal part of the target
    else if(subs->baseProp.Pc>0)  FF_GasViscTPcpLucas(T,&Vp,&subs->baseProp,&ldVisc,&ldViscCorr);
    else if(subs->gViscCorr.form>0) FF_PhysPropCorr(&subs->gViscCorr.form,subs->gViscCorr.coef,&subs->baseProp.MW,&nPoints,T,&ldVisc);
    if((ldVisc>0)&&(ldVisc0>0)&&(lVisc0>0))  result[0]=Fn*(lVisc0-ldVisc0)+ldVisc;//Calculation of target viscosity
    //printf("T:%f T0:%f f:%f h:%f Fn:%f ldVisc:%f ldVisc0:%f lVisc0:%f rVisc:%f lVisc:%f\n",*T,T0,f,h,Fn,ldVisc,ldVisc0,lVisc0,Fn*(lVisc0-ldVisc0),result[0]);

    if (ref->lThCCorr.form>0){
        FF_PhysPropCorr(&ref->lThCCorr.form,ref->lThCCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&lThC0);
        if(ref->cp0Corr.form>0){
            FF_PhysPropCorr(&ref->cp0Corr.form,ref->cp0Corr.coef,&ref->baseProp.MW,&nPoints,&T0,&Cp0);
            FF_GasLpThCondTCpChung(&T0,&Cp0,&ref->baseProp,&ldThC0);
        }
        else if(ref->gThCCorr.form>0) FF_PhysPropCorr(&ref->gThCCorr.form,ref->gThCCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&ldThC0);
        if(subs->cp0Corr.form>0){
            FF_PhysPropCorr(&subs->cp0Corr.form,subs->cp0Corr.coef,&subs->baseProp.MW,&nPoints,T,&Cp);
            FF_GasLpThCondTCpChung(T,&Cp,&subs->baseProp,&ldThC);
        }
        else if(subs->gThCCorr.form>0) FF_PhysPropCorr(&subs->gThCCorr.form,subs->gThCCorr.coef,&subs->baseProp.MW,&nPoints,T,&ldThC);
        //result[1]=Fn*(lThC0-ldThC0)+ldThC;
        result[1]=lThC0-ldThC0+ldThC;
        //printf("T:%f T0:%f f:%f h:%f Fn:%f ldThC:%f ldThC0:%f lThC0:%f rThC:%f lThC:%f\n",*T,T0,f,h,Fn,ldThC,ldThC0,lThC0,Fn*(lThC0-ldThC0),result[1]);
    }
}

//General extended correponding states calculation, using arr
void CALLCONV FF_CorrespondingStates(FF_SubstanceData *subs, FF_SubstanceData *ref, double *T, double result[2]){
    int nPoints=1;
    int eosType=FF_NoType;
    double answerL[3],answerG[3];
    char option='l';
    char state;
    double Vl,lDens,target,tau,delta,Arr,Z,T0,V0,lDens0,delta0,tau0,act,actPlus,der,error,ArrDer[6],lVisc0,ldVisc=0.0,ldVisc0=0.0,ldViscCorr=0.0,ldVisc0Corr=0.0,Vp,Vp0,lThC0=0.0,ldThC0=0.0,ldThC=0.0;
    double Cp0=0.0,Cp=0.0;
    double f=1.0,h=1.0,Fn=1.0;
    int i;
    result[0]=0.0;
    result[1]=0.0;
    if ((subs->model==FF_SWtype)&&(ref->model==FF_SWtype)) eosType=FF_SWtype;
    else if ((subs->model==FF_SAFTtype)&&(ref->model==FF_SAFTtype)) eosType=FF_SAFTtype;
    else if ((subs->model==FF_CubicType)&&(ref->model==FF_CubicType)) eosType=FF_CubicType;
    switch(eosType){

    case FF_SWtype:
        FF_VpEOS(&eosType,T,&subs->swData,&Vp);
        FF_VfromTPsw(T,&Vp,&subs->swData,&option,answerL,answerG,&state);
        tau=subs->swData.tRef/ *T;
        Vl=answerL[0];
        delta=1/(Vl*subs->swData.rhoRef);
        FF_ArrZfromTVsw(T,&Vl,&subs->swData,&Arr,&Z);//Arr and Z to match
        tau0=tau;//initial T for the reference
        T0=ref->swData.tRef/tau0;
        FF_VpEOS(&eosType,&T0,&ref->swData,&Vp0);
        FF_VfromTPsw(&T0,&Vp0,&ref->swData,&option,answerL,answerG,&state);
        V0=answerL[0];
        delta0=1/(V0*ref->swData.rhoRef);
        FF_ArrDerSW(&tau0,&delta0,&ref->swData,ArrDer);
        //printf("T:%f Arr:%f, T0:%f Arr0:%f\n",*T,Arr,T0,ArrDer[0]);
        while (fabs(Arr-ArrDer[0])>0.002){
            tau0=tau0+(Arr-ArrDer[0])/(1.5*ArrDer[3]);
            T0=ref->swData.tRef/tau0;
            FF_VpEOS(&eosType,&T0,&ref->swData,&Vp0);
            FF_VfromTPsw(&T0,&Vp0,&ref->swData,&option,answerL,answerG,&state);
            V0=answerL[0];
            delta0=1/(V0*ref->swData.rhoRef);
            FF_ArrDerSW(&tau0,&delta0,&ref->swData,ArrDer);
            //printf("Arr:%f, Z:%f T:%f T0:%f Arr0:%f Z0:%f\n",Arr,Z,*T,T0,ArrDer[0],1+delta0*ArrDer[1]);
        }
        break;
    case FF_SAFTtype:
        FF_VpEOS(&eosType,T,&subs->saftData,&Vp);
        FF_VfromTPSAFT(T,&Vp,&subs->saftData,&option,answerL,answerG,&state);
        Vl=answerL[0];
        FF_ArrZfromTVSAFT(T,&Vl,&subs->saftData,&Arr,&Z);//Arr and Z to match
        T0=*T*ref->baseProp.Tc/subs->baseProp.Tc;//Initial T for the reference
        FF_VpEOS(&eosType,&T0,&ref->saftData,&Vp0);
        FF_VfromTPSAFT(&T0,&Vp0,&ref->saftData,&option,answerL,answerG,&state);
        V0=answerL[0];
        FF_ArrDerSAFT(&T0,&V0,&ref->saftData,ArrDer);
        //printf("Arr:%f, T:%f T0:%f ArrRef:%f\n",Arr,*T,T0,ArrDer[0]);
        while (fabs(Arr-ArrDer[0])>0.002){
            T0=T0+fabs(Arr-ArrDer[0])/(1.5*ArrDer[3]);
            FF_VpEOS(&eosType,&T0,&ref->saftData,&Vp0);
            FF_VfromTPSAFT(&T0,&Vp0,&ref->saftData,&option,answerL,answerG,&state);
            V0=answerL[0];
            FF_ArrDerSAFT(&T0,&V0,&ref->saftData,ArrDer);
            //printf("Arr:%f, Z:%f T:%f T0:%f rho0:%f ArrRef:%f Zref:%f\n",Arr,Z,*T,T0,1/V0,ArrDer[0],1-V0*ArrDer[1]);
        }
        break;
    }
    //V0=ref->baseProp.MW/lDens0*1000;

    f=*T/T0;
    h=Vl/V0;
    Fn=pow((f*subs->baseProp.MW/ref->baseProp.MW),0.5)*pow(h,-0.66667);
    //Fn=pow(f,1.0)*pow(subs->baseProp.MW/ref->baseProp.MW,1.0)*pow(h,-1.2);
    if(ref->baseProp.Vc>0) FF_GasViscTVcpChung(&T0,&V0,&ref->baseProp,&ldVisc0,&ldVisc0Corr);
    else if(ref->baseProp.Pc>0) FF_GasViscTPcpLucas(&T0,&Vp0,&ref->baseProp,&ldVisc0,&ldVisc0Corr);
    else if(ref->gViscCorr.form>0) FF_PhysPropCorr(&ref->gViscCorr.form,ref->gViscCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&ldVisc0);
    else return;
    if(subs->baseProp.Vc>0) FF_GasViscTVcpChung(T,&Vl,&subs->baseProp,&ldVisc,&ldViscCorr);
    else if(subs->baseProp.Pc>0)  FF_GasViscTPcpLucas(T,&Vp,&subs->baseProp,&ldVisc,&ldViscCorr);
    else if(subs->gViscCorr.form>0) FF_PhysPropCorr(&subs->gViscCorr.form,subs->gViscCorr.coef,&subs->baseProp.MW,&nPoints,T,&ldVisc);
    else return;
    if (ref->lViscCorr.form>0){
        FF_PhysPropCorr(&ref->lViscCorr.form,ref->lViscCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&lVisc0);
        result[0]=Fn*(lVisc0-ldVisc0)+ldVisc;
        //result[0]=Fn*lVisc0+ldVisc;
        //printf("T:%f T0:%f f:%f h:%f Fn:%f ldVisc:%f ldVisc0:%f lVisc0:%f rVisc:%f lVisc:%f\n",*T,T0,f,h,Fn,ldVisc,ldVisc0,lVisc0,Fn*(lVisc0-ldVisc0),result[0]);
    }
    if (ref->lThCCorr.form>0){
        FF_PhysPropCorr(&ref->lThCCorr.form,ref->lThCCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&lThC0);
        if(ref->cp0Corr.form>0){
            FF_PhysPropCorr(&ref->cp0Corr.form,ref->cp0Corr.coef,&ref->baseProp.MW,&nPoints,&T0,&Cp0);
            FF_GasLpThCondTCpChung(&T0,&Cp0,&ref->baseProp,&ldThC0);
        }
        else if(ref->gThCCorr.form>0) FF_PhysPropCorr(&ref->gThCCorr.form,ref->gThCCorr.coef,&ref->baseProp.MW,&nPoints,&T0,&ldThC0);
        if(subs->cp0Corr.form>0){
            FF_PhysPropCorr(&subs->cp0Corr.form,subs->cp0Corr.coef,&subs->baseProp.MW,&nPoints,T,&Cp);
            FF_GasLpThCondTCpChung(T,&Cp,&subs->baseProp,&ldThC);
        }
        else if(subs->gThCCorr.form>0) FF_PhysPropCorr(&subs->gThCCorr.form,subs->gThCCorr.coef,&subs->baseProp.MW,&nPoints,T,&ldThC);
        //result[1]=Fn*(lThC0-ldThC0)+ldThC;
        result[1]=lThC0-ldThC0+ldThC;
        //printf("T:%f T0:%f f:%f h:%f Fn:%f ldThC:%f ldThC0:%f lThC0:%f rThC:%f lThC:%f\n",*T,T0,f,h,Fn,ldThC,ldThC0,lThC0,Fn*(lThC0-ldThC0),result[1]);
    }
}
