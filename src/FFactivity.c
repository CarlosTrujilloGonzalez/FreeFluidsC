/*
 * FFactivity.c
 *
 *  Created on: 10/04/2017
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

// contains mainly calculations for activity coefficients
//=======================================================

#include <math.h>
#include <stdio.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFactivity.h"

//Calculates activity coefficients according to Wilson (modified) equation, at given T and composition
void CALLCONV FF_ActivityWilson(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[],
                                const enum FF_IntParamForm *form,const double *T,const double x[],double gamma[],double *gE){
    double V[*numSubs];//will contain the saturated liquid volume for each substance
    double Ru=0;//The gas constant in the necessary units
    double RuTinv;//1/(Ru*T) to speed calculation
    double Lambda[*numSubs][*numSubs];
    double LambdaX[*numSubs];
    double lnGamma[*numSubs];
    int i,j;//the loop variables
    *gE=0;
    if (*form!=FF_Ant3) for(i=0;i<*numSubs;i++) V[i]=R*baseProp[i].Tc/baseProp[i].Pc*pow(baseProp[i].Zra,(1+pow(1-*T/baseProp[i].Tc,0.2857)));
    //printf("V:%f %f\n",V[0],V[1]);
    switch(*form){
    case FF_Pol1K://if energy supplied in K, we use R=1
    case FF_Pol2K:
        Ru=1;
        break;
    case FF_Pol1J://if energy supplied in J/mol
    case FF_Pol2J:
        Ru=8.314472;
        break;
    case FF_Pol1C://if in calories/mol
    case FF_Pol2C:
        Ru=1.98588;
        break;
    }
    RuTinv=1/(Ru * *T);
    for (i=0;i<*numSubs;i++){
        for(j=0;j<*numSubs;j++){
            switch(*form){
            case FF_Pol1K:
            case FF_Pol1J:
            case FF_Pol1C:
                Lambda[i][j]=V[j]/V[i]*exp(-(pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]* *T+
                            pintParam[(i* *numSubs+j)*6+2]* *T * *T)*RuTinv);
                break;
            case FF_Pol2K:
            case FF_Pol2J:
            case FF_Pol2C:
                Lambda[i][j]=V[j]/V[i]*exp(-(pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]* *T+
                            pintParam[(i* *numSubs+j)*6+2]/ *T)*RuTinv);
                break;
            case FF_Ant3:
                Lambda[i][j]=exp(pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]/ *T+pintParam[(i* *numSubs+j)*6+2]*log(*T)+
                             pintParam[(i* *numSubs+j)*6+3]* *T+ pintParam[(i* *numSubs+j)*6+4]/ (*T * *T));//extended Antoine form
                break;
            }
            //printf("Lambda:%f\n",Lambda[i][j]);
        }
    }
    for (i=0;i<*numSubs;i++){
        LambdaX[i]=0;
        for(j=0;j<*numSubs;j++){
            LambdaX[i]=LambdaX[i]+Lambda[i][j]*x[j];
        }
    }
    for (i=0;i<*numSubs;i++){
        lnGamma[i]=1-log(LambdaX[i]);
        for(j=0;j<*numSubs;j++){
            lnGamma[i]=lnGamma[i]-x[j]*Lambda[j][i]/LambdaX[j];
        }
        gamma[i]=exp(lnGamma[i]);
        *gE=*gE+x[i]*lnGamma[i];
    }
}

//Calculates activity coefficients according to NRTL equation, at given T and composition
//the first four parameters are for the calculation of tau, the 2 last ones for the calculation of alpha
void CALLCONV FF_ActivityNRTL(const int *numSubs,const double pintParam[],const enum FF_IntParamForm *form,const double *T,const double x[],
                              double gamma[],double *gE){
    double Ru=0;//The gas constant in different units
    double RuTinv;//1/(Ru*T) to speed calculation
    double tau[*numSubs][*numSubs];
    double alpha[*numSubs][*numSubs];
    double G[*numSubs][*numSubs];
    double xG[*numSubs],xGtau[*numSubs];
    double lnGamma[*numSubs];
    int i,j;//the loop variables
    *gE=0;
    switch(*form){
    case FF_Pol1K://if energy supplied in K, we use R=1
    case FF_Pol2K:
        Ru=1;
        break;
    case FF_Pol1J://if energy supplied in J/mol
    case FF_Pol2J:
        Ru=8.314472;
        break;
    case FF_Pol1C://if in calories/mol
    case FF_Pol2C:
        Ru=1.98588;
        break;
    }
    RuTinv=1/(Ru * *T);
    for (i=0;i<*numSubs;i++){
        for(j=0;j<*numSubs;j++){//we switch between the different calculations of tau, and alpha
            switch(*form){
            case FF_Pol1K:
            case FF_Pol1J:
            case FF_Pol1C://polynomial 1 form: a+b*T+c*T^2
                tau[i][j]=(pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]* *T+pintParam[(i* *numSubs+j)*6+2]* *T * *T)*RuTinv;
                alpha[i][j]=pintParam[(i* *numSubs+j)*6+4];
                break;
            case FF_Pol2K:
            case FF_Pol2J:
            case FF_Pol2C://polynomial 2 form: a+b*T+c/T
                tau[i][j]=(pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]* *T+pintParam[(i* *numSubs+j)*6+2]/ *T)*RuTinv;//polynomial 2 form
                alpha[i][j]=pintParam[(i* *numSubs+j)*6+4];
                break;
            case FF_Ant1://Antoine form: a +b/T+c*ln(T)+d*T
                tau[i][j]=pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]/ *T+pintParam[(i* *numSubs+j)*6+2]*log(*T)+
                                                        pintParam[(i* *numSubs+j)*6+3]* *T;//extended Antoine form
                alpha[i][j]=pintParam[(i* *numSubs+j)*6+4]+pintParam[(i* *numSubs+j)*6+5]* (*T-273.15);
                break;
            }
            G[i][j]=exp(-alpha[i][j]*tau[i][j]);
            //printf("alpha:%f tau:%f\n",alpha[i][j],tau[i][j]);
        }
    }
    for (i=0;i<*numSubs;i++){
        xG[i]=0;
        xGtau[i]=0;
        for(j=0;j<*numSubs;j++){
            xG[i]=xG[i]+x[j]*G[j][i];
            xGtau[i]=xGtau[i]+x[j]*G[j][i]*tau[j][i];
        }
    }
    for (i=0;i<*numSubs;i++){
        lnGamma[i]=xGtau[i]/xG[i];
        for(j=0;j<*numSubs;j++){
            lnGamma[i]=lnGamma[i]+x[j]*G[i][j]/xG[j]*(tau[i][j]-xGtau[j]/xG[j]);
        }
        gamma[i]=exp(lnGamma[i]);
        *gE=*gE+x[i]*lnGamma[i];
    }


}

//Calculates activity coefficients according to UNIQUAQ equation, at given T and composition
void CALLCONV FF_ActivityUNIQUAC(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[],
                                 const enum FF_IntParamForm *form,const double *T,const double x[],double gamma[],double *gE){
    double Ru=0;//The gas constant in different units
    double RuTinv;//1/(Ru*T) to speed calculation
    double R=0,Q=0,Qres=0;//mean molecular volume and surface
    double Phi[*numSubs],Theta[*numSubs],ThetaRes[*numSubs];//the volume and surface fractions of all individual substances
    double tau[*numSubs][*numSubs];
    double S[*numSubs];
    double lnGamma[*numSubs],lnGammaC[*numSubs],lnGammaR[*numSubs];
    int i,j;//the loop variables
    *gE=0;

    switch(*form){
    case FF_Pol1K://if energy supplied in K, we use R=1
    case FF_Pol2K:
        Ru=1;
        break;
    case FF_Pol1J://if energy supplied in J/mol
    case FF_Pol2J:
        Ru=8.314472;
        break;
    case FF_Pol1C://if in calories/mol
    case FF_Pol2C:
        Ru=1.98588;
        break;
    }
    RuTinv=1/(Ru * *T);
    //combinatorial part calculation
    for (i=0;i<*numSubs;i++){
        R=R+x[i]*baseProp[i].r;
        Q=Q+x[i]*baseProp[i].q;
        Qres=Qres+x[i]*baseProp[i].qRes;
    }
    for (i=0;i<*numSubs;i++){
        Phi[i]=x[i]*baseProp[i].r/R;
        Theta[i]=x[i]*baseProp[i].q/Q;
        lnGammaC[i]=log(Phi[i]/x[i])+1-Phi[i]/x[i]-5*baseProp[i].q*(log(Phi[i]/Theta[i])+1-Phi[i]/Theta[i]);
        //printf("GammaC: %f\n",exp(lnGammaC[i]));
    }
    //residual part
    for (i=0;i<*numSubs;i++){
        ThetaRes[i]=x[i]*baseProp[i].qRes/Qres;
        for(j=0;j<*numSubs;j++){
            switch(*form){
            case FF_Pol1K:
            case FF_Pol1J:
            case FF_Pol1C:
                tau[i][j]=exp(-(pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]* *T+pintParam[(i* *numSubs+j)*6+2]* *T * *T)*RuTinv);//polynomial form
                break;
            case FF_Pol2K:
            case FF_Pol2J:
            case FF_Pol2C:
                tau[i][j]=exp(-(pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]* *T+pintParam[(i* *numSubs+j)*6+2]/ *T)*RuTinv);//polynomial form 2
                break;
            case FF_Ant3:
                tau[i][j]=exp(pintParam[(i* *numSubs+j)*6]+pintParam[(i* *numSubs+j)*6+1]/ *T+pintParam[(i* *numSubs+j)*6+2]*log(*T)+
                                    pintParam[(i* *numSubs+j)*6+3]* *T+pintParam[(i* *numSubs+j)*6+4]/ (*T * *T));//extended Antoine form
                break;
            }
        }
    }
    for (i=0;i<*numSubs;i++){
        S[i]=0;
        for(j=0;j<*numSubs;j++){
            S[i]=S[i]+ThetaRes[j]*tau[j][i];
            //S[i]=S[i]+Theta[j]*tau[j][i];
        }
    }
    for (i=0;i<*numSubs;i++){
        lnGammaR[i]=1-log(S[i]);
        for(j=0;j<*numSubs;j++){
            lnGammaR[i]=lnGammaR[i]-tau[i][j]*ThetaRes[j]/S[j];
            //lnGammaR[i]=lnGammaR[i]-tau[i][j]*Theta[j]/S[j];
        }
        lnGammaR[i]=baseProp[i].q*lnGammaR[i];
        //printf("GammaR: %f\n",exp(lnGammaR[i]));
        lnGamma[i]=lnGammaC[i]+lnGammaR[i];
        gamma[i]=exp(lnGamma[i]);
        *gE=*gE+x[i]*lnGamma[i];
    }
}

//Calculates fugacity and activity coefficients, at given T and composition, from an activity model
void CALLCONV FF_PhiFromActivity(const int *numSubs,enum FF_ActModel *model,const  FF_BaseProp baseProp[],const double pintParam[],
                                 const enum FF_IntParamForm *form,const bool *useVp,const enum FF_EosType *eosType,const void *data,
                                 const double *T,const double *P,const double x[],double gamma[],double phi[],double *gE){
    //First the calculation of the activity coefficient
    switch(*model){
    case Wilson:
        FF_ActivityWilson(numSubs,baseProp,pintParam,form,T,x,gamma,gE);
        break;
    case NRTL:
        FF_ActivityNRTL(numSubs,pintParam,form,T,x,gamma,gE);
        break;
    case UNIQUAC:
        FF_ActivityUNIQUAC(numSubs,baseProp,pintParam,form,T,x,gamma,gE);
        break;
    }
    //Now the calculation of the reference fugacity at same T and P for each component
    int i;
    if (*useVp==false){
         FF_SaftEOSdata *dataP;
         FF_SWEOSdata *dataS;
         FF_CubicEOSdata *dataC;
        char option='b',state;
        double answerL[3],answerG[3];
        double phiL0[*numSubs];
        switch (*eosType)
        {
            case FF_SAFTtype:
                dataP=data;
                for(i=0;i<*numSubs;i++){
                    dataP=dataP+i;
                    FF_VfromTPeos(eosType,T,P,dataP,&option,answerL,answerG,&state);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    phi[i]=phiL0[i]*gamma[i];
                }
                break;
            case FF_SWtype:
                dataS=data;
                for(i=0;i<*numSubs;i++){
                    dataS=dataS+i;
                    FF_VfromTPeos(eosType,T,P,dataS,&option,answerL,answerG,&state);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    phi[i]=phiL0[i]*gamma[i];
                }
                break;
            default:
                dataC=data;
                for(i=0;i<*numSubs;i++){
                    dataC=dataC+i;
                    //printf("eos:%i dataC:%f %f %f\n",eos[i],dataC->Tc,dataC->Pc,dataC->w);
                    FF_VfromTPeos(eosType,T,P,dataC,&option,answerL,answerG,&state);
                    //printf("Vl,Z:%f %f\n",answerL[0],answerL[2]);
                    phiL0[i]=exp(answerL[1]+answerL[2]-1)/answerL[2];
                    //printf("phiL0:%f\n",phiL0[i]);
                    phi[i]=phiL0[i]*gamma[i];
                }
                break;
        }
    }
    //In case of using Vp, it makes no sense to calculate the reference fugacity at Vp, and later the Poynting factor, to change to actual pressure
    //We can do directly the calculation at actual pressure. It makes sense only if assuming that reference fugacity is Vp
    else{
         FF_Correlation *vpCorr;
        vpCorr=data;
        double Vp;
        int nPoints=1;
        for(i=0;i<*numSubs;i++){
            vpCorr=vpCorr+i;
            FF_PhysPropCorr(&vpCorr->form,vpCorr->coef,&baseProp[i].MW,&nPoints,T,&Vp);
            //printf("Vp:%f\n",Vp);
            phi[i]=gamma[i]*Vp/ *P;
        }
    }
}

