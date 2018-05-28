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

//Calculates the common data, independent of the molar fraction and temperature for the UNIFAC model.
//In data, we receive for each substance, its composition in number of each subgroup. So each register is: id substance, id subgroup, num of ocurrences
//In uni we answer all the calculations
void CALLCONV FF_UNIFACParams(int numData, const int data[][3], FF_UnifacData *uni){
    int i,j,k,newSubs,newSubg;
    //We need to count the different subgroups used, and to store their id
    int numSubgroups=0;
    //We need also to count the different substances, storing it original id
    int numSubs=0;
    int substance[20];//Will contain the original id of the substances
    for (i=0;i<20;i++) for(j=0;j<30;j++) uni->subsSubg[i][j]=0;//We fill with 0 the conmposition of the substances in subgroups
    //We fill the first register with the first data
    uni->subgroup[0][0]=data[0][1];//this will contain the list of subgroups and their corresponding group
    numSubgroups=1;
    substance[0]=data[0][0];
    numSubs=1;
    uni->subsSubg[0][0]=data[0][2];//this will contain the ocurrence of subgroups in the substances

    for (i=1; i<numData;i++){//We follow all the descriptors
        newSubs=1;//We say that the substance is new
        for (j=0;j<numSubs;j++){//We follow all the already defined substances
            if (data[i][0]==substance[j]){
                newSubs=0;//We say that the substance is not new
                break;
            }
        }
        if (newSubs){//If not found we add the substance
            j=numSubs;//j must contain the index of the substance
            substance[j]=data[i][0];
            numSubs++;
        }
        newSubg=1;//we say the subgroup is new
        for (k=0; k<numSubgroups;k++){//We follow the already defined subgroups
            if (data[i][1]==uni->subgroup[k][0]){
                newSubg=0;//we say that the subgroup is not new
                break;
            }
        }
        if (newSubg){//if not found we add the subgroups
            k=numSubgroups;
            uni->subgroup[k][0]=data[i][1];
            numSubgroups++;
        }
        uni->subsSubg[j][k]=data[i][2];
    }
    //we recover from the corresponding file the subgroups information, and fill the subgroup and subgData table
    FILE *f;
    unsigned sg,g,g1,g2;
    double r,q,A12,B12,C12,A21,B21,C21;
    if (uni->model==FF_UNIFACStd) f=fopen("UnifacSubgStd.txt","r");
    else if ((uni->model==FF_UNIFACPSRK)||(uni->model==FF_EntropicFV)||(uni->model==FF_UNIFACZM)) f=fopen("UnifacSubgPSRK.txt","r");
    else if (uni->model==FF_UNIFACDort) f=fopen("UnifacSubgDort.txt","r");
    else if (uni->model==FF_UNIFACNist) f=fopen("UnifacSubgNist.txt","r");
    for (i=0;i<numSubgroups;i++){
        do{
            fscanf(f,"%3lu%3lu%6lf%6lf\n",&sg,&g,&r,&q);
            if(uni->subgroup[i][0]==sg){
                uni->subgroup[i][1]=g;
                uni->subgData[i][0]=r;
                uni->subgData[i][1]=q;
                break;
            }
        } while (sg<300);
        //printf("%i %i %i %f %f\n",i,uni->subgroup[i][0],uni->subgroup[i][1],uni->subgData[i][0],uni->subgData[i][1]);
        rewind(f);
    }
    fclose(f);

    //we recover now the subgroup interaction and fill the subgroupInter table
    for (i=0;i<numSubgroups;i++){
        for (j=0;j<numSubgroups;j++){
            uni->subgInt[i][j][0]=0;
            uni->subgInt[i][j][1]=0;
            uni->subgInt[i][j][2]=0;
        }
    }
    if (uni->model==FF_UNIFACStd) f=fopen("UnifacInterStd.txt","r");
    else if ((uni->model==FF_UNIFACPSRK)||(uni->model==FF_EntropicFV)||(uni->model==FF_UNIFACZM)) f=fopen("UnifacInterPSRK.txt","r");
    else if (uni->model==FF_UNIFACDort) f=fopen("UnifacInterDort.txt","r");
    else if (uni->model==FF_UNIFACNist) f=fopen("UnifacInterNist.txt","r");
    if (f==NULL) printf("Error\n");
    if (uni->model==FF_UNIFACStd){
        for (i=0;i<numSubgroups;i++){
            for (j=0;j<numSubgroups;j++){
                if (uni->subgroup[i][1]<uni->subgroup[j][1]){
                    do{
                        fscanf(f,"%03lu%03lu%10lf%10lf\n",&g1,&g2,&A12,&A21);
                        //printf("%i %i %f %f\n",g1,g2,A12,A21);
                        if ((g1==uni->subgroup[i][1])&& (g2==uni->subgroup[j][1])){
                            uni->subgInt[i][j][0]=A12;
                            uni->subgInt[j][i][0]=A21;
                            //printf("%i %i %f %f\n",g1,g2,uni->subgInt[i][j][0],uni->subgInt[j][i][0]);
                            break;
                        }
                    } while(g1<85);
                    rewind(f);
                }
            }
        }
    }
    else if ((uni->model==FF_UNIFACPSRK)||(uni->model==FF_UNIFACDort)||(uni->model==FF_UNIFACNist)||(uni->model==FF_UNIFACZM)||(uni->model==FF_EntropicFV)){
        for (i=0;i<numSubgroups;i++){
            for (j=0;j<numSubgroups;j++){
                if (uni->subgroup[i][1]<uni->subgroup[j][1]){
                    do{
                        fscanf(f,"%03lu%03lu%10lf%10lf%10lf%10lf%10lf%10lf\n",&g1,&g2,&A12,&B12,&C12,&A21,&B21,&C21);
                        //printf("%i %i %f %f\n",g1,g2,A12,A21);
                        if ((g1==uni->subgroup[i][1])&& (g2==uni->subgroup[j][1])){
                            uni->subgInt[i][j][0]=A12;
                            uni->subgInt[i][j][1]=B12;
                            uni->subgInt[i][j][2]=C12;
                            uni->subgInt[j][i][0]=A21;
                            uni->subgInt[j][i][1]=B21;
                            uni->subgInt[j][i][2]=C21;
                            //printf("%i %i %f %f %f %f %f %f\n",g1,g2,uni->subgInt[i][j][0],uni->subgInt[i][j][1],uni->subgInt[i][j][2],
                                    //uni->subgInt[j][i][0],uni->subgInt[j][i][1],uni->subgInt[j][i][2]);
                            break;
                        }
                    } while(g1<89);
                    rewind(f);
                }
            }
        }
    }
    fclose(f);

    for (i=0; i<numSubs;i++){
        uni->subsR[i]=0;
        uni->subsQ[i]=0;
        for (j=0;j<numSubgroups;j++){
            uni->subsR[i]=uni->subsR[i]+uni->subsSubg[i][j]*uni->subgData[j][0];
            uni->subsQ[i]=uni->subsQ[i]+uni->subsSubg[i][j]*uni->subgData[j][1];
        }
    }
    uni->numSubs=numSubs;
    uni->numSubg=numSubgroups;
}

//Calculates activity coefficients according to UNIFAC models, at given T and composition.
void CALLCONV FF_ActivityUNIFAC(FF_UnifacData *data, const double *T, const double x[], double lnGammaC[], double lnGammaR[], double *gE){
    int i,j,k;
    //Combinatorial part calculation
    double meanR=0, meanQ=0, meanL=0,theta[data->numSubs],phi[data->numSubs];
    if (data->model==FF_EntropicFV){
        for (i=0;i<data->numSubs;i++){
            meanR=meanR+data->FV[i]*x[i];
            meanQ=meanQ+data->subsQ[i]*x[i];
        }
        for (i=0;i<data->numSubs;i++){
            phi[i]=data->FV[i]/meanR;//Free Volume ratio
            theta[i]=data->subsQ[i]/meanQ;//Area fraction
        }
    }
    else{
        for (i=0;i<data->numSubs;i++){
            meanR=meanR+data->subsR[i]*x[i];
            meanQ=meanQ+data->subsQ[i]*x[i];
        }
        for (i=0;i<data->numSubs;i++){
            phi[i]=data->subsR[i]/meanR;//Volume fraction/x[i]
            theta[i]=data->subsQ[i]/meanQ;//Area fraction
        }
    }
    if ((data->model==FF_UNIFACStd)||(data->model==FF_UNIFACPSRK)){
            for (i=0;i<data->numSubs;i++){
                lnGammaC[i]=log(phi[i])+1-phi[i]-5*data->subsQ[i]*(log(phi[i]/theta[i])+1-phi[i]/theta[i]);
                //lnGammaC[i]=log(phi[i])+5*data->subsQ[i]*log(theta[i]/phi[i])+data->subsL[i]-phi[i]*meanL;//alternative from Fredenslund
                //printf("lnGammaC[%i]:%f\n",i,lnGammaC[i]);
            }
    }
    else if((data->model==FF_EntropicFV)){
        for (i=0;i<data->numSubs;i++){
            lnGammaC[i]=1-phi[i]+log(phi[i]);
            //printf("lnGammaC[%i]:%f\n",i,lnGammaC[i]);
        }
    }
    else if ((data->model==FF_UNIFACDort)||(data->model==FF_UNIFACNist)){
        double meanR2=0, phi2[data->numSubs];
        for (i=0;i<data->numSubs;i++) meanR2=meanR2+pow(data->subsR[i],0.75)*x[i];
        for (i=0;i<data->numSubs;i++){
            phi2[i]=pow(data->subsR[i],0.75)/meanR2;//Volume fraction /x[i]
            lnGammaC[i]=(1-phi2[i]+log(phi2[i])-5*data->subsQ[i]*(1-phi[i]/theta[i]+log(phi[i]/theta[i])));
            //printf("lnGammaC[%i]:%f\n",i,lnGammaC[i]);
        }
    }
    else if (data->model==FF_UNIFACZM){
        double meanR2=0, phi2[data->numSubs];
        for (i=0;i<data->numSubs;i++){
            if (data->subsR[i]>65) meanR2=meanR2+data->subsR[i]*0.6583*x[i];//we consider a polymer if R is over 65
            else meanR2=meanR2+data->subsR[i]*x[i];
        }
        for (i=0;i<data->numSubs;i++){
            if (data->subsR[i]>65) phi2[i]=data->subsR[i]*0.6583/meanR2;//Volume fraction /x[i]
            else phi2[i]=data->subsR[i]/meanR2;
            lnGammaC[i]=(1-phi2[i]+log(phi2[i])-5*data->subsQ[i]*(1-phi[i]/theta[i]+log(phi[i]/theta[i])));
            //printf("lnGammaC[%i]:%f\n",i,lnGammaC[i]);
        }
    }
    //Residual part calculation
    //first we calculate the interaction between subgroups, as a temperature function
    double psi[data->numSubg][data->numSubg], substSubgrTheta[data->numSubs][data->numSubg],substSubgrLambda[data->numSubs][data->numSubg],
            substSubgrSum[data->numSubs][data->numSubg];
    if (data->model==FF_UNIFACStd){
        for (i=0;i<data->numSubg;i++) for (j=0;j<data->numSubg;j++) psi[i][j]=exp(-data->subgInt[i][j][0]/ *T);
    }
    else if ((data->model==FF_UNIFACPSRK)||(data->model==FF_UNIFACDort)||(data->model==FF_UNIFACNist)||(data->model==FF_EntropicFV)||(data->model==FF_UNIFACZM)){
        for (i=0;i<data->numSubg;i++) for (j=0;j<data->numSubg;j++) psi[i][j]=exp(-data->subgInt[i][j][0]/ *T -
               data->subgInt[i][j][1] - data->subgInt[i][j][2]* *T);
    }
    //now the calculation of theta and lambda por each subgroup in each substance
    for (i=0;i<data->numSubs;i++){
        for (j=0;j<data->numSubg;j++) substSubgrTheta[i][j]=data->subsSubg[i][j]*data->subgData[j][1]/data->subsQ[i];
    }
    for (i=0;i<data->numSubs;i++) for (j=0;j<data->numSubg;j++){
        substSubgrSum[i][j]=0;
        for (k=0;k<data->numSubg;k++) substSubgrSum[i][j]=substSubgrSum[i][j]+substSubgrTheta[i][k]*psi[k][j];
    }
    for (i=0;i<data->numSubs;i++) for (j=0;j<data->numSubg;j++){
        substSubgrLambda[i][j]=0;
        if (data->subsSubg[i][j]>0){
            for (k=0;k<data->numSubg;k++) substSubgrLambda[i][j]=substSubgrLambda[i][j]+substSubgrTheta[i][k]*psi[k][j];
            substSubgrLambda[i][j]=1-log(substSubgrLambda[i][j]);
            for (k=0;k<data->numSubg;k++) substSubgrLambda[i][j]=substSubgrLambda[i][j]-substSubgrTheta[i][k]*psi[j][k]/substSubgrSum[i][k];
             substSubgrLambda[i][j]=data->subgData[j][1]* substSubgrLambda[i][j];
        }
    }
    //Now we calculate theta and lambda of the subgroups in the whole mixture
    double subgrTheta[data->numSubg],subgrSum[data->numSubg],subgrLambda[data->numSubg];
    for (j=0;j<data->numSubg;j++){
        subgrTheta[j]=0;
        for (i=0;i<data->numSubs;i++) subgrTheta[j]=subgrTheta[j]+data->subgData[j][1]*x[i]*data->subsSubg[i][j]/meanQ;
    }
    for (j=0;j<data->numSubg;j++){
        subgrSum[j]=0;
        for (k=0;k<data->numSubg;k++) subgrSum[j]=subgrSum[j]+subgrTheta[k]*psi[k][j];
    }
    for (j=0;j<data->numSubg;j++){
        subgrLambda[j]=0;
        for (k=0;k<data->numSubg;k++) subgrLambda[j]=subgrLambda[j]+ subgrTheta[k]*psi[k][j];
        subgrLambda[j]=1-log(subgrLambda[j]);
        for (k=0;k<data->numSubg;k++) subgrLambda[j]=subgrLambda[j]-subgrTheta[k]*psi[j][k]/subgrSum[k];
        subgrLambda[j]=data->subgData[j][1]*subgrLambda[j];
    }
    *gE=0;
    for (i=0;i<data->numSubs;i++){
        lnGammaR[i]=0;
        for (j=0;j<data->numSubg;j++) lnGammaR[i]=lnGammaR[i]+data->subsSubg[i][j]*(subgrLambda[j]-substSubgrLambda[i][j]);
        //printf("lnGammaR[%i]:%f\n",i,lnGammaR[i]);
        *gE=*gE+x[i]*(lnGammaC[i]+lnGammaR[i]);
    }
}


//Calculates fugacity and activity coefficients, at given T and composition, from an activity model
void CALLCONV FF_PhiFromActivity(const int *numSubs,enum FF_ActModel *model,const  FF_BaseProp baseProp[],const double pintParam[],
                                 const enum FF_IntParamForm *form,const bool *useVp,const enum FF_EosType *eosType,const void *data,
                                 const double *T,const double *P,const double x[],double gamma[],double phi[],double *gE){
    //First the calculation of the activity coefficient
    switch(*model){
    case FF_Wilson:
        FF_ActivityWilson(numSubs,baseProp,pintParam,form,T,x,gamma,gE);
        break;
    case FF_NRTL:
        FF_ActivityNRTL(numSubs,pintParam,form,T,x,gamma,gE);
        break;
    case FF_UNIQUAC:
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

