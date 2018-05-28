/*

 * FFeosMix.c
 *
 *  Created on: 14/07/2013
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

// contains EOS calculations for pure substances
//==============================================

#include <math.h>
#include <stdio.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFeosMix.h"


//Mixture Cubic EOS calculations
//==============================

//Calculates Theta,b,dTheta/dT, d2Theta/dT2, dTheta/dX[i] and db/dX[i] for a mixture, given a cubic EOS,a mixing rule, composition, and pure substance parameters
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamCubicEOS(const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
        const double pintParam[],const double x[], FF_CubicParam *param,double dTheta_dXi[],double db_dXi[],double dc_dXi[])
{

    int i,j;
    FF_CubicParam sParam[*numSubs];
    double Pr,dPr,d2Pr;//Will be theta[i]*theta[j], and their derivatives
    double Co,dCo,d2Co;//This will be (theta[i]*theta[j])^0.5 and their derivatives regarding T
    double k,dk,d2k;//Will be the interaction coefficient for a[i][j], and their derivatives regarding T
    double kInv,dkInv,d2kInv,aInv,bInv;//Will be the corresponding [j][i] values
    double kCo,dkCo,d2kCo;//will be the combined interaction in PR and MKP rules
    double l;//the interaction coefficient for b[i][j]in VdW and PR rules, and for a2[i][j] for MKP rule
    double a,da,d2a;//Is the combined Theta[i,j] used in VdW and PR for the calculation of mixture Theta. It is also the first part of MKP Theta calculation
    double b;


    //double sumdTheta_dx=0,sumdb_dx=0;
    //First we get the parameters of the individual substances
    //later, calculation of the interaction coefficient and of its derivatives
    for (i=0;i< *numSubs;i++)//First we get the parameters of the individual substances
    {
        FF_FixedParamCubic(&data[i],&sParam[i]);
        FF_ThetaDerivCubic(T,&data[i],&sParam[i]);
        //printf("i a,b,Theta,dTheta,d2Theta: %i %f %f %f %f %f\n",i,sParam[i].a,sParam[i].b,sParam[i].Theta*1e3,sParam[i].dTheta*1e3,sParam[i].d2Theta*1e3);
    }
    param->Theta=0;
    param->dTheta=0;
    param->d2Theta=0;
    param->b=0;
    param->c=0;
    param->u=sParam[0].u;
    param->w=sParam[0].w;
    for (i=0;i< *numSubs;i++){
        if (!((sParam[i].u==param->u)&&(sParam[i].w==param->w))) return;//We check that all cubic eos are of the same type
        param->c=param->c+x[i]*sParam[i].c;//Calculation of mixture volume translation
        dc_dXi[i]=sParam[i].c;
    }
    switch (*rule)//six values are passed by every substances pair. Not all of them are used in every mixing rule.
    {
        case FF_VdW://van der Waals quadratic mixing rule
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    Pr=sParam[i].Theta*sParam[j].Theta;
                    dPr=sParam[i].dTheta*sParam[j].Theta+sParam[i].Theta*sParam[j].dTheta;
                    d2Pr=sParam[i].d2Theta*sParam[j].Theta+2*sParam[i].dTheta*sParam[j].dTheta+sParam[i].Theta*sParam[j].d2Theta;
                    Co=fabs(pow(Pr,0.5));
                    dCo=0.5*pow(Pr,-0.5)*dPr;
                    d2Co=-0.25*pow(Pr,-1.5)*dPr*dPr+0.5*pow(Pr,-0.5)*d2Pr;
                    //printf("i,j k: %i %i %f\n",i,j,pintParam[(i* *numSubs+j)*6]);
                    //Three values are used for the calculation of the mixture parameter k[i,j]=k[j,i], the forth for l[i,j]=l[j,i], used in the calculation of b[i,j]
                    k=pintParam[(i* *numSubs+j)*6] + pintParam[(i* *numSubs+j)*6+1]* *T+pintParam[(i* *numSubs+j)*6+2]/ *T ;
                    dk=pintParam[(i* *numSubs+j)*6+1]-pintParam[(i* *numSubs+j)*6+2]/(*T * *T);
                    d2k=2* pintParam[(i* *numSubs+j)*6+2]/(*T * *T * *T);
                    //printf("i,j k,dk,d2k: %i %i %f %f %f\n",i,j,k,dk,d2k);
                    a=Co*(1-k);
                    da=dCo*(1-k)-Co*dk;
                    d2a=d2Co*(1-k)-2*dCo*dk-Co*d2k;
                    param->Theta=param->Theta+x[i]*x[j]*a;
                    param->dTheta=param->dTheta+x[i]*x[j]*da;
                    param->d2Theta=param->d2Theta+x[i]*x[j]*d2a;
                    //printf("a,da,d2a[i][j]: %i %i %f %f %f\n",i,j,a*1e3,da*1e3,d2a*1e3);
                    l=pintParam[(i* *numSubs+j)*6+3];
                    b=(sParam[i].b+sParam[j].b)/2*(1-l);//this is b[i,j]=(b[i]+b[j])/2*(1-l)
                    param->b=param->b+x[i]*x[j]*b;
                    //kInv=pintParam[(j* *numSubs+i)*6] + pintParam[(j* *numSubs+i)*6+1]* *T+pintParam[(j* *numSubs+i)*6+2]/ *T ;
                    //aInv=Co*(1-kInv);
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*2*a;//This is the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    //This allows for asimetric k[i,j] and , al
                    db_dXi[i]=db_dXi[i]+x[j]*2*b;
                }
            }
            break;


        case FF_PR://Panagiotopoulos and Reid composition dependent mixing rule
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    Pr=sParam[i].Theta*sParam[j].Theta;
                    dPr=sParam[i].dTheta*sParam[j].Theta+sParam[i].Theta*sParam[j].dTheta;
                    d2Pr=sParam[i].d2Theta*sParam[j].Theta+2*sParam[i].dTheta*sParam[j].dTheta+sParam[i].Theta*sParam[j].d2Theta;
                    Co=fabs(pow(Pr,0.5));
                    dCo=0.5*pow(Pr,-0.5)*dPr;
                    d2Co=-0.25*pow(Pr,-1.5)*dPr*dPr+0.5*pow(Pr,-0.5)*d2Pr;
                    k=pintParam[(i* *numSubs+j)*6] + pintParam[(i* *numSubs+j)*6+1]* *T+pintParam[(i* *numSubs+j)*6+2]/ *T ;
                    dk=pintParam[(i* *numSubs+j)*6+1]-pintParam[(i* *numSubs+j)*6+2]/(*T * *T);
                    d2k=2* pintParam[(i* *numSubs+j)*6+2]/(*T * *T * *T);
                    kInv=pintParam[(j* *numSubs+i)*6] + pintParam[(j* *numSubs+i)*6+1]* *T+pintParam[(j* *numSubs+i)*6+2]/ *T ;
                    dkInv=pintParam[(j* *numSubs+i)*6+1]-pintParam[(j* *numSubs+i)*6+2]/(*T * *T);
                    d2kInv=2* pintParam[(j* *numSubs+i)*6+2]/(*T * *T * *T);
                    //printf("i,j k,dk,d2k: %i %i %f %f %f\n",i,j,k,dk,d2k);
                    kCo=1-k+(k-kInv)*x[i];
                    dkCo=-dk+(dk-dkInv)*x[i];
                    d2kCo=-d2k+(d2k-d2kInv)*x[i];
                    a=Co*kCo;//this is theta[i,j]=(1-k[i,j]+x[i]*(k[i,j]-k[j,i]))*(theta[i]*theta[j])^0.5
                    da=dCo*kCo+Co*dkCo;//this is the derivarive of a[i][j] regarding T
                    d2a=d2Co*kCo+2*dCo*dkCo+Co*d2kCo;//this is the second derivarive of a[i][j] regarding T
                    param->Theta=param->Theta+x[i]*x[j]*a;
                    param->dTheta=param->dTheta+x[i]*x[j]*da;
                    param->d2Theta=param->d2Theta+x[i]*x[j]*d2a;
                    l=*(pintParam+(i* *numSubs+j)*6+3);
                    b=(1-l)*(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(1-l[i,j])*(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b;
                    aInv=Co*(1-kInv+(kInv-k)*x[j]);
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a+x[i]*Co*(k-kInv)+aInv);//This is the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    db_dXi[i]=db_dXi[i]+x[j]*2*b;
                }

            }
            for (i=0;i<*numSubs;i++)//And now db/dX[i] and dTheta/dX[i]
            {
                for (j=0;j<*numSubs;j++)
                {

                }
            }
            break;
/*
        case FF_MKP://Mathias, Klotz and Prausnitz composition dependent mixing rule
        {
            double a2[*numSubs][*numSubs];//It is the secon part of combined Theta[i,j] calculation in MKP rule
            double da2[*numSubs][*numSubs];//It is the secon part of combined Theta[i,j] calculation in MKP rule
            double a3[*numSubs];//Is the summatory of second part of Theta calculation for each substance, and it derivative regarding T
            double da3[*numSubs];//The same for second part in MKP
            double d2a3[*numSubs];//The same for second part in MKP
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                a3[i]=0;
                da3[i]=0;
                d2a3[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    //Three values are used for the calculation of the interaction parameters, one is k[i,j]=k[j,i], the second lambda[i,j]=-lambda[j,i], and l[i,j]=l[j,i]
                    a[i][j]=fabs((1-*(pintParam+(i* *numSubs+j)*6)+x[i]*(*(pintParam+(i* *numSubs+j)*6)-*(pintParam+(j* *numSubs+i)*6)))*pow(sParam[i].Theta*sParam[j].Theta,0.5));
                    //this is Theta1[i,j]=(1-k[i,j]+x[i]*(k[i,j]-k[j,i]))*(theta[i]*theta[j])^0.5. It is the first part of Theta calculation
                    da[i][j]=0.5*a[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);
                    param->Theta=param->Theta+x[i]*x[j]*a[i][j];//We increase theta in the part that is as in VdW mixing rule
                    param->dTheta=param->dTheta+x[i]*x[j]*da[i][j];
                    param->d2Theta=param->d2Theta+x[i]*x[j]*(da[i][j]*da[i][j]/a[i][j]+0.5*a[i][j]*((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                            sParam[i].Theta/sParam[i].Theta+(sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta));
                    b[i][j]=(1-*(pintParam+(i* *numSubs+j)*6+2))*(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(1-l[i,j])*(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b[i][j];
                    a2[i][j]= x[j]*pow(sParam[i].Theta*sParam[j].Theta,0.1666667)* pow(fabs(*(pintParam+(j* *numSubs+i)*6+1)),0.333333);
                    if (*(pintParam+(j* *numSubs+i)*6+1)<0) a2[i][j]=-a2[i][j];
                    a3[i]=a3[i]+a2[i][j];
                    //printf("%f ",*(pintParam+(j* *numSubs+i)*3+1));
                    //printf("a2[i,j]:%f\n",a2[i][j]);
                    da2[i][j]=a2[i][j]/6*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);
                    da3[i]=da3[i]+da2[i][j];
                    if (i!=i)
                    d2a3[i]=d2a3[i]+da2[i][j]*da2[i][j]/a2[i][j]+a2[i][j]*(((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                        sParam[i].Theta/sParam[i].Theta)+((sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta))/6;
                    //printf("d2a3[i]%f\n",d2a3[i]);
                }
                param->Theta=param->Theta+x[i]*pow(a3[i],3);//We increase Theta and its derivatives in the part special for MKP
                param->dTheta=param->dTheta+x[i]*3*pow(a3[i],2)*da3[i];
                param->d2Theta=param->d2Theta+x[i]*3*(2*a3[i]*da3[i]*da3[i]+a3[i]*a3[i]*d2a3[i]);
            }
            for (i=0;i<*numSubs;i++)//And now db/dX[i] and dTheta/dX[i]
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a[i][j]+a[j][i]);//This is the first part (VdW) of the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    //This allows for asimetric k[i,j]
                    dTheta_dXi[i]=dTheta_dXi[i]+3*x[j]*a3[j]*a3[j]*a2[j][i]/x[i];//This is the derivative for the second part when [i] is not the first component of the pair
                    db_dXi[i]=db_dXi[i]+x[j]*(b[i][j]+b[j][i]);
                }
                dTheta_dXi[i]=dTheta_dXi[i]+pow(a3[i],3);//this is the derivative for the second part when [i] is the first component of the pair
            }
            break;
        }
        default:
            param->Theta=0;
            param->b=0;
            param->c=0;
            break;
            */
    }
    //printf("Theta:%f dTheta:%f d2Theta:%f b:%f\n",param->Theta,param->dTheta,param->d2Theta,param->b);
}

//Calculates Theta,b,dTheta/dT, d2Theta/dT2, dTheta/dX[i] and db/dX[i] for a mixture, given a cubic EOS,a mixing rule, composition, and pure substance parameters
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixParamCubicEOSOld(const enum FF_EOS eos[],const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
        const double pintParam[],const double x[], FF_CubicParam *param,double dTheta_dXi[],double db_dXi[],double dc_dXi[])
{
    int i,j;
     FF_CubicParam sParam[*numSubs];
    double rootTij;//This will be (theta[i]*theta[j])^0.5
    double dRootTij;//Will be the derivative of (theta[i]*theta[j])^0.5 regarding T
    double k[*numSubs][*numSubs];//the interaction coefficient for a[i][j]
    double dk[*numSubs][*numSubs];//Derivative of the interaction coefficient for a[i][j] regarding T
    double d2k[*numSubs][*numSubs];//Second derivative of the interaction coefficient for a[i][j] regarding T
    double l[*numSubs][*numSubs];//the interaction coefficient for b[i][j]in VdW and PR rules, and for b[i][j]
    double dl[*numSubs][*numSubs];//Derivative of the interaction coefficient for b[i][j]in VdW and PR rules, and of a2[i][j] for MKP rule
    double a[*numSubs][*numSubs];//Is the combined Theta[i,j] used in VdW and PR for the calculation of mixture Theta. It is also the first part of MKP Theta calculation
    double da[*numSubs][*numSubs];//Is the derivative of Theta[i,j] regarding T
    double b[*numSubs][*numSubs];

    //double sumdTheta_dx=0,sumdb_dx=0;
    //First we get the parameters of the individual substances
    //later, calculation of the interaction coefficient and of its derivatives
    for (i=0;i< *numSubs;i++)//First we get the parameters of the individual substances
    {
        FF_FixedParamCubic(&data[i],&sParam[i]);
        FF_ThetaDerivCubic(T,&data[i],&sParam[i]);
        for (j=0;j<*numSubs;j++){
            k[i][j]=pintParam[(i* *numSubs+j)*6] + pintParam[(i* *numSubs+j)*6+1]* *T+pintParam[(i* *numSubs+j)*6+2]/ *T ;
            dk[i][j]=pintParam[(i* *numSubs+j)*6+1]-pintParam[(i* *numSubs+j)*6+2]/(*T * *T);
            d2k[i][j]=2* pintParam[(i* *numSubs+j)*6+2]/(*T * *T * *T);
            l[i][j]=pintParam[(i* *numSubs+j)*6+3] + pintParam[(i* *numSubs+j)*6+4]* *T+pintParam[(i* *numSubs+j)*6+5]/ *T;
            dl[i][j]=pintParam[(i* *numSubs+j)*6+4]-pintParam[(i* *numSubs+j)*6+5]/(*T * *T);
        }
    }
    param->Theta=0;
    param->dTheta=0;
    param->d2Theta=0;
    param->b=0;
    param->c=0;
    param->u=sParam[0].u;
    param->w=sParam[0].w;
    for (i=0;i< *numSubs;i++){
        param->c=param->c+x[i]*sParam[i].c;//Calculation of mixture volume translation
        dc_dXi[i]=sParam[i].c;
    }
    switch (*rule)//six values are passed by every substances pair. Not all of them are used in every mixing rule.
    {
        case FF_VdW://van der Waals quadratic mixing rule
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                for (j=0;j<*numSubs;j++)
                {
                    //Two values are used for the calculation of the mixture parameters, one for k[i,j]=k[j,i], the other for l[i,j]=l[j,i], used in the calculation of b[i,j]
                    rootTij=fabs(pow(sParam[i].Theta*sParam[j].Theta,0.5));
                    dRootTij=0.5*rootTij*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);
                    a[i][j]=(1-k[i][j])*rootTij;//this is theta[i,j]=(1-k[i,j])*(theta[i]*theta[j])^0.5
                    param->Theta=param->Theta+x[i]*x[j]*a[i][j];
                    //da[i][j]=dRootTij*(1-k[i][j])-rootTij*dk[i][j];//this is the derivarive of a[i][j] regarding T
                    da[i][j]=0.5*a[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta)-dk[i][j]*rootTij;
                    param->dTheta=param->dTheta+x[i]*x[j]*da[i][j];
                    param->d2Theta=param->d2Theta+x[i]*x[j]*(0.5*da[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta)+0.5*a[i][j]*
                                   ((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/sParam[i].Theta/sParam[i].Theta+(sParam[j].d2Theta*
                                   sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta)-d2k[i][j]*rootTij-dk[i][j]*dRootTij);
                    //param->d2Theta=param->d2Theta+x[i]*x[j]*(da[i][j]*da[i][j]/a[i][j]+0.5*a[i][j]*((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                                //sParam[i].Theta/sParam[i].Theta+(sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta));
                    b[i][j]=(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b[i][j];
                }
            }
            for (i=0;i<*numSubs;i++)//And now db/dX[i] and dTheta/dX[i]
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a[i][j]+a[j][i]);//This is the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    //This allows for asimetric k[i,j] and , al
                    db_dXi[i]=db_dXi[i]+x[j]*(b[i][j]+b[j][i]);
                }
            }
            break;
        case FF_PR://Panagiotopoulos and Reid composition dependent mixing rule
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                for (j=0;j<*numSubs;j++)
                {
                    //Three values are used for the calculation of the interaction parameters, one is k[i,j] not equal to k[j,i], the other is l[i,j]=l[j,i]
                    a[i][j]=fabs((1-*(pintParam+(i* *numSubs+j)*6)+x[i]*(*(pintParam+(i* *numSubs+j)*6)-*(pintParam+(j* *numSubs+i)*6)))*pow(sParam[i].Theta*sParam[j].Theta,0.5));
                    //this is theta[i,j]=(1-k[i,j]+x[i]*(k[i,j]-k[j,i]))*(theta[i]*theta[j])^0.5
                    param->Theta=param->Theta+x[i]*x[j]*a[i][j];
                    da[i][j]=0.5*a[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);//this is the derivarive of a[i][j] regarding T
                    param->dTheta=param->dTheta+x[i]*x[j]*da[i][j];
                    param->d2Theta=param->d2Theta+x[i]*x[j]*(da[i][j]*da[i][j]/a[i][j]+0.5*a[i][j]*((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                                sParam[i].Theta/sParam[i].Theta)+((sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta));
                    b[i][j]=(1-*(pintParam+(i* *numSubs+j)*6+3))*(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(1-l[i,j])*(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b[i][j];
                }
            }
            for (i=0;i<*numSubs;i++)//And now db/dX[i] and dTheta/dX[i]
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a[i][j]+pow(sParam[i].Theta*sParam[j].Theta,0.5)*x[i]*(*(pintParam+(i* *numSubs+j)*6)-*(pintParam+(j* *numSubs+i)*6))+
                                                      a[j][i]);//This is the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    db_dXi[i]=db_dXi[i]+x[j]*(b[i][j]+b[j][i]);
                }
            }
            break;
        case FF_MKP://Mathias, Klotz and Prausnitz composition dependent mixing rule
        {
            double a2[*numSubs][*numSubs];//It is the secon part of combined Theta[i,j] calculation in MKP rule
            double da2[*numSubs][*numSubs];//It is the secon part of combined Theta[i,j] calculation in MKP rule
            double a3[*numSubs];//Is the summatory of second part of Theta calculation for each substance, and it derivative regarding T
            double da3[*numSubs];//The same for second part in MKP
            double d2a3[*numSubs];//The same for second part in MKP
            for (i=0;i<*numSubs;i++)//We calculate b, theta, dTheta/dT, and d2Theta/dT for the mix
            {
                a3[i]=0;
                da3[i]=0;
                d2a3[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    //Three values are used for the calculation of the interaction parameters, one is k[i,j]=k[j,i], the second lambda[i,j]=-lambda[j,i], and l[i,j]=l[j,i]
                    a[i][j]=fabs((1-*(pintParam+(i* *numSubs+j)*6)+x[i]*(*(pintParam+(i* *numSubs+j)*6)-*(pintParam+(j* *numSubs+i)*6)))*pow(sParam[i].Theta*sParam[j].Theta,0.5));
                    //this is Theta1[i,j]=(1-k[i,j]+x[i]*(k[i,j]-k[j,i]))*(theta[i]*theta[j])^0.5. It is the first part of Theta calculation
                    da[i][j]=0.5*a[i][j]*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);
                    param->Theta=param->Theta+x[i]*x[j]*a[i][j];//We increase theta in the part that is as in VdW mixing rule
                    param->dTheta=param->dTheta+x[i]*x[j]*da[i][j];
                    param->d2Theta=param->d2Theta+x[i]*x[j]*(da[i][j]*da[i][j]/a[i][j]+0.5*a[i][j]*((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                            sParam[i].Theta/sParam[i].Theta+(sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta));
                    b[i][j]=(1-*(pintParam+(i* *numSubs+j)*6+2))*(sParam[i].b+sParam[j].b)/2;//this is b[i,j]=(1-l[i,j])*(b[i]+b[j])/2
                    param->b=param->b+x[i]*x[j]*b[i][j];
                    a2[i][j]= x[j]*pow(sParam[i].Theta*sParam[j].Theta,0.1666667)* pow(fabs(*(pintParam+(j* *numSubs+i)*6+1)),0.333333);
                    if (*(pintParam+(j* *numSubs+i)*6+1)<0) a2[i][j]=-a2[i][j];
                    a3[i]=a3[i]+a2[i][j];
                    //printf("%f ",*(pintParam+(j* *numSubs+i)*3+1));
                    //printf("a2[i,j]:%f\n",a2[i][j]);
                    da2[i][j]=a2[i][j]/6*(sParam[i].dTheta/sParam[i].Theta+sParam[j].dTheta/sParam[j].Theta);
                    da3[i]=da3[i]+da2[i][j];
                    if (i!=i)
                    d2a3[i]=d2a3[i]+da2[i][j]*da2[i][j]/a2[i][j]+a2[i][j]*(((sParam[i].d2Theta*sParam[i].Theta-sParam[i].dTheta*sParam[i].dTheta)/
                        sParam[i].Theta/sParam[i].Theta)+((sParam[j].d2Theta*sParam[j].Theta-sParam[j].dTheta*sParam[j].dTheta)/sParam[j].Theta/sParam[j].Theta))/6;
                    //printf("d2a3[i]%f\n",d2a3[i]);
                }
                param->Theta=param->Theta+x[i]*pow(a3[i],3);//We increase Theta and its derivatives in the part special for MKP
                param->dTheta=param->dTheta+x[i]*3*pow(a3[i],2)*da3[i];
                param->d2Theta=param->d2Theta+x[i]*3*(2*a3[i]*da3[i]*da3[i]+a3[i]*a3[i]*d2a3[i]);
            }
            for (i=0;i<*numSubs;i++)//And now db/dX[i] and dTheta/dX[i]
            {
                dTheta_dXi[i]=0;
                db_dXi[i]=0;
                for (j=0;j<*numSubs;j++)
                {
                    dTheta_dXi[i]=dTheta_dXi[i]+x[j]*(a[i][j]+a[j][i]);//This is the first part (VdW) of the derivative of Theta regarding x[i]. Due to the fact that there is x[i]*x[j] + x[j]*x[i]
                    //This allows for asimetric k[i,j]
                    dTheta_dXi[i]=dTheta_dXi[i]+3*x[j]*a3[j]*a3[j]*a2[j][i]/x[i];//This is the derivative for the second part when [i] is not the first component of the pair
                    db_dXi[i]=db_dXi[i]+x[j]*(b[i][j]+b[j][i]);
                }
                dTheta_dXi[i]=dTheta_dXi[i]+pow(a3[i],3);//this is the derivative for the second part when [i] is the first component of the pair
            }
            break;
        }
        default:
            param->Theta=0;
            param->b=0;
            param->c=0;
            break;
    }
    //printf("Theta:%f dTheta:%f d2Theta:%f b:%f\n",param->Theta,param->dTheta,param->d2Theta,param->b);
}


/*
//Calculates Theta,b,delta and epsilon for a mixture, given a cubic EOS,a gE mixing rule, composition,pure substance parameters, and gE
void CALLCONV calcMixParamCubicEOSgE(const EOS& eos,const FF_MixingRule& rule,const double& T,const int& *numSubs,const FF_CubicParam sParam[],
        const double& gE,const double x[],FF_CubicParam& param,double dNb_dNi[])
{
    double a[numSubs][numSubs];
    double b[numSubs][numSubs];
    double sPart=0;//sum of individual substances contribution
    double q1=0,q2=0,lambda=0;//Parameters for  alpha (and Theta) calculation, lambda is for LCVM mixing rule
    param->Theta=0;
    param->b=0;
    param->c=0;
    int i;
    for (i=0;i<numSubs;i++) param->c=param->c+x[i]*sParam[i].c;//Calculation of mixture volume translation
    switch (eos)
    {
    case FF_PR76:
    case FF_PR78://Peng-Robinson.
    case FF_PRSV1://Peng-Robinson, Stryjek-Vera.
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        break;
    case FF_SRK://Soaves-Redlich-Kwong.
        param->u=1;
        param->w=0;
        break;
    }
    switch (rule)//Charge of fixed parameters for each mixing rule
    {
    case MHV1:
        if (eos==FF_PR78 or eos==FF_PRSV1) q1=-0.53;
        if (eos==FF_SRK) q1=-0.593;
        break;
    case PSRK:
        if (eos==FF_SRK) q1=-0.64663;
        break;
    case HVOS:
        if (eos==FF_PR78 or eos==FF_PRSV1) q1=-0.62323;
        break;
    case LCVM:
        if (eos==FF_PR78 or eos==FF_PRSV1)
        {
            q1=-0.52;//Am. Michelsen(zero pressure) part coefficient
            q2=-0.623;///Av. Huron-Vidal(infinite pressure) part coefficient
        }
        if (eos==FF_SRK)
        {
            q1=-0.593;//Am. Michelsen(zero pressure) part coefficient
            q2=-0.693;///Av. Huron-Vidal(infinite pressure) part coefficient
        }
        lambda=0.36;
        break;
    case MHV2:
        if (eos==FF_PR78 or eos==FF_PRSV1)
        {
            q1=-0.4347;
            q2=-0.003654;
        }
        if (eos==FF_SRK)
        {
            q1=-0.4783;
            q2=-0.0047;
        }
        break;
    }
    switch (rule)//Calculation of b and Theta, for each mixing rule
    {
    case MHV1:
    case PSRK:
    case HVOS:
        if (not(q1==0))
        {
            for (i=0;i<numSubs;i++)	param->b=param->b+x[i]*sParam[i].b;
            for (i=0;i<numSubs;i++) dNb_dNi[i]=sParam[i].b;
            for (i=0;i<numSubs;i++) sPart=sPart+x[i]*(R*T/q1*log(param->b/sParam[i].b)+sParam[i].Theta/sParam[i].b);
            param->Theta=(sPart+gE/q1)*param->b;//q1*alpha=sum(x[i]*(ln(b/b[i])+q1*alpha[i]))+gE/(R·T);With alpha=a/(bRT)
        }
        break;
    case LCVM://Good for highly asymmetric mixtures
        if (not(q1==0))
        {
            for (i=0;i<numSubs;i++)	param->b=param->b+x[i]*sParam[i].b;
            for (i=0;i<numSubs;i++) dNb_dNi[i]=sParam[i].b;
            for (i=0;i<numSubs;i++) param->Theta=param->Theta+x[i]*(sParam[i].Theta/sParam[i].b+(1-lambda)/q1*R*T*log(param->b/sParam[i].b));
        param->Theta=(param->Theta+(lambda/q2+(1-lambda)/q1)*gE)*param->b;
        }
        else
        {
            param->b=0;
            param->Theta=0;
        }
        break;
    case MHV2://Excellent for polar compounds and difficult mixtures
        if (not(q2==0))
        {
            for (i=0;i<numSubs;i++)	param->b=param->b+x[i]*sParam[i].b;
            for (i=0;i<numSubs;i++) dNb_dNi[i]=sParam[i].b;
            for (i=0;i<numSubs;i++) sPart=sPart+x[i]*(log(param->b/sParam[i].b)+(q1*sParam[i].Theta/sParam[i].b/R/T+
                        q2*pow(sParam[i].Theta/sParam[i].b/R/T,2)));
            param->Theta=param->b*R*T*(-q1-pow(pow(q1,2)-4*q2*(-sPart-gE/R/T),0.5))/2/q2;
            //q2*alpha^2+q1*alpha=sum(x[i]*(ln(b/b[i])+q2*alpha[i]^2+q1*alpha[i]))+gE/(R·T);With alpha=a/(bRT)
        }
        break;
    default:
        param->c=0;
        param->b=0;
        param->Theta=0;
        break;
    }
}
*/


//Mixture SAFT EOS calculations
//==============================

//Z and Arr calculation for a mixture, given T and V, according to FF_PCSAFT EOS
//------------------------------------------------------------------------------
CALLCONV FF_MixArrZfromTVPCSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                                        const  FF_SaftEOSdata data[],const double pintParam[],const double x[],double *Arr,double *Z)
{
    //Initial calculations: molecular volume, molecular density, and hard molecule volume fraction, etc.
    int i,j,k;
    double Vmolecular,rho,mM=0.0,substRho[*numSubs],d[*numSubs],pairSigma[*numSubs][*numSubs],pairEpsilon[*numSubs][*numSubs],
            pairD[*numSubs][*numSubs],pairKAB[*numSubs][*numSubs],pairEpsilonAB[*numSubs][*numSubs];
    Vmolecular = *V / Av * 1E+30;
    rho = 1 / Vmolecular;
    for (i=0;i<*numSubs;i++)
    {
        mM=mM+x[i]*data[i].m;
        substRho[i]=x[i]*rho;
        d[i] = data[i].sigma * (1 - 0.12 * exp(-3 * data[i].epsilon / *T));//Diameter of each segment of the chain
    }
    int combRul=2;
    for (i=0;i<*numSubs;i++)
        for (j=0;j<*numSubs;j++)
        {
            pairSigma[i][j]=(data[i].sigma+data[j].sigma)/2;
//			pairEpsilon[i][j]=pow((data[i].epsilon*data[j].epsilon),0.5);
            pairEpsilon[i][j]=(1-*(pintParam+(i**numSubs+j)*2)-*(pintParam+(i**numSubs+j)*2+1)* *T)*pow((data[i].epsilon*data[j].epsilon),0.5);
            pairD[i][j]=(d[i]+d[j])/2;
            switch (combRul)
            {
            case 1://CR-1 combining rule: (epsilonAB+epsilonAB)/2, kAB=(kAB*kAB)^0.5
                if (data[i].kAB==0) pairKAB[i][j]=data[j].kAB;//if only one substance is associating, we use its kAB
                else if (data[j].kAB==0) pairKAB[i][j]=data[i].kAB;
                else pairKAB[i][j]=pow((data[i].kAB*data[j].kAB),0.5);
                pairEpsilonAB[i][j]=(data[i].epsilonAB+data[j].epsilonAB)/2;
                break;
            case 2://CR-1 modified by Wolbach and Sandler
                if (data[i].kAB==0) pairKAB[i][j]=data[j].kAB;//if only one substance is associating, we use its kAB
                else if (data[j].kAB==0) pairKAB[i][j]=data[i].kAB;
                else pairKAB[i][j]=pow((data[i].kAB*data[j].kAB),0.5)*pow(2*pow(data[i].sigma*data[j].sigma,0.5)/(data[i].sigma+data[j].sigma),3);
                pairEpsilonAB[i][j]=(data[i].epsilonAB+data[j].epsilonAB)/2;
                break;
            case 3://combRul (epsilon*epsilon)**0.5, kAB=(kAB+kAB)/2. Not used
                pairKAB[i][j]=(data[i].kAB+data[j].kAB)/2;
                pairEpsilonAB[i][j]=pow((data[i].epsilonAB*data[j].epsilonAB),0.5);
                break;
            case 4://combRul (epsilon*epsilon)**0.5 , kAB=(...)**3. Not used
                pairKAB[i][j]=pow((pow(data[i].kAB,1/3)+pow(data[j].kAB,1/3))/2,3);
                pairEpsilonAB[i][j]=pow((data[i].epsilonAB*data[j].epsilonAB),0.5);
                break;
            }
        }
    //Contribution by hard sphere
    double dseta[4],auxDseta,Zhs,Ahs;
        for (i=0;i<4;i++)
        {
            auxDseta=0.0;
            for (j=0;j<*numSubs;j++)
                auxDseta=auxDseta+x[j]*data[j].m*pow(d[j],i)/6;
            dseta[i]=Pi*rho*auxDseta;
        }
        Zhs=mM*(1/(1-dseta[3])+3*dseta[1]*dseta[2]/dseta[0]/pow((1-dseta[3]),2)+(3*pow(dseta[2],3)-dseta[3]*pow(dseta[2],3))/dseta[0]/pow((1-dseta[3]),3)-1);
        Ahs=mM*1/dseta[0]*(3*dseta[1]*dseta[2]/(1-dseta[3])+pow(dseta[2],3)/dseta[3]/pow((1-dseta[3]),2)+(pow(dseta[2],3)/pow(dseta[3],2)-dseta[0])*log(1-dseta[3]));
    //Contribution by chain
    double pairGhs[*numSubs][*numSubs],pairDelta[*numSubs][*numSubs],dLghs_drho[*numSubs],Achain=0.0,Zchain=0.0,etaM,dLghsM_drho;
    for (i=0;i<*numSubs;i++)
    {
        for (j=0;j<*numSubs;j++)
        {
            pairGhs[i][j]=1/(1-dseta[3])+(d[i]*d[j]/(d[i]+d[j]))*3*dseta[2]/pow((1-dseta[3]),2)+
            pow((d[i]*d[j]/(d[i]+d[j])),2)*2*pow(dseta[2],2)/pow((1-dseta[3]),3);
            pairDelta[i][j]=pow(pairD[i][j],3)*pairGhs[i][j]*pairKAB[i][j]*(exp(pairEpsilonAB[i][j]/ *T)-1);//This will be used in contribution by molecular association
        }
        //Now we calculate d(ln(ghs(di))/d(rho) for each substance
        dLghs_drho[i] =(dseta[3]/pow((1-dseta[3]),2)+3*d[i]*dseta[2]/2/pow((1-dseta[3]),2)+3*d[i]*dseta[2]*dseta[3]/pow((1-dseta[3]),3)+
            pow(d[i],2)*pow(dseta[2],2)/pow((1-dseta[3]),3)+3*pow(d[i],2)*pow(dseta[2],2)*dseta[3]/2/pow((1-dseta[3]),4))*Vmolecular/pairGhs[i][i];
        Achain = Achain + x[i]*(1-data[i].m)* log(pairGhs[i][i]);
        Zchain= Zchain + x[i]*(1-data[i].m)* dLghs_drho[i];
    }
    Zchain= rho*Zchain;
    etaM=dseta[3];
    dLghsM_drho=(5*etaM/2-(pow(etaM,2)))/(1-etaM)/(1-0.5*etaM)*Vmolecular; //Derivative of ghs(mean diameter). Used later in association

    //contribution by dispersion
    double I[4]={0.0,0.0,0.0,0.0};
    FF_calcI1I2(mM,etaM,I);
    double sum1=0.0,sum2=0.0,Z1,C1,C2,Z2,Zdisp,Adisp;
    for (i=0;i<*numSubs;i++)
        for (j=0;j<*numSubs;j++)
        {	sum1=sum1+x[i]*x[j]*data[i].m*data[j].m*pairEpsilon[i][j]*pow(pairSigma[i][j],3);
            sum2=sum2+x[i]*x[j]*data[i].m*data[j].m*pow(pairEpsilon[i][j],2)*pow(pairSigma[i][j],3);
        }
    sum1=sum1/ *T;
    sum2=sum2/pow(*T,2);
    Z1 = -2 * Pi / Vmolecular * I[1] * sum1;
    C1 = 1/(1 + mM * (8 * etaM - 2 * pow(etaM,2)) / pow((1 - etaM),4) + (1 - mM) *
            (20 * etaM - 27 * pow(etaM,2)+ 12 * pow(etaM,3) - 2 * pow(etaM,4)) / pow(((1 - etaM) * (2 - etaM)),2));
    C2 = C1 * (mM * (-4 * pow(etaM,2) + 20 * etaM + 8) / pow((1 - etaM),5) + (1 - mM) * (2 * pow(etaM,3)+
            12 * pow(etaM,2) - 48 * etaM + 40) / pow(((1 - etaM) * (2 - etaM)),3));
    Z2 = -Pi / Vmolecular * mM * C1 * (I[3] - C2 * etaM * I[2])* sum2;
    Zdisp = Z1 + Z2;
    Adisp = - Pi * rho *(2* I[0] * sum1 +  mM * C1 * I[2] * sum2);

    //contribution by molecular association
    double xPos[*numSubs],xNeg[*numSubs],xAcid[*numSubs],Aassoc=0.0,Zassoc=0.0;//Fraction of non associated sites in each molecule
    for (i=0;i<*numSubs;i++)
    {
        if (data[i].nPos>0) xPos[i]=xNeg[i]=0.5; else xPos[i]=xNeg[i]=1.0;
        if (data[i].nAcid>0) xAcid[i]=0.5; else xAcid[i]=1.0;
    }
    for (k=0;k<13;k++) //This is a iteration to approximate non associated fraction
        for (i=0;i<*numSubs;i++)
        {
            if (data[i].nPos>0)
            {
                xPos[i]=0.0;
                for (j=0;j<*numSubs;j++)
                    xPos[i]=xPos[i]+substRho[j]*(data[j].nNeg*xNeg[j]+data[j].nAcid*xAcid[j])*pairDelta[i][j];
                xPos[i]=1/(1+xPos[i]);
            }
            if (data[i].nNeg>0)
            {
                xNeg[i]=0.0;
                for (j=0;j<*numSubs;j++)
                    xNeg[i]=xNeg[i]+substRho[j]*(data[j].nPos*xPos[j]+data[j].nAcid*xAcid[j])*pairDelta[i][j];
                xNeg[i]=1/(1+xNeg[i]);
            }
            if (data[i].nAcid>0)
            {
                xAcid[i]=0.0;
                for (j=0;j<*numSubs;j++)
                    xAcid[i]=xAcid[i]+substRho[j]*(data[j].nNeg*xNeg[j]+data[j].nPos*xPos[j])*pairDelta[i][j];
                xAcid[i]=1/(1+xAcid[i]);

            }
        }
    for (i=0;i<*numSubs;i++)
    {
        Aassoc=Aassoc+x[i]*(((data[i].nPos+data[i].nNeg+data[i].nAcid)/2)+data[i].nPos*(log(xPos[i])-xPos[i]/2)+data[i].nNeg*(log(xNeg[i])-xNeg[i]/2)+
                            data[i].nAcid*(log(xAcid[i])-xAcid[i]/2));
        Zassoc=Zassoc+x[i]*(data[i].nPos*(1-xPos[i])+data[i].nNeg*(1-xNeg[i])+data[i].nAcid*(1-xAcid[i]));
    }
        Zassoc=-0.5*(1+rho*dLghsM_drho)*Zassoc;

    *Arr=Ahs+Achain+Adisp+Aassoc;//Arr
    *Z=1+Zhs+Zchain+Zdisp+Zassoc;//Z
}

//Mixture P calculation given T and V, according to FF_PCSAFT EOS
//---------------------------------------------------------------
void CALLCONV FF_MixPfromTVPCSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                                        const  FF_SaftEOSdata data[],const double pintParam[],const double x[],double *P)
{
    double Arr,Z;
    FF_MixArrZfromTVPCSAFT(rule,T,V,numSubs,data,pintParam,x,&Arr,&Z);
    //printf("T:%f V:%f Z:%f\n",*T,*V,Z);
    *P=Z*R* *T/ *V;
}


//Mixture V,Arr and Z  calculation, given T and P and composition, according to FF_PCSAFT EOS
//-------------------------------------------------------------------------------------------
void CALLCONV FF_MixVfromTPPCSAFT(const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const  FF_SaftEOSdata data[],
                              const double pintParam[],const double x[],const char *option,double resultL[3],double resultG[3],char *state)
{
    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    double V,Arr,Z,Vplus,ArrPlus,Zplus,error,errorRel,dP_dV;
    double maxError=0.000001;
    int i;
    double dM=0,mM=0;
    double eta,Vmolecular;
    for (i=0;i<*numSubs;i++)
    {
        dM = dM + x[i]*data[i].sigma * (1 - 0.12 * exp(-3 * data[i].epsilon / *T));
        mM = mM + x[i]*data[i].m;
    }
    eta = 0.6;  //As initial guess, we suppose the volume fraction occupied by the hard molecule is 0.4
    Vmolecular = mM * (Pi * pow(dM,3) / 6) / eta * Av / 1E+30; //This will be the initial guess for the molecular volume in m3/mol
    if (*option!='g'){//we calculate the liquid phase if not only the gas phase has been asked for
        V = Vmolecular; //initial guess for mole volume in m3
            //Vplus=V; //Vplus will be a volume that produces a P lower than that wanted
        Vplus=V*1.0000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_MixArrZfromTVPCSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
        FF_MixArrZfromTVPCSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
        dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        while (dP_dV>=0){//Initial loop to find a negative derivative of P regarding V
            //printf("finding a liquid positive derivative. V:%f\n",V);
            V = V*1.5;
            Vplus=V*1.0000001;
            FF_MixArrZfromTVPCSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
            FF_MixArrZfromTVPCSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        }
        error =*P-R* *T*Z/V;
        errorRel=error/ *P;
        i=1;
        //printf("Initial liquid:T:%f V:%f dP/dV:%f err:%f\n",*T-273.15,V,dP_dV,error);
        while ((fabs(errorRel)>maxError)&&(dP_dV <0)&&(V>0)&&(i<51))
        {
            V=V+error/dP_dV;//Newton method for root finding
            if ((V<=0)||(dP_dV>0))//We slow the Newton method if V is negative, or dP/dV is positive
            {
               V=V-0.9*error/dP_dV;
               //printf("Hola V negativo o derivada positiva");
            }

            Vplus=V*1.0000001;//we calculate dP/dV numerically
            FF_MixArrZfromTVPCSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
            FF_MixArrZfromTVPCSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
            error =*P-R* *T*Z/V;
            errorRel=error/ *P;
            //printf("i:%d Vl:%f Vplus:%f Z:%f Zplus:%fdP_dV:%f error:%f\n",i,V,Vplus,Z,Zplus,dP_dV,error);
            i++;
        }
        if ((fabs(errorRel)<maxError)){
            resultL[0]=V;
            resultL[1]=Arr;
            resultL[2]=Z;
            *state='l';//We inform that we have done the calculation from the liquid end
        }
        else resultL[0]=resultL[1]=resultL[2]=0;
    }
    if (*option!='l')//and the gas phase if not only the liquid one has been asked for
    {
        V = R * *T / *P + Vmolecular;//initial guess for gas volume
        Vplus=V*1.000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_MixArrZfromTVPCSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
        FF_MixArrZfromTVPCSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
        error =*P-R* *T*Z/V;
        errorRel=error/ *P;
        dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        i=1;
        //printf("%d V:%f dP/dV:%f error:%f\n",i,V,dP_dV,error);
        while ((fabs(errorRel)>maxError)&&(dP_dV <0)&&(V>0)&&(i<51))
        {
            V=V+error/dP_dV;//Newton method for root finding
            if ((V<=0)||(dP_dV>0))//We slow the Newton method if V is negative, or dP/dV is positive
            {
               V=V-0.9*error/dP_dV;//Newton method for root finding, slowed
            }

            Vplus=V*1.000001;
            FF_MixArrZfromTVPCSAFT(rule,T,&V,numSubs,data,pintParam,x,&Arr,&Z);
            FF_MixArrZfromTVPCSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrPlus,&Zplus);
            error =*P-R* *T*Z/V;
            errorRel=error/ *P;
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus- V);
            i=i+1;
            //printf("i:%d V:%f dP_dV:%f error:%f\n",i,V,dP_dV,error);
        }
        if (fabs(errorRel)<maxError){
            resultG[0]=V;
            resultG[1]=Arr;
            resultG[2]=Z;
            if (*state=='l') *state='b';
            else *state='g';
            //printf("hola\n");
        }
        else resultG[0]=resultG[1]=resultG[2]=0;
    }
}


//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a mixture, given T and V, according to FF_PCSAFT EOS
//--------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_MixArrDerPCSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                              const  FF_SaftEOSdata data[],const double pintParam[],const double x[],double result[6])
{
    double dV=*V * 0.00001;//increments of V and T used to obtain dArr/dV,dArr/dT and d2Arr/dT2 in SAFT eos
    double Vplus=*V + dV;
    double dT=0.01;
    double Tplus=*T+dT;
    double Tminus=*T-dT;
    double Arr,Z,ArrVplus,ZVplus,ArrTplus,ZTplus,ArrTminus,ZTminus;
    FF_MixArrZfromTVPCSAFT(rule,T,V,numSubs,data,pintParam,x,&Arr,&Z);
    FF_MixArrZfromTVPCSAFT(rule,T,&Vplus,numSubs,data,pintParam,x,&ArrVplus,&ZVplus);
    FF_MixArrZfromTVPCSAFT(rule,&Tplus,V,numSubs,data,pintParam,x,&ArrTplus,&ZTplus);
    FF_MixArrZfromTVPCSAFT(rule,&Tminus,V,numSubs,data,pintParam,x,&ArrTminus,&ZTminus);
    result[0]=Arr;//This is Arr
    result[1]=(1- Z)/ *V;//dArr/dV at constant T
    result[2]=((1- ZVplus)/ Vplus-result[1])/dV;//d2Arr/dV2 at constant T
    result[3]=(ArrTplus-Arr)/dT;//dArr/T at constant V
    result[4]=(result[3]-(Arr-ArrTminus)/dT)/dT;//d2Arr/dT2 at constant V
    result[5]=((1- ZTplus)/ *V-result[1])/dT;//d2Arr/dVdT
}


//Mixtures common calculations
//============================

//Mixture P calculation from T and V by eos
//-----------------------------------------
void CALLCONV FF_MixPfromTVeos(const enum FF_EOS eos[],const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,const void *data,const double pintParam[],
                           const double x[],double *P)
{
     FF_CubicParam param;
    double Arr,Z;
    double dTheta_dXi[*numSubs],db_dXi[*numSubs],dc_dXi[*numSubs];
    switch (*eos)
    {
        case FF_IdealGas:
            *P=R* *T/ *V;
            break;
        case FF_PCSAFT:
        case FF_DPCSAFT_GV:
        case FF_DPCSAFT_JC:
            //*( FF_SaftEOSdata*) data;
            FF_MixArrZfromTVPCSAFT(rule,T,V,numSubs,data,pintParam,x,&Arr,&Z);
            //printf("T:%f V:%f Z:%f\n",*T,*V,Z);
            *P=Z*R* *T/ *V;
            break;
        case FF_SW:
            //FF_PfromTVsw(T,V,data,P);
            break;
        default://if we have a cubic eos the first step is to calculate its parameters
            //*( FF_CubicEOSdata*) data;
            FF_MixParamCubicEOS(rule,T,numSubs,data,pintParam,x,&param,dTheta_dXi,db_dXi,dc_dXi);
            FF_PfromTVcubic(T,V,&param,P);
            break;
    }
}

//Mixture V calculation from T and P by eos
//-----------------------------------------
void CALLCONV FF_MixVfromTPeos(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const void *data,const double pintParam[],
                           const double x[],const char *option,double resultL[3],double resultG[3],char *state)
{
     FF_CubicParam param;//if we have a cubic eos the first step is to calculate its parameters
    double dTheta_dXi[*numSubs],db_dXi[*numSubs],dc_dXi[*numSubs];
    switch (*eosType)
    {
        case FF_IdealType:
            resultL[0]=resultG[0]=R* *T/ *P;
            resultL[1]=resultG[1]=0;
            resultL[2]=resultG[2]=1;
            *state='u';
            break;
        case FF_SAFTtype:
            //( FF_SaftEOSdata*) data;
            FF_MixVfromTPPCSAFT(rule,T,P,numSubs,data,pintParam,x,option,resultL,resultG,state);
            break;
        case FF_SWtype:
            //*( FF_SWEOSdata*) data;
            //FF_VfromTPsw(eos,T,P,data,option,resultL,resultG,state);
            break;
        default://Cubic eos
            FF_MixParamCubicEOS(rule,T,numSubs,data,pintParam,x,&param,dTheta_dXi,db_dXi,dc_dXi);
            //printf("a:%f, b:%f\n",param.a,param.b);
            FF_VfromTPcubic(T,P,&param,option,resultL,resultG,state);
            break;
    }

    //printf("T:%f  P:%f Vl:%f ArrL:%f Zl:%f\n",*T,*P,resultL[0],resultL[1],resultL[2]);
    if (*option=='s'){
        if (*state=='b'){
            if (fabs((resultL[0]-resultG[0])/resultL[0])>0.001){
                if ((resultL[1]+resultL[2]-1-log(resultL[2]))<(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='L';//we compare Gdr
                else if ((resultL[1]+resultL[2]-1-log(resultL[2]))>(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='G';
                else *state='E';//if Gdr is the same we are in equilibrium
            }
            else *state='U';
        }
    }
}

//Mixture Ideal gas thermodynamic properties calculation, from a reference state, specified by T and P, where H and S are 0
void CALLCONV FF_MixFF_IdealThermoEOS(const int *numSubs,const  FF_Correlation cp0[],const double x[],double *refT,double *refP, FF_ThermoProperties *th0)
{
    //printf("Hola, soy FF_MixFF_IdealThermoEOS\n");
     FF_ThermoProperties th0S;//Ideal thermodynamic properties of the selected substance
    int i;
    th0S.T=th0->T;
    th0S.V=th0->V;
    th0->Cp=th0->H=th0->S=0;
    //printf("I am going for individual calculation. Num subs: %i\n",*numSubs);
    //printf("%i %f %f %f %f %f\n",equation[0],coef[0][0],coef[0][1],coef[0][2],coef[0][3],coef[0][4]);
    for (i=0;i<*numSubs;i++){
        FF_IdealThermoEOS(&cp0[i].form,cp0[i].coef,refT,refP,&th0S);
        th0->Cp=th0->Cp+x[i]*th0S.Cp;
        th0->H=th0->H+x[i]*th0S.H;
        th0->S=th0->S+x[i]*th0S.S;
    }
    th0->Cv=th0->Cp-R;
    th0->U=th0->H-R*(th0->T- *refT);//We need to substrat the integration of d(P*V)=d(R*T)=R*dT from reference T to actual T
    th0->A=th0->U-th0->T*th0->S;
    th0->G=th0->H-th0->T*th0->S;
}


//Mixture Residual thermodynamic properties calculation from T and V, using EOS
void CALLCONV FF_MixResidualThermoEOS(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const int *numSubs,const void *data,const double pintParam[],
                                   const double x[], FF_ThermoProperties *thR)
{
    double ArrDer[6];
     FF_CubicParam param;
    double dN2Theta_dNi[*numSubs],dNb_dNi[*numSubs],dc_dXi[*numSubs];
    double delta,tau;
    switch (*eosType)
    {
        case FF_SAFTtype:
            FF_MixArrDerPCSAFT(rule,&thR->T,&thR->V,numSubs,data,pintParam,x,ArrDer);
            break;
        case FF_SWtype:
            delta=1/(thR->V*((( FF_SWEOSdata*)data)->rhoRef));
            tau=((( FF_SWEOSdata*)data)->tRef)/thR->T;
            //printf("delta:%f  tau:%f\n",delta,tau);
            //FF_ArrDerSW(&tau,&delta,data,ArrDer);
            break;
        default://Cubic EOS
            switch(*rule){
            case FF_VdW:
            case FF_PR:
            case FF_MKP:
                 FF_MixParamCubicEOS(rule,&thR->T,numSubs,data,pintParam,x,&param,dN2Theta_dNi,dNb_dNi,dc_dXi);
                break;
            }
        //else//gE mixing rules
            FF_ArrDerCubic(&thR->T,&thR->V,&param,ArrDer);
            break;
    }
    //printf("Arr:%f dArr/dV:%f d2Arr/dV2:%f dArr/dT:%f d2Arr/dT2:%f d2Arr/dVdT:%f\n",ArrDer[0],ArrDer[1],ArrDer[2],ArrDer[3],ArrDer[4],ArrDer[5]);
    if (*eosType==FF_SWtype)
    {
        thR->P=R*thR->T*(1+delta*ArrDer[1])/thR->V;
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=R*(tau*ArrDer[3]-ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        thR->H=thR->G+thR->T*thR->S;
        thR->dP_dV=-R*thR->T*(1+2*delta*ArrDer[1]+delta*delta*ArrDer[2])/thR->V/thR->V;
        thR->dP_dT=R*(1+delta*ArrDer[1]-delta*tau*ArrDer[5])/thR->V;
        thR->Cv=-R*tau*tau*ArrDer[4];
        thR->Cp=thR->Cv+R*pow((1+delta*ArrDer[1]-delta*tau*ArrDer[5]),2)/(1+2*delta*ArrDer[1]+delta*delta*ArrDer[2]);

    }
    else
    {
        thR->P=R*thR->T*(1/thR->V-ArrDer[1]);
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=-R*(thR->T*ArrDer[3]+ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        //thR->U=-R*thR->T*thR->T*ArrDer[3];
        //printf("%f %f\n",thR->U,-R*thR->T*thR->T*ArrDer[3]);
        thR->H=thR->G+thR->T*thR->S;
        //thR->H=thR->U-R*thR->T*thR->V*ArrDer[1];
        //thR->H=-R*thR->T*(thR->V*ArrDer[1]+thR->T*ArrDer[3]);
        thR->Cv=-2*R*thR->T*ArrDer[3]-R*thR->T*thR->T*ArrDer[4];
        thR->dP_dV=R*thR->T*(-1/thR->V/thR->V-ArrDer[2]);
        thR->dP_dT=R*(1/thR->V-ArrDer[1]-thR->T*ArrDer[5]);
        //thR->Cp=thR->Cv-thR->T*pow(R*(1/thR->V-ArrDer[1]-thR->T*ArrDer[5]),2)/(R*thR->T*(-1/thR->V/thR->V-ArrDer[2]))-R;
        thR->Cp=thR->Cv-thR->T*thR->dP_dT*thR->dP_dT/thR->dP_dV-R;

    }
}

//Mixture thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0
void CALLCONV FF_MixThermoEOS(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const int *numSubs,const void *data,const double pintParam[],
                                 const double x[],const  FF_Correlation cp0[],double *refT,double *refP, FF_ThermoProperties *th)
{
    //printf("Hi, I am FF_MixThermoEOS\n");
    th->MW=0;
    int i;
     FF_ThermoProperties th0,thR;
    switch (*eosType)
    {
        case FF_SAFTtype:
            for (i=0;i<*numSubs;i++) th->MW=th->MW+x[i]*(( FF_SaftEOSdata *)data+i)->MW;
            break;
        case FF_SWtype:
            for (i=0;i<*numSubs;i++) th->MW=th->MW+x[i]*(( FF_SWEOSdata *)data+i)->MW;
            break;
        default://Cubic EOS
            for (i=0;i<*numSubs;i++) th->MW=th->MW+x[i]*(( FF_CubicEOSdata *)data+i)->MW;
            break;
    }
    th0.MW=thR.MW=th->MW;//Perhaps not necessary
    th0.T=thR.T=th->T;
    th0.V=thR.V=th->V;
    //printf("Now going to calculate Ideal part\n");
    //printf("%i %i %f\n",*numSubs,equation[0],coef[0][0]);
    FF_MixFF_IdealThermoEOS(numSubs,cp0,x,refT,refP,&th0);
    if (*eosType==FF_IdealType)
    {
        th->P=R*th->T/th->V;
        th->A=th0.A;
        th->G=th0.G;
        th->S=th0.S;
        th->U=th0.U;
        th->H=th0.H;
        th->dP_dV=-R*th->T/th->V/th->V;
        th->dP_dT=R/th->V;
        th->Cv=th0.Cv;
        th->Cp=th0.Cp;
    }
    else
    {
       // printf("Now going for residual part\n");
        FF_MixResidualThermoEOS(eosType,rule,numSubs,data,pintParam,x,&thR);
        th->P=thR.P;
        th->A=th0.A+thR.A;
        th->G=th0.G+thR.G;
        th->S=th0.S+thR.S;
        th->U=th0.U+thR.U;
        th->H=th0.H+thR.H;
        th->dP_dV=thR.dP_dV;
        th->dP_dT=thR.dP_dT;
        th->Cv=th0.Cv+thR.Cv;
        th->Cp=th0.Cp+thR.Cp;
        /*
        th->P=R*th->T/th->V;
        th->A=th0.A;
        th->G=th0.G;
        th->S=th0.S;
        th->U=th0.U;
        th->H=th0.H;
        th->dP_dV=-R*th->T/th->V/th->V;
        th->dP_dT=R/th->V;
        th->Cv=th0.Cv;
        th->Cp=th0.Cp;*/
    }
    th->SS=th->V*pow(-th->Cp*th->dP_dV/th->Cv/th->MW*1e3,0.5);
    th->JT=-(th->T*th->dP_dT/th->dP_dV+th->V)/th->Cp;//Joule Thomson coefficient =(dT/dP) at constant H
    th->IT=-th->JT*th->Cp;//Isothermal throttling coefficient = (dH/dP) at constant T
    //printf("MW:%f T:%f V:%f P:%f Cv:%f Cp:%f H:%f S:%f dP_dV:%f dP_dT:%f SS:%f JT:%f IT:%f\n",th->MW,th->T,th->V,th->P,th->Cv,th->Cp,th->H,th->S,th->dP_dV,th->dP_dT,th->SS,th->JT,th->IT);
    //printf("Ideal T:%f V:%f P:%f Cv0:%f Cp0:%f H0:%f S0:%f\n",th0.T,th0.V,th0.P,th0.Cv,th0.Cp,th0.H,th0.S);
    //printf("Residual T:%f V:%f P:%f Cvr:%f Cpr:%f Ur:%f Hr:%f Sr:%f\n",thR.T,thR.V,thR.P,thR.Cv,thR.Cp,thR.U,thR.H,thR.S);
}


EXP_IMP void CALLCONV FF_MixFugacityEOS(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,const void *data,
                                        const double pintParam[],const double x[],double phi[])
{
     FF_CubicParam param;
    double xPlus[*numSubs],ArrPlus,Zplus;
    double Arr,Z;
    double dTheta_dXi[*numSubs],db_dXi[*numSubs],dc_dXi[*numSubs],dArr_dXi[*numSubs],sumdArr_dXi=0;
    int i;
    switch (*eosType)//Mixture parameters calculation, and initial solution of the EOS for the given T and P
    {
        case FF_SAFTtype:
            for (i=0;i<*numSubs;i++) xPlus[i]=x[i];
            FF_MixArrZfromTVPCSAFT(rule,T,V,numSubs,data,pintParam,x,&Arr,&Z);
            for (i=0;i<*numSubs;i++)
            {
                xPlus[i]=x[i]*1.00001; //This is the composition increment for the calculation of the variation of red.res.Helmholtz energy
                //printf("i:%i xPlus:%f\n",i,x[i]);
                FF_MixArrZfromTVPCSAFT(rule,T,V,numSubs,data,pintParam,xPlus,&ArrPlus,&Zplus);
                dArr_dXi[i]=(ArrPlus-Arr)/(xPlus[i]-x[i]);
                sumdArr_dXi=sumdArr_dXi+x[i]*dArr_dXi[i];
                xPlus[i]=x[i];
            }
            break;
        default:
            FF_MixParamCubicEOS(rule,T,numSubs,data,pintParam,x,&param,dTheta_dXi,db_dXi,dc_dXi);
            FF_ArrZfromTVcubic(T,V,&param,&Arr,&Z);
            for (i=0;i<*numSubs;i++)
            {
                //dArr_dXi[i]=param.Theta/(param.b*R* *T*(param.w-param.u))*((dTheta_dXi[i]/param.Theta-db_dXi[i]/param.b)*log((*V+param.u*param.b)/(*V+param.w*param.b))+
                //            *V*db_dXi[i]*(param.u-param.w)/(*V+param.u*param.b)/(*V+param.w*param.b))+db_dXi[i]/(*V-param.b);//This is without volume traslation
                dArr_dXi[i]=param.Theta/(param.b*R* *T*(param.w-param.u))*((dTheta_dXi[i]/param.Theta-db_dXi[i]/param.b)*log((*V+param.c+param.u*param.b)/
                           (*V+param.c+param.w*param.b))+((*V+param.c)*db_dXi[i]-param.b*dc_dXi[i])*(param.u-param.w)/
                           ((*V+param.c+param.u*param.b)*(*V+param.c+param.w*param.b)))+(db_dXi[i]-dc_dXi[i])/(*V+param.c-param.b);
                sumdArr_dXi=sumdArr_dXi+x[i]*dArr_dXi[i];
            }
            break;
    }
    for (i=0;i<*numSubs;i++)
    {
        ;
        phi[i]=exp(dArr_dXi[i]+Arr-sumdArr_dXi+Z-1-log(Z));
        //printf("phi:%f\n",phi[i]);
    }
}


/*
//Fugacity coefficients calculation for a mixture, given T and P, according to a EOS
void CALLCONV calcMixPhiEOS(const enum FF_EOS *eos,const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const void *sData,
        const double correction[],const double x[],const char *option,double fugacityCoefL[],double fugacityCoefG[],char *state)
{
    int i;
    double answer3[3],answerG[3],answer[2],Zl,Zg,ArrL,ArrG,Vl,Vg;
    double answerP[2];
     FF_CubicParam sParam[*numSubs],param,paramPlus;
    double gE=0;//excess mixture Gibbs energy
    double dN2Theta_dNi[*numSubs],dNb_dNi[*numSubs];//Derivatives of Theta and b regarding each substance number of moles
    double alpha=0,alphai[*numSubs],dNalpha_dNi[*numSubs];//alpha=Theta/(bRT), and its derivatives regarding number of moles of each substance
    double A,B;
    double xPlus[*numSubs],substdAl_dx[*numSubs],substdAg_dy[*numSubs],sumdAl_dx=0.0,sumdAg_dy=0.0;
    for (i=0;i<*numSubs;i++) xPlus[i]=x[i];
    switch (*eos)//Mixture parameters calculation, and initial solution of the EOS for the given T and P
    {
    case FF_PR78:
    case FF_PRSV1:
    case FF_SRK:
        for (i=0;i<*numSubs;i++)
        {
            FF_FixedParamCubic(&sData[i],&sParam[i]);
            FF_ThetaDerivCubic(T,&sData[i],&sParam[i]);
            //calcParamCubicEOS(eos,T,sData[i],sParam[i]);//We obtain the individual substances parameters
        }
        switch (*rule)//And now the mixture parameters using the selected mixing rule
        {
            case FF_VdW:
                FF_MixParamCubicEOS(eos,rule,T,numSubs,sData,correction,x,&param,dN2Theta_dNi,dNb_dNi);
                printf("%f\n",param.Theta);
                //calcMixParamCubicEOS(eos,rule,T,numSubs,sParam,correction,x,param,dN2Theta_dNi,dNb_dNi);
                break;
            case MHV1://gE EOS
            case PSRK:
            case HVOS:
            case LCVM:
            case MHV2:
                for (i=0;i<numSubs;i++) gE=gE+x[i]*log(correction[i]);
                //gE=gE*R*T;
                //calcMixParamCubicEOSgE(eos,rule,T,numSubs,sParam,gE,x,param,dNb_dNi);
                break;
        }
        FF_VfromTPeos(eos,T,P,sData,option,answer3,answerG,state);
        //FF_VfromTPcubic(T,P,param,answer3,optLiquid);//Here we obtain the volumes that match the requested pressure, at given temperature
        //FF_VfromTPcubic(T,P,param,answerG,optGas);
        alpha=param.Theta/(param.b*R* *T);
        for (i=0;i<numSubs;i++)
        {
            alphai[i]=sParam[i].Theta/(sParam[i].b*R* *T);

        }

        //A = param.Theta* *P/ pow((R* *T),2);
        //B = param.b * *P / (R * *T);
        break;
    case FF_PCSAFT:
    case FF_DPCSAFT_JC:
        //calcMixVolFF_PCSAFT(eos,T,P,numSubs,sData,correction,x,answer3,optLiquid);//First we obtain the volumes that match the requested pressure
        //calcMixVolFF_PCSAFT(eos,T,P,numSubs,sData,correction,x,answerG,not(optGas));
        A=A;
        break;
    }

    //Charge of the fixed parameters for each gE EOS
    double q1=0,q2=0,lambda=0;//Coefficients for  alpha calculation. lambda is for LCVM mixing rule
    switch (*rule)
    {
        case MHV1:
            if ((*eos==FF_PR78) || (*eos==FF_PRSV1)) q1=-0.53;
            if (*eos==FF_SRK) q1=-0.593;
            break;
        case PSRK:
            if (*eos==FF_SRK) q1=-0.64663;
            break;
        case HVOS:
            if ((*eos==FF_PR78) || (*eos==FF_PRSV1)) q1=-0.62323;
            break;
        case LCVM:
            if ((*eos==FF_PR78) || (*eos==FF_PRSV1))
            {
                q1=-0.52;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.623;///Av. Huron-Vidal(infinite pressure) part coefficient
            }
            if (eos==FF_SRK)
            {
                q1=-0.593;//Am. Michelsen(zero pressure) part coefficient
                q2=-0.693;///Av. Huron-Vidal(infinite pressure) part coefficient
            }
            lambda=0.36;
            break;
        case MHV2:
            if ((eos==FF_PR78) || (eos==FF_PRSV1))
            {
                q1=-0.4347;
                q2=-0.003654;
            }
            if (eos==FF_SRK)
            {
                q1=-0.4783;
                q2=-0.0047;
            }
            break;
    }

    //Calculation of N*alpha derivatives regarding each component number of moles (dNalpha_dNi)
    switch (*rule)
    {
    case FF_VdW:
        for (i=0;i< *numSubs;i++) dNalpha_dNi[i]=alpha*(dN2Theta_dNi[i]/param.Theta-dNb_dNi[i]/param.b);
        break;
    case MHV1:
    case PSRK:
    case HVOS:
        for (i=0;i< *numSubs;i++) dNalpha_dNi[i]=alphai[i]+(log(correction[i])+log(param.b/sParam[i].b)+sParam[i].b/param.b-1)/q1;
        break;
    case LCVM:
        for (i=0;i< *numSubs;i++) dNalpha_dNi[i]=alphai[i]+(lambda/q2+(1-lambda)/q1)*log(correction[i])-(1-lambda)/0.52*
        (log(param.b/sParam[i].b)+sParam[i].b/param.b-1);
        break;
    case MHV2:
        for (i=0;i< *numSubs;i++)  dNalpha_dNi[i]=(q1*alphai[i]+q2*(pow(alpha,2)+pow(alphai[i],2))+log(correction[i])+log(param.b/sParam[i].b)+sParam[i].b/param.b-1)/
                (q1+2*q2*alpha);
        break;
    default:
        for (i=0;i< *numSubs;i++) dNalpha_dNi[i]=0;
        break;
    }

    //And now the calculation of activity coefficient for each substance
    if ((*option!='g') && (answer3[0]>0)) //optional calculation of the liquid phase
    {
        Zl=answer3[2];//Zl. dArr/dV=(1-Z)/V
        ArrL=answer3[1];//Reduced residual Helmholtz energy of the liquid ArrL
        Vl=answer3[0];//Vl
        switch (*eos)
        {
        case FF_PCSAFT:
        case FF_DPCSAFT_JC:
            for (i=0;i<numSubs;i++)
            {
                xPlus[i]=x[i]*1.00001; //This is the composition increment for the calculation of the variation of red.res.Helmholtz energy
                //calcMixPresFF_PCSAFT(eos,T,Vl,numSubs,sData,correction,xPlus,answerP);
                substdAl_dx[i]=(answerP[1]-ArrL)/(xPlus[i]-x[i]);
                sumdAl_dx=sumdAl_dx+x[i]*substdAl_dx[i];
                xPlus[i]=x[i];
            }
            for (i=0;i<numSubs;i++) fugacityCoefL[i]=exp(substdAl_dx[i]-sumdAl_dx+ArrL+Zl-1-log(Zl));//Liquid fugacity coefficients
            break;
        case FF_PR78:
        case FF_PRSV1:
        case FF_SRK:
            for (i=0;i< *numSubs;i++)	fugacityCoefL[i]=exp(dNb_dNi[i]/param.b*(Zl-1)-log((Vl-param.b)*Zl/Vl)+1/(param.w-param.u)*
                    dNalpha_dNi[i]*log((Vl+param.u*param.b)/(Vl+param.w*param.b)));
                //Is equal for both FF_PR78 and FF_SRK
            break;

        default:
            for (i=0;i< *numSubs;i++) fugacityCoefG[i]=0;
            break;
        }

    }

    else for (i=0;i<numSubs;i++) fugacityCoefL[i]=0;

    if ((*option!='l') && (answerG[0]>0)) //Optional calculation for the gas phase
    {
        Zg=answerG[2];//Zg
        ArrG=answerG[1];//Reduced residual Helmholtz energy of the gas ArrG
        Vg=answerG[0];//Vg
        switch (*eos)
        {
        case FF_PCSAFT:
        case FF_DPCSAFT_JC:
            for (i=0;i<numSubs;i++)
            {
                xPlus[i]=x[i]*1.0001; //This is the composition increment for the calculation of the variation of red.res.Helmholtz energy
                    //calcMixPresFF_PCSAFT(eos,T,Vg,numSubs,sData,correction,xPlus,answerP);
                    substdAg_dy[i]=(answerP[1]-ArrG)/(xPlus[i]-x[i]);
                    sumdAg_dy=sumdAg_dy+x[i]*substdAg_dy[i];
                xPlus[i]=x[i];
            }
            for (i=0;i<numSubs;i++) fugacityCoefG[i]=exp(substdAg_dy[i]-sumdAg_dy+ArrG+Zg-1-log(Zg));//Gas fugacity coefficients
            break;
        case FF_PR78:
        case FF_PRSV1:
        case FF_SRK:
            for (i=0;i< *numSubs;i++)	fugacityCoefG[i]=exp(dNb_dNi[i]/param.b*(Zg-1)-log((Vg-param.b)*Zg/Vg)+1/(param.w-param.u)*
                                dNalpha_dNi[i]*log((Vg+param.u*param.b)/(Vg+param.w*param.b)));
            break;
        default:
            for (i=0;i<numSubs;i++) fugacityCoefG[i]=0;
            break;
        }
    }
    else for (i=0;i<numSubs;i++) fugacityCoefG[i]=0;

}
    */

//Mixture bubble temperature calculation, given P, composition, eos and mixing rule
void CALLCONV FF_BubbleT(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *P,const int *numSubs,const void *data,
                     const double pintParam[],const double x[],const double *bTguess, double *bT,double y[],double substPhiL[],double substPhiG[])
{
    //printf("%i %i %f %i %f %f %f %f\n",*eosType,*rule,*P,*numSubs,pintParam[6],pintParam[12],x[0],x[1]);
    *bT=0;
    double Tmax=800;
    double T;
    //double substPhiL[*numSubs],substPhiG[*numSubs];
    double newY[*numSubs],yTotal,dyTotal;
    double answerL[3],answerG[3];
    int i,n=4,counter=0;
    char option,state;
    if (*bTguess==0) T=Tmax*0.5;
    else T=*bTguess;
    option='b';
    FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
    //printf("BubbleT initial guess, T:%f state:%c Vl:%f Zl:%f Vg:%f Zg:%f\n",T,state,answerL[0],answerL[2],answerG[0],answerG[2]);

    if (*bTguess==0){//If no guess given, we approximate the temperature
        //Till we find a temperature with different possitive liquid and gas solutions, we approximate a liquid phase close to boil
        while (((state=='l')||(state=='g')||(fabs(answerL[2]-answerG[2])<0.01))&&(counter<7)){
            if ((state=='l')||(answerL[2]<=0.2)) T=T+Tmax/n;
            else if((state=='g')||(answerL[2]>=0.25)) T=T-Tmax/n;
            else if (state=='f'){
                *bT=0;
                return;
            }
            FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
            n=n*2;
            counter++;
            //printf("Approximating T, counter:%i n:%i T:%f state:%c  Vl:%f  Zl:%f  Vg%f  Zg%f\n",counter,n/2,T,state,answerL[0],answerL[2],answerG[0],answerG[2]);
        }
    }

    yTotal=0;
    FF_MixFugacityEOS(eosType,rule,&T,&answerL[0],numSubs,data,pintParam,x,substPhiL);//phi as liquid phase
    if ((state=='b')&&(fabs(answerL[2]-answerG[2])>0.01)){
        FF_MixFugacityEOS(eosType,rule,&T,&answerG[0],numSubs,data,pintParam,x,substPhiG);//phi as gas phase
        for (i=0;i<*numSubs;i++){
            y[i]=x[i]*substPhiL[i]/substPhiG[i];//Initial composition of gas phase
            yTotal=yTotal+y[i];
        }
    }
    else{
        for (i=0;i<*numSubs;i++){
            y[i]=x[i]*substPhiL[i];//we assume that phi of gas is 1
            yTotal=yTotal+y[i];
        }
    }


    //printf("yTotal inicial: %f\n",yTotal);
    //T=T/pow(yTotal,0.06);
    for (i=0;i<*numSubs;i++) y[i]=y[i]/yTotal;//We obtain the normalized gas phase composition
    //printf("initial liquid phase: %f  %f  gas phase: %f  %f  phis: %f  %f  %f  %f\n",x[0],x[1],y[0],y[1],substPhiL[0],substPhiL[1],substPhiG [0],substPhiG[1]);
    do
    {
        T=T/pow(yTotal,0.06);
        option='l';
        FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
        FF_MixFugacityEOS(eosType,rule,&T,&answerL[0],numSubs,data,pintParam,x,substPhiL);
        option='g';
        FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
        FF_MixFugacityEOS(eosType,rule,&T,&answerG[0],numSubs,data,pintParam,y,substPhiG);
        yTotal=0;
        for (i=0;i<*numSubs;i++)
        {
            y[i]=x[i]*substPhiL[i]/substPhiG[i];
            yTotal=yTotal+y[i];
        }
        for (i=0;i<*numSubs;i++) y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
        dyTotal=1.0;
        //printf("T:%f substPhiL:%f,%f substPhiG:%f,%f\n",T,substPhiL[0],substPhiL[1],substPhiG[0],substPhiG[1]);
        while (dyTotal > 0.001)
        {
            FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
            FF_MixFugacityEOS(eosType,rule,&T,&answerG[0],numSubs,data,pintParam,y,substPhiG);
            yTotal=0;
            for (i=0;i<*numSubs;i++)
            {
                newY[i]=x[i]*substPhiL[i]/substPhiG[i];
                yTotal=yTotal+newY[i];
            }
            for (i=0;i<*numSubs;i++) newY[i]=newY[i]/yTotal; //we normalize y in order to obtain molar fractions
            dyTotal=0;
            for (i=0;i<*numSubs;i++)
            {
                dyTotal=dyTotal+fabs(newY[i]-y[i]);
                y[i]=newY[i];
            }
            counter=counter+1;
            if (counter>2000) return;
        }
        //T=T/pow(yTotal,0.06);
        //printf("counter:%i yTotal:%f\n",counter,yTotal);
    } while (fabs(yTotal - 1) > 0.001);

    *bT=T;
    //printf("\nBubbleT:%f x[0]:%f x[1]:%f y[0]:%f y[1]:%f n:%i\n\n",T,x[0],x[1],y[0],y[1],counter);
}


//Mixture dew temperature calculation, given P, comoposition, eos and mixing rule
void CALLCONV FF_DewT(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *P,const int *numSubs,const void *data,
                     const double pintParam[],const double y[],const double *dTguess,double *dT,double x[],double substPhiL[],double substPhiG[])
{
    //printf("%i %i %f %i %f %f %f %f %f %f\n",*eosType,*rule,*P,*numSubs,pintParam[0],pintParam[6],pintParam[12],pintParam[18],y[0],y[1]);
    // FF_CubicEOSdata *cd=( FF_CubicEOSdata*)data;
    //printf("%f %f %f %f %f %f %f %f\n",cd[0].Tc,cd[0].Pc,cd[0].w,cd[0].c,cd[1].Tc,cd[1].Pc,cd[1].w,cd[1].c);
    *dT=0;
    double Tmax=800;
    double T;
    double newX[*numSubs],xTotal,dxTotal;
    double answerL[3],answerG[3];
    int i,n=4,counter=0;
    char option,state;
    if (*dTguess==0) T=Tmax*0.5;
    else T=*dTguess;
    option='b';
    FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
    //printf("DewT initial guess, P:%f T:%f state:%c Vl:%f Zl:%f Vg:%f Zg:%f\n",*P,T,state,answerL[0],answerL[2],answerG[0],answerG[2]);

    if (*dTguess==0){//If no guess given, we approximate the temperature
        //Till we find a temperature with different possitive liquid and gas solutions, we approximate a liquid phase close to boil
        while (((state=='l')||(state=='g')||(fabs(answerL[2]-answerG[2])<0.01))&&(counter<7)){
               if ((state=='l')||(answerL[2]<=0.2)) T=T+Tmax/n;
               else if((state=='g')||(answerL[2]>=0.25)) T=T-Tmax/n;
               else if (state=='f'){
                   *dT=0;
                   return;
               }
               FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
               n=n*2;
               counter++;
           //printf("Approximating T, counter:%i n:%i T:%f state:%c  Vl:%f  Zl:%f  Vg:%f  Zg:%f\n",counter,n/2,T,state,answerL[0],answerL[2],answerG[0],answerG[2]);
           }
    }

    xTotal=0;
    FF_MixFugacityEOS(eosType,rule,&T,&answerL[0],numSubs,data,pintParam,y,substPhiL);//phi as liquid phase
    if ((state=='b')&&(fabs(answerL[2]-answerG[2])>0.01)){
        FF_MixFugacityEOS(eosType,rule,&T,&answerG[0],numSubs,data,pintParam,y,substPhiG);//phi as gas phase
        for (i=0;i<*numSubs;i++){
            x[i]=y[i]*substPhiG[i]/substPhiL[i];//Initial composition of liquid phase
            xTotal=xTotal+x[i];
        }
    }
    else{
        for (i=0;i<*numSubs;i++){
            x[i]=y[i]/substPhiL[i];//we assume that phi of gas is 1
            xTotal=xTotal+x[i];
        }
    }
    //printf("xTotal inicial: %f\n",xTotal);
    for (i=0;i<*numSubs;i++) x[i]=x[i]/xTotal;//normalized initial liquid phase composition
    //printf("initial liquid phase: %f  %f  gas phase: %f  %f phis: %f %f %f %f\n",x[0],x[1],y[0],y[1],substPhiL[0],substPhiL[1],substPhiG [0],substPhiG[1]);
    //T=T*pow(xTotal,0.06);
    //Now we begin a double loop. The external one modifies T in order to obtain sum of liquid fractions equal to 1
    //The internal one equals fugacities of liquid phase to those of the gas,modifying liquid concentrations
    do
    {
        T=T*pow(xTotal,0.06);
        option='g';
        FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);//If it seems liquid phase, we shoould increase T. To do
        FF_MixFugacityEOS(eosType,rule,&T,&answerG[0],numSubs,data,pintParam,y,substPhiG);
        //printf("Finding T counter:%i n:%i T:%f state:%c  Vl:%f  Zl:%f  Vg:%f  Zg:%f\n",counter,n/2,T,state,answerL[0],answerL[2],answerG[0],answerG[2]);
        option='l';
        FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);//If it seems gas phase, we shoould increase T. To do
        FF_MixFugacityEOS(eosType,rule,&T,&answerL[0],numSubs,data,pintParam,x,substPhiL);
        //printf("Finding T counter:%i n:%i T:%f state:%c  Vl:%f  Zl:%f  Vg:%f  Zg:%f\n",counter,n/2,T,state,answerL[0],answerL[2],answerG[0],answerG[2]);
        xTotal=0;
        for (i=0;i<*numSubs;i++)
        {
            x[i]=y[i]*substPhiG[i]/substPhiL[i];
            xTotal=xTotal+x[i];
        }
        for (i=0;i<*numSubs;i++) x[i]=x[i]/xTotal; //we normalize y in order to obtain molar fractions
        dxTotal=1.0;
        //printf("T:%f substPhiL:%f,%f substPhiG:%f,%f\n",T,substPhiL[0],substPhiL[1],substPhiG[0],substPhiG[1]);
        while (dxTotal > 0.001)
        {
            FF_MixVfromTPeos(eosType,rule,&T,P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
            FF_MixFugacityEOS(eosType,rule,&T,&answerL[0],numSubs,data,pintParam,x,substPhiL);
            xTotal=0;
            for (i=0;i<*numSubs;i++)
            {
                newX[i]=y[i]*substPhiG[i]/substPhiL[i];
                xTotal=xTotal+newX[i];
            }
            for (i=0;i<*numSubs;i++) newX[i]=newX[i]/xTotal; //we normalize y in order to obtain molar fractions
            dxTotal=0;
            for (i=0;i<*numSubs;i++)
            {
                dxTotal=dxTotal+fabs(newX[i]-x[i]);
                x[i]=newX[i];
            }
            counter=counter+1;
            //printf("counter:%i dxTotal:%f\n",counter,dxTotal);
            if (counter>2000) return;
        }
        //T = T * (1 - (1 - xTotal) * 0.06);
        //T=T*pow(xTotal,0.06);
        //printf("counter:%i T:%f xTotal:%f\n",counter,T,xTotal);
    }while (fabs(xTotal - 1) > 0.001);
    *dT=T;
    //printf("DewT:%f y[0]:%f y[1]:%f x[0]:%f x[1]:%f n:%i\n",T,y[0],y[1],x[0],x[1],counter);
}


//Mixture bubble pressure calculation, given T, comoposition, eos and mixing rule
void CALLCONV FF_BubbleP(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const int *numSubs,const void *data,
                     const double pintParam[],const double x[],const double *bPguess, double *bP,double y[],double substPhiL[],double substPhiG[]){
    *bP=0;
    double Pmax=200e5;
    double P;
    //double substPhiL[*numSubs],substPhiG[*numSubs];
    double newY[*numSubs],yTotal,dyTotal;
    double answerL[3],answerG[3];
    int i,n=4,counter=0;
    char option,state;
    if (*bPguess==0) P=Pmax*0.5;
    else P=*bPguess;
    option='b';
    FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
    //printf("BubbleP initial guess, T:%f P:%f state:%c Vl:%f Zl:%f Vg:%f Zg:%f\n",*T,P,state,answerL[0],answerL[2],answerG[0],answerG[2]);

    if (*bPguess==0){//If no guess given, we approximate the temperature
        //Till we find a temperature with different possitive liquid and gas solutions, we approximate a liquid phase close to boil
        while (((state=='l')||(state=='g')||(fabs(answerL[2]-answerG[2])<0.01))&&(counter<7)){
            if ((state=='l')||(answerL[2]<=0.2)) P=P-Pmax/n;
            else if((state=='g')||(answerL[2]>=0.25)) P=P+Pmax/n;
            else if (state=='f'){
                *bP=0;
                return;
            }
            FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
            n=n*2;
            counter++;
            //printf("Approximating P, counter:%i n:%i T:%f state:%c  Vl:%f  Zl:%f  Vg%f  Zg%f\n",counter,n/2,P,state,answerL[0],answerL[2],answerG[0],answerG[2]);
        }
    }

    yTotal=0;
    FF_MixFugacityEOS(eosType,rule,T,&answerL[0],numSubs,data,pintParam,x,substPhiL);//phi as liquid phase
    if ((state=='b')&&(fabs(answerL[2]-answerG[2])>0.01)){
        FF_MixFugacityEOS(eosType,rule,T,&answerG[0],numSubs,data,pintParam,x,substPhiG);//phi as gas phase
        for (i=0;i<*numSubs;i++){
            y[i]=x[i]*substPhiL[i]/substPhiG[i];//Initial composition of gas phase
            yTotal=yTotal+y[i];
        }
    }
    else{
        for (i=0;i<*numSubs;i++){
            y[i]=x[i]*substPhiL[i];//we assume that phi of gas is 1
            yTotal=yTotal+y[i];
        }
    }


    //printf("yTotal inicial: %f\n",yTotal);
    //P=P*yTotal;
    for (i=0;i<*numSubs;i++) y[i]=y[i]/yTotal;//We obtain the normalized gas phase composition
    //printf("initial liquid phase: %f  %f  gas phase: %f  %f  phis: %f  %f  %f  %f\n",x[0],x[1],y[0],y[1],substPhiL[0],substPhiL[1],substPhiG [0],substPhiG[1]);

    do
    {
        P = P * yTotal;//(1 - (1 - yTotal) * 0.05);
        option='l';
        FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
        FF_MixFugacityEOS(eosType,rule,T,&answerL[0],numSubs,data,pintParam,x,substPhiL);
        option='g';
        FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
        //printf("y[i]:%f,%f Vg:%f\n",y[0],y[1],answerG[0]);
        FF_MixFugacityEOS(eosType,rule,T,&answerG[0],numSubs,data,pintParam,y,substPhiG);
        //printf("P:%f substPhiL:%f,%f substPhiG:%f,%f\n",P,substPhiL[0],substPhiL[1],substPhiG[0],substPhiG[1]);
        yTotal=0;
        for (i=0;i<*numSubs;i++)
        {
            y[i]=x[i]*substPhiL[i]/substPhiG[i];
            yTotal=yTotal+y[i];
        }
        for (i=0;i<*numSubs;i++) y[i]=y[i]/yTotal; //we normalize y in order to obtain molar fractions
        dyTotal=1.0;

       // printf("Counter:%i P:%f yTotal:%f\n",counter,P,yTotal);
        while (dyTotal > 0.001)
        {
            FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
            FF_MixFugacityEOS(eosType,rule,T,&answerG[0],numSubs,data,pintParam,y,substPhiG);
            yTotal=0;
            for (i=0;i<*numSubs;i++)
            {
                newY[i]=x[i]*substPhiL[i]/substPhiG[i];
                yTotal=yTotal+newY[i];
            }
            for (i=0;i<*numSubs;i++) newY[i]=newY[i]/yTotal; //we normalize y in order to obtain molar fractions
            dyTotal=0;
            for (i=0;i<*numSubs;i++)
            {
                dyTotal=dyTotal+fabs(newY[i]-y[i]);
                y[i]=newY[i];
            }
            counter=counter+1;
            if (counter>2000) return;
        }
        //P = P * yTotal;//(1 - (1 - yTotal) * 0.05);
        //printf("counter:%i dyTotal:%f\n",counter,dyTotal);
    } while (fabs(yTotal - 1) > 0.001);

    *bP=P;
    //printf("BubbleP:%f n:%i\n",P,counter);
}

//void CALLCONV FF_DewT(const enum FF_EOS eosType[],const enum FF_MixingRule *rule,const double *P,const int *numSubs,const void *data,
//                     const double pintParam[],const double y[],double *dT,double x[])

//Mixture dew pressure calculation, given T, comoposition, eos and mixing rule
void CALLCONV FF_DewP(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const int *numSubs,const void *data,
                     const double pintParam[],const double y[],const double *dPguess,double *dP,double x[],double substPhiL[],double substPhiG[])
{
    *dP=0;
    double Pmax=300e5;
    double P;
    //double substPhiL[*numSubs],substPhiG[*numSubs];
    double newX[*numSubs],xTotal=0,dxTotal;
    double answerL[3],answerG[3];
    int i,n=4,counter=0;
    char option,state;
    if (*dPguess==0) P=Pmax*0.5;
    else P=*dPguess;
    option='b';
    FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
    //printf("DewP initial guess, T:%f P:%f state:%c Vl:%f Zl:%f Vg:%f Zg:%f\n",*T,P,state,answerL[0],answerL[2],answerG[0],answerG[2]);

    if (*dPguess==0){//If no guess given, we approximate the temperature
        //Till we find a temperature with different possitive liquid and gas solutions, we approximate a liquid phase close to boil
        while (((state=='l')||(state=='g')||(fabs(answerL[2]-answerG[2])<0.01))&&(counter<7)){
               if ((state=='l')||(answerL[2]<=0.2)) P=P-Pmax/n;
               else if((state=='g')||(answerL[2]>=0.25)) P=P+Pmax/n;
               else if (state=='f'){
                   *dP=0;
                   return;
               }
               FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
               n=n*2;
               counter++;
           //printf("Approximating P, counter:%i n:%i T:%f state:%c  Vl:%f  Zl:%f  Vg:%f  Zg:%f\n",counter,n/2,P,state,answerL[0],answerL[2],answerG[0],answerG[2]);
           }
    }

    xTotal=0;
    FF_MixFugacityEOS(eosType,rule,T,&answerL[0],numSubs,data,pintParam,y,substPhiL);//phi as liquid phase
    if ((state=='b')&&(fabs(answerL[2]-answerG[2])>0.01)){
        FF_MixFugacityEOS(eosType,rule,T,&answerG[0],numSubs,data,pintParam,y,substPhiG);//phi as gas phase
        for (i=0;i<*numSubs;i++){
            x[i]=y[i]*substPhiG[i]/substPhiL[i];//Initial composition of liquid phase
            xTotal=xTotal+x[i];
        }
    }
    else{
        for (i=0;i<*numSubs;i++){
            x[i]=y[i]/substPhiL[i];//we assume that phi of gas is 1
            xTotal=xTotal+x[i];
        }
    }
    //printf("xTotal inicial: %f\n",xTotal);
    for (i=0;i<*numSubs;i++) x[i]=x[i]/xTotal;//normalized initial liquid phase composition
    //printf("initial liquid phase: %f  %f  gas phase: %f  %f phis: %f %f %f %f\n",x[0],x[1],y[0],y[1],substPhiL[0],substPhiL[1],substPhiG [0],substPhiG[1]);
    //P=P/xTotal;
    //Now we begin a double loop. The external one modifies P in order to obtain sum of liquid fractions equal to 1
    //The internal one equals fugacities of liquid phase to those of the gas,modifying liquid concentrations
    do
    {
        P = P / xTotal;//* (1 + (1 - xTotal) * 0.06);
        option='g';
        FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,y,&option,answerL,answerG,&state);
        FF_MixFugacityEOS(eosType,rule,T,&answerG[0],numSubs,data,pintParam,y,substPhiG);
        option='l';
        FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
        FF_MixFugacityEOS(eosType,rule,T,&answerL[0],numSubs,data,pintParam,x,substPhiL);
        xTotal=0;
        for (i=0;i<*numSubs;i++)
        {
            x[i]=y[i]*substPhiG[i]/substPhiL[i];
            xTotal=xTotal+x[i];
        }
        for (i=0;i<*numSubs;i++) x[i]=x[i]/xTotal; //we normalize y in order to obtain molar fractions
        dxTotal=1.0;
        //printf("P:%f substPhiL:%f,%f substPhiG:%f,%f\n",P,substPhiL[0],substPhiL[1],substPhiG[0],substPhiG[1]);
        while (dxTotal > 0.001)
        {
            FF_MixVfromTPeos(eosType,rule,T,&P,numSubs,data,pintParam,x,&option,answerL,answerG,&state);
            FF_MixFugacityEOS(eosType,rule,T,&answerL[0],numSubs,data,pintParam,x,substPhiL);
            xTotal=0;
            for (i=0;i<*numSubs;i++)
            {
                newX[i]=y[i]*substPhiG[i]/substPhiL[i];
                xTotal=xTotal+newX[i];
            }
            for (i=0;i<*numSubs;i++) newX[i]=newX[i]/xTotal; //we normalize y in order to obtain molar fractions
            dxTotal=0;
            for (i=0;i<*numSubs;i++)
            {
                dxTotal=dxTotal+fabs(newX[i]-x[i]);
                x[i]=newX[i];
            }
            counter=counter+1;
            //printf("counter:%i dxTotal:%f\n",counter,dxTotal);
            if (counter>2000) return;
        }
        //P = P / xTotal;//* (1 + (1 - xTotal) * 0.06);
        //printf("counter:%i P:%f xTotal:%f\n",counter,P,xTotal);
    }while (fabs(xTotal - 1) > 0.001);
    *dP=P;
    //printf("DewP:%f n:%i\n",P,counter);
}

//Pressure envelope of a binary mixture
void CALLCONV FF_PressureEnvelope(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const void *data,
                                  const double pintParam[],const int *numPoints,double c[],double x[],double y[],double bP[],double dP[]){
    double pGuess=0;
    int numSubs=2;
    int i;
    float interval;
    double Vp;
    double in[2],out[2];
    double substPhiL[2];
    double substPhiG[2];
    interval=1.0 /(*numPoints - 1);
    for(i=0;i<*numPoints;i++){
        c[i]=i*interval;
    }
    FF_VpEOS(eosType,T,&data[0],&Vp);
    bP[*numPoints - 1]=dP[*numPoints - 1]=Vp;
    //printf("bP:%f\n",bP[*numPoints - 1]);
    for(i=*numPoints-2;i>-1;i--){
        in[0]=c[i];
        in[1]=1-c[i];
        FF_BubbleP(eosType,rule,T,&numSubs,data,pintParam,in,&pGuess,&bP[i],out,substPhiL,substPhiG);
        y[i]=out[0];
        FF_DewP(eosType,rule,T,&numSubs,data,pintParam,in,&pGuess,&dP[i],out,substPhiL,substPhiG);
        x[i]=out[0];
        //printf("c:%f bP:%f y:%f dP:%f x:%f\n",in[0],bP[i],y[i],dP[i],x[i]);
    }
}



//VL flash calculation, given T, P, feed composition, eos and mixing rule
void CALLCONV FF_VLflashTP(const enum FF_EOS *eos,const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const void *data,
                     const double pintParam[],const double f[],double x[],double y[],double *beta)
{
    //first we will try to avoid bubble and dew point calculation, using the wilson formula for K[i]
    //We recover Yc,Pc and w for each substance
    double Tc[*numSubs],Pc[*numSubs],w[*numSubs];
    double k[*numSubs];
    char option,state;
    int i;
    double answerL[3],answerG[3],substPhiL[*numSubs],substPhiG[*numSubs];
    option='b';
    switch (*eos)
    {
        case FF_PCSAFT:
        case FF_DPCSAFT_GV:
        case FF_DPCSAFT_JC:
            for (i=0;i<*numSubs;i++)
            {
                Tc[i]=(( FF_SaftEOSdata*)data)[i].Tc;
                Pc[i]=(( FF_SaftEOSdata*)data)[i].Pc;
                w[i]=(( FF_SaftEOSdata*)data)[i].w;
            }
            break;
        case FF_SW:
            for (i=0;i<*numSubs;i++)
            {
                Tc[i]=(( FF_SWEOSdata*)data)[i].Tc;
                Pc[i]=(( FF_SWEOSdata*)data)[i].Pc;
                w[i]=(( FF_SWEOSdata*)data)[i].w;
            }
            break;
        case IF97:
            break;
        default://Cubic eos
            for (i=0;i<*numSubs;i++)
            {
                Tc[i]=(( FF_CubicEOSdata*)data)[i].Tc;
                Pc[i]=(( FF_CubicEOSdata*)data)[i].Pc;
                w[i]=(( FF_CubicEOSdata*)data)[i].w;
            }
            break;
    }
    for (i=0;i<*numSubs;i++)
    {
        k[i]=exp(log(Pc[i]/ *P)+5.373*(1-w[i])*(1-Tc[i]/ *T));
        //printf("k:%f\n",k[i]);
    }
    FF_MixVfromTPeos(eos,rule,T,P,numSubs,data,pintParam,f,&option,answerL,answerG,&state);
    //printf("state:%c\n",state);
    FF_MixFugacityEOS(eos,rule,T,&answerL[0],numSubs,data,pintParam,f,substPhiL);
    FF_MixFugacityEOS(eos,rule,T,&answerG[0],numSubs,data,pintParam,f,substPhiG);
    for (i=0;i<*numSubs;i++)
    {
        k[i]=substPhiL[i]/substPhiG[i];
        //printf("k:%f\n",k[i]);
    }
}






/*
def calcIsothermFlash(self,T,P,composition):
    bubbleT=self.calcBubbleT(P,composition,initT=T)[0]
    dewT=self.calcDewT(P,composition,initT=(bubbleT+50))[0]
    if T>dewT:
        lPhase=Phase([0.0],0.0)
        gPhase=Phase(composition,1.0)
    elif T<bubbleT:
        lPhase=Phase(composition,1.0)
        gPhase=Phase([0.0],0.0)
    else: #if T between bubble and dew points it is necessary to calculate the equilibrium
        substPhi=[[],[]]
        x=[1*composition[i] for i in range(len(composition))]
        y=[1*composition[i] for i in range(len(composition))]
        if self.fLmodel==1: #fugacity by EOS calculation
            if self.eosL==4 or self.eosL==5: #SAFT based EOS
                answer=self.calcMixPropSAFT(T,P,self.eosL,composition)
                substPhi=[answer[0][2],answer[1][2]]
            else:
                answer=self.calcMixPropCubicEOS(T,P,self.eosL,self.mixRulL,composition)
                substPhi=substPhi=[answer[0][2],answer[1][2]]
        elif self.fLmodel==0: #fugacity from activity
            substPhi[0]=[self.calcFugacityFromActivity(T,P,self.aLmodel,self.fLref,self.eosL,composition)[i]
                            /composition[i]/P for i in range(len(composition))]
            substPhi[1]=self.calcMixPropCubicEOS(T,P,self.eosG,self.mixRulG,composition)[1][2]
        k=[substPhi[0][i]/substPhi[1][i] for i in range(len(substPhi[0]))]#k first aproximation
        B=sum([composition[i]*k[i]/(k[i]+1.0) for i in range(len(k))])#gas fraction first aproximation
        counter=0
        for i in range(len(k)):# first aprox. for compositions
            x[i]=composition[i]/(1+B*(k[i]-1))
            y[i]=x[i]*k[i]
        yTotal=sum([y[i] for i in range(len(y))])#A sumatory different of 1 indicates a incorrect gas fraction
        while abs(yTotal-1)>0.0001 and B>0.0001 and B<0.9999: #Loop to obtain a summatory of gas fractions=1
            B=B*yTotal
            dK=1
            while dK>0.001: #Loop to obtain a stabilized composition for the selected gas fraction
                xTotal=sum([x[i] for i in range(len(x))])
                yTotal=sum([y[i] for i in range(len(y))])
                for i in range(len(x)):
                    x[i]=x[i]/xTotal  #we normalize in order to obtain molar fractions
                    y[i]=y[i]/yTotal
                if self.fLmodel==1: #fugacity by EOS calculation
                    if self.eosL==4 or self.eosL==5: #SAFT based EOS
                        substPhi[0]=self.calcMixPropSAFT(T,P,self.eosL,x)[0][2]
                    else:
                        #print T,P,self.eosL,self.mixRulL,x
                        substPhi[0]=self.calcMixPropCubicEOS(T,P,self.eosL,self.mixRulL,x)[0][2]

                    if self.eosG==0: #Ideal gas
                        substPhi[1]=[1.0 for i in range(len(composition))]
                    elif self.eosG==4 or self.eosG==5: #SAFT type EOS
                        substPhi[1]=self.calcMixPropSAFT(T,P,self.eosG,y,optLiquid=0,optGas=1)[1][2]
                    else: #Cubic EOS
                        substPhi[1]=self.calcMixPropCubicEOS(T,P,self.eosG,self.mixRulG,y,optLiquid=1)[1][2]

                    #substPhi[1]=self.calcMixPropCubicEOS(T,P,self.eosG,self.mixRulG,y)[1][2]
                elif self.fLmodel==0: #fugacity from activity
                    substPhi[0]=[self.calcFugacityFromActivity(T,P,self.aLmodel,self.fLref,self.eosL,x)[i]
                                /x[i]/P for i in range(len(composition))]
                    substPhi[1]=self.calcMixPropCubicEOS(T,P,self.eosG,self.mixRulG,y)[1][2]
                newK=[substPhi[0][i]/substPhi[1][i] for i in range(len(substPhi[0]))]
                dK=sum([abs(k[i]-newK[i]) for i in range(len(k))])
                k=[1*newK[i] for i in range(len(k))]
                for i in range(len(k)):
                    x[i]=composition[i]/(1+B*(k[i]-1))
                    y[i]=x[i]*k[i]
                counter=counter+1
                #print counter
                if counter>2000:
                    lPhase=Phase([0.0,0.0],0.0)
                    gPhase=Phase([0.0,0.0],0.0)
                    return [lPhase,gPhase,1]
            yTotal=sum([y[i] for i in range(len(y))])
        suma=[y[i]*B+x[i]*(1-B) for i in range(len(x))]
        #print suma
        lPhase=Phase(x,(1-B))
        gPhase=Phase(y,B)
    return [lPhase,gPhase,0]*/
