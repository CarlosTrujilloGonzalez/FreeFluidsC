/*

 * FFtools.c
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

// contains parameter fitting calculations for pure substances
//============================================================
#include <stdio.h>
#include <math.h>
#include <nlopt.h>
#include "FFbasic.h"
#include "FFtools.h"

//Determines the error of a correlation, using the given coefficients, and the real value supplied inside functionData
double CALLCONV FF_CorrelationError(unsigned nVar, const double coef[], double grad[], const  FF_CorrelationData *data)
//this is the function that gives the result in return, and that gives also the gradients
//In coef receives the coefficients to optimize, in grad there is the answer of the gradients, and in data there are the id of the correlation to use and the real data used
//for the optimization
{
    double result[data->nPoints];
    double error=0;
    int i;
    if (grad)
    {
        for (i=0;i<nVar;i++) grad[i]=0;
    }
    FF_CorrelationResult(&data->eq,coef,&data->nPoints,&data->x,result);
    for (i=0;i<data->nPoints;i++)
    {
        if (!(data->y[i]==0)) error=error+fabs((result[i]-data->y[i])/data->y[i]);
        //else if (!(result[i]==0)) error=error+fabs((result[i]-data->y[i])/result[i]);
    }
    return error/data->nPoints;
}


//Optimizer function for correlations. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
void CALLCONV FF_OptCorrelation(unsigned nVar,double lb[],double ub[],char enforceLimits[],const  FF_CorrelationData *data,double var[],double *error)
{
    nlopt_algorithm alg,algL;
    nlopt_opt opt,optL;
    unsigned i;
    double rLb[nVar],rUb[nVar],rVar[nVar];
    double varAlt[nVar],errorAlt;//alternative coefficients and error in different conditions
    for (i=0;i<nVar;i++){//we store the received values for limits and initial value
        if (enforceLimits[i]=='y'){
            rLb[i]=lb[i];
            rUb[i]=ub[i];
            rVar[i]=var[i];
        }
        //printf("%c %f %f %f\n",enforceLimits[i],rLb[i],rUb[i],rVar[i]);
    }
    //Default values for bounds and local algorithm
    for (i=0;i<nVar;i++){
        lb[i]=-1e10;
        ub[i]=1e10;
        var[i]=0;
    }

    algL=NLOPT_LN_NELDERMEAD;
    //algL=NLOPT_LN_SBPLX;
    //We prepare the bounds and initial values, and the optimization algorithm to use, according to the intended use
    switch (data->eq)
    {
    case FF_DIPPR100:
        for (i=0;i<nVar;i++){
            lb[i]=-1e15;
            ub[i]=1e15;
            var[i]=0;
        }
        break;
    case FF_DIPPR100Ld:
        lb[0]=200;
        ub[0]=3000;
        var[0]=800;
        lb[1]=-1e6;
        ub[1]=1e6;
        lb[2]=-1e4;
        ub[2]=1e4;
        lb[3]=-1e2;
        ub[3]=1e2;
        lb[4]=-1;
        ub[4]=1;
        for (i=1;i<nVar;i++) var[i]=1;
        break;
    case FF_Polynomial2:
        lb[0]=-1e3;
        ub[0]=1e3;
        lb[1]=-1e2;
        ub[1]=1e2;
        lb[2]=-1e1;
        ub[2]=1e1;
        lb[3]=-1;
        ub[3]=1;
        lb[4]=-1e-1;
        ub[4]=1e-1;
        break;
    case FF_DIPPR101Lv:
        lb[0]=-5e2;
        ub[0]=1e1;
        var[0]=1;
        lb[1]=-1e4;
        ub[1]=1e4;
        var[1]=1e3;
        lb[2]=-1e1;
        ub[2]=1e2;
        var[2]=1e1;
        lb[3]=-1e-4;
        ub[3]=1e-4;
        var[3]=-2e-6;
        lb[4]=2;
        ub[4]=2;
        var[4]=2;
        break;
    case FF_DIPPR101Vp:
        lb[0]=1e1;
        ub[0]=2e2;
        var[0]=63.9;
        lb[1]=-1.5e4;
        ub[1]=-4e3;
        var[1]=-4e3;
        lb[2]=-1e2;
        ub[2]=-1e0;
        var[2]=-10.34;
        lb[3]=-1e-4;
        ub[3]=1e-4;
        var[3]=4.16e-6;
        lb[4]=0;
        ub[4]=6;
        var[4]=2;
        break;
    case FF_DIPPR102:
        lb[0]=0;
        ub[0]=1e-4;
        lb[1]=0;
        ub[1]=1e1;
        lb[2]=0;
        ub[2]=1e4;
        break;
    case FF_DIPPR106Hv://It seems that LsurfT requires very narrow limits. For Hv changing back +-30 to +-10 seems better. For Ldens better to compare
        lb[0]=0;
        ub[0]=1e8;
        var[0]=1e7;
        lb[1]=-10;
        ub[1]=10;
        var[1]=0;
        lb[2]=-30;
        ub[2]=30;
        var[2]=0;
        lb[3]=-30;
        ub[3]=30;
        var[3]=0;
        lb[4]=-10;
        ub[4]=10;
        var[4]=0;
        break;
    case FF_DIPPR106Ld:
        lb[0]=0;
        ub[0]=3000;
        var[0]=500;
        lb[1]=-10;
        ub[1]=10;
        var[1]=0;
        lb[2]=-10;
        ub[2]=10;
        var[2]=0;
        lb[3]=-10;
        ub[3]=10;
        var[3]=0;
        lb[4]=-10;
        ub[4]=10;
        var[4]=0;
        break;
    case FF_DIPPR106SurfT://It seems that LsurfT requires very narrow limits. For Hv changing back +-30 to +-10 seems better. For Ldens better to compare
        lb[0]=0;
        ub[0]=1e0;
        var[0]=0.01;
        lb[1]=-1e2;
        ub[1]=1e2;
        var[1]=0;
        lb[2]=-1e2;
        ub[2]=1e2;
        var[2]=0;
        lb[3]=-1e2;
        ub[3]=1e2;
        var[3]=0;
        lb[4]=-1e2;
        ub[4]=1e2;
        var[4]=0;
        break;
    case FF_DIPPR107Cp:
        for (i=0;i<2;i++)
        {
            lb[i]=0;
            ub[i]=1e5;
            var[i]=1e4;
        }
        /*lb[0]=1e2;
        ub[0]=3.3e5;
        var[0]=5e3;
        lb[1]=0;
        ub[1]=1e5;
        var[1]=5e3;*/
        lb[2]=300;
        ub[2]=3.1e3;
        var[2]=1000;
        lb[3]=-1e8;
        ub[3]=1e5;
        var[3]=5e3;
        lb[4]=20;
        ub[4]=3100;
        var[4]=700;
        algL=NLOPT_LN_SBPLX;
        break;
    case FF_DIPPR116Ld://For liquid density
    case FF_PPDS10:
        lb[0]=0;//coef0 is not used, instead rho critical is used
        ub[0]=0;
        var[0]=0;
        lb[1]=-1e4;
        ub[1]=1e4;
        var[1]=0;
        lb[2]=-1e4;
        ub[2]=1e4;
        var[2]=0;
        lb[3]=-1e4;
        ub[3]=1e4;
        var[3]=0;
        lb[4]=-1e4;
        ub[4]=1e4;
        var[4]=0;
        break;
    case FF_Antoine1:
    case FF_Antoine2:
        lb[0]=0;
        ub[0]=5e1;
        var[0]=10.0;
        lb[1]=0;
        ub[1]=1e4;
        var[1]=4e3;
        lb[2]=-1e2;
        ub[2]=1e3;
        var[2]=2e2;
        break;
    case FF_Wagner25://For Vp
    case FF_Wagner36:
        //coef0 is forced to Pc, and coef5 to Tc, when recovering bounds
        for (i=1;i<nVar;i++){
            lb[i]=-1e2;
            ub[i]=1e2;
            var[i]=0;
        }
        break;
    case FF_Rackett://For liquid density
        lb[0]=0;
        ub[0]=1e3;
        var[0]=3e2;
        lb[1]=0;
        ub[1]=1;
        var[1]=0.2;
        lb[2]=1e2;
        ub[2]=1e3;
        var[2]=4e2;
        lb[3]=0;
        ub[3]=1;
        var[3]=0.2;
        break;
    case FF_WagnerGd://For gas saturated density
        //coef0 is forced to Rhoc, and coef5 to Tc, when recovering bounds
        for (i=1;i<nVar-3;i++){
            lb[i]=-10;
            ub[i]=0;
            var[i]=-2;
            lb[3]=-50;
            ub[3]=0;
            var[3]=-25;
            lb[4]=-1e2;
            ub[4]=0;
            var[4]=-50;
        }
        break;
    }

    for (i=0;i<nVar;i++){//we recover the received values for limits and initial value,if needed
        if (enforceLimits[i]=='y'){
            lb[i]=rLb[i];
            ub[i]=rUb[i];
            var[i]=rVar[i];
        }
    }
    alg=NLOPT_GN_MLSL;
    opt=nlopt_create(alg,nVar);
    optL=nlopt_create(algL,nVar);
    nlopt_set_ftol_rel(optL, 0.001);
    nlopt_set_xtol_rel(optL, 0.00001);
    //nlopt_set_population(opt,16);
    nlopt_set_local_optimizer(opt, optL);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_maxeval(opt,2000000);//number of evaluations
    nlopt_set_maxtime(opt,20.0);//time in seconds
    nlopt_set_min_objective(opt,FF_CorrelationError , data);
    int code=nlopt_optimize(opt, var, error);
    if (*error>0.1){
        nlopt_set_maxtime(opt,5.0);//time in seconds
        //We change de default values for bounds and local algorithm
        for (i=0;i<nVar;i++){
            if (enforceLimits[i]!='y'){
                lb[i]=-1e6;
                ub[i]=1e6;
                var[i]=0;
            }
        }
        int code=nlopt_optimize(opt, varAlt, &errorAlt);
        if (errorAlt<*error){
            *error=errorAlt;
            for (i=0;i<nVar;i++){
                var[i]=varAlt[i];
            }
        }
    }
    printf("Code:%i\n",code);
    if (code < 0)
    {
        printf("nlopt failed!\n");

        //printf("found minimum after %d evaluations\n", count);
        printf("found minimum at f(%g,%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],var[3],var[4],*error);
    }
    else
    {
        //printf("found minimum after %d evaluations\n", count);
        printf("found minimum at A:%g, B:%g, C:%g, D:%g, E:%g, Error:%0.10g\n", var[0], var[1], var[2],var[3],var[4],*error);
    }
    //nlopt_destroy(opt);
}

//Determines the error of a SAFT EOS, using the given coefficients, and the real value supplied inside data
double CALLCONV FF_SaftPvRhoError(unsigned nVar, const double coef[], double grad[],  FF_EOSPvRhoData *data)
{
    //printf("eos:%u \n",data->eos);
    static double minError=2;//used to actualize partial errors when there is an improvement in the total error found
    static unsigned int eosUsed=9999;//used to actualize partial errors data when changing the eos
     FF_SaftEOSdata saftData;
    enum FF_EosType eos=FF_SAFTtype;
    saftData.MW=data->MW;
    saftData.Tc=data->Tc;
    saftData.Pc=data->Pc;
    saftData.Zc=data->Zc;
    saftData.w=data->w;
    saftData.sigma=coef[0];
    saftData.m=coef[1];
    saftData.epsilon=coef[2];
    switch (data->eos){
    case FF_PCSAFT:
        saftData.kAB=0;
        saftData.epsilonAB=0;
        saftData.nPos=0;
        saftData.nNeg=0;
        saftData.nAcid=0;
        //printf("sigma:%f m:%f epsilon:%f\n",saftData.sigma,saftData.m,saftData.epsilon);
        break;
    case FF_PCSAFT1:
        saftData.kAB=coef[3];
        saftData.epsilonAB=coef[4];
        saftData.nPos=0;
        saftData.nNeg=0;
        saftData.nAcid=1;
        //printf("sigma:%f m:%f epsilon:%f Kab:%f epsilonAB:%f\n",saftData.sigma,saftData.m,saftData.epsilon,saftData.kAB,saftData.epsilonAB);
        break;
    case FF_PCSAFT2B:
        saftData.kAB=coef[3];
        saftData.epsilonAB=coef[4];
        saftData.nPos=1;
        saftData.nNeg=1;
        saftData.nAcid=0;
        //printf("PCSAFT 2B");
        break;
    case FF_PCSAFT4C:
        //printf("sigma:%f m:%f epsilon:%f\n",saftData.sigma,saftData.m,saftData.epsilon);
        saftData.kAB=coef[3];
        saftData.epsilonAB=coef[4];
        saftData.nPos=2;
        saftData.nNeg=2;
        saftData.nAcid=0;
        //printf("PCSAFT 4C");
        break;
    }
    //printf("MW:%f Tc:%f Pc:%f Zc:%f\n",saftData.MW,saftData.Tc,saftData.Pc,saftData.Zc,saftData.nAcid,saftData.nPos,saftData.nNeg);
    //printf("sigma:%f m:%f epsilon:%f Kab:%f epsilonAB:%f\n",saftData.sigma,saftData.m,saftData.epsilon,saftData.kAB,saftData.epsilonAB);



    unsigned i;
    //Zc calculated value constraint
    double ZcDiff=0;

    double Vc,Arr,Z;//Vc is the critical volume calculated from Tc,Pc and Zc. Arr and Z are the saft calculated values at Tc,Vc point
    Vc=saftData.Zc*R*saftData.Tc/saftData.Pc;
    FF_ArrZfromTVPCSAFT(&saftData.Tc,&Vc,&saftData,&Arr,&Z);
    ZcDiff=fabs((Z-saftData.Zc)/saftData.Zc);
    //printf("Vc:%f Zc:%f Zsaft:%f Zcdiff:%f\n",Vc,saftData.Zc,Z,ZcDiff);
    if (!(Z>0)) return +HUGE_VALF;
    if (ZcDiff>data->zcFilter) return 4.0;


    //Liquid density value constraint
    char option='l',state;
    double answerL[3],answerG[3];
    double densDiff=0,error;
    for (i=0;i<data->nPoints;i++){
        FF_VfromTPPCSAFT(&data->points[i][0],&data->points[i][1],&saftData,&option,answerL,answerG,&state);
        densDiff=densDiff+fabs((saftData.MW*1e-3/answerL[0]-data->points[i][2])/data->points[i][2]);
        //printf("T:%f d:%f Vfound:%f Dfound:%f\n",data->points[i][0],data->points[i][2],answerL[0],saftData.MW*1e-3/answerL[0]);
        //printf("densDiff:%f\n",densDiff);
    }
    densDiff=densDiff/data->nPoints;
    //printf("densDiff:%f\n",densDiff);
    if (densDiff>data->ldensFilter) return 2.0;

    //If Zc and liquid density OK we conduct the optimization
    double Vp,VpDiff=0;
    for (i=0;i<data->nPoints;i++)
    {
        FF_VpEOS(&eos,&data->points[i][0],&saftData,&Vp);
        VpDiff=VpDiff+fabs((Vp-data->points[i][1])/data->points[i][1]);
        //printf("Vp found:%f Vp in:%f\n",Vp,data->points[i][1]);
    }
    VpDiff=VpDiff/data->nPoints;
    error=0.5*VpDiff+0.5*densDiff;
    if ((minError>error)||!(eosUsed==data->eos)){
        minError=error;
        eosUsed=data->eos;
        data->vpError=VpDiff;
        data->ldensError=densDiff;
        data->zcError=ZcDiff;
        //printf("sigma:%f m:%f epsilon:%f Kab:%f epsilonAB:%f\n",saftData.sigma,saftData.m,saftData.epsilon,saftData.kAB,saftData.epsilonAB);
        //printf("ZcDiff:%f desDiff;%f VpDiff:%f\n",ZcDiff,densDiff,VpDiff);
    }
    //printf("Error:%f\n",error/data->nPoints);
    return error;
}


//Determines the error of a cubic EOS, using the given coefficients, and the real value supplied inside data
EXP_IMP double CALLCONV FF_CubicPvRhoError(unsigned nVar, const double coef[], double grad[],  FF_EOSPvRhoData *data)
{
    unsigned i;
    double Vp,VpDiff=0,error=0;
    FF_CubicEOSdata cubicData;
    FF_CubicParam param;
    enum FF_EosType eos=FF_CubicType;
    cubicData.eos=data->eos;
    cubicData.MW=data->MW;
    cubicData.Zc=data->Zc;
    switch (data->eos){
    case FF_PR76:
    case FF_PR78:
        cubicData.Tc=coef[0];
        cubicData.Pc=coef[1];
        cubicData.w=coef[2];
        cubicData.c=0;
        cubicData.k1=0;
        cubicData.k2=0;
        cubicData.k3=0;
        cubicData.k4=0;
        break;
    case FF_PRFIT3:
        cubicData.Tc=data->Tc;
        cubicData.Pc=data->Pc;
        cubicData.w=data->w;
        cubicData.c=0;
        cubicData.k1=coef[0];
        cubicData.k2=coef[1];
        cubicData.k3=coef[2];
        cubicData.k4=0;
        break;
    case FF_PRFIT4:
        cubicData.Tc=data->Tc;
        cubicData.Pc=data->Pc;
        cubicData.w=data->w;
        cubicData.c=-1;
        cubicData.k1=coef[0];
        cubicData.k2=coef[1];
        cubicData.k3=coef[2];
        cubicData.k4=coef[3];
        break;
    default:
        cubicData.Tc=data->Tc;
        cubicData.Pc=data->Pc;
        cubicData.w=data->w;
        cubicData.VdWV=data->VdWV;
        cubicData.c=0;
        cubicData.k1=coef[0];
        cubicData.k2=coef[1];
        cubicData.k3=coef[2];
        cubicData.k4=0;
        //printf("k1:%f k2:%f k3:%f\n",coef[0],coef[1],coef[2]);
        break;
    }

    //Liquid density value constraint
    char option='l',state;
    double answerL[3],answerG[3];
    double densDiff=0;
    if ((data->eos==FF_PR78)||(data->eos==FF_PRFIT3)||(data->eos==FF_PRFIT4)){
        for (i=0;i<data->nPoints;i++){
            FF_FixedParamCubic(&cubicData,&param);
            FF_ThetaDerivCubic(&data->points[i][0],&cubicData,&param);
            FF_VfromTPcubic(&data->points[i][0],&data->points[i][1],&param,&option,answerL,answerG,&state);
            densDiff=densDiff+fabs((cubicData.MW*1e-3/answerL[0]-data->points[i][2])/data->points[i][2]);
        }
        densDiff=densDiff/data->nPoints;
        if (densDiff>data->ldensFilter) return 2.0;//if density error is higher than the constraint we return a 200% error in order to reject the coefficients
    }
    static double minError=2;//used to actualize partial errors when there is an improvement in the total error found
    static unsigned int eosUsed=9999;//used to actualize partial errors when changing the used eos

    for (i=0;i<data->nPoints;i++){
        FF_VpEOS(&eos,&data->points[i][0],&cubicData,&Vp);
        VpDiff=VpDiff+fabs((Vp-data->points[i][1])/data->points[i][1]);
    }
    VpDiff=VpDiff/data->nPoints;
    if ((data->eos==FF_PR78)||(data->eos==FF_PRFIT3)||(data->eos==FF_PRFIT4)) error=0.6*VpDiff+0.4*densDiff;
    else error=VpDiff;
    if (isnan(error)) return 4.0;
    if ((minError>error)||!(eosUsed==data->eos)){
        minError=error;
        eosUsed=data->eos;
        data->vpError=VpDiff;
        if ((data->eos==FF_PR78)||(data->eos==FF_PRFIT3)||(data->eos==FF_PRFIT4)) data->ldensError=densDiff;
        else data->ldensError=+HUGE_VALF;
        switch (data->eos){//We print the new optimal values found
        case FF_PR76:
        case FF_PR78:
            printf("Tc:%f Pc:%f w:%f\n",cubicData.Tc,cubicData.Pc,cubicData.w);
            break;
        case FF_PRFIT3:
            printf("a:%f Tc:%f k1:%f\n",cubicData.k1,cubicData.k2,cubicData.k3);
            break;
        case FF_PRFIT4:
            printf("b:%f a:%f w:%f Tc:%f\n",cubicData.k1,cubicData.k2,cubicData.k3,cubicData.k4);
            break;
        case FF_PRSV1:
        case FF_PRBM:
            printf("k1:%f\n",cubicData.k1);
            break;
        case FF_PRMELHEM:
        case FF_PRSOF:
        case FF_SRKSOF:
            printf("k1:%f k2:%f\n",cubicData.k1,cubicData.k2);
            break;
        default:
            printf("k1:%f k2:%f k3:%f\n",cubicData.k1,cubicData.k2,cubicData.k3);
            break;
        }
        printf("densDiff:%f VpDiff:%f\n",densDiff,VpDiff);
    }
    return error;
}


//Constraints the parameters of a cubic EOS to a generated Apha= straight EOS Alpha +-5%
EXP_IMP double CALLCONV FF_CubicEOSAlphaConstraint(unsigned nVar, const double coef[], double grad[], const  FF_EOSPvRhoData *data)
{
    double Tr,Tx,alphaBase,alphaDiff=0;
    int i;
    double diff=0;
    for (i=0;i<data->nPoints;i++)
    {
        Tr=data->points[i][0]/data->Tc;
        Tx = 1-pow(Tr,0.5);
        alphaBase=pow(1+(0.37464+1.54226*data->w-0.26992*pow(data->w,2))*Tx,2);
        switch (data->eos){
        case FF_PRSOF:
            alphaBase=pow(1+(0.37464+1.54226*data->w-0.26992*pow(data->w,2))*Tx,2);
            alphaDiff=fabs((Tr*(1+coef[0]/(coef[1]-1)*(1-pow(Tr,(1-coef[1]))))-alphaBase)/alphaBase);
            break;
        case FF_PRMC:
            alphaBase=pow(1+(0.37464+1.54226*data->w-0.26992*pow(data->w,2))*Tx,2);
            alphaDiff=fabs((pow(1+coef[0]*Tx+coef[1]*pow(Tx,2)+coef[3]*pow(Tx,3),2)-alphaBase)/alphaBase);
            break;
        /*case FF_PRSV1:
            alphaBase=pow(1+(0.37464+1.54226*data->w-0.26992*pow(data->w,2))*Tx,2);
            alphaDiff=fabs((pow((1+(0.378893+1.4897153*data->w-0.17131848*pow(data->w,2)+0.0196554*pow(data->w,3))*Tx+coef[0]*(1-Tr)*(0.7-Tr)),2)-alphaBase)/alphaBase);
            break;*/
        case FF_PRTWU91:
            alphaBase=pow(1+(0.37464+1.54226*data->w-0.26992*pow(data->w,2))*Tx,2);
            alphaDiff=fabs((pow(Tr,coef[2]*(coef[1]-1))*exp(coef[0]*(1-pow(Tr,(coef[1]*coef[2]))))-alphaBase)/alphaBase);
            break;
        //case FF_PRFIT4:
        //    alphaBase=pow(1+(0.37464+1.54226*data->w-0.26992*pow(data->w,2))*Tx,2);
        //    alphaDiff=fabs((pow(1+coef[2]*(1-pow(data->points[i][0]/coef[3],0.5)),2)-alphaBase)/alphaBase);
        //    break;
        case FF_SRKTWU91:
            alphaBase=pow(1+(0.48 + 1.574 * data->w - 0.176 * pow(data->w,2))*Tx,2);
            alphaDiff=fabs((pow(Tr,coef[2]*(coef[1]-1))*exp(coef[0]*(1-pow(Tr,(coef[1]*coef[2]))))-alphaBase)/alphaBase);
            break;
        }
        if (alphaDiff > diff) diff=alphaDiff;

    }
    //printf("Cubic constraint diff:%f\n",diff);
    return diff-0.05;
}

//Optimizer function for EOS parameters. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
void CALLCONV FF_OptEOSparam(unsigned optTime,unsigned nVar,double lb[],double ub[],char enforceLimits[], FF_EOSPvRhoData*data,double var[],double *error)
{
    double rLb[nVar],rUb[nVar],rVar[nVar];
    nlopt_algorithm alg,algL;
    nlopt_opt opt,optL;
    int i;
    for (i=0;i<nVar;i++){//we store the received values for limits and initial value
        if (enforceLimits[i]=='y'){
            rLb[i]=lb[i];
            rUb[i]=ub[i];
            rVar[i]=var[i];
        }
    }
    switch (data->eos){
    case FF_PR76://optimal Tc,Pc and w will be calculated
    case FF_PR78:
        nVar=3;
        lb[0]=0.9*data->Tc;//will replace Tc
        ub[0]=1.1*data->Tc;
        lb[1]=0.9*data->Pc;//will replace Pc
        ub[1]=1.1*data->Pc;
        lb[2]=0.9*data->w;//will replace w
        ub[2]=1.1*data->w;
        var[0]=data->Tc;
        var[1]=data->Pc;
        var[2]=data->w;
        break;
    case FF_PRSV1:
        nVar=1;
        lb[0]=-1;
        ub[0]=1;
        var[0]=0;
        break;
    case FF_PRBM:
        nVar=1;
        lb[0]=-1;
        ub[0]=1;
        var[0]=0;
        break;
    case FF_PRMELHEM:
        nVar=2;
        lb[0]=0.3;
        ub[0]=2;
        lb[1]=-11;
        ub[1]=1.5;
        var[0]=0.7;
        var[1]=0.3;
        break;
    case FF_PRSOF:
    case FF_SRKSOF:
        nVar=2;
        lb[0]=-4;
        ub[0]=-1;
        lb[1]=1.9;
        ub[1]=3;
        var[0]=-2;
        var[1]=2.5;
        break;
    case FF_PRALMEIDA:
        nVar=3;
        lb[0]=0.2;
        ub[0]=1.5;
        lb[1]=-0.15;
        ub[1]=0.5;
        lb[2]=0.3;
        ub[2]=5;
        var[0]=0.7;
        var[1]=0.15;
        var[2]=0.6;
        break;
    case FF_PRMC:
    case FF_SRKMC:
        nVar=3;
        lb[0]=0.3;
        ub[0]=1.6;
        lb[1]=-2;
        ub[1]=4;
        lb[2]=-7;
        ub[2]=6;
        var[0]=1;
        var[1]=1;
        var[2]=0;
        break;
    case FF_PRTWU91:
    case FF_SRKTWU91:
        nVar=3;
        lb[0]=0.1;
        ub[0]=2.2;
        lb[1]=0.8;
        ub[1]=1.5;
        lb[2]=0.3;
        ub[2]=3;
        var[0]=0.6;
        var[1]=0.9;
        var[2]=1.6;
        break;
    case FF_PRFIT3://will optimize a,Tc and k1, being k1 the parameter used in the Peng-Robinson Stryjeck-Vera Theta calculation
        nVar=3;
        lb[0]=0.7*0.45724 * pow(R*data->Tc,2)/ data->Pc;;//equivalent to a
        ub[0]=1.3 * pow(R*data->Tc,2)/ data->Pc;
        lb[1]=0.9*data->Tc;//equivalent to Tc
        ub[1]=1.1*data->Tc;
        lb[2]=-1;//equivalent to k1 of FF_PRSV1
        ub[2]=1;
        var[0]=0.45724 * pow(R*data->Tc,2)/ data->Pc;
        var[1]=data->Tc;
        var[2]=0;
        break;
    case FF_PRFIT4://Will optimize a,b,Tc and w regarding liquid density and vapor pressure
        nVar=4;
        lb[0]=0.7*0.0778 * R * data->Tc / data->Pc;//equivalent to b
        ub[0]=1.2*0.0778 * R * data->Tc / data->Pc;
        lb[1]=0.7*0.45724 * pow(R*data->Tc,2)/ data->Pc;//equivalent to a
        ub[1]=1.3*0.45724 * pow(R*data->Tc,2)/ data->Pc;
        lb[2]=0.8*data->w;//equivalent to w
        ub[2]=1.2*data->w;
        lb[3]=0.9*data->Tc;//equivalent to Tc
        ub[3]=1.1*data->Tc;
        var[0]=0.85*0.0778 * R * data->Tc / data->Pc;//equivalent to b
        var[1]=0.45724 * pow(R*data->Tc,2)/ data->Pc;//equivalent to a
        var[2]=data->w;//equivalent to w
        var[3]=data->Tc;//equivalent to Tc
        break;
    case FF_PRFIT4B:
        nVar=4;
        lb[0]=-1;
        ub[0]=1;
        lb[1]=0.95*data->Tc;
        ub[1]=1.05*data->Tc;
        lb[2]=0.95*data->Pc;
        ub[2]=1.05*data->Pc;
        lb[3]=0.95*data->w;
        ub[3]=1.05*data->w;
        var[0]=0;
        var[1]=data->Tc;
        var[2]=data->Pc;
        var[3]=data->w;
        break;
    case FF_PRvTWU91:
        nVar=3;
        lb[0]=0.1;
        ub[0]=2.2;
        lb[1]=0.8;
        ub[1]=1.5;
        lb[2]=0.3;
        ub[2]=3;
        //lb[3]=-0.2*data->VdWV;
        //ub[3]=0.2*data->VdWV;
        var[0]=0.6;
        var[1]=0.9;
        var[2]=1.6;
        //var[3]=0;
        break;
    case FF_PCSAFT://use directly Nelder-Mead algorithm
        nVar=3;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=5;
        lb[2]=150;
        ub[2]=350;
        var[0]=2.5;
        var[1]=1;
        var[2]=350;
        break;
    case FF_PCSAFT1://PCSAFT 1A assoc. schema Use MLSL_LDS plus Nealder-Mead
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=5;
        lb[2]=150;
        ub[2]=350;
        lb[3]=0.001;
        ub[3]=0.5;
        lb[4]=1000;
        ub[4]=5000;
        var[0]=3;
        var[1]=1.5;
        var[2]=250;
        var[3]=0.001;
        var[4]=2300;
        break;
    case FF_PCSAFT2B://Use MLSL_LDS plus Nealder-Mead
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=5;
        lb[2]=150;
        ub[2]=350;
        lb[3]=0.005;
        ub[3]=0.5;
        lb[4]=1000;
        ub[4]=2300;
        var[0]=3;
        var[1]=1.5;
        var[2]=250;
        var[3]=0.02;
        var[4]=2300;
        break;
    case FF_PCSAFT4C://Use MLSL_LDS plus Nealder-Mead
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=5;
        lb[2]=100;
        ub[2]=400;
        lb[3]=0.005;
        ub[3]=0.5;
        lb[4]=1000;
        ub[4]=3200;
        var[0]=2.0;
        var[1]=3.5;
        var[2]=300.0;
        var[3]=0.02;
        var[4]=2300;
        break;
    }

    for (i=0;i<nVar;i++){//we recover the received values for limits and initial value,if needed
        if (enforceLimits[i]=='y'){
            lb[i]=rLb[i];
            ub[i]=rUb[i];
            var[i]=rVar[i];
        }
    }

    switch (data->eos){
    case FF_PCSAFT:
    case FF_PCSAFT1:
    case FF_PCSAFT2B:
    case FF_PCSAFT4C:
        alg=NLOPT_GN_MLSL_LDS;
        algL=NLOPT_LN_NELDERMEAD;
        opt=nlopt_create(alg,nVar);
        optL=nlopt_create(algL,nVar);
        nlopt_set_ftol_rel(optL, 0.0001);
        nlopt_set_xtol_rel(optL, 0.0001);
        nlopt_set_local_optimizer(opt, optL);
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);
        nlopt_set_maxeval(opt,1e8);//number of evaluations
        nlopt_set_maxtime(opt,optTime);//Max.time in seconds
        nlopt_set_min_objective(opt, FF_SaftPvRhoError, data);
        break;
    default://cubic EOS
        alg=NLOPT_GN_MLSL_LDS;
        if (data->eos==FF_PR76 || data->eos==FF_PR78)algL=NLOPT_LN_SBPLX;
        else algL=NLOPT_LN_NELDERMEAD;
        opt=nlopt_create(alg,nVar);
        optL=nlopt_create(algL,nVar);
        nlopt_set_ftol_rel(optL, 0.0001);
        nlopt_set_xtol_rel(optL, 0.0001);
        nlopt_set_local_optimizer(opt, optL);
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);
        nlopt_set_maxeval(opt,1e6);//max. number of evaluations
        nlopt_set_maxtime(opt,optTime);//Max.time in seconds
        nlopt_set_min_objective(opt, FF_CubicPvRhoError, data);
        //nlopt_add_inequality_constraint(opt,FF_CubicEOSAlphaConstraint,data,1e-2);//this makes that the calculated alpha is not too different from PR or FF_SRK original alpha
        break;
    }
    int code=nlopt_optimize(opt, var, error);
    printf("Return code:%i\n",code);
    if (code < 0){
        printf("nlopt failed!\n");
    }
    else{
        switch(data->eos){
        case FF_PCSAFT1:
        case FF_PCSAFT2B:
        case FF_PCSAFT4C:
            printf("found minimum at f(%g,%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],var[3],var[4],*error);
            break;
        case FF_PRFIT4:
            printf("found minimum at f(%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],var[3],*error);
            break;
        case FF_PR76:
        case FF_PR78:
        case FF_PRALMEIDA:
        case FF_PRMC:
        case FF_PRTWU91:
        case FF_PRFIT3:
        case FF_SRKMC:
        case FF_SRKTWU91:
            printf("found minimum at f(%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],*error);
            break;
        case FF_PRMELHEM:
        case FF_PRSOF:
        case FF_SRKSOF:
            printf("found minimum at f(%g,%g) = %0.10g\n", var[0], var[1],*error);
            break;
        default:
           printf("found minimum at f(%g) = %0.10g\n", var[0],*error);
            break;
        }
    }

    //nlopt_destroy(opt);
}

//Other auxiliary functions
//-------------------------
//Determines the mass and molar fractions, given quatities in mass or moles
EXP_IMP double CALLCONV FF_FractionsCalculation(unsigned nSubs, const double MW[], const double q[], const bool mass, double massFrac[], double molarFrac[]){
    int i;
    double totalMass=0,totalMole=0;
    for (i=0;i<nSubs;i++){
        if (mass==true){
            totalMass=totalMass+q[i];
            totalMole=totalMole+q[i]/MW[i];
        }
        else{
            totalMass=totalMass+q[i]*MW[i];
            totalMole=totalMole+q[i];
        }
    }
    for (i=0;i<nSubs;i++){
        if (mass==true){
            massFrac[i]=q[i]/totalMass;
            molarFrac[i]=q[i]/MW[i]/totalMole;
        }
        else{
            massFrac[i]=q[i]*MW[i]/totalMass;
            molarFrac[i]=q[i]/totalMole;
        }
    }
}
