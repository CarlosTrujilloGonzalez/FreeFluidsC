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

// contains parameter fitting calculations for pure substances and mixtures
//=========================================================================
#include <stdio.h>
#include <math.h>
#include <nlopt.h>
#include <time.h>
#include<stdlib.h>
#include "FFbasic.h"
#include "FFtools.h"


//General Regula Falsi solver, Anderson-Bjork modification. Returns the solution. Uses auxiliary data for the function
double FF_solverRegulaBase(double y, double (*f)(double, void *), void *data, double a, double b, double ftol, int niter){//xmin and xmax define the interval
    double x=+HUGE_VALF;//The solution
    int n=0;
    int side=0;
    double ex,ea,eb;//errors at x, a and b points
    ea=(*f)(a,data)-y;
    eb=(*f)(b,data)-y;
    if (ea*eb>0){//if both functions have the same sign, it is not sure that a root exists
        double xmed,exmed;
        xmed=a+0.5*(b-a);//we try in the middle point
        exmed=(*f)(xmed,data)-y;
        if (ea*exmed<0){
            a=xmed;
            ea=exmed;
        }
        else return x;
    }
    if ((ea<ftol)&&(ea>(-ftol))) x=a;//if one of the given value was already the solution
    else if ((eb<ftol)&&(eb>(-ftol))) x=b;
    else{
        while (n<niter){
            x=a+(b-a)*ea/(ea-eb);
            ex=(*f)(x,data)-y;
            //printf("n:%i y:%f a:%f ea:%f b:%f eb:%f x:%f fx:%e ex:%f\n",n+1,y,a,ea,b,eb,x,ex+y,ex);
            if ((ex<ftol)&&(ex>(-ftol))) break;
            if ((eb*ex)>0){//x and b are at the same side
                if (side==1){//if it happened also in the previous loop
                    if ((ex/eb)<1) ea=ea-ea*(ex/eb);//we decrease y-f(a) for the next loop calculation
                    else ea*=0.5;
                }
                side = 1;
                b=x;
                eb=ex;
            }
            else{
                if (side==-1){
                    if ((ex/ea)<1) eb=eb-eb*(ex/ea);
                    else eb *= 0.5;
                }
                side = -1;
                a=x;
                ea=ex;
            }
            n=n+1;
        }
        if (n>=niter) x=HUGE_VALF;
    }
    return x;
}

//Diferential evolution optimizer
void CALLCONV FF_OptimizerDE(unsigned nVar,double lb[],double ub[],double (*f)(unsigned,const double *,double *,const void *),void *data,double var[],double *error){
    int ver=1;
    int i,j,k,l,m;
    int nPop;//number of population used
    int best;//index of the best solution
    float rm;//randon value between 0 and 1
    int rm1,rm2,rm3;//integer random values
    double errBest;
    double F,Cr;
    double auxG,auxComp[nVar];
    if (nVar<4) nPop=pow(10,nVar);
    else nPop=1200*nVar;
    double actPop[nPop][nVar],test[nVar],err[nPop],errTest;//actual points and new points, with its number of moles composition
    double grad[nVar];

    srand(time(NULL));
    F=0.3;//mutation scaling factot 0-1
    Cr=0.75;//crossover rate 0-1
    errBest=+HUGE_VALF;


    if(nVar==2){//Creation of an uniform grid of points
        for(i=0;i<10;i++){
            for(j=0;j<10;j++){
                    actPop[i*10+j][0]=lb[0]+i*(ub[0]-lb[0])/9;
                    actPop[i*10+j][1]=lb[1]+j*(ub[1]-lb[1])/9;
                }
            }
        }

    else if(nVar==3){//Creation of an uniform grid of points
        for(i=0;i<10;i++){
            for(j=0;j<10;j++){
                for(k=0;k<10;k++){
                    actPop[i*100+j*10+k][0]=lb[0]+i*(ub[0]-lb[0])/9;
                    actPop[i*100+j*10+k][1]=lb[1]+j*(ub[1]-lb[1])/9;
                    actPop[i*100+j*10+k][2]=lb[2]+k*(ub[1]-lb[2])/9;
                }
            }
        }
    }
    else{//Creation of a random grid of points
        for(i=0;i<nPop;i++){//Initial random population generation, with its mole content
            for(j=0;j<nVar;j++){
                rm=(float) rand()/RAND_MAX;//randon number between 0 and 1
                actPop[i][j]=lb[j]+rm*(ub[j]-lb[j]);
            }
        }
    }

    //Calculation of error for the initial population
    for(j=0;j<nPop;j++){
        err[j]=(*f)(nVar,actPop[j],grad,data);
        //printf("Initial population element:%i a:%f b:%f c:%f error:%f\n",j,actPop[j][0],actPop[j][1],actPop[j][2],err[j]);
        if(err[j]<errBest){
            best=j;
            errBest=err[j];
            for(i=0;i<nVar;i++){
                var[i]=actPop[j][i];//stores best variables
            }
        }
    }
    if (ver==1) printf("Initial population best element:%i error:%f\n",best,errBest);
    //Initialization of the evolution
    for(k=0;k<500;k++){
        for(j=0;j<nPop;j++){
            rm1=nPop*rand()/RAND_MAX;//randon number
            rm2=nPop*rand()/RAND_MAX;
            rm3=nPop*rand()/RAND_MAX;
            for(i=0;i<nVar;i++){//Crossover plus mutation
                rm=(float) rand()/RAND_MAX;//randon number between 0 and 1
                if(rm>Cr) test[i]=actPop[j][i];
                else{
                    test[i]=actPop[rm1][i]+F*(actPop[rm2][i]-actPop[rm3][i]);
                    if(test[i]<lb[i]) test[i]=lb[i];
                    else if(test[i]>ub[i]) test[i]=ub[i];
                }
            }
            //Calculation of error of the test
            errTest=(*f)(nVar,test,grad,data);

            if(errTest<err[j]){
                for(i=0;i<nVar;i++)actPop[j][i]=test[i];
                err[j]=errTest;
                if(err[j]<errBest){
                    best=j;
                    errBest=err[j];
                    for(i=0;i<nVar;i++){
                        var[i]=actPop[j][i];//stores best variables
                    }
                    if (ver==1) printf("generation:%i best element:%i error:%f\n",k,best,errBest);
                }
            }
        }
        if(nPop>50){
            if((k==10)||(k==20)||(k==50)||(k==100)||(k==150)||(k==200)||(k==250)){//Population reduction after n evolutions
                for(i=0;i<nPop;i++){
                    for(l=i+1;l<nPop;l++){
                        if(err[i]>err[l]){
                            auxG=err[i];
                            err[i]=err[l];
                            err[l]=auxG;
                            for(m=0;m<nVar;m++){
                                auxComp[m]=actPop[i][m];
                                actPop[i][m]=actPop[l][m];
                                actPop[l][m]=auxComp[m];
                            }
                        }
                    }
                }
                nPop=nPop/2;
            }
        }

    }
    *error=errBest;
}


//Recibe la ecuación a usar y establece los límites para las variables
void CALLCONV FF_CorrelationBounds(unsigned nVar, int eq, double lb[], double ub[]){
    unsigned i;
    switch (eq)
    {
    case FF_DIPPR100:
        for (i=0;i<nVar;i++){
            lb[i]=-1e15;
            ub[i]=1e15;
        }
        break;
    case FF_Polynomial:
        for (i=0;i<nVar;i++){
            lb[i]=-1e15;
            ub[i]=1e15;
        }
        break;
    case FF_DIPPR100Ld:
        lb[0]=200;
        ub[0]=3000;
        lb[1]=-1e6;
        ub[1]=1e6;
        lb[2]=-1e4;
        ub[2]=1e4;
        lb[3]=-1e2;
        ub[3]=1e2;
        lb[4]=-1;
        ub[4]=1;
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
        lb[1]=-1e4;
        ub[1]=1e4;
        lb[2]=-1e1;
        ub[2]=1e2;
        lb[3]=-1e-4;
        ub[3]=1e-4;
        lb[4]=2;
        ub[4]=2;
        break;
    case FF_DIPPR101Vp:
        lb[0]=1e1;
        ub[0]=2e2;
        lb[1]=-1.5e4;
        ub[1]=-4e3;
        lb[2]=-1e2;
        ub[2]=-1e0;
        lb[3]=-1e-4;
        ub[3]=1e-4;
        lb[4]=0;
        ub[4]=6;
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
        lb[1]=-10;
        ub[1]=10;
        lb[2]=-30;
        ub[2]=30;
        lb[3]=-30;
        ub[3]=30;
        lb[4]=-10;
        ub[4]=10;
        break;
    case FF_DIPPR106Ld:
        lb[0]=0;
        ub[0]=3000;
        lb[1]=-10;
        ub[1]=10;
        lb[2]=-10;
        ub[2]=10;
        lb[3]=-10;
        ub[3]=10;
        lb[4]=-10;
        ub[4]=10;
        break;
    case FF_DIPPR106SurfT://It seems that LsurfT requires very narrow limits. For Hv changing back +-30 to +-10 seems better. For Ldens better to compare
        lb[0]=0;
        ub[0]=1e0;
        lb[1]=-1e2;
        ub[1]=1e2;
        lb[2]=-1e2;
        ub[2]=1e2;
        lb[3]=-1e2;
        ub[3]=1e2;
        lb[4]=-1e2;
        ub[4]=1e2;
        break;
    case FF_DIPPR107Cp:
        for (i=0;i<2;i++)
        {
            lb[i]=0;
            ub[i]=1e5;
        }
        /*lb[0]=1e2;
        ub[0]=3.3e5;
        var[0]=5e3;
        lb[1]=0;
        ub[1]=1e5;
        var[1]=5e3;*/
        lb[2]=300;
        ub[2]=3.1e3;
        lb[3]=-1e8;
        ub[3]=1e5;
        lb[4]=20;
        ub[4]=3100;
        break;
    case FF_DIPPR116Ld://For liquid density
    case FF_PPDS10:
        lb[0]=0;//coef0 is not used, instead rho critical is used
        ub[0]=0;
        lb[1]=-1e4;
        ub[1]=1e4;
        lb[2]=-1e4;
        ub[2]=1e4;
        lb[3]=-1e4;
        ub[3]=1e4;
        lb[4]=-1e4;
        ub[4]=1e4;
        break;
    case FF_PPDS9://Liquid Viscosity
        lb[0]=-5.0;
        ub[0]=1e1;
        lb[1]=-1.0;
        ub[1]=1e1;
        lb[2]=1.0e2;
        ub[2]=2.0e3;
        lb[3]=-1.0e4;
        ub[3]=5.0e2;
        lb[4]=1.0e-7;
        ub[4]=1.0e-3;
        break;
    case FF_Antoine1:
    case FF_Antoine2:
        lb[0]=0;
        ub[0]=5e1;
        lb[1]=0;
        ub[1]=1e4;
        lb[2]=-1e2;
        ub[2]=1e3;
        break;
    case FF_Wagner25://For Vp
    case FF_Wagner36:
        //coef0 is forced to Pc, and coef5 to Tc, when recovering bounds
        for (i=1;i<nVar;i++){
            lb[i]=-1e2;
            ub[i]=1e2;
        }
        break;
    case FF_Rackett://For liquid density
        lb[0]=0;
        ub[0]=1e3;
        lb[1]=0;
        ub[1]=1;
        lb[2]=1e2;
        ub[2]=1e3;
        lb[3]=0;
        ub[3]=1;
        break;
    case FF_WagnerGd://For gas saturated density
        //coef0 is forced to Rhoc, and coef5 to Tc, when recovering bounds
        for (i=1;i<nVar-3;i++){
            lb[i]=-10;
            ub[i]=0;
        }
            lb[3]=-50;
            ub[3]=0;
            lb[4]=-1e2;
            ub[4]=0;
        break;
    default:
        for (i=0;i<nVar;i++){
            lb[i]=-1e10;
            ub[i]=1e10;
        }
    }
}


//Determines the error of a correlation, using the given coefficients, and the real value supplied inside functionData
double CALLCONV FF_CorrelationError(unsigned nVar, const double coef[], double grad[], const void *data1)
//this is the function that gives the result in return, and that gives also the gradients
//In coef receives the coefficients to optimize, in grad there is the answer of the gradients, and in data there are the id of the correlation to use and the real data used
//for the optimization
{
    FF_CorrelationData *data=data1;
    double result[data->nPoints];
    double error=0;
    int i;
    if (grad)
    {
        for (i=0;i<nVar;i++) grad[i]=0;
    }
    FF_CorrelationResult(data->eq,coef,data->nPoints,&data->x,result);
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

    algL=NLOPT_LN_NELDERMEAD;
    //algL=NLOPT_LN_SBPLX;
    //We prepare the initial values, and the optimization algorithm to use, according to the intended use
    switch (data->eq)
    {
    case FF_DIPPR100:
        for (i=0;i<nVar;i++){
            var[i]=0;
        }
        break;
    case FF_Polynomial:
        for (i=0;i<nVar;i++){
            var[i]=0;
        }
        break;
    case FF_DIPPR100Ld:
        var[0]=800;
        for (i=1;i<nVar;i++) var[i]=1;
        break;
    case FF_Polynomial2:
        break;
    case FF_DIPPR101Lv:
        var[0]=1;
        var[1]=1e3;
        var[2]=1e1;
        var[3]=-2e-6;
        var[4]=2;
        break;
    case FF_DIPPR101Vp:
        var[0]=63.9;
        var[1]=-4e3;
        var[2]=-10.34;
        var[3]=4.16e-6;
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
        var[0]=1e7;
        var[1]=0;
        var[2]=0;
        var[3]=0;
        var[4]=0;
        break;
    case FF_DIPPR106Ld:
        var[0]=500;
        var[1]=0;
        var[2]=0;
        var[3]=0;
        var[4]=0;
        break;
    case FF_DIPPR106SurfT://It seems that LsurfT requires very narrow limits. For Hv changing back +-30 to +-10 seems better. For Ldens better to compare
        var[0]=0.01;
        var[1]=0;
        var[2]=0;
        var[3]=0;
        var[4]=0;
        break;
    case FF_DIPPR107Cp:
        for (i=0;i<2;i++)
        {
            var[i]=1e4;
        }
        var[2]=1000;
        var[3]=5e3;
        var[4]=700;
        algL=NLOPT_LN_SBPLX;
        break;
    case FF_DIPPR116Ld://For liquid density
    case FF_PPDS10:
        var[1]=0;
        var[2]=0;
        var[3]=0;
        var[4]=0;
        break;
    case FF_PPDS9://Liquid Viscosity
        var[0]=3.0;
        var[1]=5.0;
        var[2]=4.0e2;
        var[3]=-1.0e2;
        var[4]=1.0e-5;
        break;
    case FF_Antoine1:
    case FF_Antoine2:
        var[0]=10.0;
        var[1]=4e3;
        var[2]=2e2;
        break;
    case FF_Wagner25://For Vp
    case FF_Wagner36:
        //coef0 is forced to Pc, and coef5 to Tc, when recovering bounds
        for (i=1;i<nVar-1;i++){
            var[i]=0;
        }
        break;
    case FF_Rackett://For liquid density
        var[0]=3e2;
        var[1]=0.2;
        var[2]=4e2;
        var[3]=0.2;
        break;
    case FF_WagnerGd://For gas saturated density
        //coef0 is forced to Rhoc, and coef5 to Tc, when recovering bounds
        for (i=1;i<nVar-3;i++){
            var[i]=-2;
        }
            var[3]=-25;
            var[4]=-50;
        break;
      case FF_PPDS15://For liquid Cp
        //coef5 is forced to Tc
        for (i=1;i<nVar-1;i++){
            var[i]=0;
        }
        break;
    default:
        //Default values for local algorithm
        for (i=0;i<nVar;i++){
            var[i]=0;
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
    nlopt_destroy(optL);
    nlopt_destroy(opt);
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
    //printf("Code:%i\n",code);
    if (code < 0)
    {
        //printf("nlopt failed!\n");

        //printf("found minimum after %d evaluations\n", count);
        //printf("found minimum at f(%g,%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],var[3],var[4],*error);
    }
    else
    {
        //printf("found minimum after %d evaluations\n", count);
        //printf("found minimum at A:%g, B:%g, C:%g, D:%g, E:%g, Error:%0.10g\n", var[0], var[1], var[2],var[3],var[4],*error);
    }
    //nlopt_destroy(opt);
}

//Recibe la EOS a usar y establece los límites para las variables
void CALLCONV FF_SaftBounds(unsigned nVar, int eos, double lb[], double ub[]){
    switch (eos){
    case FF_PCSAFT:
    case FF_PPCSAFT_GV:
    case FF_PPCSAFT_JC:
        nVar=3;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=12;
        lb[2]=150;
        ub[2]=350;
        //printf("Hola estoy en PCSAF normal. eosType:%i eos:%i\n",data->eosType,data->eos.eos);
        break;
    case FF_PCSAFT1A:
    case FF_PPCSAFT1A_GV:
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=7;
        lb[2]=150;
        ub[2]=350;
        lb[3]=0.001;
        ub[3]=0.5;
        lb[4]=1000;
        ub[4]=5000;
        break;
    case FF_PCSAFT2B:
    case FF_PPCSAFT2B_GV:
    case FF_PPCSAFT2B_JC:
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=7;
        lb[2]=150;
        ub[2]=350;
        lb[3]=0.005;
        ub[3]=0.5;
        lb[4]=1000;
        ub[4]=6000;
        break;
    case FF_PCSAFT3B:
    case FF_PPCSAFT3B_GV:
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=7;
        lb[2]=150;
        ub[2]=350;
        lb[3]=0.001;
        ub[3]=0.5;
        lb[4]=200;
        ub[4]=5500;
        break;
    case FF_PCSAFT4C:
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
        ub[4]=5000;
        break;
    case FF_SAFTVRMie:
    case FF_PSAFTVRMie_GV:
    case FF_PSAFTVRMie_JC:
        nVar=6;
        lb[0]=1.5;
        ub[0]=6;
        lb[1]=1;
        ub[1]=6;
        lb[2]=130;
        ub[2]=500;
        lb[5]=8;
        ub[5]=20;
        break;
    case FF_SAFTVRMie2B:
    case FF_SAFTVRMie2B_GV:
        nVar=6;
        lb[0]=1.5;
        ub[0]=6;
        lb[1]=1;
        ub[1]=6;
        lb[2]=130;
        ub[2]=500;
        lb[3]=0.005;
        ub[3]=0.5;
        lb[4]=1000;
        ub[4]=3000;
        lb[5]=8;
        ub[5]=20;
        break;
    }
}


//Determines the error of a SAFT EOS, using the given coefficients, and the real value supplied inside data
double CALLCONV FF_SaftFitError(unsigned nVar, const double coef[], double grad[],  const void *data1){
    FF_SAFTFitData *data=data1;
    data->eos->sigma=coef[0];
    data->eos->m=coef[1];
    data->eos->epsilon=coef[2];
    if((data->eos->nAcid>0)||(data->eos->nPos>0)){
        data->eos->kAB=coef[3];
        data->eos->epsilonAB=coef[4];
    }
    if((data->eos->eos==FF_PPCSAFT_JC)||(data->eos->eos==FF_PPCSAFT2B_JC)||(data->eos->eos==FF_PSAFTVRMie_JC)) data->eos->xp=data->xp/coef[1];
    if((data->eos->la>5)&&(data->eos->la<7)) data->eos->lr=coef[5];
    unsigned i;
    //Zc calculated value constraint
    double ZcDiff=0;
    //printf("eosType:%i eos:%i MW:%f nPos:%i nNeg:%i nAc:%i sigma:%f m:%f epsilon:%f epsilonAB:%f kAB:%f xp:%f la:%f\n",data->eosType,data->eos->eos,data->eos->MW,
    //       data->eos->nPos,data->eos->nNeg,data->eos->nAcid,data->eos->sigma,data->eos->m,data->eos->epsilon,data->eos->epsilonAB, data->eos->kAB, data->eos->xp,data->eos->la);
    double Vc,Arr,Z,aux;//Vc is the critical volume calculated from Tc,Pc and Zc. Arr and Z are the saft calculated values at Tc,Vc point
    Vc=data->eos->Zc*R*data->eos->Tc/data->eos->Pc;
    FF_ArrZfromTVSAFT(data->eos->Tc,Vc,data->eos,&Arr,&Z);
    ZcDiff=fabs((Z-data->eos->Zc)/data->eos->Zc);
    //printf("Tc:%f Vc:%f Zc:%f Zsaft:%f Zcdiff:%f\n",data->eos->Tc,Vc,data->eos->Zc,Z,ZcDiff);
    if (!(Z>0)) return +HUGE_VALF;
    if (ZcDiff>data->zcFilter) return 4.0;

    //Liquid density value constraint
    //printf("eos:%i\n",data->eos->eos);
    char option='l',state;
    double answerL[3],answerG[3];
    double densDiff=0,error;
    for (i=0;i<data->nLdPoints;i++){
        FF_VfromTPsaft(data->ldPoints[i][0],data->ldPoints[i][2],data->eos,option,answerL,answerG,&state);
        aux=fabs((data->eos->MW*1e-3/answerL[0]-data->ldPoints[i][1])/data->ldPoints[i][1]);
        densDiff=densDiff+aux;
        //printf("T:%f d:%f P:%f Vfound:%f Dfound:%f\n",data->ldPoints[i][0],data->ldPoints[i][1],data->ldPoints[i][2],answerL[0],data->eos->MW*1e-3/answerL[0]);
        //printf("densDiff:%f nLdPoints:%i\n",densDiff,data->nLdPoints);
    }

    densDiff=densDiff/data->nLdPoints;
    //printf("densDiff:%f\n",densDiff);
    if (densDiff>data->ldensFilter) return 2.0;

    //If Zc and liquid density OK we conduct the optimization
    double Vp,VpDiff=0;
    for (i=0;i<data->nVpPoints;i++)
    {
        FF_VpEos(data->eosType,data->vpPoints[i][0],data->eos,&Vp);
        aux=fabs((Vp-data->vpPoints[i][1])/data->vpPoints[i][1]);
        VpDiff=VpDiff+aux;
        //printf("Vp found:%f Vp in:%f\n",Vp,data->points[i][1]);
    }

    VpDiff=VpDiff/data->nVpPoints;
    error=0.3*VpDiff+0.7*densDiff;
    if (data->error>error){
        data->error=error;
        data->vpError=VpDiff;
        data->ldensError=densDiff;
        data->zcError=ZcDiff;
        //printf("eos:%i sigma:%f m:%f epsilon:%f Kab:%f epsilonAB:%f,lr:%f\n",data->eos->eos,data->eos->sigma,data->eos->m,data->eos->epsilon,data->eos->kAB,data->eos->epsilonAB,data->eos->lr);
        //printf("ZcDiff:%f desDiff;%f VpDiff:%f\n",ZcDiff,densDiff,VpDiff);
    }
    //printf("Error:%f\n",error/data->nPoints);
    return error;
}



//Determines the error of a cubic EOS, using the given coefficients, and the real value supplied inside data
EXP_IMP double CALLCONV FF_CubicFitError(unsigned nVar, const double coef[], double grad[],  FF_CubicFitData *data)
{
    unsigned i;
    double Vp,VpDiff=0,error=0;
    FF_CubicParam param;
    enum FF_EosType eos=FF_CubicType;
    switch (data->eos->eos){
    case FF_PR76:
    case FF_PR78:
        data->eos->Tc=coef[0];
        data->eos->Pc=coef[1];
        data->eos->w=coef[2];
        break;
    case FF_PRFIT3:
        data->eos->k1=coef[0];
        data->eos->k2=coef[1];
        data->eos->k3=coef[2];
        break;
    case FF_PRFIT4:
        data->eos->k1=coef[0];
        data->eos->k2=coef[1];
        data->eos->k3=coef[2];
        data->eos->k4=coef[3];
        break;
    default:
        data->eos->k1=coef[0];
        data->eos->k2=coef[1];
        data->eos->k3=coef[2];
        //printf("k1:%f k2:%f k3:%f\n",coef[0],coef[1],coef[2]);
        break;
    }

    //Liquid density value constraint
    char option='l',state;
    double answerL[3],answerG[3];
    double densDiff=0;
    if ((data->eos->eos==FF_PR78)||(data->eos->eos==FF_PRFIT3)||(data->eos->eos==FF_PRFIT4)){
        for (i=0;i<data->nPoints;i++){
            FF_FixedParamCubic(data->eos,&param);
            FF_ThetaDerivCubic(&data->points[i][0],data->eos,&param);
            FF_VfromTPcubic(data->points[i][0],data->points[i][1],&param,option,answerL,answerG,&state);
            densDiff=densDiff+fabs((data->eos->MW*1e-3/answerL[0]-data->points[i][2])/data->points[i][2]);
        }
        densDiff=densDiff/data->nPoints;
        if (densDiff>data->ldensFilter) return 2.0;//if density error is higher than the constraint we return a 200% error in order to reject the coefficients
    }

    for (i=0;i<data->nPoints;i++){
        FF_VpEos(eos,data->points[i][0],data->eos,&Vp);
        VpDiff=VpDiff+fabs((Vp-data->points[i][1])/data->points[i][1]);
    }
    VpDiff=VpDiff/data->nPoints;

    if ((data->eos->eos==FF_PR78)||(data->eos->eos==FF_PRFIT3)||(data->eos->eos==FF_PRFIT4)) error=0.6*VpDiff+0.4*densDiff;
    else error=VpDiff;
    if (isnan(error)) return 4.0;


    if (data->error>error){
        data->error=error;
        data->vpError=VpDiff;
        if ((data->eos->eos==FF_PR78)||(data->eos->eos==FF_PRFIT3)||(data->eos->eos==FF_PRFIT4)) data->ldensError=densDiff;

        /*
        switch (data->eos){//We print the new optimal values found
        case FF_PR76:
        case FF_PR78:
            //printf("Tc:%f Pc:%f w:%f\n",cubicData.Tc,cubicData.Pc,cubicData.w);
            break;
        case FF_PRFIT3:
            //printf("a:%f Tc:%f k1:%f\n",cubicData.k1,cubicData.k2,cubicData.k3);
            break;
        case FF_PRFIT4:
            //printf("b:%f a:%f w:%f Tc:%f\n",cubicData.k1,cubicData.k2,cubicData.k3,cubicData.k4);
            break;
        case FF_PRSV1:
        case FF_PRBM:
            //printf("k1:%f\n",cubicData.k1);
            break;
        case FF_PRMELHEM:
        case FF_PRSOF:
        case FF_SRKSOF:
            //printf("k1:%f k2:%f\n",cubicData.k1,cubicData.k2);
            break;
        default:
            //printf("k1:%f k2:%f k3:%f\n",cubicData.k1,cubicData.k2,cubicData.k3);
            break;
        }
        //printf("densDiff:%f VpDiff:%f\n",densDiff,VpDiff);
    */
    }

    return error;
}

//Determines the error of a binary viscosity mixing rule, using the given coefficients, and the actual value supplied inside data
EXP_IMP double CALLCONV FF_ViscFitError(unsigned nVar, const double coef[], double grad[],  FF_IntParamFitData *data)
{
    unsigned i;
    double x[2],visc,viscDiff=0,error=0;

    switch (data->mix->viscMixRule){
    case FF_Grunberg:
    case FF_Teja:
        data->mix->viscIntParam[0][1][0]=data->mix->viscIntParam[1][0][0]=coef[0];
        data->mix->viscIntParam[0][1][1]=data->mix->viscIntParam[1][0][1]=coef[1];
        break;
    case FF_Andrade:
        data->mix->viscIntParam[0][1][0]=data->mix->viscIntParam[1][0][0]=coef[0];
        data->mix->viscIntParam[0][1][1]=data->mix->viscIntParam[1][0][1]=coef[1];
        data->mix->viscIntParam[0][1][2]=data->mix->viscIntParam[1][0][2]=coef[2];
        data->mix->viscIntParam[0][1][3]=data->mix->viscIntParam[1][0][3]=coef[3];
        break;
    case FF_AspenMod:
        data->mix->viscIntParam[0][1][0]=data->mix->viscIntParam[1][0][0]=coef[0];
        data->mix->viscIntParam[0][1][1]=data->mix->viscIntParam[1][0][1]=coef[1];
        data->mix->viscIntParam[0][1][2]=data->mix->viscIntParam[1][0][2]=coef[2];
        data->mix->viscIntParam[0][1][3]=coef[3];
        data->mix->viscIntParam[1][0][3]=-coef[3];
        data->mix->viscIntParam[0][1][4]=coef[4];
        data->mix->viscIntParam[1][0][4]=-coef[4];
        data->mix->viscIntParam[0][1][5]=coef[5];
        data->mix->viscIntParam[1][0][5]=-coef[5];
        break;
    case FF_McAllister3:
        data->mix->viscIntParam[0][1][0]=data->mix->viscIntParam[1][0][2]=coef[0];
        data->mix->viscIntParam[0][1][1]=data->mix->viscIntParam[1][0][3]=coef[1];
        data->mix->viscIntParam[0][1][2]=data->mix->viscIntParam[1][0][0]=coef[2];
        data->mix->viscIntParam[0][1][3]=data->mix->viscIntParam[1][0][1]=coef[3];
        break;
    case FF_McAllister4:
        data->mix->viscIntParam[0][1][0]=data->mix->viscIntParam[1][0][4]=coef[0];
        data->mix->viscIntParam[0][1][1]=data->mix->viscIntParam[1][0][5]=coef[1];
        data->mix->viscIntParam[0][1][2]=data->mix->viscIntParam[1][0][2]=coef[2];
        data->mix->viscIntParam[0][1][3]=data->mix->viscIntParam[1][0][3]=coef[3];
        data->mix->viscIntParam[0][1][4]=data->mix->viscIntParam[1][0][0]=coef[4];
        data->mix->viscIntParam[0][1][5]=data->mix->viscIntParam[1][0][1]=coef[5];
        break;
    default:
        //printf("k1:%f k2:%f k3:%f\n",coef[0],coef[1],coef[2]);
        break;
    }
    for(i=0;i<data->nPoints;i++){
        x[0]=data->points[i][0];
        x[1]=1.0-x[0];
        switch (data->mix->viscMixRule){
        case FF_Grunberg:
            FF_MixLiqViscGrunberg(data->mix,data->points[i][1],data->points[i][2],x,&visc);
            break;
        case FF_Teja:
            FF_MixLiqViscTeja(data->mix,&data->points[i][1],&data->points[i][2],x,&visc);
            break;
        case FF_Andrade:
            FF_MixLiqViscAndrade(data->mix,data->points[i][1],data->points[i][2],x,&visc);
            break;
        case FF_AspenMod:
            FF_MixLiqViscAspenMod(data->mix,data->points[i][1],data->points[i][2],x,&visc);
            break;
        case FF_McAllister3:
            FF_MixLiqViscMcAllister3(data->mix,data->points[i][1],data->points[i][2],x,&visc);
            break;
        case FF_McAllister4:
            FF_MixLiqViscMcAllister4(data->mix,data->points[i][1],data->points[i][2],x,&visc);
            break;
        default:
            //printf("k1:%f k2:%f k3:%f\n",coef[0],coef[1],coef[2]);
            break;
        }
        viscDiff=viscDiff+fabs((visc-data->points[i][3])/data->points[i][3]);
    }
    error=viscDiff/data->nPoints;
    if (isnan(error)) return 4.0;
    if (data->error>error) data->error=error;
    return error;
}

//Determines the error for a mixing rule, using the given coefficients, and the actual values supplied inside data
EXP_IMP double CALLCONV FF_IntParamFitError(unsigned nVar, const double coef[], double grad[],  FF_IntParamFitData *data)
{
    unsigned i;
    double x[2],y[2],phiL[2],phiG[2],bPguess=0,bP=0,fDiff=0,error=0;
    char option;

    switch (data->mix->mixRule){
    case FF_VdW:
        data->mix->intParam[0][1][0]=data->mix->intParam[1][0][0]=coef[0];
        data->mix->intParam[0][1][1]=data->mix->intParam[1][0][1]=coef[1];
        data->mix->intParam[0][1][2]=data->mix->intParam[1][0][2]=coef[2];
        data->mix->intParam[0][1][3]=data->mix->intParam[1][0][3]=coef[3];
        break;
    case FF_PR:
        data->mix->intParam[0][1][0]=coef[0];
        data->mix->intParam[0][1][1]=coef[1];
        data->mix->intParam[0][1][2]=coef[2];
        data->mix->intParam[1][0][0]=coef[3];
        data->mix->intParam[1][0][1]=coef[4];
        data->mix->intParam[1][0][2]=coef[5];
        data->mix->intParam[0][1][3]=data->mix->intParam[1][0][3]=coef[6];
        break;
    case FF_MHV2:
        data->mix->intParam[0][1][0]=coef[0];
        data->mix->intParam[0][1][1]=coef[1];
        data->mix->intParam[0][1][2]=coef[2];
        data->mix->intParam[1][0][0]=coef[3];
        data->mix->intParam[1][0][1]=coef[4];
        data->mix->intParam[1][0][2]=coef[5];
        if(data->mix->actModel==FF_NRTL) data->mix->intParam[0][1][4]=data->mix->intParam[1][0][4]=coef[6];
        break;
    case FF_BL:
    case FF_IndAssoc:
        data->mix->intParam[0][1][0]=data->mix->intParam[1][0][0]=coef[0];
        data->mix->intParam[0][1][1]=data->mix->intParam[1][0][1]=coef[1];
        data->mix->intParam[0][1][2]=data->mix->intParam[1][0][2]=coef[2];

        break;
    default:
        //printf("k1:%f k2:%f k3:%f\n",coef[0],coef[1],coef[2]);
        break;
    }
    for(i=0;i<data->nPoints;i++){
        x[0]=data->points[i][0];
        x[1]=1.0-x[0];
        if(data->points[i][3]>0){
            option='l';
            FF_MixPhiEOS(data->mix,&data->points[i][1],&data->points[i][2],x,&option,phiL);
            y[0]=data->points[i][3];
            y[1]=1.0-y[0];
            option='g';
            FF_MixPhiEOS(data->mix,&data->points[i][1],&data->points[i][2],y,&option,phiG);
            //fDiff=fDiff+0.5*(fabs((y[0]*phiG[0]-x[0]*phiL[0])/(x[0]*phiL[0]))+fabs((y[1]*phiG[1]-x[1]*phiL[1])/(x[1]*phiL[1])));
            fDiff=fDiff+fabs(x[0]*phiL[0]/phiG[0]-y[0]);
            //printf("phiL[0]:%f phiG[0]:%f\n",phiL[0],phiG[0]);
        }
        else{
            FF_BubbleP(data->mix,&data->points[i][1],x,&bPguess,&bP,y,phiL,phiG);
            fDiff=fDiff+fabs(bP-data->points[i][2])/data->points[i][2];
        }

    }
    error=fDiff/data->nPoints;
    if (isnan(error)) return 4.0;
    if (data->error>error) data->error=error;
    return error;
}

//Optimizer function for EOS parameters. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
void CALLCONV FF_OptSAFTparam(unsigned optTime,unsigned nVar,double lb[],double ub[],char enforceLimits[], FF_SAFTFitData *data,double var[],double *error)
{
    if (data->eos->Vc==0) data->eos->Vc=data->eos->Zc*R*data->eos->Tc/data->eos->Pc;
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
    switch (data->eos->eos){
    case FF_PCSAFT:
    case FF_PPCSAFT_GV:
    case FF_PPCSAFT_JC:
        nVar=3;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=12;
        lb[2]=150;
        ub[2]=350;
        var[0]=2.5;
        var[1]=1;
        var[2]=350;
        //printf("Hola estoy en PCSAF normal. eosType:%i eos:%i\n",data->eosType,data->eos.eos);
        break;
    case FF_PCSAFT1A:
    case FF_PPCSAFT1A_GV:
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=7;
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
    case FF_PCSAFT2B:
    case FF_PPCSAFT2B_GV:
    case FF_PPCSAFT2B_JC:
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=7;
        lb[2]=150;
        ub[2]=350;
        lb[3]=0.005;
        ub[3]=0.5;
        lb[4]=1000;
        ub[4]=6000;
        var[0]=3;
        var[1]=1.5;
        var[2]=250;
        var[3]=0.02;
        var[4]=2300;
        break;
    case FF_PCSAFT3B:
    case FF_PPCSAFT3B_GV:
        nVar=5;
        lb[0]=2.0;
        ub[0]=4.1;
        lb[1]=1;
        ub[1]=7;
        lb[2]=150;
        ub[2]=350;
        lb[3]=0.001;
        ub[3]=0.5;
        lb[4]=200;
        ub[4]=5500;
        var[0]=3;
        var[1]=1.5;
        var[2]=250;
        var[3]=0.02;
        var[4]=2300;
        break;
    case FF_PCSAFT4C:
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
        ub[4]=5000;
        var[0]=2.0;
        var[1]=3.5;
        var[2]=300.0;
        var[3]=0.02;
        var[4]=2300;
        break;
    case FF_SAFTVRMie:
    case FF_PSAFTVRMie_GV:
    case FF_PSAFTVRMie_JC:
        nVar=6;
        lb[0]=1.5;
        ub[0]=6;
        lb[1]=1;
        ub[1]=6;
        lb[2]=130;
        ub[2]=500;
        var[0]=2.5;
        var[1]=1;
        var[2]=350;
        var[3]=lb[3]=ub[3]=0;
        var[4]=lb[4]=ub[4]=0;
        lb[5]=8;
        ub[5]=20;
        var[5]=12;
        break;
    case FF_SAFTVRMie2B:
    case FF_SAFTVRMie2B_GV:
        nVar=6;
        lb[0]=1.5;
        ub[0]=6;
        lb[1]=1;
        ub[1]=6;
        lb[2]=130;
        ub[2]=500;
        lb[3]=0.005;
        ub[3]=0.5;
        lb[4]=1000;
        ub[4]=3000;
        lb[5]=8;
        ub[5]=20;
        var[0]=2.5;
        var[1]=1;
        var[2]=350;
        var[3]=0.02;
        var[4]=2300;
        var[5]=12;
        break;
    }
    for (i=0;i<nVar;i++){//we recover the received values for limits and initial value,if needed
        if (enforceLimits[i]=='y'){
            lb[i]=rLb[i];
            ub[i]=rUb[i];
            var[i]=rVar[i];
        }
    }
    alg=NLOPT_GN_MLSL_LDS;
    algL=NLOPT_LN_NELDERMEAD;
    opt=nlopt_create(alg,nVar);
    optL=nlopt_create(algL,nVar);
    //nlopt_set_population(opt,10);
    nlopt_set_ftol_rel(optL, 0.001);
    nlopt_set_xtol_rel(optL, 0.0001);
    nlopt_set_local_optimizer(opt, optL);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_maxeval(opt,1e8);//number of evaluations
    nlopt_set_maxtime(opt,optTime);//Max.time in seconds
    nlopt_set_min_objective(opt, FF_SaftFitError, data);
    int code=nlopt_optimize(opt, var, error);
    nlopt_destroy(optL);
    nlopt_destroy(opt);
    printf("Return code:%i\n",code);
    if (code < 0){
        printf("nlopt failed!\n");
    }
    /*
    else{
        switch(data->eos.eos){
        case FF_SAFTVRMie:
            printf("found minimum at f(%g,%g,%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2], var[3], var[4],var[5],*error);
            break;
        case FF_PCSAFT1A:
        case FF_PCSAFT2B:
        case FF_PCSAFT4C:
            printf("found minimum at f(%g,%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],var[3],var[4],*error);
            break;
        default:
           printf("found minimum at f(%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],*error);
            break;
        }
    }
    */
    //nlopt_destroy(opt);
}


//Optimizer function for cubic EOS parameters. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
void CALLCONV FF_OptCubicParam(unsigned optTime,unsigned nVar,double lb[],double ub[],char enforceLimits[], FF_CubicFitData *data,double var[],double *error)
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
    switch (data->eos->eos){
    case FF_PR76://optimal Tc,Pc and w will be calculated
    case FF_PR78:
        nVar=3;
        lb[0]=0.9*data->eos->Tc;//will replace Tc
        ub[0]=1.1*data->eos->Tc;
        lb[1]=0.9*data->eos->Pc;//will replace Pc
        ub[1]=1.1*data->eos->Pc;
        lb[2]=0.9*data->eos->w;//will replace w
        ub[2]=1.1*data->eos->w;
        var[0]=data->eos->Tc;
        var[1]=data->eos->Pc;
        var[2]=data->eos->w;
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
        lb[0]=0.7*0.45724 * pow(R*data->eos->Tc,2)/ data->eos->Pc;;//equivalent to a
        ub[0]=1.3 * pow(R*data->eos->Tc,2)/ data->eos->Pc;
        lb[1]=0.9*data->eos->Tc;//equivalent to Tc
        ub[1]=1.1*data->eos->Tc;
        lb[2]=-1;//equivalent to k1 of FF_PRSV1
        ub[2]=1;
        var[0]=0.45724 * pow(R*data->eos->Tc,2)/ data->eos->Pc;
        var[1]=data->eos->Tc;
        var[2]=0;
        break;
    case FF_PRFIT4://Will optimize a,b,Tc and w regarding liquid density and vapor pressure
        nVar=4;
        lb[0]=0.7*0.0778 * R * data->eos->Tc / data->eos->Pc;//equivalent to b
        ub[0]=1.2*0.0778 * R * data->eos->Tc / data->eos->Pc;
        lb[1]=0.7*0.45724 * pow(R*data->eos->Tc,2)/ data->eos->Pc;//equivalent to a
        ub[1]=1.3*0.45724 * pow(R*data->eos->Tc,2)/ data->eos->Pc;
        lb[2]=0.8*data->eos->w;//equivalent to w
        ub[2]=1.2*data->eos->w;
        lb[3]=0.9*data->eos->Tc;//equivalent to Tc
        ub[3]=1.1*data->eos->Tc;
        var[0]=0.85*0.0778 * R * data->eos->Tc / data->eos->Pc;//equivalent to b
        var[1]=0.45724 * pow(R*data->eos->Tc,2)/ data->eos->Pc;//equivalent to a
        var[2]=data->eos->w;//equivalent to w
        var[3]=data->eos->Tc;//equivalent to Tc
        break;
    case FF_PRFIT4B:
        nVar=4;
        lb[0]=-1;
        ub[0]=1;
        lb[1]=0.95*data->eos->Tc;
        ub[1]=1.05*data->eos->Tc;
        lb[2]=0.95*data->eos->Pc;
        ub[2]=1.05*data->eos->Pc;
        lb[3]=0.95*data->eos->w;
        ub[3]=1.05*data->eos->w;
        var[0]=0;
        var[1]=data->eos->Tc;
        var[2]=data->eos->Pc;
        var[3]=data->eos->w;
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
    }

    for (i=0;i<nVar;i++){//we recover the received values for limits and initial value,if needed
        if (enforceLimits[i]=='y'){
            lb[i]=rLb[i];
            ub[i]=rUb[i];
            var[i]=rVar[i];
        }
    }

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
    nlopt_set_min_objective(opt, FF_CubicFitError, data);
    //nlopt_add_inequality_constraint(opt,FF_CubicEOSAlphaConstraint,data,1e-2);//this makes that the calculated alpha is not too different from PR or FF_SRK original alpha

    int code=nlopt_optimize(opt, var, error);
    nlopt_destroy(optL);
    nlopt_destroy(opt);
    printf("Return code:%i\n",code);
    if (code < 0){
        printf("nlopt failed!\n");
    }
    else{
        switch(data->eos->eos){
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

//Optimizer function for viscosity BIPs. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
void CALLCONV FF_OptViscBIPs(unsigned optTime,unsigned nVar,double lb[],double ub[],char Tindependent, FF_IntParamFitData *data,double var[],double *error)
{
    nlopt_algorithm alg,algL;
    nlopt_opt opt,optL;

    switch (data->mix->viscMixRule){
    case FF_Grunberg:
    case FF_Teja:
        nVar=2;
        lb[0]=-1.0;
        ub[0]=3.0;
        var[0]=1.0;
        if(Tindependent=='y')lb[1]=ub[1]=var[1]=0;
        else{
            lb[1]=-300.0;
            ub[1]=900.0;
            var[1]=0.0;
        }
        break;
    case FF_Andrade:
        nVar=4;
        lb[0]=-1.0;
        ub[0]=3.0;
        var[0]=1.0;
        lb[2]=-1.0;
        ub[2]=3.0;
        var[2]=1.0;
        if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[3]=ub[3]=var[3]=0;
        else{
            lb[1]=-100.0;
            ub[1]=300.0;
            var[1]=0.0;
            lb[3]=-400.0;
            ub[3]=1200.0;
            var[3]=0.0;
        }
        break;
    case FF_AspenMod:
        nVar=6;
        lb[0]=-1.0;
        ub[0]=1.0;
        var[0]=0.002;
        lb[3]=-100.0;
        ub[3]=100.0;
        var[3]=0.002;

        if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[2]=ub[2]=var[2]=lb[4]=ub[4]=var[4]=lb[5]=ub[5]=var[5]=0;
        else{
            lb[1]=-100.0;
            ub[1]=300.0;
            var[1]=0.0;
            lb[2]=-1.0;
            ub[2]=1.0;
            var[2]=0.0;
            lb[4]=-100.0;
            ub[4]=300.0;
            var[4]=0.0;
            lb[5]=-1.0;
            ub[5]=1.0;
            var[5]=0.0;
        }
        break;
    case FF_McAllister3:
        nVar=4;
        lb[0]=-1.0;
        ub[0]=1.0;
        var[0]=0.002;
        lb[2]=-1.0;
        ub[2]=1.0;
        var[2]=0.002;
        if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[3]=ub[3]=var[3]=0;
        else{
            lb[1]=-100.0;
            ub[1]=300.0;
            var[1]=0.0;
            lb[3]=-100.0;
            ub[3]=300.0;
            var[3]=0.0;
        }
        break;
    case FF_McAllister4:
        nVar=6;
        lb[0]=-1.0;
        ub[0]=1.0;
        var[0]=0.002;
        lb[2]=-1.0;
        ub[2]=1.0;
        var[2]=0.002;
        lb[4]=-1.0;
        ub[4]=1.0;
        var[4]=0.002;
        if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[3]=ub[3]=var[3]=lb[5]=ub[5]=var[5]=0;
        else{
            lb[1]=-100.0;
            ub[1]=300.0;
            var[1]=0.0;
            lb[3]=-100.0;
            ub[3]=300.0;
            var[3]=0.0;
            lb[5]=-100.0;
            ub[5]=300.0;
            var[5]=0.0;
        }
        break;
    }

    alg=NLOPT_GN_MLSL_LDS;
    algL=NLOPT_LN_NELDERMEAD;
    opt=nlopt_create(alg,nVar);
    optL=nlopt_create(algL,nVar);
    nlopt_set_ftol_rel(optL, 0.0001);
    nlopt_set_xtol_rel(optL, 0.0001);
    nlopt_set_local_optimizer(opt, optL);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_maxeval(opt,1e6);//max. number of evaluations
    nlopt_set_maxtime(opt,optTime);//Max.time in seconds
    nlopt_set_min_objective(opt, FF_ViscFitError, data);
    //nlopt_add_inequality_constraint(opt,FF_CubicEOSAlphaConstraint,data,1e-2);//this makes that the calculated alpha is not too different from PR or FF_SRK original alpha

    int code=nlopt_optimize(opt, var, error);
    nlopt_destroy(optL);
    nlopt_destroy(opt);
    printf("Return code:%i\n",code);
    if (code < 0){
        printf("nlopt failed!\n");
    }
    else{
        switch(data->mix->viscMixRule){
        case FF_AspenMod:
        case FF_McAllister4:
            printf("found minimum at f(%g,%g,%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],var[3],var[4],var[5],*error);
            break;
        case FF_Andrade:
        case FF_McAllister3:
            printf("found minimum at f(%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],var[3],*error);
            break;
        case FF_Grunberg:
        case FF_Teja:
            printf("found minimum at f(%g,%g) = %0.10g\n", var[0], var[1],*error);
            break;
        default:
           printf("found minimum at f(%g) = %0.10g\n", var[0],*error);
            break;
        }
    }

    //nlopt_destroy(opt);
}

//Optimizer function for fugacity BIPs. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
void CALLCONV FF_OptPhiBIPs(unsigned optTime,unsigned nVar,double lb[],double ub[],char Tindependent, FF_IntParamFitData *data,double var[],double *error)
{
    nlopt_algorithm alg,algL;
    nlopt_opt opt,optL;

    switch (data->mix->mixRule){
    case FF_VdW:
        nVar=4;
        lb[0]=-0.6;
        ub[0]=0.6;
        var[0]=0.0;
        lb[3]=-0.3;
        ub[3]=0.5;
        var[3]=0.0;
        if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[2]=ub[2]=var[2]=0;
        else{
            lb[1]=-0.002;
            ub[1]=0.002;
            var[1]=0.0;
            lb[2]=-200.0;
            ub[2]=200.0;
            var[2]=0.0;
        }
        break;
    case FF_PR:
        nVar=7;
        lb[0]=-0.6;
        ub[0]=0.6;
        var[0]=0.0;
        lb[3]=-0.6;
        ub[3]=0.6;
        var[3]=0.0;
        lb[6]=-0.3;
        ub[6]=0.5;
        var[6]=0.0;
        if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[2]=ub[2]=var[2]=lb[4]=ub[4]=var[4]=lb[5]=ub[5]=var[5]=0;
        else{
            lb[1]=-0.002;
            ub[1]=0.002;
            var[1]=0.0;
            lb[2]=-200.0;
            ub[2]=200.0;
            var[2]=0.0;
            lb[4]=-0.002;
            ub[4]=0.002;
            var[4]=0.0;
            lb[5]=-200.0;
            ub[5]=200.0;
            var[5]=0.0;
        }
        break;
    case FF_MHV2:
    case FF_MHV1:
    case FF_PSRK:
    case FF_LCVM:
    case FF_HV:
        if(data->mix->actModel==FF_NRTL){
            nVar=7;
            lb[0]=-5e4;
            ub[0]=6e4;
            var[0]=0.0;
            lb[3]=-5e4;
            ub[3]=5e4;
            var[3]=0.0;
            lb[6]=-10.0;//0.1;//-10
            ub[6]=100.0;//0.9//100
            var[6]=0.3;
            if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[2]=ub[2]=var[2]=lb[4]=ub[4]=var[4]=lb[5]=ub[5]=var[5]=0;
            else{
                lb[1]=-170;
                ub[1]=170;
                var[1]=0.0;
                lb[2]=-0.6;
                ub[2]=0.6;
                var[2]=0.0;
                lb[4]=-170;
                ub[4]=170;
                var[4]=0.0;
                lb[5]=-0.3;
                ub[5]=0.3;
                var[5]=0.0;
            }
        }
        else if(data->mix->actModel==FF_UNIQUAC){
            nVar=6;
            lb[0]=-1.5e3;
            ub[0]=2.5e3;
            var[0]=0.0;
            lb[3]=-1.5e3;
            ub[3]=2.5e3;
            var[3]=0.0;
            if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[2]=ub[2]=var[2]=lb[4]=ub[4]=var[4]=lb[5]=ub[5]=var[5]=0;
            else{
                lb[1]=-5;
                ub[1]=9;
                var[1]=0.0;
                lb[2]=-0.02;
                ub[2]=0.03;
                var[2]=0.0;
                lb[4]=-5;
                ub[4]=9;
                var[4]=0.0;
                lb[5]=-0.02;
                ub[5]=0.03;
                var[5]=0.0;
            }
        }
        break;
    case FF_BL:
    case FF_IndAssoc:
        nVar=3;
        lb[0]=-0.15;
        ub[0]=0.15;
        var[0]=0.0;
        if(Tindependent=='y')lb[1]=ub[1]=var[1]=lb[2]=ub[2]=var[2]=0;
        else{
            lb[1]=-0.0005;
            ub[1]=0.0005;
            var[1]=0.0;
            lb[2]=-65.0;
            ub[2]=45.0;
            var[2]=0.0;
        }
        break;
    }

    alg=NLOPT_GN_MLSL_LDS;
    algL=NLOPT_LN_NELDERMEAD;
    opt=nlopt_create(alg,nVar);
    optL=nlopt_create(algL,nVar);
    nlopt_set_ftol_rel(optL, 0.0001);
    nlopt_set_xtol_rel(optL, 0.0001);
    nlopt_set_local_optimizer(opt, optL);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_maxeval(opt,1e6);//max. number of evaluations
    nlopt_set_maxtime(opt,optTime);//Max.time in seconds
    nlopt_set_min_objective(opt, FF_IntParamFitError, data);
    //nlopt_add_inequality_constraint(opt,FF_CubicEOSAlphaConstraint,data,1e-2);//this makes that the calculated alpha is not too different from PR or FF_SRK original alpha

    int code=nlopt_optimize(opt, var, error);
    nlopt_destroy(optL);
    nlopt_destroy(opt);
    printf("Return code:%i\n",code);
    if (code < 0){
        printf("nlopt failed!\n");
    }
    else{
        switch(data->mix->mixRule){
        case FF_VdW:
            printf("found minimum at f(%g,%g,%g,%g) = %0.10g\n", var[0], var[1], var[2],var[3],*error);
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
EXP_IMP void CALLCONV FF_FractionsCalculation(unsigned nSubs, const double MW[], const double q[], const bool mass, double massFrac[], double molarFrac[]){
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

double CALLCONV FF_PresDerError(unsigned nVar, const double coef[], double grad[], const  FF_SubstanceData *data){
    if (data->model==FF_SAFTtype){
        double T=coef[0];
        double V=coef[1];
        double Vl,Vh;
        double Arr,Z;
        double Pl,P,Ph;
        double der1,der2;//first and second derivatives of pressure regarding volume, absolute value
        Vl=V*0.9999;
        Vh=V*1.0001;
        FF_ArrZfromTVSAFT(T,Vl,&data->saftData,&Arr,&Z);
        Pl=Z*R* T/ Vl;
        FF_ArrZfromTVSAFT(T,V,&data->saftData,&Arr,&Z);
        P=Z*R* T/ V;
        FF_ArrZfromTVSAFT(T,Vh,&data->saftData,&Arr,&Z);
        Ph=Z*R* T/ Vh;
        der1=fabs((Ph-P)/(Vh-V));
        der2=fabs(((der1-(P-Pl)/(V-Vl)))/(V-Vl));
        return der1+der2*1e-4;
    }
}

//Finds critical points reported by SAFT EOS, using optimizer
void CALLCONV FF_FindCriticalPoint(const FF_SubstanceData *subs, double *Tc, double *Pc, double *Vc){
    unsigned nVar=2;
    double var[2]={subs->baseProp.Tc,subs->baseProp.Vc};
    double lb[2]={subs->baseProp.Tc-5.0,subs->baseProp.Vc*0.85};
    double ub[2]={subs->baseProp.Tc+40.0,subs->baseProp.Vc*1.15};
    double error;
    double Arr,Z;
    nlopt_algorithm alg,algL;
    nlopt_opt opt,optL;
    algL=NLOPT_LN_NELDERMEAD;
    alg=NLOPT_GN_MLSL;
    opt=nlopt_create(alg,nVar);
    optL=nlopt_create(algL,nVar);
    nlopt_set_ftol_rel(optL, 0.0001);
    nlopt_set_xtol_rel(optL, 0.001);
    //nlopt_set_population(opt,16);
    nlopt_set_local_optimizer(opt, optL);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_maxeval(opt,2000000);//number of evaluations
    nlopt_set_maxtime(opt,15.0);//time in seconds
    nlopt_set_min_objective(opt,FF_PresDerError,subs);
    int code=nlopt_optimize(opt, var, &error);
    nlopt_destroy(optL);
    nlopt_destroy(opt);
    *Tc=var[0];
    *Vc=var[1];
    FF_ArrZfromTVSAFT(*Tc,*Vc,&subs->saftData,&Arr,&Z);
    *Pc=Z*R* *Tc/ *Vc;
    if (code < 0)
    {
        printf("nlopt failed!\n");

        //printf("found minimum after %d evaluations\n", count);
        printf("found minimum at f(%g,%g,%g) = %0.10g\n", var[0], var[1], error);
    }
    else
    {
        //printf("found minimum after %d evaluations\n", count);
        printf("found minimum at T:%g, V:%e, P:%f Error:%0.10g\n", *Tc, *Vc, *Pc, error);
    }
}
