/*
 * FFtools.h
 *
 *  Created on: 10/12/2015
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

// contains mainly definitions for pure substances parameters fitting
#ifndef FFTOOLS_H
#define FFTOOLS_H

#if (defined(_WIN32) || defined(__WIN32__))
  #define CALLCONV __cdecl//Definition of the calling convention, option is __stdcall
  /*** You should define EXPORTS at the gcc command line only when building the shared library ***/
  #ifdef EXPORTS
    #define EXP_IMP __declspec(dllexport)/*used to declare functions as export*/
  #else
    #ifdef IMPORTS
      #define EXP_IMP __declspec(dllimport)/*used to declare functions as import*/
    #else
      #define EXP_IMP
    #endif
  #endif
#else/*If we are not in Windows, we give EXP_IMP and CALLCONV null value*/
  #define EXP_IMP
  #define CALLCONV
#endif

#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFphysprop.h"

enum FF_OptModel{corr,cubicParam,SAFTParam};
typedef struct {enum FF_CorrEquation eq ;double Tc,Pc,rhoC;unsigned nPoints;double x[40],y[40];}FF_CorrelationData;//The dimension of arrays must be declared in structures
typedef struct {enum FF_EOS eos ;double MW,Tc,Pc,Zc,w,VdWV,mu,xp,m,chi,ldensFilter,zcFilter,error,vpError,ldensError,zcError;unsigned nPoints;double points[40][3];}FF_EOSPvRhoData;//The dimension of arrays must be declared in structures
typedef struct{int eosType;FF_SaftEOSdata *eos;double xp,ldensFilter,zcFilter,error,vpError,ldensError,zcError;unsigned nPoints,nVpPoints,nLdPoints;double vpPoints[30][2],ldPoints[30][3];}FF_SAFTFitData;
typedef struct{int eosType;FF_CubicEOSdata *eos;double ldensFilter,zcFilter,error,vpError,ldensError,zcError;unsigned nPoints;double points[40][3];}FF_CubicFitData;
#ifdef __cplusplus
extern "C"
{
#endif



//Error calculation functions
//---------------------------
//Determines the error of a correlation, using the given coefficients, and the real value supplied inside functionData
EXP_IMP double CALLCONV FF_CorrelationError(unsigned nVar, const double coef[], double grad[], const  FF_CorrelationData *data);
//Optimizer function for correlations. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
EXP_IMP void CALLCONV FF_OptCorrelation(unsigned nVar,double lb[],double ub[],char enforceLimits[],const  FF_CorrelationData *data,double var[],double *error);

//Determines the error of a  cubic EOS, using the given coefficients, and the real value supplied inside data
EXP_IMP double CALLCONV FF_CubicPvRhoError(unsigned nVar, const double coef[], double grad[], FF_EOSPvRhoData *data);
//Determines the error of a SAFT EOS, using the given coefficients, and the real value supplied inside data
double CALLCONV FF_SaftPvRhoError(unsigned nVar, const double coef[], double grad[],  FF_EOSPvRhoData *data);
//Optimizer function for SAFT parameters. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
void CALLCONV FF_OptSAFTparam(unsigned optTime,unsigned nVar,double lb[],double ub[],char enforceLimits[], FF_SAFTFitData *data,double var[],double *error);
//Optimizer function for cubic EOS parameters. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
EXP_IMP void CALLCONV FF_OptCubicParam(unsigned optTime,unsigned nVar,double lb[],double ub[],char enforceLimits[], FF_CubicFitData *data,double var[],double *error);

//Other auxiliary functions
//-------------------------
//Determines the mass and molar fractions, given quatities in mass or moles
EXP_IMP double CALLCONV FF_FractionsCalculation(unsigned nSubs, const double MW[], const double q[], const bool mass, double massFrac[], double molarFrac[]);

#ifdef __cplusplus
}
#endif
#endif // FFTOOLS_H
