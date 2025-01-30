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
#include "FFeosMix.h"
#include "FFequilibrium.h"

#ifdef __cplusplus
extern "C"
{
#endif

//Solvers
//-------
//General Regula Falsi, Anderson-Bjork modified, solver for one variant equations. The supplied interval must comprise a root of the function
EXP_IMP double FF_solverRegulaBase(double y, double (*f)(double, void *), void *data, double a, double b, double ftol, int niter);//a and b define the interval
//Diferential evolution optimizer
EXP_IMP void CALLCONV FF_OptimizerDE(unsigned nVar,double lb[],double ub[],double (*f)(unsigned,const double *,double *,const void *),void *data,double var[],double *error);

//Error calculation functions
//---------------------------
//Recibe la ecuación a usar y establece los límites para las variables
EXP_IMP void CALLCONV FF_CorrelationBounds(unsigned nVar, int eq, double lb[], double ub[]);
//Determines the error of a correlation, using the given coefficients, and the real value supplied inside functionData
EXP_IMP double CALLCONV FF_CorrelationError(unsigned nVar, const double coef[], double grad[], const  void *data1);
//Optimizer function for correlations. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
EXP_IMP void CALLCONV FF_OptCorrelation(unsigned nVar,double lb[],double ub[],char enforceLimits[],const  FF_CorrelationData *data,double var[],double *error);
//Determines the error of a SAFT EOS, using the given coefficients, and the real value supplied inside data
EXP_IMP double CALLCONV FF_SaftFitError(unsigned nVar, const double coef[], double grad[],  const void *data1);
//Determines the error of a cubic EOS, using the given coefficients, and the real value supplied inside data
EXP_IMP double CALLCONV FF_CubicFitError(unsigned nVar, const double coef[], double grad[],  FF_CubicFitData *data);
//Determines the error of a binary viscosity mixing rule, using the given coefficients, and the actual value supplied inside data
EXP_IMP double CALLCONV FF_ViscFitError(unsigned nVar, const double coef[], double grad[],  FF_IntParamFitData *data);
//Determines the error for a mixing rule, using the given coefficients, and the actual values supplied inside data
EXP_IMP double CALLCONV FF_IntParamFitError(unsigned nVar, const double coef[], double grad[],  FF_IntParamFitData *data);


//Determines the error of a  cubic EOS, using the given coefficients, and the real value supplied inside data
EXP_IMP double CALLCONV FF_CubicPvRhoError(unsigned nVar, const double coef[], double grad[], FF_EOSPvRhoData *data);
//Recibe la EOS a usar y establece los límites para las variables
EXP_IMP void CALLCONV FF_SaftBounds(unsigned nVar, int eos, double lb[], double ub[]);
//Determines the error of a SAFT EOS, using the given coefficients, and the real value supplied inside data
EXP_IMP double CALLCONV FF_SaftPvRhoError(unsigned nVar, const double coef[], double grad[],  FF_EOSPvRhoData *data);
//Optimizer function for SAFT parameters. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
EXP_IMP void CALLCONV FF_OptSAFTparam(unsigned optTime,unsigned nVar,double lb[],double ub[],char enforceLimits[], FF_SAFTFitData *data,double var[],double *error);
//Optimizer function for cubic EOS parameters. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
EXP_IMP void CALLCONV FF_OptCubicParam(unsigned optTime,unsigned nVar,double lb[],double ub[],char enforceLimits[], FF_CubicFitData *data,double var[],double *error);
//Optimizer function for viscosity BIPs. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
EXP_IMP void CALLCONV FF_OptViscBIPs(unsigned optTime,unsigned nVar,double lb[],double ub[],char Tindependent, FF_IntParamFitData *data,double var[],double *error);
//Optimizer function for fugacity BIPs. Recibes: number of variables to optimize, bounds, data to pass to the error calculation function. Returns optimized values and error
EXP_IMP void CALLCONV FF_OptPhiBIPs(unsigned optTime,unsigned nVar,double lb[],double ub[],char Tindependent, FF_IntParamFitData *data,double var[],double *error);
EXP_IMP double CALLCONV FF_PresDerError(unsigned nVar, const double coef[], double grad[], const  FF_SubstanceData *data);
EXP_IMP void CALLCONV FF_FindCriticalPoint(const FF_SubstanceData *subs, double *Tc, double *Pc, double *Vc);
//Other auxiliary functions
//-------------------------
//Determines the mass and molar fractions, given quatities in mass or moles
EXP_IMP void CALLCONV FF_FractionsCalculation(unsigned nSubs, const double MW[], const double q[], const bool mass, double massFrac[], double molarFrac[]);

#ifdef __cplusplus
}
#endif
#endif // FFTOOLS_H
