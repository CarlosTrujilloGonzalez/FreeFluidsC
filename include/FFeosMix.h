/*
 * FFeosMix.h
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

// contains EOS definitions for pure substances
#ifndef FFEOSMIX_H
#define FFEOSMIX_H

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

#include <stdbool.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFactivity.h"


#ifdef __cplusplus
extern "C"
{
#endif

//Mixture calculations
//====================

//Mixture cubic EOS calculations
//-------------------------------
//Calculates Theta,b,delta and epsilon for a mixture, given a cubic EOS,a mixing rule, composition, and pure substance parameters
EXP_IMP void CALLCONV FF_MixParamCubicEOS(const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
                                     const double pintParam[],const double x[], FF_CubicParam *param,double dTheta_dXi[],double db_dXi[],double dc_dXi[]);

//Calculates Theta,b,dTheta/dT, d2Theta/dT2, dTheta/dX[i] and db/dX[i] for a mixture, given a cubic EOS,a mixing rule, composition, and pure substance parameters
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_MixParamCubicEOSOld(const enum FF_EOS eos[],const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicEOSdata data[],
        const double pintParam[],const double x[], FF_CubicParam *param,double dTheta_dXi[],double db_dXi[],double dc_dXi[]);

//Calculates Theta,b,delta and epsilon for a mixture, given a cubic EOS,a gE mixing rule, composition,pure substance parameters, and gE
//EXP_IMP void CALLCONV calcMixParamCubicEOSgE(const enum FF_EOS *eos,const enum FF_MixingRule *rule,const double *T,const int *numSubs,const  FF_CubicParam sParam[],
//        const double *gE,const double x[], FF_CubicParam *param,double dNb_dNi[]);


//Mixture FF_PCSAFT EOS calculation
//------------------------------
//Mixture Z and Arr calculation for a mixture, given T and V, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_MixArrZfromTVPCSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                                        const  FF_SaftEOSdata data[],const double pintParam[],const double x[],double *Arr,double *Z);
//Mixture P calculation given T, V, and composition according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_MixPfromTVPCSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                                        const  FF_SaftEOSdata data[],const double pintParam[],const double x[],double *P);
//Mixture V,Arr and Z  calculation, given T and P and composition, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_MixVfromTPPCSAFT(const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const  FF_SaftEOSdata data[],
                                     const double pintParam[],const double x[],const char *option,double resultL[3],double resultG[3],char *state);
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a mixture, given T and V, according to FF_PCSAFT EOS
//--------------------------------------------------------------------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_MixArrDerPCSAFT(const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,
                              const  FF_SaftEOSdata data[],const double pintParam[],const double x[],double result[6]);

//Mixture common calculations
//---------------------------
//Mixture P calculation from T and V by eos
EXP_IMP void CALLCONV FF_MixPfromTVeos(const enum FF_EOS eos[],const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,const void *data,
                                  const double pintParam[],const double x[],double *P);
//Mixture V calculation from T and P by eos
EXP_IMP void CALLCONV FF_MixVfromTPeos(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const void *data,
                                  const double pintParam[],const double x[],const char *option,double resultL[3],double resultG[3],char *state);
//Mixture Ideal gas thermodynamic properties calculation, from a reference state, specified by T and P, where H and S are 0
EXP_IMP void CALLCONV FF_MixFF_IdealThermoEOS(const int *numSubs,const  FF_Correlation cp0[],const double x[],double *refT,double *refP, FF_ThermoProperties *th0);

//Mixture Residual thermodynamic properties calculation from T and V, using EOS
EXP_IMP void CALLCONV FF_MixResidualThermoEOS(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const int *numSubs,const void *data,const double pintParam[],
                                         const double x[], FF_ThermoProperties *thR);
//Mixture thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0
EXP_IMP void CALLCONV FF_MixThermoEOS(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const int *numSubs,const void *data,const double pintParam[],
                                 const double x[],const  FF_Correlation cp0[],double *refT,double *refP, FF_ThermoProperties *th);
//Mixture fugacity coeff. calculation from T and P by eos
EXP_IMP void CALLCONV FF_MixFugacityEOS(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const double *V,const int *numSubs,const void *data,
                                        const double pintParam[],const double x[],double phi[]);
//Mixture bubble temperature calculation, given P, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_BubbleT(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *P,const int *numSubs,const void *data,
                     const double pintParam[],const double x[],const double *bTguess,double *bT,double y[],double substPhiL[],double substPhiG[]);
//Mixture dew temperature calculation, given P, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_DewT(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *P,const int *numSubs,const void *data,
                           const double pintParam[],const double y[],const double *dTguess,double *dT,double x[],double substPhiL[],double substPhiG[]);
//Mixture bubble pressure calculation, given T, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_BubbleP(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const int *numSubs,const void *data,
                              const double pintParam[],const double x[],const double *bPguess, double *bP,double y[],double substPhiL[],double substPhiG[]);
//Mixture dew pressure calculation, given T, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_DewP(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const int *numSubs,const void *data,
                           const double pintParam[],const double y[],const double *dPguess,double *dP,double x[],double substPhiL[],double substPhiG[]);
//Pressure envelope of a binary mixture
EXP_IMP void CALLCONV FF_PressureEnvelope(const enum FF_EosType *eosType,const enum FF_MixingRule *rule,const double *T,const void *data,
                                  const double pintParam[],const int *numPoints,double c[],double x[],double y[],double bP[],double dP[]);
//VL flash calculation, given T, P, feed composition, eos and mixing rule
EXP_IMP void CALLCONV FF_VLflashTP(const enum FF_EOS *eos,const enum FF_MixingRule *rule,const double *T,const double *P,const int *numSubs,const void *data,
                     const double pintParam[],const double f[],double x[],double y[],double *beta);

#ifdef __cplusplus
}
#endif

#endif // FFEOSMIX_H

