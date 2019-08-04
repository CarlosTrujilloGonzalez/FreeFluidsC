/*
 * FFphysprop.h
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

// contains mainly definitions for pure substances physical properties calculations, by non-eos methods
#ifndef FFPHYSPROP
#define FFPHYSPROP
#if (defined(_WIN32) || defined(__WIN32__))
  #define CALLCONV __cdecl//Definition of the calling convention, option is __stdcall
  /*** You should define EXPORTS at the gcc command line only when building the shared library ***/
  #ifdef EXPORTS
    #define EXP_IMP __declspec(dllexport)/*used to declare functions as export*/
  #else
    #ifdef IMPORTS
      #defineEXP_IMP __declspec(dllimport)/*used to declare functions as import*/
    #else
      #define EXP_IMP
    #endif
  #endif
#else/*If we are not in Windows, we give EXP_IMP and CALLCONV null value*/
  #define EXP_IMP
  #define CALLCONV
#endif

#include "FFbasic.h"

#ifdef __cplusplus
extern "C"
{
#endif
//Single substance, physical properties correlations
//--------------------------------------------------
//Calculates the result of the given equation
EXP_IMP void CALLCONV FF_CorrelationResult(const int *eq,const double coef[],const int *nPoints,double x[],double y[]);

//Calculates the result of the given physical property in SI units
EXP_IMP void CALLCONV FF_PhysPropCorr(const int *cor,const double coef[],const double *MW,const int *nPoints,double x[],double y[]);

//Calculates specific enthalpy and entropy from a Cp correlation, with reference T=0 K
void FF_SpecificEnthalpyEntropyCorr(int *cor,const double coef[],double *MW,const int *nPoints,double x[],double H[],double S[]);

//Vapor pressure related calculations
//-----------------------------------
//Acentric factor calculation by the definition equation
EXP_IMP void CALLCONV FF_WfromDefinition(const double *Pc,const double *Vp,double *w);

//Acentric factor calculation from one vapor pressure
EXP_IMP void CALLCONV FF_WfromOneVap(const double *Tc,const double *Pc,const double *T,const double *Vp,double *w);

//Vapor pressure using Ambrose-Walton equation. Needs only Tc,Pc and w
EXP_IMP void CALLCONV FF_VpAmbroseWalton(const FF_BaseProp *baseProp,const double *T,double *Vp);

//Vapor pressure using Riedel-Vetere equation. Needs only Tc,Pc and one boiling point
EXP_IMP void CALLCONV FF_VpRiedelVetere(const  FF_BaseProp *baseProp,const double *Tref,const double *VpRef,const double *T,double *Vp);

//Vapor pressure calculation by correlations
EXP_IMP void CALLCONV FF_Vp(double *T,FF_SubstanceData *data,double *Vp);

//Liquid density calculations
//---------------------------
//Rackett equation for satured liquid density. If you supply ref rho, Vc is not used. It is better to use w than Zra, and Zra than Zc
EXP_IMP void CALLCONV FF_LiqDensSatRackett(const  FF_BaseProp *baseProp,const double *Tref,const double *rhoRef,const double *T,double *rho);

//Chueh-Prausnitz pressure correction for liquid density
EXP_IMP void CALLCONV FF_LiqDensChuehPrausnitz(const  FF_BaseProp *baseProp, const double *T,const double *P,const double *Pin,double *rhoIn, double *rho);

//Tait equation for polymer density, with T and P dependence
EXP_IMP void CALLCONV FF_LiqDensTait(const int *eq,const double coef[],const  FF_BaseProp *baseProp, const int *nPoints,const double T[],const double P[], double rho[]);

//Liquid density from T,P
EXP_IMP void CALLCONV FF_LiqDensTP(double *T,double *P,FF_SubstanceData *data,double *liqDens);

//Others thermodynamic properties
//-------------------------------
//Liquid Cp. Bondi method
void CALLCONV FF_LiqCpBondi(const  FF_SubstanceData *data,const double *T,double *Cp);

//Transport properties of liquids
//-------------------------------
//Lucas liquid viscosity pressure correction.
EXP_IMP void CALLCONV FF_LiqViscPcorLucas(double *T,double *P,double *Pref,FF_BaseProp *data,double *rpVisc,double *visc);

//Liquid viscosity from T,P
EXP_IMP void CALLCONV FF_LiqViscTP(double *T,double *P,FF_SubstanceData *data,double *liqVisc);

//Teja-Rice method for mixture liquid viscosity
//BIP is missing
EXP_IMP void CALLCONV FF_MixLiqViscTeja(FF_MixData *mix,double *T,double *P,double x[],double *visc);

//Thermal conductivity of liquids. Latini method
EXP_IMP void CALLCONV FF_LiquidThCondLatini(double *T,FF_BaseProp *data,double *thCond);

//Liquid thermal conductivity from T
EXP_IMP void CALLCONV FF_LiqThCondT(double *T,FF_SubstanceData *data,double *liqThCond);

//Li method for mixture liquid thermal conductivity
EXP_IMP void CALLCONV FF_MixLiqThCondLi(FF_MixData *mix,double *T,double *P,double x[],double *thCond);

//SurfaceTension, MacLeod-Sugden method. Very sensible to Parachor value
EXP_IMP void CALLCONV FF_SurfTensMcLeod(double *T,FF_SubstanceData *data,double *surfTens);

//Surface tension,Sastri-Rao method
EXP_IMP void CALLCONV FF_SurfTensSastri(double *T,FF_BaseProp *data,double *surfTens);

//Surface tension from T
EXP_IMP void CALLCONV FF_SurfTensT(double *T,FF_SubstanceData *data,double *surfTens);

//Linear method for mixture surface tension
EXP_IMP void CALLCONV FF_MixLiqSurfTensLinear(FF_MixData *mix,double *T,double *P,double x[],double *surfTens);

//Winterfeld method for mixture surface tension.
void CALLCONV FF_MixLiqSurfTensWinterfeld(FF_MixData *mix,double *T,double *P,double x[],double *surfTens);

//McLeod-Sugden method for mixture surface tension.
void CALLCONV FF_MixLiqSurfTensMcLeod(FF_MixData *mix,double *rhoL,double *rhoG,double x[],double y[],double *surfTens);

//Transport properties of gases
//-----------------------------
//Gas viscosity pressure prediction/correction. Chung method
EXP_IMP void CALLCONV FF_GasViscTVcpChung(double *T,double *V,FF_BaseProp *data,double *ldVisc,double *visc);

//Gas viscosity pressure prediction/correction. Lucas method
EXP_IMP void CALLCONV FF_GasViscTPcpLucas(const double *T,const double *P,const FF_BaseProp *data,double *lpVisc,double *visc);

//Gas viscosity from T,P
void CALLCONV FF_GasViscTP(double *T,double *P,FF_SubstanceData *data,double *gasVisc);

//Viscosity of gas mixtures. Wilke method
EXP_IMP void CALLCONV FF_MixGasViscTPcpWilke(FF_MixData *mix,double *T,double *P,double y[],double *gVisc);

//Viscosity of gas mixtures. Lucas method
EXP_IMP void CALLCONV FF_MixGasViscTPcpLucas(FF_MixData *mix,double *T,double *P,double y[],double *gVisc);

//Gas low pressure thermal conductivity prediction
EXP_IMP void CALLCONV FF_GasLpThCondTCpChung(double *T,double *Cp0,FF_BaseProp *data,double *ldThCond);

//Gas thermal conductivity pressure correction
EXP_IMP void CALLCONV FF_GasThCondTVcorChung(double *T,double *V,FF_BaseProp *data,double *ldThCond,double *thCond);

//Gas thermal conductivity from T,V
EXP_IMP void CALLCONV FF_GasThCondTV(double *T,double *V,FF_SubstanceData *data,double *gasThCond);

//Thermal conductivity of low pressure gas mixtures. Mason and Saxena method
EXP_IMP void CALLCONV FF_MixLpGasThCondTpMason(FF_MixData *mix,double *T,double y[],double *gThCond);
#ifdef __cplusplus
}
#endif
#endif // FFPHYSPROP

