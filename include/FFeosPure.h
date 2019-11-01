/*
 * FFeosPure.h
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

#ifndef FFEOSPURE_H
#define FFEOSPURE_H

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

#ifdef __cplusplus
extern "C"
{
#endif

//Single substance calculations
//=============================

//Get substance data from an exported file
EXP_IMP FF_SubstanceData * CALLCONV FF_SubsDataFromFile(const char *name);
//Write a substance data to a file. Adds ".sd" extension
EXP_IMP void CALLCONV FF_SubsDataToFile(const char *name,FF_SubstanceData *subsData);

//Single substance, cubic EOS calculations
//----------------------------------------
//Calculates u,w,b a cubic EOS. Using critical and others constants
EXP_IMP void CALLCONV FF_FixedParamCubic(const  FF_CubicEOSdata *data, FF_CubicParam *param);
//Calculates Theta and its derivatives, given a cubic EOS and T. Using critical and others constants
EXP_IMP void CALLCONV FF_ThetaDerivCubic(const double *T,const  FF_CubicEOSdata *data, FF_CubicParam *param);
//Arr and dArr/dV at constant T calculation for a pure substance, given T and V, according to cubic EOS
EXP_IMP void CALLCONV FF_ArrZfromTVcubic(const double *T,const double *V,const  FF_CubicParam *param,double *Arr,double *Z);
//P calculation from T and V using cubic eos
EXP_IMP void CALLCONV FF_PfromTVcubic(const double *T,const double *V,const  FF_CubicParam *param,double *P);
//V calculation for a pure substance, given T and P, according to cubic EOS
EXP_IMP void CALLCONV FF_VfromTPcubic(const double *T,const double *P,const  FF_CubicParam *param,const char *option,double resultL[3],double resultG[3],char *state);
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given T and V, according to cubic EOS
EXP_IMP void CALLCONV FF_ArrDerCubic(const double *T,const double *V,const  FF_CubicParam *param,double result[6]);
//Arr and dArr/dT as for cubic EOS
EXP_IMP void CALLCONV FF_ArrDerCubic0T(const double *T,const double *V,const  FF_CubicParam *param,double result[2]);

//Single substance FF_PCSAFT EOS calculation
//---------------------------------------
//Auxiliary calculation for ArrZfromTVPCSAFT
EXP_IMP void CALLCONV FF_calcI1I2(double m,double eta,double I[4]);
//Arr and Z calculation for a pure substance, given T and V, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_ArrZfromTVSAFT(const double *T,const double *V,const  FF_SaftEOSdata *data,double *Arr,double *Z);
//P calculation from T and V using according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_PfromTVSAFT(const double *T,const double *V,const  FF_SaftEOSdata *data,double *P);
//V,Arr and Z calculation for a pure substance, given T and P, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_VfromTPSAFT(const double *T,const double *P,const  FF_SaftEOSdata *data,
                                  const char *option,double resultL[3],double resultG[3],char *state);
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given T and V, according to FF_PCSAFT EOS
EXP_IMP void CALLCONV FF_ArrDerSAFT(const double *T,const double *V,const  FF_SaftEOSdata *data,double result[6]);
//Arr and dArr/dT as per SAFT EOS
EXP_IMP void CALLCONV FF_ArrDerSAFT0T(const double *T,const double *V,const  FF_SaftEOSdata *data,double result[2]);

//Single substance, Span and Wagner EOS calculation
//----------------------------------------
//Arr and Z calculation for a pure substance, given T and V, according to Span and Wagner EOS
EXP_IMP void CALLCONV FF_ArrZfromTVsw(const double *T,const double *V,const  FF_SWEOSdata *data,double *Arr,double *Z);
//P calculation from T and V using according to Span and Wagner EOS
EXP_IMP void CALLCONV FF_PfromTVsw(const double *T,const double *V,const  FF_SWEOSdata *data,double *P);
//P and dP_ddelta calculation for a pure substance, given reduced T and reduced rho, according to SW EOS
EXP_IMP void CALLCONV FF_PresDerSW(const double *tau,const double *delta,const  FF_SWEOSdata *data,double result[2]);
//V,Arr and Z calculation for a pure substance, given T and P, according to SW EOS
EXP_IMP void CALLCONV FF_VfromTPsw(const double *T,const double *P,const  FF_SWEOSdata *data,const char *option
                              ,double resultL[3],double resultG[3],char *state);
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given tau and delta, according to SW EOS
EXP_IMP void CALLCONV FF_ArrDerSW(const double *tau,const double *delta,const  FF_SWEOSdata *data,double result[6]);
//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given T and V, according to SW EOS
EXP_IMP void CALLCONV FF_ArrDerSWTV(const double *T,const double *V,const  FF_SWEOSdata *data,double result[6]);
//Arr and dArr/dT for a pure substance, given tau and delta, according to SW EOS
EXP_IMP void CALLCONV FF_ArrDerSW0T(const double *tau,const double *delta,const  FF_SWEOSdata *data,double result[2]);

//Single substance common calculations
//------------------------------------
//Arr and dArr/dV at constant T calculation for a pure substance, given T and V by eos
EXP_IMP void CALLCONV FF_ArrZfromTVeos(const int *eosType,const double *T,const double *V,const void *data,double *Arr,double *Z);
//P calculation from T and V by eos
EXP_IMP void CALLCONV FF_PfromTVeos(const int *eosType,const double *T,const double *V,const void *data,double *P);
//P calculation from T and V by eos
EXP_IMP void CALLCONV FF_PfromTVeosS(const double *T,const double *V,const FF_SubstanceData *data,double *P);
//V,Arr and Z calculation for a pure substance, given T and P by eos
EXP_IMP void CALLCONV FF_VfromTPeos(const int *eosType,const double *T,const double *P,const void *data,const char *option,double resultL[3],double resultG[3],char *state);
//V,Arr and Z calculation for a pure substance, given T and P by eos.
EXP_IMP void CALLCONV FF_VfromTPeosS(const double *T,const double *P,const FF_SubstanceData *data,const char *option,double resultL[3],double resultG[3],char *state);
//Boiling point calculation
EXP_IMP void CALLCONV FF_TbEOS(const int *eosType,const double *P,const void *data,double *Tb);
//Boiling point calculation
EXP_IMP void CALLCONV FF_TbEOSs(const double *P,const FF_SubstanceData *data,double *Tb);
//Vapor pressure calculation
EXP_IMP void CALLCONV FF_VpEOS(const int *eosType,const double *T,const void *data,double *Vp);
//Vapor pressure calculation
EXP_IMP void CALLCONV FF_VpEOSs(const double *T,const FF_SubstanceData *data,double *Vp);
//Thermodynamic properties calculation for a ideal gas at same T and V, from a reference state, specified by refT and refP, where H and S are 0
EXP_IMP void CALLCONV FF_IdealThermoEOS(const int *equation,const double coef[],double *refT,double *refP, FF_ThermoProperties *th0);
//Ideal gas thermodynamic properties of water calculation , from a reference state specified by the triple point where H and S are 0
EXP_IMP void CALLCONV FF_IdealThermoWater( FF_ThermoProperties *th0);
//Enthalpy and entropy calculation from T,V and P using EOS
EXP_IMP void CALLCONV FF_HSfromTVPeosS(double *T, double *V, double *P, const FF_SubstanceData *data, double *H, double *S);
//Residual extended thermodynamic properties calculation from T and V, using EOS
EXP_IMP void CALLCONV FF_ExtResidualThermoEOS(const int *eosType,const void *data, FF_ThermoProperties *thR);
//Residual extended thermodynamic properties calculation from T and V, using EOS
EXP_IMP void CALLCONV FF_ExtResidualThermoEOSs(const FF_SubstanceData *data, FF_ThermoProperties *thR);
//Thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0
EXP_IMP void CALLCONV FF_ThermoEOS(const int *eosType,const void *data,const int *equation,const double coef[],double *refT,double *refP, FF_ThermoProperties *th);
//Thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0
EXP_IMP void CALLCONV FF_ThermoEOSs(const FF_SubstanceData *data, FF_ThermoProperties *th);
//Calculation of thermo properties and liquid/gas fractions from P and H(or U or S)
EXP_IMP void CALLCONV FF_ThermoEOSfromPX(const int *eosType,const void *data,const int *equation,const double coef[],double *refT,double *refP,char *var,
                                     FF_ThermoProperties *th,double *liqFraction);
//Calculation of thermo properties and liquid/gas fractions from V and T(or U or S)
EXP_IMP void CALLCONV FF_ThermoEOSfromVX(const int *eosType,const void *data,const int *equation,const double coef[],double *refT,double *refP,char *var,
                                     FF_ThermoProperties *th,double *liqFraction);

//Calculation of T,V and liquid/gas fractions from S and P
EXP_IMP void CALLCONV TVfromVHeos(const int *eosType,const void *data,const  FF_Correlation *cp,double *refT,double *refP, FF_ThermoProperties *th,double *liqFraction);

#ifdef __cplusplus
}
#endif
#endif /* FFEOSPURE_H */
