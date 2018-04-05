/*
 * FFbasic.h
 *
 *  Created on: 11/03/2018
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

// contains basic definitions for FreeFluids

#ifndef FFBASIC_H
#define FFBASIC_H

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

//Global variables
//----------------
extern const double Av; //molecules/mol
extern const double kb;//J/K
extern const double Pi;
extern const double R; //Pa·m3/(K·mol)
extern const double FF_PCSAFTap[7][3];
extern const double FF_PCSAFTbp[7][3];

enum FF_EOS{FF_IdealGas,FF_PR76,FF_PR78,FF_PRSV1,FF_PRBM,FF_PRMELHEM,FF_PRSOF,FF_PRALMEIDA,FF_PRMC,FF_PRTWU91,FF_PRvTWU91,FF_PRTWU95,FF_PRPOL1,FF_PRFIT3,FF_PRFIT4,FF_PRFIT4B,FF_SRK,
         FF_SRKSOF,FF_SRKMC,FF_SRKTWU91,FF_SRKPOL2,FF_PCSAFT,FF_PCSAFT1,FF_PCSAFT2B,FF_PCSAFT4C,FF_DPCSAFT_GV,FF_DPCSAFT_JC,FF_PCSAFTPOL1,FF_SW,FF_IAPWS95,IF97};//OK. Used
enum FF_EosType{FF_NoType,FF_IdealType,FF_CubicType,FF_CubicPR,FF_CubicSRK,FF_SAFTtype,FF_SWtype};
enum FF_MixingRule{FF_NoMixRul,FF_VdW,FF_PR,FF_MKP,MHV1,PSRK,HVOS,LCVM,MHV2,OPTgE,PSRKnew,VTPR};
enum FF_IntParamForm{FF_NoForm,FF_Pol1,FF_Pol1K,FF_Pol1J,FF_Pol1C,FF_Pol2,FF_Pol2K,FF_Pol2J,FF_Pol2C,FF_Pol3,FF_Pol3K,FF_Pol3J,FF_Pol3C,FF_Ant1,FF_Ant2,FF_Ant3};
//FF_Pol1 corresponds to:a+b*T
//FF_Pol2: a+b*T+c*T^2
//FF_Pol3: a+b/T+c*T
//FF_Ant1: a+b/T+c*ln(T)+d*T
//FF_Ant2: a+b/T+c*ln(T)+d*T+e*T^2
//FF_Ant3: a+b/T+c*ln(T)+d*T+e/T^2.
//In case the result contains energy, K,J,C indicates if the energy is expresed in K, Joules o Calories

typedef struct {int id;double MW,Tc,Pc,Vc,Zc,w,Zra,r,q,qRes,mu,Tb;}FF_BaseProp;
typedef struct {double x,y;}FF_SinglePointData;
typedef struct {double a,Theta,b,c,u,w,dTheta,d2Theta;}FF_CubicParam;
typedef struct {int id;enum FF_EOS eos;double MW,Tc,Pc,Zc,w,VdWV,c,k1,k2,k3,k4;}FF_CubicEOSdata;
typedef struct {int id;enum FF_EOS eos;double MW,Tc,Pc,Zc,w,sigma,m,epsilon,kAB,epsilonAB,mu,xp;int nPos,nNeg,nAcid;}FF_SaftEOSdata;
typedef struct {int id;enum FF_EOS eos;double MW,Tc,Pc,Zc,w,tRef,rhoRef,n[24],t[24];int d[24],c[24],nPol,nExp;}SWSEOSdata;
typedef struct {int id;enum FF_EOS eos;double MW,Tc,Pc,Zc,w,tRef,rhoRef,n[60],t[60],a[8],e[8],b[8],g[8],af[5],bf[5],
                Af[5],Bf[5],Cf[5],Df[5],betaf[5];int d[60],c[55],nPol,nExp,nSpec,nFinal;}FF_SWEOSdata;
typedef struct {double MW,T,P,V,A,G,S,U,H,dP_dT,dP_dV,Cv,Cp,SS,JT,IT;}FF_ThermoProperties;
typedef struct {double MW,T,P,V,Vt,rho,Z,Arr,dP_dV,dP_dT,dArr_dT,d2Arr_dT2,Ar,Gd,Phi,Ur,Hr,Sr,Cvr,Cpr,U,H,S,Cv,Cp,SS,JT;}ThermoProp;//deprecated
typedef struct {double uE,hE,gE;}ExcessProp;
//typedef struct FF_Correlation{int form;double A,B,C,D,E,F,G,H,I,J,K,L,M,limI,limS;};//deprecated
typedef struct {int id,form;double coef[13],limI,limS;}FF_Correlation;

enum FF_CorrEquation{FF_DIPPR100,FF_Polynomial,FF_Polynomial2,FF_DIPPR100Ld,FF_expDIPPR100,FF_DIPPR101,FF_DIPPR101Vp,FF_DIPPR101Lv,FF_logDIPPR101,
                  FF_DIPPR102,FF_DIPPR103,FF_DIPPR104,FF_DIPPR105,FF_DIPPR106,FF_DIPPR106Hv,FF_DIPPR106Ld,FF_DIPPR106SurfT,FF_DIPPR107,
                  FF_DIPPR107Cp,FF_DIPPR114,FF_DIPPR115,FF_DIPPR116,FF_DIPPR116Ld,FF_Wilhoit,FF_Cooper,FF_Jaechske,FF_ChemSep16,FF_Antoine1,
                  FF_Antoine2,FF_Wagner25,FF_Wagner36,FF_PPDS10,FF_PCWIN,FF_Rackett,FF_ExtAndrade1,FF_ExtAndrade2,FF_ChericVisc,FF_WagnerGd};

enum FF_ActModel{FF_NoModel,Wilson,NRTL,UNIQUAC};

typedef struct {int id,model;double refT,refP;FF_BaseProp baseProp;FF_SinglePointData cp0,vp,hVsat,lCp,lDens,lVisc,lThC,lSurfT,gVisc,gThC,sDens,sCp;FF_CubicEOSdata cubicData;FF_SaftEOSdata saftData;FF_SWEOSdata swData;FF_Correlation cp0Corr,vpCorr,btCorr,hVsatCorr,lCpCorr,lDensCorr,lViscCorr,
                lThCCorr,lSurfTCorr,gDensCorr,gViscCorr,gThCCorr,sDensCorr,sCpCorr;}FF_SubstanceData;

#endif /* FFBASIC_H */
