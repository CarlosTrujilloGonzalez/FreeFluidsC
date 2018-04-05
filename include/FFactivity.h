/*
 * FFactivity.h
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


#ifndef FFACTIVITY_H
#define FFACTIVITY_H

#include <math.h>
#include <stdio.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFphysprop.h"


#ifdef __cplusplus
extern "C"
{
#endif
//Calculates activity coefficients according to Wilson (modified) equation, at given T and composition
EXP_IMP void CALLCONV FF_ActivityWilson(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[],
                                        const enum FF_IntParamForm *form,const double *T,const double x[],double gamma[],double *gE);

//Calculates activity coefficients according to NRTL equation, at given T and composition
EXP_IMP void CALLCONV FF_ActivityNRTL(const int *numSubs,const double pintParam[],const enum FF_IntParamForm *form,const double *T,const double x[],double gamma[],double *gE);

//Calculates activity coefficients according to UNIQUAQ equation, at given T and composition
void CALLCONV FF_ActivityUNIQUAC(const int *numSubs,const  FF_BaseProp baseProp[],const double pintParam[],
                                 const enum FF_IntParamForm *form,const double *T,const double x[],double gamma[],double *gE);

//Calculates fugacity and activity coefficients, at given T and composition, from an activity model
void CALLCONV FF_PhiFromActivity(const int *numSubs,enum FF_ActModel *model,const  FF_BaseProp baseProp[],const double pintParam[],
                                 const enum FF_IntParamForm *form,const bool *useVp,const enum FF_EosType *eosType,const void *data,
                                 const double *T,const double *P,const double x[],double gamma[],double phi[],double *gE);
#ifdef __cplusplus
}
#endif

#endif // FFACTIVITY_H
