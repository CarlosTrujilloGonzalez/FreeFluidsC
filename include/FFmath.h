/*
 * FFmath.h
 *
 *  Created on: 23/02/2017
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

// contains auxiliary mathematical functions
//==========================================

#ifndef FFMATH_H
#define FFMATH_H

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

#ifdef __cplusplus
extern "C"
{
#endif

void FFnormalize(int n, double x[]);
//Bisection solver for one variant equations. The supplied interval must comprise a root of the function
void FFsolverBisection(double (*f)(double), double *x, double *y, int *n, double xmin, double xmax, double *ftol, int *niter);//xmin and xmax define the interval
//Regula falsi solver for one variant equations. The supplied interval must comprise a root of the function
void FFsolverRegula(double (*f)(double), double *x, double *y, int *n, double xmin, double xmax, double *ftol, int *niter);//xmin and xmax define the interval
//Regula falsi, Illinois modified, solver for one variant equations. The supplied interval must comprise a root of the function
void FFsolverRegulaI(double (*f)(double), double *x, double *y, int *n, double xmin, double xmax, double *ftol, int *niter);//xmin and xmax define the interval
//Regula falsi, Anderson-Bjork modified, solver for one variant equations. The supplied interval must comprise a root of the function
void FFsolverRegulaA(double (*f)(double), double *x, double *y, int *n, double xmin, double xmax, double *ftol, int *niter);//xmin and xmax define the interval
void FFsolverNewtonND(double (*f)(double), double *x, double *y, int *n, double xini, double *ftol, int *niter);//xini is the initial guess
void FFsolverNewton(double (*f)(double),double (*d)(double), double *x, double *y, int *n, double xini, double *ftol, int *niter);//xini is the initial guess
void FFsolverIntervalNewtonND(double (*f)(double), double x[], double y[], int *n, double xmin, double xmax, double *ftol, int *niter);//xmin and xmax define the interval
void FFsolverIntervalNewtonND2(double (*f)(double), double x[], double y[], int *n, double xmin, double xmax, double *ftol, int *niter);//xmin and xmax define the interval

#ifdef __cplusplus
}
#endif
#endif /* FFMATH_H */
