/*

 * FFmath.c
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

#include <math.h>
#include <stdio.h>
#include "FFmath.h"
void FFnormalize(int n, double x[])
{   double sum=0;
    int i;
    for (i=0;i<n;i++) sum+=x[i];
    for ( i = 0; i < n; i++) x[i] /= sum;
}
void FFsolverBisection(double (*f)(double), double *x, double *y, int *n, double xmin, double xmax, double *ftol, int *niter)//xmin and xmax define the interval
{
    *n=0;
    double fxmin,fxmax;
    fxmin=(*f)(xmin);
    fxmax=(*f)(xmax);
    if ((fxmin* fxmax)>0){//if both functions have the same sign, it is not sure that a root exists
        double interval=xmax-xmin;
        double xmed,fxmed;
        xmed=xmin+interval*0.5;//we try in the middle point
        fxmed=(*f)(xmed);
        if (fxmin*fxmed<0)xmin=xmed;
        else{
            *x=*y=HUGE_VAL;
            return;
        }
    }
    else if ((fxmin* fxmax)==0){//if one of the given value was already the solution
        if (fxmin==0){
            *x=xmin;
            *y=fxmin;
        }
        else{
            *x=xmax;
            *y=fxmax;
        }
        return;
    }
    while (*n<*niter){
        *x=(xmin+xmax)*0.5;
        *y=(*f)(*x);
        if ((*y<*ftol)&&(*y>(-*ftol))) break;
        if ((fxmax* *y)>0){
            xmax=*x;
            fxmax=*y;
        }
        else{
            xmin=*x;
        }
        *n=*n+1;
    }
    return;
}

void FFsolverRegula(double (*f)(double), double *x, double *y, int *n, double xmin, double xmax, double *ftol, int *niter)//xmin and xmax define the interval
{
    *n=0;
    double fxmin,fxmax;
    fxmin=(*f)(xmin);
    fxmax=(*f)(xmax);
    if ((fxmin* fxmax)>0){//if both functions have the same sign, it is not sure that a root exists
        double interval=xmax-xmin;
        double xmed,fxmed;
        xmed=xmin+interval*0.5;//we try in the middle point
        fxmed=(*f)(xmed);
        if (fxmin*fxmed<0){
            xmin=xmed;
            fxmin=fxmed;
        }
        else{
            *x=*y=HUGE_VAL;
            return;
        }
    }
    else if ((fxmin* fxmax)==0){//if one of the given value was already the solution
        if (fxmin==0){
            *x=xmin;
            *y=fxmin;
        }
        else{
            *x=xmax;
            *y=fxmax;
        }
        return;
    }
    while (*n<*niter){
        *x=(xmin*fxmax-xmax*fxmin)/(fxmax-fxmin);
        *y=(*f)(*x);
        if ((*y<*ftol)&&(*y>(-*ftol))) break;
        if ((fxmax* *y)>0){
            xmax=*x;
            fxmax=*y;
        }
        else{
            xmin=*x;
            fxmin=*y;
        }
        *n=*n+1;
    }
    return;
}

//Regula Falsi, Illinois modification
void FFsolverRegulaI(double (*f)(double), double *x, double *y, int *n, double xmin, double xmax, double *ftol, int *niter)//xmin and xmax define the interval
{
    *n=0;
    int side=0;
    double fxmin,fxmax;
    fxmin=(*f)(xmin);
    fxmax=(*f)(xmax);
    if ((fxmin* fxmax)>0){//if both functions have the same sign, it is not sure that a root exists
        double interval=xmax-xmin;
        double xmed,fxmed;
        xmed=xmin+interval*0.5;//we try in the middle point
        fxmed=(*f)(xmed);
        if (fxmin*fxmed<0){
            xmin=xmed;
            fxmin=fxmed;
        }
        else{
            *x=*y=HUGE_VAL;
            return;
        }
    }
    else if ((fxmin* fxmax)==0){//if one of the given value was already the solution
        if (fxmin==0){
            *x=xmin;
            *y=fxmin;
        }
        else{
            *x=xmax;
            *y=fxmax;
        }
        return;
    }
    while (*n<*niter){
        *x=(xmin*fxmax-xmax*fxmin)/(fxmax-fxmin);
        *y=(*f)(*x);
        if ((*y<*ftol)&&(*y>(-*ftol))) break;
        if ((fxmax* *y)>0){
            xmax=*x;
            fxmax=*y;
            if (side==-1) fxmin *= 0.5;
            side = -1;
        }
        else{
            xmin=*x;
            fxmin=*y;
            if (side==1) fxmax *= 0.5;
            side = 1;
        }
        *n=*n+1;
    }
    return;
}

//Regula Falsi, Anderson-Bjork modification
void FFsolverRegulaA(double (*f)(double), double *x, double *y, int *n, double xmin, double xmax, double *ftol, int *niter)//xmin and xmax define the interval
{
    *n=0;
    int side=0;
    double fxmin,fxmax;
    fxmin=(*f)(xmin);
    fxmax=(*f)(xmax);
    if ((fxmin* fxmax)>0){//if both functions have the same sign, it is not sure that a root exists
        double interval=xmax-xmin;
        double xmed,fxmed;
        xmed=xmin+interval*0.5;//we try in the middle point
        fxmed=(*f)(xmed);
        if (fxmin*fxmed<0){
            xmin=xmed;
            fxmin=fxmed;
        }
        else{
            *x=*y=HUGE_VAL;
            return;
        }
    }
    else if ((fxmin* fxmax)==0){//if one of the given value was already the solution
        if (fxmin==0){
            *x=xmin;
            *y=fxmin;
        }
        else{
            *x=xmax;
            *y=fxmax;
        }
        return;
    }
    while (*n<*niter){//we begin to seek for the root
        *n=*n+1;
        *x=(xmin*fxmax-xmax*fxmin)/(fxmax-fxmin);//This is the regula falsi method
        *y=(*f)(*x);
        if ((*y<*ftol)&&(*y>(-*ftol))) break;//if we have arrived to the solution exit
        if ((fxmax * *y)>0){//if the proposed solution is of the same sign than f(xmin)
            if (side==-1){//if it happened also in the previous loop
                if ((1-*y/fxmax)>0) fxmin *= (1-*y/fxmax);//we decrease f(xmax) for the next loop calculation
                else fxmin *= 0.5;
            }
            xmax=*x;
            fxmax=*y;
            side = -1;//we register than the solution was of the same sign than previous xmin
        }
        else{//If the prosed solution is of the same sign than f(smax) we apply the same technic
            if (side==1){
                if ((1-*y/fxmin)>0) fxmax *= (1-*y/fxmin);
                else fxmax *= 0.5;
            }
            xmin=*x;
            fxmin=*y;
            side = 1;
        }
    }
    return;
}

void FFsolverNewtonND(double (*f)(double), double *x, double *y, int *n, double xini, double *ftol, int *niter)//xini is the initial guess
{
    *n=0;//loop counter
    double dx,dfx;//differential of independent variable to use, differential of funtion found
    *x=xini;
    *y=(*f)(*x);
    while ((*n<*niter)&&((*y>*ftol)||(*y<(-*ftol)))){
        dx=*x * *ftol / *y;
        dfx=((*f)(*x + dx)- *y);
        *x=*x - *y/dfx*dx;
        *y=(*f)(*x);
        *n=*n+1;
    }
    return;
}


void FFsolverNewton(double (*f)(double),double (*d)(double), double *x, double *y, int *n, double xini, double *ftol, int *niter)//xini is the initial guess
{
    *n=0;//loop counter
    *x=xini;
    *y=(*f)(*x);
    while ((*n<*niter)&&((*y>*ftol)||(*y<(-*ftol)))){
        *x=*x - *y/(*d)(*x);
        *y=(*f)(*x);
        *n=*n+1;
    }
    return;
}


//Interval-Newton monovariable solver.
void FFsolverIntervalNewtonND(double (*f)(double), double x[], double y[], int *n, double xmin, double xmax, double *ftol, int *niter)//xmin and xmax define the interval
{
    *n=0;
    double interval[10][2];//wher to store the intervals
    interval[0][0]=xmin;//we charge the first interval with the initial values
    interval[0][1]=xmax;
    int intCreated=1,intSolved=0;//will count the intervals created and solved
    double xmed;//the middle point of the interval
    double fxmin,fxmax,fxmed;//the value of the function
    double dxmin,dxmax;//the derivatives
    double change;//auxiliar for ordering intervals
    double min,max;//we will store here the interval image
    if (interval[0][0]>interval[0][1]){//we order the initial interval, if necessary
        change=interval[0][0];
        interval[0][0]=interval[0][1];
        interval[0][1]=change;
    }
    xmed=(interval[0][0]+interval[0][1])*0.5;
    fxmed=(*f)(xmed);
    while ((intCreated>intSolved)&&(*n<*niter)){
        printf("\nInitial interval: %f %f\n",interval[intSolved][0],interval[intSolved][1]);
        //printf("hola\n");
        *n=*n+1;
        fxmin=(*f)(interval[intSolved][0]);//We calculate the value of the function and its derivates at extremes
        dxmin=((*f)(interval[intSolved][0]+*ftol)-fxmin)/ *ftol;
        fxmax=(*f)(interval[intSolved][1]);
        dxmax=((*f)(interval[intSolved][1]+*ftol)-fxmax)/ *ftol;
        if (dxmin>dxmax){//we order the interval derivarives
            change=dxmin;
            dxmin=dxmax;
            dxmax=change;
        }
        //printf("dxmin%f dxmax:%f\n",dxmin,dxmax);
        //Here we begin to compute the interval constraint
        if ((dxmin==0)&&(dxmax>0)&&(fxmed>0)){
            min=-HUGE_VAL;
            max=xmed-fxmed/dxmax;
            printf("dxmin=0, fxmed>0\n");
        }
        if ((dxmin==0)&&(dxmax>0)&&(fxmed<0)){
            min=xmed-fxmed/dxmax;
            max=+HUGE_VAL;
            printf("dxmin=0, fxmed<0\n");
        }
        else if ((dxmin<0)&&(dxmax>0)&&(fxmed>0)){//We need to add here a new interval also
            min=xmed-fxmed/dxmin;
            max=HUGE_VAL;
            interval[intCreated][0]=interval[intSolved][0];//It should be -HUGE_VAL, but we make simultaneously the intersection with the original;
            if ((xmed-fxmed/dxmax) < interval[intSolved][1]) interval[intCreated][1]=xmed-fxmed/dxmax;//It should be xmed-fxmed/dxmax
            else interval[intCreated][1]=interval[intSolved][1];
            printf("\n0 between derivatives, fxmed>0\n");
            printf("New interval: %f %f\n\n",interval[intCreated][0],interval[intCreated][1]);
            //printf("%f %f\n",interval[intCreated][0],interval[intCreated][1]);
            intCreated++;
        }
        else if ((dxmin<0)&&(dxmax>0)&&(fxmed<0)){//We need to add here a new interval also
            min=-HUGE_VAL;
            max=xmed-fxmed/dxmin;
            if ((xmed-fxmed/dxmax) > interval[intSolved][0]) interval[intCreated][0]=xmed-fxmed/dxmax;//It should be xmed-fxmed/dxmax
            else interval[intCreated][0]=interval[intSolved][0];
            interval[intCreated][1]=interval[intSolved][1];//It should be HUGE_VAL, but we make simultaneously the intersection with the original;
            printf("\n0 between derivatives, fxmed<0\n");
            printf("New interval: %f %f\n\n",interval[intCreated][0],interval[intCreated][1]);
            //printf("%f %f\n",interval[intCreated][0],interval[intCreated][1]);
            intCreated++;
        }
        else if ((dxmin<0)&&(dxmax==0)&&(fxmed>0)){
            min=xmed-fxmed/dxmin;
            max=-HUGE_VAL;
            printf("dxmax=0, fxmed>0\n");
        }
        else if ((dxmin<0)&&(dxmax==0)&&(fxmed>0)){
            min=-HUGE_VAL;
            max=xmed-fxmed/dxmin;
            printf("dxmax=0, fxmed<0\n");
        }
        else{
            min=xmed-fxmed/dxmax;
            max=xmed-fxmed/dxmin;
        }
        if (min>max){//we order the constraint interval
            change=min;
            min=max;
            max=change;
        }
        printf("xmed:%f fxmed:%f dxmin:%f dxmax:%f\n",xmed,fxmed,dxmin,dxmax);
        printf("Constraint:%f %f\n",min,max);
        if ((max<interval[intSolved][0])||(interval[intSolved][1]<min)){//if there is no intersection, there is no solution
            x[intSolved]=y[intSolved]=HUGE_VAL;
            intSolved++;
            printf("No hay solution en este intervalo\n");
        }
        /*if((interval[intSolved][0]>min)&&(interval[intSolved][1]>max)){
            interval[intCreated][1]=interval[intSolved][1];
            interval[intSolved][1]=(interval[intSolved][1]-interval[intSolved][0])*0.5;//half the interval
            interval[intCreated][0]=interval[intSolved][1];
            intCreated++;
            printf("\nSpliting the interval\n");
        }*/
        else{
            if (min>interval[intSolved][0]) interval[intSolved][0]=min;//we intersect the intervals, to get the new working interval
            if (max<interval[intSolved][1]) interval[intSolved][1]=max;
        }

        //printf("intersection: %f %f\n",interval[intSolved][0],interval[intSolved][1]);
        xmed=(interval[intSolved][0]+interval[intSolved][1])*0.5;
        fxmed=(*f)(xmed);
        if (fabs(fxmed)<*ftol){
            x[intSolved]=xmed;
            y[intSolved]=fxmed;
            //printf("n:%i xmed:%f fxmed;%f fxmin:%f fxmax:%f dxmin:%f dxmax:%f min:%f max:%f xmin:%f xmax:%f\n",*n,xmed,fxmed,fxmin, fxmax,
            //       dxmin,dxmax,min,max,interval[intSolved][0],interval[intSolved][1]);
            intSolved++;
            if (intCreated>intSolved){//If there are more intervals to test, we need to prepare for the next
                xmed=(interval[intSolved][0]+interval[intSolved][1])*0.5;
                fxmed=(*f)(xmed);
            }
        }
    }
    //*n=intSolved;

}
//Interval-Newton monovariable solver.
void FFsolverIntervalNewtonND2(double (*f)(double), double x[], double y[], int *n, double xmin, double xmax, double *ftol, int *niter)//xmin and xmax define the interval
{
    *n=0;
    double interval[20][2];//where to store the intervals
    interval[0][0]=xmin;//we charge the first interval with the initial values
    interval[0][1]=xmax;
    int intCreated=1,intSolved=0;//will count the intervals created and solved
    double xmed;//the middle point of the interval
    double fxmin,fxmax,fxmed;//the value of the function
    double dxmin,dxmax;//the derivatives
    double change;//auxiliar for ordering intervals
    double min,max;//we will store here the interval constraint
    if (interval[0][0]>interval[0][1]){//we order the initial interval, if necessary
        change=interval[0][0];
        interval[0][0]=interval[0][1];
        interval[0][1]=change;
    }
    xmed=(interval[0][0]+interval[0][1])*0.5;
    fxmed=(*f)(xmed);
    while ((intCreated>intSolved)&&(*n<*niter)){
        printf("\nIteration: %i  Initial interval:%i %f %f\n",*n+1,intSolved,interval[intSolved][0],interval[intSolved][1]);
        //printf("hola\n");
        *n=*n+1;
        fxmin=(*f)(interval[intSolved][0]);//We calculate the value of the function and its derivates at extremes
        dxmin=((*f)(interval[intSolved][0]+*ftol)-fxmin)/ *ftol;
        fxmax=(*f)(interval[intSolved][1]);
        dxmax=((*f)(interval[intSolved][1]+*ftol)-fxmax)/ *ftol;
        if (dxmin>dxmax){//we order the interval derivarives
            change=dxmin;
            dxmin=dxmax;
            dxmax=change;
        }
            //Here we begin to compute the interval constraint
        if ((dxmin==0)||(dxmax==0)){
            interval[intSolved][0]=interval[intSolved][0]+*ftol;
            interval[intSolved][1]=interval[intSolved][1]-*ftol;
            printf("dxmin=0, intermal moved\n");
        }
        else if ((dxmin<0)&&(dxmax>0)){//We need to split the interval
            interval[intCreated][0]=xmed;
            interval[intCreated][1]=interval[intSolved][1];
            interval[intSolved][1]=xmed;
            printf("(dxmin<0)&&(dxmax>0) Interval splitted\n");
            intCreated++;
        }
        else {
            min=xmed-fxmed/dxmax;
            max=xmed-fxmed/dxmin;
            if (min>max){//we order the constraint interval
                change=min;
                min=max;
                max=change;
            }
            if (*n==49){
                printf("%i min:%f max:%f solved[0]:%f solved[1]:%f\n",*n,min,max,interval[intSolved][0],interval[intSolved][1]);
            }
            if ((min>interval[intSolved][0])&& (max<interval[intSolved][1])){//If the constraint is inside there is solution
                interval[intSolved][0]=min;
                interval[intSolved][1]=max;
                printf("We are advancing, contrained interval:%f %f\n",min,max);
            }
            else if((min>=interval[intSolved][1])||(max<=interval[intSolved][0])){//if there is no intersection there is no solution
                 intSolved++;
                 printf("No intersection, interval removed\n");
            }
            else{
                interval[intCreated][0]=xmed;
                interval[intCreated][1]=interval[intSolved][1];
                interval[intSolved][1]=xmed;
                printf("No advance, interval splitted\n");
                intCreated++;
            }

        }

        //printf("xmed:%f fxmed:%f dxmin:%f dxmax:%f\n",xmed,fxmed,dxmin,dxmax);
        xmed=(interval[intSolved][0]+interval[intSolved][1])*0.5;
        fxmed=(*f)(xmed);
        if (fabs(fxmed)<*ftol){
            //x[intSolved]=xmed;
            //y[intSolved]=fxmed;
            printf("n:%i xmed:%f fxmed;%f fxmin:%f fxmax:%f dxmin:%f dxmax:%f min:%f max:%f xmin:%f xmax:%f\n",*n,xmed,fxmed,fxmin, fxmax,
                   dxmin,dxmax,min,max,interval[intSolved][0],interval[intSolved][1]);
            intSolved++;
            if (intCreated>intSolved){//If there are more intervals to test, we need to prepare for the next
                xmed=(interval[intSolved][0]+interval[intSolved][1])*0.5;
                fxmed=(*f)(xmed);
            }
        }
    }
    //*n=intSolved;

}
