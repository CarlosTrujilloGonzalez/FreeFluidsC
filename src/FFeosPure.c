/*

 * FFeosPure.c
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

// contains EOS calculations for pure substances
//==============================================

#include <math.h>
#include <stdio.h>
#include "FFbasic.h"
#include "FFeosPure.h"

const double Av = 6.02214E+23; //molecules/mol
const double kb = 1.3806504E-023;//J/K
const double Pi = 3.141592;
const double R = 8.314472; //Pa·m3/(K·mol)
const double FF_PCSAFTap[7][3]={{0.910563145, -0.308401692, -0.090614835},{0.636128145, 0.186053116, 0.452784281},
        {2.686134789, -2.503004726, 0.596270073},{-26.54736249, 21.41979363, -1.724182913},
        {97.75920878, -65.25588533, -4.130211253},{-159.5915409, 83.31868048, 13.77663187},
        {91.29777408, -33.74692293, -8.672847037}};
const double FF_PCSAFTbp[7][3]={{0.724094694, -0.575549808, 0.097688312},{2.238279186, 0.699509552, -0.255757498},
        {-4.002584948, 3.892567339, -9.155856153},{-21.00357681, -17.21547165, 20.64207597},
        {26.85564136, 192.6722645, -38.80443005},{206.5513384, -161.8264616, 93.62677408},
        {-355.6023561, -165.2076935, -29.66690559}};

//Single substance calculations
//=============================
//=============================


//Single substance, cubic EOS calculations
//========================================

//Calculates a,b,u,w for a cubic EOS. Using critical constants
//------------------------------------------------------------
void CALLCONV FF_FixedParamCubic(const  FF_CubicEOSdata *data, FF_CubicParam *param){
    switch (data->eos)
    {
    case FF_PR76:
    case FF_PR78://Peng-Robinson.
    case FF_PRSV1://Peng-Robinson Stryjek-Vera with 1 extra parameter.
    case FF_PRBM://Peng-Robinson Boston-Mathias with 1 extra parameter. Very similar to Stryjeck-Vera
    case FF_PRMELHEM://Peng-robinson Melhem with 2 extra parameters.
    case FF_PRSOF://Peng-Robinson Stamateria-Olivera-Fuentes with 2 extra parameters
    case FF_PRALMEIDA://Peng-Robinson Almeida with 3 extra parameters.
    case FF_PRMC://Peng-Robinson Mathias-Copeman with 3 extra parameters.
    case FF_PRTWU91://Peng-Robinson Twu with 3 extra parameters.
    case FF_PRTWU95://Peng-robinson Twu with parameters from Tc and w
    case FF_PRFIT4B://Peng Robinson Stryjeck- Vera with Tc,Pc,w,k1 fitted
        //printf("Tc:%f %f %f\n",data->Tc,data->Pc,data->w);
        param->a = 0.457235 * pow(R*data->Tc,2)/ data->Pc;
        param->b = 0.077796 * R * data->Tc / data->Pc;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        if (data->c>0) param->c=data->c;//if we supply a volume correction it is used
        else if ((data->c==0)&&(data->Zc>0)&&(data->Zc<0.45)) param->c=R*data->Tc/data->Pc*(0.1014048-0.3892896*data->Zc);//if volume correction is 0, we apply that of Peneloux
        else param->c=0.0;//a negative number will indicate not to use volume correction
        break;
    case FF_PRvTWU91:
        param->a = param->a = 0.42748 * pow(R*data->Tc,2)/ data->Pc*1.08;//0.599877;
        param->b = data->VdWV*1.55;//data->r*15.17*1.6*1e-6;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        //if (data->c>0) param->c=data->c;
        //else if ((data->c==0)&&(data->Zc>0)&&(data->Zc<0.5)) param->c=R*data->Tc/data->Pc*(0.1014048-0.3892896*data->Zc);
        //else param->c=0;
        param->c=data->c;
        break;
    case FF_PRFIT3://similar to FF_PRSV1, some parameters have no relationship with critical properties. There are 3 adjustable parameters:a, fw, k1
        param->a = data->k1;
        param->b = 0.9*0.0778 * R * data->Tc / data->Pc;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        if (data->c>0) param->c=data->c;
        else if ((data->c==0)&&(data->Zc>0)&&(data->Zc<0.5)) param->c=R*data->Tc/data->Pc*(0.1014048-0.3892896*data->Zc);
        else param->c=0;
        break;
    case FF_PRFIT4://parameters have no relationship with critical properties. There are 4 adjustable parameters:b,a,fw and Tc
        param->a = data->k2;
        param->b = data->k1;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        if (data->c>0) param->c=data->c;
        else if ((data->c==0)&&(data->Zc>0)&&(data->Zc<0.5)) param->c=R*data->Tc/data->Pc*(0.1014048-0.3892896*data->Zc);
        else param->c=0;
        break;
    case FF_SRK://Soaves-Redlich-Kwong.
    case FF_SRKSOF://Soaves-Redlich-Kwong Stamateria-Olivera-Fuentes with 2 extra parameters
    case FF_SRKMC://Soaves-Redlich-Kwong Mathias-Copeman, with 3 extra parameters.
    case FF_SRKTWU91:
        param->a = 0.42748 * pow(R*data->Tc,2)/ data->Pc;
        param->b = 0.08664 * R * data->Tc / data->Pc;
        param->u=1;
        param->w=0;
        if (data->c>0) param->c=data->c;
        else if ((data->c==0)&&(data->w>0)&&(data->w<1)) param->c=R*data->Tc/data->Pc*(0.001569568+0.03577392*data->w);
        else param->c=0;
        break;
    case FF_PRPOL1://Peng-robinson for polymers b=k1*MW, a=1,Theta=k2*MW (a=1)
        param->a=1;
        param->b=data->MW*data->k1;
        param->u=1+pow(2,0.5);//1+2^0.5
        param->w=1-pow(2,0.5);//1-2^0.5
        break;
    case FF_SRKPOL2://Soave-Redlich-Kwong for polymers b=k1*MW, a=1
        param->a=1;
        param->b=data->MW*data->k1;
        param->u=1;
        param->w=0;
        break;

    }
}

//Calculates Theta and its derivatives, given a cubic EOS and T. Using critical and others constants
//--------------------------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_ThetaDerivCubic(const double *T,const  FF_CubicEOSdata *data, FF_CubicParam *param)
{
    double Tr,Tx,dTx,d2Tx,fw,Alpha,dAlpha,d2Alpha;
    Tr = *T / data->Tc;
    Tx = 1-pow(Tr,0.5);//As it is a common calculation, it is better to do it once
    dTx=-0.5*pow(Tr,-0.5)/data->Tc;//This is the 1sr derivative redarding T
    d2Tx=0.25*pow(Tr,-1.5)/pow(data->Tc,2);//and this the 2nd
    double d;//used in FF_PRBM
    switch (data->eos)
    {
    case FF_PR76:
        fw = 0.37464 + 1.54226 * data->w - 0.26992 * data->w * data->w;
        Alpha = pow((1 + fw * Tx),2);
        dAlpha=2*(1 + fw * Tx)*fw*dTx;
        d2Alpha=2*fw*(d2Tx+fw*pow(dTx,2)+fw*Tx*d2Tx);
        break;
    case FF_PR78://Peng-Robinson.
        if (data->w<=0.491) fw = 0.37464 + 1.54226 * data->w - 0.26992 * pow(data->w,2);//0.491
        else fw = 0.379642 + 1.487503 * data->w - 0.164423 * pow(data->w,2)+ 0.016666 * pow(data->w,3);
        Alpha = pow((1 + fw * Tx),2);
        dAlpha=2*(1 + fw * Tx)*fw*dTx;
        d2Alpha=2*fw*(d2Tx+fw*pow(dTx,2)+fw*Tx*d2Tx);
        break;
    case FF_PRSV1://Peng-Robinson Stryjek-Vera with 1 extra parameter.
    case FF_PRFIT4B://Peng Robinson Stryjeck- Vera with Tc,Pc,w,k1 fitted
        fw = 0.378893 + 1.4897153 * data->w - 0.17131848 * pow(data->w,2) + 0.0196554 * pow(data->w,3);
        if (*T < data->Tc){
            Alpha = pow((1 + fw * Tx+ data->k1 * (1 - Tr) * (0.7 - Tr)),2);
            dAlpha=2*pow(Alpha,0.5)*(fw*dTx+data->k1*(2*Tr-1.7)/data->Tc);
            d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*(fw*d2Tx+data->k1*2/pow(data->Tc,2));
        }
        else{
            Alpha = pow((1 + fw * Tx),2);
            dAlpha=2*(1 + fw * Tx)*fw*dTx;
            d2Alpha=2*fw*(d2Tx+fw*pow(dTx,2)+fw*Tx*d2Tx);
        }
        break;
    case FF_PRBM:
        fw = 0.37464 + 1.54226 * data->w - 0.26992 * data->w * data->w;
        if (*T < data->Tc){
            Alpha = pow((1 + fw * Tx- data->k1 * (1 - Tr) * (0.7 - Tr)),2);
            dAlpha=2*pow(Alpha,0.5)*(fw*dTx+data->k1*(2*Tr-1.7)/data->Tc);
            d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*(fw*d2Tx+data->k1*2/pow(data->Tc,2));
        }
        else{
            d=1+fw/2+0.3*data->k1;
            Alpha = pow(exp((1-1/d)*(1-pow(Tr,d))),2);
            dAlpha=-2*Alpha*(d-1)*pow(Tr,d-1)/data->Tc;
            d2Alpha=-2*(d-1)/data->Tc*(dAlpha*pow(Tr,d-1)+Alpha*(d-1)*pow(Tr,d-2)/data->Tc);
        }

        break;
    case FF_PRMELHEM://Peng-robinson Melhem with 2 extra parameters.
        Alpha =exp(data->k1*(1-Tr)+data->k2*pow(Tx,2));
        dAlpha=Alpha*(-data->k1/data->Tc+2*data->k2*Tx*dTx);
        d2Alpha=dAlpha*(-data->k1/data->Tc+2*data->k2*Tx*dTx)+Alpha*2*data->k2*(dTx*dTx+Tx*d2Tx);
        break;
    case FF_PRSOF://Stamateria-Olivera-Fuentes
    case FF_SRKSOF:
        Alpha=Tr*(1+data->k1/(data->k2-1))-data->k1/(data->k2-1)*pow(Tr,2-data->k2);
        dAlpha=(1+data->k1/(data->k2-1))/data->Tc-(2-data->k2)*data->k1/(data->k2-1)*pow(Tr,1-data->k2)/data->Tc;
        d2Alpha=-(2-data->k2)*data->k1*(1-data->k2)/(data->k2-1)*pow(Tr,-data->k2)/data->Tc/data->Tc;
        break;
    case FF_PRMC://Peng-Robinson Mathias-Copeman with 3 extra parameters.
    case FF_SRKMC://Soaves-Redlich-Kwong Mathias-Copeman, with 3 extra parameters.
        if (*T < data->Tc)
        {
            Alpha = pow(1+data->k1*Tx+data->k2*pow(Tx,2)+data->k3*pow(Tx,3),2);
            dAlpha=2*pow(Alpha,0.5)*(data->k1*dTx+data->k2*2*Tx*dTx+data->k3*3*pow(Tx,2)*dTx);
            d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*(data->k1*d2Tx+data->k2*2*(pow(dTx,2)+Tx*d2Tx)+data->k3*3*(2*Tx*pow(dTx,2)+pow(Tx,2)*d2Tx));
        }
        else
        {
            Alpha =pow(1+data->k1*Tx,2);
            dAlpha=2*pow(Alpha,0.5)*data->k1*dTx;
            d2Alpha=0.5*dAlpha/Alpha+2*pow(Alpha,0.5)*data->k1*d2Tx;
        }
        break;
    case FF_PRALMEIDA:
        Alpha=exp(data->k1*(1-Tr)*pow(fabs(1-Tr),data->k3-1)+data->k2*(1/Tr-1));
        dAlpha=Alpha*((-data->k1*(pow(fabs(1-Tr),data->k3-1))-(data->k1*(1-Tr))*(data->k3-1)*pow(fabs(1-Tr),data->k3-2)-data->k2/Tr/Tr)/data->Tc);
        d2Alpha=dAlpha*dAlpha/Alpha+Alpha*(2*data->k1*(data->k3-1)*pow(fabs(1-Tr),data->k3-2)-data->k1*(1-Tr)*(data->k3-1)*(data->k3-2)*pow(fabs(1-Tr),data->k3-3)+
                2*data->k2/Tr/Tr/Tr)/data->Tc/data->Tc;
        break;
    case FF_PRTWU91://Peng-Robinson Twu with 3 extra parameters.
    case FF_PRvTWU91:
    case FF_SRKTWU91://Idem SRK
        Alpha =pow(Tr,data->k3*(data->k2-1))*exp(data->k1*(1-pow(Tr,(data->k2*data->k3))));
        dAlpha=Alpha*(data->k3*(data->k2-1)/Tr-data->k1*data->k2*data->k3*pow(Tr,data->k2*data->k3-1))/data->Tc;
        d2Alpha=dAlpha*dAlpha/Alpha+Alpha*(-data->k3*(data->k2-1)/Tr/Tr-data->k1*data->k2*data->k3*(data->k2*data->k3-1)*pow(Tr,data->k2*data->k3-2))/data->Tc/data->Tc;
        break;
    case FF_PRTWU95://Peng-robinson Twu with parameters from Tc and w
        if (*T <= data->Tc)
        {
            //double Alpha0=pow(Tr,1.19764*(0.924779-1))*exp(0.272838*(1-pow(Tr,1.19764*0.924779)));
            //double Alpha1=pow(Tr,2.46022*(0.792014-1))*exp(0.625701*(1-pow(Tr,2.46022*0.792014)));
            double Alpha0=pow(Tr,-0.0900876784)*exp(0.272838*(1-pow(Tr,1.1075523216)));
            double Alpha1=pow(Tr,-0.5116913169)*exp(0.625701*(1-pow(Tr,1.9485286831)));
            double dAlpha0=Alpha0*(-0.0900876784/Tr-0.3021823603*pow(Tr,0.1075523216))/data->Tc;
            double dAlpha1=Alpha1*(-0.5116913169/Tr-1.2191963455*pow(Tr,0.9485286831))/data->Tc;
            double d2Alpha0=dAlpha0*dAlpha0/Alpha0+Alpha0*(0.0900876784/Tr/Tr-0.0325004144*pow(Tr,-0.8924476784))/data->Tc/data->Tc;
            double d2Alpha1=dAlpha1*dAlpha1/Alpha1+Alpha1*(0.5116913169/Tr/Tr-1.1564427041*pow(Tr,-0.8924476784))/data->Tc/data->Tc;
            Alpha=Alpha0+data->w*(Alpha1-Alpha0);
            dAlpha=dAlpha0+data->w*(dAlpha1-dAlpha0);
            d2Alpha=d2Alpha0+data->w*(d2Alpha1-d2Alpha0);
        }
        else
        {
            double Alpha0=pow(Tr,-0.2*(4.73020-1))*exp(0.373949*(1-pow(Tr,-0.2*4.7302)));
            double Alpha1=pow(Tr,-8.0*(1.24615-1))*exp(0.0239035*(1-pow(Tr,-8.0*1.24615)));
            double dAlpha0=Alpha0*(-0.2*(4.73020-1)/Tr-0.373949*-0.2*4.7302*pow(Tr,-0.2*4.7302-1))/data->Tc;
            double dAlpha1=Alpha1*(-8.0*(1.24615-1)/Tr-0.0239035*-8.0*1.24615*pow(Tr,-8.0*1.24615-1))/data->Tc;
            double d2Alpha0=dAlpha0*dAlpha0/Alpha0+Alpha0*(0.2*(4.73020-1)/Tr/Tr-0.373949*-0.2*4.7302*(-0.2*4.7302-1)*pow(Tr,-0.2*4.7302-2))/data->Tc/data->Tc;
            double d2Alpha1=dAlpha1*dAlpha1/Alpha1+Alpha1*(8.0*(1.24615-1)/Tr/Tr-0.0239035*-8.0*1.24615*(-8.0*1.24615-1)*pow(Tr,-8.0*1.24615-2))/data->Tc/data->Tc;
            Alpha=Alpha0+data->w*(Alpha1-Alpha0);
            dAlpha=dAlpha0+data->w*(dAlpha1-dAlpha0);
            d2Alpha=d2Alpha0+data->w*(d2Alpha1-d2Alpha0);
        }
        break;
    case FF_PRFIT3://FF_PRSV1 with some parameters no related to critical properties, 3 adjustable parameters instead of a, fw and k1 (now k3)
        fw = 0.378893 + 1.4897153 * data->w - 0.17131848 * pow(data->w,2) + 0.0196554 * pow(data->w,3);
        Tr = *T / data->k2;
        Tx = 1-pow(Tr,0.5);//As it is a common calculation, it is better to do it once
        dTx=-0.5*pow(Tr,-0.5)/data->k2;//This is the 1sr derivative redarding T
        d2Tx=0.25*pow(Tr,-1.5)/data->k2/data->k2;//and this the 2ndAlpha = pow((1 + data->k2 * (1-pow(*T/data->k3,0.5)),2);
        Alpha = pow((1 + fw * Tx+ data->k3 * (1 - Tr) * (0.7 - Tr)),2);
        dAlpha=2*pow(Alpha,0.5)*(fw*dTx+data->k3*(2*Tr-1.7)/data->k2);
        d2Alpha=0.5*dAlpha*dAlpha/Alpha+2*pow(Alpha,0.5)*(fw*d2Tx+data->k3*2/pow(data->k2,2));
        break;
    case FF_PRFIT4://PR with parameters no related to critical properties. There are 4 adjustable parameters:b,a,fw and Tc
        fw = 0.37464 + 1.54226 * data->k3 - 0.26992 * data->k3 * data->k3;
        Tr = *T / data->k4;
        Tx = 1-pow(Tr,0.5);//As it is a common calculation, it is better to do it once
        dTx=-0.5*pow(Tr,-0.5)/data->k4;//This is the 1sr derivative redarding T
        d2Tx=0.25*pow(Tr,-1.5)/data->k4/data->k4;//and this the 2ndAlpha = pow((1 + data->k2 * (1-pow(*T/data->k3,0.5)),2);
        Alpha = pow((1 + fw * Tx),2);
        dAlpha=2*(1 + fw * Tx)*fw*dTx;
        d2Alpha=2*fw*(d2Tx+fw*pow(dTx,2)+fw*Tx*d2Tx);
        break;
    case FF_SRK://Soaves-Redlich-Kwong.
        fw = 0.48 + 1.574 * data->w - 0.176 * pow(data->w,2);
        //alfa = pow((1 + fw * (1 - pow(Tr,0.5))),2);
        Alpha = pow((1 + fw * Tx),2);
        dAlpha=2*(1 + fw * Tx)*fw*dTx;
        d2Alpha=2*fw*(fw*pow(dTx,2)+fw*Tx*d2Tx);
        break;
    case FF_PRPOL1://Peng-robinson for polymers b=k1*MW, a=1,Theta=k2*MW (a=1)
        Alpha = data->MW*data->k2;
        dAlpha=0;
        d2Alpha=0;
        break;
    case FF_SRKPOL2://Soave-Redlich-Kwong for polymers b=k1*MW, a=1, Theta=k2*MW*exp(k3*T)
        Alpha = data->MW*data->k2*exp(data->k3* *T);
        dAlpha=Alpha*data->k3;
        d2Alpha=dAlpha*data->k3;
        break;
    }
    param->Theta=param->a*Alpha;
    param->dTheta=param->a*dAlpha;
    param->d2Theta=param->a*d2Alpha;
}

//Arr(reduced residual Helmholtz free energy)and Z calculation for a pure substance, given T and V, according to cubic EOS
//------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrZfromTVcubic(const double *T,const double *V,const  FF_CubicParam *param,double *Arr,double *Z)
{
    double ub=param->u*param->b;
    double wb=param->w*param->b;
    double Veos=*V+param->c;
    *Arr=param->Theta/(param->b*R* *T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(*V/(Veos-param->b));//This is Arr
    *Z=*V/(Veos-param->b)-param->Theta* *V/(R * *T * (Veos + ub)*(Veos + wb));//Z
}

//P calculation from T and V using cubic eos
//-------------------------------------------
EXP_IMP void CALLCONV FF_PfromTVcubic(const double *T,const double *V,const  FF_CubicParam *param,double *P)
{
    double Veos=*V+param->c;
    *P=R* *T/(Veos-param->b)-param->Theta/((Veos+param->u*param->b)*(Veos+param->w*param->b));
}


//V calculation for a pure substance, given T and P, according to cubic EOS. Arr and Z are also given
//----------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_VfromTPcubic(const double *T,const double *P,const  FF_CubicParam *param,const char *option,double resultL[3],double resultG[3],char *state)
{

    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    double A,B,uw,a2,a1,a0,L,M,N,Nsqr,m,phi1,Root[4],Z[4],Veos;
    double ub=param->u*param->b;
    double wb=param->w*param->b;
    A = param->Theta * *P / (R * R * *T * *T);
    B = param->b * *P / (R * *T);
    uw=param->u*param->w;
    a2=(param->u+param->w-1)*B-1;
    a1=A+B*(B*(uw-param->u-param->w)-param->u-param->w);
    a0=-B*(A+B*uw*(1+B));
    //printf("a2:%f a1:%f a0:%f\n",a2,a1,a0);
    L = (3 * a1 - a2*a2) / 3;//P
    M = (2 * a2*a2*a2 - 9 * a2 * a1 + 27 * a0) / 27;//Q
    N = M*M / 4 + L*L*L / 27;//R
    //printf("P:%f Q:%f R:%f\n",L,M,N);
    if (N > 0){
        Nsqr=pow(N,0.5);
        Root[0] = (-M / 2 + Nsqr)/fabs((-M / 2 + Nsqr))*pow(fabs((-M / 2 + Nsqr)),0.333333) + (-M / 2 - Nsqr)/fabs(-M / 2 - Nsqr)*pow((fabs(-M / 2 - Nsqr)),0.333333);//complexity is due to grant a positive base
        Z[0] = Root[0] - a2 / 3;
        //printf("Z:%f\n",Z[0]);
        Veos=Z[0]* R * *T / *P;
        resultL[0]=resultG[0]=Veos-param->c;//this is V
        resultL[1]=resultG[1]=param->Theta/(param->b*R* *T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(resultL[0]/(Veos-param->b));//This is Arr
        resultL[2]=resultG[2]=Z[0]*resultL[0]/Veos;//Z
    }
    else{
        m = 2 * pow(-L / 3,0.5);
        phi1 = acos(3 * M / L / m) / 3;
        Root[1] = m * cos(phi1);
        Root[2] = m * cos(phi1 + 4*M_PI/3);
        Root[3] = m * cos(phi1 + 2*M_PI/3);
        Z[1] = Root[1] - (a2 / 3);
        Z[2] = Root[2] - a2 / 3;
        Z[3] = Root[3] - a2 / 3;
        //printf("Z:%f %f %f\n",Z[1],Z[2],Z[3]);
        if ((Z[1] > Z[2])&& (Z[1] > Z[3])) resultG[2] = Z[1];
        else if (Z[2] > Z[3]) resultG[2] = Z[2];
        else resultG[2] = Z[3];//this is Z gas as P*Veos/(R*T)
        Veos=resultG[2]* R * *T / *P;
        resultG[0]=Veos-param->c;//this is V
        resultG[1]=param->Theta/(param->b*R* *T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(resultG[0]/(Veos-param->b));//This is Arr
        resultG[2]=resultG[2]*resultG[0]/Veos;//This is Z gas as P*V/(R*T)
        if ((Z[1] < Z[2]) && (Z[1] < Z[3]) && (Z[1]>0)) resultL[2] = Z[1];
        else if ((Z[2] < Z[3]) && (Z[2]>0)) resultL[2] = Z[2];
        else if (Z[3]>0) resultL[2] = Z[3];
        else resultL[2]=Z[1];
        Veos=resultL[2]* R * *T / *P;
        //printf("Veos: %f c:%f\n",Veos,param->c);
        resultL[0]=Veos-param->c;//this is V
        resultL[1]=param->Theta/(param->b*R* *T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(resultL[0]/(Veos-param->b));//This is Arr
        resultL[2]=resultL[2]*resultL[0]/Veos;
    }
    *state='b';
}


//V calculation for a pure substance, given T and P, according to cubic EOS. Arr and Z are also given
//----------------------------------------------------------------------------------------------------------------------
void CALLCONV VfromTPcubicOriginal(const double *T,const double *P,const  FF_CubicParam *param,const char *option,double resultL[3],double resultG[3],char *state)
{


    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    double Veos,dP_dVeos,error,dP_dVeosPrev;
    double ub=param->u*param->b;
    double wb=param->w*param->b;
    double maxError=0.00001;
    int i;
    if (*option!='g')//we calculate the liquid phase if not only the gas phase has been asked
    {
        Veos=1.01*param->b;//first approximation for liquid volume
        while ((R* *T/(Veos-param->b)-param->Theta/((Veos+ub)*(Veos+wb))<0)||
               ((-R* *T/pow((Veos-param->b),2)+param->Theta*(Veos+ub+Veos+wb)/pow((Veos+ub)*(Veos+wb),2))>0)) Veos=Veos*1.02;
        error =*P-R* *T/(Veos-param->b)+param->Theta/((Veos+ub)*(Veos+wb));
        dP_dVeos=-R* *T/pow((Veos-param->b),2)+param->Theta*(Veos+ub+Veos+wb)/pow((Veos+ub)*(Veos+wb),2);
        i=1;
        printf("Liquido Inicial:T:%f V:%f dP/dV:%f err:%f\n",*T-273.15,Veos, dP_dVeos,error);
        while ((fabs(error/ *P)> maxError) &&(dP_dVeos <0)&&(Veos>0)&&(i<51))//till error is minimum or we have pass a maximum or minimum
        {
            Veos=Veos+error/dP_dVeos;//Newton method for root finding
            dP_dVeosPrev=dP_dVeos;
            dP_dVeos=-R* *T/pow((Veos-param->b),2)+param->Theta*(Veos+ub+Veos+wb)/pow((Veos+ub)*(Veos+wb),2);
            if ((Veos<=0)||(dP_dVeos>=0))
            {
                Veos=Veos-0.8*error/dP_dVeosPrev;//Newton method for root finding, slowed
                dP_dVeos=-R* *T/pow((Veos-param->b),2)+param->Theta*(Veos+ub+Veos+wb)/pow((Veos+ub)*(Veos+wb),2);
            }

            error =*P-R* *T/(Veos-param->b)+param->Theta/((Veos+ub)*(Veos+wb));
            i++;
            printf("Bucle liquido:%d %f %f %f\n",i,Veos,dP_dVeos,error);
        }
        if ((Veos<=0)||(dP_dVeos>0)||(i>=51))
        {
            resultL[0]=resultL[1]=resultL[2]=0;
        }
        else
        {
            resultL[0]=Veos-param->c;//this is V
            resultL[1]=param->Theta/(param->b*R* *T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(resultL[0]/(Veos-param->b));//This is Arr
            resultL[2]=resultL[0]/(Veos-param->b)-param->Theta* resultL[0]/(R * *T * (Veos + ub)*(Veos + wb));//Z
            *state='l';
        }
    }
    if (*option!='l')//and the gas phase if not only the liquid one has been asked for
    {
        Veos=param->b+R * *T / *P;//first approximation for gas volume
        dP_dVeos=-R* *T/pow((Veos-param->b),2)+param->Theta*(Veos+ub+Veos+wb)/pow((Veos+ub)*(Veos+wb),2);
        error =*P-R* *T/(Veos-param->b)+param->Theta/((Veos+ub)*(Veos+wb));
        printf("Gas Inicial:T:%f V:%f dP/dV:%f err:%f\n",*T-273.15,Veos, dP_dVeos,error);
        i=1;
        while ((fabs(error/ *P)>maxError) &&(dP_dVeos <0)&&(Veos>0)&&(i<51))//till error is minimum or we have pass a maximum or minimum
        {
            Veos=Veos+error/dP_dVeos;//Newton method for root finding
            dP_dVeosPrev=dP_dVeos;
            dP_dVeos=-R* *T/pow((Veos-param->b),2)+param->Theta*(Veos+ub+Veos+wb)/pow((Veos+ub)*(Veos+wb),2);
            if ((Veos<=0)||(dP_dVeos>=0))
            {
                Veos=Veos-0.8*error/dP_dVeosPrev;//Newton method for root finding, slowed
                dP_dVeos=-R* *T/pow((Veos-param->b),2)+param->Theta*(Veos+ub+Veos+wb)/pow((Veos+ub)*(Veos+wb),2);
            }
            dP_dVeos=-R* *T/pow((Veos-param->b),2)+param->Theta*(Veos+ub+Veos+wb)/pow((Veos+ub)*(Veos+wb),2);
            error =*P-R* *T/(Veos-param->b)+param->Theta/((Veos+ub)*(Veos+wb));
            i++;
            printf("Bucle gas:%d %f %f %f\n",i,Veos,dP_dVeos,error);
        }
        if ((dP_dVeos>0)||(Veos<=0))
        {
            resultG[0]=resultG[1]=resultG[2]=0;
            //printf("hola\n");
        }
        else
        {
            resultG[0]=Veos-param->c;
            resultG[1]=param->Theta/(param->b*R* *T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(resultG[0]/(Veos-param->b));//This is Arr
            resultG[2]=resultG[0]/(Veos-param->b)-param->Theta* resultG[0]/(R * *T * (Veos + ub)*(Veos + wb));//Z
            if (*state=='l') *state='b';
            else *state='g';
        }
    }
}


void CALLCONV FF_ArrDerCubic(const double *T,const double *V,const  FF_CubicParam *param,double result[6])
{
    double ub=param->u*param->b;
    double wb=param->w*param->b;
    double Veos=*V+param->c;
    double T2=*T * *T;
    result[0]=param->Theta/(param->b*R* *T*(param->w-param->u))*log((Veos+ub)/(Veos+wb))+log(*V/(Veos-param->b));//This is Arr
    result[1]=1/ *V- 1/(Veos-param->b)+param->Theta/(R * *T * (Veos + ub)*(Veos + wb));//dArr/dV
    result[2]=-1/pow(*V,2)+1/pow(Veos-param->b,2)-param->Theta * (Veos * 2+ub+wb)/(R * *T *pow(Veos+ub,2)*pow(Veos+wb,2));//d2Arr/dV2
    result[3]=log((Veos+ub)/(Veos+wb))*(*T * param->dTheta-param->Theta)/(R*param->b*(param->w-param->u)*T2);//dArr/dT
    result[4]=log((Veos+ub)/(Veos+wb))* (param->d2Theta*T2-2*(*T * param->dTheta-param->Theta))/(R*param->b*(param->w-param->u)*T2* *T);//d2Arr/dT2
    result[5]=(*T * param->dTheta-param->Theta)/(R*T2*(Veos+ub)*(Veos+wb));//d2Arr/dVdT
}


//Single substance FF_PCSAFT EOS calculation
//=======================================

//Auxiliary calculation for FF_ArrZfromTVPCSAFT and calcMixPresFF_PCSAFT
//--------------------------------------------------------------
void CALLCONV FF_calcI1I2(double m,double eta,double I[4])
{
    int i;
    double Aeta[7],Beta[7];
    for (i=0;i<7;i++)
    {
        Aeta[i]=(FF_PCSAFTap[i][0] + (m - 1) / m * FF_PCSAFTap[i][1] + (m - 1) * (m - 2) / pow(m,2) * FF_PCSAFTap[i][2])* pow(eta,i);
        I[0]=I[0]+Aeta[i];
        I[1]=I[1]+Aeta[i]*(i+1);
        Beta[i]=(FF_PCSAFTbp[i][0] + (m - 1) / m * FF_PCSAFTbp[i][1] + (m - 1) * (m - 2) / pow(m,2) * FF_PCSAFTbp[i][2])* pow(eta,i);
        I[2]=I[2]+Beta[i];
        I[3]=I[3]+Beta[i]*(i+1);
    }
}


//Z and Arr calculation for a pure substance, given T and V, according to FF_PCSAFT EOS
//-----------------------------------------------------------------------------------
void CALLCONV FF_ArrZfromTVPCSAFT(const double *T,const double *V,const  FF_SaftEOSdata *data,double *Arr,double *Z)
{
    //static int counter=0;
    //counter++;
    double Vmolecular,rho,d,eta,Zhs,Ahs,ghs,dLghs_drho,Zchain,Achain;
    Vmolecular = *V / Av * 1E+30;//molecular volume in A3
    rho = 1 / Vmolecular;//molecules/A3
    d = data->sigma * (1 - 0.12 * exp(-3 * data->epsilon / *T)); //Hard sphere diameter in A, at given T
    eta = data->m * (Pi * pow(d,3) / 6) / Vmolecular; //Volume fraction filled with hard sphere. The terms in Angstrom cancel, so not necessary to pass to SI

    //Contribution by hard spheres chain.
    Zhs = data->m * (4 * eta - 2 * pow(eta,2)) / pow((1 - eta),3);
    Ahs = data->m*(4 * eta - 3 * pow(eta,2)) / pow((1 - eta),2);
    ghs = (2 - eta) / 2 / pow((1 - eta),3);
    dLghs_drho = (5 * eta / 2 - pow(eta,2)) / (1 - eta) / (1 - 0.5 * eta) * Vmolecular;
    Zchain = (1 - data->m) * dLghs_drho / Vmolecular;
    Achain = (1 - data->m) * log(ghs);

    //contribution by dispersion (attraction between non associated molecules)
    double Z1,C1,C2,Z2,Zdisp,Adisp,I[4]={0.0,0.0,0.0,0.0};;
    FF_calcI1I2(data->m,eta,I);
    Z1 = -2 * Pi / Vmolecular * I[1] * pow(data->m,2) * data->epsilon / *T * pow(data->sigma,3);
    C1 = 1/(1 + data->m * (8 * eta - 2 * pow(eta,2)) / pow((1 - eta),4) + (1 - data->m) * (20 * eta - 27 * pow(eta,2)
            + 12 * pow(eta,3) - 2 * pow(eta,4)) / pow(((1 - eta) * (2 - eta)),2));
    C2 = C1 * (data->m * (-4 * pow(eta,2) + 20 * eta + 8) / pow((1 - eta),5) + (1 - data->m) * (2 * pow(eta,3)
            + 12 * pow(eta,2) - 48 * eta + 40) / pow(((1 - eta) * (2 - eta)),3));
    Z2 = -Pi / Vmolecular * data->m * C1 * (I[3] - C2 * eta * I[2])* pow(data->m,2) * pow((data->epsilon / *T),2) * pow(data->sigma,3);
    Zdisp = Z1 + Z2;
    Adisp = -2 * Pi / Vmolecular * I[0] * pow(data->m,2) * data->epsilon * pow(data->sigma,3) / *T - Pi / Vmolecular * data->m * C1
            * I[2] * pow(data->m,2) * pow((data->epsilon / *T),2) * pow(data->sigma,3);

    //Contribution by molecular association
    double DeltaAB,X[data->nPos+data->nNeg+data->nAcid],Zassoc,Aassoc; //X=[] is fraction of molecules not associated at site i
    double sum;
    if ((data->kAB > 0) && (data->epsilonAB > 0)) //If the molecule has association parameters
    {   //DeltaAB = pow(d,3) * ghs * data->kAB * (exp(data->epsilonAB / *T) - 1);
        DeltaAB = pow(data->sigma,3) * ghs * data->kAB * (exp(data->epsilonAB / *T) - 1);
        //Calculation taking account of number of association sites of the molecule(1=acids,2=alcohol,4=water or diols)
        if (data->nAcid==1){//1A
            X[0]=(-1 + pow((1 + 4 * rho * DeltaAB),0.5)) / (2 * rho * DeltaAB);
            //printf("Delta:%f\n",DeltaAB);
        }
        else if (data->nPos==1 && data->nNeg==1){//2B
            X[0]=X[1]=(-1 + pow((1 + 4 * rho * DeltaAB),0.5)) / (2 * rho * DeltaAB);
        }
        else if (data->nPos==2 && data->nNeg==2)//4C
        {
            X[0]=X[1]=X[2]=X[3]=(-1 + pow((1 + 8 * rho * DeltaAB),0.5)) / (4 * rho * DeltaAB);
        }
        else if ((data->nPos==2 && data->nNeg==1)||(data->nPos==1 && data->nNeg==2)){//3B
            X[0]=X[1]=(-(1 - rho * DeltaAB) + pow((pow((1 + rho * DeltaAB),2) + 4 * rho * DeltaAB),0.5)) / (4 * rho * DeltaAB);
            X[2]=(2 * X[0] - 1);
        }
        else if (data->nAcid==2){//2A
            X[0]=X[1]=(-1 + pow((1 + 8 * rho * DeltaAB),0.5)) / (4 * rho * DeltaAB);
            //printf("Delta:%f\n",DeltaAB);
        }
        else if ((data->nPos==3 && data->nNeg==1)||(data->nPos==1 && data->nNeg==3)){//4B
            X[0]=X[1]=X[2]=(-(1 - 2*rho * DeltaAB) + pow((pow((1 + 2*rho * DeltaAB),2) + 4 * rho * DeltaAB),0.5)) / (6 * rho * DeltaAB);
            X[3]=(3 * X[0] - 2);
        }
        else if (data->nAcid==3){//3A
            X[0]=X[1]=X[2]=(-1 + pow((1 + 12 * rho * DeltaAB),0.5)) / (6 * rho * DeltaAB);
            //printf("Delta:%f\n",DeltaAB);
        }
        else if (data->nAcid==4){//4A
            X[0]=X[1]=X[3]=X[4]=(-1 + pow((1 + 16 * rho * DeltaAB),0.5)) / (8 * rho * DeltaAB);
            //printf("Delta:%f\n",DeltaAB);
        }
        else{
            int i;
            for (i==0;i<data->nPos+data->nNeg+data->nAcid;i++) X[i]=1;
        }
        sum=0;
        int i;
        for (i=0;i<(data->nPos+data->nNeg+data->nAcid);i++){
            //printf("i:%i Xi:%f\n",i,X[i]);
            sum = sum + 1 -X[i];
        }
        //printf("sum:%f dLghs_drho:%f\n",sum,dLghs_drho);
        Zassoc=-0.5*(1+rho*dLghs_drho)*sum;
        Aassoc = (data->nPos+data->nNeg+data->nAcid)/ 2;
        for (i=0; i<(data->nPos+data->nNeg+data->nAcid);i++)
            Aassoc = Aassoc + (log(X[i]) - X[i] / 2);
        //printf("Aassoc:%f Zassoc:%f\n",Aassoc,Zassoc);
    }
    else
    {
        Zassoc=0.0;
        Aassoc=0.0;
    }

    //contribution by polar forces
    double Add2,Add3,Add=0,Zdd=0;
    double Vplus,rhoPlus,Add2Plus,Add3Plus,AddPlus;
    if (data->eos== FF_DPCSAFT_GV){//Gross and Vrabeck model. To be studied in the future. Results show higher effect than expected
        double mEfec,muRed2;
        double ad[3][5]={{0.3043504,-0.1358588,1.4493329,0.3556977,-2.0653308},{0.9534641,-1.8396383,2.0131180,-7.3724958,8.2374135},{-1.1610080,4.5258607,0.9751222,-12.281038,5.9397575}};
        double bd[3][5]={{0.2187939,-1.1896431,1.1626889,0.0,0.0},{-0.5873164,1.2489132,-0.5085280,0.0,0.0},{3.4869576,-14.915974,15.372022,0.0,0.0}};
        double cd[3][5]={{-0.0646774,0.1975882,-0.8087562,0.6902849,0.0},{-0.9520876,2.9924258,-2.3802636,-0.2701261,0.0},{-0.6260979,1.2924686,1.6542783,-3.4396744,0.0}};

        double a[5],b[5],c[5];
        double Jdd2,Jdd3;
        double etaPlus,Jdd2Plus,Jdd3Plus;
        Jdd2=Jdd3=Jdd2Plus=Jdd3Plus=0;
        if (data->m>2) mEfec=2;
        else mEfec=data->m;
        muRed2=pow(data->mu,2)/(data->m*pow(data->sigma,3)*data->epsilon)*1e4/1.3807;
        //cout<<"muRed "<<pow(muRed2,0.5)<<endl;
        int n;
        for (n=0;n<5;n++)
        {
            a[n]=ad[0][n]+ad[1][n]*(mEfec-1)/mEfec+ad[2][n]*(mEfec-1)*(mEfec-2)/mEfec/mEfec;
            b[n]=bd[0][n]+bd[1][n]*(mEfec-1)/mEfec+bd[2][n]*(mEfec-1)*(mEfec-2)/mEfec/mEfec;
            c[n]=cd[0][n]+cd[1][n]*(mEfec-1)/mEfec+cd[2][n]*(mEfec-1)*(mEfec-2)/mEfec/mEfec;
            Jdd2=Jdd2+(a[n]+b[n]*data->epsilon/ *T)*pow(eta,n);
            Jdd3=Jdd3+c[n]*pow(eta,n);
        }
        //Add2=-Pi*rho*1e30*pow(muSI2,2)/(pow((kb*T*data->m),2)*pow((data->sigma*1e-10),3))*Jdd2;//Alternative, gives the same result
        //Add3=-4/3*pow((Pi*rho*1e30),2)*pow(muSI2,3)/pow((kb*T*data->m*data->sigma*1e-10),3)*Jdd3;//Alternative, gives the same result
        Add2=-Pi*rho*1e30*pow(muRed2,2)*(pow((data->epsilon/ *T),2)*pow((data->sigma*1e-10),3))*Jdd2;
        Add3=-4/3*pow((Pi*rho*1e30),2)*pow(muRed2,3)*pow((data->epsilon/ *T),3)*pow((data->sigma*1e-10),6)*Jdd3;
        Add=Add2/(1-Add3/Add2);
        //cout<<"Add2 "<<Add2<<" Add3 "<<Add3<<endl;
        //cout<<"Add GV " <<Add<<endl;
        Vplus=*V*1.000001;
        rhoPlus=Av/Vplus*1e-30;
        etaPlus=data->m * (Pi * pow(d,3) / 6) *rhoPlus;
        for (n=0;n<5;n++)
        {
            Jdd2Plus=Jdd2Plus+(a[n]+b[n]*data->epsilon/ *T)*pow(etaPlus,n);
            Jdd3Plus=Jdd3Plus+c[n]*pow(etaPlus,n);
        }
        //Add2Plus=-Pi*rhoPlus*1e30*pow(muSI2,2)/(pow((kb*T*data->m),2)*pow((data->sigma*1e-10),3))*Jdd2Plus;
        //Add3Plus=-4/3*pow((Pi*rhoPlus*1e30),2)*pow(muSI2,3)/pow((kb*T*data->m*data->sigma*1e-10),3)*Jdd3Plus;
        Add2Plus=-Pi*rhoPlus*1e30*pow(muRed2,2)*(pow((data->epsilon/ *T),2)*pow((data->sigma*1e-10),3))*Jdd2Plus;
        Add3Plus=-4/3*pow((Pi*rhoPlus*1e30),2)*pow(muRed2,3)*pow((data->epsilon/ *T),3)*pow((data->sigma*1e-10),6)*Jdd3Plus;
        AddPlus=Add2Plus/(1-Add3Plus/Add2Plus);
        Zdd=-*V*(AddPlus-Add)/(Vplus- *V);
        //cout<<"Zdd GV "<<Zdd<<endl;
    }
    else if (data->eos==FF_DPCSAFT_JC){	//Jog and Chapman model
        double muSI2,rhoRed,I2,I3,rhoRedPlus,I2Plus,I3Plus;
        if (data->mu>0 && data->mu<10 && data->xp>0 && data->xp<1)
        {
            muSI2=pow((data->mu*3.33564e-30),2)/(4*Pi*8.854e-12);
            Vplus=*V*1.000001;
            rhoPlus=Av/Vplus*1e-30;
            rhoRed=rho*data->m*pow(d,3);
            I2=(1-0.3618*rhoRed-0.3205*rhoRed*rhoRed+0.1078*rhoRed*rhoRed*rhoRed)/pow((1-0.5236*rhoRed),2);
            I3=(1+0.62378*rhoRed-0.11658*rhoRed*rhoRed)/(1-0.59056*rhoRed+0.20059*rhoRed*rhoRed);
            Add2=-2*Pi*rho*1e30*pow((muSI2*data->m*data->xp/kb/ *T),2)/(9*pow((d*1e-10),3))*I2;
            Add3=5*pow((Pi*rho*1e30),2)*pow((muSI2*data->m*data->xp/kb/ *T),3)/(162*pow((d*1e-10),3))*I3;
            Add=Add2/(1-Add3/Add2);
            rhoRedPlus=rhoPlus*data->m*pow(d,3);
            I2Plus=(1-0.3618*rhoRedPlus-0.3205*rhoRedPlus*rhoRedPlus+0.1078*rhoRedPlus*rhoRedPlus*rhoRedPlus)/pow((1-0.5236*rhoRedPlus),2);
            I3Plus=(1+0.62378*rhoRedPlus-0.11658*rhoRedPlus*rhoRedPlus)/(1-0.59056*rhoRedPlus+0.20059*rhoRedPlus*rhoRedPlus);
            Add2Plus=-2*Pi*rhoPlus*1e30*pow((muSI2*data->m*data->xp/kb/ *T),2)/(9*pow((d*1e-10),3))*I2Plus;
            Add3Plus=5*pow((Pi*rhoPlus*1e30),2)*pow((muSI2*data->m*data->xp/kb/ *T),3)/(162*pow((d*1e-10),3))*I3Plus;
            AddPlus=Add2Plus/(1-Add3Plus/Add2Plus);
            Zdd=- *V*(AddPlus-Add)/(Vplus- *V);
        }
    }
    *Arr= Ahs + Achain + Adisp + Aassoc+Add;//Reduced residual Helmholz free energy
    *Z = 1 + Zhs + Zchain + Zdisp + Zassoc+Zdd;//Z
}

//P calculation from T and V using according to FF_PCSAFT EOS
//--------------------------------------------------------
EXP_IMP void CALLCONV FF_PfromTVPCSAFT(const double *T,const double *V,const  FF_SaftEOSdata *data,double *P)
{
    double Arr,Z;
    FF_ArrZfromTVPCSAFT(T,V,data,&Arr,&Z);
    *P=Z*R* *T/ *V;
}


//Z,Arr,V calculation for a pure substance, given T and P, according to FF_PCSAFT EOS
//--------------------------------------------------------------------------------
void CALLCONV FF_VfromTPPCSAFT(const double *T,const double *P,const  FF_SaftEOSdata *data,const char *option,double resultL[3],double resultG[3],char *state)
{
    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    double V,Arr,Z,Vplus,ArrPlus,Zplus,error,errorRel,dP_dV;
    int i;
    double maxError=0.000001;
    double d,eta,Vmolecular;
    d = data->sigma * (1 - 0.12 * exp(-3 * data->epsilon / *T)); //Hard sphere diameter, at given T
    eta = 0.6;  //As initial guess, we suppose the volume fraction occupied by the hard molecule is 0.6
    Vmolecular = data->m * (Pi * pow(d,3) / 6) / eta * Av / 1E+30; //This will be the initial guess for the molecular volume in m3/mol
    if (*option!='g')//we calculate the liquid phase if not only the gas phase has been asked
    {
        V = Vmolecular; //initial guess for liquid mole volume in m3
        Vplus=V*1.0000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_ArrZfromTVPCSAFT(T,&V,data,&Arr,&Z);
        FF_ArrZfromTVPCSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
        dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        while (dP_dV>=0){//Initial loop to find a negative derivative of P regarding V
            //printf("finding a liquid positive derivative\n");
            V = V*1.5;
            Vplus=V*1.0000001;
            FF_ArrZfromTVPCSAFT(T,&V,data,&Arr,&Z);
            FF_ArrZfromTVPCSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        }
        error =*P-R* *T*Z/V;
        errorRel=error/ *P;
        i=1;
        //printf("Liquido Inicial:T:%f V:%f dP/dV:%f err:%f\n",*T-273.15,V,dP_dV,error);
        while ((fabs(errorRel)>maxError)&&(dP_dV <0)&&(V>0)&&(i<51))
        {
            V=V+error/dP_dV;//Newton method for root finding
            if ((V<=0)||(dP_dV>0))//We slow the Newton method if V is negative, or dP/dV is positive
            {
               V=V-0.9*error/dP_dV;
            }

            Vplus=V*1.0000001;//we calculate dP/dV numerically
            FF_ArrZfromTVPCSAFT(T,&V,data,&Arr,&Z);
            FF_ArrZfromTVPCSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus- V);
            error =*P-R* *T*Z/V;
            errorRel=error/ *P;
            //printf("i:%d Vl:%f dP_dV:%f error:%f\n",i,V,dP_dV,error);
            i++;
        }
        if ((fabs(errorRel)<maxError)){
            resultL[0]=V;
            resultL[1]=Arr;
            resultL[2]=Z;
            *state='l';//We inform that we have done the calculation from the liquid end

        }
        else resultL[0]=resultL[1]=resultL[2]=0;
    }
    if (*option!='l')//and the gas phase if not only the liquid one has been asked for
    {
        V = R * *T / *P + Vmolecular;//initial guess for gas volume
        Vplus=V*1.000001; //Vplus will mind a volume which corresponding pressure is lower than target pressure.
        FF_ArrZfromTVPCSAFT(T,&V,data,&Arr,&Z);
        FF_ArrZfromTVPCSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
        dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus-V);
        error =*P-R* *T*Z/V;
        errorRel=error/ *P;
        i=1;
        //printf("Gas Inicial:T:%f V:%f dP/dV:%f err:%f\n",*T-273.15,V,dP_dV,error);
        while ((fabs(errorRel)>maxError)&&(dP_dV <0)&&(V>0)&&(i<51))
        {
            V=V+error/dP_dV;//Newton method for root finding
            if ((V<=0)||(dP_dV>0))//We slow the Newton method if V is negative, or dP/dV is positive
            {
               V=V-0.9*error/dP_dV;//Newton method for root finding, slowed
            }

            Vplus=V*1.000001;
            FF_ArrZfromTVPCSAFT(T,&V,data,&Arr,&Z);
            FF_ArrZfromTVPCSAFT(T,&Vplus,data,&ArrPlus,&Zplus);
            dP_dV=R* *T*(Zplus/Vplus-Z/ V)/(Vplus- V);
            error =*P-R* *T*Z/V;
            errorRel=error/ *P;
            //printf("Bucle gas: i:%i V:%f dP/dV:%f err:%f\n",i,V,dP_dV,error);
            i=i+1;
        }
        if (fabs(errorRel)<maxError){
            resultG[0]=V;
            resultG[1]=Arr;
            resultG[2]=Z;
            if (*state=='l') *state='b';
            else *state='g';
            //printf("hola\n");
        }
        else resultG[0]=resultG[1]=resultG[2]=0;
    }
    //printf("state:%c\n",*state);
}



//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given T and V, according to FF_PCSAFT EOS
//--------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrDerPCSAFT(const double *T,const double *V,const  FF_SaftEOSdata *data,double result[6])
{
    double dV=*V * 0.00001;//increments of V and T used to obtain dArr/dV,dArr/dT and d2Arr/dT2 in SAFT eos
    double Vplus=*V + dV;
    double dT=0.01;
    double Tplus=*T+dT;
    double Tminus=*T-dT;
    double Arr,Z,ArrVplus,ZVplus,ArrTplus,ZTplus,ArrTminus,ZTminus;
    FF_ArrZfromTVPCSAFT(T,V,data,&Arr,&Z);
    FF_ArrZfromTVPCSAFT(T,&Vplus,data,&ArrVplus,&ZVplus);
    FF_ArrZfromTVPCSAFT(&Tplus,V,data,&ArrTplus,&ZTplus);
    FF_ArrZfromTVPCSAFT(&Tminus,V,data,&ArrTminus,&ZTminus);
    result[0]=Arr;//This is Arr
    result[1]=(1- Z)/ *V;//dArr/dV at constant T
    result[2]=((1- ZVplus)/ Vplus-result[1])/dV;//d2Arr/dV2 at constant T
    result[3]=(ArrTplus-Arr)/dT;//dArr/T at constant V
    result[4]=(result[3]-(Arr-ArrTminus)/dT)/dT;//d2Arr/dT2 at constant V
    result[5]=((1- ZTplus)/ *V-result[1])/dT;//d2Arr/dVdT
}

//Single substance, SW EOS calculation
//====================================

//Arr and dArr/dV at constant T calculation for a pure substance, given T and V, according to Span and Wagner EOS
//---------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrZfromTVsw(const double *T,const double *V,const  FF_SWEOSdata *data,double *Arr,double *Z)
{
    double dArr_ddelta=0;
    *Arr=0;
    double tau=data->tRef/ *T;
    double delta=1/(*V*data->rhoRef);//rhoRef en mol/m3
    double partial;
    int i,j;
    double d12,D,z,p;
    double dp_ddelta,dD_ddelta,dDb_ddelta;
    for (i=0;i<data->nPol;i++)
    {
        *Arr=*Arr+data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i]);//This will be Arr
        dArr_ddelta=dArr_ddelta+data->d[i]*data->n[i]*pow(delta,data->d[i]-1)*pow(tau,data->t[i]);//And this dArr/ddelta at constant T
    }
    for (i=data->nPol;i<(data->nPol+data->nExp);i++)
    {
        *Arr=*Arr+data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-pow(delta,data->c[i]));
        dArr_ddelta=dArr_ddelta+data->n[i]*pow(delta,data->d[i]-1)*pow(tau,data->t[i])*exp(-pow(delta,data->c[i]))*(data->d[i]-data->c[i]*pow(delta,data->c[i]));
    }
    for (i=(data->nPol+data->nExp);i<(data->nPol+data->nExp+data->nSpec);i++)
    {
        j=i-data->nPol-data->nExp;
        //printf("j:%i\n",j);
        partial=data->n[i]*pow(delta,data->d[i])*pow(tau,data->t[i])*exp(-data->a[j]*pow((delta-data->e[j]),2)-data->b[j]*pow(tau-data->g[j],2));
        *Arr=*Arr+partial;
        dArr_ddelta=dArr_ddelta+(data->d[i]-2* delta*data->a[j]*(delta-data->e[j]))*partial/delta;
    }
    for (i=(data->nPol+data->nExp+data->nSpec);i<(data->nPol+data->nExp+data->nSpec+data->nFinal);i++)//additional terms
    {
        j=i-data->nPol-data->nExp-data->nSpec;
        d12=pow((delta-1),2);
        //Arr[j]=n*D^bf*delta*p with D=z^2+Bf*(delta-1)^(2*af) ; z=(1-tau)+Af*(delta-1)^(1/betaf) ; p=exp(-Cf*(delta-1)^2-Df*(tau-1)^2)
        p=exp(-data->Cf[j]*d12-data->Df[j]*pow((tau-1),2));
        z=(1-tau)+data->Af[j]*pow(d12,(1/2/data->betaf[j]));
        D=pow(z,2)+data->Bf[j]*pow(d12,data->af[j]);
        *Arr=*Arr+data->n[i]*pow(D,data->bf[j])* delta*p;//Arr
        dp_ddelta=-2*data->Cf[j]*(delta-1)*p;
        dD_ddelta=(delta-1)*(data->Af[j]*z*2/data->betaf[j]*pow(d12,(1/2/data->betaf[j]-1))+2*data->Bf[j]*data->af[j]*pow(d12,(data->af[j]-1)));
        dDb_ddelta=data->bf[j]*pow(D,(data->bf[j]-1))*dD_ddelta;
        dArr_ddelta=dArr_ddelta+data->n[i]*(pow(D,data->bf[j])*(p+ delta*dp_ddelta)+dDb_ddelta* delta*p);//dArr/ddelta
    }
    *Z=1+delta*dArr_ddelta;//this is Z
}

//P calculation from T and V using according to Span and Wagner EOS
//-----------------------------------------------------------------
EXP_IMP void CALLCONV FF_PfromTVsw(const double *T,const double *V,const  FF_SWEOSdata *data,double *P)
{
    double Arr,Z;
    FF_ArrZfromTVsw(T,V,data,&Arr,&Z);
    *P=R* *T*Z / *V;//This is P in Pa
}

//P and dP/drho calculation for a pure substance, given reduced T and rho, according to SW EOS
//--------------------------------------------------------------------------------------------
void CALLCONV FF_PresDerSW(const double *tau,const double *delta,const  FF_SWEOSdata *data,double result[2])
{
    double partial,der[2]={0,0};
    int i,j;
    double d12,D,z,p;
    double dp_ddelta,dD_ddelta,dDb_ddelta;
    double d2p_ddelta2,d2D_ddelta2,d2Db_ddelta2;
    for (i=0;i<data->nPol;i++)
    {
        partial=data->n[i]*pow(*delta,data->d[i])*pow(*tau,data->t[i]);
        der[0]=der[0]+data->d[i]*partial;//This will be delta*dArr/ddelta at constant tau
        der[1]=der[1]+data->d[i]*(data->d[i]-1)*partial;//This will be delta*delta*d2Arr/ddelta2 at constant tau
    }
    for (i=data->nPol;i<(data->nPol+data->nExp);i++)
    {
        partial=data->n[i]*pow(*delta,data->d[i])*pow(*tau,data->t[i])*exp(-pow(*delta,data->c[i]));
        der[0]=der[0]+(data->d[i]-data->c[i]*pow(*delta,data->c[i]))*partial;
        der[1]=der[1]+((data->d[i]-data->c[i]*pow(*delta,data->c[i]))*(data->d[i]-1-data->c[i]*pow(*delta,data->c[i]))-data->c[i]*data->c[i]*pow(*delta,data->c[i]))
                *partial;
    }
    for (i=(data->nPol+data->nExp);i<(data->nPol+data->nExp+data->nSpec);i++)
    {
        j=i-data->nPol-data->nExp;
        partial=data->n[i]*pow(*delta,data->d[i])*pow(*tau,data->t[i])*exp(-data->a[j]*pow((*delta-data->e[j]),2)-data->b[j]*pow(*tau-data->g[j],2));
        der[0]=der[0]+(data->d[i]-2* *delta*data->a[j]*(*delta-data->e[j]))*partial;
        der[1]=der[1]+(pow(data->d[i]-2* *delta*data->a[j]*(*delta-data->e[j]),2)-data->d[i]-2* *delta* *delta*data->a[j])*partial;
    }
    for (i=(data->nPol+data->nExp+data->nSpec);i<(data->nPol+data->nExp+data->nSpec+data->nFinal);i++)//additional terms
    {
        j=i-data->nPol-data->nExp-data->nSpec;
        d12=pow((*delta-1),2);
        //Arr[j]=n*D^bf*delta*p with D=z^2+Bf*(delta-1)^(2*af) ; z=(1-tau)+Af*(delta-1)^(1/betaf) ; p=exp(-Cf*(delta-1)^2-Df*(tau-1)^2)
        p=exp(-data->Cf[j]*d12-data->Df[j]*pow((*tau-1),2));
        z=(1-*tau)+data->Af[j]*pow(d12,(1/2/data->betaf[j]));
        D=pow(z,2)+data->Bf[j]*pow(d12,data->af[j]);
        //result[0]=result[0]+data->n[i]*pow(D,data->bf[j])* *delta*p;//Arr
        dp_ddelta=-2*data->Cf[j]*(*delta-1)*p;
        dD_ddelta=(*delta-1)*(data->Af[j]*z*2/data->betaf[j]*pow(d12,(1/2/data->betaf[j]-1))+2*data->Bf[j]*data->af[j]*pow(d12,(data->af[j]-1)));
        dDb_ddelta=data->bf[j]*pow(D,(data->bf[j]-1))*dD_ddelta;
        der[0]=der[0]+ *delta*(data->n[i]*(pow(D,data->bf[j])*(p+ *delta*dp_ddelta)+dDb_ddelta* *delta*p));//delta*dArr/ddelta
        d2p_ddelta2=2*data->Cf[j]*p*(2*data->Cf[j]*d12-1);
        d2D_ddelta2=dD_ddelta/(*delta-1)+d12*(4*data->Bf[j]*data->af[j]*(data->af[j]-1)*pow(d12,(data->af[j]-2))+2*pow((data->Af[j]/data->betaf[j]),2)*
                    pow(d12,(1/data->betaf[j]-2))+data->Af[j]*z*4/data->betaf[j]*(1/2/data->betaf[j]-1)*pow(d12,(1/2/data->betaf[j]-2)));
        d2Db_ddelta2=data->bf[j]*(pow(D,(data->bf[j]-1))*d2D_ddelta2+(data->bf[j]-1)*pow(D,(data->bf[j]-2))*pow(dD_ddelta,2));
        der[1]=der[1]+ *delta * *delta*(data->n[i]*(pow(D,data->bf[j])*(2*dp_ddelta+ *delta*d2p_ddelta2)+2*dDb_ddelta*(p+ *delta*dp_ddelta)+
               d2Db_ddelta2* *delta*p));//delta*delta*d2Arr/ddelta2
    }
    result[0]=R * (data->tRef/ *tau)*(1+der[0])* *delta*data->rhoRef;//This is P in Pa
    result[1]=R* (data->tRef/ *tau)*(1+2*der[0]+der[1])*data->rhoRef;//this is dP/ddelta
    //printf("Der:%f Der2:%f\n",der[0],der[1]);
    //printf("P:%f dP/drho:%f\n",result[0],result[1]);

}


//Z,Arr,V and dP/dV at constant T calculation for a pure substance, given T and P, according to SW EOS
//----------------------------------------------------------------------------------------------------
void CALLCONV FF_VfromTPsw(const double *T,const double *P,const  FF_SWEOSdata *data,const char *option,double resultL[3],double resultG[3],char *state)
{
    *state='f';//We beging puting calculation state information to fail. If calculation finish OK we will change this information
    double tau=data->tRef/ *T;
    //printf("tau: %f \n",tau);
    double delta,error;
    double answer[2];
    int i;
    double maxError;
    if (*T>0.95*data->Tc) maxError=0.0001;
    else maxError=0.00001;
    if (*option!='g')//we calculate the liquid phase if not only the gas phase has been asked
    {
        delta=1/(data->rhoRef*0.0778*R*data->Tc/data->Pc);//we use b calculated as per Peng-Robinson to initiate V
        FF_PresDerSW(&tau,&delta,data,answer);
        error =*P-answer[0];

        i=1;
        //printf("i= %d  error= %f  dP_ddelta= %f\n",i,error,answer[1]);
        while ((fabs(error/ *P)>maxError) && (answer[1] >  0)&&(delta>1))
        {
            if (delta<1.3) delta=delta+error/answer[1]*0.5;
            else delta=delta+error/answer[1];//Newton method for root finding
            FF_PresDerSW(&tau,&delta,data,answer);
            error =error =*P-answer[0];
            i=i+1;
            //printf("i= %d  deltaL= %f error= %f  P= %f  dP_ddelta= %f\n",i,delta,error,answer[0],answer[1]);
        }
        if ((answer[1]>=0)&&(delta>1))//if dP/ddelta >=0
        {
            resultL[0]=1/(delta*data->rhoRef); //This is V
            FF_ArrZfromTVsw(T,&resultL[0],data,&resultL[1],&resultL[2]);
            *state='l';
        }
        else
        {
            resultL[0]=resultL[1]=resultL[2]=0;
        }
    }

    if (*option!='l')//and the gas phase if not only the liquid one has been asked for
    {
        delta=*P * tau/(R* data->tRef*data->rhoRef);//initial guess for gas volume
        //printf("tau= %f  delta= %f \n",tau,delta);
        FF_PresDerSW(&tau,&delta,data,answer);
        error =*P-answer[0];
        i=1;
        //printf("i= %d  error= %f  dP_ddelta= %f\n",i,error,answer[1]);
        while ((fabs(error/ *P)>maxError) && (answer[1] >  0)&&(delta<1))
        {
            if (delta>0.75) delta=delta+error/answer[1]*0.5;
            //else if (delta>0.6) delta=delta+error/answer[1]*0.5;
            else delta=delta+error/answer[1];//Newton method for root finding
            FF_PresDerSW(&tau,&delta,data,answer);
            error =error =*P-answer[0];
            i=i+1;
            //printf("i= %d  deltaG= %f error= %f  P= %f  dP_ddelta= %f\n",i,delta,error,answer[0],answer[1]);
        }
        if ((answer[1]>=0)&&(delta<1))
        {
            resultG[0]=1/(delta*data->rhoRef); //This is V
            FF_ArrZfromTVsw(T,&resultG[0],data,&resultG[1],&resultG[2]);
            if (*state=='l') *state='b';
            else *state='g';
        }
        else
        {
            resultG[0]=resultG[1]=resultG[2]=0;
        }
    }

}

//Arr (reduced residual Helmholtz energy) and its partial derivatives calculation for a pure substance, given tau and delta, according to SW EOS
//----------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ArrDerSW(const double *tau,const double *delta,const  FF_SWEOSdata *data,double result[6])
{
    int i,j;
    double partial;
    double d12,D,z,p;
    double dp_ddelta,dD_ddelta,dDb_ddelta;
    double d2p_ddelta2,d2D_ddelta2,d2Db_ddelta2;
    double dp_dtau,dDb_dtau,d2p_dtau2,d2Db_dtau2;
    double d2p_ddelta_dtau,d2Db_ddelta_dtau;
    for (i=0;i<6;i++) result[i]=0;
    for (i=0;i<data->nPol;i++)//Polinomial terms
    {
        partial=data->n[i]*pow(*delta,data->d[i])*pow(*tau,data->t[i]);//n*delta^d*tau^t
        result[0]=result[0]+partial;//This will be Arr
        result[1]=result[1]+data->d[i]*partial/ *delta;//This will be dArr/ddelta at constant tau
        result[2]=result[2]+data->d[i]*(data->d[i]-1)*partial/ *delta / *delta;//This will be d2Arr/ddelta2 at constant tau
        result[3]=result[3]+data->t[i]*partial/ *tau;//This will be dArr/dtau at constant delta
        result[4]=result[4]+data->t[i]*(data->t[i]-1)*partial/ *tau/ *tau;//This will be d2Arr/dtau2 at constant delta
        result[5]=result[5]+data->t[i]*data->d[i]*partial/ *delta / *tau;//This will be d2Arr/ddelta_dteta
        //result[6]=result[6]+data->d[i]*(data->d[i]-1)*(data->d[i]-2)*partial/ *delta / *delta/ *delta;//This will be d3Arr/ddelta3 at constant tau
    }
    for (i=data->nPol;i<(data->nPol+data->nExp);i++)//exponential terms
    {
        partial=data->n[i]*pow(*delta,data->d[i])*pow(*tau,data->t[i])*exp(-pow(*delta,data->c[i]));//n*delta^d*tau^t*e^(-delta^c)
        result[0]=result[0]+partial;
        result[1]=result[1]+(data->d[i]-data->c[i]*pow(*delta,data->c[i]))*partial/ *delta;
        result[2]=result[2]+((data->d[i]-data->c[i]*pow(*delta,data->c[i]))*(data->d[i]-1-data->c[i]*pow(*delta,data->c[i]))-data->c[i]*data->c[i]*pow(*delta,data->c[i]))
                *partial/ *delta/ *delta;
        result[3]=result[3]+data->t[i]*partial/ *tau;
        result[4]=result[4]+data->t[i]*(data->t[i]-1)*partial/ *tau/ *tau;
        result[5]=result[5]+data->t[i]*(data->d[i]-data->c[i]*pow(*delta,data->c[i]))*partial/ *delta/ *tau;
        //result[6]=result[6]+
    }
    for (i=(data->nPol+data->nExp);i<(data->nPol+data->nExp+data->nSpec);i++)//special exponential terms
    {
        j=i-data->nPol-data->nExp;
        partial=data->n[i]*pow(*delta,data->d[i])*pow(*tau,data->t[i])*exp(-data->a[j]*pow((*delta-data->e[j]),2)-data->b[j]*pow(*tau-data->g[j],2));
        //n*delta^d*tau^t*e^(-a*(delta-e)^2-b*(tau-g)^2)
        result[0]=result[0]+partial;
        result[1]=result[1]+(data->d[i]-2* *delta*data->a[j]*(*delta-data->e[j]))*partial/ *delta;
        result[2]=result[2]+(pow(data->d[i]-2* *delta*data->a[j]*(*delta-data->e[j]),2)-data->d[i]-2* *delta* *delta*data->a[j])*partial/ *delta/ *delta;
        result[3]=result[3]+(data->t[i]-2* *tau*data->b[j]*(*tau-data->g[j]))*partial/ *tau;
        result[4]=result[4]+(pow(data->t[i]-2* *tau*data->b[j]*(*tau-data->g[j]),2)-data->t[i]-2* *tau* *tau*data->b[j])*partial/ *tau/ *tau;
        result[5]=result[5]+(data->d[i]-2* *delta*data->a[j]*(*delta-data->e[j]))*(data->t[i]-2* *tau*data->b[j]*(*tau-data->g[j]))*partial/ *delta/ *tau;
    }
    for (i=(data->nPol+data->nExp+data->nSpec);i<(data->nPol+data->nExp+data->nSpec+data->nFinal);i++)//additional terms
    {
        j=i-data->nPol-data->nExp-data->nSpec;
        d12=pow((*delta-1),2);
        //Arr[j]=n*D^bf*delta*p with D=z^2+Bf*(delta-1)^(2*af) ; z=(1-tau)+Af*(delta-1)^(1/betaf) ; p=exp(-Cf*(delta-1)^2-Df*(tau-1)^2)
        p=exp(-data->Cf[j]*d12-data->Df[j]*pow((*tau-1),2));
        z=(1-*tau)+data->Af[j]*pow(d12,(1/2/data->betaf[j]));
        D=pow(z,2)+data->Bf[j]*pow(d12,data->af[j]);
        result[0]=result[0]+data->n[i]*pow(D,data->bf[j])* *delta*p;//Arr
        dp_ddelta=-2*data->Cf[j]*(*delta-1)*p;
        dD_ddelta=(*delta-1)*(data->Af[j]*z*2/data->betaf[j]*pow(d12,(1/2/data->betaf[j]-1))+2*data->Bf[j]*data->af[j]*pow(d12,(data->af[j]-1)));
        dDb_ddelta=data->bf[j]*pow(D,(data->bf[j]-1))*dD_ddelta;
        result[1]=result[1]+data->n[i]*(pow(D,data->bf[j])*(p+ *delta*dp_ddelta)+dDb_ddelta* *delta*p);//dArr/ddelta
        d2p_ddelta2=2*data->Cf[j]*p*(2*data->Cf[j]*d12-1);
        d2D_ddelta2=dD_ddelta/(*delta-1)+d12*(4*data->Bf[j]*data->af[j]*(data->af[j]-1)*pow(d12,(data->af[j]-2))+2*pow((data->Af[j]/data->betaf[j]),2)*
                    pow(d12,(1/data->betaf[j]-2))+data->Af[j]*z*4/data->betaf[j]*(1/2/data->betaf[j]-1)*pow(d12,(1/2/data->betaf[j]-2)));
        d2Db_ddelta2=data->bf[j]*(pow(D,(data->bf[j]-1))*d2D_ddelta2+(data->bf[j]-1)*pow(D,(data->bf[j]-2))*pow(dD_ddelta,2));
        result[2]=result[2]+data->n[i]*(pow(D,data->bf[j])*(2*dp_ddelta+ *delta*d2p_ddelta2)+2*dDb_ddelta*(p+ *delta*dp_ddelta)+d2Db_ddelta2* *delta*p);//dArr/ddelta
        dp_dtau=-2*data->Df[j]*(*tau-1)*p;
        dDb_dtau=-2*z*data->bf[j]*pow(D,(data->bf[j]-1));
        result[3]=result[3]+data->n[i]* *delta*(dDb_dtau*p+pow(D,data->bf[j])*dp_dtau);//dArr/dtau
        d2p_dtau2=2*data->Df[j]*p*(2*data->Df[j]*pow((* tau-1),2)-1);
        d2Db_dtau2=2*data->bf[j]*pow(D,(data->bf[j]-1))+4*pow(z,2)*data->bf[j]*(data->bf[j]-1)*pow(D,(data->bf[j]-2));
        result[4]=result[4]+data->n[i]* *delta*(d2Db_dtau2*p+2*dDb_dtau*dp_dtau+pow(D,data->bf[j])*d2p_dtau2);//d2Arr/dtau2
    }
        //printf("%f  %f  %f %f %f %f\n",result[0],result[1],result[2],result[3],result[4],result[5]);
}

void CALLCONV FF_ArrDerSWTV(const double *T,const double *V,const  FF_SWEOSdata *data,double result[6])
{
    double delta,tau;
    double ArrDer[6];
    delta=1/(*V * data->rhoRef);
    tau=data->tRef/ *T;
    FF_ArrDerSW(&tau,&delta,data,ArrDer);
    result[0]=ArrDer[0];
    result[1]=-ArrDer[1]/(data->rhoRef* pow(*V,2));
    result[2]=ArrDer[2]/(data->rhoRef* pow(*V,2))/(data->rhoRef* pow(*V,2))+ArrDer[1]*2/(data->rhoRef* pow(*V,3));
    result[3]=-ArrDer[3]*data->tRef/pow(*T,2);
    result[4]=ArrDer[4]*data->tRef/pow(*T,2)*data->tRef/pow(*T,2)+ArrDer[3]*data->tRef*2/pow(*T,3);
    result[5]=ArrDer[5]/(data->rhoRef* pow(*V,2))*data->tRef/pow(*T,2);
}

//Single substance common calculations
//====================================

//Arr, and dArr/dV at constant T, calculation for a pure substance, given T and V by eos
//--------------------------------------------------------------------------------------
void CALLCONV FF_ArrZfromTVeos(const enum FF_EosType *eosType,const double *T,const double *V,const void *data,double *Arr,double *Z)
{
     FF_CubicParam param;
    switch (*eosType)//if we have a cubic eos the first step is to calculate its parameters
    {
        case FF_SAFTtype:
            //*( FF_SaftEOSdata*) data;
            FF_ArrZfromTVPCSAFT(T,V,data,Arr,Z);
            break;
        case FF_SWtype:
            FF_ArrZfromTVsw(T,V,data,Arr,Z);
            break;
        default:
            //*( FF_CubicEOSdata*) data;
            FF_FixedParamCubic(data, &param);
            FF_ThetaDerivCubic(T,data, &param);
            FF_ArrZfromTVcubic(T,V,&param,Arr,Z);
            break;
    }
    //printf("Arr:%f Z:%f")
}

//P calculation from T and V by eos
//---------------------------------
void CALLCONV FF_PfromTVeos(const enum FF_EosType *eosType,const double *T,const double *V,const void *data,double *P)
{
     FF_CubicParam param;
    switch (*eosType)
    {
        case FF_IdealType:
            *P=R* *T/ *V;
            break;
        case FF_SAFTtype:
            //*( FF_SaftEOSdata*) data;
            FF_PfromTVPCSAFT(T,V,data,P);
            break;
        case FF_SWtype:
            FF_PfromTVsw(T,V,data,P);
            break;
        default://if we have a cubic eos the first step is to calculate its parameters
            //*( FF_CubicEOSdata*) data;
            FF_FixedParamCubic(data, &param);
            FF_ThetaDerivCubic(T,data, &param);
            FF_PfromTVcubic(T,V,&param,P);
            break;
    }
}

//V,Arr and Z calculation for a pure substance, given T and P by eos. *state= L,only liquid, G:only gas, U:unique, E:equilibrium, F:fail
//--------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_VfromTPeos(const enum FF_EosType *eosType,const double *T,const double *P,const void *data,const char *option,double resultL[3],double resultG[3],char *state)
{
    //printf("P:%f",*P);
     FF_CubicParam param;
    switch (*eosType)//if we have a cubic eos the first step is to calculate its parameters
    {
        case FF_IdealType:
            resultL[0]=resultG[0]=R* *T/ *P;
            resultL[1]=resultG[1]=0;
            resultL[2]=resultG[2]=1;
            *state='u';
            break;
        case FF_SAFTtype:
            //( FF_SaftEOSdata*) data;
            FF_VfromTPPCSAFT(T,P,data,option,resultL,resultG,state);
            break;
        case FF_SWtype:
            //*( FF_SWEOSdata*) data;
            FF_VfromTPsw(T,P,data,option,resultL,resultG,state);
            break;
        default:
            //( FF_CubicEOSdata*) data;
            FF_FixedParamCubic(data, &param);
            FF_ThetaDerivCubic(T,data, &param);
            //printf("Hola soy V from TP eos: c:%f b:%f a:%f Theta:%f dTheta:%f d2Theta:%f\n",param.c,param.b,param.a,param.Theta,param.dTheta,param.d2Theta);
            FF_VfromTPcubic(T,P,&param,option,resultL,resultG,state);
            break;
    }
    //printf("T:%f  P:%f Vl:%f ArrL:%f Zl:%f\n",*T,*P,resultL[0],resultL[1],resultL[2]);
    if (*option=='s'){
        if ((*state=='b')&&(fabs((resultL[0]-resultG[0])/resultL[0])>0.001)){
                if ((resultL[1]+resultL[2]-1-log(resultL[2]))<(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='L';//we compare Gdr
                else if ((resultL[1]+resultL[2]-1-log(resultL[2]))>(resultG[1]+resultG[2]-1-log(resultG[2]))) *state='G';
                else *state='E';//if Gdr is the same we are in equilibrium
        }
    }
}

//Boiling point calculation
//-------------------------
void CALLCONV FF_TbEOS(const enum FF_EosType *eosType,const double *P,const void *data,double *Tb)
{
    int n=0;//number of calculations don
    //printf("hola soy Tb\n");
    double Tc,Pc;
    switch (*eosType)
    {
        case FF_SAFTtype:
            if ((( FF_SaftEOSdata*)data)->eos==FF_PCSAFTPOL1){
                *Tb=0;
                return;
            }
            Tc=(( FF_SaftEOSdata*)data)->Tc;
            Pc=(( FF_SaftEOSdata*)data)->Pc;
            break;
        case FF_SWtype:
            if (((( FF_SWEOSdata*)data)->eos==FF_IAPWS95)||((( FF_SWEOSdata*)data)->eos==IF97)){
                Pc=22.064e6;
                if (*P>=Pc)
                {
                    *Tb=647.096;
                }
                else
                {
                    double beta=pow(*P * 1e-6,0.25);
                    double B2=beta*beta;
                    double E=B2-0.17073846940092e2*beta+0.1491510861353e2;
                    double F=0.11670521452767e4*B2+0.1202082470247e5*beta-0.48232657361591e4;
                    double G=-0.72421316703206e6*B2-0.32325550322333e7*beta+0.40511340542057e6;
                    double D=2*G/(-F-pow(F*F-4*E*G,0.5));
                    *Tb=(0.65017534844798e3+D-pow((0.65017534844798e3+D)*(0.65017534844798e3+D)-4*(-0.23855557567849+0.65017534844798e3*D),0.5))/2;
                }
                return;
            }
            Tc=(( FF_SWEOSdata*)data)->Tc;
            Pc=(( FF_SWEOSdata*)data)->Pc;
            break;
        default://Cubic eos
            if (((( FF_CubicEOSdata*)data)->eos==FF_PRPOL1)||((( FF_CubicEOSdata*)data)->eos==FF_SRKPOL2)){
                *Tb=0;
                return;
            }
            Tc=(( FF_CubicEOSdata*)data)->Tc;
            Pc=(( FF_CubicEOSdata*)data)->Pc;
            break;
    }
    if ((*P>=Pc)&&(!(Pc==0)))
    {
        *Tb=Tc;
    }
    else
    {
        if (Tc==0) Tc=800.0;
        double T,phiL,phiG;//Temperature and fugacity coef.
        double answerL[3],answerG[3];
        char option='b',state;
        T=Tc/2;//Initial temperature guess
        FF_VfromTPeos(eosType,&T,P,data,&option,answerL,answerG,&state);
        //printf("T:%f  Vliq:%f  Vgas:%f\n",T,answerL[0],answerG[0]);
        int i=4;//i will be used in the calculation of the new temperature

        while (!((answerL[0]>0) && (answerG[0]>0) && ((fabs(answerL[0]-answerG[0])/answerL[0])>0.002)))//Till we find a temperature with different possitive liquid and gas solutions
        {
            if (state=='g')T=T-Tc/i;
            else if (state=='l') T=T+Tc/i;
            if ((state=='b')&&(((fabs(answerL[0]-answerG[0])/answerL[0])<=0.002))){//If we have found only one result
                if (answerL[2]<0.2) T=T+Tc/i;
                else if (answerL[2]>0.5) T=T-Tc/i;
                else//if we are not sure about if it is liquid or gas
                {
                    double Vprov=answerL[0]*0.7,ArrProv,Zprov;
                    FF_ArrZfromTVeos(eosType,&T,&Vprov,data,&ArrProv,&Zprov);
                    //printf("Finding by Z n:%i i:%i  Vl:%f  Zl:%f  Vprov:%f  Zprov:%f\n",n,i/2,answerL[0],answerL[2],Vprov,Zprov);
                    if (Zprov<(answerL[2]-0.06)) T=T-Tc/i;
                    else if (Zprov>(answerL[2]-0.03)) T=T+Tc/i;

                    else {
                        *Tb=T;
                        return;
                    }
                }
            }

            FF_VfromTPeos(eosType,&T,P,data,&option,answerL,answerG,&state);
            i=i*2;
            n=n+1;
            //printf("Finding two phases n:%i i:%i T:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i/2,T,answerL[0],answerL[2],answerG[0],answerG[2]);
            if (n>15){
                *Tb=0;
                return;
            //printf("T:%f  Vliq:%f  Vgas:%f\n",T,answerL[0],answerG[0]);
            }
        }

        phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
        phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
        //printf("T:%f  phiL:%f  phiG:%f\n",T,phiL,phiG);
        while (fabs(phiL-phiG)>0.0001)
        {
            T=T* pow((phiG / phiL),0.06);
            FF_VfromTPeos(eosType,&T,P,data,&option,answerL,answerG,&state);
            phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
            phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
            //printf("%f  %f  %f\n",T,phiL,phiG);
            n=n+1;
            //printf("%i Tb:%f  phiL:%f  phiG:%f\n",n-1,T,phiL,phiG);
            if (n>100){
                *Tb=0;
                return;
            }
        }
        *Tb=T;
    }

}

//Vapor pressure calculation
//--------------------------
void CALLCONV FF_VpEOS(const enum FF_EosType *eosType,const double *T,const void *data,double *Vp)
{
    //printf("Hola, soy Vp, EOS type:%i\n",*eosType);
    int n=0;//number of calculations done
    double phiMaxError;
    double Tc,Pc;
    switch (*eosType)
    {
        case FF_SAFTtype:
            if ((( FF_SaftEOSdata*)data)->eos==FF_PCSAFTPOL1){
                *Vp=0;
                return;
            }
            Tc=(( FF_SaftEOSdata*)data)->Tc;
            Pc=(( FF_SaftEOSdata*)data)->Pc;
            phiMaxError=0.01;
            //printf("Tc:%f Pc:%f sigma:%f m:%f epsilon:%f \n",Tc,Pc,(( FF_SaftEOSdata*)data)->sigma,(( FF_SaftEOSdata*)data)->m,(( FF_SaftEOSdata*)data)->epsilon);
            break;
        case FF_SWtype:
            if (((( FF_SWEOSdata*)data)->eos==FF_IAPWS95)||((( FF_SWEOSdata*)data)->eos==IF97)){
                Tc=647.096;
                if (*T>=Tc)
                {
                    *Vp=1e12;
                }
                else
                {
                    double tau=*T-0.23855557567849/(*T-0.65017534844798e3);
                    double tau2=tau*tau;
                    double A=tau2+0.11670521452767e4*tau-0.72421316703206e6;
                    double B=-0.17073846940092e2*tau2+0.1202082470247e5*tau-0.32325550322333e7;
                    double C=0.1491510861353e2*tau2-0.48232657361591e4*tau+0.40511340542057e6;
                    *Vp=pow(2*C/(-B+pow(B*B-4*A*C,0.5)),4)*1e6;
                }
                return;
            }
            Tc=(( FF_SWEOSdata*)data)->Tc;
            Pc=(( FF_SWEOSdata*)data)->Pc;
            phiMaxError=0.001;
            break;
        //case FF_PRPOL1:
            //Tc=1000;
            //Pc=1e6;
            //break;
        default://Cubic eos
        {   enum FF_EOS tipo=(( FF_CubicEOSdata*)data)->eos;
            //printf("estoy en Vp, a ver que esos es:%i\n",tipo);
            if (((( FF_CubicEOSdata*)data)->eos==FF_PRPOL1)||((( FF_CubicEOSdata*)data)->eos==FF_SRKPOL2)){
                *Vp=0;
                return;
            }
            Tc=(( FF_CubicEOSdata*)data)->Tc;
            Pc=(( FF_CubicEOSdata*)data)->Pc;
            phiMaxError=0.001;
            //printf("Hola, soy Vp cubic eos:%i T:%f Tc:%f Pc:%f w:%f k1:%f\n",*eosType,*T,Tc,Pc,(( FF_CubicEOSdata*)data)->w,(( FF_CubicEOSdata*)data)->k1);
            break;
        }
    }
    if (((*T>=Tc)&&(!(Tc==0)))||((*T>0.999*Tc)&&(*eosType==FF_SW)))//If T> supplied Tc no calculation is made
    {
        *Vp=1e12;
        //printf("T:%f Tc:%f\n",*T,Tc);
    }
    else
    {
        double P,phiL,phiG;//Pressure and fugacity coef.
        double answerL[3],answerG[3];
        char option='b',state;
        if (Pc==0) Pc=200e5;
        else P=Pc/2;//Initial pressure guess
        //printf("Initial pressure guess:%f\n",P);
        FF_VfromTPeos(eosType,T,&P,data,&option,answerL,answerG,&state);
        //printf("T:%f P:%f  Vl:%f  Zl:%f  Vg:%f  Zg:%f\n",*T,P,answerL[0],answerL[2],answerG[0],answerG[2]);
        int i=4;//i will be used to help in the calculation of the new pressure

        while (!((answerL[0]>0) && (answerG[0]>0) && ((fabs(answerL[0]-answerG[0])/answerL[0])>0.002)))//Till we find a pressure with different possitive liquid and gas solutions
        {
            if (state=='g')P=P+Pc/i;
            else if (state=='l')P=P-Pc/i;
            if ((state=='b')&&(((fabs(answerL[0]-answerG[0])/answerL[0])<=0.002))){//If there is no specific gas phase we decrease pressure 0.35 antes
                if (answerL[2]<0.2) P=P-Pc/i;
                else if (answerL[2]>0.5) P=P+Pc/i;
                else//if we are not sure about if it is liquid or gas
                {
                    double Vprov=answerL[0]*0.7,ArrProv,Zprov;
                    FF_ArrZfromTVeos(eosType,T,&Vprov,data,&ArrProv,&Zprov);
                    //printf("Finding by Z n:%i i:%i  Vl:%f  Zl:%f  Vprov:%f  Zprov:%f\n",n,i/2,answerL[0],answerL[2],Vprov,Zprov);
                    if (Zprov<(answerL[2]-0.06)) P=P+Pc/i;
                    else if (Zprov>(answerL[2]-0.03)) P=P-Pc/i;

                    else {
                        *Vp=P;
                        return;
                    }
                }
            }

            FF_VfromTPeos(eosType,T,&P,data,&option,answerL,answerG,&state);
            i=i*2;
            n=n+1;
            //printf("Finding two phases n:%i i:%i P:%f  Vl:%f  Zl:%f  Vg%f  Zg%f\n",n,i/2,P,answerL[0],answerL[2],answerG[0],answerG[2]);
            if (n>15){
                *Vp=0;
                return;
            }
        }
        phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
        phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
        //printf("Initial phi calculation T:%f P:%f  phiL:%f  phiG:%f\n",*T,P,phiL,phiG);
        while ((phiL<0.3)){
            P=0.5*P;
            FF_VfromTPeos(eosType,T,&P,data,&option,answerL,answerG,&state);
            phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
            phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
            n=n+1;
            //printf("%i Pa:%f  phiL:%f  phiG:%f\n",n-1,P,phiL,phiG);
            if (n>50){
                *Vp=0;
                return;
            }
        }
        while (fabs(phiL-phiG)>phiMaxError)
        {
            P=P* pow((phiL / phiG),0.3);
            FF_VfromTPeos(eosType,T,&P,data,&option,answerL,answerG,&state);
            phiL=exp(answerL[1]+answerL[2]-1)/answerL[2];
            phiG=exp(answerG[1]+answerG[2]-1)/answerG[2];
            n=n+1;
            //printf("%i Pb:%f  phiL:%f  phiG:%f\n",n-1,P,phiL,phiG);
            if (n>100){
                *Vp=0;
                return;
            }
        }
        *Vp=P;
    }
}

//Thermodynamic properties calculation for a ideal gas at same T and V, from a reference state, specified by refT and refP, where H and S are 0
//---------------------------------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_IdealThermoEOS(const int *equation,const double coef[],double *refT,double *refP, FF_ThermoProperties *th0)
{   //the result delivered is always in J/mol·K
    //printf("Hi, I am individual FF_IdealThermoEOS\n");
    int i;
    th0->P=R*th0->T/th0->V;//this is the pressure of a ideal gas with the same volume
    //printf("T:%f P:%f refT:%f refP:%f\n",th0->T,th0->P,*refT,*refP);
    switch (*equation)
    {
        case 1://DIPPR 100 in KJ/kgr·K
        case 2://DIPPR 100 in J/mol·K
        case 6://DIPPR 100 in J/kgr·K
            th0->Cp=coef[0]+coef[1]*th0->T+coef[2]*pow(th0->T,2)+coef[3]*pow(th0->T,3)+coef[4]*pow(th0->T,4);
            th0->H=coef[0]*(th0->T-*refT)+coef[1]*(th0->T*th0->T-*refT* *refT)/2+coef[2]*(th0->T*th0->T*th0->T-*refT* *refT* *refT)/3+
                    coef[3]*(th0->T*th0->T*th0->T*th0->T-*refT* *refT* *refT* *refT)/4+
                    coef[4]*(th0->T*th0->T*th0->T*th0->T*th0->T-*refT* *refT* *refT* *refT* *refT)/5;//This is the integration from reference T to actual T at constant pressure
                    //as the derivative (dH/dP) at constant T is 0 for an ideal gas, this is all needed
            th0->S=coef[0]*log(th0->T/ *refT)+coef[1]*(th0->T- *refT)+coef[2]*(th0->T*th0->T- *refT* *refT)/2+coef[3]*(th0->T*th0->T*th0->T-
                    *refT* *refT* *refT)/3+coef[4]*(th0->T*th0->T*th0->T*th0->T-*refT* *refT* *refT* *refT)/4;
            break;
        case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 in cal/(mol·K)
            th0->Cp=coef[0]+coef[1]*th0->T+coef[2]*pow(th0->T,2)+coef[3]*pow(th0->T,3)+coef[4]*pow(th0->T,4)+coef[5]*pow(th0->T,5);
            th0->H=coef[0]*(th0->T-*refT)+coef[1]*(th0->T*th0->T-*refT* *refT)/2+coef[2]*(th0->T*th0->T*th0->T-*refT* *refT* *refT)/3+
                coef[3]*(th0->T*th0->T*th0->T*th0->T-*refT* *refT* *refT* *refT)/4+
                coef[4]*(th0->T*th0->T*th0->T*th0->T*th0->T-*refT* *refT* *refT* *refT* *refT)/5+
                coef[5]*(th0->T*th0->T*th0->T*th0->T*th0->T*th0->T-*refT* *refT* *refT* *refT* *refT* *refT)/6;//This is the integration from reference T to actual T at constant pressure
                //as the derivative (dH/dP) at constant T is 0 for an ideal gas, this is all needed
            th0->S=coef[0]*log(th0->T/ *refT)+coef[1]*(th0->T- *refT)+coef[2]*(th0->T*th0->T- *refT* *refT)/2+coef[3]*(th0->T*th0->T*th0->T-
                *refT* *refT* *refT)/3+coef[4]*(th0->T*th0->T*th0->T*th0->T-*refT* *refT* *refT* *refT)/4+
                coef[5]*(th0->T*th0->T*th0->T*th0->T*th0->T-*refT* *refT* *refT* *refT* *refT)/5;
            break;
        case 3://DIPPR 107 correlation in calories/mol·K
        case 4://DIPPR 107 correlation in J/Kmol·K
            th0->Cp=(coef[0]+coef[1]*pow(coef[2]/th0->T/sinh(coef[2]/th0->T),2)+coef[3]*pow(coef[4]/th0->T/cosh(coef[4]/th0->T),2));
            th0->H=(coef[0]*(th0->T-*refT)+coef[1]*coef[2]*(1/tanh(coef[2]/th0->T)-1/tanh(coef[2]/ *refT))+coef[3]*coef[4]*(tanh(coef[4]/ *refT)-tanh(coef[4]/th0->T)));
            th0->S=(coef[0]*log(th0->T/ *refT)+coef[1]*(coef[2]/th0->T/tanh(coef[2]/th0->T)-log(sinh(coef[2]/th0->T)))-coef[1]*(coef[2]/
                    *refT/tanh(coef[2]/ *refT)-log(sinh(coef[2]/ *refT)))-coef[3]*(coef[4]/th0->T*tanh(coef[4]/th0->T)-log(cosh(coef[4]/th0->T)))+
                    coef[3]*(coef[4]/ *refT*tanh(coef[4]/ *refT)-log(cosh(coef[4]/ *refT))));
            break;
        case 9:{//ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K. Integration is done numerically
            th0->Cp=coef[0]+exp(coef[1]/th0->T+coef[2]+coef[3]*th0->T+coef[4]*pow(th0->T,2));
            th0->H=coef[0]*(th0->T- *refT);
            th0->S=coef[0]*log(th0->T/ *refT);
            int i=20;
            double interval=(th0->T- *refT)/i;
            double T,Cp1,Cp2;
            Cp1=exp(coef[1]/ *refT+coef[2]+coef[3]* *refT+coef[4]*pow(*refT,2));
            for (i=1;i<21;i++){
                T=*refT+i*interval;
                Cp2=exp(coef[1]/T+coef[2]+coef[3]*T+coef[4]*pow(T,2));
                th0->H=th0->H+(Cp1+Cp2)/2*interval;
                th0->S=th0->S+(Cp1+Cp2)/(T+T-interval)*interval;
                Cp1=Cp2;
                }
            }
            break;
        case 5:{//Wilhoit equation J/mol·K (8 coefficients)
            double y,y2,y4,h,a7_a6,a7_a6_2,a7_a6_4,x1,z,w,s;
            a7_a6=coef[7]/coef[6];
            a7_a6_2=a7_a6*a7_a6;
            a7_a6_4=a7_a6_2*a7_a6_2;
            x1=(coef[4]*coef[7]*coef[7] - coef[5])/(coef[6]*coef[6]);

            if (th0->T<=coef[7]) y=0;
            else y=(th0->T-coef[7])/(th0->T+coef[6]);
            y2=y*y;
            y4=y2*y2;
            th0->Cp=R*(coef[0] + coef[1]/pow(th0->T,2)*exp(-coef[2]/th0->T) + coef[3]*y2 + (coef[4] - coef[5] /pow((th0->T -coef[7]),2))*y4*y4);

            if (th0->T<=coef[7]) h=0;
            else h=(coef[6]+coef[7])*((2*coef[3]+8*coef[4])*log(1-y)+ (coef[3]*(1+1/(1-y))+coef[4]*(7+1/(1-y)))*y+
                    coef[4]*(3*y2+5*y*y2/3+y4+0.6*y4*y+y4*y2/3)+ (coef[4]-coef[5]/pow((coef[6]+coef[7]),2))*y4*y2*y/7);
            th0->H= R*th0->T*(coef[0]+coef[1]*exp(-coef[2]/th0->T)/(coef[2]*th0->T))+R*h;

            if (th0->T<=coef[7]) s=0;
            else{
                z = th0->T*(coef[7] + coef[6])/((th0->T + coef[6])*coef[7]);
                w=0;
                for (i=1;i<8;i++) w=w+(x1*pow(-a7_a6,6-i) - coef[4])*pow(y,i)/i;
                s=(coef[3] + ((coef[4]*coef[7]*coef[7]-coef[5])*a7_a6_4/(coef[6]*coef[6])))*a7_a6_2*log(z)+
                    (coef[3] + coef[4])*log((th0->T + coef[6])/(coef[6] + coef[7]))-
                    (coef[3]*(coef[6] + coef[7])/coef[6] + coef[5]*y4*y2/(7.*coef[7]*(coef[6] + coef[7])))*y+w;
            }
            th0->S=R*(coef[0]*log(th0->T)+coef[1]*(1+coef[2]/th0->T)*exp(-coef[2]/th0->T)/(coef[2]*coef[2])+s);

            if (*refT<coef[7]) y=0;
            else y=(*refT-coef[7])/(*refT+coef[6]);
            y2=y*y;
            y4=y2*y2;
            if (*refT<=coef[7]) h=0;
            else h=(coef[6]+coef[7])*((2*coef[3]+8*coef[4])*log(1-y)+ (coef[3]*(1+1/(1-y))+coef[4]*(7+1/(1-y)))*y+
                    coef[4]*(3*y2+5*y*y2/3+y4+0.6*y4*y+y4*y2/3)+ (coef[4]-coef[5]/pow((coef[6]+coef[7]),2))*y4*y2*y/7);
            th0->H=th0->H-R* *refT*(coef[0]+coef[1]*exp(-coef[2]/ *refT)/(coef[2]*th0->T))-R*h;

            if (*refT<=coef[7]) s=0;
            else{
                z = *refT*(coef[7] + coef[6])/((*refT + coef[6])*coef[7]);
                w=0;
                for (i=1;i<8;i++) w=w+(x1*pow(-a7_a6,6-i) - coef[4])*pow(y,i)/i;
                s=(coef[3] + ((coef[4]*coef[7]*coef[7]-coef[5])*a7_a6_4/(coef[6]*coef[6])))*a7_a6_2*log(z)+
                        (coef[3] + coef[4])*log((*refT + coef[6])/(coef[6] + coef[7]))-
                        (coef[3]*(coef[6] + coef[7])/coef[6] + coef[5]*y4*y2/(7.*coef[7]*(coef[6] + coef[7])))*y+w;
            }
            th0->S=th0->S-R*(coef[0]*log(*refT)+coef[1]*(1+coef[2]/ *refT)*exp(-coef[2]/ *refT)/(coef[2]*coef[2])+s);
            }
            break;
        case 7://Cooper (11 coefficients used in IAPWS95 and CO2) plus potential term  (used in short fundamental equations with 11 coefficients also,lacks last exp terms)
            th0->Cp=coef[0]+coef[1]*pow(th0->T,coef[2]);
            for (i=3;i<13;i=i+2){
                if (coef[i]>0) th0->Cp=th0->Cp+coef[i]*pow((coef[i+1]/th0->T),2)*exp(coef[i+1]/th0->T)/pow((exp(coef[i+1]/th0->T)-1),2);
            }
            th0->Cp=th0->Cp*R;
            th0->H=coef[0]*(th0->T- *refT)+coef[1]/(coef[2]+1)+(pow(th0->T,(coef[2]+1))-pow(*refT,(coef[2]+1)));
            for (i=3;i<13;i=i+2) if (coef[i]>0) th0->H=th0->H+coef[i]*coef[i+1]*(1/(exp(coef[i+1]/th0->T)-1)-1/(exp(coef[i+1]/ *refT)-1));
            th0->H=th0->H*R;
            if (coef[1]>0) th0->S=coef[0]*log(th0->T/ *refT)+coef[1]/coef[2]*(pow(th0->T,coef[2])-pow(*refT,coef[2]));
            else th0->S=coef[0]*log(th0->T/ *refT);
            //printf("T:%f S0(0):%f\n",th0->T,th0->S*R/th0->MW);
            for (i=3;i<13;i=i+2){
                if (coef[i]>0) th0->S=th0->S+coef[i]*coef[i+1]*(exp(coef[i+1]/th0->T)/(th0->T*(exp(coef[i+1]/th0->T)-1))-
                    log(exp(coef[i+1]/th0->T)-1)/coef[i+1]-exp(coef[i+1]/ *refT)/(*refT*(exp(coef[i+1]/ *refT)-1))+log(exp(coef[i+1]/ *refT)-1)/coef[i+1]);
                //printf("S0(%i):%f\n",i,th0->S*R/th0->MW);
            }
            th0->S=th0->S*R;
            break;
        case 8://Jaeschke and Schley equation (9 coefficients). Used by GERG2004
            th0->Cp=1+coef[0];
            for (i=1;i<9;i=i+4) if (coef[i]>0) th0->Cp=th0->Cp+coef[i]*pow(coef[i+1]/th0->T/sinh(coef[i+1]/th0->T),2)+coef[i+2]*pow(coef[i+3]/th0->T/cosh(coef[i+3]/th0->T),2);
            th0->Cp=th0->Cp*R;
            th0->H=(1+coef[0])*(th0->T- *refT);
            for (i=1;i<9;i=i+4) if (coef[i]>0) th0->H=th0->H+2*coef[i]*coef[i+1]*(1/(exp(2*coef[i+1]/th0->T)-1)-1/(exp(2*coef[i+1]/ *refT)-1))+
                    +2*coef[i+2]*coef[i+3]*(1/(exp(2*coef[i+3]/th0->T)+1)-1/(exp(2*coef[i+3]/ *refT)+1));
            th0->H=th0->H*R;
            th0->S=(1+coef[0])*log(th0->T/ *refT);
            for (i=1;i<9;i=i+4) if (coef[i]>0) th0->S=th0->S+coef[i]*(coef[i+1]/th0->T/tanh(coef[i+1]/th0->T)-log(sinh(coef[i+1]/th0->T)))-coef[i]*(coef[i+1]/
                    *refT/tanh(coef[i+1]/ *refT)-log(sinh(coef[i+1]/ *refT)))-coef[i+2]*(coef[i+3]/th0->T*tanh(coef[i+3]/th0->T)-log(cosh(coef[i+3]/th0->T)))+
                    coef[i+2]*(coef[i+3]/ *refT*tanh(coef[i+3]/ *refT)-log(cosh(coef[i+3]/ *refT)));
            th0->S=th0->S*R;
            break;
        default:
        th0->Cp=th0->H=th0->S=0;
    }
    switch (*equation)//units conversion if necessary
    {
    case 1://DIPPR 100 in KJ/kgr·K
        th0->Cp=th0->Cp*th0->MW;
        th0->H=th0->H*th0->MW;
        th0->S=th0->S*th0->MW;
        break;
    case 6://DIPPR 100 in J/kgr·K
        th0->Cp=th0->Cp*th0->MW*1e-3;
        th0->H=th0->H*th0->MW*1e-3;
        th0->S=th0->S*th0->MW*1e-3;
        break;
    case 3://DIPPR 107 in calories/mol·K
    case 10://Polynomial a+b*T+c*T^2+d*T^3+e*T^4+f*T^5 in cal/(mol·K)
        th0->Cp=th0->Cp*4.1868;
        th0->H=th0->H*4.1868;
        th0->S=th0->S*4.1868;
        break;
    case 4://DIPPR 107 correlation in J/Kmol·K
    case 9://ChemSep16 a + exp( b/T+ c + d*T + e*T^2 ) en J/Kmol·K
        th0->Cp=th0->Cp*1e-3;
        th0->H=th0->H*1e-3;
        th0->S=th0->S*1e-3;
        break;
    }

    if (th0->Cp>0){
        th0->S=th0->S-R*log(th0->P/ *refP);//Now we change the reference of entropy from the same P to the reference P
        th0->Cv=th0->Cp-R;
        th0->U=th0->H-R*th0->T;
        th0->A=th0->U-th0->T*th0->S;
        th0->G=th0->H-th0->T*th0->S;
    }
    else{
        th0->Cv=th0->U=th0->A=th0->G=0;
    }


    //printf("T:%f MW:%f Cp0:%f  G0:%f  H0:%f  S0:%f\n",th0->T,th0->MW,th0->Cp,th0->G,th0->H,th0->S);
}


//Ideal gas thermodynamic properties of water calculation , from a reference state specified by the triple point where H and S are 0
//----------------------------------------------------------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_IdealThermoWater( FF_ThermoProperties *th0)
{
    double delta=1/(th0->V*1.787371609e4);
    double tau=647.096/th0->T;
    double A0r,dt,dt2;//A0r and derivatives
    //double dd,dd2;
    //double n[8];//={–8.3204464837497,6.6832105275932,3.00632,0.012436,0.97315,1.27950,0.96956,0.24873};

    A0r=(log(delta)-8.3204464837497+6.6832105275932*tau+3.00632*log(tau)+0.012436*log(1-exp(-1.28728967*tau))+0.97315*log(1-exp(-3.53734222*tau))+
            1.27950*log(1-exp(-7.74073708*tau))+0.96956*log(1-exp(-9.24437796*tau))+0.24873*log(1-exp(-27.5075105*tau)));
    //dd=1/delta;
    //dd2=-1/pow(delta,2);
    dt=6.6832105275932+3.00632/tau+0.012436*1.28728967*(pow((1-exp(-1.28728967*tau)),-1)-1)+0.97315*3.53734222*(pow((1-exp(-3.53734222*tau)),-1)-1)+
            1.27950*7.74073708*(pow((1-exp(-7.74073708*tau)),-1)-1)+0.96956*9.24437796*(pow((1-exp(-9.24437796*tau)),-1)-1)+
            0.24873*27.5075105*(pow((1-exp(-27.5075105*tau)),-1)-1);
    dt2=-3.00632/pow(tau,2)-0.012436*pow(1.28728967,2)*exp(-1.28728967*tau)*pow((1-exp(-1.28728967*tau)),-2)-
            0.97315*pow(3.53734222,2)*exp(-3.53734222*tau)*pow((1-exp(-3.53734222*tau)),-2)-
            1.27950*pow(7.74073708,2)*exp(-7.74073708*tau)*pow((1-exp(-7.74073708*tau)),-2)-
            0.96956*pow(9.24437796,2)*exp(-9.24437796*tau)*pow((1-exp(-9.24437796*tau)),-2)-
            0.24873*pow(27.5075105,2)*exp(-27.5075105*tau)*pow((1-exp(-27.5075105*tau)),-2);
    th0->A=R*th0->T*A0r;
    th0->S=R*(tau*dt-A0r);
    th0->U=th0->A+th0->S*th0->T;
    th0->G=th0->A+R*th0->T;
    th0->H=th0->U+R*th0->T;
    th0->Cv=-R*pow(tau,2)*dt2;
    th0->Cp=th0->Cv+R;
}

//Residual thermodynamic properties calculation from T and V, using EOS
//---------------------------------------------------------------------
void CALLCONV FF_ResidualThermoEOS(const enum FF_EosType *eosType,const void *data, FF_ThermoProperties *thR)
{
    double ArrDer[6];
     FF_CubicParam param;
    double delta,tau;
    switch (*eosType)
    {
        case FF_SAFTtype:
            FF_ArrDerPCSAFT(&thR->T,&thR->V,data,ArrDer);
            break;
        case FF_SWtype:
            delta=1/(thR->V*((( FF_SWEOSdata*)data)->rhoRef));
            tau=((( FF_SWEOSdata*)data)->tRef)/thR->T;
            //printf("delta:%f  tau:%f\n",delta,tau);
            FF_ArrDerSW(&tau,&delta,data,ArrDer);
            break;
        default:
            FF_FixedParamCubic(data, &param);
            FF_ThetaDerivCubic(&thR->T,data, &param);
            FF_ArrDerCubic(&thR->T,&thR->V,&param,ArrDer);
            break;
    }
    if (*eosType==FF_SWtype)
    {
        thR->P=R*thR->T*(1+delta*ArrDer[1])/thR->V;
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=R*(tau*ArrDer[3]-ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        thR->H=thR->G+thR->T*thR->S;
    }
    else
    {
        thR->P=R*thR->T*(1/thR->V-ArrDer[1]);
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=-R*(thR->T*ArrDer[3]+ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        //thR->U=-R*thR->T*thR->T*ArrDer[3];
        thR->H=thR->G+thR->T*thR->S;
        //thR->H=thR->U-R*thR->T*thR->V*ArrDer[1];
    }
}

//Residual extended thermodynamic properties calculation from T and V, using EOS
//------------------------------------------------------------------------------
EXP_IMP void CALLCONV FF_ExtResidualThermoEOS(const enum FF_EosType *eosType,const void *data, FF_ThermoProperties *thR)
{
    double ArrDer[6];
     FF_CubicParam param;
    double delta,tau;
    switch (*eosType)
    {
        case FF_SAFTtype:
            FF_ArrDerPCSAFT(&thR->T,&thR->V,data,ArrDer);
            break;
        case FF_SWtype:
            delta=1/(thR->V*((( FF_SWEOSdata*)data)->rhoRef));
            tau=((( FF_SWEOSdata*)data)->tRef)/thR->T;
            //printf("delta:%f  tau:%f\n",delta,tau);
            FF_ArrDerSW(&tau,&delta,data,ArrDer);
            break;
        default:
            FF_FixedParamCubic(data, &param);
            FF_ThetaDerivCubic(&thR->T,data, &param);
            FF_ArrDerCubic(&thR->T,&thR->V,&param,ArrDer);
            break;
    }
    //printf("Arr:%f dArr/dV:%f d2Arr/dV2:%f dArr/dT:%f d2Arr/dT2:%f d2Arr/dVdT:%f\n",ArrDer[0],ArrDer[1],ArrDer[2],ArrDer[3],ArrDer[4],ArrDer[5]);
    if (*eosType==FF_SWtype)
    {
        thR->P=R*thR->T*(1+delta*ArrDer[1])/thR->V;
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=R*(tau*ArrDer[3]-ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        thR->H=thR->G+thR->T*thR->S;
        thR->dP_dV=-R*thR->T*(1+2*delta*ArrDer[1]+delta*delta*ArrDer[2])/thR->V/thR->V;
        thR->dP_dT=R*(1+delta*ArrDer[1]-delta*tau*ArrDer[5])/thR->V;
        thR->Cv=-R*tau*tau*ArrDer[4];
        thR->Cp=thR->Cv+R*pow((1+delta*ArrDer[1]-delta*tau*ArrDer[5]),2)/(1+2*delta*ArrDer[1]+delta*delta*ArrDer[2])-R;
        //thR->Cp=thR->Cv-thR->T*thR->dP_dT*thR->dP_dT/thR->dP_dV-R;
    }
    else
    {
        thR->P=R*thR->T*(1/thR->V-ArrDer[1]);
        thR->A=R*thR->T*ArrDer[0];
        thR->G=thR->A+thR->P*thR->V-R*thR->T;
        thR->S=-R*(thR->T*ArrDer[3]+ArrDer[0]);
        thR->U=thR->A+thR->T*thR->S;
        //thR->U=-R*thR->T*thR->T*ArrDer[3];
        //printf("%f %f\n",thR->U,-R*thR->T*thR->T*ArrDer[3]);
        thR->H=thR->G+thR->T*thR->S;
        //thR->H=thR->U-R*thR->T*thR->V*ArrDer[1];
        //thR->H=-R*thR->T*(thR->V*ArrDer[1]+thR->T*ArrDer[3]);
        thR->Cv=-2*R*thR->T*ArrDer[3]-R*thR->T*thR->T*ArrDer[4];
        thR->dP_dV=R*thR->T*(-1/thR->V/thR->V-ArrDer[2]);
        thR->dP_dT=R*(1/thR->V-ArrDer[1]-thR->T*ArrDer[5]);
        //thR->Cp=thR->Cv-thR->T*pow(R*(1/thR->V-ArrDer[1]-thR->T*ArrDer[5]),2)/(R*thR->T*(-1/thR->V/thR->V-ArrDer[2]))-R;
        thR->Cp=thR->Cv-thR->T*thR->dP_dT*thR->dP_dT/thR->dP_dV-R;

    }
}

//Thermodynamic properties calculation from T and V, from a reference state (specified by T and P) where H and S are 0
//--------------------------------------------------------------------------------------------------------------------
void CALLCONV FF_ThermoEOS(const enum FF_EosType *eosType,const void *data,const int *equation,const double coef[],double *refT,double *refP, FF_ThermoProperties *th)
{
     FF_ThermoProperties th0,thR;
    bool water=false;
    switch (*eosType)
    {
        case FF_SAFTtype:
            th->MW=(( FF_SaftEOSdata *)data)->MW;
            break;
        case FF_SWtype:
            th->MW=(( FF_SWEOSdata *)data)->MW;
            if ((( FF_SWEOSdata*)data)->eos==FF_IAPWS95)water=true;
            break;
        default://Cubic EOS
            th->MW=(( FF_CubicEOSdata *)data)->MW;
            break;
    }
    th0.MW=thR.MW=th->MW;//Perhaps not necessary
    th0.T=thR.T=th->T;
    th0.V=thR.V=th->V;
    if (water==true){
        FF_IdealThermoWater(&th0);
    }
    else{
        FF_IdealThermoEOS(equation,coef,refT,refP,&th0);
    }
    if (*eosType==FF_IdealType)
    {
        th->P=R*th->T/th->V;
        th->A=th0.A;
        th->G=th0.G;
        th->S=th0.S;
        th->U=th0.U;
        th->H=th0.H;
        th->dP_dV=-R*th->T/th->V/th->V;
        th->dP_dT=R/th->V;
        th->Cv=th0.Cv;
        th->Cp=th0.Cp;
    }
    else
    {
        FF_ExtResidualThermoEOS(eosType,data,&thR);
        th->P=thR.P;
        th->A=th0.A+thR.A;
        th->G=th0.G+thR.G;
        th->S=th0.S+thR.S;
        th->U=th0.U+thR.U;
        th->H=th0.H+thR.H;
        th->dP_dV=thR.dP_dV;
        th->dP_dT=thR.dP_dT;
        th->Cv=th0.Cv+thR.Cv;
        th->Cp=th0.Cp+thR.Cp;
    }
    th->SS=th->V*pow(-th->Cp*th->dP_dV/th->Cv/th->MW*1e3,0.5);
    th->JT=-(th->T*th->dP_dT/th->dP_dV+th->V)/th->Cp;//Joule Thomson coefficient =(dT/dP) at constant H
    th->IT=-th->JT*th->Cp;//Isothermal throttling coefficient = (dH/dP) at constant T
    //printf("T:%f V:%f P:%f Cv:%f Cp:%f H:%f U:%f S:%f G:%f A:%f dP_dV:%f dP_dT:%f SS:%f JT:%f IT:%f\n",th->T,th->V,th->P,th->Cv,th->Cp,th->H,th->U,th->S,th->G,th->A,th->dP_dV,th->dP_dT,th->SS,th->JT,th->IT);
    //printf("Ideal T:%f V:%f P:%f Cv0:%f Cp0:%f H0:%f U0:%f S0:%f G:%f A:%f\n",th0.T,th0.V,th0.P,th0.Cv,th0.Cp,th0.H,th0.U,th0.S,th0.G,th0.A);
    //printf("Residual T:%f V:%f P:%f Cvr:%f Cpr:%f Hr:%f Ur:%f Sr:%f Gr:%f Ar:%f\n",thR.T,thR.V,thR.P,thR.Cv,thR.Cp,thR.H,thR.U,thR.S,thR.G,thR.A);
}



//Calculation of thermo properties and liquid/gas fraction from P and (H or U or S)
//---------------------------------------------------------------------------------
void CALLCONV FF_ThermoEOSfromPX(const enum FF_EosType *eosType,const void *data,const int *equation,const double coef[],double *refT,double *refP,char *var, FF_ThermoProperties *th,double *liqFraction)
{
    //First we find the boiling point
    double Tb;
    FF_TbEOS(eosType,&th->P,data,&Tb);
    //Now we calculate the thermo properties of the liquid and gas phases at boiling point
    char option='b';
    char state;
    double answerL[3],answerG[3];
     FF_ThermoProperties thG,thL;
    FF_VfromTPeos(eosType,&Tb,&th->P,data,&option,answerL,answerG,&state);//we calculate gas and liquid volumes at boiling point
    thG.T=thL.T=Tb;
    thG.V=answerG[0];
    thL.V=answerL[0];
    //printf("Tb:%f Vg:%f Vl:%f\n",Tb,thG.V,thL.V);
    FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thG);//and now the gas and liquid thermo properties at boiling point
    FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);
    th->MW=thL.MW=thG.MW;
    int i=1;//only for printing
    if ((*var=='H')||(*var=='h')){
        if ((th->H > thL.H)&&(th->H < thG.H)) *liqFraction=(th->H-thG.H)/(thL.H-thG.H);
        else if (th->H<=thL.H)
        {
            option='l';
            while (fabs((th->H-thL.H)/th->H)>0.00001)
            {
                thL.T=thL.T+(th->H-thL.H)/thL.Cp;
                FF_VfromTPeos(eosType,&thL.T,&thL.P,data,&option,answerL,answerG,&state);//we need to recalculate the volume
                thL.V=answerL[0];
                FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);
                //printf("i:%i T:%f V:%f H:%f\n",i,thL.T,thL.V,thL.H);
                i=i+1;
            }
            *liqFraction=1;
        }
        else if (th->H>=thG.H)
        {
            //printf("Hola, miro el gas");
            option='g';
            while (fabs((th->H-thG.H)/th->H)>0.00001)
            {
                thG.T=thG.T+(th->H-thG.H)/thG.Cp;
                FF_VfromTPeos(eosType,&thG.T,&thG.P,data,&option,answerL,answerG,&state);//we need to recalculate the volume
                thG.V=answerG[0];
                FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thG);
                //printf("i:%i T:%f V:%f Hbuscada:%f Hencontrada:%f\n",i,thG.T,thG.V,th->H,thG.H);
                i=i+1;
            }
            *liqFraction=0;
        }
    }
    if ((*var=='U')||(*var=='u')){
        if ((th->U > thL.U)&&(th->U < thG.U)) *liqFraction=(th->U-thG.U)/(thL.U-thG.U);
        else if (th->U<=thL.U)
        {
            option='l';
            while (fabs((th->U-thL.U)/th->U)>0.00001)
            {
                thL.T=thL.T+(th->U-thL.U)/thL.Cv;
                FF_VfromTPeos(eosType,&thL.T,&thL.P,data,&option,answerL,answerG,&state);//we need to recalculate the volume
                thL.V=answerL[0];
                FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);
                //printf("i:%i T:%f V:%f H:%f\n",i,thL.T,thL.V,thL.U);
                i=i+1;
            }
            *liqFraction=1;
        }
        else if (th->U>=thG.U)
        {
            option='g';
            while (fabs((th->U-thG.U)/th->U)>0.00001)
            {
                thG.T=thG.T+(th->U-thG.U)/thG.Cv;
                FF_VfromTPeos(eosType,&thG.T,&thG.P,data,&option,answerL,answerG,&state);//we need to recalculate the volume
                thG.V=answerG[0];
                FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thG);
                //printf("i:%i T:%f V:%f H:%f\n",i,thG.T,thG.V,thG.U);
                i=i+1;
            }
            *liqFraction=0;
        }
    }
    if ((*var=='S')||(*var=='s')){
        if ((th->S > thL.S)&&(th->S < thG.S)) *liqFraction=(th->S-thG.S)/(thL.S-thG.S);
        else if (th->S<=thL.S)
        {
            option='l';
            while (fabs((th->S-thL.S)/th->S)>0.00001)
            {
                thL.T=thL.T+(th->S-thL.S)/thL.Cp*thL.T;
                FF_VfromTPeos(eosType,&thL.T,&thL.P,data,&option,answerL,answerG,&state);//we need to recalculate the volume
                thL.V=answerL[0];
                FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);
                //printf("i:%i T:%f V:%f H:%f\n",i,thL.T,thL.V,thL.S);
                i=i+1;
            }
            *liqFraction=1;
        }
        else if (th->S>=thG.S)
        {
            option='g';
            while (fabs((th->S-thG.S)/th->S)>0.00001)
            {
                thG.T=thG.T+(th->S-thG.S)/thG.Cp*thG.T;
                FF_VfromTPeos(eosType,&thG.T,&thG.P,data,&option,answerL,answerG,&state);//we need to recalculate the volume
                thG.V=answerG[0];
                FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thG);
                //printf("i:%i T:%f V:%f H:%f\n",i,thG.T,thG.V,thG.S);
                i=i+1;
            }
            *liqFraction=0;
        }
    }
    if (*liqFraction==0) *th=thG;
    else if (*liqFraction==1) *th=thL;
    else {
        th->T=thL.T;
        th->P=thL.P;
        th->V=thL.V* *liqFraction+thG.V*(1-*liqFraction);
        th->A=thL.A* *liqFraction+thG.A*(1-*liqFraction);
        th->G=thL.G* *liqFraction+thG.G*(1-*liqFraction);
        th->S=thL.S* *liqFraction+thG.S*(1-*liqFraction);
        th->U=thL.U* *liqFraction+thG.U*(1-*liqFraction);
        th->H=thL.H* *liqFraction+thG.H*(1-*liqFraction);
        th->Cp=thL.Cp* *liqFraction+thG.Cp*(1-*liqFraction);
        th->Cv=thL.Cv* *liqFraction+thG.Cv*(1-*liqFraction);
        th->dP_dT=0;
        th->dP_dV=0;
        th->SS=0;
        th->JT=0;
        th->IT=0;
    }
}


//Calculation of thermo properties and liquid/gas fraction from P and (H or U or S)
//---------------------------------------------------------------------------------
void CALLCONV FF_ThermoEOSfromVX(const enum FF_EosType *eosType,const void *data,const int *equation,const double coef[],double *refT,double *refP,char *var, FF_ThermoProperties *th,double *liqFraction)
{
    char option='s';
    char state;
    double answerL[3],answerG[3];
     FF_ThermoProperties thG,thL;
    double Vp;
    double Tc;
    double Hmix,Cmix;
    int i;
    switch (*eosType)
    {
        case FF_SAFTtype:
            Tc=(( FF_SaftEOSdata*)data)->Tc;
            break;
        case FF_SWtype:
            Tc=(( FF_SWEOSdata*)data)->Tc;
            break;
        default://Cubic eos
        {
            Tc=(( FF_CubicEOSdata*)data)->Tc;
            break;
        }
    }

    *liqFraction=0.5;
    thL.V=th->V;
    if ((*var=='T')||(*var=='t')){
    thL.T=th->T;
    }
    else if ((*var=='H')||(*var=='h')){
        thL.T=273.15;//we begin at 273.15
        FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);//we use thL to calculate H
        //th->MW=thG.MW=thL.MW;//we recover the molecular weight, that has been calculated
        //printf("MW:%f V:%f H:%f\n",th->MW,th->V*1e6,th->H);
        //printf("V:%f T:%f P:%f H:%f\n",thL.V*1e6,thL.T,thL.P,thL.H);
        i=1;//only for pinting and maximum loop iteration
        while (fabs((th->H-thL.H)/th->H)>0.000001)//we go to a temperature with the correct H
        {
            thL.T=thL.T+(th->H-thL.H)/thL.Cp*0.5;
            FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);
            //printf("i:%i V:%f T:%f P:%f H:%f Hobjetivo:%f\n",i,thL.V*1e6,thL.T,thL.P,thL.H,th->H);
            i=i+1;
            if (i>100) break;
        }
    }
    //Now that we have a candidate T and V, first we check that this gives a possitive pressure
    FF_PfromTVeos(eosType,&thL.T,&thL.V,data,&thL.P);//we use thL for recovering the pressure
    if (thL.P>0){
        FF_VfromTPeos(eosType,&thL.T,&thL.P,data,&option,answerL,answerG,&state);//if P>0 we check the possible volumes
        if (((state=='L')||(state=='U'))&&(fabs((th->V-answerL[0])/th->V)<0.0001)){
            *liqFraction=1;
            FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);
        }
        else if ((state=='G')&&(fabs((th->V-answerG[0])/th->V)<0.0001)){
            *liqFraction=0;
            thG.T=thL.T;
            thG.V=thL.V;
            FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thG);
        }
    }
    if ((*liqFraction>0)&&(*liqFraction<1)){// if equilibrium composition is needed
        //printf("Hola liqFraction var:%c\n",*var);
        if ((*var=='H')||(*var=='h')){//here is not clear that T is good, it is necessary a loop
            //printf("Hola H\n");
            thL.T=thG.T=0.7*Tc;//We begin at 0.7Tc
            Hmix=0;
            option='b';
            i=1;
            do {
                FF_VpEOS(eosType,&thL.T,data,&Vp);
                FF_VfromTPeos(eosType,&thL.T,&Vp,data,&option,answerL,answerG,&state);
                thL.V=answerL[0];
                thG.V=answerG[0];
                FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);
                FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thG);
                *liqFraction=(th->V-thG.V)/(thL.V-thG.V);
                Hmix=*liqFraction*thL.H+(1-*liqFraction)*thG.H;
                Cmix=*liqFraction*thL.Cp+(1-*liqFraction)*thG.Cp;
                thL.T=thG.T=thL.T+(th->H-Hmix)/Cmix*0.1;
                printf("i:%i T(C):%f Hmix:%f Hobjetivo:%f\n",i,thL.T-273.15,Hmix,th->H);
                i++;
                if ((i>20)||(thL.T>=(0.965*Tc))) return;
            }while (fabs((th->H-Hmix)/th->H)>0.000001);
        }
        //printf("Hola antes de Vp T:%f\n",thL.T-273.15);
        FF_VpEOS(eosType,&thL.T,data,&Vp);//as we know that the temperature is good, we can calculate the pressure for 2 phases
        if (Vp>0){
            option='b';
            //printf("Hola Vp T:%f Vp:%f\n",th->T-273.15,Vp*1e-5);
            FF_VfromTPeos(eosType,&thL.T,&Vp,data,&option,answerL,answerG,&state);
            thL.V=answerL[0];
            thG.V=answerG[0];
            thG.T=thL.T;
            //printf("state:%c Vl:%f Vg:%f",state,thL.V,thG.V);
            FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thL);
            FF_ThermoEOS(eosType,data,equation,coef,refT,refP,&thG);
            *liqFraction=(th->V-thG.V)/(thL.V-thG.V);
        }

    }
    th->T=thL.T;
    th->P=thL.P;
    th->A=thL.A* *liqFraction+thG.A*(1-*liqFraction);
    th->G=thL.G* *liqFraction+thG.G*(1-*liqFraction);
    th->S=thL.S* *liqFraction+thG.S*(1-*liqFraction);
    th->U=thL.U* *liqFraction+thG.U*(1-*liqFraction);
    th->H=thL.H* *liqFraction+thG.H*(1-*liqFraction);
    th->Cp=thL.Cp* *liqFraction+thG.Cp*(1-*liqFraction);
    th->Cv=thL.Cv* *liqFraction+thG.Cv*(1-*liqFraction);
    th->dP_dT=0;
    th->dP_dV=0;
    th->SS=0;
    th->JT=0;
    th->IT=0;
    //printf("P:%f\n",th->P);
}
/*
void CALLCONV TVfromVHeos(const enum FF_EosType *eosType,const void *data,const  FF_Correlation *cp,double *refT,double *refP, FF_ThermoProperties *th,double *liqFraction)
{
    char option='s';
    char state;
    double answerL[3],answerG[3];
     FF_ThermoProperties thG,thL;
    double Vp;
    double Hmix;
    double Vmix;
    //we begin at 273.15 K and calculate the enthalpy
    thL.T=273.15;
    thL.V=thG.V=th->V;
    FF_ThermoEOS(eos,data,cp,refT,refP,&thL);
    th->MW=thG.MW=thL.MW;//we recover the molecular weight, that has been calculated
    printf("MW:%f V:%f H:%f\n",th->MW,th->V*1e6,th->H);
    printf("V:%f T:%f P:%f H:%f\n",thL.V*1e6,thL.T,thL.P,thL.H);
    int i=1;//only for pinting
    while (fabs((th->H-thL.H)/th->H)>0.000001)//we go to a temperature with the correct H
    {
        thL.T=thL.T+(th->H-thL.H)/thL.Cp*0.5;
        FF_ThermoEOS(eos,data,cp,refT,refP,&thL);
        printf("i:%i V:%f T:%f P:%f H:%f\n",i,thL.V*1e6,thL.T,thL.P,thL.H);
        i=i+1;
    }
    //Once we have a solution, we need to know if it is a stable phase
    FF_VfromTPeos(eos,&thL.T,&thL.P,data,&option,answerL,answerG,&state);//we obtain the volumes, Arr, and Z for the possible phases
    if (((state=='L') || (state=="U")) && (fabs((th->V-answerL[0])/th->V)<0.0001) && (thL.P>0)){
        *liqFraction=1;
        th->T=thL.T;
        printf("H:%f T:%f,P:%f\n",thL.H,thL.T,thL.P);
    }
    else if ((state=='G') && (fabs((th->V-answerG[0])/th->V)<0.0001) && (thL.P>0) ){
        *liqFraction=0;
        th->T=thL.T;
        printf("H:%f T:%f,P:%f\n",thL.H,thL.T,thL.P);
    }
    else{//if no solution with only liquid or gas phase
        thL.T=thG.T=403.15;
        FF_VpEOS(eos,&thL.T,data,&Vp);
        FF_VfromTPeos(eos,&thL.T,&Vp,data,&option,answerL,answerG,&state);
        thL.V=answerL[0];
        thG.V=answerG[0];
        FF_ThermoEOS(eos,data,cp,refT,refP,&thL);
        FF_ThermoEOS(eos,data,cp,refT,refP,&thG);
        *liqFraction=(th->V-thG.V)/(thL.V-thG.V);
        Vmix=*liqFraction*thL.V+(1-*liqFraction)*thG.V;
        Hmix=*liqFraction*thL.H+(1-*liqFraction)*thG.H;
        printf("T:%f P:%f Hl:%f Hg:%f Vl:%f Vg:%f LiqFraction:%f Vmix:%f Hmix:%f\n",thL.T,thL.P,thL.H,thG.H,thL.V*1e6,thG.V*1e6,*liqFraction,Vmix*1e6,Hmix);
    }


}*/
