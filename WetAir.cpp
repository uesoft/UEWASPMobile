/*library UEWASP; 
注：
     本DLL 根据IAPWS-IF97 和IFC67
     (1985IAPS热力学性质国际骨架表)、湿空气
	 性质(《空气调节》1994年11月第三版，
	 中国建筑工业出版社)编写。
     本DLL中的各个参数的单位均是国际单位，
     其中压力的单位为MPa，温度的单位为℃，
	 焓的单位kJ/kg，熵的单位kJ/kg.K。
     水和蒸汽温度范围273.16 K~1073.15 K (0.01℃~800℃)
     压力达0~100 MPa。
	 湿空气温度范围273.16 K~473.15 K (0.01℃~200℃)，
     压力达40000Pa~400000Pa。
     
*/ 

//************************************************************
//     作    者：                                            *
//              长沙优易软件开发有限公司(UESoft Corp.) 李国斌*
//     文件名称：                                            * 
//                DLL外部接口实现(水蒸汽、湿空气性质计算)    *
//     完成时间：                                            *
//                2006年5月                                  *
//************************************************************ 

#include "wetair.h"
#include "UEwasp.h"
#include "UEwasp67.h"
#include "UEwasp97.h"
#include <fstream>

#define MAXNUM_GAS 12
#define T0 273.15

short AIR=5;
bool WetAir_BLIF97;
double R=8.3143;//通用气体常数,J/mol.K,8314.3kJ/kmol.K
double fPSA=101325E-6;//标准大气压,MPa,fPSA=101325Pa
double fPLA=101325E-6;//当地大气压,MPa
double fVWBSA=2.5;//湿球表面空气流速,一般取>=2.5m/s
struct GasData mArrGasData[MAXNUM_GAS]=
{0,"H2",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
1,"O2",32.00,25.48,15.200,5.0620,1.31200,273,1800,1.19,\
2,"N2",28.02,28.90,-1.570,8.0810,-28.730,273,1800,0.59,\
3,"CO",28.01,28.16,1.6750,5.3720,-2.2220,273,1800,0.89,\
4,"CO2",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
5,"AIR",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
6,"H2O_S",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
7,"CH4",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
8,"C2H4",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
9,"C2H6",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
10,"C3H6",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
11,"C3H8",2.016,29.21,-1.916,-4.004,-0.8705,273,1800,1.01,\
};

int *RANGE;
double *H;
double *S;
double *P;
double *T;
double *V;
double *X;

 void _stdcall  WetAir_SETSTD_WetBulbSurface_AirVelocity(double fWetBulbSurface_AirVelocity)
{//湿球表面空气流速,m/s,fWetBulbSurface_AirVelocity
   fVWBSA=fWetBulbSurface_AirVelocity;
}

 void _stdcall  WetAir_SETSTD_atmosphere_pressure(double fLocalAtmospherePressure)
{//当地大气压fLocalAtmospherePressure
   fPLA=fLocalAtmospherePressure;
}

 void _stdcall  WetAir_GETSTD_atmosphere_pressure(double fLocalAtmospherePressure)
{
	fLocalAtmospherePressure=fPLA;
}

 void _stdcall  WetAir_SETSTD_WASP(int  STDID)
{
     if (STDID==67)
          WetAir_BLIF97=false;
      else
          WetAir_BLIF97=true ;    
}


 void _stdcall  WetAir_GETSTD_WASP(int *STDID)
{
    if (WetAir_BLIF97==false)
        *STDID=67;
    else
        *STDID=97;
}

 short _stdcall  WetAir_GetIndex(char* cName)
{//获得气体索引
	short i=0;
	 for( i=0;i<MAXNUM_GAS;i++)
	 {
		 if(cName==mArrGasData[i].cName)
			 return i;
	 }
	 return i;//此时返回最大索引加1	 
}

  double _stdcall  WetAir_GetRva()
{//获得水/空气分子量之比Rva=0.6219
	 return WetAir_GetIndex("H2O_S")/WetAir_GetIndex("AIR");
}

  double _stdcall  WetAir_GetCp(short i,double T)
{//获得气体比热Cp
	T=T+T0;
	double fCp=mArrGasData[i].a0+mArrGasData[i].a1*1e-3*T+mArrGasData[i].a2*1e-6*T*T+mArrGasData[i].a3*1e-9*T*T*T;
	fCp=fCp/mArrGasData[i].fMolecular_Weigh;//unit:kJ/kg.K
	return fCp;
}

  double _stdcall  WetAir_GetEnthalpy(short i,double T)
{//获得气体焓值Enthalpy
	return WetAir_GetCp(i,T+T0)*T;//kJ/kg
}

  double _stdcall  WetAir_D2P(double d)
{//获得水汽分压P--unit:MPa;绝对湿度D--kg/kg干空气
	return fPLA*d/(d+WetAir_GetRva());
}

  double _stdcall  WetAir_P2D(double P)
{//获得绝对湿度D--kg/kg干空气;水汽分压P--unit:MPa;
	return WetAir_GetRva()*P/(fPSA-P);
}

  double _stdcall  WetAir_GetEnthalpy_wetair(double T,double W)
{//获得湿空气焓值Enthalpy--unit:kJ/kg,=干空气的焓ha+水蒸汽的焓hv
    if (WetAir_BLIF97==true ) 
    {
        PT2H97(WetAir_D2P(W),T,H,RANGE);
    }
    else
    {
        PT2H67(WetAir_D2P(W),T,H,RANGE);
    }
	return WetAir_GetEnthalpy(AIR,T)+W*(*H);
}

//已知干球温度t(℃)和湿球温度tw(℃)，求水汽分压P(MPa)
 double _stdcall  WetAir_TT2P(double T, double Tw, double phi1)
 {//以下见《空气调节(第三版)》P13，1994年11月第三版，2006年1月第21次印刷，中国建筑工业出版社
	 //A=alpha/(r*beta*fPSA)
	 //经验公式A=(65+6.75/v)*1e-5,一般空气流速v>=2.5m/s
	 //Pq=湿球周围空气的水蒸汽压力,Pa
	 //Pqb1=湿球表面水温下的饱和水蒸汽压力,Pa，
	 //Pqb1也相当于水表面一个饱和空气薄层的水蒸汽压力,Pa
	 double Ps1;//Ps1=干球温度T对应的饱和水蒸汽压力,MPa
	 //phi1=入口相对湿度,0~1
	double A=(65+6.75/fVWBSA)*1e-5;//in general A as 6.67e-4
	double *Pqb1=NULL;
    if ( WetAir_BLIF97==true ) 
    {
        T2P97(Tw,Pqb1,RANGE);
		T2P97(T,&Ps1,RANGE);
    }
    else
    {
        T2P67(Tw,Pqb1,RANGE);
		T2P67(T,&Ps1,RANGE);
    }
	double p=*Pqb1-A*(T-Tw)*fPLA;
	phi1=p/Ps1;
	return p;
}

 //已知干球温度t(℃)和绝热饱和温度tw(℃)，求绝对湿度D(kg/kg.dry)
 double _stdcall  WetAir_TT2D(double T, double Tw, double phi1)
 {//以下见《供暖、通风及空气调节--分析与设计(原著第6版)》P43，2005年5月第1版，2005年5月第1次印刷，化学工业出版社
	 //根据稳定流动状态，列出能量守恒方程，化简得:
	 //W1*(iv1-iw)=(ia2-ia1)+Ws2*(iv2-iw)
	 //W1=绝热饱和器入口湿空气绝对湿度,kg/kg干
	 //Ws2=绝热饱和器出口饱和湿空气绝对湿度,kg/kg干
	 //ia1=绝热饱和器入口干空气焓,kJ/kg
	 //ia2=绝热饱和器出口干空气焓,kJ/kg
	 double iv1;//iv1=绝热饱和器入口水蒸汽焓,kJ/kg
	 double iv2;//iv2=绝热饱和器出口水蒸汽焓,kJ/kg
	 double iw;//iw=绝热饱和器中液体水焓,kJ/kg
	 double Pv2;//Pv2=湿球表面水温下的饱和水蒸汽压力,MPa
	 //Pv2 也相当于水表面一个饱和空气薄层的水蒸汽压力,MPa
	 double Pv1;//Pv1=绝热饱和器入口水蒸汽压力,MPa
	 double Ps1;//Ps1=绝热饱和器入口温度T对应的饱和水蒸汽压力,MPa
	 //phi1=绝热饱和器入口相对湿度,0~1
	 double fTmp=0;
	Pv1=WetAir_TT2P(T,Tw,fTmp);
	
    if ( WetAir_BLIF97==true ) 
    {
		PT2H97(Pv1,T,&iv1,RANGE);
        T2HL97(Tw,&iw,RANGE);
		T2HG97(Tw,&iv2,RANGE);
		T2P97(Tw,&Pv2,RANGE);
		T2P97(T,&Ps1,RANGE);
    }
    else
    {
		PT2H67(Pv1,T,&iv1,RANGE);
        T2HL67(Tw,&iw,RANGE);
		T2HG67(Tw,&iv2,RANGE);
		T2P67(Tw,&Pv2,RANGE);
		T2P67(T,&Ps1,RANGE);
    }
	phi1=Pv1/Ps1;//Relative Humidity(0~1)
	return (WetAir_GetEnthalpy(AIR,Tw)-WetAir_GetEnthalpy(AIR,T)+WetAir_P2D(Pv2)*(iv2-iw))/(iv1-iw);
}

