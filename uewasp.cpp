//http://www.iciba.com/搜索(ip->192.168.0.24，用完恢复老ip)
//
//动力粘度
//dynamic viscosity
//
//运动粘度
//kinematic viscosity
//
//导热系数
//thermal conductivity
//
//定熵指数
//isentropic exponent
//
//普朗特数
//Prandtl number
//
//介电常数
//dielectric constant
//
//波长
//wavelength

/*library UEWASP; 
注：
     本DLL 根据IAPWS-IF97 和IFC67
     (1985IAPS热力学性质国际骨架表)编写。
     本DLL中的各个参数的单位均是国际单位，
     其中压力的单位为MPa，温度的单位为℃。
     温度范围273.16 K~1073.15 K (0.01℃~800℃)
     压力达0~100 MPa
     
*/ 

//************************************************************
//     作    者：                                            *
//              长沙优易软件开发有限公司(UESoft Corp.) 邝祝芳*
//     文件名称：                                            * 
//                DLL外部接口实现(水和蒸汽性质计算动态连接库)*
//     完成时间：                                            *
//                2005年5月                                  *
//************************************************************ 

#include "UEwasp.h"
#include "UEwasp67.h"
#include "UEwasp97.h"
#include <math.h>
bool BLIF97;

 void _stdcall SETSTD_WASP(int  STDID)
{
     if (STDID==67)
          BLIF97=false;
      else
          BLIF97=true ;    
}


 void _stdcall GETSTD_WASP(int * STDID)
{//C++中要得到返回值需用指针,by fordao ligb on 2007.02.03
    if (BLIF97==false)
        *STDID=67;
    else
        *STDID=97;
}

 void _stdcall P2T(double P, double * T, int * RANGE )
{
    if (BLIF97==true ) 
    {
        P2T97(P,T,RANGE);
    }
    else
    {
        P2T67(P,T,RANGE);
    }
}

 void _stdcall P2HL(double P, double *H, int * RANGE)
{
    if  (BLIF97==true )
    {
        P2HL97(P,H,RANGE);
    }
    else
    {
        P2HL67(P,H,RANGE);
    }
}

 void _stdcall P2HG(double P , double * H , int * RANGE)
{
    if  ( BLIF97==true ) 
    {
        P2HG97(P,H,RANGE);
    }
    else
    {
        P2HG67(P,H,RANGE);
    }
}

 void _stdcall P2SL(double P  ,double * S , int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2SL97(P,S,RANGE);
    }
    else
    {
        P2SL67(P,S,RANGE);
    }
}

 void _stdcall P2SG(double P  ,double * S , int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2SG97(P,S,RANGE);
    }
    else
    {
        P2SG67(P,S,RANGE);
    }
}

 void _stdcall P2VL(double P  ,double * V , int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2VL97(P,V,RANGE);
    }
    else
    {
        P2VL67(P,V,RANGE);
    }
}

 void _stdcall P2VG(double P  ,double * V , int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2VG97(P,V,RANGE);
    }
    else
    {
        P2VG67(P,V,RANGE);
    }
}

 void _stdcall P2CPL(double P ,double *CP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2CPL97(P,CP,RANGE);
    }
    else
    {
        P2CPL67(P,CP,RANGE);
    }
}

 void _stdcall P2CPG(double P ,double *CP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2CPG97(P,CP,RANGE);
    }
    else
    {
        P2CPG67(P,CP,RANGE);
    }
}

 void _stdcall P2CVL(double P ,double *CV,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2CVL97(P,CV,RANGE);
    }
    else
    {
        P2CVL67(P,CV,RANGE);
    }
}

 void _stdcall P2CVG(double P ,double *CV, int* RANGE)
{
    if ( BLIF97==true ) 
    {
        P2CVG97(P,CV,RANGE);
    }
    else
    {
        P2CVG67(P,CV,RANGE);
    }
}

 void _stdcall P2EL(double P  ,double * E,int    * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2EL97(P,E,RANGE);
    }
    else
    {
        P2EL67(P,E,RANGE);
    }
}

 void _stdcall P2EG(double P  ,double * E,int    * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2EG97(P,E,RANGE);
    }
    else
    {
        P2EG67(P,E,RANGE);
    }
}

 void _stdcall P2SSPL(double P,double *SSP,int   * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2SSPL97(P,SSP,RANGE);
    }
    else
    {
        P2SSPL67(P,SSP,RANGE);
    }
}

 void _stdcall P2SSPG(double P,double *SSP,int   * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2SSPG97(P,SSP,RANGE);
    }
    else
    {
        P2SSPG67(P,SSP,RANGE);
    }
}

 void _stdcall P2KSL(double P ,double *KS ,int   * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2KSL97(P,KS,RANGE);
    }
    else
    {
        P2KSL67(P,KS,RANGE);
    }
}

 void _stdcall P2KSG(double P ,double *KS ,int   * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2KSG97(P,KS,RANGE);
    }
    else
    {
        P2KSG67(P,KS,RANGE);
    }
}

 void _stdcall P2ETAL(double P,double *ETA,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2ETAL97(P,ETA,RANGE);
    }
    else
    {
        P2ETAL67(P,ETA,RANGE);
    }
}

 void _stdcall P2ETAG(double P,double *ETA,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2ETAG97(P,ETA,RANGE);
    }
    else
    {
        P2ETAG67(P,ETA,RANGE);
    }
}


 void _stdcall P2UL(double P  ,double *U,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2UL97(P,U,RANGE);
    }
    else
    {
        P2UL67(P,U,RANGE);
    }
}



 void _stdcall P2UG(double P  ,double *U,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2UG97(P,U,RANGE);
    }
    else
    {
        P2UG67(P,U,RANGE);
    }
}

 void _stdcall P2RAMDL(double P ,double *RAMD,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2RAMDL97(P,RAMD,RANGE);
    }
    else
    {
        P2RAMDL67(P,RAMD,RANGE);
    }
}

 void _stdcall P2RAMDG(double P ,double *RAMD,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2RAMDG97(P,RAMD,RANGE);
    }
    else
    {
        P2RAMDG67(P,RAMD,RANGE);
    }
}

 void _stdcall P2PRNL(double P  ,double *PRN,  int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2PRNL97(P,PRN,RANGE);
    }
    else
    {
        P2PRNL67(P,PRN,RANGE);
    }
}

 void _stdcall P2PRNG(double P  ,double *PRN,  int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2PRNG97(P,PRN,RANGE);
    }
    else
    {
        P2PRNG67(P,PRN,RANGE);
    }
}

 void _stdcall P2EPSL(double P  ,double *EPS,  int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2EPSL97(P,EPS,RANGE);
    }
    else
    {
        P2EPSL67(P,EPS,RANGE);
    }
}

 void _stdcall P2EPSG(double P  ,double *EPS,  int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2EPSG97(P,EPS,RANGE);
    }
    else
    {
        P2EPSG67(P,EPS,RANGE);
    }
}


 void _stdcall P2NL(double P, double Lamd,double *N,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2NL97(P,Lamd,N,RANGE);
    }
    else
    {
        P2NL67(P,Lamd,N,RANGE);
    }
}

 void _stdcall P2NG(double P, double Lamd,double *N,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2NG97(P,Lamd,N,RANGE);
    }
    else
    {
        P2NG67(P,Lamd,N,RANGE);
    }
}



 void _stdcall P2L(double P   ,double * T,double * H, double *S, double *V, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2L97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        P2L67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall P2G(double P   ,double * T,double * H, double *S, double *V, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        P2G97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        P2G67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall PT2H(double P, double T,double *H,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2H97(P,T,H,RANGE);
    }
    else
    {
        PT2H67(P,T,H,RANGE);
    }
}

 void _stdcall PT2S(double P, double T,double *S,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2S97(P,T,S,RANGE);
    }
    else
    {
        PT2S67(P,T,S,RANGE);
    }
}

 void _stdcall PT2V(double P, double T,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2V97(P,T,V,RANGE);
    }
    else
    {
        PT2V67(P,T,V,RANGE);
    }
}

 void _stdcall PT2X(double P, double T,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2X97(P,T,X,RANGE);
    }
    else
    {
        PT2X67(P,T,X,RANGE);
    }
}


 void _stdcall PT2CP(double P, double T,double *CP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2CP97(P,T,CP,RANGE);
    }
    else
    {
        PT2CP67(P,T,CP,RANGE);
    }
}

 void _stdcall PT2CV(double P, double T,double *CV,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2CV97(P,T,CV,RANGE);
    }
    else
    {
        PT2CV67(P,T,CV,RANGE);
    }
}

 void _stdcall PT2E(double P, double T,double *E,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2E97(P,T,E,RANGE);
    }
    else
    {
        PT2E67(P,T,E,RANGE);
    }
}

 void _stdcall PT2SSP(double P, double T,double *SSP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2SSP97(P,T,SSP,RANGE);
    }
    else
    {
        PT2SSP67(P,T,SSP,RANGE);
    }
}


 void _stdcall PT2KS(double P, double T,double *KS,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2KS97(P,T,KS,RANGE);
    }
    else
    {
        PT2KS67(P,T,KS,RANGE);
    }
}

 void _stdcall PT2ETA(double P, double T,double *ETA,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2ETA97(P,T,ETA,RANGE);
    }
    else
    {
        PT2ETA67(P,T,ETA,RANGE);
    }
}

 void _stdcall PT2U(double P, double T,double *U,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2U97(P,T,U,RANGE);
    }
    else
    {
        PT2U67(P,T,U,RANGE);
    }
}

 void _stdcall PT2RAMD(double P, double T,double *RAMD,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2RAMD97(P,T,RAMD,RANGE);
    }
    else
    {
        PT2RAMD67(P,T,RAMD,RANGE);
    }
}


 void _stdcall PT2PRN(double P, double T,double *PRN,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2PRN97(P,T,PRN,RANGE);
    }
    else
    {
        PT2PRN67(P,T,PRN,RANGE);
    }
}

 void _stdcall PT2EPS(double P, double T,double *EPS,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2EPS97(P,T,EPS,RANGE);
    }
    else
    {
        PT2EPS67(P,T,EPS,RANGE);
    }
}

 void _stdcall PT2N(double P, double T,double LAMD, double *N,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT2N97(P,T,LAMD,N,RANGE);
    }
    else
    {
        PT2N67(P,T,LAMD,N,RANGE);
    }
}

 void _stdcall PT(double P, double T, double * H, double *S, double *V, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PT97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        PT67(P,T,H,S,V,X,RANGE);
    }
}
//6.2 Supplementary Equation for the Metastable-Vapor Region by ligb on 2007.08.14
void _stdcall PT2MV(double P,double T,double *H,double *S,double *V,double *E,double *CP,double *CV,double *SSP,int *Range)
{
//	Equation (18) is valid in the metastable-vapor region from the saturated vapor line to the
//5 % equilibrium moisture line (determined from the equilibrium h￠ and h￠￠ values) at pressures
//from the triple-point pressure, see Eq. (9), up to 10 MPa.
	PX2H(P,0.95,H,Range);
	double H1=*H;
	double *T1,t;
	T1 = &t;
	P2T(P,T1,Range);
	*H=PT2HReg2MV(P,T);
	if(P<=10.0 && *H>=H1 && *T1>=T)
		*Range=128;//存在亚稳态参数
	else
		*Range=0;//不存在亚稳态参数
// 	if(fabs(*T1-T)<10e-6)
// 	{
// 
// 		if(P<=10.0&&T<=310)
// 			*Range=128;//存在亚稳态参数
// 		else
// 			*Range=0;//不存在亚稳态参数
// 	}
// 	else 
// 	{
// 		T2P(T,P1,Range);
// 		if(fabs(*P1-P)<10e-6&&P<=10&&T<=310)
// 			*Range=128;
// 		else if(*T1!=T&&*P1!=P&&P<=10&&T<=310)
// 			*Range=128;//存在亚稳态参数
// 		else
// 			*Range=0;
// 	}
// 	*V=PT2VReg2MV(P,T);
// 	*E=PT2EReg2MV(P,T);
// 	*S=PT2SReg2MV(P,T);
// 	*CP=PT2CPReg2MV(P,T);
// 	*SSP=PT2SSPReg2MV(P,T);
// 	*CV=PT2CVReg2MV(P,T);
// 
//	CString strPT;
//	CString strMsg1,strMsg;
//	strPT.Format("\r\nP,T= %11.9e, %11.9e\n\r",P,T);
//	strMsg1.Format("V,H,E,S,CP,SSP,CV= %11.9e, %11.9e, %11.9e, %11.9e, %11.9e, %11.9e, %11.9e\n",*V,*H,*E,*S,*CP,*SSP,*CV);
//	strMsg=strPT+strMsg1;
//	AfxMessageBox(strPT+strMsg1);
	return;
}
 void _stdcall PT2HMV(double P,double T,double *H)
 {
	*H=PT2HReg2MV(P,T);

 }
 void _stdcall PT2VMV(double P,double T,double *V)
 {
	*V=PT2VReg2MV(P,T);

 }
 void _stdcall PT2EMV(double P,double T,double *E)
 {
	*E=PT2EReg2MV(P,T);
 }
 void _stdcall PT2SMV(double P,double T,double *S)
 {
	*S=PT2SReg2MV(P,T);
 }
 void _stdcall PT2CPMV(double P,double T,double *CP)
 {
     *CP=PT2CPReg2MV(P,T);
 }
 void _stdcall PT2SSPMV(double P,double T,double *SSP)
 {
	 *SSP=PT2SSPReg2MV(P,T);
 }
 void _stdcall PT2CVMV(double P,double T,double *CV)
 {
	 *CV=PT2CVReg2MV(P,T);
 }

 void _stdcall PH2T(double P, double H,double *T,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PH2T97(P,H,T,RANGE);
    }
    else
    {
        PH2T67(P,H,T,RANGE);
    }
}

 void _stdcall PH2S(double P, double H,double *S,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PH2S97(P,H,S,RANGE);
    }
    else
    {
        PH2S67(P,H,S,RANGE);
    }
}

 void _stdcall PH2V(double P, double H,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PH2V97(P,H,V,RANGE);
    }
    else
    {
        PH2V67(P,H,V,RANGE);
    }
}

 void _stdcall PH2X(double P, double H,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PH2X97(P,H,X,RANGE);
    }
    else
    {
        PH2X67(P,H,X,RANGE);
    }
}
 //Add 2007.08.28
void _stdcall PH2SSP(double P, double H,double *SSP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PH2SSP97(P,H,SSP,RANGE);
    }
    else
    {
        PH2SSP67(P,H,SSP,RANGE);
    }

}

 void _stdcall PH(double P  ,double *T,double H, double * S,double *V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PH97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        PH67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall PS2T(double P, double S,double *T,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PS2T97(P,S,T,RANGE);
    }
    else
    {
        PS2T67(P,S,T,RANGE);
    }
}

 void _stdcall PS2H(double P, double S,double *H,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PS2H97(P,S,H,RANGE);
    }
    else
    {
        PS2H67(P,S,H,RANGE);
    }
}

 void _stdcall PS2V(double P, double S,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PS2V97(P,S,V,RANGE);
    }
    else
    {
        PS2V67(P,S,V,RANGE);
    }
}

 void _stdcall PS2X(double P, double S,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PS2X97(P,S,X,RANGE);
    }
    else
    {
        PS2X67(P,S,X,RANGE);
    }
}

 void _stdcall PS(double P  ,double *T, double *H,double S,double *V, double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PS97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        PS67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall PV2T(double P, double V,double *T,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PV2T97(P,V,T,RANGE);
    }
    else
    {
        PV2T67(P,V,T,RANGE);
    }
}

 void _stdcall PV2H(double P, double V,double *H,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PV2H97(P,V,H,RANGE);
    }
    else
    {
        PV2H67(P,V,H,RANGE);
    }
}

 void _stdcall PV2S(double P, double V,double *S,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PV2S97(P,V,S,RANGE);
    }
    else
    {
        PV2S67(P,V,S,RANGE);
    }
}

 void _stdcall PV2X(double P, double V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PV2X97(P,V,X,RANGE);
    }
    else
    {
        PV2X67(P,V,X,RANGE);
    }
}

 void _stdcall PV(double P  ,double *T, double *H, double *S,double V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PV97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        PV67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall PX2T(double P, double X,double *T,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PX2T97(P,X,T,RANGE);
    }
    else
    {
        PX2T67(P,X,T,RANGE);
    }
}

 void _stdcall PX2H(double P, double X,double *H,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PX2H97(P,X,H,RANGE);
    }
    else
    {
        PX2H67(P,X,H,RANGE);
    }
}

 void _stdcall PX2S(double P, double X,double *S,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PX2S97(P,X,S,RANGE);
    }
    else
    {
        PX2S67(P,X,S,RANGE);
    }
}

 void _stdcall PX2V(double P, double X,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PX2V97(P,X,V,RANGE);
    }
    else
    {
        PX2V67(P,X,V,RANGE);
    }
}

 void _stdcall PX(double P  ,double *T, double *H, double *S, double *V,double X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        PX97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        PX67(P,T,H,S,V,X,RANGE);
    }
}



 void _stdcall T2P(double T  ,double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2P97(T,P,RANGE);
    }
    else
    {
        T2P67(T,P,RANGE);
    }
}

 void _stdcall T2HL(double T  ,double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2HL97(T,H,RANGE);
    }
    else
    {
        T2HL67(T,H,RANGE);
    }
}

 void _stdcall T2HG(double T  ,double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2HG97(T,H,RANGE);
    }
    else
    {
        T2HG67(T,H,RANGE);
    }
}

 void _stdcall T2SL(double T  ,double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2SL97(T,S,RANGE);
    }
    else
    {
        T2SL67(T,S,RANGE);
    }
}

 void _stdcall T2SG(double T  ,double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2SG97(T,S,RANGE);
    }
    else
    {
        T2SG67(T,S,RANGE);
    }
}

 void _stdcall T2VL(double T  ,double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2VL97(T,V,RANGE);
    }
    else
    {
        T2VL67(T,V,RANGE);
    }
}

 void _stdcall T2VG(double T  ,double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2VG97(T,V,RANGE);
    }
    else
    {
        T2VG67(T,V,RANGE);
    }
}

 void _stdcall T2CPL(double T  ,double *CP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2CPL97(T,CP,RANGE);
    }
    else
    {
        T2CPL67(T,CP,RANGE);
    }
}

 void _stdcall T2CPG(double T ,double *CP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2CPG97(T,CP,RANGE);
    }
    else
    {
        T2CPG67(T,CP,RANGE);
    }
}

 void _stdcall T2CVL(double T  ,double *CV,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2CVL97(T,CV,RANGE);
    }
    else
    {
        T2CVL67(T,CV,RANGE);
    }
}

 void _stdcall T2CVG(double T  ,double *CV,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2CVG97(T,CV,RANGE);
    }
    else
    {
        T2CVG67(T,CV,RANGE);
    }
}

 void _stdcall T2EL(double T  ,double *E,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2EL97(T,E,RANGE);
    }
    else
    {
        T2EL67(T,E,RANGE);
    }
}

 void _stdcall T2EG(double T  ,double *E,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2EG97(T,E,RANGE);
    }
    else
    {
        T2EG67(T,E,RANGE);
    }
}

 void _stdcall T2SSPL(double T  ,double *SSP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2SSPL97(T,SSP,RANGE);
    }
    else
    {
        T2SSPL67(T,SSP,RANGE);
    }
}

 void _stdcall T2SSPG(double T  ,double *SSP,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2SSPG97(T,SSP,RANGE);
    }
    else
    {
        T2SSPG67(T,SSP,RANGE);
    }
}

 void _stdcall T2KSL(double T  ,double *KS,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2KSL97(T,KS,RANGE);
    }
    else
    {
        T2KSL67(T,KS,RANGE);
    }
}

 void _stdcall T2KSG(double T  ,double *KS,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2KSG97(T,KS,RANGE);
    }
    else
    {
        T2KSG67(T,KS,RANGE);
    }
}

 void _stdcall T2ETAL(double T  ,double *ETA,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2ETAL97(T,ETA,RANGE);
    }
    else
    {
        T2ETAL67(T,ETA,RANGE);
    }
}

 void _stdcall T2ETAG(double T  ,double *ETA,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2ETAG97(T,ETA,RANGE);
    }
    else
    {
        T2ETAG67(T,ETA,RANGE);
    }
}

 void _stdcall T2UL(double T  ,double *U,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2UL97(T,U,RANGE);
    }
    else
    {
        T2UL67(T,U,RANGE);
    }
}

 void _stdcall T2UG(double T  ,double *U,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2UG97(T,U,RANGE);
    }
    else
    {
        T2UG67(T,U,RANGE);
    }
}

 void _stdcall T2RAMDL(double T  ,double *RAMD,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2RAMDL97(T,RAMD,RANGE);
    }
    else
    {
        T2RAMDL67(T,RAMD,RANGE);
    }
}

 void _stdcall T2RAMDG(double T  ,double *RAMD,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2RAMDG97(T,RAMD,RANGE);
    }
    else
    {
        T2RAMDG67(T,RAMD,RANGE);
    }
}

 void _stdcall T2PRNL(double T  ,double *PRN,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2PRNL97(T,PRN,RANGE);
    }
    else
    {
        T2PRNL67(T,PRN,RANGE);
    }
}

 void _stdcall T2PRNG(double T  ,double *PRN,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2PRNG97(T,PRN,RANGE);
    }
    else
    {
        T2PRNG67(T,PRN,RANGE);
    }
}

 void _stdcall T2EPSL(double T  ,double *EPS,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2EPSL97(T,EPS,RANGE);
    }
    else
    {
        T2EPSL67(T,EPS,RANGE);
    }
}

 void _stdcall T2EPSG(double T  ,double *EPS,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2EPSG97(T,EPS,RANGE);
    }
    else
    {
        T2EPSG67(T,EPS,RANGE);
    }
}

 void _stdcall T2NL(double T, double Lamd,double *N,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2NL97(T,Lamd,N,RANGE);
    }
    else
    {
        T2NL67(T,Lamd,N,RANGE);
    }
}

 void _stdcall T2NG(double T, double Lamd,double *N,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2NG97(T,Lamd,N,RANGE);
    }
    else
    {
        T2NG67(T,Lamd,N,RANGE);
    }
}

 void _stdcall T2SURFT(double T  ,double *SURFT,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2SURFT97(T,SURFT,RANGE);
    }
    else
    {
        T2SURFT67(T,SURFT,RANGE);
    }
}

 void _stdcall T2L(double *P,double T  ,double *H,double *S,double *V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2L97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        T2L67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall T2G(double *P,double T  ,double *H,double *S,double *V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        T2G97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        T2G67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall TH2PHP(double T, double H,double *P,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2PHP97(T,H,P,RANGE);
    }
    else
    {
        TH2PHP67(T,H,P,RANGE);
    }
}

 void _stdcall TH2PLP(double T, double H,double *P,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2PLP97(T,H,P,RANGE);
    }
    else
    {
        TH2PLP67(T,H,P,RANGE);
    }
}

 void _stdcall TH2P(double T, double H,double *P,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2P97(T,H,P,RANGE);
    }
    else
    {
        TH2P67(T,H,P,RANGE);
    }
}

 void _stdcall TH2SHP(double T, double H,double *S,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2SHP97(T,H,S,RANGE);
    }
    else
    {
        TH2SHP67(T,H,S,RANGE);
    }
}

 void _stdcall TH2SLP(double T, double H,double *S,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2SLP97(T,H,S,RANGE);
    }
    else
    {
        TH2SLP67(T,H,S,RANGE);
    }
}

 void _stdcall TH2S(double T, double H,double *S,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2S97(T,H,S,RANGE);
    }
    else
    {
        TH2S67(T,H,S,RANGE);
    }
}

 void _stdcall TH2VHP(double T, double H,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2VHP97(T,H,V,RANGE);
    }
    else
    {
        TH2VHP67(T,H,V,RANGE);
    }
}

 void _stdcall TH2VLP(double T, double H,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2VLP97(T,H,V,RANGE);
    }
    else
    {
        TH2VLP67(T,H,V,RANGE);
    }
}

 void _stdcall TH2V(double T, double H,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2V97(T,H,V,RANGE);
    }
    else
    {
        TH2V67(T,H,V,RANGE);
    }
}

 void _stdcall TH2XHP(double T, double H,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2XHP97(T,H,X,RANGE);
    }
    else
    {
        TH2XHP67(T,H,X,RANGE);
    }
}

 void _stdcall TH2XLP(double T, double H,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2XLP97(T,H,X,RANGE);
    }
    else
    {
        TH2XLP67(T,H,X,RANGE);
    }
}

 void _stdcall TH2X(double T, double H,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH2X97(T,H,X,RANGE);
    }
    else
    {
        TH2X67(T,H,X,RANGE);
    }
}

 void _stdcall TH(double *P,double T, double H,double * S,double *V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TH97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        TH67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall THHP(double *P,double T, double H,double * S,double *V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        THHP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        THHP67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall THLP(double *P,double T, double H,double * S,double *V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        THLP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        THLP67(P,T,H,S,V,X,RANGE);
    }
}


 void _stdcall TS2PHP(double T, double S,double *P,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2PHP97(T,S,P,RANGE);
    }
    else
    {
        TS2PHP67(T,S,P,RANGE);
    }
}

 void _stdcall TS2PLP(double T, double S,double *P,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2PLP97(T,S,P,RANGE);
    }
    else
    {
        TS2PLP67(T,S,P,RANGE);
    }
}

 void _stdcall TS2P(double T, double S,double *P,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2P97(T,S,P,RANGE);
    }
    else
    {
        TS2P67(T,S,P,RANGE);
    }
}

 void _stdcall TS2HHP(double T, double S,double *H,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2HHP97(T,S,H,RANGE);
    }
    else
    {
        TS2HHP67(T,S,H,RANGE);
    }
}

 void _stdcall TS2HLP(double T, double S,double *H,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2HLP97(T,S,H,RANGE);
    }
    else
    {
        TS2HLP67(T,S,H,RANGE);
    }
}



 void _stdcall TS2H(double T, double S,double *H,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2H97(T,S,H,RANGE);
    }
    else
    {
        TS2H67(T,S,H,RANGE);
    }
}

 void _stdcall TS2VHP(double T, double S,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2VHP97(T,S,V,RANGE);
    }
    else
    {
        TS2VHP67(T,S,V,RANGE);
    }
}

 void _stdcall TS2VLP(double T, double S,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2VLP97(T,S,V,RANGE);
    }
    else
    {
        TS2VLP67(T,S,V,RANGE);
    }
}



 void _stdcall TS2V(double T, double S,double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2V97(T,S,V,RANGE);
    }
    else
    {
        TS2V67(T,S,V,RANGE);
    }
}

 void _stdcall TS2X(double T, double S,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS2X97(T,S,X,RANGE);
    }
    else
    {
        TS2X67(T,S,X,RANGE);
    }
}


 void _stdcall TS(double *P,double T  ,double *H,double S, double *V, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TS97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        TS67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall TSHP(double *P, double T, double *H, double S, double *V, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TSHP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        TSHP67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall TSLP(double *P, double T,  double *H, double S, double *V,  double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TSLP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        TSLP67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall TV2P(double T, double V, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TV2P97(T,V,P,RANGE);
    }
    else
    {
        TV2P67(T,V,P,RANGE);
    }
}

 void _stdcall TV2H(double T, double V, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TV2H97(T,V,H,RANGE);
    }
    else
    {
        TV2H67(T,V,H,RANGE);
    }
}

 void _stdcall TV2S(double T, double V, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TV2S97(T,V,S,RANGE);
    }
    else
    {
        TV2S67(T,V,S,RANGE);
    }
}

 void _stdcall TV2X(double T, double V,double *X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TV2X97(T,V,X,RANGE);
    }
    else
    {
        TV2X67(T,V,X,RANGE);
    }
}

 void _stdcall TV(double *P,double T, double *H, double *S , double V , double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TV97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        TV67(P,T,H,S,V,X,RANGE);
    }
}


 void _stdcall TX2P(double T, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TX2P97(T,X,P,RANGE);
    }
    else
    {
        TX2P67(T,X,P,RANGE);
    }
}

 void _stdcall TX2H(double T, double X, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TX2H97(T,X,H,RANGE);
    }
    else
    {
        TX2H67(T,X,H,RANGE);
    }
}

 void _stdcall TX2S(double T, double X, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TX2S97(T,X,S,RANGE);
    }
    else
    {
        TX2S67(T,X,S,RANGE);
    }
}

 void _stdcall TX2V(double T, double X, double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TX2V97(T,X,V,RANGE);
    }
    else
    {
        TX2V67(T,X,V,RANGE);
    }
}

 void _stdcall TX(double *P, double T, double *H, double *S, double *V, double X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        TX97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        TX67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall H2TL(double H, double *T,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        H2TL97(H,T,RANGE);
    }
    else
    {
        H2TL67(H,T,RANGE);
    }
}

 void _stdcall HS2P(double H, double S, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HS2P97(H,S,P,RANGE);
    }
    else
    {
        HS2P67(H,S,P,RANGE);
    }
}

 void _stdcall HS2T(double H, double S, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HS2T97(H,S,T,RANGE);
    }
    else
    {
        HS2T67(H,S,T,RANGE);
    }
}

 void _stdcall HS2V(double H, double S, double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HS2V97(H,S,V,RANGE);
    }
    else
    {
        HS2V67(H,S,V,RANGE);
    }
}

 void _stdcall HS2X(double H, double S, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HS2X97(H,S,X,RANGE);
    }
    else
    {
        HS2X67(H,S,X,RANGE);
    }
}

 void _stdcall HS(double *P, double *T, double H, double S,double *V,  double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HS97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        HS67(P,T,H,S,V,X,RANGE);
    }
}



 void _stdcall HV2P(double H, double V, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HV2P97(H,V,P,RANGE);
    }
    else
    {
        HV2P67(H,V,P,RANGE);
    }
}

 void _stdcall HV2T(double H, double V, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HV2T97(H,V,T,RANGE);
    }
    else
    {
        HV2T67(H,V,T,RANGE);
    }
}

 void _stdcall HV2S(double H, double V, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HV2S97(H,V,S,RANGE);
    }
    else
    {
        HV2S67(H,V,S,RANGE);
    }
}

 void _stdcall HV2X(double H, double V, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HV2X97(H,V,X,RANGE);
    }
    else
    {
        HV2X67(H,V,X,RANGE);
    }
}

 void _stdcall HV(double *P, double *T, double H, double *S, double V, double *X,  int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HV97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        HV67(P,T,H,S,V,X,RANGE);
    }
}


 void _stdcall HX2T(double H, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2T97(H,X,T,RANGE);
    }
    else
    {
        HX2T67(H,X,T,RANGE);
    }
}

 void _stdcall HX2P(double H, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2P97(H,X,P,RANGE);
    }
    else
    {
        HX2P67(H,X,P,RANGE);
    }
}


 void _stdcall HX2S(double H, double X, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2S97(H,X,S,RANGE);
    }
    else
    {
        HX2S67(H,X,S,RANGE);
    }
}

 void _stdcall HX2V(double H, double X, double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2V97(H,X,V,RANGE);
    }
    else
    {
        HX2V67(H,X,V,RANGE);
    }
}

 void _stdcall HX(double *P, double *T, double H, double *S, double *V, double X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        HX67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall HX2THP(double H, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2THP97(H,X,T,RANGE);
    }
    else
    {
        HX2THP67(H,X,T,RANGE);
    }
}

 void _stdcall HX2PHP(double H, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2PHP97(H,X,P,RANGE);
    }
    else
    {
        HX2PHP67(H,X,P,RANGE);
    }
}


 void _stdcall HX2SHP(double H, double X, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2SHP97(H,X,S,RANGE);
    }
    else
    {
        HX2SHP67(H,X,S,RANGE);
    }
}

 void _stdcall HX2VHP(double H, double X, double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2VHP97(H,X,V,RANGE);
    }
    else
    {
        HX2VHP67(H,X,V,RANGE);
    }
}

 void _stdcall HXHP(double *P, double *T, double H, double *S, double *V, double X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HXHP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        HXHP67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall HX2TLP(double H, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2TLP97(H,X,T,RANGE);
    }
    else
    {
        HX2TLP67(H,X,T,RANGE);
    }
}

 void _stdcall HX2PLP(double H, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2PLP97(H,X,P,RANGE);
    }
    else
    {
        HX2PLP67(H,X,P,RANGE);
    }
}


 void _stdcall HX2SLP(double H, double X, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2SLP97(H,X,S,RANGE);
    }
    else
    {
        HX2SLP67(H,X,S,RANGE);
    }
}

 void _stdcall HX2VLP(double H, double X, double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HX2VLP97(H,X,V,RANGE);
    }
    else
    {
        HX2VLP67(H,X,V,RANGE);
    }
}

 void _stdcall HXLP(double *P, double *T, double H, double *S,  double *V, double X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        HXLP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        HXLP67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall S2TG(double S, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        S2TG97(S,T,RANGE);
    }
    else
    {
        S2TG67(S,T,RANGE);
    }
}


 void _stdcall SV2P(double S, double V, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SV2P97(S,V,P,RANGE);
    }
    else
    {
        SV2P67(S,V,P,RANGE);
    }
}

 void _stdcall SV2T(double S, double V, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SV2T97(S,V,T,RANGE);
    }
    else
    {
        SV2T67(S,V,T,RANGE);
    }
}

 void _stdcall SV2H(double S, double V, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SV2H97(S,V,H,RANGE);
    }
    else
    {
        SV2H67(S,V,H,RANGE);
    }
}

 void _stdcall SV2X(double S, double V, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SV2X97(S,V,X,RANGE);
    }
    else
    {
        SV2X67(S,V,X,RANGE);
    }
}

 void _stdcall SV(double *P, double *T, double *H, double S, double V, double *X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SV97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        SV67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall SX2T(double S, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2T97(S,X,T,RANGE);
    }
    else
    {
        SX2T67(S,X,T,RANGE);
    }
}

 void _stdcall SX2P(double S, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2P97(S,X,P,RANGE);
    }
    else
    {
        SX2P67(S,X,P,RANGE);
    }
}

 void _stdcall SX2H(double S, double X, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2H97(S,X,H,RANGE);
    }
    else
    {
        SX2H67(S,X,H,RANGE);
    }
}

 void _stdcall SX2V(double S, double X, double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2V97(S,X,V,RANGE);
    }
    else
    {
        SX2V67(S,X,V,RANGE);
    }
}

 void _stdcall SX(double *P, double *T, double *H, double S, double *V, double X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        SX67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall SX2THP(double S, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2THP97(S,X,T,RANGE);
    }
    else
    {
        SX2THP67(S,X,T,RANGE);
    }
}

 void _stdcall SX2PHP(double S, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2PHP97(S,X,P,RANGE);
    }
    else
    {
        SX2PHP67(S,X,P,RANGE);
    }
}

 void _stdcall SX2HHP(double S, double X, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2HHP97(S,X,H,RANGE);
    }
    else
    {
        SX2HHP67(S,X,H,RANGE);
    }
}

 void _stdcall SX2VHP(double S, double X, double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2VHP97(S,X,V,RANGE);
    }
    else
    {
        SX2VHP67(S,X,V,RANGE);
    }
}

 void _stdcall SXHP(double *P, double *T, double *H, double S, double *V, double X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SXHP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        SXHP67(P,T,H,S,V,X,RANGE);
    }
}


 void _stdcall SX2TMP(double S, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2TMP97(S,X,T,RANGE);
    }
    else
    {
        SX2TMP67(S,X,T,RANGE);
    }
}

 void _stdcall SX2PMP(double S, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2PMP97(S,X,P,RANGE);
    }
    else
    {
        SX2PMP67(S,X,P,RANGE);
    }
}

 void _stdcall SX2HMP(double S, double X, double *H,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2HMP97(S,X,H,RANGE);
    }
    else
    {
        SX2HMP67(S,X,H,RANGE);
    }
}

 void _stdcall SX2VMP(double S, double X, double *V,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2VMP97(S,X,V,RANGE);
    }
    else
    {
        SX2VMP67(S,X,V,RANGE);
    }
}

 void _stdcall SXMP(double *P, double *T, double *H, double S, double *V, double X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SXMP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        SXMP67(P,T,H,S,V,X,RANGE);
    }
}


 void _stdcall SX2TLP(double S, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2TLP97(S,X,T,RANGE);
    }
    else
    {
        SX2TLP67(S,X,T,RANGE);
    }
}

 void _stdcall SX2PLP(double S, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2PLP97(S,X,P,RANGE);
    }
    else
    {
        SX2PLP67(S,X,P,RANGE);
    }
}

 void _stdcall SX2HLP(double S, double X, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2HLP97(S,X,H,RANGE);
    }
    else
    {
        SX2HLP67(S,X,H,RANGE);
    }
}

 void _stdcall SX2VLP(double S, double X, double *V, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SX2VLP97(S,X,V,RANGE);
    }
    else
    {
        SX2VLP67(S,X,V,RANGE);
    }
}

 void _stdcall SXLP(double *P, double *T, double *H, double S, double *V, double X, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        SXLP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        SXLP67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall V2TG(double V , double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        V2TG97(V,T,RANGE);
    }
    else
    {
        V2TG67(V,T,RANGE);
    }
}

 void _stdcall VX2T(double V, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2T97(V,X,T,RANGE);
    }
    else
    {
        VX2T67(V,X,T,RANGE);
    }
}

 void _stdcall VX2P(double V, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2P97(V,X,P,RANGE);
    }
    else
    {
        VX2P67(V,X,P,RANGE);
    }
}

 void _stdcall VX2H(double V, double X, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2H97(V,X,H,RANGE);
    }
    else
    {
        VX2H67(V,X,H,RANGE);
    }
}

 void _stdcall VX2S(double V, double X, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2S97(V,X,S,RANGE);
    }
    else
    {
        VX2S67(V,X,S,RANGE);
    }
}

 void _stdcall VX(double *P, double *T, double *H, double *S, double V, double X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        VX67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall VX2THP(double V, double X, double *T,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2THP97(V,X,T,RANGE);
    }
    else
    {
        VX2THP67(V,X,T,RANGE);
    }
}

 void _stdcall VX2PHP(double V, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2PHP97(V,X,P,RANGE);
    }
    else
    {
        VX2PHP67(V,X,P,RANGE);
    }
}

 void _stdcall VX2HHP(double V, double X, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2HHP97(V,X,H,RANGE);
    }
    else
    {
        VX2HHP67(V,X,H,RANGE);
    }
}

 void _stdcall VX2SHP(double V, double X, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2SHP97(V,X,S,RANGE);
    }
    else
    {
        VX2SHP67(V,X,S,RANGE);
    }
}

 void _stdcall VXHP(double *P, double *T, double *H, double *S, double V, double X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VXHP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        VXHP67(P,T,H,S,V,X,RANGE);
    }
}

 void _stdcall VX2TLP(double V, double X, double *T, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2TLP97(V,X,T,RANGE);
    }
    else
    {
        VX2TLP67(V,X,T,RANGE);
    }
}

 void _stdcall VX2PLP(double V, double X, double *P, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2PLP97(V,X,P,RANGE);
    }
    else
    {
        VX2PLP67(V,X,P,RANGE);
    }
}

 void _stdcall VX2HLP(double V, double X, double *H, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2HLP97(V,X,H,RANGE);
    }
    else
    {
        VX2HLP67(V,X,H,RANGE);
    }
}

 void _stdcall VX2SLP(double V, double X, double *S, int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VX2SLP97(V,X,S,RANGE);
    }
    else
    {
        VX2SLP67(V,X,S,RANGE);
    }
}

 void _stdcall VXLP(double *P, double *T, double *H, double *S,double V, double X,int * RANGE)
{
    if ( BLIF97==true ) 
    {
        VXLP97(P,T,H,S,V,X,RANGE);
    }
    else
    {
        VXLP67(P,T,H,S,V,X,RANGE);
    }
}

