//UEwasp67.cpp

//************************************************************
//     作    者：                                            *
//              长沙优易软件开发有限公司(UESoft Corp.) 邝祝芳*
//     文件名称：                                            * 
//                IFC67标准实现(水和蒸汽性质计算动态连接库)  *
//     完成时间：                                            *
//                2005年5月                                  *
//************************************************************ 

#include "math.h"
//#include "stdlib.h"
#include "UEwasp67.h"

const double    Pai = 3.1415926;
const double    DeltaVal=1E-8;

const double    RGas_Water= 461.51;
const double    TC_Water  = 647.30-273.15;
const double    PC_Water  = 22.12;
const double    SC_Water  = 4.4429;
const double    HC_Water  = 2107.4;
const double    DC_Water  = 1/0.00317;
const double    VC_Water  = 0.00317;
const double    Pt_Water = 611.213e-6;
const double    T350C     = 350.0;
const double    T000C     = 0.0;
const double    T590C     = 590.0;
const double    P350C     = 16.535;

double    AREG1[23] =
{
      6.824687741E3,      -5.422063673E2,       -2.096666205E4,               3.941286787E4,       -6.733277739E4,
      9.902381028E4,      -1.093911774E5,        8.590841667E4,              -4.511168742E4,        1.418138926E4,
     -2.017271113E3,       7.982692717,         -2.616571843E-2,              1.52241179E-3,        2.284279054E-2,
      2.421647003E2,       1.269716088E-10,      2.074838328E-7,              2.17402035E-8,        1.105710498E-9,
      1.293441934E1,       1.308119072E-5,       6.047626338E-14
};
    
double   alREG1[13] =
{
     0.0, 8.438375405E-1,       5.362162162E-4,          1.72,       7.342278489E-2,    4.97585887E-2,
     6.5371543E-1,         1.15E-6,                 1.5108E-5,  1.4188E-1,         7.002753165,
     2.995284926E-4,       2.04E-1
};

double   B9REG2 [9] =
{
     0.0, 193.6587558,-1388.522425,4126.607219,-6508.211677,
     5745.984054,-2693.088365,523.5718623, 0.0
};

double   B0REG2[9] =
{
     0.0, 28.56067796,  -54.38923329,0.4330662834,-0.6547711697,
     0.08565182058, 0.0,        0.0,          0.0
};

double    BUVREG2[9][4] =
{   {0.0,         0.0,             0.0,       0.0},
	{0.0, 0.06670375918,  1.388983801,        0.0},
	{0.0, 0.08390104328,  0.02614670893,      -0.03373439453},
     {0.0, 0.4520918904,   0.1069036614,       0.0},
     {0.0, -0.5975336707, -0.08847535803999999,0.0},
     {0.0, 0.5958051609,  -0.5159303373,       0.2075021122},
     {0.0, 0.1190610271,  -0.09867174132000001,0.0},
     {0.0, 0.1683998803,  -0.05809438001,      0.0},
     {0.0, 0.006552390126,0.0005710218649,     0.0}
    };

double    ZUVREG2[9][4] =
    {
		{0, 0, 0, 0},
     {0,13, 3, 0},
     {0,18, 2, 1},
     {0,18,10, 0},
     {0,25,14, 0},
     {0,32,28,24},
     {0,12,11, 0},
     {0,24,18, 0},
     {0,24,14, 0}
    };

double    XULREG2[4][4] =
    {
		{0,0,0,0},
     {0,14, 0,0},
     {0,19, 0,0},
     {0,54,27,0}
    };

double    BULREG2[4][4]=
    {
		{0.0,  0.0,  0.0,  0.0},
     {0.0, 0.4006073948, 0.0,         0.0},
     {0.0, 0.08636081627,0.0,         0.0},
     {0.0, -0.8532322921,0.3460208861,0.0}
    };



double    CREG3[56] =
    {
	 0.0,
     -6.8399,        -1.7226042E-2,  -7.77175039,    4.20460752,      -2.76807038,
      2.10419707,    -1.14649588,     0.223138085,   0.116250363,     -0.0820900544,
      1.94129239E-2, -1.69470576E-3, -4.311577033,   0.708636085,      1.23679455E1,
     -1.20389004E1,   5.40437422,    -0.993865043,   0.0627523182,    -7.74743016,
     -4.29885092,     4.31430538E1,  -1.41619313E1,  4.04172459,       1.55546326,
     -1.66568935,     0.324881158,    2.93655325E1,  7.94841842E-6,    8.08859747E1,
     -8.36153380E1,   3.58636517E1,   7.51895954,   -1.2616064E1,      1.09717462,
      2.12145492,    -0.546529566,    8.32875413,    2.75971776E-6,   -5.09073985E-4,
      2.10636332E2,   0.05528935335, -0.2336365955,  0.369707142,     -0.259641547,
      0.06828087013, -2.571600553E2, -1.518783715E2, 2.220723208E1,   -1.80203957E2,
      2.35709622E3,  -1.462335698E4,  4.54291663E4, -7.053556432E4,    4.381571428E4
    };
double    DREG4[14] =
    {
	 0.0,
     -1.717616747,   3.526389875,    -2.690899373,0.9070982605, -0.1138791156,
     1.301023613,    -2.642777743,   1.996765362, -0.6661557013, 8.270860589E-2,
     3.426663535E-4,-1.236521258E-3, 1.155018309E-3
    };

double    EtaH[4] =
    {
      1.0,   0.978197,   0.579829,  -0.202354
    };

double    EtaHH[6][7]=
    {
		
	{  0.5132047, 0.2151778, -0.2818107, 0.1778064, -0.0417661, 0.0, 0.0},
    {  0.3205656, 0.7317883, -1.070786, 0.460504, 0.0, -0.01578386, 0.0},
    {  0.0, 1.241044, -1.263184, 0.2340379, 0.0, 0.0, 0.0},
    {  0.0, 1.476783, 0.0, -0.4924179, 0.1600435, 0.0, -0.003629481},
    {  -0.7782567,0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {  0.1885447, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };


double    RamdAL[4]=
    {
      0.0102811, 0.0299621, 0.0156146,-0.00422464
    };

double    RamdBL[3] =
    {
     -0.397070, 0.400302, 1.06
    };

double    RamdB[3] =
    {
    0.0, -0.171587, 2.39219
    };

double    RamdCL[7] =
    {
     0.0, 0.642857,-4.11717,-6.17937, 0.00308976, 0.0822994, 10.0932
    };

double    RamdDL[5] =
    {
     0.0, 0.0701309, 0.011852, 0.00169937,-1.02
    };


double ZBRENT1_KK(double X1, double X2,double TOL, double PrmA, int Nbr);
double ZBRENT2_KK(double X1, double X2,double TOL, double PrmA, double PrmB ,int Nbr);

void   T_HMAXXREG56(double Temp_MIN,double Temp_MAX, double X ,double* T);
double PHMINTREG1(double T);

double PT2HREG4(double P,double T);
double PT2SREG4(double P,double T);
double PT2VREG4(double P,double T);

double PT2HREG3(double P,double T);
double PT2SREG3(double P,double T);
double PT2VREG3(double P, double T);

double FUNC_PS2HREG1(double P, double S, double H);
double FUNC_PV2SREG1(double P, double V, double S);

double PT2VREG3NT(double P, double T);
double PT2VREG4NT(double P, double T);

void PS2V_KK(double P, double S , double*  V, int * RANGE);
void PT2S_KK(double P, double T, double* S, int * RANGE);

void P2PRNL(double P, double* PRN , int * RANGE);

void PT2E_KK(double P, double T,double* E,int * RANGE);
void PT2V_KK(double P, double T, double* V ,int * RANGE);

void TV2P_KK(double T, double V, double * P, int * RANGE);

void T2EL(double T, double *E, int * RANGE);
void T2EG(double  T, double *E,int * RANGE);
void T2EPSL(double T , double * EPS , int * RANGE );
void T2NL(double T, double Lamd,  double *N , int * RANGE);
double TVLAMD2N_KK(  double T, double V, double Lamd );    
 


double T2P_KK(double T)
{

    double ZT,Z1,Z2,Z3,Z4;
    ZT=(T+273.15)/647.3;
    Z1=1-ZT;
    Z2=((((64.23285504-118.9646225*Z1)*Z1-168.1706546)*Z1-26.08023696)*Z1-7.691234564)*Z1;
    Z3=(20.9750676*Z1+4.16711732)*Z1+1;
    Z4=Z1/(1E9*Z1*Z1+6);
    return exp((Z2/Z3)/ZT-Z4)*PC_Water;
};

double T2PBOUND23(double T)
{
    double ZT,PL;
    ZT = (T + 273.15) / 647.3;
    PL = (15.74373327 + (19.31380707 * ZT - 34.17061978) * ZT) * 221.2;
    return PL/10.0;
};



double PT2VREG1 (double P, double T)
{
    double Beta , Seta, Kesi,Beta2,Beta3;
    double Seta2,Seta6,Seta11,Seta19,Y,Z;

    Beta=P/22.12;
    Seta=(T+273.15)/(TC_Water+273.15);


    Seta2=Seta*Seta;
    Seta6=Seta2*Seta2*Seta2;
    Seta11=pow(Seta,11);
    Seta19=pow(Seta,19);

    Beta2=Beta*Beta;
    Beta3=Beta2*Beta;
    Y=1.0-alREG1[1]*Seta2-alREG1[2]/Seta6;
    Z=Y+pow(alREG1[3]*Y*Y-2.0*alREG1[4]*Seta+2.0*alREG1[5]*Beta,0.5);

    Kesi=AREG1[11]*alREG1[5]*pow(Z,-5.0/17.0);
    Kesi=Kesi+AREG1[12]+AREG1[13]*Seta+AREG1[14]*Seta2+AREG1[15]*pow(alREG1[6]-Seta,10)+AREG1[16]/(alREG1[7]+Seta19);
    Kesi=Kesi-(AREG1[17]+Beta*(2.0*AREG1[18]+Beta*3.0*AREG1[19]))/(alREG1[8]+Seta11);
    Kesi=Kesi-AREG1[20]*(Seta19/Seta)*(alREG1[9]+Seta2)*(alREG1[11]-3.0*pow(alREG1[10]+Beta,-4));
    Kesi=Kesi+3.0*AREG1[21]*(alREG1[12]-Seta)*Beta2;
    Kesi=Kesi+4.0*AREG1[22]*Beta3/(Seta19*Seta);
    return Kesi*VC_Water;
};

double PT2SREG1(double P, double T)
{
    double Beta,Seta,Fai,Seta2,Seta6,Seta12,Seta18,Seta21 ;
    double Beta2;
    double Y,Y1,Z ;
    int IFor ;

    Beta=P/PC_Water;
    Seta=(T+273.15)/(TC_Water+273.15);

    Beta2=Beta*Beta;

    Seta2=Seta*Seta;
    Seta6=Seta2*Seta2*Seta2;
    Seta12=Seta6*Seta6;
    Seta18=Seta12*Seta6;
    Seta21=Seta18*Seta2*Seta;


    Y=1.0-alREG1[1]*Seta2-alREG1[2]/Seta6;
    Z=Y+pow(alREG1[3]*Y*Y-2.0*alREG1[4]*Seta+2.0*alREG1[5]*Beta,0.5);
    Y1=-2*alREG1[1]*Seta+6*alREG1[2]/(Seta6*Seta);

    Fai=0.0;
    for(IFor=2;IFor<=10;IFor++)
        Fai=Fai+(IFor-1)*AREG1[IFor]*pow(Seta,IFor-2);
    Fai=AREG1[0]*log(Seta)-Fai;
    Fai=Fai+AREG1[11]*((5.0*Z/12.0-(alREG1[3]-1.0)*Y)*Y1+alREG1[4])*pow(Z,-5.0/17.0);
    // ??
	Fai=Fai+Beta*(0-AREG1[13]-2.0*AREG1[14]*Seta+10.0*AREG1[15]*pow(alREG1[6]-Seta,9)+19.0*AREG1[16]/pow(alREG1[7]+Seta18*Seta,2)*Seta18 );
    Fai=Fai-11.0*pow(alREG1[8]+Seta12/Seta,-2)*Seta12/Seta2*(Beta*(AREG1[17]+Beta*(AREG1[18]+Beta*AREG1[19])));
    Fai=Fai+AREG1[20]*Seta18/Seta*(18.0*alREG1[9]+20.0*Seta2)*(pow(alREG1[10]+Beta,-3)+alREG1[11]*Beta);
    Fai=Fai+AREG1[22]*Beta2*Beta;
    Fai=Fai+20.0*AREG1[22]*Beta2*Beta2/Seta21;
    return Fai*(PC_Water*1E3*VC_Water/(TC_Water+273.15));

};

double  PT2HREG1 (double  P, double T)
{
    double Beta ,Seta, Epsil  ;
    double Seta2,Seta6,Seta12,Seta18,Beta2,Beta4  ;
    double Y,Y1,Z;
    int IFor ;

    Beta=P/PC_Water; 	
    Seta=(T+273.15)/(TC_Water+273.15);

    Seta2=(Seta*Seta);
    Seta6=(Seta2*Seta2*Seta2);
    Seta12=(Seta6*Seta6);
    Seta18=(Seta12*Seta6);
    Beta2=(Beta*Beta);
    Beta4=(Beta2*Beta2);

    Y=1.0-alREG1[1]*Seta2-alREG1[2]/Seta6;
    Z=Y+pow(alREG1[3]*Y*Y-2.0*alREG1[4]*Seta+2.0*alREG1[5]*Beta,0.5);
    Y1=-2*alREG1[1]*Seta+6*alREG1[2]/(Seta6*Seta);

    Epsil=0.0;
    for(IFor=1;IFor<=10;IFor++)
        Epsil=Epsil+(IFor-2)*AREG1[IFor]*pow(Seta,IFor-1);

    Epsil=AREG1[0]*Seta-Epsil;
    Epsil=Epsil+AREG1[11]*(Z*(17.0*(Z/29.0-Y/12.0)+5.0*Seta*Y1/12.0)+alREG1[4]*Seta-(alREG1[3]-1)*Seta*Y*Y1)*pow(Z,-5.0/17.0);
    Epsil=Epsil+(AREG1[12]-AREG1[14]*Seta2+AREG1[15]*(9*Seta+alREG1[6])*pow(alREG1[6]-Seta,9)+AREG1[16]*(20*Seta18*Seta+alREG1[7])*pow(alREG1[7]+Seta18*Seta,-2))*Beta;
    Epsil=Epsil-(12*Seta12/Seta+alREG1[8])*pow(alREG1[8]+Seta12/Seta,-2)*Beta*(AREG1[17]+Beta*(AREG1[18]+Beta*AREG1[19]));
    Epsil=Epsil+AREG1[20]*Seta18*(17.0*alREG1[9]+19.0*Seta2)*(pow(alREG1[10]+Beta,-3)+alREG1[11]*Beta);
    Epsil=Epsil+AREG1[21]*alREG1[12]*Beta2*Beta;
    Epsil=Epsil+21.0*AREG1[22]*Beta4/(Seta18*Seta2);

    return Epsil*(PC_Water*1E3*VC_Water);
};

void  PTREG1 (double  P, double T ,double * H, double *  S, double *  V, double *  X)
{
    double Beta ,Seta, Kesi,Fai,Epsil ;
    double Seta2,Seta6,Seta12,Seta18,Beta2,Beta4  ;
    double Y,Y1,Z;
    int IFor ;

    Beta=P/22.12;
    Seta=(T+273.15)/(TC_Water+273.15);

    Seta2=Seta*Seta;
    Seta6=Seta2*Seta2*Seta2;
    Seta12=Seta6*Seta6;
    Seta18=Seta12*Seta6;

    Beta2=Beta*Beta;
    Beta4=Beta2*Beta2;

    Y=1.0-alREG1[1]*Seta2-alREG1[2]/Seta6;
    Z=Y+pow(alREG1[3]*Y*Y-2.0*alREG1[4]*Seta+2.0*alREG1[5]*Beta,0.5);
    Y1=-2*alREG1[1]*Seta+6*alREG1[2]/(Seta6*Seta);

    Kesi=AREG1[11]*alREG1[5]*pow(Z,-5.0/17.0);
    Kesi=Kesi+AREG1[12]+AREG1[13]*Seta+AREG1[14]*Seta2+AREG1[15]*pow(alREG1[6]-Seta,10)+AREG1[16]/(alREG1[7]+Seta18*Seta);
    Kesi=Kesi-(AREG1[17]+Beta*(2.0*AREG1[18]+Beta*3.0*AREG1[19]))/(alREG1[8]+Seta12/Seta);
    Kesi=Kesi-AREG1[20]*(Seta18)*(alREG1[9]+Seta2)*(alREG1[11]-3.0*pow(alREG1[10]+Beta,-4));
    Kesi=Kesi+3.0*AREG1[21]*(alREG1[12]-Seta)*Beta2;
    Kesi=Kesi+4.0*AREG1[22]*Beta2*Beta/(Seta18*Seta2);

    Fai=0.0;
    for(IFor=2;IFor<=10;IFor++)
        Fai=Fai+(IFor-1)*AREG1[IFor]*pow(Seta,IFor-2);
    Fai=AREG1[0]*log(Seta)-Fai;
    Fai=Fai+AREG1[11]*((5.0*Z/12.0-(alREG1[3]-1.0)*Y)*Y1+alREG1[4])*pow(Z,-5.0/17.0);
    Fai=Fai+Beta*(0-AREG1[13]-2.0*AREG1[14]*Seta+10.0*AREG1[15]*pow(alREG1[6]-Seta,9)+19.0*AREG1[16]/pow(alREG1[7]+Seta18*Seta,2)*Seta18 );
    Fai=Fai-11.0*pow(alREG1[8]+Seta12/Seta,-2)*Seta12/Seta2*(Beta*(AREG1[17]+Beta*(AREG1[18]+Beta*AREG1[19])));
    Fai=Fai+AREG1[20]*Seta18/Seta*(18.0*alREG1[9]+20.0*Seta2)*(pow(alREG1[10]+Beta,-3)+alREG1[11]*Beta);
    Fai=Fai+AREG1[22]*Beta2*Beta;
    Fai=Fai+20.0*AREG1[22]*Beta2*Beta2/(Seta18*Seta2*Seta);


    Epsil=0.0;
    for(IFor=1;IFor<=10;IFor++)
        Epsil=Epsil+(IFor-2)*AREG1[IFor]*pow(Seta,IFor-1);

    Epsil=AREG1[0]*Seta-Epsil;
    Epsil=Epsil+AREG1[11]*(Z*(17.0*(Z/29.0-Y/12.0)+5.0*Seta*Y1/12.0)+alREG1[4]*Seta-(alREG1[3]-1)*Seta*Y*Y1)*pow(Z,-5.0/17.0);
    Epsil=Epsil+(AREG1[12]-AREG1[14]*Seta2+AREG1[15]*(9*Seta+alREG1[6])*pow(alREG1[6]-Seta,9)+AREG1[16]*(20*Seta18*Seta+alREG1[7])*pow(alREG1[7]+Seta18*Seta,-2))*Beta;
    Epsil=Epsil-(12*Seta12/Seta+alREG1[8])*pow(alREG1[8]+Seta12/Seta,-2)*Beta*(AREG1[17]+Beta*(AREG1[18]+Beta*AREG1[19]));
    Epsil=Epsil+AREG1[20]*Seta18*(17.0*alREG1[9]+19.0*Seta2)*(pow(alREG1[10]+Beta,-3)+alREG1[11]*Beta);
    Epsil=Epsil+AREG1[21]*alREG1[12]*Beta2*Beta;
    Epsil=Epsil+21.0*AREG1[22]*Beta4/(Seta18*Seta2);

    *V=Kesi*VC_Water;
    *S=Fai*(PC_Water*1E3*VC_Water/(TC_Water+273.15));
    *H=Epsil*(PC_Water*1E3*VC_Water);
    *X=0.0;
};

double PT2VREG2 (double P, double T)
{
	
const double   B = 0.7633333333;
const double   B1= 16.83599274;
const double   I1= 4.260321148;

    double Beta , Seta ,Kesi;
    double BetaL , X  ;
    double Z1 , Z2 , Z3, Z4 ,Z6, Z7  ;
    int I , J ;


    Beta=P/PC_Water;
    Seta=(T+273.15)/(TC_Water+273.15);
    BetaL=T2PBOUND23(T)/PC_Water;

    X=exp(B*(1.0-Seta));

    Z1=0.0;
    for(I=1; I<=5; I++)
    {
        Z2=0.0;
        for(J=1;J<=3;J++)
            Z2=Z2+BUVREG2[I][J]*pow(X,ZUVREG2[I][J]);
        Z1=Z1+I*pow(Beta,I-1)*Z2;
    };

    Z4=0.0;
    for (I=1;I<=3;I++)
    {
        Z2=0.0;
        Z3=0.0;
        for (J=1;J<=3;J++)
        {
            Z2=Z2+BUVREG2[I+5][J]*pow(X,ZUVREG2[I+5][J]);
            Z3=Z3+BULREG2[I][J]*pow(X,XULREG2[I][J]);
        };
        Z7=pow(Beta,I+2);
        Z4=Z4+Z2*(I+3)*Z7/((1.0+Z3*Beta*Z7)*(1.0+Z3*Beta*Z7));
    };

    Z6=0.0;
    if (Beta>=0.1) 
    {
        for (J=1;J<=7;J++)
            Z6=Z6*X+B9REG2[8-J];
        Z6=11.0*Z6*pow(Beta/BetaL,10);
    };
    Kesi=I1*Seta/Beta-Z1-Z4+Z6;
    return Kesi*VC_Water;
};

double PT2SREG2 (double P, double T)
{
	
const  double    B = 0.7633333333;
const  double    B0= 16.83599274;
const  double    I1= 4.260321148;

    double Beta , Seta, Fai  ;
    double BetaL ,  BetaLP , X  ;
    double Z0 , Z1 , Z2 , Z3, Z4 , Z5 , Z6, Z7  ;
    double ZZ0 ,Z70  ;
    int I , J , K ;

    Beta=P/PC_Water;
    Seta=(T+273.15)/(TC_Water+273.15);
    BetaL=T2PBOUND23(T)/PC_Water;
    BetaLP= -34.17061978+2.0*19.31380707*Seta;

    X=exp(B*(1.0-Seta));
    ZZ0=Beta*pow(Beta/BetaL,10);

    Z0=0.0;
   	for(I=1;I<=5;I++)
    {
        J=6-I;
        Z0=Z0*Seta+(J-1)*B0REG2[J];
    };
    Z0=Z0/Seta;

    Z1=0.0;    
	for (I=1;I<=5;I++)
		for(J=1;J<=3;J++)        
            Z1=Z1+B*pow(Beta,I)*ZUVREG2[I][J]*BUVREG2[I][J]*pow(X,ZUVREG2[I][J]);
    Z5=0.0;
    if (Beta>=0.005)
	{
		for (I=1;I<=3;I++)
        {
            Z4=0.0;
			for(J=1;J<=3;J++)
            {
                Z2=0.0;
                Z3=0.0;                
				for (K=1; K<=3;K++)
                {
                    Z70=pow(X,XULREG2[I][K])*BULREG2[I][K];
                    Z2=Z2+Z70*XULREG2[I][K];
                    Z3=Z3+Z70;
                };
                Z7=1.0/pow(Beta,I+3)+Z3;
                Z4=Z4+(ZUVREG2[I+5][J]-Z2/Z7)*BUVREG2[I+5][J]*pow(X,ZUVREG2[I+5][J]);
            };
            Z5=Z5+Z4/Z7;
        }
	}
    Z5=Z5*B;

    Z6=0.0;
    if (Beta>0.1) 
    {
        for (I=1;I<=7;I++)
        {
            J=8-I;
            Z6=Z6*X+(10.0*BetaLP/BetaL+B*(J-1))*B9REG2[J];
        };
        Z6=Z6*ZZ0;
    };
    Fai=B0*log(Seta)-I1*log(Beta)-Z0-Z1-Z5+Z6;

    return Fai*(PC_Water*1E3*VC_Water/(TC_Water+273.15));
};


double PT2HREG2 (double P, double T)
{
const  double    B = 0.7633333333;
const  double    B0= 16.83599274;

   double Beta , Seta, Epsil  ;
   double BetaL, BetaLP, X  ;
   double Z0,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8 ;
   double ZZ0 ;
   int I,J,K ;

    Beta=P/PC_Water;
    Seta=(T+273.15)/(TC_Water+273.15);
    BetaL=T2PBOUND23(T)/PC_Water;
    BetaLP= -34.17061978+2.0*19.31380707*Seta;
    X=exp(B*(1.0-Seta));


    Z0=0.0;
    for (I=1;I<=5;I++)
    {
        J=6-I;
        Z0=Z0*Seta+(J-2)*B0REG2[J];
    };
    Z1=0.0;
    
	for (I=1;I<=5;I++)
    {
        Z2=0.0;
        for (J=1;J<=3;J++)
            Z2=Z2+BUVREG2[I][J]*(1.0+B*ZUVREG2[I][J]*Seta)*pow(X,ZUVREG2[I][J]);
        Z1=Z1+Z2*pow(Beta,I);
    };
    Z5=0.0;
    if (Beta>=0.005) 
    {
        for (I=1;I<=3;I++)
        {
            Z4=0.0;
            for (J=1;J<=3;J++)
            {
                Z2=0.0;
                Z3=0.0;
                for (K=1;K<=3;K++)
                {
                    Z2=Z2+XULREG2[I][K]*BULREG2[I][K]*pow(X,XULREG2[I][K]);
                    Z3=Z3+BULREG2[I][K]*pow(X,XULREG2[I][K]);
                };
                Z7=pow(Beta,I+3);
                Z8=pow(X,ZUVREG2[I+5][J])*(1.0+(ZUVREG2[I+5][J]-Z2/(Z3+1.0/Z7))*B*Seta);
                Z4=Z4+Z8*BUVREG2[I+5][J];
            };
            Z5=Z5+Z4/(Z3+1.0/pow(Beta,I+3));
        };
    };
    Z6=0.0;
    if (Beta>=0.1) 
	{
     for (I=1;I<=7;I++)
	 {
        J=8-I;
        Z6=Z6*X+B9REG2[J]*(1.0+Seta*(10.0*BetaLP/BetaL+B*(J-1)));
	 }
	}
    ZZ0=Beta*pow(Beta/BetaL,10);

    Epsil=B0*Seta-Z0-Z1-Z5+Z6*ZZ0;
    return Epsil*(PC_Water*1E3*VC_Water);
};


void PTREG2 (double  P, double T ,double * H,double *S,double *V,double *X)
{
	
const double    B = 0.7633333333;
const double    B0= 16.83599274;
const double    I1= 4.260321148;

double    Beta , Seta ,Kesi,Epsil, Fai ;
double    BetaL,BetaLP ;
double    Z0,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8 ;
double    ZZ0,Z70  ;
int       I,J,K ;


    Beta=P/PC_Water;
    Seta=(T+273.15)/(TC_Water+273.15);
    BetaL=T2PBOUND23(T)/PC_Water;

    *X=exp(B*(1.0-Seta));

    BetaLP= -34.17061978+2.0*19.31380707*Seta;
    ZZ0=Beta*pow(Beta/BetaL,10);

    Z1=0.0;
    for (I=1;I<=5;I++)
    {
        Z2=0.0;
        for (J=1;J<=3;J++)
            Z2=Z2+BUVREG2[I][J]*pow(*X,ZUVREG2[I][J]);
        Z1=Z1+I*pow(Beta,I-1)*Z2;
    }


    Z4=0.0;
    for (I=1;I<=3;I++)
    {
        Z2=0.0;
        Z3=0.0;
        for (J=1;J<=3;J++)
        {
            Z2=Z2+BUVREG2[I+5][J]*pow(*X,ZUVREG2[I+5][J]);
            Z3=Z3+BULREG2[I][J]*pow(*X,XULREG2[I][J]);
        }
        Z7=pow(Beta,I+2);
        Z4=Z4+Z2*(I+3)*Z7/((1.0+Z3*Beta*Z7)*(1.0+Z3*Beta*Z7));
	}
    Z6=0.0;
    if (Beta>=0.1) 
    {
        for (I=1;I<=7;I++)
            Z6=Z6*(*X)+B9REG2[8-I];
        Z6=11.0*Z6*pow(Beta/BetaL,10);
    }
    Kesi=I1*Seta/Beta-Z1-Z4+Z6;

    Z0=0.0;    
    for (I=1;I<=5;I++)
    {
        J=6-I;
        Z0=Z0*Seta+(J-1)*B0REG2[J];
    };
    Z0=Z0/Seta;
    Z1=0.0;    
	for (I=1;I<=5;I++)
         for(J=1;J<=3;J++)
            Z1=Z1+B*pow(Beta,I)*ZUVREG2[I][J]*BUVREG2[I][J]*pow(*X,ZUVREG2[I][J]);
    Z5=0.0;
    if (Beta>=0.005)        
		for (I=1;I<=3;I++)
        {
            Z4=0.0;            
			for(J=1;J<=3;J++)
            {
                Z2=0.0;
                Z3=0.0;                
				for (K=1;K<=3;K++)
                {
                    Z70=pow(*X,XULREG2[I][K])*BULREG2[I][K];
                    Z2=Z2+Z70*XULREG2[I][K];
                    Z3=Z3+Z70;
                }
                Z7=1.0/pow(Beta,I+3)+Z3;
                Z4=Z4+(ZUVREG2[I+5][J]-Z2/Z7)*BUVREG2[I+5][J]*pow(*X,ZUVREG2[I+5][J]);
            }
            Z5=Z5+Z4/Z7;
        }
    Z5=Z5*B;

    Z6=0.0;
    if (Beta>0.1) 
    {
        for (I=1;I<=7;I++)
        {
            J=8-I;
            Z6=Z6*(*X) +(10.0*BetaLP/BetaL+B*(J-1))*B9REG2[J];
		}
        Z6=Z6*ZZ0;
    }
    Fai=B0*log(Seta)-I1*log(Beta)-Z0-Z1-Z5+Z6;

    Z0=0.0;    
	for (I=1;I<=5;I++)
    {
        J=6-I;
        Z0=Z0*Seta+(J-2)*B0REG2[J];
    }

    Z1=0.0;    
	for (I=1;I<=5;I++)
    {
        Z2=0.0;        
		for (J=1;J<=3;J++)
            Z2=Z2+BUVREG2[I][J]*(1.0+B*ZUVREG2[I][J]*Seta)*pow(*X,ZUVREG2[I][J]);
        Z1=Z1+Z2*pow(Beta,I);
    }

    Z5=0.0;
    if  (Beta>=0.005)
    {
        for(I=1;I<=3;I++)
        {
            Z4=0.0;            
			for(J=1;J<=3;J++)
            {
                Z2=0.0;
                Z3=0.0;                
				for (K=1;K<=3;K++)
                {
                    Z2=Z2+XULREG2[I][K]*BULREG2[I][K]*pow(*X,XULREG2[I][K]);
                    Z3=Z3+BULREG2[I][K]*pow(*X,XULREG2[I][K]);
                }
                Z7=pow(Beta,I+3);
                Z8=pow(*X,ZUVREG2[I+5][J])*(1.0+(ZUVREG2[I+5][J]-Z2/(Z3+1.0/Z7))*B*Seta);
                Z4=Z4+Z8*BUVREG2[I+5][J];
            }
            Z5=Z5+Z4/(Z3+1.0/pow(Beta,I+3));
        }
    }

    Z6=0.0;
    if (Beta>=0.1) 
    {
      for (I=1;I<=7;I++)
	  {
        J=8-I;
        Z6=Z6*(*X)+B9REG2[J]*(1.0+Seta*(10.0*BetaLP/BetaL+B*(J-1)));
	  };
    }
      ZZ0=Beta*pow(Beta/BetaL,10);

      Epsil=B0*Seta-Z0-Z1-Z5+Z6*ZZ0;

      *V=Kesi*VC_Water;
      *S=Fai*(PC_Water*1E3*VC_Water/(TC_Water+273.15));
      *H=Epsil*(PC_Water*1E3*VC_Water);
      *X=1.0;	
};


double TV2PsiREG3(double T, double V )
{
double    Seta;
double    CHI,CHI2,CHI5,CHI6;
double    Seta23,TT;
double    LN_Chi;
double    PS11,PS12, PS13;
double    PS1,PS2,PS3,PS4,PS5,PS6,PS7;
double    PSIC;

    Seta=(T+273.15)/(TC_Water+273.15);
    CHI= V/VC_Water;
    CHI2 =CHI *CHI;
    CHI5 =CHI2 *CHI2 *CHI;
    CHI6 =CHI5 *CHI;
    LN_Chi=log(CHI);

    Seta23 =pow(Seta,23);
    TT =Seta -1;

    PS11 =CREG3[0+1]+CREG3[1+1]*CHI;
    PS12 =((((CREG3[11+1]/CHI+CREG3[10+1])/CHI+CREG3[9+1])/CHI+CREG3[8+1])/CHI+CREG3[7+1])/CHI;
    PS12 =(((((PS12+CREG3[6+1])/CHI+CREG3[5+1])/CHI+CREG3[4+1])/CHI+CREG3[3+1])/CHI+CREG3[2+1])/CHI;
    PS13 =CREG3[12+1]*LN_Chi;
    PS1 =PS11+PS12+PS13;

    PS2 =CREG3[13+1]*CHI+((((CREG3[19]/CHI+CREG3[18])/CHI+CREG3[17])/CHI+CREG3[16])/CHI+CREG3[15])/CHI;
    PS2 =(PS2+CREG3[20]*LN_Chi)*TT;
    PS3 =CREG3[21]*CHI+(((((CREG3[27]/CHI+CREG3[26])/CHI+CREG3[25])/CHI+CREG3[24])/CHI+CREG3[23])/CHI+CREG3[22])/CHI;
    PS3 =(PS3+CREG3[28]*LN_Chi)*(TT*TT);

    PS4 =((((CREG3[37]/CHI+CREG3[36])/CHI+CREG3[35])/CHI +CREG3[34])/CHI +CREG3[33])/CHI;
    PS4 =(CREG3[29]*CHI+(((PS4+CREG3[32])/CHI+CREG3[31])/CHI+CREG3[30])/CHI+CREG3[38]*LN_Chi)*(TT*TT*TT);

    PS5 =(CREG3[39]+CREG3[40]/CHI5)/Seta23 *TT +CREG3[41]*Seta*log(Seta);
    PS6 =CHI6*((((CREG3[46]/Seta+CREG3[45])/Seta+CREG3[44])/Seta+CREG3[43])/Seta+CREG3[42])/Seta/Seta;

    PS7 =((((CREG3[55]*TT +CREG3[54])*TT +CREG3[53])*TT +CREG3[52])*TT +CREG3[51])*TT;
    PS7 =((((PS7 +CREG3[50])*TT +CREG3[49])*TT +CREG3[48])*TT +CREG3[47])*TT;

    PSIC=PS1 +PS2 +PS3 +PS4 +PS5 +PS6 +PS7;
    return PSIC+0;
}



double TV2PREG3(double T, double V)
{
double    Seta;
double    CHI,CHI2,CHI5,CHI6;
double    Seta23,TT;
double    DPS11CH,DPS12CH,DPS13CH;
double    DPS1CH,DPS2CH,DPS3CH,DPS4CH,DPS5CH,DPS6CH;
double    BETA3;

    Seta=(T+273.15)/(TC_Water+273.15);
    CHI= V/VC_Water;
    CHI2 =CHI *CHI;
    CHI5 =CHI2 *CHI2 *CHI;
    CHI6 =CHI5 *CHI;

    Seta23 =pow(Seta,23);
    TT =Seta -1;

    DPS11CH =(((10 *CREG3[12]/CHI+9*CREG3[11])/CHI+8*CREG3[10])/CHI+7*CREG3[9])/CHI;
    DPS12CH =((DPS11CH +6 *CREG3[8])/CHI +5 *CREG3[7])/CHI;
    DPS13CH =((((DPS12CH+4*CREG3[6])/CHI+3*CREG3[5])/CHI+2*CREG3[4])/CHI+CREG3[3])/(CHI*CHI);
    DPS1CH =CREG3[2]-DPS13CH +CREG3[13]/CHI;

    DPS2CH =((((5*CREG3[19]/CHI+4*CREG3[18])/CHI+3*CREG3[17])/CHI+2*CREG3[16])/CHI+CREG3[15])/(CHI*CHI);
    DPS2CH =(CREG3[14]-DPS2CH+CREG3[20]/CHI)*TT;

    DPS3CH =((((6*CREG3[27]/CHI+5*CREG3[26])/CHI+4*CREG3[25])/CHI+3*CREG3[24])/CHI+2*CREG3[23])/CHI;
    DPS3CH =(CREG3[21]-(DPS3CH +CREG3[22])/CHI /CHI +CREG3[28]/CHI)*TT *TT;

    DPS4CH =((((8*CREG3[37]/CHI+7*CREG3[36])/CHI+6*CREG3[35])/CHI+5*CREG3[34])/CHI+4*CREG3[33])/CHI;
    DPS4CH =(CREG3[29]-(((DPS4CH+3*CREG3[32])/CHI+2*CREG3[31])/CHI+CREG3[30])/CHI/CHI+CREG3[38]/CHI)*(TT*TT*TT);

    DPS5CH =-5 *CREG3[40]/CHI6 /Seta23 *TT;
    DPS6CH =6*CHI5*((((CREG3[46]/Seta+CREG3[45])/Seta+CREG3[44])/Seta+CREG3[43])/Seta+CREG3[42])/(Seta*Seta);

    BETA3 =-DPS1CH -DPS2CH -DPS3CH -DPS4CH -DPS5CH -DPS6CH;

    return BETA3*PC_Water;
}

double TV2SREG3(double T, double V )
{
double    Seta;
double    CHI, CHI2,CHI5, CHI6;
double    LN_Chi;
double    Seta23,TT;
double    DPS2TH,DPS3TH,DPS4TH,DPS5TH,DPS6TH,DPS7TH;
double    SEGMA3;

    Seta=(T+273.15)/(TC_Water+273.15);
    CHI= V/VC_Water;
    CHI2=CHI*CHI;
    CHI5 =CHI2 *CHI2 *CHI;
    CHI6 =CHI5 *CHI;
    Seta23 =pow(Seta,23);
    LN_Chi=log(CHI);
    TT =Seta -1;


    DPS2TH =CREG3[14]*CHI+((((CREG3[19]/CHI+CREG3[18])/CHI+CREG3[17])/CHI+CREG3[16])/CHI+CREG3[15])/CHI+CREG3[20]*LN_Chi;

    DPS3TH =((((CREG3[27]/CHI+CREG3[26])/CHI+CREG3[25])/CHI+CREG3[24])/CHI+CREG3[23])/CHI;
    DPS3TH =2*(CREG3[21]*CHI+(DPS3TH+CREG3[22])/CHI +CREG3[28]*LN_Chi)*TT;

    DPS4TH =(((((CREG3[37]/CHI+CREG3[36])/CHI+CREG3[35])/CHI+CREG3[34])/CHI+CREG3[33])/CHI+CREG3[32])/CHI;
    DPS4TH =3*(CREG3[29]*CHI+((DPS4TH+CREG3[31])/CHI+CREG3[30])/CHI+CREG3[38]*LN_Chi)*TT *TT;

    DPS5TH=(CREG3[39]+CREG3[40]/CHI5)*(-22/Seta23+23/Seta23/Seta)+CREG3[41]*(log(Seta)+1);

    DPS6TH=(((6*CREG3[46]/Seta+5*CREG3[45])/Seta+4*CREG3[44])/Seta+3*CREG3[43])/Seta;
    DPS6TH=-CHI6*(DPS6TH+2*CREG3[42])/(Seta*Seta*Seta);

    DPS7TH=((((9*CREG3[55]*TT+8*CREG3[54])*TT+7*CREG3[53])*TT+6*CREG3[52])*TT+5*CREG3[51])*TT;
    DPS7TH=(((DPS7TH+4*CREG3[50])*TT+3*CREG3[49])*TT+2*CREG3[48])*TT+CREG3[47];

    SEGMA3=-(DPS2TH+DPS3TH+DPS4TH+DPS5TH+DPS6TH+DPS7TH)+0.0000015929;

    return SEGMA3*(PC_Water*1E3*VC_Water/(TC_Water+273.15));
}



double TV2HREG3(double T, double V )
{
double    Seta;
double    CHI, CHI2,CHI5, CHI6;
double    LN_Chi;
double    Seta23,TT;
double    PS11,PS12, PS13;
double    PS1,PS2,PS3,PS4,PS5,PS6,PS7;
double    PSIC;
double    DPS11CH,DPS12CH,DPS13CH;
double    DPS1CH,DPS2CH,DPS3CH,DPS4CH,DPS5CH,DPS6CH;
double    DPS2TH,DPS3TH,DPS4TH,DPS5TH,DPS6TH,DPS7TH;
double    BETA3,SEGMA3,EPSILON3;


    Seta=(T+273.15)/(TC_Water+273.15);
    CHI=V/VC_Water;
    CHI2 =CHI *CHI;
    CHI5 =CHI2 *CHI2 *CHI;
    CHI6 =CHI5 *CHI;
    LN_Chi=log(CHI);

    Seta23 =pow(Seta,23);
    TT =Seta -1;

    PS11 =CREG3[1]+CREG3[2]*CHI;
    PS12 =((((CREG3[12]/CHI+CREG3[11])/CHI+CREG3[10])/CHI+CREG3[9])/CHI+CREG3[8])/CHI;
    PS12 =(((((PS12+CREG3[7])/CHI+CREG3[6])/CHI+CREG3[5])/CHI+CREG3[4])/CHI+CREG3[3])/CHI;
    PS13 =CREG3[13]*LN_Chi;
    PS1 =PS11+PS12+PS13;

    PS2 =CREG3[14]*CHI+((((CREG3[19]/CHI+CREG3[18])/CHI+CREG3[17])/CHI+CREG3[16])/CHI+CREG3[15])/CHI;
    PS2 =(PS2+CREG3[20]*LN_Chi)*TT;
    PS3 =CREG3[21]*CHI+(((((CREG3[27]/CHI+CREG3[26])/CHI+CREG3[25])/CHI+CREG3[24])/CHI+CREG3[23])/CHI+CREG3[22])/CHI;
    PS3 =(PS3+CREG3[28]*LN_Chi)*(TT*TT);

    PS4 =((((CREG3[37]/CHI+CREG3[36])/CHI+CREG3[35])/CHI +CREG3[34])/CHI +CREG3[33])/CHI;
    PS4 =(CREG3[29]*CHI+(((PS4+CREG3[32])/CHI+CREG3[31])/CHI+CREG3[30])/CHI+CREG3[38]*LN_Chi)*(TT*TT*TT);

    PS5 =(CREG3[39]+CREG3[40]/CHI5)/Seta23 *TT +CREG3[41]*Seta*log(Seta);
    PS6 =CHI6*((((CREG3[46]/Seta+CREG3[45])/Seta+CREG3[44])/Seta+CREG3[43])/Seta+CREG3[42])/Seta/Seta;

    PS7 =((((CREG3[55]*TT +CREG3[54])*TT +CREG3[53])*TT +CREG3[52])*TT +CREG3[51])*TT;
    PS7 =((((PS7 +CREG3[50])*TT +CREG3[49])*TT +CREG3[48])*TT +CREG3[47])*TT;

    PSIC =PS1 +PS2 +PS3 +PS4 +PS5 +PS6 +PS7;

    DPS11CH =(((10 *CREG3[12]/CHI+9*CREG3[11])/CHI+8*CREG3[10])/CHI+7*CREG3[9])/CHI;
    DPS12CH =((DPS11CH +6 *CREG3[8])/CHI +5 *CREG3[7])/CHI;
    DPS13CH =((((DPS12CH+4*CREG3[6])/CHI+3*CREG3[5])/CHI+2*CREG3[4])/CHI+CREG3[3])/(CHI*CHI);
    DPS1CH =CREG3[2]-DPS13CH +CREG3[13]/CHI;

    DPS2CH =((((5*CREG3[19]/CHI+4*CREG3[18])/CHI+3*CREG3[17])/CHI+2*CREG3[16])/CHI+CREG3[15])/(CHI*CHI);
    DPS2CH =(CREG3[14]-DPS2CH+CREG3[20]/CHI)*TT;

    DPS3CH =((((6*CREG3[27]/CHI+5*CREG3[26])/CHI+4*CREG3[25])/CHI+3*CREG3[24])/CHI+2*CREG3[23])/CHI;
    DPS3CH =(CREG3[21]-(DPS3CH +CREG3[22])/CHI /CHI +CREG3[28]/CHI)*TT *TT;

    DPS4CH =((((8*CREG3[37]/CHI+7*CREG3[36])/CHI+6*CREG3[35])/CHI+5*CREG3[34])/CHI+4*CREG3[33])/CHI;
    DPS4CH =(CREG3[29]-(((DPS4CH+3*CREG3[32])/CHI+2*CREG3[31])/CHI+CREG3[30])/CHI/CHI+CREG3[38]/CHI)*(TT*TT*TT);

    DPS5CH =-5 *CREG3[40]/CHI6 /Seta23 *TT;
    DPS6CH =6*CHI5*((((CREG3[46]/Seta+CREG3[45])/Seta+CREG3[44])/Seta+CREG3[43])/Seta+CREG3[42])/(Seta*Seta);

    BETA3 =-DPS1CH -DPS2CH -DPS3CH -DPS4CH -DPS5CH -DPS6CH;

    DPS2TH =CREG3[14]*CHI+((((CREG3[19]/CHI+CREG3[18])/CHI+CREG3[17])/CHI+CREG3[16])/CHI+CREG3[15])/CHI+CREG3[20]*LN_Chi;

    DPS3TH =((((CREG3[27]/CHI+CREG3[26])/CHI+CREG3[25])/CHI+CREG3[24])/CHI+CREG3[23])/CHI;
    DPS3TH =2*(CREG3[21]*CHI+(DPS3TH+CREG3[22])/CHI +CREG3[28]*LN_Chi)*TT;

    DPS4TH =(((((CREG3[37]/CHI+CREG3[36])/CHI+CREG3[35])/CHI+CREG3[34])/CHI+CREG3[33])/CHI+CREG3[32])/CHI;
    DPS4TH =3*(CREG3[29]*CHI+((DPS4TH+CREG3[31])/CHI+CREG3[30])/CHI+CREG3[38]*LN_Chi)*TT *TT;

    DPS5TH=(CREG3[39]+CREG3[40]/CHI5)*(-22/Seta23+23/Seta23/Seta)+CREG3[41]*(log(Seta)+1);

    DPS6TH=(((6*CREG3[46]/Seta+5*CREG3[45])/Seta+4*CREG3[44])/Seta+3*CREG3[43])/Seta;
    DPS6TH=-CHI6*(DPS6TH+2*CREG3[42])/(Seta*Seta*Seta);

    DPS7TH=((((9*CREG3[55]*TT+8*CREG3[54])*TT+7*CREG3[53])*TT+6*CREG3[52])*TT+5*CREG3[51])*TT;
    DPS7TH=(((DPS7TH+4*CREG3[50])*TT+3*CREG3[49])*TT+2*CREG3[48])*TT+CREG3[47];

    SEGMA3=-(DPS2TH+DPS3TH+DPS4TH+DPS5TH+DPS6TH+DPS7TH)+0.0000015929;


    EPSILON3=PSIC+Seta*SEGMA3+CHI*BETA3+0.000000010418-0.0000015929*Seta;


    return EPSILON3*(PC_Water*1E3*VC_Water);

}


void TVREG3(double* P, double T,double* H , double *S, double V ,double *X)
{
double    Seta;
double    CHI, CHI2,CHI5, CHI6;
double    LN_Chi;
double    Seta23,TT;
double    PS11,PS12, PS13;
double    PS1,PS2,PS3,PS4,PS5,PS6,PS7;
double    PSIC;
double    DPS11CH,DPS12CH,DPS13CH;
double    DPS1CH,DPS2CH,DPS3CH,DPS4CH,DPS5CH,DPS6CH;
double    DPS2TH,DPS3TH,DPS4TH,DPS5TH,DPS6TH,DPS7TH;
double    BETA3,SEGMA3,EPSILON3;


    Seta=(T+273.15)/(TC_Water+273.15);
    CHI= V/VC_Water;
    CHI2 =CHI *CHI;
    CHI5 =CHI2 *CHI2 *CHI;
    CHI6 =CHI5 *CHI;

    Seta23 =pow(Seta,23);
    TT =Seta -1;
    LN_Chi=log(CHI);


    PS11 =CREG3[1]+CREG3[2]*CHI;
    PS12 =((((CREG3[12]/CHI+CREG3[11])/CHI+CREG3[10])/CHI+CREG3[9])/CHI+CREG3[8])/CHI;
    PS12 =(((((PS12+CREG3[7])/CHI+CREG3[6])/CHI+CREG3[5])/CHI+CREG3[4])/CHI+CREG3[3])/CHI;
    PS13 =CREG3[13]*LN_Chi;
    PS1 =PS11+PS12+PS13;

    PS2 =CREG3[14]*CHI+((((CREG3[19]/CHI+CREG3[18])/CHI+CREG3[17])/CHI+CREG3[16])/CHI+CREG3[15])/CHI;
    PS2 =(PS2+CREG3[20]*LN_Chi)*TT;
    PS3 =CREG3[21]*CHI+(((((CREG3[27]/CHI+CREG3[26])/CHI+CREG3[25])/CHI+CREG3[24])/CHI+CREG3[23])/CHI+CREG3[22])/CHI;
    PS3 =(PS3+CREG3[28]*LN_Chi)*(TT*TT);

    PS4 =((((CREG3[37]/CHI+CREG3[36])/CHI+CREG3[35])/CHI +CREG3[34])/CHI +CREG3[33])/CHI;
    PS4 =(CREG3[29]*CHI+(((PS4+CREG3[32])/CHI+CREG3[31])/CHI+CREG3[30])/CHI+CREG3[38]*LN_Chi)*(TT*TT*TT);

    PS5 =(CREG3[39]+CREG3[40]/CHI5)/Seta23 *TT +CREG3[41]*Seta*log(Seta);
    PS6 =CHI6*((((CREG3[46]/Seta+CREG3[45])/Seta+CREG3[44])/Seta+CREG3[43])/Seta+CREG3[42])/Seta/Seta;

    PS7 =((((CREG3[55]*TT +CREG3[54])*TT +CREG3[53])*TT +CREG3[52])*TT +CREG3[51])*TT;
    PS7 =((((PS7 +CREG3[50])*TT +CREG3[49])*TT +CREG3[48])*TT +CREG3[47])*TT;

    PSIC =PS1 +PS2 +PS3 +PS4 +PS5 +PS6 +PS7;


    DPS11CH =(((10 *CREG3[12]/CHI+9*CREG3[11])/CHI+8*CREG3[10])/CHI+7*CREG3[9])/CHI;
    DPS12CH =((DPS11CH +6 *CREG3[8])/CHI +5 *CREG3[7])/CHI;
    DPS13CH =((((DPS12CH+4*CREG3[6])/CHI+3*CREG3[5])/CHI+2*CREG3[4])/CHI+CREG3[3])/(CHI*CHI);
    DPS1CH =CREG3[2]-DPS13CH +CREG3[13]/CHI;

    DPS2CH =((((5*CREG3[19]/CHI+4*CREG3[18])/CHI+3*CREG3[17])/CHI+2*CREG3[16])/CHI+CREG3[15])/(CHI*CHI);
    DPS2CH =(CREG3[14]-DPS2CH+CREG3[20]/CHI)*TT;

    DPS3CH =((((6*CREG3[27]/CHI+5*CREG3[26])/CHI+4*CREG3[25])/CHI+3*CREG3[24])/CHI+2*CREG3[23])/CHI;
    DPS3CH =(CREG3[21]-(DPS3CH +CREG3[22])/CHI /CHI +CREG3[28]/CHI)*TT *TT;

    DPS4CH =((((8*CREG3[37]/CHI+7*CREG3[36])/CHI+6*CREG3[35])/CHI+5*CREG3[34])/CHI+4*CREG3[33])/CHI;
    DPS4CH =(CREG3[29]-(((DPS4CH+3*CREG3[32])/CHI+2*CREG3[31])/CHI+CREG3[30])/CHI/CHI+CREG3[38]/CHI)*(TT*TT*TT);

    DPS5CH =-5 *CREG3[40]/CHI6 /Seta23 *TT;
    DPS6CH =6*CHI5*((((CREG3[46]/Seta+CREG3[45])/Seta+CREG3[44])/Seta+CREG3[43])/Seta+CREG3[42])/(Seta*Seta);

    BETA3 =-DPS1CH -DPS2CH -DPS3CH -DPS4CH -DPS5CH -DPS6CH;


    DPS2TH =CREG3[14]*CHI+((((CREG3[19]/CHI+CREG3[18])/CHI+CREG3[17])/CHI+CREG3[16])/CHI+CREG3[15])/CHI+CREG3[20]*LN_Chi;

    DPS3TH =((((CREG3[27]/CHI+CREG3[26])/CHI+CREG3[25])/CHI+CREG3[24])/CHI+CREG3[23])/CHI;
    DPS3TH =2*(CREG3[21]*CHI+(DPS3TH+CREG3[22])/CHI +CREG3[28]*LN_Chi)*TT;

    DPS4TH =(((((CREG3[37]/CHI+CREG3[36])/CHI+CREG3[35])/CHI+CREG3[34])/CHI+CREG3[33])/CHI+CREG3[32])/CHI;
    DPS4TH =3*(CREG3[29]*CHI+((DPS4TH+CREG3[31])/CHI+CREG3[30])/CHI+CREG3[38]*LN_Chi)*TT *TT;

    DPS5TH=(CREG3[39]+CREG3[40]/CHI5)*(-22/Seta23+23/Seta23/Seta)+CREG3[41]*(log(Seta)+1);

    DPS6TH=(((6*CREG3[46]/Seta+5*CREG3[45])/Seta+4*CREG3[44])/Seta+3*CREG3[43])/Seta;
    DPS6TH=-CHI6*(DPS6TH+2*CREG3[42])/(Seta*Seta*Seta);

    DPS7TH=((((9*CREG3[55]*TT+8*CREG3[54])*TT+7*CREG3[53])*TT+6*CREG3[52])*TT+5*CREG3[51])*TT;
    DPS7TH=(((DPS7TH+4*CREG3[50])*TT+3*CREG3[49])*TT+2*CREG3[48])*TT+CREG3[47];

    SEGMA3=-(DPS2TH+DPS3TH+DPS4TH+DPS5TH+DPS6TH+DPS7TH)+0.0000015929;



    EPSILON3=PSIC+Seta*SEGMA3+CHI*BETA3+0.000000010418-0.0000015929*Seta;


    *P=BETA3*PC_Water;
    *S=SEGMA3*(PC_Water*1E3*VC_Water/(TC_Water+273.15));
    *H=EPSILON3*(PC_Water*1E3*VC_Water);
    *X=1.0;
}


double TV2PsiREG4(double T, double V )
{
double    Seta , CHI;
double    Y , Y3 , Y31 ;
double    PSID1 , PSID2 , PSID  ;

    Seta=(T+273.15)/(TC_Water+273.15);
    CHI=V/VC_Water;

    Y=(1-Seta)/(1-62315.0/64730.0);
    Y3=Y*Y*Y;
    Y31=pow(Y,31);

    PSID1=((((DREG4[10]/CHI+DREG4[9])/CHI+DREG4[8])/CHI+DREG4[7])/CHI+DREG4[6])*Y;
    PSID2=(((DREG4[5]/CHI+DREG4[4])/CHI+DREG4[3])/CHI+DREG4[2])/CHI;
    PSID=(PSID1+PSID2+DREG4[1])*Y3+((DREG4[13]*CHI+DREG4[12])*CHI+DREG4[11])*Y31 *Y;

    return TV2PsiREG3(T,V) +PSID;
}




double TV2PREG4 (double T, double V)
{
double    Seta , CHI , BETA3;
double    Y , Y3 , Y31;
double    BETA41 , BETA42 , BETA4  ;
double    P;


    Seta=(T+273.15)/(TC_Water+273.15);
    CHI=V/VC_Water;
    P=TV2PREG3(T,V);
    BETA3=P/PC_Water;

    Y=(1-Seta)/(1-62315.0/64730.0);
    Y3=Y*Y*Y;
    Y31=pow(Y,31);

    BETA41=(((4*DREG4[10]/CHI+3*DREG4[9])/CHI+2*DREG4[8])/CHI+DREG4[7])/CHI/CHI *Y;
    BETA42=(((4*DREG4[5]/CHI+3*DREG4[4])/CHI+2*DREG4[3])/CHI+DREG4[2])/CHI/CHI;
    BETA4=BETA3+(BETA41+BETA42)*Y3-(2*DREG4[13]*CHI+DREG4[12])*Y31*Y;

    return BETA4*PC_Water;
}


double TV2SREG4 (double T, double V )
{
double    Seta , CHI ;
double    SEGMA3;
double    Y ,Y31,S ;
double    SEGMA41 , SEGMA42 , SEGMA4  ;


    Seta=(T+273.15)/(TC_Water+273.15);
    CHI=V/VC_Water;
    S=TV2SREG3(T,V);
    SEGMA3=S/(PC_Water*1E3*VC_Water/(TC_Water+273.15));

    Y=(1-Seta)/(1-62315.0/64730.0);
    Y31=pow(Y,31);

    SEGMA41=4*((((DREG4[10]/CHI+DREG4[9])/CHI+DREG4[8])/CHI+DREG4[7])/CHI+DREG4[6])*Y;
    SEGMA42=3*((((DREG4[5]/CHI+DREG4[4])/CHI+DREG4[3])/CHI+DREG4[2])/CHI+DREG4[1]);
    SEGMA4=SEGMA3+((SEGMA41+SEGMA42)*Y*Y+32*((DREG4[13]*CHI+DREG4[12])*CHI+DREG4[11])*Y31)/(1-6.2315/6.473);

    return SEGMA4*(PC_Water*1E3*VC_Water/(TC_Water+273.15));
}



double TV2HREG4 (double T, double V )
{
double    Seta , CHI , BETA3;
double    SEGMA3 , EPSILON3;
double    Y , Y3 , Y31 ;
double    PSID1 , PSID2 , PSID  ;
double    BETA41 , BETA42 , BETA4  ;
double    SEGMA41 , SEGMA42 , SEGMA4  ;
double    EPSILON4  ;
double    *P,*S,*H,*X, p,s,h,x;

    P=&p;
	S=&s;
	H=&h;
	X=&x;
    Seta=(T+273.15)/(TC_Water+273.15);
    CHI=V/VC_Water;
    TVREG3(P,T,H,S,V,X);
    BETA3=*P/PC_Water;
    SEGMA3=*S/(PC_Water*1E3*VC_Water/(TC_Water+273.15));
    EPSILON3=*H/(PC_Water*1E3*VC_Water);

    Y=(1-Seta)/(1-62315.0/64730.0);
    Y3=Y*Y*Y;
    Y31=pow(Y,31);


    PSID1=((((DREG4[10]/CHI+DREG4[9])/CHI+DREG4[8])/CHI+DREG4[7])/CHI+DREG4[6])*Y;
    PSID2=(((DREG4[5]/CHI+DREG4[4])/CHI+DREG4[3])/CHI+DREG4[2])/CHI;
    PSID=(PSID1+PSID2+DREG4[1])*Y3+((DREG4[13]*CHI+DREG4[12])*CHI+DREG4[11])*Y31 *Y;


    BETA41=(((4*DREG4[10]/CHI+3*DREG4[9])/CHI+2*DREG4[8])/CHI+DREG4[7])/CHI/CHI *Y;
    BETA42=(((4*DREG4[5]/CHI+3*DREG4[4])/CHI+2*DREG4[3])/CHI+DREG4[2])/CHI/CHI;
    BETA4=BETA3+(BETA41+BETA42)*Y3-(2*DREG4[13]*CHI+DREG4[12])*Y31*Y;


    SEGMA41=4*((((DREG4[10]/CHI+DREG4[9])/CHI+DREG4[8])/CHI+DREG4[7])/CHI+DREG4[6])*Y;
    SEGMA42=3*((((DREG4[5]/CHI+DREG4[4])/CHI+DREG4[3])/CHI+DREG4[2])/CHI+DREG4[1]);
    SEGMA4=SEGMA3+((SEGMA41+SEGMA42)*Y*Y+32*((DREG4[13]*CHI+DREG4[12])*CHI+DREG4[11])*Y31)/(1-6.2315/6.473);


    EPSILON4=EPSILON3+PSID+(SEGMA4-SEGMA3)*Seta+(BETA4-BETA3)*CHI;

    return EPSILON4*(PC_Water*1E3*VC_Water);
}



void  TVREG4 (double *P, double T, double *H, double *S, double V, double * X)
{
double    Seta , CHI , BETA3;
double    SEGMA3 , EPSILON3;
double    Y , Y3 , Y31 ;
double    PSID1 , PSID2 , PSID  ;
double    BETA41 , BETA42 , BETA4  ;
double    SEGMA41 , SEGMA42 , SEGMA4  ;
double    EPSILON4  ;


    Seta=(T+273.15)/(TC_Water+273.15);
    CHI=V/VC_Water;
    TVREG3(P,T,H,S,V,X);
    BETA3=*P/PC_Water;
    SEGMA3=*S/(PC_Water*1E3*VC_Water/(TC_Water+273.15));
    EPSILON3=*H/(PC_Water*1E3*VC_Water);

    Y=(1-Seta)/(1-62315.0/64730.0);
    Y3=Y*Y*Y;
    Y31=pow(Y,31);


    PSID1=((((DREG4[10]/CHI+DREG4[9])/CHI+DREG4[8])/CHI+DREG4[7])/CHI+DREG4[6])*Y;
    PSID2=(((DREG4[5]/CHI+DREG4[4])/CHI+DREG4[3])/CHI+DREG4[2])/CHI;
    PSID=(PSID1+PSID2+DREG4[1])*Y3+((DREG4[13]*CHI+DREG4[12])*CHI+DREG4[11])*Y31 *Y;


    BETA41=(((4*DREG4[10]/CHI+3*DREG4[9])/CHI+2*DREG4[8])/CHI+DREG4[7])/CHI/CHI *Y;
    BETA42=(((4*DREG4[5]/CHI+3*DREG4[4])/CHI+2*DREG4[3])/CHI+DREG4[2])/CHI/CHI;
    BETA4=BETA3+(BETA41+BETA42)*Y3-(2*DREG4[13]*CHI+DREG4[12])*Y31*Y;


    SEGMA41=4*((((DREG4[10]/CHI+DREG4[9])/CHI+DREG4[8])/CHI+DREG4[7])/CHI+DREG4[6])*Y;
    SEGMA42=3*((((DREG4[5]/CHI+DREG4[4])/CHI+DREG4[3])/CHI+DREG4[2])/CHI+DREG4[1]);
    SEGMA4=SEGMA3+((SEGMA41+SEGMA42)*Y*Y+32*((DREG4[13]*CHI+DREG4[12])*CHI+DREG4[11])*Y31)/(1-6.2315/6.473);


    EPSILON4=EPSILON3+PSID+(SEGMA4-SEGMA3)*Seta+(BETA4-BETA3)*CHI;

    *P=BETA4*PC_Water;
    *S=SEGMA4*(PC_Water*1E3*VC_Water/(TC_Water+273.15));
    *H=EPSILON4*(PC_Water*1E3*VC_Water);

    *X=0.0;
}


double P2T_KK(double P)
{
double   Temp_MIN,Temp_MAX,TOL;

    if ( P<0.0 )
    {
        return 0.0;        
    }

    Temp_MIN=T000C;
    Temp_MAX=374.15; //(float)TC_Water;
    TOL=(1E-8)*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, P, 5606);

}

///////T2
double T2HBOUND23(double T)
{
    return PT2HREG2(T2PBOUND23(T),T);
}

double T2SBOUND23(double T)
{
    return PT2SREG2(T2PBOUND23(T),T);
}

double T2VBOUND23(double T)
{
    return PT2VREG2(T2PBOUND23(T),T);
}


///////P2
double P2TBOUND23(double P)
{
double    Temp_MIN,Temp_MAX,TOL;

    Temp_MIN=T350C;
    Temp_MAX=590.0;
    TOL=(1E-8)*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, P, 3200);
}

double P2HBOUND23(double P)
{
    return PT2HREG2(P,P2TBOUND23(P));
}

double P2SBOUND23(double P)
{
    return PT2SREG2(P,P2TBOUND23(P));
}

double P2VBOUND23(double P)
{
    return PT2VREG2(P,P2TBOUND23(P));
}

//////H2
double H2TLPBOUND23(double  H)  
{
double    Temp_MIN , Temp_MAX , TOL;

    if (H<2610.68079 )
    {
        Temp_MIN = 350;
        Temp_MAX = 365.12;
    }
    else if (H>2629.04908 )
    {
        Temp_MIN = 490.32;
        Temp_MAX = 590.0;
    }
    else
    {
        Temp_MIN = 365.11;
        Temp_MAX = 393.20;
    }
    TOL = 1E-8*(Temp_MIN + Temp_MAX) / 2.0;
    return  ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, H, 3201);
}



double H2TMPBOUND23(double  H) 
{
double    Temp_MIN , Temp_MAX , TOL ;


    if (H<2610.68079) 
    {
        Temp_MIN = 350;
        Temp_MAX = 365.12;
    }
    else if (H>2629.04908) 
    {
        Temp_MIN = 490.32;
        Temp_MAX = 590.0;
    }
    else
    {
        Temp_MIN = 393.20;
        Temp_MAX = 454.33;
    }
    TOL = 1E-8*(Temp_MIN + Temp_MAX) / 2.0;

    return  ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, H, 3201);
}


double H2THPBOUND23(double H)  
{
double    Temp_MIN , Temp_MAX , TOL ;


    if (H<2610.68079) 
    {
        Temp_MIN = 350;
        Temp_MAX = 365.12;
    }
    else if (H>2629.04908) 
    {
        Temp_MIN = 490.32;
        Temp_MAX = 590.0;
    }
    else
    {
        Temp_MIN = 454.33;
        Temp_MAX = 490.33;
    }
    TOL = 1E-8*(Temp_MIN + Temp_MAX) / 2.0;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, H, 3201);
}


double H2TBOUND23(double  H)  
{
    return H2THPBOUND23(H);
}


double H2PLPBOUND23(double  H) 
{
    return T2PBOUND23(H2TLPBOUND23(H));
}


double H2PMPBOUND23(double  H)  
{
    return T2PBOUND23(H2TMPBOUND23(H));
}


double H2PHPBOUND23(double H)  
{
    return T2PBOUND23(H2THPBOUND23(H));
}


double H2PBOUND23(double H)  
{
    return H2PHPBOUND23(H);
}

double H2SLPBOUND23(double H )  
{
double    T;

    T=H2TLPBOUND23(H);
    return PT2SREG2(T2PBOUND23(T),T );
}

double H2SMPBOUND23(double  H )  
{
double    T;    

    T=H2TMPBOUND23(H);
    return PT2SREG2(T2PBOUND23(T),T );
}


double H2SHPBOUND23(double  H )  
{
double    T;

    T=H2THPBOUND23(H);
    return PT2SREG2(T2PBOUND23(T),T );
}

double H2SBOUND23(double  H )  
{
    return H2SHPBOUND23(H);
}

/////////V
double V2TBOUND23(double  V )  
{
double    Temp_MIN , Temp_MAX , TOL;

    Temp_MIN = 350;
    Temp_MAX = 590;

    TOL = (1E-8)*(Temp_MIN + Temp_MAX) / 2.0;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, V, 3202);
}

double V2PBOUND23(double V)  
{
    return T2PBOUND23(V2TBOUND23(V));
}

double V2HBOUND23(double V)
{
    return T2HBOUND23(V2TBOUND23(V));
}


double V2SBOUND23(double V)
{
    return T2SBOUND23(V2TBOUND23(V));
}

//////////P2
double P2HLREG56(double  P)
{
double    T;

    T=P2T_KK(P);
    if (fabs(1-T/TC_Water)<1E-8) 
        return 107.4;
    else if (T>=T350C+DeltaVal )
        return PT2HREG4(P,T);
    else
        return PT2HREG1(P,T);
}

double P2SLREG56(double  P)
{
double    T;

    T=P2T_KK(P);
    if (fabs(1-T/TC_Water)<1E-8 )
        return 4.4429;
    else if (T>=T350C+DeltaVal )
        return PT2SREG4(P,T);
    else
        return PT2SREG1(P,T) ;
}

double P2VLREG56(double P)
{
double    T;

    T=P2T_KK(P);
    if (fabs(1-T/TC_Water)<1E-8) 
        return 0.00317;
    else if (T>=T350C+DeltaVal) 
        return PT2VREG4NT(P,T);
    else
        return PT2VREG1(P,T);
}

double P2HGREG56(double P)
{
double    T;

    T=P2T_KK(P);
    if (fabs(1-T/TC_Water)<1E-8) 
        return 2107.4;
    else if (T>=T350C+DeltaVal) 
        return PT2HREG3(P,T);
    else
        return PT2HREG2(P,T) ;
}

double P2SGREG56(double  P)
{
double    T;

    T=P2T_KK(P);
    if (fabs(1-T/TC_Water)<1E-8) 
        return 4.4429;
    else if (T>=T350C+DeltaVal) 
        return PT2SREG3(P,T);
    else
        return PT2SREG2(P,T) ;
}

double P2VGREG56(double P)
{
double    T;

    T=P2T_KK(P);
    if (fabs(1-T/TC_Water)<1E-8 )
        return 0.00317;
    else if (T>=T350C+DeltaVal )
        return PT2VREG3NT(P,T);
    else
        return PT2VREG2(P,T);
}

////////T2XLREG56
double T2HLREG56(double T)  
{
    if (fabs(1-T/TC_Water)<1E-8 )
        return 2107.4;
    else if (T>=T350C+DeltaVal) 
    {
       return PT2HREG4(T2P_KK(T),T);
    }
    else
        return PT2HREG1(T2P_KK(T),T) ;
}

double T2SLREG56(double T)  
{
    if (fabs(1-T/TC_Water)<1E-8 )
        return 4.4429;
    else if (T>=T350C+DeltaVal )
        return PT2SREG4(T2P_KK(T),T );
    else
        return PT2SREG1(T2P_KK(T),T );
}

double PT2VREG4NT(double P, double T)
{
double    ZVMAX , ZVMIN , ZF ;
double    VTEMP , PTEMP ;
double    ZPCT , ZDV , ZPCTP ;
double   ZVPREV , ZPPREV ;
int     I ;

    if ((fabs(1-P/PC_Water)<1.0E-6) && (fabs(1-T/TC_Water)<1.0E-6) )
    {
        return VC_Water;   
		//exit(1);
    }

    if ( (fabs(1-P/P350C)<1.0E-6) && (fabs(1-T/T350C)<1.0E-6) )
    {
        return 0.0017411; 
		//exit(1);
    }

    ZPCTP =0.0;
    ZVPREV=0.0;
    ZPPREV=0.0;
    ZVMAX = 0.0089272;
    ZVMIN = 0.00129226;
    ZF = 0.25;

    VTEMP = (0.004807 + 0.3994 /(10*P)) * (T - 250.0) / (T + 50.0);

    for (I=1;I<=100;I++)
    {
        if (VTEMP > ZVMAX) 
        {
            VTEMP = ZVMAX;
            ZF = 1.0;
        }
        else 
            if (VTEMP < ZVMIN) 
            {
                VTEMP = ZVMIN;
                ZF = 1.0;
            }

        PTEMP=TV2PREG4(T,VTEMP)*10.0; 
        ZPCT = PTEMP /(10*P) - 1.0;

        if (fabs(ZPCT) < 0.0000001) 
        {
            return VTEMP;
            //exit(1);
        }

        if (I < 2) 
		{
			ZDV = VTEMP * ZPCT * ZF;
		}
		else
        {
            if (ZPCTP * ZPCT > 0.0) 
			{
                if (fabs(ZPCT) <= fabs(0.3 * ZPCTP)) 
				{
                    ZDV = (VTEMP - ZVPREV) * ((10*P) - PTEMP) / (PTEMP - ZPPREV);
                }
				else
                {
                    ZF = 1.5 * ZF;
                    ZDV = VTEMP * ZPCT * ZF;
                }
			}
            else
            {
                ZF = 0.67 * ZF;
                ZDV = (VTEMP - ZVPREV) * ((10*P) - PTEMP) / (PTEMP - ZPPREV);
            }

        ZVPREV = VTEMP;
        ZPPREV = PTEMP;         
        ZPCTP = ZPCT;
        VTEMP = VTEMP + ZDV;
		}
	}
    return VTEMP;
}


double T2VLREG56(double  T)  
{
    if (T>=(T350C+DeltaVal)) 
    {
        return PT2VREG4NT(T2P_KK(T),T);
    }
    else
        return PT2VREG1(T2P_KK(T),T);
}

double T2HGREG56(double  T) 
{
    if (fabs(1-T/TC_Water)<1E-8) 
        return 2107.4;
    else 
	if (T>=(T350C+DeltaVal)) 
        return PT2HREG3(T2P_KK(T),T);
    else
        return PT2HREG2(T2P_KK(T),T) ;
}


double T2SGREG56(double  T)  
{
    if (fabs(1-T/TC_Water)<1E-8) 
        return 4.4429;
    else if (T>=(T350C+DeltaVal)) 
        return PT2SREG3(T2P_KK(T),T);
    else
        return PT2SREG2(T2P_KK(T),T);
}

double PT2VREG3NT(double P, double T)
{
double    ZVMAX , ZVMIN , ZF ;
double    VTEMP , PTEMP ;
double    ZPCT , ZDV , ZPCTP ;
double    ZVPREV , ZPPREV ;
int     I ;

    if ((fabs(1-P/PC_Water)<1.0E-6) && (fabs(1-T/TC_Water)<1.0E-6) )
    {
        return 0.00317;  
		//exit(1);
	}

    if ((fabs(1-P/P350C)<1.0E-6) && (fabs(1-T/T350C)<1.0E-6) )
    {
        return 0.008799;    
		//exit(1);
    }

    ZPCTP =0.0;
    ZVPREV=0.0;
    ZPPREV=0.0;
    ZVMAX = 0.0089272;
    ZVMIN = 0.00129226;
    ZF = 0.25;

    VTEMP = 0.272 * (T + 273.15) /(10*P);

    
	for (I=1;I<=100;I++)
    {
        if (VTEMP > ZVMAX) 
        {
            VTEMP = ZVMAX;
            ZF = 1.0;
        }
        else
            if (VTEMP < ZVMIN) 
            {
                VTEMP = ZVMIN;
                ZF = 1.0;
            }

        PTEMP=TV2PREG3(T,VTEMP)*10;
        ZPCT = PTEMP /(10*P) - 1.0;

        if (fabs(ZPCT) < 0.0000001) 
        {
            return VTEMP; 
			//exit(1);
        }

        if (I < 2) 
		{
			ZDV = VTEMP * ZPCT * ZF;
		}
        else
		{
            if (ZPCTP * ZPCT > 0.0) 
                if (fabs(ZPCT) <= fabs(0.3 * ZPCTP))                     
				{
					ZDV = (VTEMP - ZVPREV) * ((10*P) - PTEMP) / (PTEMP - ZPPREV);
				}
                else
                {
                    ZF = 1.5 * ZF;
                    ZDV = VTEMP * ZPCT * ZF;
                }
            else
            {
                ZF = 0.67 * ZF;
                ZDV = (VTEMP - ZVPREV) * ((10*P) - PTEMP) / (PTEMP - ZPPREV);
            }

        ZVPREV = VTEMP;
        ZPPREV = PTEMP;
        ZPCTP = ZPCT;
        VTEMP = VTEMP + ZDV;
		}
	}  
    return VTEMP;
}



double T2VGREG56(double T)  
{
    if (T>=(T350C+DeltaVal)) 
    {
        return PT2VREG3NT(T2P_KK(T),T);
    }
    else
        return PT2VREG2(T2P_KK(T),T);
}

///////////////H2XLREG56
double H2TLREG56(double H)
{
double    Temp_MIN,Temp_MAX,TOL;

    Temp_MIN=T000C;
    Temp_MAX=TC_Water;
    TOL=1E-8*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, H, 5600);
}


double H2PLREG56(double H)
{
    return T2P_KK(H2TLREG56(H));
}


double H2SLREG56(double H)
{
double    T;

    T=H2TLREG56(H);
    if (T>=(T350C+DeltaVal)) 
        return PT2SREG3(T2P_KK(T),T);
    else
        return PT2SREG2(T2P_KK(T),T);
}


double H2VLREG56(double H)
{
double    T;

    T=H2TLREG56(H);
    if (T>=(T350C+DeltaVal)) 
        return PT2VREG3NT(T2P_KK(T),T);
    else
        return PT2VREG2(T2P_KK(T),T);
}


double H2TGLREG56(double H)
{
double    Temp_MIN,Temp_MAX,TOL,*T,t;

    T=&t;
    Temp_MIN=T000C;
    Temp_MAX=TC_Water;
    T_HMAXXREG56(Temp_MIN,Temp_MAX,1.0,T);
    Temp_MAX=*T;
    TOL=1E-8*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, H, 5601);
}

double H2TGHREG56(double H )
{
double    Temp_MIN,Temp_MAX,TOL,*T,t;

    T=&t;
    Temp_MIN=T000C;
    Temp_MAX=TC_Water;
    T_HMAXXREG56(Temp_MIN,Temp_MAX,1.0,T);
    Temp_MIN=*T;
    TOL=1E-8*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, H, 5601);
}

double H2PGLREG56(double H)
{
    return T2P_KK(H2TGLREG56(H));
}

double H2PGHREG56(double H )
{
    return T2P_KK(H2TGHREG56(H));
}

double H2SGLREG56(double H)
{
double    T;

    T=H2TGLREG56(H);
    if (T>=(T350C+DeltaVal)) 
        return PT2SREG3(T2P_KK(T),T);
    else
        return PT2SREG2(T2P_KK(T),T);
}

double H2SGHREG56(double H)
{
double    T;

    T=H2TGHREG56(H);
    if (T>=(T350C+DeltaVal)) 
        return PT2SREG3(T2P_KK(T),T);
    else
        return PT2SREG2(T2P_KK(T),T);
}

double H2VGLREG56(double H)
{
double    T;

    T=H2TGLREG56(H);
    if (T>=(T350C+DeltaVal)) 
        return PT2VREG3NT(T2P_KK(T),T);
    else
        return PT2VREG2(T2P_KK(T),T);
}

double H2VGHREG56(double H)
{
double    T;

    T=H2TGHREG56(H);
    if (T>=T350C+DeltaVal) 
        return PT2VREG3NT(T2P_KK(T),T);
    else
        return PT2VREG2(T2P_KK(T),T);
}

double S2TLREG56(double S)
{
double    Temp_MIN,Temp_MAX,TOL;

    Temp_MIN=T000C;
    Temp_MAX=TC_Water;
    TOL=1E-8*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, S, 5602);
}

double S2PLREG56(double S)
{
    return T2P_KK(S2TLREG56(S));
}


double S2HLREG56(double S)
{
double    T;

    T=S2TLREG56(S);
    if (T>=(T350C+DeltaVal)) 
        return PT2HREG4(T2P_KK(T),T);
    else
        return PT2HREG1(T2P_KK(T),T);
}

double S2VLREG56(double S)
{
double    T;

    T=S2TLREG56(S);
    if (T>=(T350C+DeltaVal)) 
       return PT2VREG4NT(T2P_KK(T),T);
    else
       return PT2VREG1(T2P_KK(T),T);
}

double S2TGREG56(double S)
{
double    Temp_MIN,Temp_MAX,TOL;

    Temp_MIN=T000C;
    Temp_MAX=TC_Water;
    TOL=1E-8*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, S, 5603);
}

double S2PGREG56(double S)
{
    return T2P_KK(S2TGREG56(S));
}


double S2HGREG56(double S)
{
double    T;

    T=S2TGREG56(S);
    if (T>=(T350C+DeltaVal)) 
       return PT2HREG3(T2P_KK(T),T);
    else
       return PT2HREG2(T2P_KK(T),T);
}

double S2VGREG56(double S)
{
double    T;

    T=S2TGREG56(S);
    if (T>=(T350C+DeltaVal ))
        return PT2VREG3NT(T2P_KK(T),T);
    else
        return PT2VREG2(T2P_KK(T),T);
}


/////////////V2XLREG56
double V2TLREG56(double V)
{
double    Temp_MIN,Temp_MAX,TOL;

    Temp_MIN=T000C;
    Temp_MAX=TC_Water;
    TOL=1E-8*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, V, 5604);
}

double V2PLREG56(double V)
{
    return T2P_KK(V2TLREG56(V));
}

double V2HLREG56(double V)
{
double    T;

    T=V2TLREG56(V);
    if (T>=(T350C+DeltaVal)) 
        return PT2HREG4(T2P_KK(T),T);
    else
        return PT2HREG1(T2P_KK(T),T);
}

double V2SLREG56(double V)
{
double    T;

    T=V2TLREG56(V);
    if (T>=(T350C+DeltaVal)) 
        return PT2SREG4(T2P_KK(T),T);
    else
        return PT2SREG1(T2P_KK(T),T);
}

double V2TGREG56(double V)
{
double    Temp_MIN,Temp_MAX,TOL;

    Temp_MIN=T000C;
    Temp_MAX=TC_Water;
    TOL=1E-8*(Temp_MIN+Temp_MAX)/2;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, V, 5605);
}

double V2PGREG56(double V)
{
    return T2P_KK(V2TGREG56(V));
}

double V2HGREG56(double V)
{
double    T;

    T=V2TGREG56(V);
    if (T>=(T350C+DeltaVal)) 
        return PT2HREG3(T2P_KK(T),T);
    else
        return PT2HREG2(T2P_KK(T),T);
}

double V2SGREG56(double V)
{
double    T;

    T=V2TGREG56(V);
    if (T>=(T350C+DeltaVal)) 
        return PT2SREG3(T2P_KK(T),T);
    else
        return PT2SREG2(T2P_KK(T),T);
}

/////////PX2XREG56
double PX2HREG56(double  P, double X ) 
{
double    HL , HG ;

    HL = P2HLREG56(P);
    HG = P2HGREG56(P);
    return HL + (HG - HL) * X ;
}

double PX2SREG56(double P ,double X )  
{
double    SL , SG ;

    SL = P2SLREG56(P);
    SG = P2SGREG56(P);
    return SL + (SG - SL) * X;
}

double PX2VREG56(double  P ,double X )  
{
double    VL , VG ;

    VL = P2VLREG56(P);
    VG = P2VGREG56(P);
    return VL + (VG - VL) * X ;
}

void PXREG56(double P ,double * T,double *H,double *S,double *V,double X)
{
    *T=P2T_KK(P);
    *H=PX2HREG56(P,X);
    *S=PX2SREG56(P,X);
    *V=PX2VREG56(P,X);
}

//////////////P2XREG56
double PH2XREG56(double  Pressure , double Enthalpy ) 
{
double    HL , HG ;

    HL = P2HLREG56(Pressure);
    HG = P2HGREG56(Pressure);
    return (Enthalpy - HL) / (HG - HL);
}

double PH2SREG56(double  Pressure , double Enthalpy )  
{
double    X ;

    X = PH2XREG56(Pressure, Enthalpy);
    return PX2SREG56(Pressure, X);
}


double PH2VREG56(double  Pressure , double Enthalpy )  
{
double    X ;

    X = PH2XREG56(Pressure, Enthalpy);
    return PX2VREG56(Pressure, X);
}

void  PHREG56(double P,double * T, double H, double * S,double *V,double *X)
{
    *T=P2T_KK(P);
    *X=PH2XREG56(P,H);
    *S=PX2SREG56(P,*X);
    *V=PX2VREG56(P,*X);
}

////////////PS   REG56
double PS2XREG56(double  Pressure ,double Entropy )  
{
double    SL , SG ;

    SL = P2SLREG56(Pressure);
    SG = P2SGREG56(Pressure);
    return (Entropy - SL) / (SG - SL);
}


double PS2HREG56(double  Pressure , double Entropy )  
{
double    X ;

    X = PS2XREG56(Pressure, Entropy);
    return PX2HREG56(Pressure, X);
}



double PS2VREG56(double  Pressure , double Entropy )  
{
double    X ;

    X = PS2XREG56(Pressure, Entropy);
    return PX2VREG56(Pressure, X);
}


void  PSREG56(double P, double * T,double *H, double S,double*  V,double* X)
{
    *T=P2T_KK(P);
    *X=PS2XREG56(P,S);
    *H=PX2HREG56(P,*X);
    *V=PX2VREG56(P,*X);
}


////////////PV  REG56
double PV2XREG56(double  Pressure ,double Volume )  
{
double    VL , VG ;

    VL = P2VLREG56(Pressure);
    VG = P2VGREG56(Pressure);
    return (Volume - VL) / (VG - VL);
}


double PV2HREG56(double  Pressure ,double Volume )  
{
double    X ;

    X = PV2XREG56(Pressure, Volume);
    return PX2HREG56(Pressure, X);
}


double PV2SREG56(double  Pressure ,double Volume )  
{
double    X ;

    X = PV2XREG56(Pressure, Volume);
    return PX2SREG56(Pressure, X);
}

void  PVREG56(double P, double * T,double *H,double *S, double V, double* X)
{
    *T=P2T_KK(P);
    *X=PV2XREG56(P,V);
    *H=PX2HREG56(P,*X);
    *S=PX2SREG56(P,*X);
}

///////TX2  REG56
double TX2HREG56(double  T,double X)  
{
double    HL , HG ;

    HL = T2HLREG56(T);
    HG = T2HGREG56(T);
    return HL + (HG - HL) * X;
}


double TX2SREG56(double  T ,double X )  
{
double    SL , SG ;

    SL = T2SLREG56(T);
    SG = T2SGREG56(T);
    return SL + (SG - SL) * X;
}



double TX2VREG56(double  T ,  double X )  
{
double    VL , VG ;

    VL = T2VLREG56(T);
    VG = T2VGREG56(T);
    return VL + (VG - VL) * X;
}


void  TXREG56(double* P,double T, double * H,double *S,double *V, double X)
{
    *P=T2P_KK(T);
    *H=TX2HREG56(T,X);
    *S=TX2SREG56(T,X);
    *V=TX2VREG56(T,X);
}

//////////////////
void  T_HMAXXREG56(double Temp_MIN,double Temp_MAX,double X ,double * T)
{
double    Temp_Mid  ;
double    Temp_Mid_Plus , H_Mid_Plus ;
double    Temp_Mid_MINus , H_Mid_MINus ;
double    Temp_MIN2 , Temp_MAX2;

    if (fabs(Temp_MAX - Temp_MIN) < 0.02 )
    {
        *T = (Temp_MAX + Temp_MIN) / 2;
        return ;
    }

    Temp_Mid = (Temp_MAX + Temp_MIN) / 2;

    Temp_Mid_Plus = Temp_Mid + 0.01;
    H_Mid_Plus = TX2HREG56(Temp_Mid_Plus, X);

    Temp_Mid_MINus = Temp_Mid - 0.01;
    H_Mid_MINus = TX2HREG56(Temp_Mid_MINus, X);


    Temp_MIN2=Temp_MIN;
    Temp_MAX2=Temp_MAX;

    if (H_Mid_Plus > H_Mid_MINus )
        Temp_MIN2 = Temp_Mid;
    else
        Temp_MAX2 = Temp_Mid;

    T_HMAXXREG56(Temp_MIN2, Temp_MAX2, X, T);
}


double HMAXX2TREG56(double X)
{
double    Temp_MIN , Temp_MAX ;
double    *T,t ;

    T=&t;
    if (X <= 0.442) 
    {
        return TC_Water;
        
    }
    Temp_MIN = T000C;
    Temp_MAX = TC_Water;
    T_HMAXXREG56(Temp_MIN,Temp_MAX,X,T);
    return *T;
}


double HMAXXREG56(double X)
{
    return TX2HREG56(HMAXX2TREG56(X),X);
}


//////////////TH  REG56
double TH2XREG56(double  Temperature , double Enthalpy )  
{
double    HL , HG ;

    HL = T2HLREG56(Temperature);
    HG = T2HGREG56(Temperature);
    return (Enthalpy - HL) / (HG - HL);
}


double TH2SREG56(double  Temperature ,double  Enthalpy )  
{
double    X ;

    X = TH2XREG56(Temperature, Enthalpy);
    return TX2SREG56(Temperature, X);
}



double TH2VREG56(double  Temperature , double Enthalpy )  
{
double    X ;

    X = TH2XREG56(Temperature, Enthalpy);
    return TX2VREG56(Temperature, X);
}


void  THREG56(double*  P, double T,double H,double*  S,double *V,double *X)
{
    *P=T2P_KK(T);
    *X=TH2XREG56(T,H);
    *S=TX2SREG56(T,*X);
    *V=TX2VREG56(T,*X);
}

/////////////TS   REG56
double TS2XREG56(double  Temperature ,double Entropy )  
{
double    SL , SG;

    SL = T2SLREG56(Temperature);
    SG = T2SGREG56(Temperature);
    return (Entropy - SL) / (SG - SL);
}


double TS2HREG56(double  Temperature ,double Entropy ) 
{
double    X ;

    X = TS2XREG56(Temperature, Entropy);
    return TX2HREG56(Temperature, X);
}



double TS2VREG56(double  Temperature , double Entropy )  
{
double    X ;

    X = TS2XREG56(Temperature, Entropy);
    return TX2VREG56(Temperature, X);
}

void TSREG56(double * P,double T,double * H, double S, double *V,double *X)
{
    *P=T2P_KK(T);
    *X=TS2XREG56(T,S);
    *H=TX2HREG56(T,*X);
    *V=TX2VREG56(T,*X);
}


////////TV  REG56
double TV2XREG56(double  Temperature , double Volume) 
{
double    VL , VG ;

    VL = T2VLREG56(Temperature);
    VG = T2VGREG56(Temperature);
    return (Volume - VL) / (VG - VL);
}


double TV2HREG56(double  Temperature , double Volume )  
{
double    X ;

    X = TV2XREG56(Temperature, Volume);
    return TX2HREG56(Temperature, X);
}



double TV2SREG56(double  Temperature , double Volume )  
{
double    X ;

    X = TV2XREG56(Temperature, Volume);
    return TX2SREG56(Temperature, X);
}

void TVREG56(double * P, double T,double * H,double *S, double V,double * X)
{
    *P=T2P_KK(T);
    *X=TV2XREG56(T,V);
    *H=TX2HREG56(T,*X);
    *S=TX2SREG56(T,*X);
}

////////HS REG56
double HS2PREG56(double  Enthalpy , double Entropy )  
{
double    Press_MIN , Press_MAX , TOL ;

    Press_MIN = 0.0006108;
    if (Entropy <= SC_Water )
        Press_MAX = S2PLREG56(Entropy);
    else
        Press_MAX = S2PGREG56(Entropy);

    TOL = 1E-8 * (Press_MIN + Press_MAX) / 2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, Entropy, Enthalpy, 5600);
}


double HS2TREG56(double  Enthalpy ,double Entropy ) 
{
double    Pressure ;

    Pressure = HS2PREG56(Enthalpy, Entropy);
    return P2T_KK(Pressure);
}


double HS2VREG56(double  Enthalpy ,double Entropy )  
{
double    Pressure ;

    Pressure = HS2PREG56(Enthalpy, Entropy);
    return PH2VREG56(Pressure, Enthalpy);
}


double HS2XREG56(double  Enthalpy , double Entropy ) 
{
double    Pressure ;

    Pressure = HS2PREG56(Enthalpy, Entropy);
    return PS2XREG56(Pressure, Entropy);
}

void HSREG56(double*P,double *T,double  H,double S, double* V,double *X)
{
    *P=HS2PREG56(H,S);
    *T=P2T_KK(*P);
    *V=PH2VREG56(*P,H);
    *X=PS2XREG56(*P,S);
}

/////////HV REG56
double HV2PREG56(double  Enthalpy ,double Volume )  
{
double    Press_MIN , Press_MAX , TOL ;


    Press_MIN = 0.0006108;

    if (Volume < VC_Water) 
        Press_MAX = V2PLREG56(Volume);
    else
        Press_MAX = V2PGREG56(Volume);


    TOL = 1E-8 * (Press_MIN + Press_MAX) / 2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, Volume, Enthalpy, 5601);
}



double HV2TREG56(double  Enthalpy , double Volume )  
{
double    Pressure ;

    Pressure = HV2PREG56(Enthalpy, Volume);
    return P2T_KK(Pressure);
}


double HV2SREG56(double  Enthalpy , double Volume )  
{
double    Pressure ;

    Pressure = HV2PREG56(Enthalpy, Volume);
    return PH2SREG56(Pressure, Enthalpy);
}


double HV2XREG56(double  Enthalpy ,double Volume ) 
{
double    Pressure ;

    Pressure = HV2PREG56(Enthalpy, Volume);
    return PH2XREG56(Pressure, Enthalpy);
}


void HVREG56(double * P,double *T, double H, double *S, double V, double * X)
{
    *P=HV2PREG56(H,V);
    *T=P2T_KK(*P);
    *S=PH2SREG56(*P,H);
    *X=PH2XREG56(*P,H);
}



void GETTRANGEHX(double* Temp_MIN , double *Temp_MAX , double Enthalpy ,double X )
{
double    FMIN , FMid , FMidMIN;
double    Temp_Midmin , Temp_Mid;

    if (fabs(*Temp_MIN - *Temp_MAX) < 0.00001)
		return;

    FMIN =TX2HREG56(*Temp_MIN, X) - Enthalpy;

    Temp_Mid = (*Temp_MIN + *Temp_MAX) / 2;
    FMid =TX2HREG56(Temp_Mid, X) - Enthalpy;

    if (FMIN * FMid < 0 )
        *Temp_MAX = Temp_Mid;
    else
    {
        Temp_Midmin = (*Temp_MIN + Temp_Mid) / 2;
        FMidMIN =TX2HREG56(Temp_Midmin, X) - Enthalpy;

        if (FMid * FMidMIN < 0 )
        {
            *Temp_MIN = Temp_Midmin;
            *Temp_MAX = Temp_Mid;
        }
        else
        {
            *Temp_MIN = Temp_Mid;
            GETTRANGEHX(Temp_MIN, Temp_MAX, Enthalpy, X);
        }
    }
}

///////////HX   REG56
double HX2TREG56(double  Enthalpy , double X )  
{
double    *Temp_MIN , *Temp_MAX , TOL,temp_min,temp_max ;
double    FMIN , FMAX ;

    Temp_MIN=&temp_min;
	Temp_MAX=&temp_max;
    *Temp_MIN = T000C;
    *Temp_MAX = TC_Water;
    
    FMIN =TX2HREG56(*Temp_MIN, X) - Enthalpy;
    FMAX =TX2HREG56(*Temp_MAX, X) - Enthalpy;

    if (FMIN * FMAX > 0) 
        GETTRANGEHX(Temp_MIN, Temp_MAX, Enthalpy, X);

    TOL = 1E-8 * (*Temp_MIN + *Temp_MAX) / 2.0;
    return ZBRENT2_KK(*Temp_MIN, *Temp_MAX, TOL, X, Enthalpy, 5602);
}


double HX2THPREG56(double  Enthalpy , double X )  
{
double    Temp_MIN , Temp_MAX , TOL;
double    FMAX , FMIN ;
double    T_Low ;


    T_Low = HX2TREG56(Enthalpy, X);
    Temp_MAX = TC_Water;
    Temp_MIN = (T_Low + Temp_MAX) / 2;
    
    FMIN = TX2HREG56(Temp_MIN, X) - Enthalpy;
    FMAX = TX2HREG56(Temp_MAX, X) - Enthalpy;
    if (FMIN * FMAX > 0) 
    {
        return T_Low;
        //exit(1);
    }
    
    
    TOL = 1E-8 * (Temp_MIN + Temp_MAX) / 2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, X, Enthalpy, 5602);
}


double HX2PREG56(double  Enthalpy , double X )  
{
double    Temperature;

    Temperature = HX2TREG56(Enthalpy, X);
    return T2P_KK(Temperature);
}


double HX2PLPREG56(double  Enthalpy , double X )  
{
double    Temperature;

    Temperature = HX2TREG56(Enthalpy, X);
    return T2P_KK(Temperature);
}



double HX2PHPREG56(double  Enthalpy ,double X)  
{
double    Temperature;

    Temperature = HX2THPREG56(Enthalpy, X);
    return T2P_KK(Temperature);
}


double HX2SREG56(double  Enthalpy , double X )  
{
double    Temperature;

    Temperature = HX2TREG56(Enthalpy, X);
    return TX2SREG56(Temperature, X);
}


double HX2SHPREG56(double  Enthalpy , double X )  
{
double    Temperature;

    Temperature = HX2THPREG56(Enthalpy, X);
    return TX2SREG56(Temperature, X);
}


double HX2VREG56(double  Enthalpy , double X )  
{
double    Temperature;

    Temperature = HX2TREG56(Enthalpy, X);
    return TX2VREG56(Temperature, X);
}


double HX2VHPREG56(double  Enthalpy , double X )  
{
double    Temperature;

    Temperature = HX2THPREG56(Enthalpy, X);
    return TX2VREG56(Temperature, X);
}

void  HXREG56(double *P,double*T, double H, double * S,double *V, double X)
{
    *T=HX2TREG56(H,X);
    *P=T2P_KK(*T);
    *S=TX2SREG56(*T,X);
    *V=TX2VREG56(*T,X);
}

void HXHPREG56(double * P,double* T, double H, double * S,double *V,double X)
{
    *T=HX2THPREG56(H,X);
    *P=T2P_KK(*T);
    *S=TX2SREG56(*T,X);
    *V=TX2VREG56(*T,X);
}

//////////SV  REG56
double SV2PREG56(double  Entropy , double Volume ) 
{
double    Press_MIN , Press_MAX , TOL ;

    Press_MIN = 0.0006108;

    if (Volume < VC_Water) 
    {
        Press_MAX = V2PLREG56(Volume);
        if (fabs(Entropy)<1E-8) 
        {
            if (fabs(P2SLREG56(Press_MAX)-Entropy)<1E-6) 
            {
                return Press_MAX; 
				//exit(1);
            }
        }
        else
		{
            if (fabs(1-P2SLREG56(Press_MAX)/Entropy)<1E-6) 
            {
                return Press_MAX; 
				//exit(1);
            }
        }

    }
    else
    {
        Press_MAX = V2PGREG56(Volume);
        if (fabs(Entropy)<1E-8 )
        {
            if (fabs(P2SGREG56(Press_MAX)-Entropy)<1E-6 )
            {
                return Press_MAX;  
				//exit(1);
            }
        }
        else
        {
            if (fabs(1-P2SGREG56(Press_MAX)/Entropy)<1E-6 )
            {
                return Press_MAX;  
				//exit(1);
			}
        }
	}


    TOL = 1E-8 * (Press_MIN + Press_MAX) / 2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, Entropy, Volume, 5603);
}

double SV2TREG56(double  Entropy ,double  Volume )  
{
double    Pressure ;

    Pressure = SV2PREG56(Entropy, Volume);
    return P2T_KK(Pressure);
}

double SV2HREG56(double  Entropy , double Volume )  
{
double   Pressure ;

    Pressure = SV2PREG56(Entropy, Volume);
    return PS2HREG56(Pressure, Entropy);
}


double SV2XREG56(double  Entropy ,double Volume )  
{
double    Pressure ;

    Pressure = SV2PREG56(Entropy, Volume);
    return PS2XREG56(Pressure, Entropy);
}

void SVREG56(double * P,double* T,double* H, double S,double V, double * X)
{
    *P=SV2PREG56(S,V);
    *T=P2T_KK(*P);
    *H=PS2HREG56(*P,S);
    *X=PS2XREG56(*P,S);
}



void  T_SMINXREG56(double Temp_MIN,double Temp_MAX,double X, double * T)
{
double    Temp_Mid  ;
double    Temp_Mid_Plus , S_Mid_Plus ;
double    Temp_Mid_MINus , S_Mid_MINus;
double    Temp_MIN2 , Temp_MAX2;

    if (fabs(Temp_MAX-Temp_MIN) < 0.02) 
    {
        *T = (Temp_MAX + Temp_MIN) / 2;
        return;
    }

    Temp_Mid = (Temp_MAX + Temp_MIN) / 2;

    Temp_Mid_Plus = Temp_Mid + 0.01;
    S_Mid_Plus = TX2SREG56(Temp_Mid_Plus, X);

    Temp_Mid_MINus = Temp_Mid - 0.01;
    S_Mid_MINus = TX2SREG56(Temp_Mid_MINus, X);


    Temp_MIN2=Temp_MIN;
    Temp_MAX2=Temp_MAX;

    if (S_Mid_Plus < S_Mid_MINus )
        Temp_MIN2 = Temp_Mid;
    else
        Temp_MAX2 = Temp_Mid;

    T_SMINXREG56(Temp_MIN2, Temp_MAX2, X, T);
}


double SMINX2TREG56(double X)
{
double    Temp_MIN , Temp_MAX ;
double    *T,SMIN,TMIN,STX,t ;
int    IFor;

    T=&t;
    if (X >= 0.55 )
    {
        return TC_Water;
        //exit(1) ;
    }
    SMIN=100;
    TMIN=TC_Water;    
	for (IFor=0; IFor<=375; IFor++)
    {
        *T=T000C+(TC_Water-T000C)/375*IFor;
        if (*T>TC_Water)
			break;
	
        STX=TX2SREG56(*T,X);
        if (STX<SMIN )
        {
            SMIN=STX;
            TMIN=*T;
        }
    }

    Temp_MIN = TMIN-1;
    if (Temp_MIN<T000C)  Temp_MIN=T000C;
    Temp_MAX = TMIN+1;
    if (Temp_MAX>TC_Water)  Temp_MAX=TC_Water;

    T_SMINXREG56(Temp_MIN,Temp_MAX,X,T);
    return *T;
}


double SMINXREG56(double X)
{
    return TX2SREG56(SMINX2TREG56(X),X);
}


void  T_SMAXXREG56(double Temp_MIN,double Temp_MAX,double X, double* T)
{
double    Temp_Mid  ;
double    Temp_Mid_Plus , S_Mid_Plus ;
double    Temp_Mid_MINus , S_Mid_MINus ;
double    Temp_MIN2 , Temp_MAX2;

    if (fabs(Temp_MAX-Temp_MIN) < 0.02 )
    {
        *T = (Temp_MAX + Temp_MIN) / 2;
        return;
    }

    Temp_Mid = (Temp_MAX + Temp_MIN) / 2;

    Temp_Mid_Plus = Temp_Mid + 0.01;
    S_Mid_Plus = TX2SREG56(Temp_Mid_Plus, X);

    Temp_Mid_MINus = Temp_Mid - 0.01;
    S_Mid_MINus = TX2SREG56(Temp_Mid_MINus, X);


    Temp_MIN2=Temp_MIN;
    Temp_MAX2=Temp_MAX;

    if (S_Mid_Plus > S_Mid_MINus) 
        Temp_MIN2 = Temp_Mid;
    else
        Temp_MAX2 = Temp_Mid;

    T_SMAXXREG56(Temp_MIN2, Temp_MAX2, X, T);
}


double SMAXX2TREG56(double X)
{
double    Temp_MIN , Temp_MAX ;
double    *T,SMAX,TMAX,STX,t ;
int     IFor;
     
    T=&t;
    if (X >= 0.55 )
    {
        return TX2SREG56(T000C,X);
        //exit(1) ;
    }
    SMAX=-100;
    TMAX=TC_Water;    
	for(IFor=0;IFor<=375;IFor++)
    {
        *T=T000C+(TC_Water-T000C)/375*IFor;
        if (*T>TC_Water)
			break;
        STX=TX2SREG56(*T,X);
        if (STX>SMAX )
        {
            SMAX=STX;
            TMAX=*T;
        }
    }

    Temp_MIN = TMAX-1;
    if (Temp_MIN<T000C)  Temp_MIN=T000C;
    Temp_MAX = TMAX+1;
    if (Temp_MAX>TC_Water)  Temp_MAX=TC_Water;

    T_SMAXXREG56(Temp_MIN,Temp_MAX,X,T);
    return *T;
}


double SMAXXREG56(double X)
{
    return TX2SREG56(SMAXX2TREG56(X),X);
}


void  GETTRANGESX(double * Temp_MIN , double* Temp_MAX , double Entropy ,double X )
{
double    FMIN ,  FMid , FMidMIN;
double    Temp_Midmin, Temp_Mid ;

    if (fabs(*Temp_MIN-*Temp_MAX)<0.00001)  
		return;
		//exit(1);

    FMIN =TX2SREG56(*Temp_MIN, X) - Entropy;

    Temp_Mid = (*Temp_MIN + *Temp_MAX) / 2;
    FMid = TX2SREG56(Temp_Mid, X) - Entropy;

    if (FMIN * FMid < 0 )
        *Temp_MAX = Temp_Mid;
    else
    {
        Temp_Midmin = (*Temp_MIN + Temp_Mid) / 2;
        FMidMIN = TX2SREG56(Temp_Midmin, X) - Entropy;

        if (FMid * FMidMIN < 0 )
        {
            *Temp_MIN = Temp_Midmin;
            *Temp_MAX = Temp_Mid;
		}
        else
        {
            *Temp_MIN = Temp_Mid;
            GETTRANGESX(Temp_MIN, Temp_MAX, Entropy, X);
        }
    }
}


void SX2TREG56_HML(double Entropy,double X ,double * TH_KK,double *TM, double*TL,int* NBRS )
{
double    *Temp_MIN , *Temp_MAX , TOL,temp_max,temp_min ;
double    Temp ;
double    FMIN , FMAX ;

    Temp_MIN=&temp_min;
	Temp_MAX=&temp_max;
    *Temp_MIN = T000C;
    if (fabs(1 - Entropy / SC_Water) <= 1E-8 )
        *Temp_MAX = TC_Water - 0.005;
    else
        *Temp_MAX = TC_Water;
    FMIN =TX2SREG56(*Temp_MIN, X) - Entropy;
    FMAX =TX2SREG56(*Temp_MAX, X) - Entropy;

    if (FMIN * FMAX > 0 )
        GETTRANGESX(Temp_MIN, Temp_MAX, Entropy, X);

    TOL = 1E-8 * (*Temp_MIN + *Temp_MAX) / 2.0;
    *TM = ZBRENT2_KK(*Temp_MIN, *Temp_MAX, TOL, X, Entropy, 5604);


    *Temp_MAX = (*TM + T000C) / 2;
    *Temp_MIN = T000C;
    
    FMIN =TX2SREG56(*Temp_MIN, X) - Entropy;
    FMAX =TX2SREG56(*Temp_MAX, X) - Entropy;
    if (FMIN * FMAX < 0) 
    {
        TOL = (0.00000001) * (*Temp_MIN + *Temp_MAX) / 2.0;
        *TL = ZBRENT2_KK(*Temp_MIN, *Temp_MAX, TOL, X, Entropy, 5604);
    }
    else
        *TL = -1;

    if (fabs(1 - Entropy / SC_Water) <= 1E-8 )
        *TH_KK = TC_Water;
    else
    {
        *Temp_MAX = TC_Water;
        *Temp_MIN = (*TM + TC_Water) / 2;
    
        FMIN =TX2SREG56(*Temp_MIN, X) - Entropy;
        FMAX =TX2SREG56(*Temp_MAX, X) - Entropy;
        if (FMIN * FMAX < 0 )
        {
            TOL = (0.00000001) * (*Temp_MIN + *Temp_MAX) / 2.0;
            *TH_KK = ZBRENT2_KK(*Temp_MIN, *Temp_MAX, TOL, X, Entropy, 5604);
        }
        else
            *TH_KK = -1;
    }
    
    if ((*TM > 0) && (*TL > 0) && (*TH_KK < 0) )
    {
        if (*TL > *TM )
        {
            Temp = *TL;
            *TL = *TM;
            *TM = Temp;
        }


        *Temp_MAX = *TM - 0.005;
        *Temp_MIN = *TL + 0.005;

        FMIN =TX2SREG56(*Temp_MIN, X) - Entropy;
        FMAX =TX2SREG56(*Temp_MAX, X) - Entropy;
        if (FMIN * FMAX < 0) 
        {
            TOL = 1E-8 * (*Temp_MIN + *Temp_MAX) / 2.0;
            *TH_KK = ZBRENT2_KK(*Temp_MIN, *Temp_MAX, TOL, X, Entropy, 5604);
        }
    }
    else if ((*TM > 0) && (*TH_KK > 0) && (*TL < 0)) 
    {
        if (*TH_KK < *TM )
        {
            Temp = *TH_KK;
            *TH_KK = *TM;
            *TM = Temp;
        }


        *Temp_MAX = *TH_KK - 0.005;
        *Temp_MIN = *TM + 0.005;

        FMIN =TX2SREG56(*Temp_MIN, X) - Entropy;
        FMAX =TX2SREG56(*Temp_MAX, X) - Entropy;
        if (FMIN * FMAX < 0 )
        {
            TOL = 1E-8 * (*Temp_MIN + *Temp_MAX) / 2.0;
            *TL = ZBRENT2_KK(*Temp_MIN, *Temp_MAX, TOL, X, Entropy, 5604);
        }
    }


    if (*TL > *TM )
    {
      Temp = *TL;
      *TL = *TM;
      *TM = Temp;
    }
    
    if (*TM > *TH_KK) 
    {
      Temp = *TM;
      *TM = *TH_KK;
      *TH_KK = Temp;
    }
   
    if (*TL > *TM )
    {
      Temp = *TL;
      *TL = *TM;
      *TM = Temp;
    }
   
    *NBRS = 3;
    if (*TL < 0)  *NBRS = 2;
    if (*TM < 0)  *NBRS = 1;
}

//////////SX  REG56
double SX2TREG56(double  Entropy , double X )  
{
double    *TH_KK , *TM , *TL,th,tm, tl ;
int    *NBRS,nbrs ;

    TH_KK=&th;
	TM=&tm;
	TL=&tl;
	NBRS=&nbrs;
    SX2TREG56_HML(Entropy, X, TH_KK, TM, TL, NBRS);
    if (*TL > 0 )
    {
        return *TL; 
		//exit(1);
    }
    if (*TM > 0) 
    {
        return *TM;  
		//exit(1);
	}
    return *TH_KK;
}

double SX2TLPREG56(double  Entropy , double X )  
{
double    *TH_KK , *TM , *TL,th,tm,tl ;
int    *NBRS , nbrs;
    
    TH_KK=&th;
	TM=&tm;
	TL=&tl;
	NBRS=&nbrs; 
    SX2TREG56_HML(Entropy, X, TH_KK, TM, TL, NBRS);
    if (*TL > 0 )
    {
        return *TL;
        //exit(1);
    }
    if (*TM > 0) 
    {
        return *TM;  
		//exit(1);
    }
    return *TH_KK;
}

double SX2TMPREG56(double  Entropy , double X )  
{
double    *TH_KK , *TM , *TL, th,tm,tl ;
int    *NBRS ,nbrs;

    TH_KK=&th;
	TM=&tm;
	TL=&tl;
	NBRS=&nbrs;
    SX2TREG56_HML(Entropy, X, TH_KK, TM, TL, NBRS);
    if (*TM > 0 )
    {
        return *TM;        
    }
    return *TH_KK;
}


double SX2THPREG56(double  Entropy , double X )  
{
double    *TH_KK , *TM , *TL,th,tm,tl ;
int    *NBRS,nbrs ;

    TH_KK=&th;
	TM=&tm;
	TL=&tl;
	NBRS=&nbrs;
    SX2TREG56_HML(Entropy, X, TH_KK, TM, TL, NBRS);
    return *TH_KK;
}


double SX2TREG56NBRS(double  Entropy ,double X ) 
{
double    *TH_KK , *TM , *TL,th,tm,tl ;
int    *NBRS,nbrs ;

    TH_KK=&th;
	TM=&tm;
	TL=&tl;
	NBRS=&nbrs; 
    SX2TREG56_HML(Entropy, X, TH_KK, TM, TL, NBRS);
    return *NBRS;
}


double SX2PREG56(double  Entropy , double X ) 
{
double    Temperature;

    Temperature = SX2TREG56(Entropy, X);
    return T2P_KK(Temperature);
}


double SX2PLPREG56(double  Entropy ,double X )  
{
double    Temperature;

    Temperature = SX2TREG56(Entropy, X);
    return T2P_KK(Temperature);
}


double SX2PMPREG56(double  Entropy ,double X )  
{
double    Temperature;

    Temperature = SX2TMPREG56(Entropy, X);
    return T2P_KK(Temperature);
} //


double SX2PHPREG56(double  Entropy , double X )  
{
double    Temperature;

    Temperature = SX2THPREG56(Entropy, X);
    return T2P_KK(Temperature);
}


double SX2HREG56(double  Entropy ,double X )  
{
double    Temperature;

    Temperature = SX2TREG56(Entropy, X);
    return TX2HREG56(Temperature, X);
}


double SX2HLPREG56(double  Entropy , double X )  
{
double    Temperature;

    Temperature = SX2TREG56(Entropy, X);
    return TX2HREG56(Temperature, X);
}


double SX2HMPREG56(double  Entropy ,double X )  
{
double    Temperature;

    Temperature = SX2THPREG56(Entropy, X);
    return TX2HREG56(Temperature, X);
}


double SX2HHPREG56(double  Entropy , double X ) 
{
double    Temperature;

    Temperature = SX2THPREG56(Entropy, X);
    return TX2HREG56(Temperature, X);
}



double SX2VREG56(double  Entropy , double X )  
{
double    Temperature;

    Temperature = SX2TREG56(Entropy, X);
    return TX2VREG56(Temperature, X);
}


double SX2VLPREG56(double  Entropy , double X )  
{
double    Temperature;

    Temperature = SX2TREG56(Entropy, X);
    return TX2VREG56(Temperature, X);
}


double SX2VMPREG56(double  Entropy , double X)  
{
double    Temperature;

    Temperature = SX2TMPREG56(Entropy, X);
    return TX2VREG56(Temperature, X);
}


double SX2VHPREG56(double  Entropy ,double X)  
{
double    Temperature;

    Temperature = SX2THPREG56(Entropy, X);
    return TX2VREG56(Temperature, X);
}

void SXREG56(double* P,double*T,double*H,double S,double* V,double X)
{
    *T=SX2TREG56(S,X);
    *P=T2P_KK(*T);
    *H=TX2HREG56(*T,X);
    *V=TX2VREG56(*T,X);
}

void SXLPREG56(double* P,double*T,double*H,double S,double* V,double X)
{
    *T=SX2TLPREG56(S,X);
    *P=T2P_KK(*T);
    *H=TX2HREG56(*T,X);
    *V=TX2VREG56(*T,X);
}

void SXMPREG56(double* P,double*T,double*H, double S,  double* V, double X)
{
    *T=SX2TMPREG56(S,X);
    *P=T2P_KK(*T);
    *H=TX2HREG56(*T,X);
    *V=TX2VREG56(*T,X);
}

void  SXHPREG56(double* P,double*T,double*H, double S,  double* V, double X)
{
    *T=SX2THPREG56(S,X);
    *P=T2P_KK(*T);
    *H=TX2HREG56(*T,X);
    *V=TX2VREG56(*T,X);
}


void T_VMINXREG56(double Temp_MIN,double Temp_MAX,double X, double* T)
{
double    Temp_Mid  ;
double    Temp_Mid_Plus , V_Mid_Plus ;
double    Temp_Mid_MINus , V_Mid_MINus ;
double    Temp_MIN2 , Temp_MAX2;

    if (fabs(Temp_MAX - Temp_MIN) < 0.02 )
    {
        *T = (Temp_MAX + Temp_MIN) / 2;
        return;
		//exit(1);
    }

    Temp_Mid = (Temp_MAX + Temp_MIN) / 2;

    Temp_Mid_Plus = Temp_Mid + 0.01;
    V_Mid_Plus = TX2VREG56(Temp_Mid_Plus, X);

    Temp_Mid_MINus = Temp_Mid - 0.01;
    V_Mid_MINus = TX2VREG56(Temp_Mid_MINus, X);


    Temp_MIN2=Temp_MIN;
    Temp_MAX2=Temp_MAX;

    if (V_Mid_Plus < V_Mid_MINus) 
        Temp_MIN2 = Temp_Mid;
    else
        Temp_MAX2 = Temp_Mid;

    T_VMINXREG56(Temp_MIN2, Temp_MAX2, X, T);
}


double VMINX2TREG56(double X)
{
double    Temp_MIN , Temp_MAX ;
double    *T ,t;

    T=&t;
    if (X >= 0.442 )
    {
        return TC_Water;        
    }
    Temp_MIN = T000C;
    Temp_MAX = TC_Water;
    T_VMINXREG56(Temp_MIN,Temp_MAX,X,T);
    return *T;
}


double VMINXREG56(double X)
{
    return TX2VREG56(VMINX2TREG56(X),X);
}



void  GETTRANGEVX(double* Temp_MIN,double*Temp_MAX, double Volume,double  X)
{
double    FMIN , FMid , FMidMIN ;
double    Temp_MidMin , Temp_Mid ;

    if (fabs(*Temp_MIN-*Temp_MAX)<0.00001)  
		return;
		//exit(1);

    FMIN =TX2VREG56(*Temp_MIN, X)-Volume;

    Temp_Mid = (*Temp_MIN + *Temp_MAX) / 2;
    FMid =TX2VREG56(Temp_Mid, X) - Volume;

    if (FMIN * FMid < 0 )
        *Temp_MAX = Temp_Mid;
    else
    {
        Temp_MidMin = (*Temp_MIN + Temp_Mid) / 2;
        FMidMIN =TX2VREG56(Temp_MidMin, X) - Volume;

        if (FMid * FMidMIN < 0) 
        {
            *Temp_MIN = Temp_MidMin;
            *Temp_MAX = Temp_Mid;
		}
        else
        {
            *Temp_MIN = Temp_Mid;
            GETTRANGEVX(Temp_MIN, Temp_MAX, Volume, X);
		}
    }
}

/////////////VX    REG56
double VX2TREG56(double  Volume , double X )  
{
double    *Temp_MIN , *Temp_MAX , TOL , temp_min,temp_max;
double    FMIN , FMAX ;

    Temp_MIN=&temp_min;
	Temp_MAX=&temp_max;
    *Temp_MIN = T000C;
    *Temp_MAX = TC_Water;
    
    FMIN =TX2VREG56(*Temp_MIN, X) - Volume;
    FMAX =TX2VREG56(*Temp_MAX, X) - Volume;

    if (FMIN * FMAX > 0 )
        GETTRANGEVX(Temp_MIN, Temp_MAX, Volume, X);

    TOL = 1E-8 * (*Temp_MIN + *Temp_MAX) / 2.0;
    return ZBRENT2_KK(*Temp_MIN, *Temp_MAX, TOL, X, Volume, 5605);
}


double VX2TLPREG56(double  Volume ,double  X )  
{
double    *Temp_MIN , *Temp_MAX , TOL,temp_min,temp_max ;
double    FMIN , FMAX ;
    
    Temp_MIN=&temp_min;
	Temp_MAX=&temp_max;
    *Temp_MIN = T000C;
    *Temp_MAX = TC_Water;
    
    FMIN =TX2VREG56(*Temp_MIN, X) - Volume;
    FMAX =TX2VREG56(*Temp_MAX, X) - Volume;

    if (FMIN * FMAX > 0) 
        GETTRANGEVX(Temp_MIN, Temp_MAX, Volume, X);

    TOL = 1E-8 * (*Temp_MIN + *Temp_MAX) / 2.0;
    return ZBRENT2_KK(*Temp_MIN, *Temp_MAX, TOL, X, Volume, 5605);
}


double VX2THPREG56(double  Volume , double X )  
{
double    Temp_MIN , Temp_MAX , TOL ;
double    FMAX , FMIN ;
double    T_Low ;


    T_Low = VX2TREG56(Volume, X);


    Temp_MAX = TC_Water;
    Temp_MIN = (T_Low + Temp_MAX) / 2;
    
    FMIN =TX2VREG56(Temp_MIN, X) - Volume;
    FMAX =TX2VREG56(Temp_MAX, X) - Volume;
    if (FMIN * FMAX > 0 )
    {
        return T_Low; 
		//exit(1);
    }

    TOL = 1E-8 * (Temp_MIN + Temp_MAX) / 2.0;

    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, X, Volume, 5605);
}


double VX2PREG56(double  Volume , double X ) 
{
double   Temperature;

    Temperature = VX2TREG56(Volume, X);
    return T2P_KK(Temperature);
}


double VX2PLPREG56(double  Volume ,double X)  
{
double    Temperature;

    Temperature = VX2TLPREG56(Volume, X);
    return T2P_KK(Temperature);
}



double VX2PHPREG56(double  Volume ,double  X )  
{
double    Temperature;

    Temperature = VX2THPREG56(Volume, X);
    return T2P_KK(Temperature);
}


double VX2HREG56(double  Volume ,double X )  
{
double    Temperature;

    Temperature = VX2TREG56(Volume, X);
    return TX2HREG56(Temperature, X);
}


double VX2HLPREG56(double  Volume ,double X)  
{
double    Temperature;

    Temperature = VX2TLPREG56(Volume, X);
    return TX2HREG56(Temperature, X);
}


double VX2HHPREG56(double  Volume ,double X )  
{
double    Temperature;

    Temperature = VX2THPREG56(Volume, X);
    return TX2HREG56(Temperature, X);
}


double VX2SREG56(double  Volume , double X )  
{
double    Temperature;

    Temperature = VX2TREG56(Volume, X);
    return TX2SREG56(Temperature, X);
}


double VX2SLPREG56(double  Volume ,double X )  
{
double    Temperature;

    Temperature = VX2TLPREG56(Volume, X);
    return TX2SREG56(Temperature, X);
}


double VX2SHPREG56(double  Volume ,double  X )  
{
double    Temperature;

    Temperature = VX2THPREG56(Volume, X);
    return TX2SREG56(Temperature, X);
}

void VXREG56(double* P,double*T,double*H,double*S,double V,double X)
{
    *T=VX2TREG56(V,X);
    *P=T2P_KK(*T);
    *H=TX2HREG56(*T,X);
    *S=TX2SREG56(*T,X);
}

void VXLPREG56(double * P,double*T,double*H,double*S, double V,double X)
{
    *T=VX2TLPREG56(V,X);
    *P=T2P_KK(*T);
    *H=TX2HREG56(*T,X);
    *S=TX2SREG56(*T,X);
}

void VXHPREG56(double* P,double *T,double *H,double*S, double V,double X)
{
    *T=VX2THPREG56(V,X);
    *P=T2P_KK(*T);
    *H=TX2HREG56(*T,X);
    *S=TX2SREG56(*T,X);
}

///////PH2XREG1
double PH2TREG1(double P,double H)
{
double    Temp_MIN , Temp_MAX, TOL ;


    if ((P < 6.108E-4-DeltaVal) && (P > 100+DeltaVal) )
    {
        return -1.0;
		//exit(1);
    }

    Temp_MIN = T000C;

    if (P<16.535+DeltaVal) 
        Temp_MAX=P2T_KK(P);
    else
        Temp_MAX = T350C;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, H, 100);
}


double PH2SREG1(double P,double H)
{
    return PT2SREG1(P,PH2TREG1(P,H));
}


double PH2VREG1(double P,double H)
{
    return PT2VREG1(P,PH2TREG1(P,H));
}

void PHREG1(double P, double* T, double H, double* S,double*V,double*X)
{
    *T=PH2TREG1(P,H);
    *S=PT2SREG1(P,*T);
    *V=PT2VREG1(P,*T);
    *X=0.0;
}

/////PS2  REG1
double PS2TREG1(double P,double S)
{
double    Temp_MIN , Temp_MAX, TOL ;

    if ((P < 6.108E-4-DeltaVal) && (P > 100+DeltaVal)) 
    {
        return -1.0; 
		//exit(1);
	}

    Temp_MIN = T000C;
    if (P<16.535+DeltaVal) 
        Temp_MAX=P2T_KK(P);
    else
        Temp_MAX = T350C;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, S, 101);
}


double PS2HREG1(double P,double S)
{
    return PT2HREG1(P,PS2TREG1(P,S));
}


double PS2VREG1(double P,double S)
{
    return PT2VREG1(P,PS2TREG1(P,S));
}

void PSREG1(double P, double* T,double*H, double S, double* V,double*X)
{
    *T=PS2TREG1(P,S);
    *H=PT2HREG1(P,*T);
    *V=PT2VREG1(P,*T);
    *X=0.0;
}

/////////////////PV2 REG1
double PV2TREG1(double P,double V)
{
double    Temp_MIN , Temp_MAX, TOL ;


    if ( (P < 6.108E-4-DeltaVal) && (P > 100+DeltaVal) )
    {
        return -1.0; 
		//exit(1);
	}

    Temp_MIN = T000C;

    if ( P<16.535+DeltaVal) 
        Temp_MAX=P2T_KK(P);
    else
        Temp_MAX = T350C;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, V, 102);
}


double PV2HREG1(double P,double V)
{
    return PT2HREG1(P,PV2TREG1(P,V));
}


double PV2SREG1(double P,double V)
{
    return PT2SREG1(P,PV2TREG1(P,V));
}

void PVREG1(double P, double* T,double*H,double*S, double V, double* X)
{
    *T=PV2TREG1(P,V);
    *H=PT2HREG1(P,*T);
    *S=PT2SREG1(P,*T);
    *X=0.0;
}

///////TH2  REG1
double TH2PHPREG1(double  T,double H)
{
double    Press_MIN , Press_MAX , TOL ;


    Press_MIN = PHMINTREG1(T);
    Press_MAX = 100;

    if ((PT2HREG1(Press_MIN,T)-H)*(PT2HREG1(Press_MAX, T)-H)> 0 )
    {
        Press_MIN = T2P_KK(T);
        Press_MAX = PHMINTREG1(T);

        if (fabs(Press_MIN-Press_MAX)<10*DeltaVal) 
        {
            return (Press_MIN+Press_MAX)/2;            
		}

        if ((PT2HREG1(Press_MIN, T)-H)*(PT2HREG1(Press_MAX, T)-H)> 0) 
        {
            return -1.0;		
        }

        TOL = 1E-8*(Press_MIN+Press_MAX)/2.0;
        return ZBRENT2_KK(Press_MIN,Press_MAX,TOL,T,H, 103);        
    }

    TOL = 1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, T,H, 103);
}



double TH2PLPREG1(double  T,double H ) 
{
double    Press_MIN , Press_MAX , TOL;

    Press_MIN = T2P_KK(T);
    Press_MAX = PHMINTREG1(T);

    if ((PT2HREG1(Press_MIN, T)-H)*(PT2HREG1(Press_MAX, T)-H)> 0 )
    {
        return TH2PHPREG1(T,H);        
	}

    TOL = 1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN,Press_MAX,TOL,T,H, 103);
}



double TH2PREG1(double T,double H)
{
   return TH2PLPREG1(T,H);
}


double TH2SREG1(double T,double H)
{
    return PT2SREG1(TH2PREG1(T,H),T);
}


double TH2SLPREG1(double T,double H)
{
    return PT2SREG1(TH2PLPREG1(T,H),T);
}


double TH2SHPREG1(double T,double H)
{
    return PT2SREG1(TH2PHPREG1(T,H),T);
}



double TH2VREG1(double T,double H)
{
    return PT2VREG1(TH2PREG1(T,H),T);
}

double TH2VLPREG1(double T,double H)
{
    return PT2VREG1(TH2PLPREG1(T,H),T);
}

double TH2VHPREG1(double T,double H)
{
    return PT2VREG1(TH2PHPREG1(T,H),T);
}


void THREG1(double* P, double T,double H, double* S,double*V,double*X)
{
    *P=TH2PREG1(T,H);
    *S=PT2SREG1(*P,T);
    *V=PT2VREG1(*P,T);
    *X=0.0;
}

void THLPREG1(double* P,double T,double H,double*  S,double*V,double*X)
{
    *P=TH2PLPREG1(T,H);
    *S=PT2SREG1(*P,T);
    *V=PT2VREG1(*P,T);
    *X=0.0;
}

void THHPREG1(double* P, double T, double H, double* S,double*V,double*X)
{
    *P=TH2PHPREG1(T,H);
    *S=PT2SREG1(*P,T);
    *V=PT2VREG1(*P,T);
    *X=0.0;
}


double P_SMAXTREG1(double PMIN,double PMAX,double T)
{
double    P_MIN,P_MAX,P_Mid;

    if (fabs(PMAX-PMIN)<10*DeltaVal) 
    {
        return (PMIN+PMAX)/2;  
		//exit(1);
    }

    P_MIN=PMIN;
    P_MAX=PMAX;
    P_Mid=(P_MIN+P_MAX)/2;
    if (PT2SREG1(P_Mid-2*DeltaVal,T)<PT2SREG1(P_Mid+2*DeltaVal,T)) 
    {
        P_MIN=P_Mid-2*DeltaVal;
    }
    else
    {
        P_MAX=P_Mid-2*DeltaVal;
    }

    return P_SMAXTREG1(P_MIN,P_MAX,T);
}


double PSMAXTREG1(double T)
{
double    P_MIN,P_MAX;

    P_MIN=T2P_KK(T);
    P_MAX=100;
    return P_SMAXTREG1(P_MIN,P_MAX,T);
}


double SMAXTREG1(double T)
{
    return PT2SREG1(PSMAXTREG1(T),T);
}

///////TS2  REG1
double TS2PHPREG1(double T,double S)
{
double    Press_MIN , Press_MAX, TOL ;
double    PSMAX,SMIN2;

    if ((T < T000C-DeltaVal) && (T > T350C+DeltaVal)) 
    {
        return -1.0; 
		//exit(1);
    }
    if (T>5 )
    {
        Press_MIN=T2P_KK(T);
        Press_MAX=100;
	}
    else
    {
        PSMAX=PSMAXTREG1(T);
        SMIN2=T2SLREG56(T);
        if (S<SMIN2) 
        {
            Press_MIN=PSMAX;
            Press_MAX=100;
		}
        else
        {
            Press_MIN=PSMAX;
            Press_MAX=100;
        }
    }
    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, T, S, 104);
}

double TS2PLPREG1(double T,double S)
{
double   Press_MIN , Press_MAX, TOL ;
double    PSMAX,SMIN2;

    if ((T < T000C-DeltaVal) && (T > T350C+DeltaVal) )
    {
        return -1.0;  
		//exit(1);
    }

    if (T>5) 
    {
        Press_MIN=T2P_KK(T);
        Press_MAX=100;
    }
    else
    {
        PSMAX=PSMAXTREG1(T);
        SMIN2=T2SLREG56(T);
        if (S<SMIN2 )
        {
            Press_MIN=PSMAX;
            Press_MAX=100;
        }
        else
        {
            Press_MIN=T2P_KK(T);
            Press_MAX=PSMAX;
        }
    }
    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, T, S, 104);
}

double TS2PREG1(double T,double S)
{
    return TS2PLPREG1(T,S);
}

double TS2HHPREG1(double T,double S)
{
    return PT2HREG1(TS2PHPREG1(T,S),T);
}

double TS2HLPREG1(double T,double S)
{
    return PT2HREG1(TS2PLPREG1(T,S),T);
}

double TS2HREG1(double T,double S)
{
    return PT2HREG1(TS2PREG1(T,S),T);
}

double TS2VHPREG1(double T,double S)
{
    return PT2VREG1(TS2PHPREG1(T,S),T);
}

double TS2VLPREG1(double T,double S)
{
    return PT2VREG1(TS2PLPREG1(T,S),T);
}

double TS2VREG1(double T,double S)
{
    return PT2VREG1(TS2PREG1(T,S),T);
}

void TSHPREG1(double* P,double T, double* H, double S, double* V,double*X)
{
    *P=TS2PHPREG1(T,S);
    *H=PT2HREG1(*P,T);
    *V=PT2VREG1(*P,T);
    *X=0.0;
}

void TSLPREG1(double* P,double T, double* H, double S, double* V,double*X)
{
    *P=TS2PLPREG1(T,S);
    *H=PT2HREG1(*P,T);
    *V=PT2VREG1(*P,T);
    *X=0.0;
}

void TSREG1(double* P,double T, double* H, double S, double* V,double*X)
{
    *P=TS2PREG1(T,S);
    *H=PT2HREG1(*P,T);
    *V=PT2VREG1(*P,T);
    *X=0.0;
}

///////////TV2 REG1
double TV2PREG1(double T,double V)
{
double    Press_MIN , Press_MAX, TOL ;

    if ((T < T000C-DeltaVal) && (T > T350C+DeltaVal)) 
    {
        return -1.0;  
		//exit(1);
    }

    Press_MIN=T2P_KK(T);
    Press_MAX=100;

    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, T, V, 105);
}

double TV2HREG1(double T,double V)
{
    return PT2HREG1(TV2PREG1(T,V),T);
}

double TV2SREG1(double T,double V)
{
    return PT2SREG1(TV2PREG1(T,V),T);
}

void TVREG1(double*P, double T,double* H,double*S, double V, double* X)
{
    *P=TV2PREG1(T,V);
    *H=PT2HREG1(*P,T);
    *S=PT2SREG1(*P,T);
    *X=0.0;
}

//////HS2  REG1
double HS2PREG1(double H,double S)
{
double    Press_MIN , Press_MAX, TOL ;


    if (S<=-0.0002-DeltaVal) 
        Press_MIN=TS2PREG1(T000C,S);
    else
        Press_MIN=S2PLREG56(S);

    if (S<=3.3922+DeltaVal )
        Press_MAX=100;
    else
        Press_MAX=TS2PREG1(T350C,S);

    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, S, H, 106);
}

double HS2TREG1(double H,double S)
{
    return PH2TREG1(HS2PREG1(H,S),H);
}

double HS2VREG1(double H,double S)
{
    return PH2VREG1(HS2PREG1(H,S),H);
}


void HSREG1(double* P,double*T,double H,double S, double* V,double*X)
{
    *P=HS2PREG1(H,S);
    *T=PH2TREG1(*P,H);
    *V=PH2VREG1(*P,H);
    *X=0.0;
}

///////HV2 REG1
double HV2PREG1(double H,double V)
{
double    Press_MIN , Press_MAX, TOL ;

    if (V<=0.0010001*(1-DeltaVal) )
        Press_MIN=TV2PREG1(T000C,V);
    else
        Press_MIN=V2PLREG56(V);

    if (V<=0.0013148*(1+DeltaVal) )
        Press_MAX=100;
    else
        Press_MAX=TV2PREG1(T350C,V);


    if (FUNC_PS2HREG1(Press_MIN,V,H)*FUNC_PS2HREG1(Press_MIN,V,H)<0 )
    {
        TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
        return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, V, H, 107);
    }
    else
    {
        Press_MIN=H2PLREG56(H);
        if (H>=95.95*(1+DeltaVal) )
            Press_MAX=100;
        else
            Press_MAX=TH2PREG1(T000C,H);

        TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
        return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, H, V, 109);
    }
}

double HV2TREG1(double H,double V)
{
    return PH2TREG1(HV2PREG1(H,V),H);
}

double HV2SREG1(double H,double V)
{
    return PH2SREG1(HV2PREG1(H,V),H);
}

void HVREG1(double* P,double*T, double H, double* S, double V, double* X)
{
    *P=HV2PREG1(H,V);
    *T=PH2TREG1(*P,H);
    *S=PH2SREG1(*P,H);
    *X=0.0;
}


///////SV2  REG1
double SV2PREG1(double S,double V)
{
double    Press_MIN , Press_MAX, TOL;


    if (V<=0.0010001*(1-DeltaVal)) 
        Press_MIN=TV2PREG1(T000C,V);
    else
        Press_MIN=V2PLREG56(V);

    if (V<=0.0013148*(1+DeltaVal) )
        Press_MAX=100;
    else
        Press_MAX=TV2PREG1(T350C,V);

    if (FUNC_PV2SREG1(Press_MIN,V,S)*FUNC_PV2SREG1(Press_MAX,V,S)<0 )
    {
        TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
        return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, V, S, 108);
    }
    else
    {
        if (S<=-0.0002*(1+DeltaVal)) 
            Press_MIN=TS2PREG1(T000C,S);
        else
            Press_MIN=S2PLREG56(S);

        Press_MAX=100;

        TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
        return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, S, V, 110);
    }
}

double SV2TREG1(double S,double V)
{
    return PS2TREG1(SV2PREG1(S,V),S);
}

double SV2HREG1(double S,double V)
{
    return PS2HREG1(SV2PREG1(S,V),S);
}

void SVREG1(double *P,double *T,double *H, double S,double V,double* X)
{
    *P=SV2PREG1(S,V);
    *T=PS2TREG1(*P,S);
    *H=PS2HREG1(*P,S);
    *X=0.0;
}

//////PH2  REG2
double PH2TREG2(double P,double H)
{
double    Temp_MIN , Temp_MAX, TOL ;

    if (P<=16.535+DeltaVal) 
        Temp_MIN = P2T_KK(P);
    else
        Temp_MIN =P2TBOUND23(P);
    Temp_MAX = 800;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, H, 200);
}


double PH2SREG2(double P,double H)
{
	return PT2SREG2(P,PH2TREG2(P,H));
}


double PH2VREG2(double P,double H)
{
    return PT2VREG2(P,PH2TREG2(P,H));
}

void PHREG2(double P, double* T, double  H, double* S,double*V,double*X)
{
    *T=PH2TREG2(P,H);
    *S=PT2SREG2(P,*T);
    *V=PT2VREG2(P,*T);
    *X=1.0;
}

//////PS2 REG2
double PS2TREG2(double P,double S)
{
double    Temp_MIN , Temp_MAX, TOL ;


    if (P<=16.535+DeltaVal )
        Temp_MIN = P2T_KK(P);
    else
        Temp_MIN =P2TBOUND23(P);
    Temp_MAX = 800;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, S, 201);
}


double PS2HREG2(double P,double S)
{
    return PT2HREG2(P,PS2TREG2(P,S));
}

double PS2VREG2(double P,double S)
{
    return PT2VREG2(P,PS2TREG2(P,S));
}


void PSREG2(double P, double* T,double*H, double S, double* V,double*X)
{
    *T=PS2TREG2(P,S);
    *H=PT2HREG2(P,*T);
    *V=PT2VREG2(P,*T);
    *X=1.0;
}

double PV2TREG2(double P,double V)
{
double    Temp_MIN , Temp_MAX, TOL;

    if (P<=16.535+DeltaVal) 
        Temp_MIN = P2T_KK(P);
    else
        Temp_MIN =P2TBOUND23(P);
    Temp_MAX = 800;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, V, 202);
}


double PV2HREG2(double P,double V)
{
    return PT2HREG2(P,PV2TREG2(P,V));
}

double PV2SREG2(double P,double V)
{
    return PT2SREG2(P,PV2TREG2(P,V));
}

void PVREG2(double P, double* T,double*H,double*S, double V,double* X)
{
    *T=PV2TREG2(P,V);
    *H=PT2HREG2(P,*T);
    *S=PT2SREG2(P,*T);
    *X=1.0;
}

double TH2PREG2(double T,double H)
{
double    Press_MIN , Press_MAX, TOL;

    Press_MIN = 6.108E-4;

    if (T<T350C+DeltaVal) 
        Press_MAX=T2P_KK(T);
    else if (T<T590C+DeltaVal) 
        Press_MAX=T2PBOUND23(T);
    else
        Press_MAX=100;

    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, T, H, 203);
}

double TH2SREG2(double T,double H)
{
    return PT2SREG2(TH2PREG2(T,H),T);
}

double TH2VREG2(double T,double H)
{
    return PT2VREG2(TH2PREG2(T,H),T);
}

void THREG2(double*P, double T,double H, double* S,double*V,double*X)
{
    *P=TH2PREG2(T,H);
    *S=PT2SREG2(*P,T);
    *V=PT2VREG2(*P,T);
    *X=1.0;
}

double TS2PREG2(double T,double S)
{
double    Press_MIN , Press_MAX, TOL ;

    Press_MIN = 6.108E-4;

    if (T<T350C+DeltaVal) 
        Press_MAX=T2P_KK(T);
    else if (T<T590C+DeltaVal) 
        Press_MAX=T2PBOUND23(T);
    else
        Press_MAX=100;

    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, T, S, 204);
}

double TS2HREG2(double T,double S)
{
    return PT2HREG2(TS2PREG2(T,S),T);
}

double TS2VREG2(double T,double S)
{
    return PT2VREG2(TS2PREG2(T,S),T);
}

void TSREG2(double* P, double T, double* H, double S, double* V,double*X)
{
    *P=TS2PREG2(T,S);
    *H=PT2HREG2(*P,T);
    *V=PT2VREG2(*P,T);
    *X=1.0;
}

double TV2PREG2(double T,double V)
{
double    Press_MIN , Press_MAX, TOL ;


    Press_MIN = 6.108E-4;

    if (T<T350C+DeltaVal) 
        Press_MAX=T2P_KK(T);
    else if ( T<T590C+DeltaVal )
        Press_MAX=T2PBOUND23(T);
    else
        Press_MAX=100;

    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, T, V, 205);
}

double TV2HREG2(double T,double V)
{
    return PT2HREG2(TV2PREG2(T,V),T);
}

double TV2SREG2(double T,double V)
{
    return PT2SREG2(TV2PREG2(T,V),T);
}


void TVREG2(double* P, double T, double* H,double*S, double V,double* X)
{
    *P=TV2PREG2(T,V);
    *H=PT2HREG2(*P,T);
    *S=PT2SREG2(*P,T);
    *X=1.0;
}

double HS2PREG2(double H,double S)
{
double    Press_MIN , Press_MAX, TOL,T ;


    T=HMAXX2TREG56(1.0);
    if (H>T2HGREG56(T)-DeltaVal )
	{
        Press_MIN=6.108E-4;
        if (H<PT2HREG2(100,590) )
            Press_MAX=H2PBOUND23(H);
        else if (H<PT2HREG2(100,800)) 
            Press_MAX=100;
        else
            Press_MAX=TH2PREG2(800,H);
    }
    else
    {
        if (S>T2SGREG56(T)) 
        {
            Press_MIN=6.108E-4;
            Press_MAX=H2PGLREG56(H);
        }
        else
        {
            Press_MIN=H2PGHREG56(H);
            Press_MAX=H2PBOUND23(H);
        }
    }

    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, H, S, 206);
}

double HS2TREG2(double H,double S)
{
    return PH2TREG2(HS2PREG2(H,S),H);
}

double HS2VREG2(double H,double S)
{
    return PH2VREG2(HS2PREG2(H,S),H);
}

void HSREG2(double* P,double*T, double H,double S, double* V,double*X)
{
    *P=HS2PREG2(H,S);
    *T=PH2TREG2(*P,H);
    *V=PH2VREG2(*P,H);
    *X=1.0;
}

double HV2PREG2(double H,double V)
{
double    Press_MIN , Press_MAX, TOL ,T;

    T=HMAXX2TREG56(1.0);
    if (H>T2HGREG56(T)-DeltaVal )
    {
        Press_MIN=6.108E-4;
        if (H<PT2HREG2(100,590)) 
            Press_MAX=H2PBOUND23(H);
        else if (H<PT2HREG2(100,800) )
            Press_MAX=100;
        else
            Press_MAX=TH2PREG2(800,H);
    }
    else
    {
        if ( V>T2VGREG56(T) )
        {
            Press_MIN=6.108E-4;
            Press_MAX=H2PGLREG56(H);
        }
        else
        {
            Press_MIN=H2PGHREG56(H);

            if ((H>=2610.68079) && (H<=2629.04908)) 
            {
                Press_MAX=T2PBOUND23(490.32);
            }
            else
                Press_MAX=H2PBOUND23(H);


        }
    }

    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, H, V, 207);
}

double HV2TREG2(double H,double V)
{
    return PH2TREG2(HV2PREG2(H,V),H);
}

double HV2SREG2(double H,double V)
{
    return PH2SREG2(HV2PREG2(H,V),H);
}

void HVREG2(double *P,double* T, double  H, double* S, double V, double* X)
{
    *P=HV2PREG2(H,V);
    *T=PH2TREG2(*P,H);
    *S=PH2SREG2(*P,H);
    *X=1.0;
}

double SV2PREG2(double S,double V)
{
double    Press_MIN , Press_MAX, TOL ;


    if (V>T2VGREG56(T000C)+DeltaVal )
    {
        Press_MIN=6.108E-4;
    }
    else if (V>T2VGREG56(T350C)-DeltaVal) 
    {
        Press_MIN=V2PGREG56(V);
    }
    else
    {
        Press_MIN=V2PBOUND23(V);
    }

    if (V>PT2VREG2(100,800) )
        Press_MAX=TV2PREG2(800,V);
    else
        Press_MAX=100;

    TOL =1E-8*(Press_MIN+Press_MAX)/2.0;
    return ZBRENT2_KK(Press_MIN, Press_MAX, TOL, V, S, 208);
}

double SV2TREG2(double S,double V)
{
    return PS2TREG2(SV2PREG2(S,V),S);
}

double SV2HREG2(double S,double V)
{
    return PS2HREG2(SV2PREG2(S,V),S);
}

void SVREG2(double* P,double*T,double*H, double S,double  V, double* X)
{
    *P=SV2PREG2(S,V);
    *T=PS2TREG2(*P,S);
    *H=PS2HREG2(*P,S);
    *X=1.0;
}

double P100T2VREG34(double T)
{
double    Volume_MIN, Volume_MAX,TOL;

    Volume_MIN=PT2VREG1(100,350);
    Volume_MAX=PT2VREG2(100,590);

    TOL =1E-8*(Volume_MIN+Volume_MAX)/2.0;
    return ZBRENT1_KK(Volume_MIN, Volume_MAX, TOL, T,3400);
}

double P100V2TREG4(double V)
{
double    Temp_MIN, Temp_MAX,TOL;

    Temp_MIN=T350C;
    Temp_MAX=TC_Water;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, V, 3401);
}



double P100V2TREG3(double V)
{
double    Temp_MIN, Temp_MAX,TOL;

    Temp_MIN=TC_Water;
    Temp_MAX=590;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, V, 3402);
}


double P100S2TREG4(double S)
{
double    Temp_MIN, Temp_MAX,TOL;

    Temp_MIN=T350C;
    Temp_MAX=TC_Water;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, S, 3404);
}

double P100S2TREG3(double S)
{
double    Temp_MIN, Temp_MAX,TOL;

    Temp_MIN=TC_Water;
    Temp_MAX=590;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, S, 3403);
}

double P100H2TREG3(double H)
{
double    Temp_MIN, Temp_MAX,TOL;

    Temp_MIN=TC_Water;
    Temp_MAX=590;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT1_KK(Temp_MIN, Temp_MAX, TOL, H, 3405);
}

double PT2VREG4(double P,double T)
{
double    Volume_MIN, Volume_MAX,TOL;


    Volume_MIN=P100T2VREG34(T);
    Volume_MAX=T2VLREG56(T);

    TOL =1E-8*(Volume_MIN+Volume_MAX)/2.0;
    return ZBRENT2_KK(Volume_MIN, Volume_MAX, TOL, T, P, 400);
}

double PT2HREG4(double P,double T)
{
    return TV2HREG4(T,PT2VREG4NT(P,T));
}

double PT2SREG4(double P,double T)
{
    return  TV2SREG4(T,PT2VREG4NT(P,T));
}



void  PTREG4 (double P,double T,  double* H,double*S,double*V,double*X)
{
    *V=PT2VREG4NT(P,T);
    *H=TV2HREG4(T,*V);
    *S=TV2SREG4(T,*V);
    *X=0.0;
}

double TH2VREG4(double T,double H)
{
double    Volume_MIN, Volume_MAX,TOL;

    Volume_MIN=P100T2VREG34(T);
    Volume_MAX=T2VLREG56(T);

    TOL =1E-8*(Volume_MIN+Volume_MAX)/2.0;
    return ZBRENT2_KK(Volume_MIN, Volume_MAX, TOL, T, H, 401);
}


double TH2PREG4(double T,double H)
{
    return TV2PREG4(T,TH2VREG4(T,H));
}

double TH2SREG4(double T,double H)
{
    return TV2SREG4(T,TH2VREG4(T,H));
}

void THREG4(double* P, double T,double H, double* S,double*V,double*X)
{
    *V=TH2VREG4(T,H);
    *P=TV2PREG4(T,*V);
    *S=TV2SREG4(T,*V);
    *X=0.0;
}

double TS2VREG4(double T,double S)
{
double    Volume_MIN, Volume_MAX,TOL;

    Volume_MIN=P100T2VREG34(T);
    Volume_MAX=T2VLREG56(T);

    TOL =1E-8*(Volume_MIN+Volume_MAX)/2.0;
    return ZBRENT2_KK(Volume_MIN, Volume_MAX, TOL, T, S, 402);
}

double TS2PREG4(double T,double S)
{
    return  TV2PREG4(T,TS2VREG4(T,S));
}


double TS2HREG4(double T,double S)
{
    return TV2HREG4(T,TS2VREG4(T,S));
}


void TSREG4(double* P, double T, double* H, double S, double*V,double*X)
{
    *V=TS2VREG4(T,S);
    *P=TV2PREG4(T,*V);
    *H=TV2HREG4(T,*V);
    *X=1.0;
}

double PV2TREG4(double P,double V)
{
double    Temp_MIN, Temp_MAX,TOL;

    if (V<T2VLREG56(T350C)) 
        Temp_MIN=T350C;
    else
        Temp_MIN=V2TLREG56(V);

    if (V<PT2VREG4NT(100,T350C)) 
        Temp_MAX=P100V2TREG4(V);
    else
        Temp_MAX=TC_Water;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, V, P, 403);
}

double PV2HREG4(double P,double V)
{
    return TV2HREG4(PV2TREG4(P,V),V);
}

double PV2SREG4(double P,double V)
{
    return TV2SREG4(PV2TREG4(P,V),V);
}

void PVREG4(double P, double* T,double*H,double*S, double V,double* X)
{
    *T=PV2TREG4(P,V);
    *H=PT2HREG4(P,*T);
    *S=PT2SREG4(P,*T);
    *X=0.0;
}

double HV2TREG4(double H,double V)
{
double    Temp_MIN, Temp_MAX,TOL;

    if (V<T2VLREG56(T350C) )
        Temp_MIN=T350C;
    else
        Temp_MIN=V2TLREG56(V);

    if (V<PT2VREG4NT(100,T350C) )
        Temp_MAX=P100V2TREG4(V);
    else
        Temp_MAX=TC_Water;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, V, H, 404);
}

double HV2PREG4(double H,double V)
{
    return TV2PREG4(HV2TREG4(H,V),V);
}

double HV2SREG4(double H,double V)
{
    return TV2SREG4(HV2TREG4(H,V),V);
}

void HVREG4(double* P,double*T, double H, double* S, double V, double* X)
{
    *T=HV2TREG4(H,V);
    *P=TV2PREG4(*T,V);
    *S=TV2SREG4(*T,V);
    *X=0.0;
}

double SV2TREG4(double S,double V)
{
double    Temp_MIN, Temp_MAX,TOL;

    if (V<T2VLREG56(T350C)) 
        Temp_MIN=T350C;
    else
        Temp_MIN=V2TLREG56(V);

    if (V<PT2VREG4NT(100,T350C)) 
        Temp_MAX=P100V2TREG4(V);
    else
        Temp_MAX=TC_Water;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, V, S, 405);
}

double SV2PREG4(double S,double V)
{
    return TV2PREG4(SV2TREG4(S,V),V);
}

double SV2HREG4(double S,double V)
{
    return TV2HREG4(SV2TREG4(S,V),V);
}

void SVREG4(double* P,double* T,double*H, double S,double V, double* X)
{
    *T=SV2TREG4(S,V);
    *P=TV2PREG4(*T,V);
    *H=TV2HREG4(*T,V);
    *X=0.0;
}

double PH2TREG4(double P,double H)
{
double    Temp_MIN, Temp_MAX,TOL;

    Temp_MIN=T350C;

    if (P<PC_Water+DeltaVal) 
        Temp_MAX=P2T_KK(P);
    else
        Temp_MAX=TC_Water;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, H, 406);
}

double PH2VREG4(double P,double H)
{
    return PT2VREG4NT(P,PH2TREG4(P,H));
}

double PH2SREG4(double P,double H)
{
    return PT2SREG4(P,PH2TREG4(P,H));
}

void PHREG4(double P, double* T, double H,  double* S,double*V,double*X)
{
    *T=PH2TREG4(P,H);
    *S=PT2SREG4(P,*T);
    *V=PT2VREG4NT(P,*T);
    *X=0.0;
}

double PS2TREG4(double P,double S)
{
double    Temp_MIN, Temp_MAX,TOL;

    Temp_MIN=T350C;

    if (P<PC_Water+DeltaVal) 
        Temp_MAX=P2T_KK(P);
    else
        Temp_MAX=TC_Water;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, S, 407);
}

double PS2VREG4(double P,double S)
{
    return PT2VREG4NT(P,PS2TREG4(P,S));
}

double PS2HREG4(double P,double S)
{
    return PT2HREG4(P,PS2TREG4(P,S));
}

void PSREG4(double P, double* T,double*H,  double S,  double* V,double*X)
{
    *T=PS2TREG4(P,S);
    *H=PT2HREG4(P,*T);
    *V=PT2VREG4NT(P,*T);
    *X=0.0;
}

double HS2TREG4(double H,double S)
{
double    Temp_MIN, Temp_MAX,TOL;

    if (S<T2SLREG56(T350C)+DeltaVal )
        Temp_MIN=T350C;
    else
        Temp_MIN=S2TLREG56(S);

    if (S<PT2SREG4(100,TC_Water)+DeltaVal )
        Temp_MAX=P100S2TREG4(S);
    else
        Temp_MAX=TC_Water;

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, S, H, 408);
}


double HS2PREG4(double H,double S)
{
    return TH2PREG4(HS2TREG4(H,S),H);
}

double HS2VREG4(double H,double S)
{
    return TH2VREG4(HS2TREG4(H,S),H);
}

void HSREG4(double* P,double* T, double H,double S, double* V,double*X)
{
    *T=HS2TREG4(H,S);
    *P=TH2PREG4(*T,H);
    *V=TH2VREG4(*T,H);
    *X=0.0;
}


double PT2VREG3(double P,double T)
{
double    Volume_MIN, Volume_MAX,TOL;

    if (T<TC_Water-DeltaVal )
        Volume_MIN=T2VGREG56(T);
    else
        Volume_MIN=P100T2VREG34(T);

    Volume_MAX=PT2VREG2(T2PBOUND23(T),T);

    TOL =1E-8*(Volume_MIN+Volume_MAX)/2.0;
    return ZBRENT2_KK(Volume_MIN, Volume_MAX, TOL, T, P, 300);
}


double PT2HREG3(double P,double T)
{
    return TV2HREG3(T,PT2VREG3NT(P,T));
}

double PT2SREG3(double P,double T)
{
    return TV2SREG3(T,PT2VREG3NT(P,T));
}

void  PTREG3(double P,double T, double* H,double* S,double*V,double*X)
{
    *V=PT2VREG3NT(P,T);
    *H=TV2HREG3(T,*V);
    *S=TV2SREG3(T,*V);
    *X=1.0;
}

double TH2VREG3(double T,double H)
{
double    Volume_MIN, Volume_MAX,TOL;

    if (T<TC_Water-DeltaVal) 
        Volume_MIN=T2VGREG56(T);
    else
        Volume_MIN=P100T2VREG34(T);

    Volume_MAX=PT2VREG2(T2PBOUND23(T),T);

    TOL =1E-8*(Volume_MIN+Volume_MAX)/2.0;
    return ZBRENT2_KK(Volume_MIN, Volume_MAX, TOL, T, H, 301);
}


double TH2PREG3(double T,double H)
{
    return TV2PREG3(T,TH2VREG3(T,H));
}

double TH2SREG3(double T,double H)
{
    return TV2SREG3(T,TH2VREG3(T,H));
}

void THREG3(double* P, double T,double H, double* S,double*V,double*X)
{
    *V=TH2VREG3(T,H);
    *P=TV2PREG3(T,*V);
    *S=TV2SREG3(T,*V);
    *X=1.0;
}

double TS2VREG3(double T,double S)
{
double    Volume_MIN, Volume_MAX,TOL;

    if (T<TC_Water-DeltaVal )
        Volume_MIN=T2VGREG56(T);
    else
        Volume_MIN=P100T2VREG34(T);

    Volume_MAX=PT2VREG2(T2PBOUND23(T),T);

    TOL =1E-8*(Volume_MIN+Volume_MAX)/2.0;
    return ZBRENT2_KK(Volume_MIN, Volume_MAX, TOL, T, S, 302);
}

double TS2PREG3(double T,double S)
{
    return TV2PREG3(T,TS2VREG3(T,S));
}


double TS2HREG3(double T,double S)
{
    return TV2HREG3(T,TS2VREG3(T,S));
}

void TSREG3(double* P, double  T,double* H, double S, double* V,double*X)
{
    *V=TS2VREG3(T,S);
    *P=TV2PREG3(T,*V);
    *H=TV2HREG3(T,*V);
    *X=1.0;
}

double PV2TREG3(double P,double V)
{
double    Temp_MIN, Temp_MAX,TOL;

    if (V<VC_Water )
        Temp_MIN=TC_Water;
    else
        Temp_MIN=V2TGREG56(V);

    if (V<PT2VREG2(100,590) )
        Temp_MAX=P100V2TREG3(V);
    else
        Temp_MAX=V2TBOUND23(V);

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, V, P, 303);
}

double PV2HREG3(double P,double V)
{
    return TV2HREG3(PV2TREG3(P,V),V);
}

double PV2SREG3(double P,double V)
{
    return TV2SREG3(PV2TREG3(P,V),V);
}

void PVREG3(double P,double* T,double*H,double*S,double V, double* X)
{
    *T=PV2TREG3(P,V);
    *H=PT2HREG3(P,*T);
    *S=PT2SREG3(P,*T);
    *X=1.0;
}

double HV2TREG3(double H,double V)
{
double    Temp_MIN, Temp_MAX,TOL;


    if (V<VC_Water) 
        Temp_MIN=TC_Water;
    else
        Temp_MIN=V2TGREG56(V);

    if (V<PT2VREG2(100,590)) 
        Temp_MAX=P100V2TREG3(V);
    else
        Temp_MAX=V2TBOUND23(V);

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, V, H, 304);
}

double HV2PREG3(double H,double V)
{
    return TV2PREG3(HV2TREG3(H,V),V);
}

double HV2SREG3(double H,double V)
{
    return TV2SREG3(HV2TREG3(H,V),V);
}

void HVREG3(double* P,double*T, double H,double* S, double  V, double* X)
{
    *T=HV2TREG3(H,V);
    *P=TV2PREG3(*T,V);
    *S=TV2SREG3(*T,V);
    *X=1.0;
}

double SV2TREG3(double S,double V)
{
double    Temp_MIN, Temp_MAX,TOL;


    if (V<VC_Water )
        Temp_MIN=TC_Water;
    else
        Temp_MIN=V2TGREG56(V);

    if (V<PT2VREG2(100,590) )
        Temp_MAX=P100V2TREG3(V);
    else
        Temp_MAX=V2TBOUND23(V);

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, V, S, 305);
}

double SV2PREG3(double S,double V)
{
    return TV2PREG3(SV2TREG3(S,V),V);
}

double SV2HREG3(double S,double V)
{
    return TV2HREG3(SV2TREG3(S,V),V);
}

void SVREG3(double* P,double*T,double*H, double S,double V, double* X)
{
    *T=SV2TREG3(S,V);
    *P=TV2PREG3(*T,V);
    *H=TV2HREG3(*T,V);
    *X=1.0;
}

double PH2TREG3(double P,double H)
{
double    Temp_MIN, Temp_MAX,TOL;


    if (P>PC_Water+DeltaVal) 
        Temp_MIN=TC_Water;
    else
        Temp_MIN=P2T_KK(P);

    Temp_MAX=P2TBOUND23(P);

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, H, 306);
}

double PH2VREG3(double P,double H)
{
    return PT2VREG3NT(P,PH2TREG3(P,H));
}

double PH2SREG3(double P,double H)
{
    return PT2SREG3(P,PH2TREG3(P,H));
}

void PHREG3(double P, double* T, double H, double* S,double*V,double*X)
{
    *T=PH2TREG3(P,H);
    *S=PT2SREG3(P,*T);
    *V=PT2VREG3NT(P,*T);
    *X=1.0;
}

double PS2TREG3(double P,double S)
{
double   Temp_MIN, Temp_MAX,TOL;

    if (P>PC_Water+DeltaVal) 
        Temp_MIN=TC_Water;
    else
        Temp_MIN=P2T_KK(P);

    Temp_MAX=P2TBOUND23(P);

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, P, S, 307);
}

double PS2VREG3(double P,double S)
{
    return PT2VREG3NT(P,PS2TREG3(P,S));
}

double PS2HREG3(double P,double S)
{
    return PT2HREG3(P,PS2TREG3(P,S));
}

void PSREG3(double P, double* T,double*H, double S, double* V,double*X)
{
    *T=PS2TREG3(P,S);
    *H=PT2HREG3(P,*T);
    *V=PT2VREG3NT(P,*T);
    *X=1.0;
}


double HS2TREG3(double H,double S)
{
double  Temp_MIN, Temp_MAX,TOL;
double    SH,SM,SL;

    Temp_MAX = PH2TREG3(100, H);

    if (H <= HC_Water )
        Temp_MIN = H2TLREG56(H);
    else if (H<= T2HGREG56(T350C) )
        Temp_MIN = H2TGHREG56(H);
    else
    {
        if ((H>=2610.68079) && (H<=2629.04908) )
        {
            SH=H2SHPBOUND23(H);
            SM=H2SMPBOUND23(H);
            SL=H2SLPBOUND23(H);

            if ((S>=SM*(1-10*DeltaVal)) && (S<=SL*(1-10*DeltaVal)) )
            {
                Temp_MAX=H2TMPBOUND23(H);
                Temp_MIN=H2TLPBOUND23(H);
            }
            else if (S<=SH*(1-10*DeltaVal)) 
            {
                Temp_MIN=H2THPBOUND23(H);
			}
            else
                Temp_MIN=H2TBOUND23(H);
        }
        else 
            Temp_MIN=H2TBOUND23(H);
    }

    TOL =1E-8*(Temp_MIN+Temp_MAX)/2.0;
    return ZBRENT2_KK(Temp_MIN, Temp_MAX, TOL, H, S, 308);
}


double HS2PREG3(double H,double S)
{
    return TH2PREG3(HS2TREG3(H,S),H);
}

double HS2VREG3(double  H,double S)
{
    return TH2VREG3(HS2TREG3(H,S),H);
}

void HSREG3(double* P,double*T, double H,double S, double* V,double*X)
{
    *T=HS2TREG3(H,S);
    *P=TH2PREG3(*T,H);
    *V=TH2VREG3(*T,H);
    *X=1.0;
}

/////////////////////////////////////////////SubRangeBy _X

void SUBRANGEBYP(double  Pressure, int * SUBRANGE)
{
    if ((Pressure >= 6.108E-4-DeltaVal) && (Pressure<=16.535+DeltaVal) )
        *SUBRANGE = 6;
    else if  (Pressure <= PC_Water+DeltaVal) 
        *SUBRANGE = 5;
    else
        *SUBRANGE = 0;
}

void SUBRANGEBYT(double  Temperature , int*  SUBRANGE)
{
    if ((Temperature >= T000C-DeltaVal) && (Temperature<=T350C+DeltaVal) )
        *SUBRANGE = 6;
    else if (Temperature <= TC_Water+DeltaVal) 
        *SUBRANGE = 5;
    else
        *SUBRANGE = 0;
}

void SUBRANGEBYPT(double  Pressure,double  Temperature, int*  SUBRANGE)
{

    *SUBRANGE=0;

    if ((Temperature<T000C-DeltaVal)  || (Temperature>800+DeltaVal) || (Pressure<=6.108E-4-DeltaVal) || (Pressure>100.0+DeltaVal)  )
	{
		//MessageBox("参数越界");
		return ;
		//exit(1);		
	}

    if (Temperature <= TC_Water+DeltaVal) 
    {
        if (fabs(1-T2P_KK(Temperature) / Pressure) < 1E-5 )
        {

            if (Temperature<T350C+DeltaVal) 
                *SUBRANGE = 6;
            else
                *SUBRANGE = 5;
            //exit(1);
			return;
        }

        if (Temperature <= T350C+DeltaVal)  
        {
            if  (Pressure>=T2P_KK(Temperature)-DeltaVal) 
                *SUBRANGE = 1;
            else
                *SUBRANGE = 2;
        }
        else
        {
             if (Pressure>=T2P_KK(Temperature)-DeltaVal) 
                 *SUBRANGE = 4;
             else if (Pressure>T2PBOUND23(Temperature)+DeltaVal) 
                 *SUBRANGE = 3;
             else
                 *SUBRANGE = 2;
        }
    }
    else if (Temperature <= 590) 
    {
        if (Pressure >= T2PBOUND23(Temperature)+DeltaVal) 
            *SUBRANGE = 3;
        else
            *SUBRANGE = 2;
    }
    else
        *SUBRANGE = 2;
}



void SUBRANGEBYPH(double  Pressure , double Enthalpy , int* SUBRANGE )
{

    *SUBRANGE = 0;
    
    if ((Pressure < 6.108E-4-DeltaVal) || (Pressure > 100.0+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ( (Enthalpy < -0.05-DeltaVal) || (Enthalpy > 4158.8+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}


    if (Pressure >= PC_Water )
    {
        if (Enthalpy > (PT2HREG2(Pressure, 800)+DeltaVal) )  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Enthalpy < (PT2HREG1(Pressure, T000C)-DeltaVal) ) 
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Enthalpy <= (PT2HREG1(Pressure, T350C)+DeltaVal) )
            *SUBRANGE = 1;
        else 
		{
		if (Enthalpy >= PT2HREG2(Pressure, P2TBOUND23(Pressure))-DeltaVal )
            *SUBRANGE = 2;
        else if (Enthalpy >= PT2HREG4(Pressure, TC_Water)-DeltaVal )
            *SUBRANGE = 3;
        else
            *SUBRANGE = 4;
		}
    }
    else if (Pressure > T2P_KK(T350C) )
    {
        if (Enthalpy > PT2HREG2(Pressure, 800)+DeltaVal  )  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Enthalpy < PT2HREG1(Pressure, T000C)-DeltaVal )  
		{
		  //MessageBox(NULL,"参数越界","ABC",MB_OK);
		  return;
		  //exit(1);
		}
        if (Enthalpy <= PT2HREG1(Pressure, T350C)+DeltaVal) 
            *SUBRANGE = 1;
        else 
		{
		if (Enthalpy >= PT2HREG2(Pressure, P2TBOUND23(Pressure))-DeltaVal )
            *SUBRANGE = 2;
        else if (Enthalpy >= P2HGREG56(Pressure)-DeltaVal)  
            *SUBRANGE = 3;
        else if (Enthalpy <= P2HLREG56(Pressure)+ DeltaVal)  
            *SUBRANGE = 4;
        else
            *SUBRANGE = 5;
		}
    }
    else
    {
        if (Enthalpy > PT2HREG2(Pressure, 800)+DeltaVal ) 
		{
		   //MessageBox("参数越界");
		   return;
		  //exit(1);
		}
        if (Enthalpy < PT2HREG1(Pressure, T000C)-DeltaVal ) 
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}

        if (Enthalpy >= P2HGREG56(Pressure)-DeltaVal )
            *SUBRANGE = 2;
        else if (Enthalpy < P2HLREG56(Pressure)+DeltaVal )
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
    }
}

//

void SUBRANGEBYPS(double  Pressure ,double  Entropy , int * SUBRANGE )
{

    *SUBRANGE = 0;

    if ((Pressure < 6.108E-4-DeltaVal) || (Pressure > 100.0+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ((Entropy < -0.0069-DeltaVal) || (Entropy > 11.70+DeltaVal))
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}


    if (Pressure >= PC_Water+DeltaVal) 
    {
        if (Entropy > PT2SREG2(Pressure, 800)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Entropy < PT2SREG1(Pressure, T000C)-DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Entropy <= PT2SREG1(Pressure, T350C)+DeltaVal) 
            *SUBRANGE = 1;
        else if (Entropy >= PT2SREG2(Pressure, P2TBOUND23(Pressure))-DeltaVal) 
            *SUBRANGE = 2;
        else if (Entropy <= PT2SREG3(Pressure, TC_Water)-DeltaVal) 
            *SUBRANGE = 4;
        else
            *SUBRANGE = 3;
    }
    else if (Pressure > T2P_KK(T350C)+DeltaVal) 
    {
        if (Entropy > PT2SREG2(Pressure, 800)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Entropy < PT2SREG1(Pressure, T000C)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Entropy <= PT2SREG1(Pressure, T350C)+DeltaVal) 
            *SUBRANGE = 1;
        else 
		{
		if (Entropy >= PT2SREG2(Pressure, P2TBOUND23(Pressure))-DeltaVal) 
            *SUBRANGE = 2;
        else if (Entropy >=P2SGREG56(Pressure)-DeltaVal)
            *SUBRANGE = 3;
        else if (Entropy <= P2SLREG56(Pressure)+DeltaVal) 
            *SUBRANGE = 4;
        else
            *SUBRANGE = 5;
		}
    }
    else
    {
        if (Entropy > PT2SREG2(Pressure, 800)+DeltaVal) 
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Entropy < PT2SREG1(Pressure, T000C)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Entropy >= P2SGREG56(Pressure)-DeltaVal) 
            *SUBRANGE = 2;
        else if (Entropy < P2SLREG56(Pressure)+DeltaVal) 
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
    }
}


void SUBRANGEBYPV(double  Pressure ,double  Volume ,int * SUBRANGE )
{

    *SUBRANGE = 0;

    if ((Pressure < 6.108E-4-DeltaVal) || (Pressure > 100.0+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ((Volume < 0.0009564-DeltaVal) || (Volume > 495.269+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}




    if (Pressure >= PC_Water+DeltaVal) 
    {
        if (Volume > PT2VREG2(Pressure, 800)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Volume < PT2VREG1(Pressure, T000C)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Volume <= PT2VREG1(Pressure, T350C)+DeltaVal )
            *SUBRANGE = 1;
        else 
		{
		if (Volume >= PT2VREG2(Pressure, P2TBOUND23(Pressure))-DeltaVal )
            *SUBRANGE = 2;
        else 
		{if (Volume >= PT2VREG3NT(Pressure, TC_Water)+DeltaVal) 
            *SUBRANGE = 3;
        else
            *SUBRANGE = 4;
		}
		}
    }
    else if (Pressure >T2P_KK(T350C) )
    {
        if (Volume > PT2VREG2(Pressure, 800)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Volume < PT2VREG1(Pressure, T000C)-DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Volume <= PT2VREG1(Pressure, T350C)+DeltaVal) 
            *SUBRANGE = 1;
        else 
		{if (Volume >= PT2VREG2(Pressure, P2TBOUND23(Pressure))-DeltaVal) 
            *SUBRANGE = 2;
        else 
		{
		if (Volume >= P2VGREG56(Pressure)-DeltaVal) 
            *SUBRANGE = 3;
        else 
		{if (Volume <= P2VLREG56(Pressure)+DeltaVal) 
            *SUBRANGE = 4;
        else
            *SUBRANGE = 5;
		}
		}
		}
    }
    else
    {
        if (Volume > PT2VREG2(Pressure, 800)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Volume < PT2VREG1(Pressure, T000C)*(1-100000*DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if (fabs(1-P2VLREG56(Pressure)/Volume)<1E-8) 
        {
            *SUBRANGE = 6;
            return;
        }

        if (Volume >= P2VGREG56(Pressure)-DeltaVal) 
            *SUBRANGE = 2;
        else if (Volume < P2VLREG56(Pressure)*(1+100000*DeltaVal)) 
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
    }
}

void SUBRANGEBYPX(double  P,double X, int * SUBRANGE )
{
    if ((X<0.0-DeltaVal) || (X>1.0+DeltaVal) )
    {
        *SUBRANGE=0;
        return;
    }
    if (P<0.0006108-DeltaVal ) 
        *SUBRANGE=0;
    else if (P<T2P_KK(T350C)+DeltaVal) 
        *SUBRANGE=6;
    else if (P<PC_Water+DeltaVal) 
        *SUBRANGE=5;
    else
        *SUBRANGE=0;
}

void SUBRANGEBYTX(double  T,double X, int * SUBRANGE )
{
    if ((X<0.0-DeltaVal) || (X>1.0+DeltaVal) )
    {
        *SUBRANGE=0;
        return;
    }
    if (T<0.0-DeltaVal)  
        SUBRANGE=0;
    else if (T<T350C+DeltaVal) 
        *SUBRANGE=6;
    else if (T<TC_Water+DeltaVal) 
        *SUBRANGE=5;
    else
        *SUBRANGE=0;
}

void SUBRANGEBYHX(double  H,double X, int * SUBRANGE )
{
double    HMIN,HMAX;

    if ((X<0.0-DeltaVal) || (X>1.0+DeltaVal) )
    {
        *SUBRANGE=0;
        return;
    }

    HMIN=TX2HREG56(T000C,X);
    if (HMIN>HC_Water)  HMIN=HC_Water;

    if (H<HMIN-DeltaVal) 
    {
        *SUBRANGE=0;
        return;
    }

    HMAX=HMAXXREG56(X);
    if (H>HMAX+DeltaVal) 
    {
        *SUBRANGE=0;
        return;
    }

    if (HMAX>HC_Water )
    {
        if (H<TX2VREG56(T350C,X)+DeltaVal )
            *SUBRANGE=5;
        else
           *SUBRANGE=6;
    }
    else
    {
        if (H<TX2VREG56(T350C,X)+DeltaVal )
            *SUBRANGE=6;
        else
            *SUBRANGE=5;
    }
}

void SUBRANGEBYSX(double  S,double X,int * SUBRANGE )
{
double   SMIN,SMAX,S350,S349;

    if ((X<0.0-DeltaVal) || (X>1.0+DeltaVal)) 
    {
        *SUBRANGE=0;
        return;
    }
	
    SMIN=SMINXREG56(X);
    
	if (S<SMIN-DeltaVal) 
    {
        *SUBRANGE=0;
        return;
    }
    SMAX=SMAXXREG56(X);
    if (S>SMAX+DeltaVal )
    {
        *SUBRANGE=0;
        return;
    }


    S350=TX2SREG56(T350C,X);
    S349=TX2SREG56(349,X);

    if (S349<S350 )
    {
        if (S<S350-DeltaVal) 
        {
            *SUBRANGE=6;
			return;
        }
        else
        {
            *SUBRANGE=5;
			return;
        }
    }
    else
    {
        if (S<S350-DeltaVal) 
        {
            *SUBRANGE=5;
			return;
        }
        else
        {
            *SUBRANGE=6;
			return;
        }
    }

}

void SUBRANGEBYVG_KK(double  VG ,  int *SUBRANGE )
{
    *SUBRANGE = 0;

    if ((VG < VC_Water-DeltaVal) || (VG > 206.305+DeltaVal)) 
         return;

    if ( VG<0.008799 )
        *SUBRANGE=5;
    else
        *SUBRANGE=6;

}

void SUBRANGEBYVX(double  V,double X,int * SUBRANGE )
{
 double   VMIN,VMAX;

    if ((X<0.0-DeltaVal) || (X>1.0+DeltaVal)) 
    {
        *SUBRANGE=0;
        return;
    }

    VMIN=VMINXREG56(X);
    if (V<VMIN-DeltaVal) 
    {
        *SUBRANGE=0;
        return;
    }

    VMAX=TX2VREG56(T000C,X);
    if (VMAX<VC_Water)  VMAX=VC_Water;
    if (V>VMAX+DeltaVal) 
    {
        *SUBRANGE=0;
        return;
    }

    if (V<TX2VREG56(T350C,X)-DeltaVal )
    {
        *SUBRANGE=5;return;
    }
    else
    {
        *SUBRANGE=6;return;
    }
}


double P_HMINTREG1(double PMIN,double PMAX,double T)
{
double    P_MIN,P_MAX,P_Mid;

    if (fabs(PMAX-PMIN)<10*DeltaVal) 
    {
        return (PMIN+PMAX)/2;
        //exit(1);
    }

    P_MIN=PMIN;
    P_MAX=PMAX;
    P_Mid=(P_MIN+P_MAX)/2;
    if (PT2HREG1(P_Mid-2*DeltaVal,T)<PT2HREG1(P_Mid+2*DeltaVal,T)) 
    {
        P_MAX=P_Mid-2*DeltaVal;
    }
    else
    {
        P_MIN=P_Mid-2*DeltaVal;
    }

    return P_HMINTREG1(P_MIN,P_MAX,T);
}


double PHMINTREG1(double T)
{
double    P_MIN,P_MAX;

    P_MIN=T2P_KK(T);
    P_MAX=100;
    return P_HMINTREG1(P_MIN,P_MAX,T);
}


double HMINTREG1(double T)
{
    return PT2HREG1(PHMINTREG1(T),T);
}

void HMINTREG167(double T, double H)
{
    H=HMINTREG1(T);
}


double P_HMAXTREG1(double PMIN,double PMAX,double T)
{
 double   P_MIN,P_MAX,P_Mid;

    if (fabs(PMAX-PMIN)<10*DeltaVal) 
    {
        return (PMIN+PMAX)/2;
        //exit(1);
    }

    P_MIN=PMIN;
    P_MAX=PMAX;
    P_Mid=(P_MIN+P_MAX)/2;
    if (PT2HREG1(P_Mid-2*DeltaVal,T)<PT2HREG1(P_Mid+2*DeltaVal,T)) 
    {
        P_MIN=P_Mid-2*DeltaVal;

    }
    else
    {
        P_MAX=P_Mid-2*DeltaVal;
    }

    return P_HMAXTREG1(P_MIN,P_MAX,T);
}


double PHMAXTREG1(double T)
{
double    P_MIN,P_MAX;

    P_MIN=T2P_KK(T);
    P_MAX=100;
    return P_HMAXTREG1(P_MIN,P_MAX,T);
}


double HMAXTREG1(double T)
{
    return PT2HREG1(PHMAXTREG1(T),T);
}


void SUBRANGEBYTH(double  Temperature ,double Enthalpy , int * SUBRANGE ,bool* B16)
{
 double   H_MINSAT , H_MINP100,H_MIN ;


    *SUBRANGE = 0;
    *B16=false;

    if ((Temperature < T000C-DeltaVal) || (Temperature > 800+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ((Enthalpy < -0.05-DeltaVal) || (Enthalpy > 4158.8+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}


    if (Temperature >= 590-DeltaVal) 
    {
        if (Enthalpy > PT2HREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy < PT2HREG2(100, Temperature)-DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        *SUBRANGE = 2;
    }
    else if (Temperature >= TC_Water-DeltaVal) 
    {
        if (Enthalpy > PT2HREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy < PT2HREG3(100, Temperature)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy <= PT2HREG2(T2PBOUND23(Temperature), Temperature)-DeltaVal) 
            *SUBRANGE = 3;
        else
            *SUBRANGE = 2;
    }
    else if (Temperature > T350C+DeltaVal) 
    {
        if (Enthalpy > PT2HREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Enthalpy < PT2HREG4(100, Temperature)-DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}

        if (Enthalpy >= PT2HREG2(T2PBOUND23(Temperature), Temperature)-DeltaVal) 
            *SUBRANGE = 2;
        else 
		{
		if (Enthalpy >= T2HGREG56(Temperature)) 
            *SUBRANGE = 3;
        else if (Enthalpy < T2HLREG56(Temperature))
            *SUBRANGE = 4;
        else
            *SUBRANGE = 5;
		}
    }
    else if (Temperature > 336+DeltaVal) 
    {

        if (Enthalpy > PT2HREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Enthalpy < PT2HREG4(100, Temperature)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if (Enthalpy >= T2HGREG56(Temperature)*(1+DeltaVal)) 
            *SUBRANGE = 2;
        else if (Enthalpy < T2HLREG56(Temperature)) 
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
    }
    else if (Temperature > 247+DeltaVal )
    {

        if (Enthalpy > PT2HREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy >= T2HGREG56(Temperature) )
        {
            *SUBRANGE = 2;
            return;
			//exit(1);
        }

        H_MIN=HMINTREG1(Temperature);
        if (Enthalpy < H_MIN*(1-10*DeltaVal))  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}


        H_MINSAT = T2HLREG56(Temperature);
        H_MINP100 = PT2HREG1(100, Temperature);

        if (H_MINP100>H_MINSAT) 
        {

            if (Enthalpy>H_MINP100) 
            {
                *SUBRANGE = 6;
            }
            else if (Enthalpy>H_MINSAT) 
            {
                *SUBRANGE = 6;
                *B16=true;
            }
            else
            {
                *SUBRANGE = 1;
            }
        }
        else
        {

            if (Enthalpy>H_MINSAT) 
            {
                *SUBRANGE = 6;
            }
            else if (Enthalpy>H_MINP100) 
            {
                *SUBRANGE = 1;
            }
            else
            {
                *SUBRANGE = 1;
            }
        }
    }
    else
    {
        if (Enthalpy > PT2HREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy >= T2HGREG56(Temperature)*(1+DeltaVal)) 
        {
            *SUBRANGE = 2;
            return;
			//exit(1);
        }

        H_MINSAT = T2HLREG56(Temperature);
        H_MIN=H_MINSAT;
        if (Enthalpy < H_MIN*(1-10*DeltaVal) )  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        H_MINP100 = PT2HREG1(100, Temperature);

        if (Enthalpy < H_MINP100*(1+DeltaVal)) 
        {
            *SUBRANGE = 6;
            *B16=true;
        }
        else
        {
            *SUBRANGE = 6;
        }
    }
}



void SUBRANGEBYTS(double  Temperature , double Entropy ,  int * SUBRANGE )
{
double    SMAX;


    *SUBRANGE = 0;

    if ((Temperature < T000C-DeltaVal) || (Temperature > 800+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ((Entropy < -0.0069-DeltaVal) || (Entropy > 11.70+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}


    if (Temperature >= 590-DeltaVal )
    {
        if (Entropy > PT2SREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Entropy < PT2SREG2(100, Temperature)-DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        *SUBRANGE = 2;
    }
    else if (Temperature >= TC_Water+DeltaVal )
    {
        if (Entropy > PT2SREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Entropy < PT2SREG3(100, Temperature)-DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Entropy <= PT2SREG2(T2PBOUND23(Temperature), Temperature)-DeltaVal) 
            *SUBRANGE = 3;
        else
            *SUBRANGE = 2;
    }
    else if (Temperature > T350C+DeltaVal )
    {
        if (Entropy > PT2SREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Entropy < PT2SREG4(100, Temperature)-DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}

        if (Entropy >= PT2SREG2(T2PBOUND23(Temperature), Temperature)-DeltaVal) 
            *SUBRANGE = 2;
        else if (Entropy >= T2SGREG56(Temperature)) 
            *SUBRANGE = 3;
        else if (Entropy < T2SLREG56(Temperature)) 
            *SUBRANGE = 4;
        else
            *SUBRANGE = 5;
    }
    else
    {
        if (Entropy < PT2SREG1(100, Temperature)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Entropy > PT2SREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}


        if (Temperature>5 )
        {
            if (Entropy >= T2SGREG56(Temperature)-DeltaVal) 
                *SUBRANGE = 2;
            else if (Entropy <T2SLREG56(Temperature) +DeltaVal) 
                *SUBRANGE = 1;
            else
                *SUBRANGE = 6;
        }
        else
        {
            SMAX=PT2SREG1(PSMAXTREG1(Temperature),Temperature);
            if (Entropy >= T2SGREG56(Temperature)-DeltaVal )
                *SUBRANGE = 2;
            else if (Entropy < SMAX+DeltaVal )
                *SUBRANGE = 1;
            else
                *SUBRANGE = 6;
        }
    }
}


void SUBRANGEBYTV(double  Temperature , double Volume , int * SUBRANGE )
{

    *SUBRANGE = 0;

    if ((Temperature < T000C-DeltaVal) || (Temperature > 800+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ((Volume < 0.0009564-DeltaVal) || (Volume > 495.269+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if (Temperature >= 590-DeltaVal) 
    {
        if (Volume > PT2VREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Volume < PT2VREG2(100, Temperature)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        *SUBRANGE = 2;
    }
    else if (Temperature >= TC_Water+DeltaVal )
    {
        if (Volume > PT2VREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Volume < PT2VREG3NT(100, Temperature)-DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Volume <= PT2VREG2(T2PBOUND23(Temperature), Temperature)-DeltaVal )
            *SUBRANGE = 3;
        else
            *SUBRANGE = 2;
    }
    else if (Temperature > T350C+DeltaVal) 
    {
        if (Volume > PT2VREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Volume < PT2VREG4NT(100, Temperature)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if (Volume >= PT2VREG2(T2PBOUND23(Temperature), Temperature)-DeltaVal )
            *SUBRANGE = 2;
        else if (Volume >= T2VGREG56(Temperature) )
            *SUBRANGE = 3;
        else if (Volume < T2VLREG56(Temperature) )
            *SUBRANGE = 4;
        else
            *SUBRANGE = 5;
    }
    else
    {
        if (Volume > PT2VREG2(6.108E-4, Temperature)+DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Volume < PT2VREG1(100, Temperature)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if (Volume >= T2VGREG56(Temperature)-DeltaVal) 
            *SUBRANGE = 2;
        else if (Volume < T2VLREG56(Temperature)+DeltaVal )
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
    }
}


void SUBRANGEBYHL_KK(double  HL , int * SUBRANGE )
{
    *SUBRANGE = 0;

    if ((HL < -0.04-DeltaVal) || (HL > HC_Water+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if  (HL<1672.0) 
        *SUBRANGE=6;
    else
        *SUBRANGE=5;

}


void SUBRANGEBYHS(double  Enthalpy , double Entropy ,  int * SUBRANGE )
{
 double   SH,SM,SL;


    *SUBRANGE = 0;

    if ((Enthalpy < -0.05-DeltaVal) || (Enthalpy > 4158.8+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ((Entropy < -0.0069-DeltaVal) || (Entropy > 11.70+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}



    if (Entropy > P2SGREG56(6.108E-4)*(1-DeltaVal) )
    {
        if (Enthalpy < PS2HREG2(6.108E-4, Entropy)*(1-100*DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy > TS2HREG2(800, Entropy)*(1+DeltaVal))  
		{
		    //MessageBox("参数越界");
		    return;
		    //exit(1);
		}
        *SUBRANGE = 2;
    }
    else if (Entropy >= PT2SREG2(100, 800)*(1+DeltaVal) )
    {
        if (Enthalpy < PS2HREG56(6.108E-4, Entropy)*(1-DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy > TS2HREG2(800, Entropy)*(1+DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if (Enthalpy > S2HGREG56(Entropy)*(1-DeltaVal)) 
            *SUBRANGE = 2;
        else
            *SUBRANGE = 6;
    }
    else if (Entropy >= SC_Water )
    {
        if (Enthalpy < PS2HREG56(6.108E-4, Entropy)*(1-DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy > PT2HREG2(100,800)*(1+DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if (Enthalpy <= S2HGREG56(Entropy)*(1-DeltaVal)) 
        {
            if (Entropy>5.2176-DeltaVal) 
                *SUBRANGE = 6;
            else if (Enthalpy <= TS2HREG56(T350C,Entropy)*(1+DeltaVal) )
                *SUBRANGE = 6;
            else
                *SUBRANGE = 5;
        }
        else
        {

            if (Enthalpy <= PS2HREG3(100, SC_Water)*(1-DeltaVal) )
                *SUBRANGE = 3;
            else if (Enthalpy <= T2HGREG56(T350C)*(1-DeltaVal)) 
            {
                if (Entropy < PH2SREG3(100, Enthalpy)*(1-10*DeltaVal))  
				{
		          //MessageBox("参数越界");
		          return;
		          //exit(1);
				}
                *SUBRANGE = 3;
            }                      
            else if (Enthalpy <=T2HBOUND23(590) *(1-DeltaVal)) 
            {
                if (Entropy < PH2SREG3(100, Enthalpy)*(1-10*DeltaVal))  
				{
		          //MessageBox("参数越界");
		          return;
		          //exit(1);
				}

                if ((Enthalpy<=2610.68079) || (Enthalpy>=2629.04908)) 
                {
                    if (Entropy >= H2SBOUND23(Enthalpy)*(1-10*DeltaVal)) 
                        *SUBRANGE = 2;
                    else
                        *SUBRANGE = 3;
                }
                else
                {
                    SH=H2SHPBOUND23(Enthalpy);
                    SM=H2SMPBOUND23(Enthalpy);
                    SL=H2SLPBOUND23(Enthalpy);

                    if ((Entropy>=SH*(1-10*DeltaVal)) && (Entropy<=SM*(1-10*DeltaVal)) || (Entropy>=SL*(1-10*DeltaVal)) )
                        *SUBRANGE = 2;
                    else
                        *SUBRANGE = 3;
				}
            }
            else
            {
                if (Entropy < PH2SREG2(100, Enthalpy)*(1-10*DeltaVal))  
				{
		          //MessageBox("参数越界");
		          return;
		          //exit(1);
				}

                *SUBRANGE = 2;
            }
        }
    }
    else if (Entropy >= T2SLREG56(T350C)+DeltaVal )
    {
        if (Enthalpy < PS2HREG56(6.108E-4, Entropy)*(1-DeltaVal))  
		{
		    //MessageBox("参数越界");
		    return;
		    //exit(1);
		}
        if (Enthalpy > PS2HREG3(100, Entropy)*(1+DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		  //exit(1);
		}

        if (Enthalpy > TS2HREG3(TC_Water,Entropy)*(1-DeltaVal) )
            *SUBRANGE = 3;
        else 
		{
		if (Enthalpy > S2HLREG56(Entropy)*(1+DeltaVal)) 
            *SUBRANGE = 4;
        else 
		{ if (Enthalpy > TS2HREG56(T350C,Entropy)*(1+DeltaVal)) 
            *SUBRANGE = 5;
         else
            *SUBRANGE = 6;
		}
		}
    }
    else if (Entropy >= PT2SREG3(100, TC_Water) )
    {
        if (Enthalpy < PS2HREG56(6.108E-4, Entropy)*(1-DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Enthalpy > PS2HREG3(100, Entropy)*(1+DeltaVal))  
		{
		    //MessageBox("参数越界");
		    return;
		    //exit(1);
		}

        if (Enthalpy > TS2HREG3(TC_Water, Entropy)*(1-DeltaVal)) 
            *SUBRANGE = 3;
        else 
		{
		if (Enthalpy > TS2HREG1(T350C, Entropy)*(1+DeltaVal)) 
            *SUBRANGE = 4;
        else 
		{if (Enthalpy > S2HLREG56(Entropy)*(1-DeltaVal)) 
            *SUBRANGE = 1;
         else
            *SUBRANGE = 6;
		}
		}
    }
    else if (Entropy >= PT2SREG1(100, T350C)) 
    {
        if (Enthalpy < PS2HREG56(6.108E-4, Entropy)*(1-DeltaVal))  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Enthalpy > PS2HREG4(100, Entropy)*(1+DeltaVal))  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if (Enthalpy > TS2HREG1(T350C, Entropy)*(1+DeltaVal)) 
            *SUBRANGE = 4;
        else
		{
		if (Enthalpy > S2HLREG56(Entropy)*(1-DeltaVal) )
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
		}
    }
    else
    {
        if (Enthalpy < PS2HREG56(6.108E-4, Entropy)-10*DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        if (Enthalpy > PS2HREG1(100, Entropy)+10*DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}

        if (Enthalpy > S2HLREG56(Entropy)*(1-DeltaVal))
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
    }
}


void SUBRANGEBYHV(double  Enthalpy ,double  Volume , int *SUBRANGE )
{

    *SUBRANGE = 0;

    if ((Enthalpy < -0.05-DeltaVal) || (Enthalpy > 4158.8+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ((Volume < 0.0009564-DeltaVal) || (Volume > 495.269+DeltaVal)) 
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}


    if (Enthalpy > PT2HREG2(100, 800) )
    {
        if (Volume < TH2VREG2(800, Enthalpy)*(1-1000*DeltaVal))  
		{
		  //MessageBox("参数越界");
		  return;
		   //exit(1);
		}
        if (Volume > PH2VREG2(6.108E-4, Enthalpy)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        *SUBRANGE = 2;
    }
    else if (Enthalpy >T2HBOUND23(590) )
    {
        if (Volume < PH2VREG2(100, Enthalpy)-DeltaVal)  
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if (Volume > PH2VREG2(6.108E-4, Enthalpy)+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}
        *SUBRANGE = 2;
    }
    else
    {


        if (Volume > PH2VREG2(6.108E-4, T2HBOUND23(590))+DeltaVal)  
		{
		  //MessageBox("参数越界");
		  return;
		  //exit(1);
		}

        if (Volume > P2VGREG56(6.108E-4) )
        {
            if (Enthalpy < PV2HREG2(6.108E-4, Volume)-DeltaVal)  
			{
		      //MessageBox("参数越界");
		      return;
		      //exit(1);
			}
            *SUBRANGE = 2;
        }
        else if (Volume > T2VGREG56(T350C)) 
        {
            if (Enthalpy < PV2HREG56(6.108E-4, Volume)-DeltaVal)  
			{
		      //MessageBox("参数越界");
		      return;
		      //exit(1);
			}
            
			if (Enthalpy > V2HGREG56(Volume)) 
                *SUBRANGE = 2;
            else
                *SUBRANGE = 6;
        }
        else if (Volume > VC_Water) 
        {
            if (Enthalpy < PV2HREG56(6.108E-4, Volume)-DeltaVal)  
			{
		      //MessageBox("参数越界");
		      return;
		      //exit(1);
			}

            if (Enthalpy > V2HBOUND23(Volume)) 
                *SUBRANGE = 2;
            else 
			{
			if (Enthalpy > V2HGREG56(Volume)) 
                *SUBRANGE = 3;
            else 
			{
				if (Enthalpy>TV2HREG56(T350C,Volume)) 
                    *SUBRANGE = 5;
                else
                   *SUBRANGE = 6;
			}
			}
        }
        else if (Volume > PT2VREG2(100, 590)) 
        {
            if (Enthalpy < PV2HREG56(6.108E-4, Volume)-DeltaVal)  
			{
		      //MessageBox("参数越界");
		      return;
		      //exit(1);
			}

            if (Enthalpy > V2HBOUND23(Volume)) 
                *SUBRANGE = 2;
            else 
			{
			if (Enthalpy >TV2HREG3(TC_Water,Volume)) 
                *SUBRANGE = 3;
            else 
			{if (Enthalpy > V2HLREG56(Volume) )
                *SUBRANGE = 4;
            else 
			{if (Enthalpy>TV2HREG56(T350C,Volume)) 
                *SUBRANGE = 5;
            else
                *SUBRANGE = 6;
			}
			}
			}
        }
        else if (Volume > T2HLREG56(T350C)) 
        {
            if (Enthalpy < PV2HREG56(6.108E-4, Volume)-DeltaVal)  
			{
		       //MessageBox("参数越界");
		       return;
		       //exit(1);
			}
            if (Enthalpy > PV2HREG3(100, Volume)*(1+100*DeltaVal))  
			{
		      //MessageBox("参数越界");
		      return;
		      //exit(1);
			}

            if (Enthalpy >TV2HREG3(TC_Water,Volume)) 
                *SUBRANGE = 3;
            else 
			{
			if (Enthalpy > V2HLREG56(Volume)) 
                *SUBRANGE = 4;
            else 
			{if (Enthalpy>TV2HREG56(T350C,Volume) )
                *SUBRANGE = 5;
            else
                *SUBRANGE = 6;
			}
			}
        }
        else if (Volume > PT2VREG1(100, TC_Water)) 
        {
            if (Enthalpy < PV2HREG56(6.108E-4, Volume)-DeltaVal)  
			{
		      //MessageBox("参数越界");
		      return;
		      //exit(1);
			}
            if (Enthalpy > PV2HREG3(100, Volume)*(1+100*DeltaVal))  
			{
		      //MessageBox("参数越界");
		      return;
		      //exit(1);
			}

            if (Enthalpy >TV2HREG3(TC_Water,Volume)+DeltaVal) 
                *SUBRANGE = 3;
            else 
			{
			if (Enthalpy > TV2HREG1(T350C, Volume)*(1+10*DeltaVal)) 
                *SUBRANGE = 4;
            else 
			{if (Enthalpy > V2HLREG56(Volume)) 
                *SUBRANGE = 1;
            else
                *SUBRANGE = 6;
			}
			}
        }
        else if (Volume > PT2VREG1(100, T350C)) 
        {
            if (Enthalpy < PV2HREG56(6.108E-4, Volume)-DeltaVal)  
			{
		       //MessageBox("参数越界");
		       return;
		       //exit(1);
			}
            if (Enthalpy > PV2HREG3(100, Volume)*(1+100*DeltaVal))  
			{
		       //MessageBox("参数越界");
		       return;
		       //exit(1);
			}

            if (Enthalpy > TV2HREG1(T350C, Volume)*(1+10*DeltaVal)) 
                *SUBRANGE = 4;
            else 
			{if (Enthalpy > V2HLREG56(Volume) )
                *SUBRANGE = 1;
            else
                *SUBRANGE = 6;
			}
        }
        else if (Volume > P2VLREG56(6.108E-4)) 
        {
            if (Enthalpy < PV2HREG56(6.108E-4, Volume)-DeltaVal)  
			{
		        //MessageBox("参数越界");
		        return;
		        //exit(1);
			}
            if (Enthalpy > PV2HREG1(100, Volume)*(1+100*DeltaVal))  
			{
		       //MessageBox("参数越界");
		       return;
		       //exit(1);
			}
            if (Enthalpy > V2HLREG56(Volume) )
                *SUBRANGE = 1;
            else
                *SUBRANGE = 6;
        }
        else
        {
            if (Enthalpy < PV2HREG1(6.108E-4, Volume)-DeltaVal)  
			{
		       //MessageBox("参数越界");
		        return;
		       //exit(1);
			}
            if (Enthalpy > PV2HREG1(100, Volume)*(1+100*DeltaVal))  
			{
		       //MessageBox("参数越界");
		       return;
		      //exit(1);
			}
            *SUBRANGE = 1;
        }
    }
}

//

void SUBRANGEBYSG_KK(double  SG ,  int * SUBRANGE )
{

    *SUBRANGE = 0;
    if ( (SG < SC_Water-DeltaVal) || (SG > 9.1577+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}
    if ( SG<=5.2177 )
        *SUBRANGE=5;
    else
        *SUBRANGE=6;
}


void SUBRANGEBYSV(double  Entropy ,double  Volume , int * SUBRANGE )
{

    *SUBRANGE = 0;

    if ( (Entropy < -0.0069-DeltaVal) || (Entropy > 11.70+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}

    if ( (Volume < 0.0009564-DeltaVal) || (Volume > 495.269+DeltaVal) )
    {
		//MessageBox("参数越界");
		return;
		//exit(1);
	}


    if ( Volume > P2VGREG56(6.108E-4) )
    {
        if ( Entropy > TV2SREG2(800, Volume)*(1+1000*DeltaVal) ) 
		{
		    //MessageBox("参数越界");
		    return;
		    //exit(1);
		}
        if ( Entropy < PV2SREG2(6.108E-4, Volume)*(1-100*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        *SUBRANGE = 2;
    }
    else if ( Volume > T2VGREG56(T350C) )
    {
        if ( Entropy > TV2SREG2(800, Volume)*(1+100*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( Entropy < PV2SREG56(6.108E-4, Volume)*(1-10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if ( Entropy > V2SGREG56(Volume)+DeltaVal )
            *SUBRANGE = 2;
        else
            *SUBRANGE = 6;
    }
    else if ( Volume > PT2VREG2(100, 800) )
    {
        if ( Entropy > TV2SREG2(800, Volume)*(1+10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( Entropy < PV2SREG56(6.108E-4, Volume)*(1-10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if ( Entropy > V2SBOUND23(Volume) )
            *SUBRANGE = 2;
        else
		{
		if ( Entropy > V2SGREG56(Volume) )
            *SUBRANGE = 3;
        else 
		{if ( Entropy >TV2SREG56(T350C,Volume) )
            *SUBRANGE = 5;
        else
            *SUBRANGE = 6;
		}
		}
    }
    else if ( Volume > VC_Water )
    {
        if ( Entropy > PV2SREG2(100, Volume)*(1+10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( Entropy < PV2SREG56(6.108E-4, Volume)-DeltaVal ) 
		{
		   //MessageBox("参数越界");
	       return;
		   //exit(1);
		}

        if ( Entropy > V2SBOUND23(Volume) )
            *SUBRANGE = 2;
        else 
		{
		if ( Entropy > V2SGREG56(Volume) )
            *SUBRANGE = 3;
        else 
		{if ( Entropy >TV2SREG56(T350C,Volume) )
            *SUBRANGE = 5;
        else
            *SUBRANGE = 6;
		}
		}
    }
    else if ( Volume > PT2VREG2(100,590)-DeltaVal )
    {
        if ( Entropy > PV2SREG2(100, Volume)*(1+10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( Entropy < PV2SREG56(6.108E-4, Volume)-DeltaVal ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if ( Entropy > V2SBOUND23(Volume) )
            *SUBRANGE = 2;
        else 
		{
		if ( Entropy > TV2SREG3(TC_Water,Volume) )
            *SUBRANGE = 3;
        else 
		{
			if ( Entropy > V2SLREG56(Volume) )
            *SUBRANGE = 4;
        else 
		{
			if ( Entropy >TV2SREG56(T350C,Volume) )
            *SUBRANGE = 5;
           else
            *SUBRANGE = 6;
		}
		}
		}
    }
    else if ( Volume > T2VLREG56(T350C)+DeltaVal )
    {
        if ( Entropy > PV2SREG3(100, Volume)*(1+10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( Entropy < PV2SREG56(6.108E-4, Volume)-DeltaVal ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if ( Entropy > TV2SREG3(TC_Water,Volume)-DeltaVal )
            *SUBRANGE = 3;
        else 
		{
		if ( Entropy > V2SLREG56(Volume)*(1+10*DeltaVal) )
            *SUBRANGE = 4;
        else 
		{
			if ( Entropy >TV2SREG56(T350C,Volume)*(1+10*DeltaVal) )
            *SUBRANGE = 5;
            else
            *SUBRANGE = 6;
		}
		}
    }
    else if ( Volume > PT2VREG1(100, TC_Water) )
    {
        if ( Entropy > PV2SREG3(100, Volume)*(1+10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( Entropy < PV2SREG56(6.108E-4, Volume)-DeltaVal ) 
		{
		   //MessageBox("参数越界");
		   return;
		    //exit(1);
		}

        if ( Entropy > TV2SREG3(TC_Water,Volume)-DeltaVal )
            *SUBRANGE = 3;
        else 
		{
		if ( Entropy > TV2SREG1(T350C,Volume)+DeltaVal )
            *SUBRANGE = 4;
        else 
		{
		if ( Entropy > V2SLREG56(Volume)*(1-10*DeltaVal) )
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
		}
		}
    }
    else if ( Volume > PT2VREG1(100, T350C) )
    {
        if ( Entropy > PV2SREG4(100, Volume)*(1+10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( Entropy < PV2SREG56(6.108E-4, Volume)-DeltaVal ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if ( Entropy > TV2SREG1(T350C,Volume)+DeltaVal )
            *SUBRANGE = 4;
        else 
		{
			if ( Entropy > V2SLREG56(Volume)*(1-10*DeltaVal) )
            *SUBRANGE = 1;
            else
            *SUBRANGE = 6;
		}
    }
    else if ( Volume > P2VLREG56(6.108E-4) )
    {
        if ( Entropy > PV2SREG1(100, Volume)*(1+10*DeltaVal) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( Entropy < PV2SREG56(6.108E-4, Volume)-DeltaVal ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}

        if ( Entropy > V2SLREG56(Volume)*(1-10*DeltaVal) )
            *SUBRANGE = 1;
        else
            *SUBRANGE = 6;
    }
    else
    {
        if ( (Entropy<-0.00683) ) 
		  return;	
			//exit(1);

        if ( (Entropy>0.0) && (Entropy > PV2SREG1(100, Volume)*(1+100*DeltaVal)) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        if ( (Entropy>0.0) && (Entropy < PV2SREG1(6.108E-4, Volume)*(1-10*DeltaVal)) ) 
		{
		   //MessageBox("参数越界");
		   return;
		   //exit(1);
		}
        *SUBRANGE = 1;
    }
}


void PT_KK (double P,double T, double * H,double *S,double* V,double*X, int * RANGE)
{

    SUBRANGEBYPT(P,T,RANGE);
   switch ( *RANGE )
   { 
    case 1:
    {
        PTREG1(P,T,H,S,V,X); break;
    }
    case 2:
    {
        PTREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        PTREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        PTREG4(P,T,H,S,V,X);break;
    }
    case  5:
    {
        PTREG4(P,T,H,S,V,X);break;
    }
    case 6:
    {
        PTREG1(P,T,H,S,V,X);break;
    }
    
   }
}


void PH_KK (double P, double* T,double H, double* S,double* V,double*X,int*  RANGE)
{
    SUBRANGEBYPH(P,H,RANGE);
    switch ( *RANGE )
    {
	case 	1:
    {
        PHREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        PHREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        PHREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        PHREG4(P,T,H,S,V,X);break;
    }
    case 5:
	{
        PHREG56(P,T,H,S,V,X);break;
    }
    case 6: 
    {
        PHREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void PS_KK (double P, double* T,double*H, double S, double* V,double*X, int* RANGE)
{
    SUBRANGEBYPS(P,S,RANGE);
    switch ( *RANGE )
    {
	case  1:
    {
        PSREG1(P,T,H,S,V,X);break;
    }
    case  2:
    {
        PSREG2(P,T,H,S,V,X);break;
    }
    case  3:
    {
        PSREG3(P,T,H,S,V,X);break;
    }
    case  4:
    {
        PSREG4(P,T,H,S,V,X);break;
    }
    case 5:
	{
        PSREG56(P,T,H,S,V,X);break;
    }
	case 6:
    {
        PSREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
	} 
}


void PV_KK (double P, double* T,double*H,double*S, double V , double* X ,int* RANGE )
{
    SUBRANGEBYPV(P,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        PVREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        PVREG2(P,T,H,S,V,X);break;
    }
    case  3:
    {
        PVREG3(P,T,H,S,V,X);break;
    }
    case  4:
    {
        PVREG4(P,T,H,S,V,X);break;
    }
    case 5:
	{
        PVREG56(P,T,H,S,V,X);break;
    }
	case 6:
    {
        PVREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void PX_KK (double P, double*  T,double*H,double*S,double*V, double X, int* RANGE)
{
    SUBRANGEBYPX(P,X,RANGE);
    switch ( *RANGE )
    {
	case	5:
	{
        PXREG56(P,T,H,S,V,X);break;
	}
	case	6:
    {
        PXREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}

/////////////gh1
void  TH_KK (double*  P, double T, double H , double* S,double*V,double*X , int * RANGE)
{
bool    *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
    {
        THREG1(P,T,H,S,V,X); break;
    }
    case  2:
    {
        THREG2(P,T,H,S,V,X); break;
    }
    case  3:
    {
        THREG3(P,T,H,S,V,X); break;
    }
    case  4:
    {
        THREG4(P,T,H,S,V,X); break;
    }
    case 5:
	case 6:
    {
        THREG56(P,T,H,S,V,X); break;
    }
    default:
		{
         *RANGE=0; break;
		}
    }
}


void  THLP_KK (double* P, double T, double H , double* S,double*V,double*X , int* RANGE)
{
bool    *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
    {
        THLPREG1(P,T,H,S,V,X); break;
    }
    case 2:
    {
        THREG2(P,T,H,S,V,X); break;
    }
    case 3:
    {
        THREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        THREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        THREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}



void  THHP_KK (double* P, double T, double H , double* S,double*V,double*X ,int * RANGE)
{
bool    *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=1;
    switch ( *RANGE )
    {
	case 1:
    {
        THHPREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        THREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        THREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        THREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        THREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}



void  TSHP_KK (double* P, double T, double* H, double S, double* V,double*X,int * RANGE)
{
    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case  1:
    {
        TSHPREG1(P,T,H,S,V,X);break;
    }
    case  2:
    {
        TSREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        TSREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        TSREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        TSREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  TSLP_KK (double* P, double T, double* H, double S, double *V,double *X, int * RANGE)
{
    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        TSLPREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        TSREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        TSREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        TSREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        TSREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}



void  TS_KK (double* P, double T, double* H, double S, double *V,double *X, int * RANGE)
{
    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        TSREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        TSREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        TSREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        TSREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        TSREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  TV_KK (double* P, double T, double* H,double*S, double V, double* X, int * RANGE)
{
    SUBRANGEBYTV(T,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        TVREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        TVREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        TVREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        TVREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        TVREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  TX_KK (double* P, double T,  double* H,double*S,double*V , double X, int * RANGE)
{
    SUBRANGEBYTX(T,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        TXREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  HS_KK (double* P,double*T, double H,double S, double* V,double*X, int * RANGE)
{
    SUBRANGEBYHS(H,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        HSREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        HSREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        HSREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        HSREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        HSREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  HV_KK(double* P,double*T, double H, double* S, double V, double* X, int* RANGE)
{
    SUBRANGEBYHV(H,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        HVREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        HVREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        HVREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        HVREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        HVREG56(P,T,H,S,V,X);break;
    }
    default:
       *RANGE=0;
    }
}

void  HX_KK(double* P,double*T, double H, double* S,double*V, double X, int* RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        HXREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  HXLP_KK(double* P,double*T, double H, double* S,double*V,double X, int* RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        HXREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  HXHP_KK(double* P,double*T, double H, double* S,double*V, double X, int* RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        HXHPREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}

void  S2TG_KK(double S, double* T, int* RANGE)
{
    SUBRANGEBYSG_KK(S,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
        *T=S2TGREG56(S);break;
		}
	default:
    {
        *RANGE=0;
        *T=-1.0;
        break;
    }
    }
}

void  SV_KK (double* P,double*T,double*H, double S,double V, double* X, int * RANGE)
{
    SUBRANGEBYSV(S,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        SVREG1(P,T,H,S,V,X);break;
    }
    case 2:
    {
        SVREG2(P,T,H,S,V,X);break;
    }
    case 3:
    {
        SVREG3(P,T,H,S,V,X);break;
    }
    case 4:
    {
        SVREG4(P,T,H,S,V,X);break;
    }
    case 5:
	case 6:
    {
        SVREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  SX_KK (double* P,double*T,double*H, double S, double* V, double X,int * RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        SXREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  SXLP_KK (double* P,double*T,double* H, double  S, double* V, double X, int * RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        SXLPREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  SXMP_KK (double* P,double*T,double*H, double S, double* V, double X, int * RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        SXMPREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  SXHP_KK (double* P,double*T,double*H, double S, double* V, double X,int * RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        SXHPREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  VX_KK (double* P,double*T,double*H,double*S, double V,double X, int* RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        VXREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  VXLP_KK (double* P,double*T,double*H,double*S, double V,double X,int* RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        VXLPREG56(P,T,H,S,V,X);break;
    }
    default:
        *RANGE=0;
    }
}


void  VXHP_KK (double*P,double*T,double*H,double*S, double V,double X, int* RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        VXHPREG56(P,T,H,S,V,X);break;
    }
    default:
       *RANGE=0;
    }
}


void  P2HL_KK(double P,  double* H, int * RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
          *H=P2HLREG56(P); break;
		}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}


void  P2HG_KK(double P,  double* H, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
         *H=P2HGREG56(P);  break;
		}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}


void  P2SL_KK(double P,double*  S, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
          *S=P2SLREG56(P); break;
		}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}


void  P2SG_KK(double P,  double* S, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
         *S=P2SGREG56(P); break;
		}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}



void  P2VL_KK(double P,  double* V, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
         *V=P2VLREG56(P);break;
		}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}


void  P2VG_KK(double P, double* V, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
         *V=P2VGREG56(P);break;
		}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  P2L_KK(double  P, double* T,double*H,double*S,double*V,double*X, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        *T=P2T_KK(P);
        *V=PT2VREG4(P,*T);
        *H=TV2HREG4(*T,*V);
        *S=TV2SREG4(*T,*V);
        *X=0.0;
		break;
    }
    case  6:
    {
        *T=P2T_KK(P);
        *H=PT2HREG1(P,*T);
        *S=PT2SREG1(P,*T);
        *V=PT2VREG1(P,*T);
        *X=0.0;
		break;
    }
    default:
		{
        *RANGE=0;
        *T=-1.0;
        *H=-1.0;
        *S=-1.0;
        *V=-1.0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  P2G_KK(double P, double*T,double*H,double*S,double*V,double*X, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        *T=P2T_KK(P);
        *V=PT2VREG3(P,*T);
        *H=TV2HREG3(*T,*V);
        *S=TV2SREG3(*T,*V);
        *X=1.0;
		break;
    }
    case 6:
    {
        *T=P2T_KK(P);
        *H=PT2HREG2(P,*T);
        *S=PT2SREG2(P,*T);
        *V=PT2VREG2(P,*T);
        *X=1.0;
		break;
    }
    default:
		{
        *RANGE=0;
        *T=-1.0;
        *H=-1.0;
        *S=-1.0;
        *V=-1.0;
        *X=-1.0;
        //exit(1);
		}
    }
}



void  P2CPL(double P,  double* CP, int* RANGE)
{

const double    DT=0.05;

double   T,T1,H,H1;

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        T=P2T_KK(P);
        H=PT2HREG4(P,T);
        T1=T-DT;
        if ( T1<T350C-DeltaVal )
        {
            T1=T+DT;
        }
        H1=PT2HREG4(P,T1);
        if ( *RANGE!=0 )
            *CP=fabs((H-H1)/DT);
        else
            *CP=0.0;
		break;
    }
    case 6:
    {
        T=P2T_KK(P);
        H=PT2HREG1(P,T);
        T1=T-DT;
        if ( T1<T000C-DeltaVal )
        {
            T1=T+DT;
        }
        H1=PT2HREG1(P,T1);
        if ( *RANGE!=0 )
            *CP=fabs((H-H1)/DT);
        else
            *CP=0.0;
		break;
    }
    default: 
    {
        *RANGE=0;
        *CP=-1.0;
    }
    }
}


void  P2CPG(double P,  double* CP, int* RANGE)
{
const double    DT=0.05;

double   T,H,H1,T1;

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        T=P2T_KK(P);
        H=PT2HREG3(P,T);
        T1=T-DT;
        if ( T1<T350C-DeltaVal )
        {
            T1=T+DT;
        }
        H1=PT2HREG3(P,T1);
        if ( *RANGE!=0 )
            *CP=fabs((H-H1)/DT);
        else
            *CP=0.0;
		break;
    }
    case  6:
    {
        T=P2T_KK(P);
        H=PT2HREG2(P,T);
        T1=T-DT;
        if ( T1<T000C-DeltaVal )
        {
            T1=T+DT;
        }
        H1=PT2HREG2(P,T1);
        if ( *RANGE!=0 )
            *CP=fabs((H-H1)/DT);
        else
            *CP=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CP=-1.0;
    }
    }
}


void  P2EL_KK(double P, double* E, int* RANGE)
{
double    T;


    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        T=P2T_KK(P);
        *E=PT2HREG4(P,T)-P*PT2VREG4(P,T);
		break;
    }
    case 6:
    {
        T=P2T_KK(P);
        *E=PT2HREG1(P,T)-P*PT2VREG1(P,T);
		break;
    }
    default:
    {
        *RANGE=0;
        *E=-1.0;
    }
    }
}


void  P2EG_KK(double P,  double* E, int* RANGE)
{
double    T;


    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case  5:
    {
        T=P2T_KK(P);
        *E=PT2HREG3(P,T)-P*PT2VREG3(P,T);
		break;
    }
    case 6:
    {
        T=P2T_KK(P);
        *E=PT2HREG2(P,T)-P*PT2VREG2(P,T);
		break;
    }
    default:
    {
        *RANGE=0;
        *E=-1.0;
    }
    }
}



void  P2CVL(double P, double* CV, int* RANGE)
{
const double    DT=0.025;

double    T,T1,*E,*E1,V,PP,e,e1;

    E=&e; 
	E1=&e1;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        if ( (1-P/PC_Water)<1E-8 )
        {
            *CV=10.0;
            return;
			//exit(1);
        }
        P2EL_KK(P,E,RANGE);
        T=P2T_KK(P);
        V=PT2VREG4(P,T);
        T1=T+DT;
        if ( T1>TC_Water+DeltaVal )
            PP=TV2PREG3(T1,V);
        else
            PP=TV2PREG4(T1,V);
        PT2E_KK(PP,T1,E1,RANGE);
        if ( *RANGE!=0 )
            *CV=fabs((*E1-*E)/DT);
        else
            *CV=0.0;
		break;
    }
    case 6:
    {
        P2EL_KK(P,E,RANGE);
        T=P2T_KK(P);
        V=PT2VREG1(P,T);
        T1=T+DT;
        if ( T1>T350C+DeltaVal )
            T1=T-DT;
        PP=TV2PREG1(T1,V);
        PT2E_KK(PP,T1,E1,RANGE);
        if ( *RANGE!=0 )
            *CV=fabs((*E1-*E)/DT);
        else
            *CV=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CV=-1.0;
    }
    }
}


void  P2CVG(double P, double* CV, int* RANGE)
{
const double    DT=0.025;

double    T,T1,*E,*E1,V,PP,e,e1;

    E=&e;
	E1=&e1;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        if ( (1-P/PC_Water)<1E-8 )
        {
            *CV=10.0;
            return;
			//exit(1);
        }
        P2EG_KK(P,E,RANGE);
        T=P2T_KK(P);
        V=PT2VREG3(P,T);
        T1=T+DT;
        PP=TV2PREG3(T1,V);
        if ((T1<=TC_Water) && (fabs(1-PP/T2P_KK(T1))<=1E-5) )
            T2EG(T1,E1,RANGE);
        else
            PT2E_KK(PP,T1,E1,RANGE);
        
		if ( *RANGE!=0 )
            *CV=fabs((*E1-*E)/DT);
        else
            *CV=0.0;
		break;
    }
    case 6:
    {
        P2EG_KK(P,E,RANGE);
        T=P2T_KK(P);
        V=PT2VREG2(P,T);
        T1=T+DT;
        if ( T1>T350C+DeltaVal )
            T1=T-DT;
        PP=TV2PREG2(T1,V);
        if ((fabs(1-PP/T2P_KK(T1))<=1E-5) )
            T2EG(T1,E1,RANGE);
        else
            PT2E_KK(PP,T1,E1,RANGE);
        
		if ( *RANGE!=0 )
            *CV=fabs((*E1-*E)/DT);
        else
            *CV=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CV=-1.0;
    }
    }
}


void  P2SSPL(double P,  double* SSP, int* RANGE)
{
const double    RDP=5E-5;

double   T,*H,P1,*V1,*S,*V,*X,h,s,v,x,v1;

    H=&h;
	S=&s;
	V=&v;
	X=&x;
	V1=&v1;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        T=P2T_KK(P);
        PTREG4(P,T,H,S,V,X);
        P1=P*(1+RDP);
        if ( (P1>PC_Water+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            *SSP=0.0;
        }
		break;
    }
    case 6:
    {
        T=P2T_KK(P);
        PTREG1(P,T,H,S,V,X);
        P1=P*(1+RDP);
        if ( (P1>T2P_KK(T350C)+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            *SSP=0.0;
        }
		break;
    }
    default:
    {
        *RANGE=0;
        *SSP=-1.0;
    }
    }
}


void  P2SSPG(double P,  double* SSP, int* RANGE)
{
const double    RDP=5E-5;

double  T,*H,P1,*V1,*S,*V,*X,h,s,v,x,v1;

    H=&h;
	S=&s;
	V=&v;
	X=&x;
	V1=&v1;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        T=P2T_KK(P);
        PTREG3(P,T,H,S,V,X);
        P1=P*(1+RDP);
        if ( (P1>PC_Water+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            *SSP=0.0;
        }
		break;
    }
    case 6:
    {
        T=P2T_KK(P);
        PTREG2(P,T,H,S,V,X);
        P1=P*(1+RDP);
        if ( (P1>T2P_KK(T350C)+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            *SSP=0.0;
        }
		break;
    }
    default:
    {
        *RANGE=0;
        *SSP=-1.0;
    }
    }
}


void  P2KSL(double P,  double* KS, int* RANGE)
{
const double    RDP=5E-5;

double    T,S,V,*V1,P1,v1;

    V1=&v1;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        T=P2T_KK(P);
        S=PT2SREG4(P,T);
        V=PT2VREG4(P,T);
        P1=P*(1+RDP);
        if ( P1>PC_Water )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,S,V1,RANGE);
        if ( *RANGE!=0 )
            *KS=fabs(-V/P*(P1-P)/((*V1)-V));
        else
            *KS=0.0;
		break;
    }
    case 6:
    {
        T=P2T_KK(P);
        S=PT2SREG1(P,T);
        V=PT2VREG1(P,T);
        P1=P*(1+RDP);
        if ( P1>T2P_KK(T350C) )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,S,V1,RANGE);
        if ( *RANGE!=0 )
            *KS=fabs(-V/P*(P1-P)/((*V1)-V));
        else
            *KS=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *KS=-1.0;
    }
    }
}


void  P2KSG(double P,  double* KS, int* RANGE)
{
const double    RDP=5E-5;

double    T,S,V,*V1,P1,v1;

    V1=&v1;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        T=P2T_KK(P);
        S=PT2SREG3(P,T);
        V=PT2VREG3(P,T);
        P1=P*(1+RDP);
        if ( P1>PC_Water )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,S,V1,RANGE);
        if ( *RANGE!=0 )
            *KS=fabs(-V/P*(P1-P)/((*V1)-V));
        else
            *KS=0.0;
		break;
    }
    case 6:
    {
        T=P2T_KK(P);
        S=PT2SREG2(P,T);
        V=PT2VREG2(P,T);
        P1=P*(1+RDP);
        if ( P1>T2P_KK(T350C) )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,S,V1,RANGE);
        if ( *RANGE!=0 )
            *KS=fabs(-V/P*(P1-P)/((*V1)-V));
        else
            *KS=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *KS=-1.0;
    }
    }
}

double NSATLIQP(double  P,double Lamd)
{
double    T;

    T=P2T_KK(P);
    return TVLAMD2N_KK(T,T2VLREG56(T),Lamd);
}

void  P2NL(double P,double Lamd,  double* N, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	 case 5:
     case 6:
    {
        *N=NSATLIQP(P,Lamd); break;
    }
	default:
    {
        *RANGE=0;
        *N=-1.0;
    }
    }
}

double NSATVAPP(double  P,double Lamd)
{
double    T;

    T=P2T_KK(P);
    return TVLAMD2N_KK(T,T2VGREG56(T),Lamd);
}

void  P2NG(double  P,double Lamd,  double* N, int* RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        *N=NSATVAPP(P,Lamd);break;
    }
    default:
    {
        *RANGE=0;
        *N=-1.0;
    }
    }
}



void  PT2H_KK(double P,double T,  double* H, int* RANGE)
{

    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        *H=PT2HREG1(P,T);break;
    }
    case 2:
    {
        *H=PT2HREG2(P,T);break;
    }
    case 3:
    {
        *H=PT2HREG3(P,T);break;
    }
    case 4:
    {
       *H=PT2HREG4(P,T);break;
    }
    case 5:
	case 6:
    {
        *H=T2HLREG56(T);break;
    }
    default:
		{
        *H=-1.0;
        *RANGE=0;
		}
    }
}

void  PT2S_KK(double P,double T,  double* S, int* RANGE)
{

    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        *S=PT2SREG1(P,T);break;
    }
    case 2:
    {
        *S=PT2SREG2(P,T);break;
    }
    case 3:
    {
        *S=PT2SREG3(P,T);break;
    }
    case 4:
    {
        *S=PT2SREG4(P,T);break;
    }
    case 5:
	case 6:
    {
        *S=T2SLREG56(T);break;
    }
    default:
		{
        *S=-1.0;
        *RANGE=0;
		}
    }
}


void  PT2V_KK(double P,double T,  double* V, int* RANGE)
{

    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        *V=PT2VREG1(P,T);break;
    }
    case 2:
    {
        *V=PT2VREG2(P,T);break;
    }
    case 3:
    {
        *V=PT2VREG3NT(P,T);break;
    }
    case 4:
    {
        *V=PT2VREG4NT(P,T);break;
    }
    case 5:
	case 6:
    {
        *V=T2VLREG56(T);break;
    }
    default:
		{
        *V=-1.0;
        *RANGE=0;
		}
    }
}

void  PT2X(double P,double T,  double* X, int* RANGE)
{

    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
    {
        *X=0.0;break;
    }
    case 2:
    {
        *X=1.0;break;
    }
    case 3:
    {
        *X=1.0;break;
    }
    case 4:
    {
        *X=0.0;break;
    }
    case 5:
	case 6:
    {
        *X=0.0;break;
    }
    default:
		{
        *X=-1.0;
        *RANGE=0;
		}
    }
}


void  PT2CP_KK(double  P,double T, double * CP, int* RANGE)
{
const double    DT=0.05;

double    *H,*H1,T1,h1,h;
 
    H1=&h1;
	H=&h;
    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
    {
        if ( (P<PC_Water) && (T<TC_Water) )
        {
            if ( fabs(1-T2P_KK(T)/P)<1E-8 )
            {
                P2CPL(P,CP,RANGE);
                return;
				//exit(1);
            }
        }

        PT2H_KK(P,T,H,RANGE);
        T1=T+DT;
        if ( T1>800+DeltaVal )
            T1=T-DT;
        PT2H_KK(P,T1,H1,RANGE);
        if ( *RANGE!=0 )
            *CP=fabs((*H-*H1)/DT);
        else
            *CP=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CP=-1.0;
    }
    }
}

void  PT2E_KK(double P,double T, double* E, int* RANGE)
{

    SUBRANGEBYPT(P,T,RANGE);
    if ( *RANGE>0 )
    {
        if ( (P<PC_Water) && (T<TC_Water) )
        {
            if ( fabs(1-T2P_KK(T)/P)<1E-5 )
            {
                P2EL_KK(P,E,RANGE);
                return;
				//exit(1);
            }
        }
    }

    switch ( *RANGE )
    {
	case 1:
    {
        *E=PT2HREG1(P,T)-P*PT2VREG1(P,T);break;
    }
    case 2:
    {
        *E=PT2HREG2(P,T)-P*PT2VREG2(P,T);break;
    }
    case 3:
    {
        *E=PT2HREG3(P,T)-P*PT2VREG3(P,T);break;
    }
    case 4:
    {
        *E=PT2HREG4(P,T)-P*PT2VREG4(P,T);break;
    }
    case 5:
	case 6:
    {
        P2EL_KK(P,E,RANGE);break;
    }
    default:
    {
        *RANGE=0;
        *E=-1.0;
    }
    }
}


void  PT2CV_KK(double  P,double T, double* CV, int* RANGE)
{
const double    DT=0.05;

double   *E,*E1,T1,*V,*PP,e,v,e1,pp;

    E=&e;
	V=&v;
	E1=&e1;
	PP=&pp;
    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
    {
        if ( (P<PC_Water) && (T<TC_Water) )
        {
            if ( fabs(1-T2P_KK(T)/P)<1E-5 )
            {
                P2CVL(P,CV,RANGE);
                return;
				//exit(1);
            }
        }

        PT2E_KK(P,T,E,RANGE);
        PT2V_KK(P,T,V,RANGE);
        T1=T+DT;
        if ( T1>800+DeltaVal )
            T1=T-DT;
        TV2P_KK(T1,*V,PP,RANGE);
        PT2E_KK(*PP,T1,E1,RANGE);
        if ( *RANGE!=0 )
            *CV=fabs((*E-*E1)/DT);
        else
            *CV=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CV=-1.0;
    }
    }
}


void  PT2SSP_KK(double P,double T, double* SSP, int* RANGE)
{
const double    RDP=5E-5;

double    *H,*S,*V,P1,*V1,*X,*S_MIN,h,s,v,x,v1,s_min;

    H=&h;
	S=&s;
	V=&v;
	X=&x;
	V1=&v1;
	S_MIN=&s_min;
    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
    {
        if ( (P<PC_Water) && (T<TC_Water) )
        {
            if ( fabs(1-T2P_KK(T)/P)<1E-8 )
            {
                P2SSPL(P,SSP,RANGE);
                return;
				//exit(1);
            }
        }

        PT_KK(P,T,H,S,V,X,RANGE);
        P1=P*(1+RDP);
        if ( (P1>100+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            if ( (P>100-DeltaVal) && (T<0.0+DeltaVal) )
            {
                PT2S_KK(P1,T,S_MIN,RANGE);
                PS2V_KK(P,*S_MIN,V,RANGE);
                PS2V_KK(P1,*S_MIN,V1,RANGE);
                if ( *RANGE!=0 )
                {
                    *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
                    return;
					//exit(1);
                }
            }
            *SSP=0.0;
        }
		break;
    }
    default:
    {
        *RANGE=0;
        *SSP=-1.0;
    }
    }
}


void  PT2KS_KK(double  P,double T, double* KS, int* RANGE)
{
const double    RDP=5E-5;

double    *S,*V,P1,*V1,s,v,v1;

    S=&s;
	V=&v;
	V1=&v1;
    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4:
    case 5:
	case 6:
    {
        if ( (P<PC_Water) && (T<TC_Water) )
        {
            if ( fabs(1-T2P_KK(T)/P)<1E-8 )
            {
                P2KSL(P,KS,RANGE);
                return;
				//exit(1);
            }
        }

        PT2S_KK(P,T,S,RANGE);
        PT2V_KK(P,T,V,RANGE);
        P1=P*(1+RDP);
        if ( P1>100 )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *KS=fabs(-(*V)/P*(P1-P)/((*V1)-(*V)));
        else
            *KS=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *KS=-1.0;
    }
    }
}

double TV2EPS_KK(double  T,double V)
{
const double    CA= 7.62571;
const double    CB= 244.003;
const double    CC=-140.569;
const double    CD= 27.7841;
const double    CE=-96.2805;
const double    CF= 41.7909;
const double    CG=-10.2099;
const double    CH=-45.2059;
const double    CI= 84.6395;
const double    CK=-35.8644;

double    TStar,DensStar,EPSr;

    TStar=(T+273.15)/289.15;
    DensStar=(1/V)*1E-3;

    EPSr=((CH/TStar+CI)/TStar+CK);
    EPSr=EPSr*DensStar+((CG*TStar+CF)*TStar+CE/TStar);
    EPSr=EPSr*DensStar+(CD*TStar+CC+CB/TStar);
    EPSr=EPSr*DensStar+CA/TStar;
    EPSr=EPSr*DensStar+1;

    return EPSr;
}

void  PT2EPS_KK(double P,double T, double* EPS, int* RANGE)
{
double    *V,v;

    V=&v;
    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
    case 3:
	case 4:
	case 5:
	case 6:
    {
        if ( (P<PC_Water) && (T<TC_Water) )
        {
            if ( fabs(1-T2P_KK(T)/P)<1E-8 )
            {
                T2EPSL(T,EPS,RANGE);
                return;
				//exit(1);
            }
        }

        PT2V_KK(P,T,V,RANGE);
        *EPS=TV2EPS_KK(T,*V);
		break;
    }
    default:
    {
        *RANGE=0;
        *EPS=-1.0;
    }
    }
}

double TVLAMD2N_KK(double  T,double V, double Lamd)
{
const double    A1 = 3.036167E-3;
const double    A2 = 0.052421;
const double    A3 = 2.117579E-1;
const double    A4 =-5.195756E-2;
const double    A5 = 5.922248E-2;
const double    A6 =-1.918429E-2;
const double    A7 = 2.58235E-3;
const double    A8 =-2.352054E-4;
const double    A9 = 3.964628E-5;
const double    A10= 3.336153E-2;
const double    A11=-4.008264E-2;
const double    A12= 8.339681E-3;
const double    A13=-1.054741E-2;
const double    A14= 9.491575E-3;

double    LamdStar,TStar,DensStar,NN;

    LamdStar=Lamd*(1E-9)/589.0;
    TStar=(T+273.15)/273.15;
    DensStar=1E-3*(1/V);

    NN=A1/(LamdStar*LamdStar-A2)+A3+((((A8*LamdStar+A7)*LamdStar+A6)*LamdStar+A5)*LamdStar+A4)*LamdStar*LamdStar;
    NN=NN+A9/DensStar+((A12*LamdStar+A11)*LamdStar+A10)*TStar+(A13+A14*LamdStar)*LamdStar*TStar*TStar;
    NN=NN*DensStar;

    return pow(3.0/(1.0-NN)-2.0,0.5);
}

void  PT2N_KK(double  P,double T,double Lamd, double* N, int * RANGE) 
{
double    *V,v;

    V=&v;
    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
    {
        if ( (P<PC_Water) && (T<TC_Water) )
        {
            if ( fabs(1-T2P_KK(T)/P)<1E-8 )
            {
                T2NL(T,Lamd,N,RANGE);
                return;
				//exit(1);
            }
        }

        PT2V_KK(P,T,V,RANGE);
        *N=TVLAMD2N_KK(T,*V,Lamd);
		break;
    }
    default:
    {
        *RANGE=0;
        *N=-1.0;
    }
    }

}



void  T2HL(double T,  double* H, int* RANGE)
{
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
        { 
			*H=T2HLREG56(T); break; 
		}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}


void  T2HG(double T, double* H, int *RANGE)
{
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
        {
			*H=T2HGREG56(T);break;
		}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}


void  T2SL(double T,  double* S, int *RANGE)
{
    

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
        {
			*S=T2SLREG56(T);break;
		}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}


void  T2SG(double T,  double* S, int *RANGE)
{


    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
        *S=T2SGREG56(T);break;
		}
    default:
		{
		*RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}


void  T2VL(double T, double* V, int *RANGE)
{
    

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
        *V=T2VLREG56(T);break;
		}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}





void  T2VG(double T, double* V, int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{
        *V=T2VGREG56(T); break;
		}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  T2L(double* P, double T, double* H,double*S,double*V,double*X, int * RANGE)
{
    SUBRANGEBYP(*P,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        *P=T2P_KK(T);
        *V=PT2VREG4(*P,T);
        *H=TV2HREG4(T,*V);
        *S=TV2SREG4(T,*V);
        *X=0.0;
		break;
    }
    case 6:
    {
        *P=T2P_KK(T);
        *H=PT2HREG1(*P,T);
        *S=PT2SREG1(*P,T);
        *V=PT2VREG1(*P,T);
        *X=0.0;
		break;
    }
    default:
		{
        *RANGE=0;
        *P=-1.0;
        *H=-1.0;
        *S=-1.0;
        *V=-1.0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  T2G(double* P, double T, double* H,double*S,double*V,double*X,int *RANGE)
{
    SUBRANGEBYP(*P,RANGE);
    switch ( *RANGE )
    {
	case  5:
    {
        *P=T2P_KK(T);
        *V=PT2VREG3(*P,T);
        *H=TV2HREG3(T,*V);
        *S=TV2SREG3(T,*V);
        *X=1.0;
		break;
    }
    case  6:
    {
        *P=T2P_KK(T);
        *H=PT2HREG2(*P,T);
        *S=PT2SREG2(*P,T);
        *V=PT2VREG2(*P,T);
        *X=1.0;
		break;
    }
    default:
		{
        *RANGE=0;
        *P=-1.0;
        *H=-1.0;
        *S=-1.0;
        *V=-1.0;
        *X=-1.0;
        //exit(1);
		}
    }
}



void  T2CPL(double T, double* CP,int *RANGE)
{
const double    DT=0.05;

double   P,H,H1,T1;
    
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        P=T2P_KK(T);
        H=PT2HREG4(P,T);
        T1=T-DT;
        if ( T1<T350C-DeltaVal )
        {
            T1=T+DT;
        }
        H1=PT2HREG4(P,T1);
        if ( *RANGE!=0 )
            *CP=fabs((H-H1)/DT);
        else
            *CP=0.0;
		break;
    }
    case  6:
    {
        P=T2P_KK(T);
        H=PT2HREG1(P,T);
        T1=T-DT;
        if ( T1<T000C-DeltaVal )
        {
            T1=T+DT;
        }
        H1=PT2HREG1(P,T1);
        if ( *RANGE!=0 )
            *CP=fabs((H-H1)/DT);
        else
            *CP=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CP=-1.0;
    }
    }
}

void  T2CPG(double T, double* CP,int *RANGE)
{
const double    DT=0.05;

double   P,H,H1,T1;

    
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        P=T2P_KK(T);
        H=PT2HREG3(P,T);
        T1=T-DT;
        if ( T1<T350C-DeltaVal )
        {
            T1=T+DT;
        }
        H1=PT2HREG3(P,T1);
        if ( *RANGE!=0 )
            *CP=fabs((H-H1)/DT);
        else
            *CP=0.0;
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        H=PT2HREG2(P,T);
        T1=T-DT;
        if ( T1<T000C-DeltaVal )
        {
            T1=T+DT;
        }
        H1=PT2HREG2(P,T1);
        if ( *RANGE!=0 )
            *CP=fabs((H-H1)/DT);
        else
            *CP=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CP=-1.0;
    }
    }
}


void  T2EL(double T, double* E,int *RANGE)
{
double   P;

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        P=T2P_KK(T);
        *E=PT2HREG4(P,T)-P*PT2VREG4(P,T);
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        *E=PT2HREG1(P,T)-P*PT2SREG1(P,T);
		break;
    }
    default:
    {
        *RANGE=0;
        *E=-1.0;
    }
    }
}

void  T2EG(double T, double* E,int *RANGE)
{
double    P;
    
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        P=T2P_KK(T);
        *E=PT2HREG3(P,T)-P*PT2VREG3(P,T);
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        *E=PT2HREG2(P,T)-P*PT2VREG2(P,T);
		break;
    }
    default:
    {
        *RANGE=0;
        *E=-1.0;
    }
    }
}


void  T2CVL(double T, double* CV,int *RANGE)
{
const double    DT=0.025;

double    P,T1,*E,*E1,V,PP,e,e1;

    E=&e;
	E1=&e1;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        if ( (1-T/TC_Water)<1E-8 )
        {
            *CV=10.0;
            return;
			//exit(1);
        }
        P=T2P_KK(T);
        P2EL_KK(P,E,RANGE);
        V=PT2VREG4(P,T);
        T1=T+DT;
        if ( T1>TC_Water+DeltaVal )
            PP=TV2PREG3(T1,V);
        else
            PP=TV2PREG4(T1,V);
        PT2E_KK(PP,T1,E1,RANGE);
        if ( *RANGE!=0 )
            *CV=fabs((*E1-*E)/DT);
        else
            *CV=0.0;
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        P2EL_KK(P,E,RANGE);
        V=PT2VREG1(P,T);
        T1=T+DT;
        if ( T1>T350C+DeltaVal )
            T1=T-DT;
        PP=TV2PREG1(T1,V);
        PT2E_KK(PP,T1,E1,RANGE);
        if ( *RANGE!=0 )
            *CV=fabs((*E1-*E)/DT);
        else
            *CV=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CV=-1.0;
    }
    }
}


void  T2CVG(double T, double* CV,int *RANGE)
{
 double   DT=0.025;

double    P,T1,*E,*E1,V,PP,e,e1;

    E=&e;
	E1=&e1;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        if ( (1-T/TC_Water)<1E-8 )
        {
            *CV=10.0;
            return;
			//exit(1);
        }
        P=T2P_KK(T);
        P2EG_KK(P,E,RANGE);
        V=PT2VREG3(P,T);
        T1=T+DT;
        PP=TV2PREG3(T1,V);
        if ((T1<=TC_Water) && (fabs(1-PP/T2P_KK(T1))<=1E-5) )
            T2EG(T1,E1,RANGE);
        else
            PT2E_KK(PP,T1,E1,RANGE);
        if ( *RANGE!=0 )
            *CV=fabs((*E1-*E)/DT);
        else
            *CV=0.0;
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        P2EG_KK(P,E,RANGE);
        V=PT2VREG2(P,T);
        T1=T+DT;
        if ( T1>T350C+DeltaVal )
            T1=T-DT;
        PP=TV2PREG2(T1,V);
        if (fabs(1-PP/T2P_KK(T1))<=1E-5) 
            T2EG(T1,E1,RANGE);
        else
            PT2E_KK(PP,T1,E1,RANGE);
        if ( *RANGE!=0 )
            *CV=fabs((*E1-*E)/DT);
        else
            *CV=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *CV=-1.0;
    }
    }
}

void  T2SSPL(double T,double* SSP,int *RANGE)
{
const double   RDP=5E-5;

double   P,*H,P1,*V1,*S,*V,*X,h,s,v,x,v1;
  
    H=&h;
	S=&s;
	V=&v;
	X=&x;
	V1=&v1;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        P=T2P_KK(T);
        PTREG4(P,T,H,S,V,X);
        P1=P*(1+RDP);
        if ( (P1>PC_Water+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            *SSP=0.0;
        }
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        PTREG1(P,T,H,S,V,X);
        P1=P*(1+RDP);
        if ( (P1>T2P_KK(T350C)+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            *SSP=0.0;
        }
		break;
    }
    default:
    {
        *RANGE=0;
        *SSP=-1.0;
    }
    }
}


void  T2SSPG(double T, double* SSP,int *RANGE)
{
const double    RDP=5E-5;

double   P,*H,P1,*V1,*S,*V,*X,h,s,v,x,v1;
  
    H=&h;
	S=&s;
	V=&v;
	X=&x;
	V1=&v1;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        P=T2P_KK(T);
        PTREG3(P,T,H,S,V,X);
        P1=P*(1+RDP);
        if ( (P1>PC_Water+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            *SSP=0.0;
        }
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        PTREG2(P,T,H,S,V,X);
        P1=P*(1+RDP);
        if ( (P1>T2P_KK(T350C)+DeltaVal) )
            P1=P*(1-RDP);
        PS2V_KK(P1,*S,V1,RANGE);
        if ( *RANGE!=0 )
            *SSP=sqrt(fabs(-(*V)*(*V)*(P1-P)*1E6/((*V1)-(*V))));
        else
        {
            *SSP=0.0;
        }
		break;
    }
    default:
    {
        *RANGE=0;
        *SSP=-1.0;
    }
    }
}



void  T2KSL(double T, double* KS,int *RANGE)
{
const double    RDP=5E-5;

double    P,S,V,*V1,P1,v1;    

    V1=&v1;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        P=T2P_KK(T);
        S=PT2SREG4(P,T);
        V=PT2VREG4(P,T);
        P1=P*(1+RDP);
        if ( P1>PC_Water )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,S,V1,RANGE);
        if ( *RANGE!=0 )
            *KS=fabs(-V/P*(P1-P)/((*V1)-V));
        else
            *KS=0.0;
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        S=PT2SREG1(P,T);
        V=PT2VREG1(P,T);
        P1=P*(1+RDP);
        if ( P1>T2P_KK(T350C) )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,S,V1,RANGE);  //Tiao Shi
        if ( *RANGE!=0 )
            *KS=fabs(-V/P*(P1-P)/((*V1)-V));
        else
            *KS=0.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *KS=-1.0;
    }
    }
}


void  T2KSG(double T,  double* KS, int *RANGE)
{
const double    RDP=5E-5;

double    P,S,V,*V1,P1,v1;

    V1=&v1;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
    {
        P=T2P_KK(T);
        S=PT2SREG3(P,T);
        V=PT2VREG3(P,T);
        P1=P*(1+RDP);
        if ( P1>PC_Water )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,S,V1,RANGE);
        if ( *RANGE!=0 )
            *KS=fabs(-V/P*(P1-P)/((*V1)-V));
        else
            *KS=0.0;
		break;
    }
    case 6:
    {
        P=T2P_KK(T);
        S=PT2SREG2(P,T);
        V=PT2VREG2(P,T);
        P1=P*(1+RDP);
        if ( P1>T2P_KK(T350C) )
        {
            P1=P*(1-RDP);
        }
        PS2V_KK(P1,S,V1,RANGE);
        if ( *RANGE!=0 )
		{
            *KS=fabs(-V/P*(P1-P)/((*V1)-V));
		}
        else
		{
            *KS=0.0;
		}
		break;
    }
    default:
    {
        *RANGE=0;
        *KS=-1.0;
    }
    }
}


double EPSSATLIQT(double  T)
{
    return TV2EPS_KK(T,T2VLREG56(T));
}

void  T2EPSL(double T,  double* EPS, int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        *EPS=EPSSATLIQT(T);break;
    }
    default:
    {
        *RANGE=0;
        *EPS=-1.0;
    }
    }
}

double EPSSATVAPT(double  T)
{
    return TV2EPS_KK(T,T2VGREG56(T));
}

void  T2EPSG(double T, double* EPS,int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        *EPS=EPSSATVAPT(T);break;
    }
    default:
    {
        *RANGE=0;
        *EPS=-1.0;
    }
    }
}

double NSATLIQT(double  T,double Lamd)
{
    return TVLAMD2N_KK(T,T2VLREG56(T),Lamd);
}

void  T2NL(double T,double Lamd,  double* N, int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        *N=NSATLIQT(T,Lamd);break;
    }
    default:
    {
        *RANGE=0;
        *N=-1.0;
    }
    }
}

double NSATVAPT(double  T,double Lamd)
{
    return TVLAMD2N_KK(T,T2VGREG56(T),Lamd);
}

void  T2NG(double T,double Lamd,  double* N, int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        *N=NSATVAPT(T,Lamd);break;
    }
    default:
    {
        *RANGE=0;
        *N=-1.0;
    }
    }
}



void  T2SurfT(double T, double* SURFT, int *RANGE)
{
double    Sta;    

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
    {
        Sta=(T+273.15)/(TC_Water+273.15);
        *SURFT=235.8*pow(1-Sta,1.256)*(1-0.625*(1-Sta))/1000.0;
		break;
    }
    default:
    {
        *RANGE=0;
        *SURFT=-1.0;
        //exit(1);
    }
    }
}


void  PH2T_KK(double P,double H,  double* T, int *RANGE)
{
    

    SUBRANGEBYPH(P,H,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*T=PH2TREG1(P,H); break;}
    case 2:
		{ *T=PH2TREG2(P,H);break;}
    case 3:
        {*T=PH2TREG3(P,H);break;}
    case 4:
        {*T=PH2TREG4(P,H);break;}
    case 5:
	case 6:
        {*T=P2T_KK(P);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  PH2S_KK(double P,double H,  double* S, int *RANGE)
{
    

    SUBRANGEBYPH(P,H,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*S=PH2SREG1(P,H);break;}
    case 2:
        {*S=PH2SREG2(P,H);break;}
    case 3:
        {*S=PH2SREG3(P,H);break;}
    case 4:
        {*S=PH2SREG4(P,H);break;}
    case 5:
	case 6:
		{ *S=PH2SREG56(P,H);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}


void  PH2V_KK(double P,double H, double* V, int *RANGE)
{
    

    SUBRANGEBYPH(P,H,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*V=PH2VREG1(P,H);break;}
    case 2:
        {*V=PH2VREG2(P,H);break;}
    case 3:
        {*V=PH2VREG3(P,H);break;}
    case 4:
        {*V=PH2VREG4(P,H);break;}
    case 5:
	case 6:
        {*V=PH2VREG56(P,H);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  PH2X_KK(double P,double H, double* X,int *RANGE)
{
    SUBRANGEBYPH(P,H,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*X=0.0;break;}
    case 2:
        {*X=1.0;break;}
    case 3:
        {*X=1.0;break;}
    case 4:
        {*X=0.0;break;}
    case 5:
	case 6:
        {*X=PH2XREG56(P,H);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}
void  PH2SSP_KK(double P,double H, double* SSP,int *RANGE)
{
 double *X,x,*V,v,k,*T,t1;
    X=&x;
	T=&t1;
	V =&v;
    SUBRANGEBYPH(P,H,RANGE);
    switch ( *RANGE )
    {
    case 5:         
	case 6:         
        {	
			PH2X_KK(P,H,X,RANGE);
			PH2V_KK(P,H,V,RANGE);
        	k = 1.035+0.1*( *X );
	        *SSP =  pow((k*(*V)*P*(10E5)),0.5);		
			break;}   
    case 1:
	case 2:
	case 3:
	case 4:
		{ 	
			PH2T_KK(P,H,T,RANGE);
			PT2SSP_KK(P,*T,SSP,RANGE);
			break;}
	default:{
		*RANGE=0;
		*SSP=-1.0;
			break ;	}
    }
}
void  PS2T_KK(double P,double S,  double* T,int *RANGE)
{
    

    SUBRANGEBYPS(P,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *T=PS2TREG1(P,S);break;}
    case 2:
		{ *T=PS2TREG2(P,S);break;}
    case 3:
		{ *T=PS2TREG3(P,S);break;}
    case 4:
		{ *T=PS2TREG4(P,S);break;}
    case 5:
	case 6:
        {*T=P2T_KK(P);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  PS2H_KK(double P,double S, double* H, int *RANGE)
{

    SUBRANGEBYPS(P,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*H=PS2HREG1(P,S);break;}
    case  2:
        {*H=PS2HREG2(P,S);break;}
    case 3:
		{ *H=PS2HREG3(P,S);break;}
    case 4:
        {*H=PS2HREG4(P,S);break;}
    case 5:
	case 6:
		{ *H=PS2HREG56(P,S);break;}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}

void  PS2V_KK(double P,double S,  double* V, int *RANGE)
{

    SUBRANGEBYPS(P,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*V=PS2VREG1(P,S);break;}
    case 2:
        {*V=PS2VREG2(P,S);break;}
    case 3:
        {*V=PS2VREG3(P,S);break;}
    case 4:
        {*V=PS2VREG4(P,S);break;}
    case 5:
	case 6:
        {*V=PS2VREG56(P,S);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}


void  PS2X_KK(double P,double S,  double* X, int *RANGE)
{

    SUBRANGEBYPS(P,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*X=0.0;break;}
    case 2:
        {*X=1.0;break;}
    case 3:
        {*X=1.0;break;}
    case 4:
        {*X=0.0;break;}
    case 5:
	case 6:
		{*X=PS2XREG56(P,S);break;}
	default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  PV2T_KK(double P,double V, double* T, int * RANGE)
{
    

    SUBRANGEBYPV(P,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*T=PV2TREG1(P,V);break;}
    case 2:
        {*T=PV2TREG2(P,V);break;}
    case 3:
        {*T=PV2TREG3(P,V);break;}
    case 4:
        {*T=PV2TREG4(P,V);break;}
    case 5:
	case 6:
		{ *T=P2T_KK(P);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  PV2H_KK(double P,double V,  double* H, int *RANGE)
{

    SUBRANGEBYPV(P,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*H=PV2HREG1(P,V);break;}
    case 2:
        {*H=PV2HREG2(P,V);break;}
    case 3:
        {*H=PV2HREG3(P,V);break;}
    case 4:
        {*H=PV2HREG4(P,V);break;}
    case 5:
	case 6:
		{ *H=PV2HREG56(P,V);break;}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}

void  PV2S_KK(double P,double V,  double* S, int* RANGE)
{

    SUBRANGEBYPV(P,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*S=PV2SREG1(P,V);break;}
    case 2:
        {*S=PV2SREG2(P,V);break;}
    case 3:
        {*S=PV2SREG3(P,V);break;}
    case 4:
        {*S=PV2SREG4(P,V);break;}
    case 5:
	case 6:
        {*S=PV2SREG56(P,V);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}


void  PV2X_KK(double P,double V,  double* X, int *RANGE)
{

    SUBRANGEBYPV(P,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*X=0.0;break;}
    case 2:
        {*X=1.0;break;}
    case 3:
        {*X=1.0;break;}
    case 4:
        {*X=0.0;break;}
    case 5:
	case 6:
        {*X=PV2XREG56(P,V);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  PX2T_KK(double P,double X, double* T, int *RANGE)
{
    SUBRANGEBYPX(P,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *T=P2T_KK(P);break;}
	default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  PX2H_KK(double P,double X, double* H, int *RANGE)
{
    SUBRANGEBYPX(P,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *H=PX2HREG56(P,X);break;}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}

void  PX2S_KK(double P,double X, double* S,int *RANGE)
{
    SUBRANGEBYPX(P,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
        {*S=PX2SREG56(P,X);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}


void  PX2V_KK(double  P,double X, double* V,int *RANGE)
{
    SUBRANGEBYPX(P,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *V=PX2VREG56(P,X);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  TH2P_KK(double T,double H,  double* P,int *RANGE)
{
bool    *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
		{ *P=TH2PREG1(T,H);break;}
    case 2:
		{ *P=TH2PREG2(T,H);break;}
    case 3:
		{ *P=TH2PREG3(T,H);break;}
    case 4:
        {*P=TH2PREG4(T,H);break;}
    case 5:
	case 6:
		{ *P=T2P_KK(T);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  TH2PLP_KK(double T,double H,  double* P,int *RANGE)
{
bool    *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
        {*P=TH2PLPREG1(T,H);break;}
    case 2:
        {*P=TH2PREG2(T,H);break;}
    case 3:
        {*P=TH2PREG3(T,H);break;}
    case 4:
        {*P=TH2PREG4(T,H);break;}
    case 5:
	case 6:
		{*P=T2P_KK(T);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}


void  TH2PHP_KK(double T,double H, double* P,int *RANGE)
{
bool    *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=1;
    switch ( *RANGE )
    {
	case 1:
        {*P=TH2PHPREG1(T,H);break;}
    case 2:
        {*P=TH2PREG2(T,H);break;}
    case 3:
        {*P=TH2PREG3(T,H);break;}
    case 4:
        {*P=TH2PREG4(T,H);break;}
    case 5:
	case 6:
		{ *P=T2P_KK(T);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}


void  TH2S_KK(double T,double H,  double* S, int *RANGE)
{
 bool   *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
        {*S=TH2SREG1(T,H);break;}
    case 2:
        {*S=TH2SREG2(T,H);break;}
    case 3:
        {*S=TH2SREG3(T,H);break;}
    case 4:
        {*S=TH2SREG4(T,H);break;}
    case 5:
	case 6:
        {*S=TH2SREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}

void  TH2SLP_KK(double T,double H, double* S,int *RANGE)
{
 bool   *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
        {*S=TH2SLPREG1(T,H);break;}
    case 2:
		{ *S=TH2SREG2(T,H);break;}
    case 3:
        {*S=TH2SREG3(T,H);break;}
    case 4:
        {*S=TH2SREG4(T,H);break;}
    case 5:
	case 6:
		{*S=TH2SREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}


void  TH2SHP_KK(double T,double H, double* S, int *RANGE)
{
 bool   *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=1;
    switch ( *RANGE )
	{
	case 1:
        {*S=TH2SHPREG1(T,H);break;}
    case 2:
        {*S=TH2SREG2(T,H);break;}
    case 3:
		{ *S=TH2SREG3(T,H);break;}
    case 4:
		{ *S=TH2SREG4(T,H);break;}
    case 5:
	case 6:
        {*S=TH2SREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}

void  TH2V_KK(double T,double H, double* V,int *RANGE)
{
 bool   *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
        {*V=TH2VREG1(T,H);break;}
    case 2:
        {*V=TH2VREG2(T,H);break;}
    case 3:
        {*V=TH2VREG3(T,H);break;}
    case 4:
        {*V=TH2VREG4(T,H);break;}
    case 5:
	case 6:
        {*V=TH2VREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  TH2VLP_KK(double T,double H, double* V,int *RANGE)
{
 bool   *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
		{ *V=TH2VLPREG1(T,H);break;}
    case 2:
        {*V=TH2VREG2(T,H);break;}
    case 3:
		{ *V=TH2VREG3(T,H);break;}
    case 4:
        {*V=TH2VREG4(T,H);break;}
    case 5:
	case 6:
        {*V=TH2VREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}


void  TH2VHP_KK(double T,double H, double* V,int *RANGE)
{
 bool   *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=1;
    switch ( *RANGE )
    {
	case 1:
		{ *V=TH2VHPREG1(T,H);break;}
    case 2:
        {*V=TH2VREG2(T,H);break;}
    case 3:
        {*V=TH2VREG3(T,H);break;}
    case 4:
        {*V=TH2VREG4(T,H);break;}
    case 5:
	case 6:
        {*V=TH2VREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  TH2XLP_KK(double T,double H, double* X, int *RANGE)
{
 bool   *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
		{ *X=0.0;break;}
    case 2:
        {*X=1.0;break;}
    case 3:
        {*X=1.0;break;}
    case 4:
        {*X=0.0;break;}
    case 5:
	case 6:
        {*X=TH2XREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  TH2XHP_KK(double T,double H, double* X, int *RANGE)
{
 bool   *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=1;
    switch ( *RANGE )
    {
	case 1:
		{ *X=0.0;break;}
    case 2:
		{ *X=1.0;break;}
    case 3:
        {*X=1.0;break;}
    case 4:
		{ *X=0.0;break;}
    case 5:
	case 6:
		{ *X=TH2XREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  TH2X_KK(double T,double H,  double *X, int *RANGE)
{
bool    *B16,b16;

    B16=&b16;
    SUBRANGEBYTH(T,H,RANGE,B16);
    if ( *B16 )
        *RANGE=6;
    switch ( *RANGE )
    {
	case 1:
        {*X=0.0;break;}
    case 2:
        {*X=1.0;break;}
    case 3:
		{ *X=1.0;break;}
    case 4:
        {*X=0.0;break;}
    case 5:
	case 6:
		{ *X=TH2XREG56(T,H);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

/////////////////////3   19:45 开始

///////////////////////////////////////////////gq修改后20:55

void  TS2PHP_KK(double T,double S,  double* P, int *RANGE)
{
    
    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*P=TS2PHPREG1(T,S);break;}
    case 2:
        {*P=TS2PREG2(T,S);break;}
    case 3:
		{*P=TS2PREG3(T,S);break;}
    case 4:
		{ *P=TS2PREG4(T,S);break;}
    case 5:
	case 6:
		{ *P=T2P_KK(T);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  TS2PLP_KK(double T,double S,  double* P, int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *P=TS2PLPREG1(T,S);break;}
    case 2:
        {*P=TS2PREG2(T,S);break;}
    case 3:
        {*P=TS2PREG3(T,S);break;}
    case 4:
        {*P=TS2PREG4(T,S);break;}
    case 5:
	case 6:
        {*P=T2P_KK(T);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  TS2P_KK(double T,double S,  double* P, int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *P=TS2PREG1(T,S);break;}
    case 2:
		{ *P=TS2PREG2(T,S);break;}
    case 3:
		{ *P=TS2PREG3(T,S);break;}
    case 4:
        {*P=TS2PREG4(T,S);break;}
    case 5:
	case 6:
        {*P=T2P_KK(T);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  TS2HHP_KK(double T,double S, double* H, int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*H=TS2HHPREG1(T,S);break;}
    case 2:
        {*H=TS2HREG2(T,S);break;}
    case 3:
		{ *H=TS2HREG3(T,S);break;}
    case 4:
        {*H=TS2HREG4(T,S);break;}
    case 5:
	case 6:
        {*H=TS2HREG56(T,S);break;}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}

void  TS2HLP_KK(double T,double S, double* H , int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*H=TS2HLPREG1(T,S);break;}
    case 2:
        {*H=TS2HREG2(T,S);break;}
    case 3:
        {*H=TS2HREG3(T,S);break;}
    case 4:
        {*H=TS2HREG4(T,S);break;}
    case 5:
	case 6:
        {*H=TS2HREG56(T,S);break;}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}

void  TS2H_KK(double T,double S,  double* H, int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *H=TS2HREG1(T,S);break;}
    case 2:
        {*H=TS2HREG2(T,S);break;}
    case 3:
		{ *H=TS2HREG3(T,S);break;}
    case 4:
		{ *H=TS2HREG4(T,S);break;}
    case 5:
	case 6:
		{ *H=TS2HREG56(T,S);break;}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}

void  TS2VHP_KK(double T,double S,  double* V, int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *V=TS2VHPREG1(T,S);break;}
    case 2:
		{ *V=TS2VREG2(T,S);break;}
    case 3:
        {*V=TS2VREG3(T,S);break;}
    case 4:
		{ *V=TS2VREG4(T,S);break;}
    case 5:
	case 6:
		{ *V=TS2VREG56(T,S);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  TS2VLP_KK(double T,double S, double* V, int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *V=TS2VLPREG1(T,S);break;}
    case 2:
		{ *V=TS2VREG2(T,S);break;}
    case 3:
		{ *V=TS2VREG3(T,S);break;}
    case 4:
		{ *V=TS2VREG4(T,S);break;}
    case 5:
	case 6:
		{ *V=TS2VREG56(T,S);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}


void  TS2V_KK(double T,double S, double* V, int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *V=TS2VREG1(T,S);break;}
    case 2:
		{ *V=TS2VREG2(T,S);break;}
    case 3:
		{ *V=TS2VREG3(T,S);break;}
    case 4:
		{ *V=TS2VREG4(T,S);break;}
    case 5:
	case 6:
		{ *V=TS2VREG56(T,S);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  TS2X_KK(double T,double S,  double* X, int *RANGE)
{
    

    SUBRANGEBYTS(T,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *X=0.0;break;}
    case 2:
		{ *X=1.0;break;}
    case 3:
		{ *X=1.0;break;}
    case 4:
		{ *X=0.0;break;}
    case 5:
	case 6:
		{ *X=TS2XREG56(T,S);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  TV2P_KK(double T,double V,  double* P, int *RANGE)
{

    SUBRANGEBYTV(T,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *P=TV2PREG1(T,V);break;}
    case 2:
		{ *P=TV2PREG2(T,V);break;}
    case 3:
		{ *P=TV2PREG3(T,V);break;}
    case 4:
		{ *P=TV2PREG4(T,V);break;}
    case 5:
	case 6:
		{ *P=T2P_KK(T);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  TV2H_KK(double T,double V, double* H, int *RANGE)
{
    
    SUBRANGEBYTV(T,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *H=TV2HREG1(T,V);break;}
    case 2:
		{ *H=TV2HREG2(T,V);break;}
    case 3:
		{ *H=TV2HREG3(T,V);break;}
    case 4:
        {*H=TV2HREG4(T,V);break;}
    case 5:
	case 6:
        {*H=TV2HREG56(T,V);break;}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}

void  TV2S_KK(double T,double V,  double* S, int *RANGE)
{
    
    SUBRANGEBYTV(T,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*S=TV2SREG1(T,V);break;}
    case 2:
		{ *S=TV2SREG2(T,V);break;}
    case 3:
		{ *S=TV2SREG3(T,V);break;}
    case 4:
		{ *S=TV2SREG4(T,V);break;}
    case 5:
	case 6:
        {*S=TV2SREG56(T,V);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}

void  TV2X_KK(double T,double V,  double* X, int *RANGE)
{
    
    SUBRANGEBYTV(T,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*X=0.0;break;}
    case 2:
		{ *X=1.0;break;}
    case 3:
        {*X=1.0;break;}
    case 4:
		{ *X=0.0;break;}
    case 5:
	case 6:
		{ *X=TV2XREG56(T,V);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  TX2P_KK(double T,double X, double* P, int *RANGE)
{
    
    SUBRANGEBYTX(T,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *P=T2P_KK(T);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  TX2H_KK(double T, double X,  double* H, int *RANGE)
{

    SUBRANGEBYTX(T,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *H=TX2HREG56(T,X);break;}
    default:
		{
        *RANGE=0;
        *H=-1.0;
        //exit(1);
		}
    }
}

void  TX2S_KK(double T, double X, double* S, int *RANGE)
{
    
    SUBRANGEBYTX(T,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *S=TX2SREG56(T,X);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}

void  TX2V_KK(double T, double X,  double* V, int *RANGE)
{
    
    SUBRANGEBYTX(T,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *V=TX2VREG56(T,X);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  H2TL_KK(double H, double* T, int *RANGE)
{
    SUBRANGEBYHL_KK(H,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
        {*T=H2TLREG56(H);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}


void  HS2P_KK(double H,double S,  double* P,int *RANGE)
{
    SUBRANGEBYHS(H,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*P=HS2PREG1(H,S);break;}
    case 2:
        {*P=HS2PREG2(H,S);break;}
    case 3:
        {*P=HS2PREG3(H,S);break;}
    case 4:
		{ *P=HS2PREG4(H,S);break;}
    case 5:
	case 6:
		{ *P=HS2PREG56(H,S);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  HS2T_KK(double H,double S, double* T,int *RANGE)
{
    SUBRANGEBYHS(H,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*T=HS2TREG1(H,S);break;}
    case 2:
        {*T=HS2TREG2(H,S);break;}
    case 3:
        {*T=HS2TREG3(H,S);break;}
    case 4:
        {*T=HS2TREG4(H,S);break;}
    case 5:
	case 6:
        {*T=HS2TREG56(H,S);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  HS2V_KK(double H,double S, double* V, int *RANGE)
{
    SUBRANGEBYHS(H,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *V=HS2VREG1(H,S);break;}
    case 2:
        {*V=HS2VREG2(H,S);break;}
    case 3:
        {*V=HS2VREG3(H,S);break;}
    case 4:
        {*V=HS2VREG4(H,S);break;}
    case 5:
	case 6:
        {*V=HS2VREG56(H,S);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  HS2X_KK(double H,double S, double* X,int *RANGE)
{
    SUBRANGEBYHS(H,S,RANGE);
    switch ( *RANGE )
    {
	case 1:
        {*X=0.0;break;}
    case 2:
		{ *X=1.0;break;}
    case 3:
        {*X=1.0;break;}
    case 4:
		{ *X=0.0;break;}
    case 5:
	case 6:
        {*X=HS2XREG56(H,S);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  HV2P_KK(double H,double V,  double* P,int *RANGE)
{
    SUBRANGEBYHV(H,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *P=HV2PREG1(H,V);break;}
    case 2:
		{ *P=HV2PREG2(H,V);break;}
    case 3:
		{ *P=HV2PREG3(H,V);break;}
    case 4:
		{ *P=HV2PREG4(H,V);break;}
    case 5:
	case 6:
		{ *P=HV2PREG56(H,V);break;}
    default:
		{
		*RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  HV2T_KK(double H,double V, double* T,int *RANGE)
{
    SUBRANGEBYHV(H,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *T=HV2TREG1(H,V);break;}
    case 2:
        {*T=HV2TREG2(H,V);break;}
    case 3:
		{ *T=HV2TREG3(H,V);break;}
    case 4:
        {*T=HV2TREG4(H,V);break;}
    case 5:
	case 6:
		{  *T=HV2TREG56(H,V);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  HV2S_KK(double H,double V, double* S, int *RANGE)
{
    SUBRANGEBYHV(H,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *S=HV2SREG1(H,V);break;}
    case 2:
        {*S=HV2SREG2(H,V);break;}
    case 3:
        {*S=HV2SREG3(H,V);break;}
    case 4:
        {*S=HV2SREG4(H,V);break;}
    case 5:
	case 6:
        {*S=HV2SREG56(H,V);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}

void  HV2X_KK(double H,double V,  double* X, int *RANGE)
{
    SUBRANGEBYHV(H,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *X=0.0;break;}
    case 2:
        {*X=1.0;break;}
    case 3:
		{ *X=1.0;break;}
    case 4:
        {*X=0.0;break;}
    case 5:
	case 6:
		{ *X=HV2XREG56(H,V);break;}
    default:
		{
        *RANGE=0;
        *X=-1.0;
        //exit(1);
		}
    }
}

void  HX2P_KK(double H, double X,  double* P,int *RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *P=HX2PREG56(H,X); break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  HX2PLP_KK(double H, double X,  double* P, int *RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *P=HX2PLPREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  HX2PHP_KK(double H, double X,  double* P, int *RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{  *P=HX2PHPREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  HX2T_KK(double H, double X,  double* T,int *RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{  *T=HX2TREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  HX2TLP_KK(double H, double X,  double* T, int *RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *T=HX2TREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  HX2THP_KK(double H,double X, double* T, int *RANGE)
{
    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *T=HX2THPREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *T=-1.0;
        //exit(1);
		}
    }
}

void  HX2S_KK(double H, double X, double* S, int *RANGE)
{

    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
        {*S=HX2SREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}
void  HX2SLP_KK(double H, double X,  double* S, int *RANGE)
{

    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *S=HX2SREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}
void  HX2SHP_KK(double H, double X, double* S, int *RANGE)
{

    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *S=HX2SHPREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *S=-1.0;
        //exit(1);
		}
    }
}

void  HX2V_KK(double H, double X,  double* V, int *RANGE)
{

    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
        {*V=HX2VREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  HX2VLP_KK(double H, double X,  double* V, int *RANGE)
{

    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *V=HX2VREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}

void  HX2VHP_KK(double H, double X,  double* V, int *RANGE)
{

    SUBRANGEBYHX(H,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *V=HX2VHPREG56(H,X);break;}
    default:
		{
        *RANGE=0;
        *V=-1.0;
        //exit(1);
		}
    }
}
void  SV2P_KK(double S,double V,  double* P, int *RANGE)
{

    SUBRANGEBYSV(S,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *P=SV2PREG1(S,V);break;}
    case 2:
        {*P=SV2PREG2(S,V);break;}
    case 3:
		{ *P=SV2PREG3(S,V);break;}
    case 4:
		{ *P=SV2PREG4(S,V);break;}
    case 5:
	case 6:
        {*P=SV2PREG56(S,V);break;}
    default:
		{
        *RANGE=0;
        *P=-1.0;
        //exit(1);
		}
    }
}

void  SV2T_KK(double S,double V,  double * T, int *RANGE)
{

    SUBRANGEBYSV(S,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *T=SV2TREG1(S,V);break;}
    case 2:
		{ *T=SV2TREG2(S,V);break;}
    case 3:
		{ *T=SV2TREG3(S,V);break;}
    case 4:
        {*T=SV2TREG4(S,V);break;}
    case 5:
	case 6:
		{ *T=SV2TREG56(S,V);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        //exit(1);
    }
    }
}

void  SV2H_KK(double S,double V, double* H, int *RANGE)
{

    SUBRANGEBYSV(S,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *H=SV2HREG1(S,V);break;}
    case 2:
		{ *H=SV2HREG2(S,V);break;}
    case 3:
		{ *H=SV2HREG3(S,V);break;}
    case 4:
		{ *H=SV2HREG4(S,V);break;}
    case 5:
	case 6:
		{ *H=SV2HREG56(S,V);break;}
    default:
    {
        *RANGE=0;
        *H=-1.0;
        //exit(1);
    }
    }
}

void  SV2X_KK(double S,double V,  double* X, int *RANGE)
{

    SUBRANGEBYSV(S,V,RANGE);
    switch ( *RANGE )
    {
	case 1:
		{ *X=0.0;break;}
    case 2:
		{  *X=1.0;break;}
    case 3:
		{ *X=1.0;break;}
    case 4:
        {*X=0.0;break;}
    case 5:
	case 6:
		{ *X=SV2XREG56(S,V);break;}
    default:
    {
        *RANGE=0;
        *X=-1.0;
        //exit(1);
    }
    }
}

void  SX2P_KK(double S, double X,  double* P, int *RANGE)
{

    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{  *P=SX2PREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *P=-1.0;
        //exit(1);
    }
    }
}

void  SX2PLP_KK(double S, double X,  double* P, int *RANGE)
{

    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{ *P=SX2PLPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *P=-1.0;
        //exit(1);
    }
    }
}

void  SX2PMP_KK(double S, double X,  double* P, int *RANGE)
{

    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{  *P=SX2PMPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *P=-1.0;
        //exit(1);
    }
    }
}


void  SX2PHP_KK(double S, double X,  double* P, int *RANGE)
{

    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    {
	case 5:
	case 6:
		{  *P=SX2PHPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *P=-1.0;
        //exit(1);
    }
    }
}


void  SX2T_KK(double S, double X,  double* T, int *RANGE)
{

    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *T=SX2TREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        //exit(1);
    }
    }
}

void  SX2TLP_KK(double S, double X,  double* T, int *RANGE)
{

    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *T=SX2TLPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        //exit(1);
    }
    }
}

void  SX2TMP_KK(double S, double X,  double* T, int *RANGE)
{

    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *T=SX2TMPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        //exit(1);
    }
    }
}

void  SX2THP_KK(double S, double X, double* T, int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *T=SX2THPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        //exit(1);
    }
    }
}

void  SX2H_KK(double S, double X,  double* H, int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{  *H=SX2HREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *H=-1.0;
        //exit(1);
    }
    }
}

void  SX2HLP_KK(double S, double X,  double* H, int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *H=SX2HLPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *H=-1.0;
        //exit(1);
    }
    }
}

void  SX2HMP_KK(double S, double X,  double* H,int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5:
	case 6:
		{ *H=SX2HMPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *H=-1.0;
        //exit(1);
    }
    }
}

void  SX2HHP_KK(double S, double X,  double* H,int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *H=SX2HHPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *H=-1.0;
        //exit(1);
    }
    }
}

void  SX2V_KK(double S, double X, double* V,int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *V=SX2VREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *V=-1.0;
        //exit(1);
    }
    }
}

void  SX2VLP_KK(double S, double X, double* V,int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *V=SX2VLPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *V=-1.0;
        //exit(1);
    }
    }
}

void  SX2VMP_KK(double S, double X, double* V,int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{  *V=SX2VMPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *V=-1.0;
        //exit(1);
    }
    }
}

void  SX2VHP_KK(double S, double X, double* V, int *RANGE)
{
    SUBRANGEBYSX(S,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *V=SX2VHPREG56(S,X);break;}
    default:
    {
        *RANGE=0;
        *V=-1.0;
        //exit(1);
    }
    }
}

void  V2TG_KK(double V, double* T,int *RANGE)
{
    SUBRANGEBYVG_KK(V,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *T=V2TGREG56(V);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        //exit(1);
    }
    }
}

void  VX2P_KK(double V, double X , double* P,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *P=VX2PREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *P=-1.0;
        return;
		//exit(1);
    }
    }
    if (*P<=T2P_KK(T350C)+DeltaVal )
        *RANGE=6;
    else
        *RANGE=5;
}

void  VX2PLP_KK(double V, double X , double* P,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{  *P=VX2PLPREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *P=-1.0;
        return;
		//exit(1);
    }
    }
    if ( *P<T2P_KK(T350C)+DeltaVal )
        *RANGE=6;
    else
        *RANGE=5;
}

void  VX2PHP_KK(double V, double X , double* P,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{  *P=VX2PHPREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *P=-1.0;
        return;
		//exit(1);
    }
    }

    if ( *P<T2P_KK(T350C)+DeltaVal )
        *RANGE=6;
    else
        *RANGE=5;
}

void  VX2T_KK(double V, double X , double* T,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *T=VX2TREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        return;
		//exit(1);
    }
    }
    if ( *T<T350C+DeltaVal )
        *RANGE=6;
    else
        *RANGE=5;
}

void  VX2TLP_KK(double V, double X , double* T,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *T=VX2TLPREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        return;
		//exit(1);
    }
    }
    if ( *T<T350C+DeltaVal )
        *RANGE=6;
    else
        *RANGE=5;
}

void  VX2THP_KK(double V, double X , double* T,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{   *T=VX2THPREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *T=-1.0;
        return;
		//exit(1);
    }
    }
    if ( *T<T350C+DeltaVal )
        *RANGE=6;
    else
        *RANGE=5;
}

void  VX2H_KK(double V, double X , double* H,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *H=VX2HREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *H=-1.0;
        //exit(1);
    }
    }
    SUBRANGEBYHV(*H,V,RANGE);
}

void  VX2HLP_KK(double V, double X , double* H,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch (*RANGE )
    { 
	case 5: 
	case 6:
		{  *H=VX2HLPREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *H=-1.0;
        return;
		//exit(1);
    }
    }
    SUBRANGEBYHV(*H,V,RANGE);
}

void  VX2HHP_KK(double V, double X , double* H,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *H=VX2HHPREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *H=-1.0;
        return;
		//exit(1);
    }
    }
    SUBRANGEBYHV(*H,V,RANGE);
}

void  VX2S_KK(double V, double X , double* S,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{  *S=VX2SREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *S=-1.0;
        return;
		//exit(1);
    }
    }
    SUBRANGEBYSV(*S,V,RANGE);
}

void  VX2SLP_KK(double V, double X , double* S,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *S=VX2SLPREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *S=-1.0;
        return;
		//exit(1);
    }
    }
    SUBRANGEBYSV(*S,V,RANGE);
}

void  VX2SHP_KK(double V, double X , double* S,int *RANGE)
{
    SUBRANGEBYVX(V,X,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *S=VX2SHPREG56(V,X);break;}
    default:
    {
        *RANGE=0;
        *S=-1.0;
        return;
		//exit(1);
    }
    }
    SUBRANGEBYSV(*S,V,RANGE);
}


void SWETAPTV(double P, double T ,double V, double* ETA, int *RANGE)
{
double    T1 , v1 ;
double    HTI , HTI1, HTI2, HTI3;
int     I , J ;
double   ETA0 , ETA1 ;


    T1 = (T + 273.15) / 647.27;
    v1 = (1 / V) / 317.763;

    HTI = 0;    
	for(I=0;I<=3;I++)
        HTI = HTI + EtaH[I] / exp(log(T1)*I);
    ETA0 = sqrt(T1) / HTI;

    HTI1 = 0;
    
	for(I=0;I<=5;I++)  
		for(J=0;J<=6;J++)        
        {
            HTI2= exp(log(fabs(1 / T1 - 1)) *I);
            if ( (I%2==1) && (1<T1) ) 
				HTI2= -HTI2;
            HTI3=exp(log(fabs(v1 - 1))*J);
            if ( (J%2==1) && (v1<1) ) 
				HTI3=-HTI3;
            HTI1 = HTI1 + EtaHH[I][J] *HTI2  *HTI3 ;
        }
    ETA1 = exp(v1 * HTI1);

    *ETA = ETA0 * ETA1 * 1 * 55.071E-6;
}

void SWETAPT(double P,double T, double* ETA, int *RANGE)
{
double    *V,v;

    V=&v;
    PT2V_KK(P,T,V,RANGE);
    SWETAPTV(P, T ,*V,ETA, RANGE);
}


void SWETAPL(double P, double* ETA,int *RANGE)
{
double    T,*V,v;

    V=&v;
    T=P2T_KK(P);
    P2VL_KK(P,V,RANGE);
    SWETAPTV(P, T ,*V,ETA, RANGE);
}

void SWETAPG(double P, double* ETA,int *RANGE)
{
double    T,*V,v;

    V=&v; 
    T=P2T_KK(P);
    P2VG_KK(P,V,RANGE);
    SWETAPTV(P, T ,*V,ETA, RANGE);
}

void SWETATL(double T, double* ETA,int *RANGE)
{
 double   P,*V,v;
 
    V=&v;
    P=T2P_KK(T);
    P2VL_KK(P,V,RANGE);
    SWETAPTV(P, T ,*V,ETA, RANGE);
}

void SWETATG(double T, double* ETA,int *RANGE)
{
 double   P,*V,v;

    V=&v;
    P=T2P_KK(T);
    P2VG_KK(P,V,RANGE);
    SWETAPTV(P, T ,*V,ETA, RANGE);
}


void  SWRAMDPTV(double P,double T,double V, double* RAMD, int *RANGE)
{
double    ZT,ZT2,ZV,ZV2,R0,RB,RD,RD1,RD2,RD3,DT,Q,R,S;
int     I;

    ZT=(T+273.15)/647.27;
    ZV=V/0.003147;
    ZT2=ZT*ZT;
    ZV2=ZV*ZV;

    R0=0.0;   
	for (I=3;I>=0;I--)
    {
        R0=R0*ZT+RamdAL[I];
    }
    R0=sqrt(ZT)*R0;

    RB=RamdBL[0]+RamdBL[1]/ZV+RamdBL[2]*exp(RamdB[1]*(1.0/ZV+RamdB[2])*(1.0/ZV+RamdB[2]));

    RD1=(RamdDL[1]/exp(5*log(ZT2))+RamdDL[2])/exp(0.9*log(ZV2))*exp(RamdCL[1]*(1-1.0 / exp(1.4*log(ZV2)) ));

    RD3=RamdDL[4]*exp(RamdCL[2]*exp(0.75*log(ZT2))+RamdCL[3]*exp(2.5*log(ZV2)));

    DT=fabs(ZT-1.0)+RamdCL[4];
    Q=2.0+RamdCL[5]/exp(0.3*log(DT*DT));
    R=Q+1.0;
    if ( ZT>=1.0 )
        S=1.0/DT;
    else
        S=RamdCL[6]/exp(0.3*log(DT*DT));
    RD2=RamdDL[3]*S/exp(Q/2.0*log(ZV2))*exp(Q/R*(1.0-1.0/exp(R/2.0*log(ZV2))));

    RD=RD1+RD2+RD3;

    *RAMD=(R0+RB+RD);
}


void  SWRAMDPT(double P,double T, double* RAMD, int *RANGE)
{
 double   *V,v;

    V=&v; 
    PT2V_KK(P,T,V,RANGE); 
    SWRAMDPTV(P,T,*V,RAMD,RANGE);
}


void  SWRAMDPL(double P, double* RAMD,int *RANGE)
{
double    T,*V,v;

    V=&v;
    T=P2T_KK(P);
    P2VL_KK(P,V,RANGE); 
    SWRAMDPTV(P,T,*V,RAMD,RANGE);
}

void  SWRAMDPG(double P, double* RAMD,int *RANGE)
{
double    T,*V,v;

    V=&v;
    T=P2T_KK(P);
    P2VG_KK(P,V,RANGE); 
    SWRAMDPTV(P,T,*V,RAMD,RANGE);
}



void  SWRAMDTL(double T, double* RAMD,int *RANGE)
{
double    P,*V,v;

    V=&v;
    P=T2P_KK(T);
    P2VL_KK(P,V,RANGE); 
    SWRAMDPTV(P,T,*V,RAMD,RANGE);
}

void  SWRAMDTG(double T, double* RAMD,int *RANGE)
{
double    P,*V,v;

    V=&v;
    P=T2P_KK(T);
    P2VG_KK(P,V,RANGE); 
    SWRAMDPTV(P,T,*V,RAMD,RANGE);
}

void  PT2ETA(double P,double T,  double* ETA, int *RANGE)
{

    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4: 
	case 5: 
	case 6:     
    {
        SWETAPT(P,T,ETA,RANGE);    break;       
    }
    default:
    {
        *RANGE=0;
        *ETA=-1.0;
    }
    }
}

void  PT2U_KK(double  P,double T, double* U, int  *RANGE)
{
 double   *Eta,*V,v,eta;

    V=&v;
	Eta=&eta;
    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4: 
	case 5: 
	case 6:     
    {
        PT2V_KK(P,T,V,RANGE);
        SWETAPT(P,T,Eta,RANGE);
        *U =*Eta*(*V);    
		break;
    }
    default:
    {
        *RANGE=0;
        *U=-1.0;
    }
    }
}



void  PT2RAMD_KK(double P,double T, double* RAMD,int *RANGE)
{

    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4: 
	case 5: 
	case 6:
    {
        SWRAMDPT(P,T,RAMD,RANGE);  break;        
    }
    default:
    {
        *RANGE=0;
        *RAMD=-1.0;
    }
    }
}


void  PT2PRN_KK(double P,double T, double* PRN, int *RANGE)
{
double    *CP,*ETA,*RAMD,cp,eta,ramd;

    CP=&cp;
	ETA=&eta;
	RAMD=&ramd;
    SUBRANGEBYPT(P,T,RANGE);
    switch ( *RANGE )
    {
	case 1:
	case 2:
	case 3:
	case 4: 
	case 5: 
	case 6:     
    {
        if ( (P<PC_Water) && (T<TC_Water) )
        {
            if ( fabs(1-T2P_KK(T)/P)<1E-8 )
            {
                P2PRNL(P,PRN,RANGE);
                return;
				//exit(1);
            }
        }

        PT2CP_KK(P,T,CP,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        PT2ETA(P,T,ETA,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        PT2RAMD_KK(P,T,RAMD,RANGE);
        if ( (*RANGE!=0) )
        {
            if ( fabs(*RAMD)<1E-10 )
            {
                *PRN=0.0;
                *RANGE=0;
                return;
				//exit(1);
            }
            else
			{
                *PRN=(*ETA)*(*CP)/(*RAMD)*1E3;
			}
        }
        else
		{
            *PRN=0.0;
		}
		break;
    }
    default:
    {
        *RANGE=0;
        *ETA=-1.0;
    }
    }
}



void  P2ETAL(double P, double* ETA,int *RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5:
	case 6:     
    {
        SWETAPL(P,ETA,RANGE);  break;         
    }
    default:
    {
        *RANGE=0;
        *ETA=-1.0;
    }
    }
}


void  P2ETAG(double P, double* ETA,int *RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        SWETAPG(P,ETA,RANGE);   break;        
    }
    default:
    {
        *RANGE=0;
        *ETA=-1.0;
    }
    }
}


void  P2UL(double P, double* U, int *RANGE)
{
 double   *ETA,V,eta;

    ETA=&eta; 
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        V=P2VLREG56(P);
        SWETAPL(P,ETA,RANGE);           
        *U=V*(*ETA);
		break;
    }
    default:
    {
        *RANGE=0;
        *U=-1.0;
    }
    }
}


void  P2UG(double P, double* U, int *RANGE)
{
double    *ETA,V,eta;

    ETA=&eta;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        V=P2VGREG56(P);
        SWETAPG(P,ETA,RANGE);           
        *U=V*(*ETA);    
		break;
    }
    default:
    {
        *RANGE=0;
        *U=-1.0;
    }
    }
}

void  P2RAMDL(double P, double* RAMD, int *RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        SWRAMDPL(P,RAMD,RANGE); break;          
    }
    default:
    {
        *RANGE=0;
        *RAMD=-1.0;
    }
    }
}


void  P2RAMDG(double P, double* RAMD,int *RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        SWRAMDPG(P,RAMD,RANGE);   break;        
    }
    default:
    {
        *RANGE=0;
        *RAMD=-1.0;
    }
    }
}


void  P2PRNL(double P, double* PRN, int *RANGE)
{
double   *CP,*ETA,*RAMD,cp,eta,ramd;

    CP=&cp;
	ETA=&eta;
	RAMD=&ramd;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        P2CPL(P,CP,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        P2ETAL(P,ETA,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        P2RAMDL(P,RAMD,RANGE);
        if ( (*RANGE!=0) )
        {
            if ( fabs(*RAMD)<1E-10 )
            {
                *PRN=0.0;
                *RANGE=0;
                return;
				//exit(1);
            }
            else
			{
                *PRN=(*ETA)*(*CP)/(*RAMD)*1E3;
			}
        }
        else
		{
            *PRN=0.0;
		}
		break;
        
    }
    default:
    {
        *RANGE=0;
        *PRN=-1.0;
    }
    }
}


void  P2PRNG(double P, double* PRN,int *RANGE)
{
double   *CP,*ETA,*RAMD,cp,eta,ramd;

    CP=&cp;
	ETA=&eta;
	RAMD=&ramd;
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        P2CPG(P,CP,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        P2ETAG(P,ETA,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        P2RAMDG(P,RAMD,RANGE);
        if ( (*RANGE!=0) )
        {
            if ( fabs(*RAMD)<1E-10 )
            {
                *PRN=0.0;
                *RANGE=0;
                return;
				//exit(1);
            }
            else
			{
                *PRN=(*ETA)*(*CP)/(*RAMD)*1E3;
			}
        }
        else
		{
            *PRN=0.0;
		}
		break;
        
    }
    default:
    {
        *RANGE=0;
        *PRN=-1.0;
    }
    }
}

double EPSSATLIQP(double  P)
{
double    T ;

    T = P2T_KK(P);                    
    return TV2EPS_KK(T,T2VLREG56(T));
}

void  P2EPSL(double P, double* EPS,int *RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        *EPS=EPSSATLIQP(P);break;
    }
    default:
    {
        *RANGE=0;
        *EPS=-1.0;
    }
    }
}

double EPSSATVAPP(double  P)
{
double    T ;

    T = P2T_KK(P);
    return TV2EPS_KK(T,T2VGREG56(T));
}

void  P2EPSG(double P, double* EPS, int *RANGE)
{

    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        *EPS=EPSSATVAPP(P);break;
    }
	default:
    {
        *RANGE=0;
        *EPS=-1.0;
    }
    }
}



void  T2ETAL(double T,  double* ETA, int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        SWETATL(T,ETA,RANGE);   break;        
    }
    default:
    {
        *RANGE=0;
        *ETA=-1.0;
    }
    }
}



void  T2ETAG(double T, double* ETA, int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        SWETATG(T,ETA,RANGE);    break;       
    }
    default:
    {
        *RANGE=0;
        *ETA=-1.0;
    }
    }
}


void  T2UL(double T,  double* U, int *RANGE)
{
double   V,*ETA,eta;

    ETA=&eta;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        V=T2VLREG56(T);
        SWETATL(T,ETA,RANGE);           
        *U=V*(*ETA);
		break;
    }
    default:
    {
        *RANGE=0;
        *U=-1.0;
    }
    }
}


void  T2UG(double T,  double* U, int *RANGE)
{
 double   V,*ETA,eta;

    ETA=&eta; 
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        V=T2VGREG56(T);
        SWETATG(T,ETA,RANGE);
        *U=V*(*ETA);
		break;
    }
    default:
    {
        *RANGE=0;
        *U=-1.0;
    }
    }
}


void  T2RAMDL(double T, double* RAMD, int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        SWRAMDTL(T,RAMD,RANGE);break;
    }
    default:
    {
        *RANGE=0;
        *RAMD=-1.0;
    }
    }
}


void  T2RAMDG(double T, double* RAMD,int *RANGE)
{

    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        SWRAMDTG(T,RAMD,RANGE);break;
    }
    default:
    {
        *RANGE=0;
        *RAMD=-1.0;
    }
    }
}


void  T2PRNL(double T,  double* PRN, int *RANGE)
{
double   *CP,*ETA,*RAMD,cp,eta,ramd;

    CP=&cp;
	ETA=&eta;
	RAMD=&ramd;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:     
    {
        T2CPL(T,CP,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        T2ETAL(T,ETA,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        T2RAMDL(T,RAMD,RANGE);
        if ( (*RANGE!=0) )
        {
            if ( fabs(*RAMD)<1E-10 )
            {
                *PRN=0.0;
                *RANGE=0;
                return;
				//exit(1);
            }
            else
			{
                *PRN=(*ETA)*(*CP)/(*RAMD)*1E3;
			}
        }
        else
		{
            *PRN=0.0;
		}
		break;
    }
    default:
    {
        *RANGE=0;
        *PRN=-1.0;
    }
    }
}


void  T2PRNG(double T,  double* PRN,int *RANGE)
{
double   *CP,*ETA,*RAMD,cp,eta,ramd;

    CP=&cp;
	ETA=&eta;
	RAMD=&ramd;
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
    {
        T2CPG(T,CP,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        T2ETAG(T,ETA,RANGE);
        if ( *RANGE==0 )
        {
            *PRN=0.0;
            return;
			//exit(1);
        }

        T2RAMDG(T,RAMD,RANGE);
        if ( (*RANGE!=0) )
        {
            if ( fabs(*RAMD)<1E-10 )
            {
                *PRN=0.0;
                *RANGE=0;
                return;
				//exit(1);
            }
            else
			{
                *PRN=(*ETA)*(*CP)/(*RAMD)*1E3;
			}
        }
        else
		{
            *PRN=0.0;
		}
		break;
    }
    default:
    {
        *RANGE=0;
        *PRN=-1.0;
    }
    }
}






double FUNC_PSAT(double Temp,double  PrmA)
{
    return T2P_KK(Temp)-PrmA;
}


double FUNC_T2PBOUND23(double Temp,double  PrmA)
{
    return T2PBOUND23(Temp)-PrmA;
}

double FUNC_T2HBOUND23(double Temp,double  PrmA)
{
    return PT2HREG2(T2PBOUND23(Temp),Temp)-PrmA;
}

double FUNC_T2VBOUND23(double Temp,double  PrmA)
{
    return PT2VREG2(T2PBOUND23(Temp),Temp)-PrmA;
}


double FUNC_TV2PREG34(double Temp,double  PrmA)
{
    if ( Temp>TC_Water-DeltaVal )
        return TV2PREG3(Temp, PrmA)-100;
    else
        return TV2PREG4(Temp, PrmA)-100;
}

double FUNC_TV2P100REG4(double Temp,double  PrmA)
{
    return TV2PREG4(Temp, PrmA)-100;
}

double FUNC_TV2P100REG3(double Temp,double  PrmA)
{
    return TV2PREG3(Temp, PrmA)-100;
}

double FUNC_TS2P100REG4(double Temp,double  PrmA)
{
    return TS2PREG4(Temp, PrmA)-100;
}

double FUNC_TS2P100REG3(double Temp,double  PrmA)
{
    return TS2PREG3(Temp, PrmA)-100;
}


double FUNC_TH2P100REG3(double Temp,double  PrmA)
{
    return PT2HREG3(100,Temp)-PrmA;
}


double FUNC_TV2PSATREG3(double Temp,double  PrmA)
{
    return TV2PREG3(Temp, PrmA)-T2P_KK(Temp);
}


double FUNC_TH2PSATREG3(double Temp,double  PrmA)
{
    return TH2PREG3(Temp, PrmA)-T2P_KK(Temp);
}


double FUNC_H2TLREG56(double Temp,double  PrmA)
{
    return T2HLREG56(Temp)-PrmA;
}


double FUNC_H2TGREG56(double Temp,double  PrmA)
{
    return T2HGREG56(Temp)-PrmA;
}



double FUNC_S2TLREG56(double Temp,double  PrmA)
{
    return T2SLREG56(Temp)-PrmA;
}



double FUNC_S2TGREG56(double Temp,double  PrmA)
{
    return T2SGREG56(Temp)-PrmA;
}



double FUNC_V2TLREG56(double Temp,double  PrmA)
{
    return T2VLREG56(Temp)-PrmA;
}


double FUNC_V2TGREG56(double Temp,double  PrmA)
{
    return T2VGREG56(Temp)-PrmA;
}


double FUNC_TV2PSATREG4(double Temp,double  PrmA)
{
    return TV2PREG4(Temp, PrmA)-T2P_KK(Temp);
}


double FUNC_TS2PSATREG4(double Temp,double  PrmA)
{
    return T2SLREG56(Temp)-PrmA;
}


double FUNC_T2SLREG6(double Temp,double  PrmA)
{
    return T2SLREG56(Temp)-PrmA;
}


double FUNC_T2VLREG6(double Temp,double  PrmA)
{
    return T2VLREG56(Temp)-PrmA;
}


double FUNC_T2HGREG6(double Temp,double  PrmA)
{
    return T2HGREG56(Temp)-PrmA;
}



double FUNC_T2VGREG6(double Temp,double  PrmA)
{
    return T2VGREG56(Temp)-PrmA;
}

double FUNC_T2SGREG6(double Temp,double  PrmA)
{
    return T2SGREG56(Temp)-PrmA;
}

double ZBRENT1_KK(double  X1 , double X2 , double TOL ,double PrmA , int Nbr )
{
double    A , B , C , D , E , FA , FB , FC ,ZZ;
double    TOL1 , XM , S , P , Q , R , AAA ;
int     ITMAX , ITER ;
double    EPS ;

    ITMAX = 100;
    EPS = 0.00000003;
    A = X1;
    B = X2;
    C =0;
    D =0;
    E =0;

 switch ( Nbr )
    {
	case 3200:
    {
        FA = FUNC_T2PBOUND23(A, PrmA);
        FB = FUNC_T2PBOUND23(B, PrmA);break; 
    }
    case 3201:
    {
        FA = FUNC_T2HBOUND23(A, PrmA);
        FB = FUNC_T2HBOUND23(B, PrmA);break; 
    }
    case 3202:
    {
        FA = FUNC_T2VBOUND23(A, PrmA);
        FB = FUNC_T2VBOUND23(B, PrmA);break; 
    }
    case 3400:
    {
        FA = FUNC_TV2PREG34(PrmA,A);
        FB = FUNC_TV2PREG34(PrmA,B);break; 
    }
    case 3401:
    {

        FA = FUNC_TV2P100REG4(A,PrmA);
        FB = FUNC_TV2P100REG4(B,PrmA);break; 
    }
    case 3402:
    {

        FA = FUNC_TV2P100REG3(A,PrmA);
        FB = FUNC_TV2P100REG3(B,PrmA);break; 
    }
    case 3403:
    {

        FA = FUNC_TS2P100REG3(A,PrmA);
        FB = FUNC_TS2P100REG3(B,PrmA);break; 
    }
    case 3404:
    {
        FA = FUNC_TS2P100REG4(A,PrmA);
        FB = FUNC_TS2P100REG4(B,PrmA);break; 
    }
    case 3405:
    {
        FA = FUNC_TH2P100REG3(A,PrmA);
        FB = FUNC_TH2P100REG3(B,PrmA);break; 
    }
    case 5300:
    {
        FA = FUNC_TV2PSATREG3(PrmA,A);
        FB = FUNC_TV2PSATREG3(PrmA,B);break; 
    }
    case 5301:
    {
        FA = FUNC_TV2PSATREG3(A,PrmA);
        FB = FUNC_TV2PSATREG3(B,PrmA);break; 
    }
    case 5302:
    {
        FA = FUNC_TH2PSATREG3(A,PrmA);
        FB = FUNC_TH2PSATREG3(B,PrmA);break; 
    }
    case 5400:
    {
        FA = FUNC_TV2PSATREG4(PrmA,A);
        FB = FUNC_TV2PSATREG4(PrmA,B);break; 
    }
    case 5401:
    {
        FA = FUNC_TV2PSATREG4(A,PrmA);
        FB = FUNC_TV2PSATREG4(B,PrmA);break; 
    }
    case 5402:
    {
        FA = FUNC_TS2PSATREG4(A,PrmA);
        FB = FUNC_TS2PSATREG4(B,PrmA);break; 
    }
    case 5600:
    {
        FA = FUNC_H2TLREG56(A,PrmA);
        FB = FUNC_H2TLREG56(B,PrmA);break; 
    }
    case 5601:
    {
        FA = FUNC_H2TGREG56(A, PrmA);
        FB = FUNC_H2TGREG56(B, PrmA);break; 
    }
    case 5602:
    {
        FA = FUNC_S2TLREG56(A, PrmA);
        FB = FUNC_S2TLREG56(B, PrmA);break; 
    }
    case 5603:
    {
        FA = FUNC_S2TGREG56(A, PrmA);
        FB = FUNC_S2TGREG56(B, PrmA);break; 
    }
    case 5604:
    {
        FA = FUNC_V2TLREG56(A, PrmA);
        FB = FUNC_V2TLREG56(B, PrmA);break; 
    }
    case 5605:
    {
        FA = FUNC_V2TGREG56(A, PrmA);
        FB = FUNC_V2TGREG56(B, PrmA);break; 
    }
    case 5606:
    {
        FA = FUNC_PSAT(A, PrmA);
        FB = FUNC_PSAT(B, PrmA);break; 
    }
    case 6100:
    {
        FA = FUNC_T2SLREG6(A, PrmA);
        FB = FUNC_T2SLREG6(B, PrmA);break; 
    }
    case 6101:
    {
        FA = FUNC_T2VLREG6(A, PrmA);
        FB = FUNC_T2VLREG6(B, PrmA);break; 
    }
    case 6200:
    {
        FA = FUNC_T2HGREG6(A, PrmA);
        FB = FUNC_T2HGREG6(B, PrmA);break; 
    }
    case 6201:
    {
        FA = FUNC_T2VGREG6(A, PrmA);
        FB = FUNC_T2VGREG6(B, PrmA);break; 
    }
    case 6202:
    {
        FA = FUNC_T2SGREG6(A, PrmA);
        FB = FUNC_T2SGREG6(B, PrmA);break; 
    }
    default:
    {
        FA =0.0;
        FB =0.0;  
    }
    } //switch ( Nbr )

    if ( FB * FA > 0.0 )
    {
    }
    
	FC = FB;   
	for(ITER=1;ITER<=ITMAX;ITER++)
    {
        if ( FB * FC > 0 )
        {
            C = A;
            FC = FA;
            D = B - A;
            E = D;
        }
        if ( fabs(FC) < fabs(FB) )
        {
            A = B;
            B = C;
            C = A;
            FA = FB;
            FB = FC;
            FC = FA;
        }
        TOL1 = 2.0 * EPS * fabs(B) + 0.5 * TOL;
        XM = 0.5 * (C - B);
        if ( (fabs(XM) <= TOL1) || (FB == 0) )
        {
            return  B;
            //exit(1);
        }
        if ( (fabs(E) >= TOL1) && (fabs(FA) > fabs(FB)) )
        {
            S = FB / FA;
            if ( A == C )
            {
                P = 2.0 * XM * S;
                Q = 1 - S;
            }
            else
            {
                Q = FA / FC;
                R = FB / FC;
                P = S * (2.0 * XM * Q * (Q - R) - (B - A) * (R - 1.0));
                Q = (Q - 1.0) * (R - 1.0) * (S - 1.0);
            }
            if ( P > 0.0 ) Q = -Q;
            P = fabs(P);
            if ( 3.0 * XM * Q - fabs(TOL1 * Q) < fabs(E * Q) )
                AAA = 3.0 * XM * Q - fabs(TOL1 * Q);
            else
                AAA = fabs(E * Q);
            if ( 2.0 * P < AAA )
            {
                E = D;
                D = P / Q;
            }
            else
            {
                D = XM;
                E = D;
            }
        }
        else
        {
            D = XM;
            E = D;
        }
        A = B;
        FA = FB;
        if ( fabs(D) > TOL1 )
            B = B + D ;
        else
        {
            if ( XM>=0 )
                ZZ=1;
            else
                ZZ=-1;
            B=B + fabs(TOL1) * ZZ;
        }

  switch ( Nbr )
	{
       case  3200:
        {
            FB = FUNC_T2PBOUND23(B, PrmA); break;
        }
       case  3201:
        {
            FB = FUNC_T2HBOUND23(B, PrmA); break;
        }
        case 3202:
        {
            FB = FUNC_T2VBOUND23(B, PrmA); break;
        }
        case 3400:
        {
            FB = FUNC_TV2PREG34(PrmA,B); break;
        }
        case 3401:
        {
            FB = FUNC_TV2P100REG4(B,PrmA); break;
        }
        case 3402:
        {
            FB = FUNC_TV2P100REG3(B,PrmA); break;
        }
        case 3403:
        {
            FB = FUNC_TS2P100REG3(B,PrmA); break;
        }
        case 3404:
        {
            FB = FUNC_TS2P100REG4(B,PrmA); break;
        }
       case  3405:
        {
            FB = FUNC_TH2P100REG3(B,PrmA); break;
        }
        case 5300:
        {
            FB = FUNC_TV2PSATREG3(PrmA,B); break;
        }
        case 5301:
        {
            FA = FUNC_TV2PSATREG3(A,PrmA);
            FB = FUNC_TV2PSATREG3(B,PrmA);
			break;
        }
       case  5302:
        {
            FB = FUNC_TH2PSATREG3(B,PrmA); break;
        }
       case  5400:
        {
            FB = FUNC_TV2PSATREG4(PrmA,B); break;
        }
      case   5401:
        {
            FB = FUNC_TV2PSATREG4(B,PrmA); break;
        }
      case   5402:
        {
            FB = FUNC_TS2PSATREG4(B,PrmA); break;
        }
      case   5600:
        {
            FB = FUNC_H2TLREG56(B,PrmA); break;
        }
      case   5601:
        {
            FB = FUNC_H2TGREG56(B, PrmA); break;
        }
       case  5602:
        {
            FB = FUNC_S2TLREG56(B, PrmA); break;
        }
        case 5603:
        {
            FB = FUNC_S2TGREG56(B, PrmA); break;
        }
        case 5604:
        {
            FB = FUNC_V2TLREG56(B, PrmA); break;
        }
       case  5605:
        {
            FB = FUNC_V2TGREG56(B, PrmA); break;
        }
       case  5606:
        {
            FB = FUNC_PSAT(B, PrmA); break;
        }
        case 6100:
        {
            FB = FUNC_T2SLREG6(B, PrmA); break;
        }
        case 6101:
        {
            FB = FUNC_T2VLREG6(B, PrmA); break;
        }
        case 6200:
        {
            FB = FUNC_T2HGREG6(B, PrmA); break;
        }
        case 6201:
        {
            FB = FUNC_T2VGREG6(B, PrmA); break;
        }
        case 6202:
        {
            FB = FUNC_T2SGREG6(B, PrmA); break;
        }
        default:
        {
            FB =0.0; break;
        }
        }//switch ( Nbr )
    }
    return  B;
}


double FUNC_PT2HREG1(double P,double T,double H)
{
    return PT2HREG1 (P,T)-H;
}

double FUNC_PT2SREG1(double P,double T,double S)
{
    return PT2SREG1 (P,T)-S;
}

double FUNC_PT2VREG1(double P,double T,double V)
{
    return PT2VREG1 (P,T)-V;
}

double FUNC_PS2HREG1(double P,double S,double H)
{
    return PS2HREG1(P,S)-H;
}

double FUNC_PV2HREG1(double P,double V,double H)
{
    return PV2HREG1(P,V)-H;
}


double FUNC_PV2SREG1(double P,double V,double S)
{
    return PV2SREG1(P,V)-S;
}


double FUNC_PH2VREG1(double P,double H,double V)
{
    return PH2VREG1(P,H)-V;
}

double FUNC_PS2VREG1(double P,double S,double V)
{
    return PS2VREG1(P,S)-V;
}


double FUNC_PT2HREG2(double P,double T,double H)
{
    return PT2HREG2(P,T)-H;
}

double FUNC_PT2SREG2(double P,double T,double S)
{
    return PT2SREG2 (P,T)-S;
}

double FUNC_PT2VREG2(double P,double T,double V)
{
    return PT2VREG2 (P,T)-V;
}

double FUNC_PH2SREG2(double P,double H,double S)
{
    return PH2SREG2(P,H)-S;
}

double FUNC_PH2VREG2(double P,double H,double V)
{
    return PH2VREG2(P,H)-V;
}


double FUNC_PV2SREG2(double P,double V,double S)
{
    return PV2SREG2(P,V)-S;
}

double FUNC_TV2PREG3(double T,double V,double P)
{
    return TV2PREG3(T,V)-P;
}

double FUNC_TV2PREG4(double T,double V,double P)
{
    return TV2PREG4(T,V)-P;
}


double FUNC_TV2HREG4(double T,double V,double H)
{
    return TV2HREG4(T,V)-H;
}

double FUNC_TV2SREG4(double T,double V,double S)
{
    return TV2SREG4(T,V)-S;
}


double FUNC_PT2HREG4(double P,double T,double H)
{
    return PT2HREG4(P,T)-H;
}

double FUNC_PT2SREG4(double P,double T,double S)
{
    return PT2SREG4(P,T)-S;
}

double FUNC_TS2HREG4(double T,double S,double H)
{
    return TS2HREG4(T,S)-H;
}

double FUNC_TV2HREG3(double T,double V,double H)
{
    return TV2HREG3(T,V)-H;
}

double FUNC_TV2SREG3(double T,double V,double S)
{
    return TV2SREG3(T,V)-S;
}


double FUNC_PT2HREG3(double P,double T,double H)
{
    return PT2HREG3(P,T)-H;
}

double FUNC_PT2SREG3(double P,double T,double S)
{
    return PT2SREG3(P,T)-S;
}

double FUNC_TH2SREG3(double T,double H,double S)
{
    return TH2SREG3(T,H)-S;
}

double FUNC_PS2HREG56(double P,double S,double H)
{
    return PS2HREG56(P,S)-H;
}

double FUNC_PV2HREG56(double P,double V,double H)
{
    return PV2HREG56(P,V)-H;
}


double FUNC_TX2HREG56(double T,double X,double H)
{
    return TX2HREG56(T,X)-H;
}



double FUNC_PS2VREG56(double P,double S,double V)
{
    return PS2VREG56(P,S)-V;
}



double FUNC_TX2SREG56(double T,double X,double S)
{
    return TX2SREG56(T,X)-S;
}

double FUNC_TX2VREG56(double T,double X,double V)
{
    return TX2VREG56(T,X)-V;
}


double ZBRENT2_KK(double X1,double X2, double TOL ,double PrmA ,double PrmB ,int Nbr )
{
double    A , B , C , D , E , FA , FB , FC ;
double    ZZ,TOL1, XM , S , P , Q , R , AAA,EPS ;
int     ITMAX , ITER ;

    ITMAX = 100;
    EPS = 0.00000003;
    A = X1;
    B = X2;
    C=0.0;
    D=0.0;
    E=0.0;

 switch ( Nbr )
{
    case 100:
    {

        FA = FUNC_PT2HREG1(PrmA, A, PrmB);
        FB = FUNC_PT2HREG1(PrmA, B, PrmB); break;
    }
    case 101:
    {

        FA = FUNC_PT2SREG1(PrmA, A, PrmB);
        FB = FUNC_PT2SREG1(PrmA, B, PrmB); break;
    }
    case 102:
    {

        FA = FUNC_PT2VREG1(PrmA, A, PrmB);
        FB = FUNC_PT2VREG1(PrmA, B, PrmB); break;
    }
    case 103:
    {

        FA = FUNC_PT2HREG1(A, PrmA, PrmB);
        FB = FUNC_PT2HREG1(B, PrmA, PrmB); break;
    }
    case 104:
    {

        FA = FUNC_PT2SREG1(A, PrmA, PrmB);
        FB = FUNC_PT2SREG1(B, PrmA, PrmB); break;
    }
    case 105:
    {
        FA = FUNC_PT2VREG1(A, PrmA, PrmB);
        FB = FUNC_PT2VREG1(B, PrmA, PrmB); break;
    }
   case  106:
    {
        FA = FUNC_PS2HREG1(A, PrmA, PrmB);
        FB = FUNC_PS2HREG1(B, PrmA, PrmB); break;
    }
    case 107:
    {
        FA = FUNC_PV2HREG1(A, PrmA, PrmB);
        FB = FUNC_PV2HREG1(B, PrmA, PrmB); break;
    }
   case  108:
    {
        FA = FUNC_PV2SREG1(A, PrmA, PrmB);
        FB = FUNC_PV2SREG1(B, PrmA, PrmB); break;
    }
   case  109:
    {
        FA = FUNC_PH2VREG1(A, PrmA, PrmB);
        FB = FUNC_PH2VREG1(B, PrmA, PrmB); break;
    }
   case  110:
    {
        FA = FUNC_PS2VREG1(A, PrmA, PrmB);
        FB = FUNC_PS2VREG1(B, PrmA, PrmB); break;
    }
   case  200:
    {
        FA = FUNC_PT2HREG2(PrmA, A, PrmB);
        FB = FUNC_PT2HREG2(PrmA, B, PrmB); break;
    }
   case  201:
    {
        FA = FUNC_PT2SREG2(PrmA, A, PrmB);
        FB = FUNC_PT2SREG2(PrmA, B, PrmB); break;
    }
   case  202:
    {
        FA = FUNC_PT2VREG2(PrmA, A, PrmB);
        FB = FUNC_PT2VREG2(PrmA, B, PrmB); break;
    }
   case  203:
    {
        FA = FUNC_PT2HREG2(A, PrmA, PrmB);
        FB = FUNC_PT2HREG2(B, PrmA, PrmB); break;
    }
   case  204:
    {
        FA = FUNC_PT2SREG2(A, PrmA, PrmB);
        FB = FUNC_PT2SREG2(B, PrmA, PrmB); break;
    }
    case 205:
    {
        FA = FUNC_PT2VREG2(A, PrmA, PrmB);
        FB = FUNC_PT2VREG2(B, PrmA, PrmB); break;
    }
   case  206:
    {
        FA = FUNC_PH2SREG2(A, PrmA, PrmB);
        FB = FUNC_PH2SREG2(B, PrmA, PrmB); break;
    }
   case  207:
    {
        FA = FUNC_PH2VREG2(A, PrmA, PrmB);
        FB = FUNC_PH2VREG2(B, PrmA, PrmB); break;
    }
   case  208:
    {
        FA = FUNC_PV2SREG2(A, PrmA, PrmB);
        FB = FUNC_PV2SREG2(B, PrmA, PrmB); break;
    }
   case  300:
    {
        FA = FUNC_TV2PREG3(PrmA, A, PrmB);
        FB = FUNC_TV2PREG3(PrmA, B, PrmB); break;
    }
   case  301:
    {
        FA = FUNC_TV2HREG3(PrmA, A, PrmB);
        FB = FUNC_TV2HREG3(PrmA, B, PrmB); break;
    }
   case  302:
    {
        FA = FUNC_TV2SREG3(PrmA, A, PrmB);
        FB = FUNC_TV2SREG3(PrmA, B, PrmB); break;
    }
    case 303:
    {
        FA = FUNC_TV2PREG3(A, PrmA, PrmB);
        FB = FUNC_TV2PREG3(B, PrmA, PrmB); break;
    }
    case 304:
    {
        FA = FUNC_TV2HREG3(A, PrmA, PrmB);
        FB = FUNC_TV2HREG3(B, PrmA, PrmB); break;
    }
    case 305:
    {
        FA = FUNC_TV2SREG3(A, PrmA, PrmB);
        FB = FUNC_TV2SREG3(B, PrmA, PrmB); break;
    }
    case 306:
    {
        FA = FUNC_PT2HREG3(PrmA, A, PrmB);
        FB = FUNC_PT2HREG3(PrmA, B, PrmB); break;
    }
   case  307:
    {
        FA = FUNC_PT2SREG3(PrmA, A, PrmB);
        FB = FUNC_PT2SREG3(PrmA, B, PrmB); break;
    }
    case 308:
    {
        FA = FUNC_TH2SREG3(A,PrmA,PrmB);
        FB = FUNC_TH2SREG3(B,PrmA,PrmB); break;
    }
   case  400:
    {
        FA = FUNC_TV2PREG3(PrmA, A, PrmB);
        FB = FUNC_TV2PREG3(PrmA, B, PrmB); break;
    }
   case  401:
    {
        FA = FUNC_TV2HREG4(PrmA, A, PrmB);
        FB = FUNC_TV2HREG4(PrmA, B, PrmB); break;
    }
   case  402:
    {
        FA = FUNC_TV2SREG4(PrmA, A, PrmB);
        FB = FUNC_TV2SREG4(PrmA, B, PrmB); break;
    }
   case  403:
    {
        FA = FUNC_TV2PREG4(A, PrmA, PrmB);
        FB = FUNC_TV2PREG4(B, PrmA, PrmB); break;
    }
   case  404:
    {
        FA = FUNC_TV2HREG4(A, PrmA, PrmB);
        FB = FUNC_TV2HREG4(B, PrmA, PrmB); break;
    }
   case  405:
    {
        FA = FUNC_TV2SREG4(A, PrmA, PrmB);
        FB = FUNC_TV2SREG4(B, PrmA, PrmB); break;
    }
  case   406:
    {
        FA = FUNC_PT2HREG4(PrmA, A, PrmB);
        FB = FUNC_PT2HREG4(PrmA, B, PrmB); break;
    }
   case  407:
    {
        FA = FUNC_PT2SREG4(PrmA, A, PrmB);
        FB = FUNC_PT2SREG4(PrmA, B, PrmB); break;
    }
   case  408:
    {
        FA = FUNC_TS2HREG4(A,PrmA,PrmB);
        FB = FUNC_TS2HREG4(B,PrmA,PrmB); break;
    }
   case  5600:
    {
        FA = FUNC_PS2HREG56(A,PrmA,PrmB); 
        FB = FUNC_PS2HREG56(B,PrmA,PrmB); break;
    }
    case 5601:
    {
        FA = FUNC_PV2HREG56(A,PrmA,PrmB); 
        FB = FUNC_PV2HREG56(B,PrmA,PrmB); break;
    }
   case  5602:
    {
        FA = FUNC_TX2HREG56(A,PrmA,PrmB);
        FB = FUNC_TX2HREG56(B,PrmA,PrmB); break;
    }
   case  5603:
    {
        FA = FUNC_PS2VREG56(A,PrmA,PrmB);
        FB = FUNC_PS2VREG56(B,PrmA,PrmB); break;
    }
   case  5604:
    {
        FA = FUNC_TX2SREG56(A,PrmA,PrmB);
        FB = FUNC_TX2SREG56(B,PrmA,PrmB); break;
    }
   case  5605:
    {
        FA = FUNC_TX2VREG56(A,PrmA,PrmB);
        FB = FUNC_TX2VREG56(B,PrmA,PrmB); break;
    }
    default:
    {
        FA=0;
        FB=0;
    }
  }//switch ( Nbr )

    if ( (FB * FA > 0.0) )
    {
    }
    FC = FB;    
	for(ITER=1;ITER<=ITMAX;ITER++)
    {
        if ( FB * FC > 0 )
        {
            C = A;
            FC = FA;
            D = B - A;
            E = D;
        }
        if ( fabs(FC) < fabs(FB) )
        {
            A = B;
            B = C;
            C = A;
            FA = FB;
            FB = FC;
            FC = FA;
        }
        TOL1 = 2.0 * EPS * fabs(B) + 0.5 * TOL;
        XM = 0.5 * (C - B);
        if ( (fabs(XM) <= TOL1) || (FB == 0) )
        {
            return B;
            //exit(1);
        }
        if ( (fabs(E) >= TOL1) && (fabs(FA) > fabs(FB)) )
        {
            S = FB / FA;
            if ( A == C )
            {
                P = 2.0 * XM * S;
                Q = 1 - S;
            }
            else
            {
                Q = FA / FC;
                R = FB / FC;
                P = S * (2.0 * XM * Q * (Q - R) - (B - A) * (R - 1.0));
                Q = (Q - 1.0) * (R - 1.0) * (S - 1.0);
            }
            if ( (P > 0.0) ) Q = -Q;
            P = fabs(P);
            if ( (3.0 * XM * Q - fabs(TOL1 * Q) < fabs(E * Q)) )
                AAA = 3.0 * XM * Q - fabs(TOL1 * Q);
            else
                AAA = fabs(E * Q);
            if ( 2.0 * P < AAA )
            {
                E = D;
                D = P / Q;
            }
            else
            {
                D = XM;
                E = D;
            }
        }
        else
        {
            D = XM;
            E = D;
        }
        A = B;
        FA = FB;
        if ( fabs(D) > TOL1 )
            B = B + D;
        else
        {
            if ( XM>=0 )
                ZZ=1;
            else
                ZZ=-1;
            B=B + fabs(TOL1) * ZZ;
        }


        switch ( Nbr )
		{
        case 100:
        {
            FB = FUNC_PT2HREG1(PrmA, B, PrmB);break;
        }
        case 101:
        {
            FB = FUNC_PT2SREG1(PrmA, B, PrmB);break;
        }
        case 102:
        {
            FB = FUNC_PT2VREG1(PrmA, B, PrmB);break;
        }
       case  103:
        {
            FB = FUNC_PT2HREG1(B, PrmA, PrmB);break;
        }
        case 104:
        {
            FB = FUNC_PT2SREG1(B, PrmA, PrmB);break;
        }
        case 105:
        {
           FB = FUNC_PT2VREG1(B, PrmA, PrmB);break;
        }
        case 106:
        {
            FB = FUNC_PS2HREG1(B, PrmA, PrmB);break;
        }
        case 107:
        {
            FB = FUNC_PV2HREG1(B, PrmA, PrmB);break;
        }
        case 108:
        {
            FB = FUNC_PV2SREG1(B, PrmA, PrmB);break;
        }
        case 109:
        {
            FB = FUNC_PH2VREG1(B, PrmA, PrmB);break;
        }
        case 110:
        {
            FB = FUNC_PS2VREG1(B, PrmA, PrmB);break;
        }
        case 200:
        {
            FB = FUNC_PT2HREG2(PrmA, B, PrmB);break;
        }
        case 201:
        {
            FB = FUNC_PT2SREG2(PrmA, B, PrmB);break;
        }
        case 202:
        {
            FB = FUNC_PT2VREG2(PrmA, B, PrmB);break;
        }
        case 203:
        {
            FB = FUNC_PT2HREG2(B, PrmA, PrmB);break;
        }
        case 204:
        {
            FB = FUNC_PT2SREG2(B, PrmA, PrmB);break;
        }
        case 205:
        {
           FB = FUNC_PT2VREG2(B, PrmA, PrmB);break;
        }
        case 206:
        {
            FB = FUNC_PH2SREG2(B, PrmA, PrmB);break;
        }
        case 207:
        {
            FB = FUNC_PH2VREG2(B, PrmA, PrmB);break;
        }
        case 208:
        {
            FB = FUNC_PV2SREG2(B, PrmA, PrmB);break;
        }
        case 300:
        {
            FB = FUNC_TV2PREG3(PrmA, B, PrmB);break;
        }
        case 301:
        {
            FB = FUNC_TV2HREG3(PrmA, B, PrmB);break;
        }
        case 302:
        {
            FB = FUNC_TV2SREG3(PrmA, B, PrmB);break;
        }
        case 303:
        {
            FB = FUNC_TV2PREG3(B, PrmA, PrmB);break;
        }
        case 304:
        {
            FB = FUNC_TV2HREG3(B, PrmA, PrmB);break;
        }
        case 305:
        {
            FB = FUNC_TV2SREG3(B, PrmA, PrmB);break;
        }
        case 306:
        {
            FB = FUNC_PT2HREG3(PrmA, B, PrmB);break;
        }
       case  307:
        {
            FB = FUNC_PT2SREG3(PrmA, B, PrmB);break;
        }
       case  308:
        {
            FB = FUNC_TH2SREG3(B,PrmA,PrmB);break;
        }
       case  400:
        {
            FB = FUNC_TV2PREG4(PrmA, B, PrmB);break;
        }
       case  401:
        {
            FB = FUNC_TV2HREG4(PrmA, B, PrmB);break;
        }
        case 402:
        {
            FB = FUNC_TV2SREG4(PrmA, B, PrmB);break;
        }
       case  403:
        {
            FB = FUNC_TV2PREG4(B, PrmA, PrmB);break;
        }
       case  404:
        {
            FB = FUNC_TV2HREG4(B, PrmA, PrmB);break;
        }
        case 405:
        {
            FB = FUNC_TV2SREG4(B, PrmA, PrmB);break;
        }
        case 406:
        {
            FB = FUNC_PT2HREG4(PrmA, B, PrmB);break;
        }
        case 407:
        {
            FB = FUNC_PT2SREG4(PrmA, B, PrmB);break;
        }
       case  408:
        {
            FB = FUNC_TS2HREG4(B, PrmA, PrmB);break;
        }
       case  5600:
        {
            FB = FUNC_PS2HREG56(B,PrmA,PrmB);break;
        }
       case  5601:
        {
            FB = FUNC_PV2HREG56(B,PrmA,PrmB);break;
        }
       case  5602:
        {
            FB = FUNC_TX2HREG56(B,PrmA,PrmB);break;
        }
       case  5603:
        {
            FB = FUNC_PS2VREG56(B,PrmA,PrmB);break;
        }
       case  5604:
        {
            FB = FUNC_TX2SREG56(B,PrmA,PrmB);break;
        }
       case  5605:
        {
            FB = FUNC_TX2VREG56(B,PrmA,PrmB); break;
        }
	   default:
        {
        }
        } //switch ( Nbr )

    }
    return B;
}
//****************************************8888  20:55



//**********************************************P
void  P2T67(double P, double* T, int *RANGE)
{
    SUBRANGEBYP(P,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
		{ *T=P2T_KK(P); break;}
    default:
		{*T=-1.0; break;}
    }
}

void  P2HL67(double P, double* H,int *RANGE)
{
    P2HL_KK(P,H,RANGE);
}

void  P2HG67(double P, double* H,int *RANGE)
{
    P2HG_KK(P,H,RANGE);
}

void  P2SL67(double P, double* S,int *RANGE)
{
    P2SL_KK(P,S,RANGE);
}

void  P2SG67(double P, double* S,int *RANGE)
{
    P2SG_KK(P,S,RANGE);
}

void  P2VL67(double P, double* V,int *RANGE)
{
    P2VL_KK(P,V,RANGE);
}

void  P2VG67(double P, double* V,int *RANGE)
{
    P2VG_KK(P,V,RANGE);
}

void  P2L67(double P, double* T,double*H,double*S,double*V,double*X,int *RANGE)
{
    P2L_KK(P,T,H,S,V,X,RANGE);
}

void  P2G67(double P, double* T,double*H,double*S,double*V,double*X,int *RANGE)
{
    P2G_KK(P,T,H,S,V,X,RANGE);
}


void  P2CPL67(double P ,double *CP,int *RANGE)
{
    P2CPL(P,CP,RANGE);
}

void  P2CPG67(double P ,double *CP,int *RANGE)
{
    P2CPG(P,CP,RANGE);
}

void  P2CVL67(double P ,double *CV,int *RANGE)
{
    P2CVL(P,CV,RANGE);
}

void  P2CVG67(double P ,double *CV,int *RANGE)
{
    P2CVG(P,CV,RANGE);
}

void  P2EL67(double P  ,double * E,int *RANGE)
{
    P2EL_KK(P,E,RANGE);
}

void  P2EG67(double P  ,double * E,int *RANGE)
{
    P2EG_KK(P,E,RANGE);
}

void  P2SSPL67(double P,double *SSP,int *RANGE)
{
    P2SSPL(P,SSP,RANGE);
}

void  P2SSPG67(double P,double *SSP,int *RANGE)
{
    P2SSPG(P,SSP,RANGE);
}

void  P2KSL67(double P ,double *KS ,int *RANGE)
{
    P2KSL(P,KS,RANGE);
}

void  P2KSG67(double P ,double *KS ,int *RANGE)
{
    P2KSG(P,KS,RANGE);
}

void  P2ETAL67(double P,double *ETA,int *RANGE)
{
    P2ETAL(P,ETA,RANGE);
}

void  P2ETAG67(double P,double *ETA,int *RANGE)
{
    P2ETAG(P,ETA,RANGE);
}


void  P2UL67(double P  ,double *U,int *RANGE)
{
    P2UL(P,U,RANGE);
}

void  P2UG67(double P  ,double *U,int *RANGE)
{
    P2UG(P,U,RANGE);
}


void  P2RAMDL67(double P ,double *RAMD,int *RANGE)
{
    P2RAMDL(P,RAMD,RANGE);
}

void  P2RAMDG67(double P ,double *RAMD,int *RANGE)
{
    P2RAMDG(P,RAMD,RANGE);
}

void  P2PRNL67(double P  ,double *PRN, int *RANGE)
{
    P2PRNL(P,PRN,RANGE);
}

void  P2PRNG67(double P  ,double *PRN, int *RANGE)
{
    P2PRNG(P,PRN,RANGE);
}

void  P2EPSL67(double P  ,double *EPS,int *RANGE)
{
    P2EPSL(P,EPS,RANGE);
}

void  P2EPSG67(double P  ,double *EPS,int *RANGE)
{
    P2EPSG(P,EPS,RANGE);
}

void  P2NL67(double P, double Lamd,double *N,int *RANGE)
{
    P2NL(P,Lamd,N,RANGE);
}

void  P2NG67(double P, double Lamd,double *N,int *RANGE)
{
    P2NG(P,Lamd,N,RANGE);
}


void  PT2H67(double P, double T,double *H,int *RANGE)
{
    PT2H_KK(P,T,H,RANGE);
}

void  PT2S67(double P, double T,double *S,int *RANGE)
{
    PT2S_KK(P,T,S,RANGE);
}

void  PT2V67(double P, double T,double *V,int *RANGE)
{
    PT2V_KK(P,T,V,RANGE);
}

void  PT2X67(double P, double T,double *X,int *RANGE)
{
    PT2X(P,T,X,RANGE);
}

void  PT67(double P, double T, double * H, double *S, double *V, double *X,int *RANGE)
{
    PT_KK(P,T,H,S,V,X,RANGE);
}

void  PT2CP67(double P, double T,double *CP,int *RANGE)
{
    PT2CP_KK(P,T,CP,RANGE);
}

void  PT2CV67(double P, double T,double *CV,int *RANGE)
{
    PT2CV_KK(P,T,CV,RANGE);
}

void  PT2E67(double P , double T,double *E,int *RANGE)
{
    PT2E_KK(P,T,E,RANGE);
}

void  PT2SSP67(double P, double T,double *SSP,int *RANGE)
{
    PT2SSP_KK(P,T,SSP,RANGE);
}

void  PT2KS67(double P , double T,double *KS,int *RANGE)
{
    PT2KS_KK(P,T,KS,RANGE);
}

void  PT2ETA67(double P, double T,double *ETA,int *RANGE)
{
    PT2ETA(P,T,ETA,RANGE);
}


void  PT2U67(double P, double T,double *U,int *RANGE)
{
    PT2U_KK(P,T,U,RANGE);
}

void  PT2RAMD67(double P,double T,double *RAMD,int *RANGE)
{
    PT2RAMD_KK(P,T,RAMD,RANGE);
}

void  PT2PRN67(double P, double T,double *PRN,int *RANGE)
{
    PT2PRN_KK(P,T,PRN,RANGE);
}


void  PT2EPS67(double P, double T,double *EPS,int *RANGE)
{
    PT2EPS_KK(P,T,EPS,RANGE);
}


void  PT2N67(double P, double T,double LAMD, double *N,int *RANGE)
{
    PT2N_KK(P,T,LAMD,N,RANGE);
}


void  PH2T67(double P, double H,double *T,int *RANGE)
{
    PH2T_KK(P,H,T,RANGE);
}

void PH2S67(double P, double H,double *S,int *RANGE)
{
    PH2S_KK(P,H,S,RANGE);
}

void  PH2V67(double P, double H,double *V,int *RANGE)
{
    PH2V_KK(P,H,V,RANGE);
}

void  PH2X67(double P, double H,double *X,int *RANGE)
{
    PH2X_KK(P,H,X,RANGE);
}
void  PH2SSP67(double P, double H,double *SSP,int *RANGE)
{
    PH2SSP_KK(P,H,SSP,RANGE);
}
void  PH67(double P  ,double *T,double H, double * S,double *V,double *X,int *RANGE)
{
    PH_KK(P,T,H,S,V,X,RANGE);
}

void  PS2T67(double P, double S,double *T,int *RANGE)
{
    PS2T_KK(P,S,T,RANGE);
}

void  PS2H67(double P, double S,double *H,int *RANGE)
{
    PS2H_KK(P,S,H,RANGE);
}

void  PS2V67(double P, double S,double *V,int *RANGE)
{
    PS2V_KK(P,S,V,RANGE);
}

void  PS2X67(double P, double S,double *X,int *RANGE)
{
    PS2X_KK(P,S,X,RANGE);
}


void  PS67(double P  , double *T, double *H,double S,double *V, double *X,int *RANGE)
{
    PS_KK(P,T,H,S,V,X,RANGE);
}

void  PV2T67(double P, double V,double *T,int *RANGE)
{
    PV2T_KK(P,V,T,RANGE);
}

void  PV2H67(double P, double V,double *H,int *RANGE)
{
    PV2H_KK(P,V,H,RANGE);
}

void  PV2S67(double P, double V,double *S,int *RANGE)
{
    PV2S_KK(P,V,S,RANGE);
}

void  PV2X67(double P, double V,double *X,int *RANGE)
{
    PV2X_KK(P,V,X,RANGE);
}


void  PV67(double P  ,double *T, double *H, double *S,double V,double *X,int *RANGE)
{
    PV_KK(P,T,H,S,V,X,RANGE);
}

void  PX2T67(double P, double X,double *T,int *RANGE)
{
    PX2T_KK(P,X,T,RANGE);
}

void  PX2H67(double P, double X,double *H,int *RANGE)
{
    PX2H_KK(P,X,H,RANGE);
}

void  PX2S67(double P, double X,double *S,int *RANGE)
{
    PX2S_KK(P,X,S,RANGE);
}

void  PX2V67(double P, double X,double *V,int *RANGE)
{
    PX2V_KK(P,X,V,RANGE);
}


void  PX67(double P  ,double *T, double *H, double *S, double *V,double X,int *RANGE)
{
    PX_KK(P,T,H,S,V,X,RANGE);
}

void  T2P67(double T ,double *P, int *RANGE)
{
    SUBRANGEBYT(T,RANGE);
    switch ( *RANGE )
    { 
	case 5: 
	case 6:
        { *P=T2P_KK(T); break;}
    default:
        { *P=-1.0; break;}
    }
}


void  T2HL67(double T,double *H, int *RANGE)
{
    T2HL(T,H,RANGE);
}

void  T2HG67(double T,double *H, int *RANGE)
{
    T2HG(T,H,RANGE);
}

void  T2SL67(double T,double *S,  int *RANGE)
{
    T2SL(T,S,RANGE);
}

void  T2SG67(double T,double *S,  int *RANGE)
{
    T2SG(T,S,RANGE);
}

void  T2VL67(double T,double *V, int *RANGE)
{
    T2VL(T,V,RANGE);
}

void  T2VG67(double T,double *V, int *RANGE)
{
    T2VG(T,V,RANGE);
}

void  T2L67(double *P,double T  ,double *H,double *S,double *V,double *X,int *RANGE)
{
    T2L(P,T,H,S,V,X,RANGE);
}

void  T2G67(double *P,double T  ,double *H,double *S,double *V,double *X,int *RANGE)
{
    T2G(P,T,H,S,V,X,RANGE);
}

void  T2CPL67(double T,double *CP,int *RANGE)
{
    T2CPL(T,CP,RANGE);
}

void  T2CPG67(double T,double *CP,int *RANGE)
{
    T2CPG(T,CP,RANGE);
}

void  T2CVL67(double T,double *CV,int *RANGE)
{
    T2CVL(T,CV,RANGE);
}

void  T2CVG67(double T,double *CV,int *RANGE)
{
    T2CVG(T,CV,RANGE);
}

void  T2EL67(double T ,double *E,int *RANGE)
{
    T2EL(T,E,RANGE);
}

void  T2EG67(double T ,double *E,int *RANGE)
{
    T2EG(T,E,RANGE);
}

void  T2SSPL67(double T ,double *SSP,int *RANGE)
{
    T2SSPL(T,SSP,RANGE);
}

void  T2SSPG67(double T ,double *SSP,int *RANGE)
{
    T2SSPG(T,SSP,RANGE);
}

void  T2KSL67(double T  ,double *KS,int *RANGE)
{
    T2KSL(T,KS,RANGE);
}

void  T2KSG67(double T  ,double *KS,int *RANGE)
{
    T2KSG(T,KS,RANGE);
}

void  T2ETAL67(double T ,double *ETA,int *RANGE)
{
    T2ETAL(T,ETA,RANGE);
}

void  T2ETAG67(double T ,double *ETA,int *RANGE)
{
    T2ETAG(T,ETA,RANGE);
}

void  T2UL67(double T  ,double *U,int *RANGE)
{
    T2UL(T,U,RANGE);
}

void  T2UG67(double T  ,double *U,int *RANGE)
{
    T2UG(T,U,RANGE);
}

void  T2RAMDL67(double T ,double *RAMD,int *RANGE)
{
    T2RAMDL(T,RAMD,RANGE);
}

void  T2RAMDG67(double T ,double *RAMD,int *RANGE)
{
    T2RAMDG(T,RAMD,RANGE);
}

void  T2PRNL67(double T  ,double *PRN,int *RANGE)
{
    T2PRNL(T,PRN,RANGE);
}

void  T2PRNG67(double T  ,double *PRN,int *RANGE)
{
    T2PRNG(T,PRN,RANGE);
}

void  T2EPSL67(double T  ,double *EPS,int *RANGE)
{
    T2EPSL(T,EPS,RANGE);
}

void  T2EPSG67(double T  ,double *EPS,int *RANGE)
{
    T2EPSG(T,EPS,RANGE);
}

void  T2NL67(double T, double Lamd,double *N,int *RANGE)
{
    T2NL(T,Lamd,N,RANGE);
}

void  T2NG67(double T, double Lamd,double *N,int *RANGE)
{
    T2NG(T,Lamd,N,RANGE);
}

void  T2SURFT67(double T ,double *SURFT,int *RANGE)
{
    T2SurfT(T,SURFT,RANGE);
}

void  TH2P67(double T, double H,double *P,int *RANGE)
{
    TH2P_KK(T,H,P,RANGE);
}

void  TH2PLP67(double T, double H,double *P,int *RANGE)
{
    TH2PLP_KK(T,H,P,RANGE);
}


void  TH2PHP67(double T, double H,double *P,int *RANGE)
{
    TH2PHP_KK(T,H,P,RANGE);
}


void  TH2S67(double T, double H,double *S,int *RANGE)
{
    TH2S_KK(T,H,S,RANGE);
}

void  TH2SLP67(double T, double H,double *S,int *RANGE)
{
    TH2SLP_KK(T,H,S,RANGE);
}

void  TH2SHP67(double T, double H,double *S,int *RANGE)
{
    TH2SHP_KK(T,H,S,RANGE);
}

void  TH2V67(double T, double H,double *V,int *RANGE)
{
    TH2V_KK(T,H,V,RANGE);
}

void  TH2VLP67(double T, double H,double *V,int *RANGE)
{
    TH2VLP_KK(T,H,V,RANGE);
}

void  TH2VHP67(double T, double H,double *V,int *RANGE)
{
    TH2VHP_KK(T,H,V,RANGE);
}

void  TH2X67(double T, double H,double *X,int *RANGE)
{
    TH2X_KK(T,H,X,RANGE);
}

void  TH2XLP67(double T, double H,double *X,int *RANGE)
{
    TH2XLP_KK(T,H,X,RANGE);
}

void  TH2XHP67(double T, double H,double *X,int *RANGE)
{
    TH2XHP_KK(T,H,X,RANGE);
}


void  TH67(double *P,double T, double H,double * S,double *V,double *X,int *RANGE)
{
    TH_KK(P,T,H,S,V,X,RANGE);
}

void  THLP67(double *P,double T, double H,double * S,double *V,double *X,int *RANGE)
{
    THLP_KK(P,T,H,S,V,X,RANGE);
}

void  THHP67(double *P,double T, double H,double * S,double *V,double *X,int *RANGE)
{
    THHP_KK(P,T,H,S,V,X,RANGE);
}

void  TS2PHP67(double T, double S,double *P,int *RANGE)
{
    TS2PHP_KK(T,S,P,RANGE);
}


void  TS2PLP67(double T, double S,double *P,int *RANGE)
{
    TS2PLP_KK(T,S,P,RANGE);
}

void  TS2P67(double T, double S,double *P,int *RANGE)
{
    TS2P_KK(T,S,P,RANGE);
}


void  TS2HHP67(double T, double S,double *H,int *RANGE)
{
    TS2HHP_KK(T,S,H,RANGE);
}


void  TS2HLP67(double T, double S,double *H,int *RANGE)
{
    TS2HLP_KK(T,S,H,RANGE);
}

void  TS2H67(double T, double S,double *H,int *RANGE)
{
    TS2H_KK(T,S,H,RANGE);
}

void  TS2VHP67(double T, double S,double *V,int *RANGE)
{
    TS2VHP_KK(T,S,V,RANGE);
}

void  TS2VLP67(double T, double S,double *V,int *RANGE)
{
    TS2VLP_KK(T,S,V,RANGE);
}

void  TS2V67(double T, double S,double *V,int *RANGE)
{
    TS2V_KK(T,S,V,RANGE);
}

void  TS2X67(double T, double S,double *X,int *RANGE)
{
    TS2X_KK(T,S,X,RANGE);
}


void  TSHP67(double *P,double T  ,double *H,double S,double *V, double *X,int *RANGE)
{
    TSHP_KK(P,T,H,S,V,X,RANGE);
}


void  TSLP67(double *P,double T  ,double *H,double S,double *V, double *X,int *RANGE)
{
    TSLP_KK(P,T,H,S,V,X,RANGE);
}

void  TS67(double *P,double T  ,double *H,double S,double *V, double *X,int *RANGE)
{
    TS_KK(P,T,H,S,V,X,RANGE);
}


void  TV2P67(double T, double V,double *P,int *RANGE)
{
    TV2P_KK(T,V,P,RANGE);
}

void  TV2H67(double T, double V,double *H,int *RANGE)
{
    TV2H_KK(T,V,H,RANGE);
}

void  TV2S67(double T, double V,double *S,int *RANGE)
{
    TV2S_KK(T,V,S,RANGE);
}

void  TV2X67(double T, double V,double *X,int *RANGE)
{
    TV2X_KK(T,V,X,RANGE);
}

void  TV67(double *P,double T  ,double *H,double *S ,double V , double *X,int *RANGE)
{
    TV_KK(P,T,H,S,V,X,RANGE);
}

void  TX2P67(double T, double X,double *P,int *RANGE)
{
    TX2P_KK(T,X,P,RANGE);
}

void  TX2H67(double T, double X,double *H,int *RANGE)
{
    TX2H_KK(T,X,H,RANGE);
}

void  TX2S67(double T, double X,double *S,int *RANGE)
{
    TX2S_KK(T,X,S,RANGE);
}

void  TX2V67(double T, double X,double *V,int *RANGE)
{
    TX2V_KK(T,X,V,RANGE);
}

void  TX67(double *P,double T  ,double *H, double *S, double *V,double X,int *RANGE)
{
    TX_KK(P,T,H,S,V,X,RANGE);
}

void  H2TL67(double H, double *T,int *RANGE)
{
    H2TL_KK(H,T,RANGE);
}

void  HS2P67(double H, double S,double *P,int *RANGE)
{
    HS2P_KK(H,S,P,RANGE);
}

void  HS2T67(double H, double S,double *T,int *RANGE)
{
    HS2T_KK(H,S,T,RANGE);
}

void  HS2V67(double H, double S,double *V,int *RANGE)
{
    HS2V_KK(H,S,V,RANGE);
}

void  HS2X67(double H, double S,double *X,int *RANGE)
{
    HS2X_KK(H,S,X,RANGE);
}

void  HS67(double *P, double *T,double H, double S,double *V, double *X,int *RANGE)
{
    HS_KK(P,T,H,S,V,X,RANGE);
}

void  HV2P67(double H, double V,double *P,int *RANGE)
{
    HV2P_KK(H,V,P,RANGE);
}

void  HV2T67(double H, double V,double *T,int *RANGE)
{
    HV2T_KK(H,V,T,RANGE);
}

void  HV2S67(double H, double V,double *S,int *RANGE)
{
    HV2S_KK(H,V,S,RANGE);
}

void  HV2X67(double H, double V,double *X,int *RANGE)
{
    HV2X_KK(H,V,X,RANGE);
}

void  HV67(double *P, double *T,double H,double *S, double V,double *X,int *RANGE)
{
    HV_KK(P,T,H,S,V,X,RANGE);
}

void  HX2P67(double H, double X,double *P,int *RANGE)
{
    HX2P_KK(H,X,P,RANGE);
}

void  HX2PLP67(double H, double X,double *P,int *RANGE)
{
    HX2PLP_KK(H,X,P,RANGE);
}

void  HX2PHP67(double H, double X,double *P,int *RANGE)
{
    HX2PHP_KK(H,X,P,RANGE);
}

void  HX2T67(double H, double X,double *T,int *RANGE)
{
    HX2T_KK(H,X,T,RANGE);
}

void  HX2TLP67(double H, double X,double *T,int *RANGE)
{
    HX2TLP_KK(H,X,T,RANGE);
}

void  HX2THP67(double H, double X,double *T,int *RANGE)
{
    HX2THP_KK(H,X,T,RANGE);
}

void  HX2S67(double H, double X,double *S,int *RANGE)
{
    HX2S_KK(H,X,S,RANGE);
}

void  HX2SLP67(double H, double X,double *S,int *RANGE)
{
    HX2SLP_KK(H,X,S,RANGE);
}

void  HX2SHP67(double H, double X,double *S,int *RANGE)
{
    HX2SHP_KK(H,X,S,RANGE);
}

void  HX2V67(double H, double X,double *V,int *RANGE)
{
    HX2V_KK(H,X,V,RANGE);
}

void  HX2VLP67(double H, double X,double *V,int *RANGE)
{
    HX2VLP_KK(H,X,V,RANGE);
}

void  HX2VHP67(double H, double X,double *V,int *RANGE)
{
    HX2VHP_KK(H,X,V,RANGE);
}

void  HX67(double *P, double *T,double H,double *S, double *V,double X,int *RANGE)
{
    HX_KK(P,T,H,S,V,X,RANGE);
}


void  HXLP67(double *P, double *T,double H,double *S, double *V,double X,int *RANGE)
{
    HXLP_KK(P,T,H,S,V,X,RANGE);
}

void  HXHP67(double *P, double *T,double H,double *S, double *V,double X,int *RANGE)
{
    HXHP_KK(P,T,H,S,V,X,RANGE);
}

void  S2TG67(double S,double *T,int *RANGE)
{
    S2TG_KK(S,T,RANGE);
}


void  SV2P67(double S, double V,double *P,int *RANGE)
{
    SV2P_KK(S,V,P,RANGE);
}

void  SV2T67(double S, double V,double *T,int *RANGE)
{
    SV2T_KK(S,V,T,RANGE);
}

void  SV2H67(double S, double V,double *H,int *RANGE)
{
    SV2H_KK(S,V,H,RANGE);
}

void  SV2X67(double S, double V,double *X,int *RANGE)
{
    SV2X_KK(S,V,X,RANGE);
}


void  SV67(double *P, double *T, double *H,double S, double V,double *X,int *RANGE)
{
    SV_KK(P,T,H,S,V,X,RANGE);
}

void  SX2P67(double S, double X,double *P,int *RANGE)
{
    SX2P_KK(S,X,P,RANGE);
}

void  SX2PLP67(double S, double X,double *P,int *RANGE)
{
    SX2PLP_KK(S,X,P,RANGE);
}

void  SX2PMP67(double S, double X,double *P,int *RANGE)
{
    SX2PMP_KK(S,X,P,RANGE);
}

void  SX2PHP67(double S, double X,double *P,int *RANGE)
{
    SX2PHP_KK(S,X,P,RANGE);
}

void  SX2TLP67(double S, double X,double *T,int *RANGE)
{
    SX2TLP_KK(S,X,T,RANGE);
}

void  SX2TMP67(double S, double X,double *T,int *RANGE)
{
    SX2TMP_KK(S,X,T,RANGE);
}

void  SX2THP67(double S, double X,double *T,int *RANGE)
{
    SX2THP_KK(S,X,T,RANGE);
}

void  SX2T67(double S, double X,double *T,int *RANGE)
{
    SX2T_KK(S,X,T,RANGE);
}

void  SX2H67(double S, double X,double *H,int *RANGE)
{
    SX2H_KK(S,X,H,RANGE);
}

void  SX2HLP67(double S, double X,double *H,int *RANGE)
{
    SX2HLP_KK(S,X,H,RANGE);
}

void  SX2HMP67(double S, double X,double *H,int *RANGE)
{
    SX2HMP_KK(S,X,H,RANGE);
}

void  SX2HHP67(double S, double X,double *H,int *RANGE)
{
    SX2HHP_KK(S,X,H,RANGE);
}

void  SX2V67(double S, double X,double *V,int *RANGE)
{
    SX2V_KK(S,X,V,RANGE);
}

void  SX2VLP67(double S, double X,double *V,int *RANGE)
{
    SX2VLP_KK(S,X,V,RANGE);
}

void  SX2VMP67(double S, double X,double *V,int *RANGE)
{
    SX2VMP_KK(S,X,V,RANGE);
}

void  SX2VHP67(double S, double X,double *V,int *RANGE)
{
    SX2VHP_KK(S,X,V,RANGE);
}

void  SX67(double *P, double *T, double *H,double S,double *V,double X,int *RANGE)
{
    SX_KK(P,T,H,S,V,X,RANGE);
}

void  SXLP67(double *P, double *T, double *H,double S,double *V,double X,int *RANGE)
{
    SXLP_KK(P,T,H,S,V,X,RANGE);
}


void  SXMP67(double *P, double *T, double *H,double S,double *V,double X,int *RANGE)
{
    SXMP_KK(P,T,H,S,V,X,RANGE);
}


void  SXHP67(double *P, double *T, double *H,double S,double *V,double X,int *RANGE)
{
    SXHP_KK(P,T,H,S,V,X,RANGE);
}

void  V2TG67(double V , double *T,int *RANGE)
{
    V2TG_KK(V,T,RANGE);
}

void  VX2P67(double V, double X,double *P,int *RANGE)
{
    VX2P_KK(V,X,P,RANGE);
}
void  VX2PLP67(double V, double X,double *P,int *RANGE)
{
    VX2PLP_KK(V,X,P,RANGE);
}

void  VX2PHP67(double V, double X,double *P,int *RANGE)
{
    VX2PHP_KK(V,X,P,RANGE);
}

void  VX2T67(double V, double X,double *T,int *RANGE)
{
    VX2T_KK(V,X,T,RANGE);
}

void  VX2TLP67(double V, double X,double *T,int *RANGE)
{
    VX2TLP_KK(V,X,T,RANGE);
}

void  VX2THP67(double V, double X,double *T,int *RANGE)
{
    VX2THP_KK(V,X,T,RANGE);
}

void  VX2H67(double V, double X,double *H,int *RANGE)
{
    VX2H_KK(V,X,H,RANGE);
}

void  VX2HLP67(double V, double X,double *H,int *RANGE)
{
    VX2HLP_KK(V,X,H,RANGE);
}

void  VX2HHP67(double V, double X,double *H,int *RANGE)
{
    VX2HHP_KK(V,X,H,RANGE);
}

void  VX2S67(double V, double X,double *S,int *RANGE)
{
    VX2S_KK(V,X,S,RANGE);
}

void  VX2SLP67(double V, double X,double *S,int *RANGE)
{
    VX2SLP_KK(V,X,S,RANGE);
}

void  VX2SHP67(double V, double X,double *S,int *RANGE)
{
    VX2SHP_KK(V,X,S,RANGE);
}

void  VX67(double *P, double *T, double *H, double *S,double V, double X,int *RANGE)
{
    VX_KK(P,T,H,S,V,X,RANGE);
}

void  VXLP67(double *P, double *T, double *H, double *S,double V, double X,int *RANGE)
{
    VXLP_KK(P,T,H,S,V,X,RANGE);
}

void  VXHP67(double *P, double *T, double *H, double *S,double V, double X,int *RANGE)
{
    VXHP_KK(P,T,H,S,V,X,RANGE);
}





