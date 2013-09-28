//UEwasp.h

//************************************************************
//     作    者：                                            *
//              长沙优易软件开发有限公司(UESoft Corp.) 邝祝芳*
//     文件名称：                                            * 
//                DLL外部接口(水和蒸汽性质计算动态连接库)    *
//     完成时间：                                            *
//                2005年5月                                  *
//************************************************************ 

extern "C" void _stdcall  SETSTD_WASP(int  STDID);
extern "C" void _stdcall  GETSTD_WASP(int  * STDID);
extern "C" void _stdcall  P2T(double P, double * T, int * RANGE );
extern "C" void _stdcall  P2HL(double P, double * H , int * RANGE);
extern "C" void _stdcall  P2HG(double P  ,double * H , int * RANGE);
extern "C" void _stdcall  P2SL(double P  ,double * S , int * RANGE);
extern "C" void _stdcall  P2SG(double P  ,double * S , int * RANGE);
extern "C" void _stdcall  P2VL(double P  ,double * V , int * RANGE);
extern "C" void _stdcall  P2VG(double P  ,double * V , int * RANGE);
extern "C" void _stdcall  P2L(double P   ,double * T,double * H, double *S, double *V, double *X, int * RANGE);
extern "C" void _stdcall  P2G(double P   ,double * T,double * H, double *S, double *V, double *X, int * RANGE);
extern "C" void _stdcall  P2CPL(double P ,double *CP,int * RANGE);
extern "C" void _stdcall  P2CPG(double P ,double *CP,int * RANGE);
extern "C" void _stdcall  P2CVL(double P ,double *CV,int * RANGE);
extern "C" void _stdcall  P2CVG(double P ,double *CV,int * RANGE);
extern "C" void _stdcall  P2EL(double P  ,double * E,int    * RANGE);
extern "C" void _stdcall  P2EG(double P  ,double * E,int    * RANGE);
extern "C" void _stdcall  P2SSPL(double P,double *SSP,int   * RANGE);
extern "C" void _stdcall  P2SSPG(double P,double *SSP,int   * RANGE);
extern "C" void _stdcall  P2KSL(double P ,double *KS ,int   * RANGE);
extern "C" void _stdcall  P2KSG(double P ,double *KS ,int   * RANGE);
extern "C" void _stdcall  P2ETAL(double P,double *ETA,int * RANGE);
extern "C" void _stdcall  P2ETAG(double P,double *ETA,int * RANGE);
extern "C" void _stdcall  P2UL(double P  ,double *U,int * RANGE);
extern "C" void _stdcall  P2UG(double P  ,double *U,int * RANGE);
extern "C" void _stdcall  P2RAMDL(double P ,double *RAMD, int * RANGE);
extern "C" void _stdcall  P2RAMDG(double P ,double *RAMD, int * RANGE);
extern "C" void _stdcall  P2PRNL(double P  ,double *PRN,  int * RANGE);
extern "C" void _stdcall  P2PRNG(double P  ,double *PRN,  int * RANGE);
extern "C" void _stdcall  P2EPSL(double P  ,double *EPS,  int * RANGE);
extern "C" void _stdcall  P2EPSG(double P  ,double *EPS,  int * RANGE);
extern "C" void _stdcall  P2NL(double P, double Lamd,double *N,int * RANGE);
extern "C" void _stdcall  P2NG(double P, double Lamd,double *N,int * RANGE);

extern "C" void _stdcall  PT2H(double P, double T,double *H,int * RANGE);
extern "C" void _stdcall  PT2S(double P, double T,double *S,int * RANGE);
extern "C" void _stdcall  PT2V(double P, double T,double *V,int * RANGE);
extern "C" void _stdcall  PT2X(double P, double T,double *X,int * RANGE);
extern "C" void _stdcall  PT(double P, double T, double * H, double *S, double *V, double *X, int * RANGE);
extern "C" void _stdcall  PT2MV(double P,double T,double *H,double *S,double *V,double *E,double *CP,double *CV,double *SSP,int *Range);
extern "C" void _stdcall  PT2CP(double P, double T,double *CP,int * RANGE);
extern "C" void _stdcall  PT2CV(double P, double T,double *CV,int * RANGE);
extern "C" void _stdcall  PT2E(double P, double T,double *E,int * RANGE);
extern "C" void _stdcall  PT2SSP(double P, double T,double *SSP,int * RANGE);
extern "C" void _stdcall  PT2KS(double P, double T,double *KS,int * RANGE);
extern "C" void _stdcall  PT2ETA(double P, double T,double *ETA,int * RANGE);
extern "C" void _stdcall  PT2U(double P, double T,double *U,int * RANGE);
extern "C" void _stdcall  PT2RAMD(double P, double T,double *RAMD,int * RANGE);
extern "C" void _stdcall  PT2PRN(double P, double T,double *PRN,int * RANGE);
extern "C" void _stdcall  PT2EPS(double P, double T,double *EPS,int * RANGE);
extern "C" void _stdcall  PT2N(double P, double T,double LAMD, double *N,int * RANGE);

extern "C" void _stdcall  PT2HMV(double P, double T,double *H);
extern "C" void _stdcall  PT2VMV(double P, double T,double *V);
extern "C" void _stdcall  PT2SMV(double P, double T,double *S);
extern "C" void _stdcall  PT2EMV(double P, double T,double *E);
extern "C" void _stdcall  PT2CPMV(double P, double T,double *CP);
extern "C" void _stdcall  PT2CVMV(double P, double T,double *CV);
extern "C" void _stdcall  PT2SSPMV(double P, double T,double *SSP);

extern "C" void _stdcall  PH2T(double P, double H,double *T,int * RANGE);
extern "C" void _stdcall  PH2S(double P, double H,double *S,int * RANGE);
extern "C" void _stdcall  PH2V(double P, double H,double *V,int * RANGE);
extern "C" void _stdcall  PH2X(double P, double H,double *X,int * RANGE);
extern "C" void _stdcall  PH2SSP(double P, double H,double *ssp,int * RANGE);
extern "C" void _stdcall  PH(double P  ,double *T,double H, double * S,double *V,double *X,int * RANGE);

extern "C" void _stdcall  PS2T(double P, double S,double *T,int * RANGE);
extern "C" void _stdcall  PS2H(double P, double S,double *H,int * RANGE);
extern "C" void _stdcall  PS2V(double P, double S,double *V,int * RANGE);
extern "C" void _stdcall  PS2X(double P, double S,double *X,int * RANGE);
extern "C" void _stdcall  PS(double P  ,double *T, double *H,double S,double *V, double *X,int * RANGE);

extern "C" void _stdcall  PV2T(double P, double V,double *T,int * RANGE);
extern "C" void _stdcall  PV2H(double P, double V,double *H,int * RANGE);
extern "C" void _stdcall  PV2S(double P, double V,double *S,int * RANGE);
extern "C" void _stdcall  PV2X(double P, double V,double *X,int * RANGE);
extern "C" void _stdcall  PV(double P  ,double *T, double *H, double *S,double V,double *X,int * RANGE);

extern "C" void _stdcall  PX2T(double P, double X,double *T,int * RANGE);
extern "C" void _stdcall  PX2H(double P, double X,double *H,int * RANGE);
extern "C" void _stdcall  PX2S(double P, double X,double *S,int * RANGE);
extern "C" void _stdcall  PX2V(double P, double X,double *V,int * RANGE);
extern "C" void _stdcall  PX(double P  ,double *T, double *H, double *S, double *V,double X,int * RANGE);

extern "C" void _stdcall  T2P(double T  ,double *P, int * RANGE);
extern "C" void _stdcall  T2HL(double T  ,double *H, int * RANGE);
extern "C" void _stdcall  T2HG(double T  ,double *H, int * RANGE);
extern "C" void _stdcall  T2SL(double T  ,double *S, int * RANGE);
extern "C" void _stdcall  T2SG(double T  ,double *S, int * RANGE);
extern "C" void _stdcall  T2VL(double T  ,double *V, int * RANGE);
extern "C" void _stdcall  T2VG(double T  ,double *V, int * RANGE);
extern "C" void _stdcall  T2L(double *P,double T  ,double *H,double *S,double *V,double *X,int * RANGE);
extern "C" void _stdcall  T2G(double *P,double T  ,double *H,double *S,double *V,double *X,int * RANGE);
extern "C" void _stdcall  T2CPL(double T  ,double *CP,int * RANGE);
extern "C" void _stdcall  T2CPG(double T  ,double *CP,int * RANGE);
extern "C" void _stdcall  T2CVL(double T  ,double *CV,int * RANGE);
extern "C" void _stdcall  T2CVG(double T  ,double *CV,int * RANGE);
extern "C" void _stdcall  T2EL(double T  ,double *E,int * RANGE);
extern "C" void _stdcall  T2EG(double T  ,double *E,int * RANGE);
extern "C" void _stdcall  T2SSPL(double T  ,double *SSP,int * RANGE);
extern "C" void _stdcall  T2SSPG(double T  ,double *SSP,int * RANGE);
extern "C" void _stdcall  T2KSL(double T  ,double *KS,int * RANGE);
extern "C" void _stdcall  T2KSG(double T  ,double *KS,int * RANGE);
extern "C" void _stdcall  T2ETAL(double T  ,double *ETA,int * RANGE);
extern "C" void _stdcall  T2ETAG(double T  ,double *ETA,int * RANGE);
extern "C" void _stdcall  T2UL(double T  ,double *U,int * RANGE);
extern "C" void _stdcall  T2UG(double T  ,double *U,int * RANGE);
extern "C" void _stdcall  T2RAMDL(double T  ,double *RAMD,int * RANGE);
extern "C" void _stdcall  T2RAMDG(double T  ,double *RAMD,int * RANGE);
extern "C" void _stdcall  T2PRNL(double T  ,double *PRN,int * RANGE);
extern "C" void _stdcall  T2PRNG(double T  ,double *PRN,int * RANGE);
extern "C" void _stdcall  T2EPSL(double T  ,double *EPS,int * RANGE);
extern "C" void _stdcall  T2EPSG(double T  ,double *EPS,int * RANGE);
extern "C" void _stdcall  T2NL(double T, double Lamd,double *N,int * RANGE);
extern "C" void _stdcall  T2NG(double T, double Lamd,double *N,int * RANGE);
extern "C" void _stdcall  T2SURFT(double T  ,double *SURFT,int * RANGE);

extern "C" void _stdcall  TH2P(double T, double H,double *P,int * RANGE);
extern "C" void _stdcall  TH2PLP(double T, double H,double *P,int * RANGE);
extern "C" void _stdcall  TH2PHP(double T, double H,double *P,int * RANGE);
extern "C" void _stdcall  TH2S(double T, double H,double *S,int * RANGE);
extern "C" void _stdcall  TH2SLP(double T, double H,double *S,int * RANGE);
extern "C" void _stdcall  TH2SHP(double T, double H,double *S,int * RANGE);
extern "C" void _stdcall  TH2V(double T, double H,double *V,int * RANGE);
extern "C" void _stdcall  TH2VLP(double T, double H,double *V,int * RANGE);
extern "C" void _stdcall  TH2VHP(double T, double H,double *V,int * RANGE);
extern "C" void _stdcall  TH2X(double T, double H,double *X,int * RANGE);
extern "C" void _stdcall  TH2XLP(double T, double H,double *X,int * RANGE);
extern "C" void _stdcall  TH2XHP(double T, double H,double *X,int * RANGE);
extern "C" void _stdcall  TH(double *P,double T, double H,double * S,double *V,double *X,int * RANGE);
extern "C" void _stdcall  THLP(double *P,double T, double H,double * S,double *V,double *X,int * RANGE);
extern "C" void _stdcall  THHP(double *P,double T, double H,double * S,double *V,double *X,int * RANGE);

extern "C" void _stdcall  TS2PHP(double T, double S,double *P,int * RANGE);
extern "C" void _stdcall  TS2PLP(double T, double S,double *P,int * RANGE);
extern "C" void _stdcall  TS2P(double T, double S,double *P,int * RANGE);
extern "C" void _stdcall  TS2HHP(double T, double S,double *H,int * RANGE);
extern "C" void _stdcall  TS2HLP(double T, double S,double *H,int * RANGE);
extern "C" void _stdcall  TS2H(double T, double S,double *H,int * RANGE);
extern "C" void _stdcall  TS2VHP(double T, double S,double *V,int * RANGE);
extern "C" void _stdcall  TS2VLP(double T, double S,double *V,int * RANGE);
extern "C" void _stdcall  TS2V(double T, double S,double *V,int * RANGE);
extern "C" void _stdcall  TS2X(double T, double S,double *X,int * RANGE);
extern "C" void _stdcall  TSHP(double *P,double T,double *H,double S,double *V, double *X,int * RANGE);
extern "C" void _stdcall  TSLP(double *P,double T,double *H,double S,double *V, double *X,int * RANGE);
extern "C" void _stdcall  TS(double *P,double T, double *H,double S,double *V, double *X,int * RANGE);

extern "C" void _stdcall  TV2P(double T, double V,double *P,int * RANGE);
extern "C" void _stdcall  TV2H(double T, double V,double *H,int * RANGE);
extern "C" void _stdcall  TV2S(double T, double V,double *S,int * RANGE);
extern "C" void _stdcall  TV2X(double T, double V,double *X,int * RANGE);
extern "C" void _stdcall  TV(double *P,double T  ,double *H,double *S ,double V , double *X,int * RANGE);

extern "C" void _stdcall  TX2P(double T, double X,double *P,int * RANGE);
extern "C" void _stdcall  TX2H(double T, double X,double *H,int * RANGE);
extern "C" void _stdcall  TX2S(double T, double X,double *S,int * RANGE);
extern "C" void _stdcall  TX2V(double T, double X,double *V,int * RANGE);
extern "C" void _stdcall  TX(double *P,double T  ,double *H, double *S, double *V,double X,int * RANGE);

extern "C" void _stdcall  H2TL(double H, double *T,int * RANGE);
extern "C" void _stdcall  HS2P(double H, double S,double *P,int * RANGE);
extern "C" void _stdcall  HS2T(double H, double S,double *T,int * RANGE);
extern "C" void _stdcall  HS2V(double H, double S,double *V,int * RANGE);
extern "C" void _stdcall  HS2X(double H, double S,double *X,int * RANGE);
extern "C" void _stdcall  HS(double *P, double *T,double H, double S,double *V, double *X,int * RANGE);

extern "C" void _stdcall  HV2P(double H, double V,double *P,int * RANGE);
extern "C" void _stdcall  HV2T(double H, double V,double *T,int * RANGE);
extern "C" void _stdcall  HV2S(double H, double V,double *S,int * RANGE);
extern "C" void _stdcall  HV2X(double H, double V,double *X,int * RANGE);
extern "C" void _stdcall  HV(double *P, double *T,double H,double *S,double V,double *X,int * RANGE);

extern "C" void _stdcall  HX2P(double H, double X,double *P,int * RANGE);
extern "C" void _stdcall  HX2PLP(double H, double X,double *P,int * RANGE);
extern "C" void _stdcall  HX2PHP(double H, double X,double *P,int * RANGE);
extern "C" void _stdcall  HX2T(double H, double X,double *T,int * RANGE);
extern "C" void _stdcall  HX2TLP(double H, double X,double *T,int * RANGE);
extern "C" void _stdcall  HX2THP(double H, double X,double *T,int * RANGE);
extern "C" void _stdcall  HX2S(double H, double X,double *S,int * RANGE);
extern "C" void _stdcall  HX2SLP(double H, double X,double *S,int * RANGE);
extern "C" void _stdcall  HX2SHP(double H, double X,double *S,int * RANGE);
extern "C" void _stdcall  HX2V(double H, double X,double *V,int * RANGE);
extern "C" void _stdcall  HX2VLP(double H, double X,double *V,int * RANGE);
extern "C" void _stdcall  HX2VHP(double H, double X,double *V,int * RANGE);
extern "C" void _stdcall  HX(double *P, double *T,double H,double *S, double *V,double X,int * RANGE);
extern "C" void _stdcall  HXLP(double *P, double *T,double H,double *S, double *V,double X,int * RANGE);
extern "C" void _stdcall  HXHP(double *P, double *T,double H,double *S, double *V,double X,int * RANGE);

extern "C" void _stdcall  S2TG(double S,double *T,int * RANGE);
extern "C" void _stdcall  SV2P(double S, double V,double *P,int * RANGE);
extern "C" void _stdcall  SV2T(double S, double V,double *T,int * RANGE);
extern "C" void _stdcall  SV2H(double S, double V,double *H,int * RANGE);
extern "C" void _stdcall  SV2X(double S, double V,double *X,int * RANGE);
extern "C" void _stdcall  SV(double *P, double *T, double *H,double S, double V,double *X,int * RANGE);

extern "C" void _stdcall  SX2P(double S, double X,double *P,int * RANGE);
extern "C" void _stdcall  SX2PLP(double S, double X,double *P,int * RANGE);
extern "C" void _stdcall  SX2PMP(double S, double X,double *P,int * RANGE);
extern "C" void _stdcall  SX2PHP(double S, double X,double *P,int * RANGE);
extern "C" void _stdcall  SX2TLP(double S, double X,double *T,int * RANGE);
extern "C" void _stdcall  SX2TMP(double S, double X,double *T,int * RANGE);
extern "C" void _stdcall  SX2THP(double S, double X,double *T,int * RANGE);
extern "C" void _stdcall  SX2T(double S, double X,double *T,int * RANGE);
extern "C" void _stdcall  SX2H(double S, double X,double *H,int * RANGE);
extern "C" void _stdcall  SX2HLP(double S, double X,double *H,int * RANGE);
extern "C" void _stdcall  SX2HMP(double S, double X,double *H,int * RANGE);
extern "C" void _stdcall  SX2HHP(double S, double X,double *H,int * RANGE);
extern "C" void _stdcall  SX2V(double S, double X,double *V,int * RANGE);
extern "C" void _stdcall  SX2VLP(double S, double X,double *V,int * RANGE);
extern "C" void _stdcall  SX2VMP(double S, double X,double *V,int * RANGE);
extern "C" void _stdcall  SX2VHP(double S, double X,double *V,int * RANGE);
extern "C" void _stdcall  SX(double *P, double *T, double *H,double S,double *V,double X,int * RANGE);
extern "C" void _stdcall  SXLP(double *P, double *T, double *H,double S,double *V,double X,int * RANGE);
extern "C" void _stdcall  SXMP(double *P, double *T, double *H,double S,double *V,double X,int * RANGE);
extern "C" void _stdcall  SXHP(double *P, double *T, double *H,double S,double *V,double X,int * RANGE);

extern "C" void _stdcall  V2TG(double V , double *T,int * RANGE);
extern "C" void _stdcall  VX2P(double V, double X,double *P,int * RANGE);
extern "C" void _stdcall  VX2PLP(double V, double X,double *P,int * RANGE);
extern "C" void _stdcall  VX2PHP(double V, double X,double *P,int * RANGE);
extern "C" void _stdcall  VX2T(double V, double X,double *T,int * RANGE);
extern "C" void _stdcall  VX2TLP(double V, double X,double *T,int * RANGE);
extern "C" void _stdcall  VX2THP(double V, double X,double *T,int * RANGE);
extern "C" void _stdcall  VX2H(double V, double X,double *H,int * RANGE);
extern "C" void _stdcall  VX2HLP(double V, double X,double *H,int * RANGE);
extern "C" void _stdcall  VX2HHP(double V, double X,double *H,int * RANGE);
extern "C" void _stdcall  VX2S(double V, double X,double *S,int * RANGE);
extern "C" void _stdcall  VX2SLP(double V, double X,double *S,int * RANGE);
extern "C" void _stdcall  VX2SHP(double V, double X,double *S,int * RANGE);
extern "C" void _stdcall  VX(double *P, double *T, double *H, double *S,double V, double X,int * RANGE);
extern "C" void _stdcall  VXLP(double *P, double *T, double *H, double *S,double V, double X,int * RANGE);
extern "C" void _stdcall  VXHP(double *P, double *T, double *H, double *S,double V, double X,int * RANGE);