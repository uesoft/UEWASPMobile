//WetAir.h

//************************************************************
//     作    者：                                            *
//              长沙优易软件开发有限公司(UESoft Corp.) 李国斌*
//     文件名称：                                            * 
//                DLL外部接口(湿空气性质计算焓湿图动态连接库)*
//                Psychrometric charts(焓湿图)计算动态连接库 *
//     完成时间：                                            *
//                2006年5月                                  *
//************************************************************ 

//各种气体的数据
typedef struct GasData
{
	short iSEQ;//gas sequence number
	char* cName;//gas name,e.g. H2,AIR,etc.
	double fMolecular_Weigh;//gas molecular weight,kg/kmol
	//Molar Specific Heats formula coefficient of true Gases
	//at constant pressure are as follows
	double a0;//temperature coefficient 
	double a1;//temperature coefficient 
	double a2;//temperature coefficient 
	double a3;//temperature coefficient 
	double fMinTemp;//lower temperature range
	double fMaxTemp;//upper temperature range
	double fMaxError;//max error of specific heat 
} mGasData;


extern "C" extern "C" double _stdcall  WetAir_GetRva();

extern "C" void _stdcall  WetAir_SETSTD_WASP(int  STDID);//为湿空气计算设定水蒸气公式IF97/IFC-67
extern "C" void _stdcall  WetAir_GETSTD_WASP(int  *STDID);

extern "C" void _stdcall  WetAir_SETSTD_WetBulbSurface_AirVelocity(double fWetBulbSurface_AirVelocity);
extern "C" void _stdcall  WetAir_SETSTD_atmosphere_pressure(double fLocalAtmospherePressure);
extern "C" void _stdcall  WetAir_GETSTD_atmosphere_pressure(double fLocalAtmospherePressure);

extern "C" void _stdcall  WetAir_SETSTD_WASP(int  STDID);
extern "C" void _stdcall  WetAir_GETSTD_WASP(int *STDID);
extern "C" short _stdcall  WetAir_GetIndex(char* cName);
extern "C" double _stdcall  WetAir_GetRva();
extern "C" double _stdcall  WetAir_GetCp(short i,double T);
extern "C" double _stdcall  WetAir_GetEnthalpy(short i,double T);
extern "C" double _stdcall  WetAir_D2P(double d);
extern "C" double _stdcall  WetAir_P2D(double P);
extern "C" double _stdcall  WetAir_GetEnthalpy_wetair(double T,double W);
//已知干球温度t(℃)和湿球温度tw(℃)，求水汽分压P(MPa)
extern "C" double _stdcall  WetAir_TT2P(double T, double Tw, double phi1);
//已知干球温度t(℃)和绝热饱和温度tw(℃)，求绝对湿度D(kg/kg.dry)
extern "C" double _stdcall  WetAir_TT2D(double T, double Tw, double phi1);

