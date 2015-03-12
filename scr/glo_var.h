/*全局变量*/
#ifndef topoft_h
#include"topoft.h"
#endif
#ifndef glo_var_h
#define glo_var_h

string TITLE_NAME;
string VERBOSE;
/*2-D网格定义*/
int   NX, NZ; //NX，NZ：模型网格大小，个数，之间距离为DH
float DX, DZ; //DX, DZ:空间步长，单位(m)
int free_depth = 1;
/*时间步长*/
float TIME, DT;// TIME：地震波传播的时间，单位(秒),DT：时间步长, 单位(秒)
int NT;//NT：时间步数，一个单位代表一个采样率

/*震源参数*/
string SOURCE_WAVELET;
string SOURCE_TYP;// SOURCE_WAVELET:震源函数;SOURCE_TYP：爆炸震源 = 0，点震源 在x = 上1，在Z上 = 2
char SOURCE_FILE[STRING_SIZE];//SOURCE_FILE：震源参数文件
float ANGLE, F0, XS, ZS = 0.0;//f0：震源主频，单位(秒);XS，YS，ZS震源位置，单位(m)
float AMP = 1.0;
/*模型*/
string READMOD; //READMOD：不从文件读取模型 = 'NO'，从文件读取模型 = 'YES'
char MFILE_VP[STRING_SIZE], MFILE_VS[STRING_SIZE], MFILE_DEN[STRING_SIZE];//MFILE：存放模型文件
string MD_FORMAT; //模型文件格式
/*自由边界*/
string FREE_SURF;// FREE_SURF：自由表面条件 ，水平 = 0， 带地形 = 1自由表面条件(不是自由边界--'NO_FREE', 镜像水平--'IMAGING_PLANER', 
                  //镜像带地形--'IMAGING_TOPOGRAPHY', 水平AEA方案--'AEA_PLANER', AEA起伏方案--'AEA_TOPOGRAPHY')

/*吸收自由边界*/
float FW;//FW:吸收边界宽度,单位(m)
float  DAMPING; //DAMPING：衰减系数

/*快照*/
string FILEDS;//SNAP：输出剖面，不输出 = 0，输出速度 = 1，输出应力场 = 2
float TSNAP1;    //TSNAP1：第一张波形快照时间，单位(秒)
float TSNAP2;    //TSNAP2：最后张波形快照时间，单位(秒)
float TSNAPINC;  //TSNAPINC：之间波场输出时间间隔，单位(秒)
string  SNAP_FORMAT; //SNAP_FORMAT：输出格式(ASCII(2);BINARY(3)) = 3
char SNAP_FILE[STRING_SIZE]; //SNAP_FILE：存取快照文件名

/*检波器布置*/
string SEISMO;// SEISMO:不输出地震记录=0，输出速度场=1，输出水下地震检波器=2 
string READREC; //READREC：不从文件读取检波器位置 = 0，从文件读取检波器位置 = 1
char REC_FILE[STRING_SIZE];//REC_FILE：检波器位置文件
char RECIN[STRING_SIZE];
string RECTYPE; //地面布置检波器为0，VSP为 1
float XREC1, XREC2, ZREC1 = 0.0, ZREC2 = 0.0; // XREC1,ZREC1：第一个检波器位置,XREC2, ZREC2：最后一个检波器位置

float   DRX, DRZ; //DRX,DRZ：X,Z方向检波器之间的距离

/*地震记录*/
int   NDT; //NDT：地震道输出重采样，单位(个数)
string SEIS_FORMAT;//SEIS_FORMAT：地震记录输出格式，SU(1); ASCII(2); BINARY(3)
char     SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VZ[STRING_SIZE];// SEIS_FILE_VX：存放VX速度的文件;SEIS_FILE_VZ：存放VZ速度的文件


#endif
