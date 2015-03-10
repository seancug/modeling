/*ȫ�ֱ���*/
#ifndef topoft_h
#include"topoft.h"
#endif
#ifndef glo_var_h
#define glo_var_h

string TITLE_NAME;
string VERBOSE;
/*2-D������*/
int   NX, NZ; //NX��NZ��ģ�������С��������֮�����ΪDH
float DX, DZ; //DX, DZ:�ռ䲽������λ(m)
int free_depth = 1;
/*ʱ�䲽��*/
float TIME, DT;// TIME�����𲨴�����ʱ�䣬��λ(��),DT��ʱ�䲽��, ��λ(��)
int NT;//NT��ʱ�䲽����һ����λ����һ��������

/*��Դ����*/
string SOURCE_WAVELET;
string SOURCE_TYP;// SOURCE_WAVELET:��Դ����;SOURCE_TYP����ը��Դ = 0������Դ ��x = ��1����Z�� = 2
char SOURCE_FILE[STRING_SIZE];//SOURCE_FILE����Դ�����ļ�
float ANGLE, F0, XS, ZS = 0.0;//f0����Դ��Ƶ����λ(��);XS��YS��ZS��Դλ�ã���λ(m)
float AMP = 1.0;
/*ģ��*/
string READMOD; //READMOD�������ļ���ȡģ�� = 'NO'�����ļ���ȡģ�� = 'YES'
char MFILE_VP[STRING_SIZE], MFILE_VS[STRING_SIZE], MFILE_DEN[STRING_SIZE];//MFILE�����ģ���ļ�
string MD_FORMAT; //ģ���ļ���ʽ
/*���ɱ߽�*/
string FREE_SURF;// FREE_SURF�����ɱ������� ��ˮƽ = 0�� ������ = 1���ɱ�������(�������ɱ߽�--'NO_FREE', ����ˮƽ--'IMAGING_PLANER', 
                  //���������--'IMAGING_TOPOGRAPHY', ˮƽAEA����--'AEA_PLANER', AEA�������--'AEA_TOPOGRAPHY')

/*�������ɱ߽�*/
float FW;//FW:���ձ߽���,��λ(m)
float  DAMPING; //DAMPING��˥��ϵ��

/*����*/
string FILEDS;//SNAP��������棬����� = 0������ٶ� = 1�����Ӧ���� = 2
float TSNAP1;    //TSNAP1����һ�Ų��ο���ʱ�䣬��λ(��)
float TSNAP2;    //TSNAP2������Ų��ο���ʱ�䣬��λ(��)
float TSNAPINC;  //TSNAPINC��֮�䲨�����ʱ��������λ(��)
string  SNAP_FORMAT; //SNAP_FORMAT�������ʽ(ASCII(2);BINARY(3)) = 3
char SNAP_FILE[STRING_SIZE]; //SNAP_FILE����ȡ�����ļ���

/*�첨������*/
string SEISMO;// SEISMO:����������¼=0������ٶȳ�=1�����ˮ�µ���첨��=2 
string READREC; //READREC�������ļ���ȡ�첨��λ�� = 0�����ļ���ȡ�첨��λ�� = 1
char REC_FILE[STRING_SIZE];//REC_FILE���첨��λ���ļ�
char RECIN[STRING_SIZE];
string RECTYPE; //���沼�ü첨��Ϊ0��VSPΪ 1
float XREC1, XREC2, ZREC1 = 0.0, ZREC2 = 0.0; // XREC1,ZREC1����һ���첨��λ��,XREC2, ZREC2�����һ���첨��λ��

float   DRX, DRZ; //DRX,DRZ��X,Z����첨��֮��ľ���

/*�����¼*/
int   NDT; //NDT�����������ز�������λ(����)
string SEIS_FORMAT;//SEIS_FORMAT�������¼�����ʽ��SU(1); ASCII(2); BINARY(3)
char     SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VZ[STRING_SIZE];// SEIS_FILE_VX�����VX�ٶȵ��ļ�;SEIS_FILE_VZ�����VZ�ٶȵ��ļ�


#endif
