#include"topoft.h"
#include "segy.h"
/*
*参数：
*    rescpos:检波器坐标，输出segy数据所用，这里还没使用到   
*/
string saveSeis(float **sectionvx, float **sectionvz, int  **recpos, int ntr, int ns)
{

	extern string SEISMO, SEIS_FORMAT;
	extern char  SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VZ[STRING_SIZE];
	string message;
	if (SEISMO == "V")
	{
		FILE *fpvx = NULL;
		FILE *fpvz = NULL;
        char temvx[124];
        char temvz[124];
        sprintf(temvx, "%s%s","../gather/", SEIS_FILE_VX);
        sprintf(temvz, "%s%s","../gather/", SEIS_FILE_VZ);
		if (SEIS_FORMAT == "SEGY" || SEIS_FORMAT == "SU")
		{
            fpvx = fopen(temvx, "wb");
            fpvz = fopen(temvz, "wb");
		}
		else if (SEIS_FORMAT == "ASCII")
		{
            fpvx = fopen(temvx, "w");
            fpvz = fopen(temvz, "w");
		}
		char temp[252] = " ";
		sprintf(temp, "\tWriting %d traces seismograms (vx) into\n\t %s \n", ntr, SEIS_FILE_VX);
		message.append(temp);
		message.append(outseis(fpvx, sectionvx, recpos, ntr, ns, SEIS_FORMAT));
		sprintf(temp, "\tWriting %d traces seismograms (vz) into\n\t %s \n", ntr, SEIS_FILE_VZ);
		message.append(temp);
		message.append(outseis(fpvz, sectionvz, recpos, ntr, ns, SEIS_FORMAT));
	}
	else if (SEISMO == "V")
	{

	}
	else if (SEISMO == "SV" || SEISMO == "VS")
	{

	}
	/*写入水听器，还没有写入*/
	/*if (SEISMO == 2){
		fprintf(fp, " PE %d is writing %d seismograms of pressure into\n\t %s \n", ntr, SEIS_FILE_VX);
		outseis(fp, fopen(SEIS_FILE_VX, "w"), sectionvx, recpos, ntr, ns, SEIS_FORMAT);
		}*/
	return message;
}
