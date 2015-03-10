#include <string.h>
#include <stdlib.h>
#include"topoft.h"

void check_par()
{
	extern float DX, DZ, TIME, DT, ANGLE, F0, XS, ZS;
	extern float TSNAP1, TSNAP2, TSNAPINC, FW;
	extern string RECTYPE;
	extern float XREC1, XREC2, ZREC1, ZREC2;
	extern float   DRX, DRZ;
	extern float  DAMPING;
	extern int    NDT;
	extern int   NX, NZ;
	extern string SOURCE_WAVELET, SOURCE_TYP, READMOD, MD_FORMAT, FREE_SURF, FILEDS, SNAP_FORMAT, READREC, SEISMO, SEIS_FORMAT;
	extern char  SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE];
	extern char  MFILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char  SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VZ[STRING_SIZE];
	string message;

	bool out = false;
	if (TIME == 0.0)
	{
		parerror("ti in TIME");
		out = true;
	}
	if (DT == 0)
	{
		parerror("dt in TIME");
		out = true;
	}
	if (XS == 0 || ZS == 0)
	{
		parerror("coorid in SOURCE");
		out = true;
	}
	if (SOURCE_WAVELET != "FILE"&&SOURCE_TYP == "POINT"&&F0 == 0)
	{
		parerror("f0 in SOURCE");
		out = true;
	}
	if (FILEDS != "NO"&&SNAP_FILE == 0)
		strcpy(SNAP_FILE, "waveFiled");
	if (SEISMO != "NO"&&SEIS_FILE_VX == 0)
		strcpy(SNAP_FILE, "vx");
	if (SEISMO != "NO"&&SEIS_FILE_VZ == 0)
		strcpy(SNAP_FILE, "vz");
	if (out)
	{
		exit(1);
	}
	
}
