#include"topoft.h"
#include "segy.h"

string outseis(FILE *fpdata, float **section, int **recpos, int ntr, int ns, string seis_form)
{

	/* declaration of extern variables */
	extern int NDT;
	extern float XS, ZS, DX, DZ, TIME, DT;

	string message;
	/* declaration of local variables */
	int i, j;
	//segy tr;
	//int tracl;
	//float xr, yr, zr, x, z;
	
	if (seis_form == "SEGY")
	{
		segy tr;
		bhed bh;
		
		Write_Dedault_Seg_Y_head(fpdata);
		fseek(fpdata, 3200L, SEEK_SET);
		/*设置400字节头*/
		bh.format = 5;       /* IEEE 32 float bit format*/
		bh.lino = 1;         /* Line number */
		bh.reno = 1;         /* Reel number */
		bh.nart = 0;         /* Number of aux traces per rec.*/
		bh.hdt = (short)((NDT*DT)*1.0e6);
		bh.tsort = 1;
		bh.schn = ntr;        /* Traces per record*/
		bh.hns = (short)ns;   /* Samples per trace          */
		bh.ntrpr = 1;
		bh.mfeet = 1;

		fwrite(&bh, 400, 1, fpdata);   //400字节道头

		float xr, zr;
		/*BINARY,二进制文件,功能不完善 */
		for (i = 1; i <= ntr; i++)
		{
			xr = recpos[1][ntr] * DX;
			zr = recpos[2][ntr] * DZ;
			tr.tracl = i;      /* trace sequence number within line */
			tr.cdp = i;
			tr.trid = 1;           /* trace identification code: 1=seismic*/
			tr.offset = iround(sqrt((XS - xr)*(XS - xr)));
			tr.sdepth = iround(ZS);

			tr.scalel = (short)-3;
			tr.scalco = (short)-3;
			tr.sx = iround(XS);  /* X source coordinate */
			tr.sy = 0;

			/* group coordinates */
			tr.gx = 0;
			tr.gy = 0;
			tr.ns = (unsigned short)ns; /* number of samples in this trace */
			tr.dt = (unsigned short)((NDT*DT)*1.0e6); /* sample interval in micro-seconds */
			tr.d1 = 0.0;        /* sample spacing for non-seismic data */
			tr.f1 = 0.0;              /* first sample location for non-seismic data */
			tr.d2 = 0.0;        /* sample spacing between traces */
			tr.f2 = 0.0;

			fwrite(&tr, 240, 1, fpdata); //道头
			fwrite(section[i], 4, ns, fpdata);//地震数据
		}
	}
	else if (seis_form == "ASCII")
	{
		for (i = 1; i <= ntr; i++)  /*ASCII ONE COLUMN*/
		{
			for (j = 1; j <= ns; j++)
				fprintf(fpdata, "%e\t", section[i][j]);
			fprintf(fpdata, "\n");
		}
	}
	else if (seis_form == "SU")
	{
	
	}
	else
	{
		message.append(" Don't know the format for the seismogramm-data !\n");
		message.append(" No output written. ");
	}
	fclose(fpdata);
	return message;
}
