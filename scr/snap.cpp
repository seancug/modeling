#include"topoft.h"
/*
������
FILE *fp log�ļ�
*/
string snap( int nt, string format, string type,
	float **vx, float **vz, int idx, int idz, int nx1, int nz1, int nx2, int nz2)
{
	/*
	different data formats of output available:
	format=GRID  :  Sufer��ʽ,����δ���
	format=ASCII  :  ASCII,�ı���ʽ
	format=BINARY  :  BINARY,������,SU,ximage����ʹ��

	different types:
	type=V : values in vx, and vz
	���¹���û��д
	type=S : -(vx+vy+vz) (pressure field)
	type=3 : divergence of vx, vy and vz (energy of compressional waves)
	and curl of vx, vy and vz (energy of shear waves)
	type=4 : both particle velocities (type=1) and energy (type=3)
	*/

	string message;
	char xfile[STRING_SIZE], zfile[STRING_SIZE];
	char ext[8];

	FILE *fpx1, *fpz1; //��������ٶ��ļ�ָ��

	extern float DT;
	extern char SNAP_FILE[STRING_SIZE];

	if (format == "ASCII")
	{
		sprintf(ext, ".mat");
	}
	else if (format == "BINARY")
	{
		sprintf(ext, ".dat");
	}
	else if (format == "GRID")
	{
		sprintf(ext, ".grd");
	}


    sprintf(xfile, "%s%s_%.2f_vx%s", "../snap/", SNAP_FILE, nt*DT, ext);
    sprintf(zfile, "%s%s_%.2f_vz%s", "../snap/", SNAP_FILE, nt*DT, ext);
	char in[252]="";
	sprintf(in, "\nWriting snapshot-data at T=%.2f(s) to \n", nt*DT);
	message.append(in);
	
	if (type == "V")
	{
		sprintf(in, "\t%s\n", xfile);
		message.append(in);
		sprintf(in, "\t%s\n\n", zfile);
		message.append(in);

		if (format == "ASCII")
		{
			fpx1 = fopen(xfile, "w");
			fpz1 = fopen(zfile, "w");
			for (int k = nz1; k <= nz2; k += idz)
			{
				/*����Ϊx������Ϊz*/
				for (int i = nx1; i <= nx2; i += idx)
				{
					writedsk(fpx1, vx[i][k], format);
					writedsk(fpz1, vz[i][k], format);
				}
				fprintf(fpx1, "\n");
				fprintf(fpz1, "\n");
			}
			fclose(fpx1);
			fclose(fpz1);
		}
		else if (format == "BINARY")
		{
			fpx1 = fopen(xfile, "wb");
			fpz1 = fopen(zfile, "wb");
			/*����Ϊx������Ϊz*/
			for (int i = nx1; i <= nx2; i += idx)
			{
				for (int k = nz1; k <= nz2; k += idz)
				{
					writedsk(fpx1, vx[i][k], format);
					writedsk(fpz1, vz[i][k], format);
				}
			}
			fclose(fpx1);
			fclose(fpz1);
		}
		else if (format == "GRID")
		{
		}
		else
		{
			cout << "\tSnap Shot format is not correct!" << endl;
		}
	}
	else if (type == "S")
	{

	}
	else if (type == "SV" || type == "VS")
	{

	}
	return message;
}
