#include"topoft.h"

void caculate_coff(float **rden, float **lam, float **u, float **vpIn, float **vsIn, float **rhoIn)
{
	extern float DX, DZ, FW;
	extern int NX, NZ;

	int ifw_x, ifw_z;
	ifw_x = iround(FW / DX);
	ifw_z = iround(FW / DZ);

	int NX_f = 1, NX_e = NX, NZ_f = 1, NZ_e = NZ;
	for (int i = NX_f + ifw_x; i <= NX_e - ifw_x; i++)
	{
		for (int j = NZ_f; j <= NZ_e - ifw_z; j++)
		{
			if (rhoIn[i][j] == 0)
				rden[i][j] = 0;
			else
				rden[i][j] = 1 / (rhoIn[i][j] * 1000.0); //密度倒数
			lam[i][j] = (pow(vpIn[i][j], 2) - 2 * pow(vsIn[i][j], 2))*(rhoIn[i][j] * 1000.0); //拉梅常数
			u[i][j] = pow(vsIn[i][j], 2)*(rhoIn[i][j] * 1000.0);
		}
	}

	/*PML 区域延伸过来*/
	int ii = 0;
	for (int i = NX_f - 1; i < NX_f + ifw_x; i++)
	{
		ii = NX_e - i + NX_f;
		for (int j = NZ_f; j <= NZ_e - ifw_z; j++)
		{
			vpIn[i][j] = vpIn[NX_f + ifw_x][j];
			vpIn[ii][j] = vpIn[NX_e - ifw_x][j];

			rden[i][j] = rden[NX_f + ifw_x][j];
			rden[ii][j] = rden[NX_e - ifw_x][j];

			lam[i][j] = lam[NX_f + ifw_x][j];
			lam[ii][j] = lam[NX_e - ifw_x][j];

			u[i][j] = u[NX_f + ifw_x][j];
			u[ii][j] = u[NX_e - ifw_x][j];
		}
	}

	for (int j = NZ_e - ifw_z + 1; j <= NZ_e + 1; j++)
	{
		for (int i = NX_f - 1; i <= NX_e + 1; i++)
		{
			vpIn[i][j] = vpIn[i][NZ_e - ifw_z];
			rden[i][j] = rden[i][NZ_e - ifw_z];
			lam[i][j] = lam[i][NZ_e - ifw_z];
			u[i][j] = u[i][NZ_e - ifw_z];
		}
	}
}