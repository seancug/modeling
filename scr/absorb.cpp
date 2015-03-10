#include"topoft.h"

void absorb(float **absorb_coeff, int NX_f, int NX_e, int NZ_f, int NZ_e)
{
	extern float DX, DZ, FW, DAMPING;
	extern int NX, NZ;

	/*中间变量*/
	float amp, ax, az, *coeff_x, *coeff_z;

	/*吸收层数*/
	int ifw_x, ifw_z;
	float damp = 8.0;
	amp = 1 - damp / 100.0;
	ifw_x = iround(FW / DX);
	ifw_z = iround(FW / DZ);

	coeff_x = dvector(NX_f, NX_f + ifw_x);
	coeff_z = dvector(NZ_f, NZ_f + ifw_z);
	ax = sqrt(-log(amp) / ((ifw_x-1)*(ifw_x-1)));
	az = sqrt(-log(amp) / ((ifw_z-1)*(ifw_z-1)));
	float b = log(amp);
	/*计算衰减系数*/
	for (int i = NX_f; i < NX_f + ifw_x; i++)
		coeff_x[i] = exp(-(ax*ax*(NX_f + ifw_x - i)*(NX_f + ifw_x - i)));

	for (int i = NZ_f; i < NZ_f + ifw_z; i++)
		coeff_z[i] = exp(-(az*az*(NZ_f + ifw_z - i)*(NZ_f + ifw_z - i)));

	for (int i = NX_f; i <= NX_e; i++)
		for (int j = NZ_f; j <= NZ_e; j++)
			absorb_coeff[i][j] = 1.0;

	int ii = 0;
	/*x 方向即左右生成吸收系数*/
	for (int i = NX_f; i < NX_f + ifw_x; i++)
	{
		ii = NX_e - i + NX_f;
		for (int j = NZ_f; j <= NZ_e; j++)
			absorb_coeff[i][j] = absorb_coeff[ii][j] = coeff_x[i];
	}
		
	/*如果不是自由地表上面加上吸收边界条件*/
	if (false)
	{
		for (int j = NZ_f; j < NZ_f + ifw_z; j++)
		{
			for (int i = NX_f; i <= NX_e; i++)
				absorb_coeff[i][j] = coeff_z[j];
		}
	}
	int jj = 0;
	/*z 方向即下生成衰减系数*/
	for (int j = NZ_f; j < NZ_f + ifw_z; j++)
	{
		jj = NZ_e - j + NZ_f;
		for (int i = NX_f; i <= NX_e; i++)
			absorb_coeff[i][jj] = coeff_z[j];
	}
	
	free_dvector(coeff_x, NX_f, NX_f + ifw_x);
	free_dvector(coeff_z, NZ_f, NZ_f + ifw_z);
}
