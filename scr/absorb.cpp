#include"topoft.h"

void absorb(float **absorb_coeff, int NX_f, int NX_e, int NZ_f, int NZ_e)
{
	extern float DX, DZ, FW, DAMPING;
	extern int NX, NZ;

	/*�м����*/
	float amp, ax, az, *coeff_x, *coeff_z;

	/*���ղ���*/
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
	/*����˥��ϵ��*/
	for (int i = NX_f; i < NX_f + ifw_x; i++)
		coeff_x[i] = exp(-(ax*ax*(NX_f + ifw_x - i)*(NX_f + ifw_x - i)));

	for (int i = NZ_f; i < NZ_f + ifw_z; i++)
		coeff_z[i] = exp(-(az*az*(NZ_f + ifw_z - i)*(NZ_f + ifw_z - i)));

	for (int i = NX_f; i <= NX_e; i++)
		for (int j = NZ_f; j <= NZ_e; j++)
			absorb_coeff[i][j] = 1.0;

	int ii = 0;
	/*x ����������������ϵ��*/
	for (int i = NX_f; i < NX_f + ifw_x; i++)
	{
		ii = NX_e - i + NX_f;
		for (int j = NZ_f; j <= NZ_e; j++)
			absorb_coeff[i][j] = absorb_coeff[ii][j] = coeff_x[i];
	}
		
	/*����������ɵر�����������ձ߽�����*/
	if (false)
	{
		for (int j = NZ_f; j < NZ_f + ifw_z; j++)
		{
			for (int i = NX_f; i <= NX_e; i++)
				absorb_coeff[i][j] = coeff_z[j];
		}
	}
	int jj = 0;
	/*z ����������˥��ϵ��*/
	for (int j = NZ_f; j < NZ_f + ifw_z; j++)
	{
		jj = NZ_e - j + NZ_f;
		for (int i = NX_f; i <= NX_e; i++)
			absorb_coeff[i][jj] = coeff_z[j];
	}
	
	free_dvector(coeff_x, NX_f, NX_f + ifw_x);
	free_dvector(coeff_z, NZ_f, NZ_f + ifw_z);
}
