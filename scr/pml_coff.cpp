#include"topoft.h"

void pml_coff(pml_c pml, float Vs, int NX_f, int NX_e, int NZ_f, int NZ_e)
{
	extern float DX, DZ, FW, DAMPING;
	extern int NX, NZ;

	/*中间变量*/
	float N;
	/*吸收层数*/
	float ifw_x, ifw_z;

	float d0 = 0.0, R = 0.0, sgam = 0.0;
	ifw_x = iround(FW / DX);
	ifw_z = iround(FW / DZ);
	N = 2;
	R = DAMPING;
	sgam = 3.5;

	d0 = log10(1 / R)*sgam*Vs / (ifw_x*DX);
	/*生成pml系数*/
	for (int i = NX_f; i < NX_f + ifw_x; i++)
	{
		pml.left[i] = d0*pow((NX_f + ifw_x - i) / ifw_x, N);
		pml.leftv[i] = d0*pow((NX_f + ifw_x - i - 0.5) / ifw_x, N);//差半个网格，速度与应力分量
		//pml.left[i] = (A*Vs*pow((NX_f + ifw_x - i) / ifw_x, N)) / (ifw_x*DX);
	}

	for (int i = NX_e - ifw_x + 1; i <= NX_e; i++)
	{
		pml.right[i] = d0*pow((i - NX_e + ifw_x - 0.5) / ifw_x, N);
		pml.rightv[i] = d0*pow((i - NX_e + ifw_x) / ifw_x, N);
		//pml.right[i] = (A*Vs*pow((i - NX_e + ifw_x) / ifw_x, N)) / (ifw_x*DX);
	}

	for (int i = NZ_e - ifw_z + 1; i <= NZ_e; i++)
	{
		pml.bottom[i] = d0*pow((i - NZ_e + ifw_z - 0.5) / ifw_z, N);
		pml.bottomv[i] = d0*pow((i - NZ_e + ifw_z) / ifw_z, N);
		//pml.bottom[i] = (A*Vs*pow((i - NZ_e + ifw_z) / ifw_z, N)) / (ifw_z*DZ);
	}
}