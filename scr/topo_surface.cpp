#include"topoft.h"

void topo_surface(int **flag, float ** vx, float **vz, float ** sxx, float **sxz, float **szz,
	float ** rden, float **lam, float **u, int depth)
{
	extern float DT, DX, DZ;
	extern int NX, NZ;

	float vx_dx = 0, vz_dz = 0;
	float sxx_dx = 0, sxz_dz = 0;
	float vp = 0;

	if (depth == -1)
		depth = NZ;

	/*避免角落网格重叠，对水平方向和垂直方向*/
	for (int i = 1; i <= NX; i++)
	{
		for (int j = 1; j <= depth; j++)
		{
			switch (flag[i][j])
			{
			case 1:/*H-point,正常镜像法*/
				/*镜像法 imageing stress*/
				szz[i][j] = 0;
				szz[i][j - 1] = -szz[i][j + 1];
				szz[i][j - 2] = -szz[i][j + 2];
				szz[i][j - 3] = -szz[i][j + 3];

				sxz[i][j - 1] = -sxz[i][j];
				sxz[i][j - 2] = -sxz[i][j + 1];
				sxz[i][j - 3] = -sxz[i][j + 2];
				sxz[i][j - 4] = -sxz[i][j + 3];

				vx_dx =
					(
					a1*(vx[i][j] - vx[i - 1][j]) +
					a2*(vx[i + 1][j] - vx[i - 2][j]) +
					a3*(vx[i + 2][j] - vx[i - 3][j]) +
					a4*(vx[i + 3][j] - vx[i - 4][j])) / DX;
				/*Dvz用Dvx代替*/
				vz_dz = (-lam[i][j] / (lam[i][j] + 2.0*u[i][j]))*vx_dx;

				/*更新sxx*/
				sxx[i][j] += DT*((lam[i][j] + 2.0 * u[i][j])*vx_dx + lam[i][j] * vz_dz);
				break;
			case 2:

				break;
			case 3:

				break;
			case 4:

				break;
			case 5:

				break;
			case 6:

				break;
			case 7:
				break;
			}

		}
	}
}
