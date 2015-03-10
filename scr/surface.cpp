#include"topoft.h"

void surface(int depth, float ** vx, float **vz, float ** sxx, float **sxz, float **szz,
	float ** rden, float **lam, float **u)
{
	extern float DT, DX, DZ;
	extern int NX;

	/*自由表面深度*/
	int j = depth;
	float vx_dx = 0, vz_dz = 0;
	float sxx_dx = 0, sxz_dz = 0;
	float vx_dz = 0, vz_dx = 0;
	for (int i = 1; i <= NX; i++)
	{
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

		/*****sxz 迭代刚刚加上去，需要重新测试*****/
		/*速度空间偏导*/
		vx_dz = (
			a1*(vx[i][j + 1] - vx[i][j]) +
			a2*(vx[i][j + 2] - vx[i][j - 1]) +
			a3*(vx[i][j + 3] - vx[i][j - 2]) +
			a4*(vx[i][j + 4] - vx[i][j - 3])) / DZ;

		vz_dx = (
			a1*(vz[i + 1][j] - vz[i][j]) +
			a2*(vz[i + 2][j] - vz[i - 1][j]) +
			a3*(vz[i + 3][j] - vz[i - 2][j]) +
			a4*(vz[i + 4][j] - vz[i - 3][j])) / DX;

		/*sxz 剪应力迭代*/
		sxz[i][j] += DT*rden[i][j] * (vx_dz + vz_dx);
	}
}
