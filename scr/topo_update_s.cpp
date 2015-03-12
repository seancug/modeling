
#include"topoft.h"
/*
����ر�Ӧ������
*/
void topo_update_s(int nt, int **flag, float ** vx, float ** vz,
	float ** sxx, float ** szz, float ** sxz,
	float ** lam, float ** u)
{
	/*ȫ�ֱ���*/
	extern float DT, DX, DZ;
	extern int   NX, NZ;

	/*�м����*/
	float vx_dx = 0, vz_dz = 0;
	float vx_dz = 0, vz_dx = 0;
	float eff_u_xz = 0;
	int m, n;
	for ( m = 1; m <= NX; m++)
	{
		for ( n = 1; n <= NZ; n++)
		{
			if (flag[m][n] == 0)
			{
				/*��Ч����ģ��������ƽ��*/
				eff_u_xz = 1 / ((1 / u[m][n] + 1 / u[m + 1][n] + 1 / u[m][n + 1] + 1 / u[m + 1][n + 1]) / 4);
			}
			else
			{
				eff_u_xz = u[m][n];
			}

			switch (flag[m][n])
			{
			case -1:
				break;
			case 0:/*�ڲ�����*/
				/*�ٶȿռ�ƫ��*/
				vx_dx = (
					a1*(vx[m][n] - vx[m - 1][n]) +
					a2*(vx[m + 1][n] - vx[m - 2][n]) +
					a3*(vx[m + 2][n] - vx[m - 3][n]) +
					a4*(vx[m + 3][n] - vx[m - 4][n])) / DX;

				vz_dz = (
					a1*(vz[m][n] - vz[m][n - 1]) +
					a2*(vz[m][n + 1] - vz[m][n - 2]) +
					a3*(vz[m][n + 2] - vz[m][n - 3]) +
					a4*(vz[m][n + 3] - vz[m][n - 4])) / DZ;

				/*sxx��szz ��Ӧ������*/
				sxx[m][n] += DT*((lam[m][n] + 2.0 * u[m][n])*vx_dx + lam[m][n] * vz_dz);
				szz[m][n] += DT*((lam[m][n] + 2.0 * u[m][n])*vz_dz + lam[m][n] * vx_dx);

				/*�ٶȿռ�ƫ��*/
				vx_dz = (
					a1*(vx[m][n + 1] - vx[m][n]) +
					a2*(vx[m][n + 2] - vx[m][n - 1]) +
					a3*(vx[m][n + 3] - vx[m][n - 2]) +
					a4*(vx[m][n + 4] - vx[m][n - 3])) / DZ;

				vz_dx = (
					a1*(vz[m + 1][n] - vz[m][n]) +
					a2*(vz[m + 2][n] - vz[m - 1][n]) +
					a3*(vz[m + 3][n] - vz[m - 2][n]) +
					a4*(vz[m + 4][n] - vz[m - 3][n])) / DX;

				/*sxz ��Ӧ������*/
				sxz[m][n] += DT*eff_u_xz*(vx_dz + vz_dx);
				break;
			case 1:/*H-point*/
				vx_dx =
					(
					a1*(vx[m][n] - vx[m - 1][n]) +
					a2*(vx[m + 1][n] - vx[m - 2][n]) +
					a3*(vx[m + 2][n] - vx[m - 3][n]) +
					a4*(vx[m + 3][n] - vx[m - 4][n])) / DX;
				/*Dvz��Dvx����*/
				vz_dz = (-lam[m][n] / (lam[m][n] + 2.0*u[m][n]))*vx_dx;

				/*����sxx*/
				sxx[m][n] += DT*((lam[m][n] + 2.0 * u[m][n])*vx_dx + lam[m][n] * vz_dz);

				/*�ٶȿռ�ƫ��*/
				vx_dz = (
					a1*(vx[m][n + 1] - vx[m][n]) +
					a2*(vx[m][n + 2] - vx[m][n - 1]) +
					a3*(vx[m][n + 3] - vx[m][n - 2]) +
					a4*(vx[m][n + 4] - vx[m][n - 3])) / DZ;

				vz_dx = (
					a1*(vz[m + 1][n] - vz[m][n]) +
					a2*(vz[m + 2][n] - vz[m - 1][n]) +
					a3*(vz[m + 3][n] - vz[m - 2][n]) +
					a4*(vz[m + 4][n] - vz[m - 3][n])) / DX;

				/*sxz ��Ӧ������*/
				sxz[m][n] += DT*u[m][n] * (vx_dz + vz_dx);
				break;
			case 4:/*IR-point*/
				/*�ٶȿռ�ƫ��*/
				vx_dz = (
					a1*(vx[m][n + 1] - vx[m][n]) +
					a2*(vx[m][n + 2] - vx[m][n - 1]) +
					a3*(vx[m][n + 3] - vx[m][n - 2]) +
					a4*(vx[m][n + 4] - vx[m][n - 3])) / DZ;

				vz_dx = (
					a1*(vz[m + 1][n] - vz[m][n]) +
					a2*(vz[m + 2][n] - vz[m - 1][n]) +
					a3*(vz[m + 3][n] - vz[m - 2][n]) +
					a4*(vz[m + 4][n] - vz[m - 3][n])) / DX;

				/*sxz ��Ӧ������*/
				sxz[m][n] += DT*eff_u_xz*(vx_dz + vz_dx);
				//printf("%d\t,%d   %d\n", flag[m][n], m, n);
				break;
			case 5:/*IL-point*/
				/*�ٶȿռ�ƫ��*/
				vx_dz = (
					a1*(vx[m][n + 1] - vx[m][n]) +
					a2*(vx[m][n + 2] - vx[m][n - 1]) +
					a3*(vx[m][n + 3] - vx[m][n - 2]) +
					a4*(vx[m][n + 4] - vx[m][n - 3])) / DZ;

				vz_dx = (
					a1*(vz[m + 1][n] - vz[m][n]) +
					a2*(vz[m + 2][n] - vz[m - 1][n]) +
					a3*(vz[m + 3][n] - vz[m - 2][n]) +
					a4*(vz[m + 4][n] - vz[m - 3][n])) / DX;

				/*sxz ��Ӧ������*/
				sxz[m][n] += DT*eff_u_xz*(vx_dz + vz_dx);
				break;
			case 6:/*VL-point*/
				/*�ٶȿռ�ƫ��*/
				vz_dz = (
					a1*(vz[m][n] - vz[m][n - 1]) +
					a2*(vz[m][n + 1] - vz[m][n - 2]) +
					a3*(vz[m][n + 2] - vz[m][n - 3]) +
					a4*(vz[m][n + 3] - vz[m][n - 4])) / DZ;

				vx_dx = (-lam[m][n] / (lam[m][n] + 2.0*u[m][n]))*vz_dz;

				/*sxx��szz ��Ӧ������*/
				szz[m][n] += DT*((lam[m][n] + 2.0 * u[m][n])*vz_dz + lam[m][n] * vx_dx);

				/*�ٶȿռ�ƫ��*/
				vx_dz = (
					a1*(vx[m][n + 1] - vx[m][n]) +
					a2*(vx[m][n + 2] - vx[m][n - 1]) +
					a3*(vx[m][n + 3] - vx[m][n - 2]) +
					a4*(vx[m][n + 4] - vx[m][n - 3])) / DZ;

				vz_dx = (
					a1*(vz[m + 1][n] - vz[m][n]) +
					a2*(vz[m + 2][n] - vz[m - 1][n]) +
					a3*(vz[m + 3][n] - vz[m - 2][n]) +
					a4*(vz[m + 4][n] - vz[m - 3][n])) / DX;

				/*sxz ��Ӧ������*/
				sxz[m][n] += DT*u[m][n] * (vx_dz + vz_dx);
				break;
			case 7:/*VR-point*/
				
				/*�ٶȿռ�ƫ��*/
				vz_dz = (
					a1*(vz[m][n] - vz[m][n - 1]) +
					a2*(vz[m][n + 1] - vz[m][n - 2]) +
					a3*(vz[m][n + 2] - vz[m][n - 3]) +
					a4*(vz[m][n + 3] - vz[m][n - 4])) / DZ;

				vx_dx = (-lam[m][n] / (lam[m][n] + 2.0*u[m][n]))*vz_dz;

				/*sxx��szz ��Ӧ������*/
				szz[m][n] += DT*((lam[m][n] + 2.0 * u[m][n])*vz_dz + lam[m][n] * vx_dx);

				/*�ٶȿռ�ƫ��*/
				vx_dz = (
					a1*(vx[m][n + 1] - vx[m][n]) +
					a2*(vx[m][n + 2] - vx[m][n - 1]) +
					a3*(vx[m][n + 3] - vx[m][n - 2]) +
					a4*(vx[m][n + 4] - vx[m][n - 3])) / DZ;
				vz_dx = (
					a1*(vz[m + 1][n] - vz[m][n]) +
					a2*(vz[m + 2][n] - vz[m - 1][n]) +
					a3*(vz[m + 3][n] - vz[m - 2][n]) +
					a4*(vz[m + 4][n] - vz[m - 3][n])) / DX;
				/*sxz ��Ӧ������*/
				sxz[m][n] += DT*u[m][n] * (vx_dz + vz_dx);
				//printf("%d", flag[m][n]);
				break;
			}
		}
	}
}
