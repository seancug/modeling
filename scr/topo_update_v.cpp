#include"topoft.h"

void topo_update_v(int nt, int **flag, float ** vx, float ** vz,
	float ** sxx, float ** szz, float ** sxz,
	float **  rden, float ** absorb_coeff)
{
	extern float DT, DX, DZ, FW;
	extern int   NX, NZ;

	/**避免角落网格重叠，对水平方向和垂直方向求导分开进行**/
	float sxx_dx = 0, sxz_dz = 0;
	float sxz_dx = 0, szz_dz = 0;
	float eff_rden_x = 0, eff_rden_z = 0;

	/*先对Z方向进行求导*/
	for (int m = 1; m <= NX; m++)
	{
		for (int n = 1; n <= NZ; n++)
		{
			/******垂向求导**********/
			/*有效介质参数*/
			if (flag[m][n] == 0)
			{
				eff_rden_x = (rden[m][n] + rden[m + 1][n]) / 2;
				eff_rden_z = (rden[m][n] + rden[m][n + 1]) / 2;
			}
			else
			{
				eff_rden_x = eff_rden_z = rden[m][n];
			}

			switch (flag[m][n])
			{
			case -1:
				break;
			case 0:
				/*z方向应力空间偏导*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z方向应力空间偏导*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 1:/*H-point,正常镜像法,Check*/
				szz[m][n] = 0;
				szz[m][n - 1] = -szz[m][n + 1];
				szz[m][n - 2] = -szz[m][n + 2];
				szz[m][n - 3] = -szz[m][n + 3];

				if (m < NX&&flag[m + 1][n] != 3)
				{
					sxz[m][n - 1] = -sxz[m][n];
					sxz[m][n - 2] = -sxz[m][n + 1];
					sxz[m][n - 3] = -sxz[m][n + 2];
					sxz[m][n - 4] = -sxz[m][n + 3];
				}
				else if (m < NX&&flag[m + 1][n] == 3)
				{
					sxz[m][n] = 0;
				}

				/*z方向应力空间偏导*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z方向应力空间偏导*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 2:/*OL-point,sxz在自由表面处设为零,szz,sxz垂直镜像*/
				sxz[m][n] = 0;
				sxz[m][n - 1] = -sxz[m][n + 1];
				sxz[m][n - 2] = -sxz[m][n + 2];
				sxz[m][n - 3] = -sxz[m][n + 3];
				sxz[m][n - 4] = -sxz[m][n + 4];

				szz[m][n] = -szz[m][n + 1];
				szz[m][n - 1] = -szz[m][n + 2];
				szz[m][n - 2] = -szz[m][n + 3];
				szz[m][n - 3] = -szz[m][n + 4];

				if (flag[m + 1][n] == 5 || flag[m][n + 1] == 5)
				{
					eff_rden_x = 0;
					eff_rden_z = 0;
				}

				/*z方向应力空间偏导*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z方向应力空间偏导*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 3:/*OR-point sxz自由表面出设为零，szz,sxz垂直镜像*/
				sxz[m - 1][n] = 0;

				szz[m][n] = -szz[m][n + 1];
				szz[m][n - 1] = -szz[m][n + 2];
				szz[m][n - 2] = -szz[m][n + 3];
				szz[m][n - 3] = -szz[m][n + 4];

				eff_rden_x = 0;
				if (flag[m - 1][n] == 4 || flag[m][n + 1] == 4)
				{
					eff_rden_z = 0;
				}
				break;

				/*z方向应力空间偏导*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z方向应力空间偏导*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
			case 4:/*IR-point sxx,szz位于自由边界上,置为0,sxz垂直镜像*/
				sxx[m][n] = 0;
				szz[m][n] = 0;

				if (m < NX&&flag[m + 1][n] != 3)
				{
					sxz[m][n - 1] = -sxz[m][n];
					sxz[m][n - 2] = -sxz[m][n + 1];
					sxz[m][n - 3] = -sxz[m][n + 2];
					sxz[m][n - 4] = -sxz[m][n + 3];	
				}
				else if (m < NX&&flag[m + 1][n] == 3)
				{
					sxz[m][n] = 0;
				}

				if (flag[m + 1][n] == 3 || flag[m][n + 1] == 3)
				{
					eff_rden_x = 0;
					eff_rden_z = 0;
				}
				
				/*z方向应力空间偏导*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z方向应力空间偏导*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 5:/*IL-point 无镜像值,sxx,szz位于自由边界上,置为0*/
				sxx[m][n] = 0;
				szz[m][n] = 0;
				
				/*z方向应力空间偏导*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z方向应力空间偏导*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 6:/*VL-point*/
				/*z方向应力空间偏导*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z方向应力空间偏导*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 7:/*VR-point*/
				/*z方向应力空间偏导*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z方向应力空间偏导*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			}
		}
	}

	/*先对x方向进行求导*/
	for (int m = 1; m <= NX; m++)
	{
		for (int n = 1; n <= NZ; n++)
		{
			/*有效介质参数*/
			if (flag[m][n] == 0)
			{
				eff_rden_x = (rden[m][n] + rden[m + 1][n]) / 2;
				eff_rden_z = (rden[m][n] + rden[m][n + 1]) / 2;
			}
			else
			{
				eff_rden_x = eff_rden_z = rden[m][n];
			}

			switch (flag[m][n])
			{
			case -1:/*空气*/
				break;
			case 0:/*内部介质*/
				/******水平求导**********/
				/*应力空间偏导*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*应力空间偏导*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 1:/*H-point*/
				/******水平求导**********/
				/*应力空间偏导*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*应力空间偏导*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 2:/*OL-point,sxz在自由表面处设为零,szz,sxz垂直镜像,前面已经设置过*/
				sxz[m - 1][n] = -sxz[m + 1][n];
				sxz[m - 2][n] = -sxz[m + 2][n];
				sxz[m - 3][n] = -sxz[m + 3][n];
				sxz[m - 4][n] = -sxz[m + 4][n];

				sxx[m][n] = -sxx[m + 1][n];
				sxx[m - 1][n] = -sxx[m + 2][n];
				sxx[m - 2][n] = -sxx[m + 3][n];
				sxx[m - 3][n] = -sxx[m + 4][n];
				if (flag[m + 1][n] == 5 || flag[m][n + 1] == 5)
				{
					eff_rden_x = 0;
					eff_rden_z = 0;
				}
				
				/******水平求导**********/
				/*应力空间偏导*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*应力空间偏导*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 3:/*OR-point sxz自由表面出设为零,szz,sxz垂直镜像*/
				sxz[m][n] = -sxz[m - 2][n];
				sxz[m + 1][n] = -sxz[m - 3][n];
				sxz[m + 2][n] = -sxz[m - 4][n];
				sxz[m + 3][n] = -sxz[m - 5][n];

				sxx[m][n] = -sxx[m - 1][n];
				sxx[m + 1][n] = -sxx[m - 2][n];
				sxx[m + 2][n] = -sxx[m - 3][n];
				sxx[m + 3][n] = -sxx[m - 4][n];

				eff_rden_x = 0;
				if (flag[m - 1][n] == 4 || flag[m][n + 1] == 4)
				{
					eff_rden_z = 0;
				}
				
				/******水平求导**********/
				/*应力空间偏导*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*应力空间偏导*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 4:/*IR-point*/
				/******水平求导**********/
				/*应力空间偏导*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				if (flag[m + 1][n] == 3 || flag[m][n + 1] == 3)
				{
					eff_rden_x = 0;
					eff_rden_z = 0;
				}

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*应力空间偏导*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 5:/*IL-point*/
								/******水平求导**********/
				/*应力空间偏导*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*应力空间偏导*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 6:/*VL-point 水平镜像*/
				sxx[m][n] = 0;
				sxx[m - 1][n] = -sxx[m + 1][n];
				sxx[m - 2][n] = -sxx[m + 2][n];
				sxx[m - 3][n] = -sxx[m + 3][n];

				sxz[m - 1][n] = -sxz[m][n];
				sxz[m - 2][n] = -sxz[m + 1][n];
				sxz[m - 3][n] = -sxz[m + 2][n];
				sxz[m - 4][n] = -sxz[m + 3][n];

				/******水平求导**********/
				/*应力空间偏导*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*应力空间偏导*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 7:/*VR-point 水平镜像*/
				sxx[m][n] = 0;
				sxx[m + 1][n] = -sxx[m - 1][n];
				sxx[m + 2][n] = -sxx[m - 2][n];
				sxx[m + 3][n] = -sxx[m - 3][n];
				sxx[m + 4][n] = -sxx[m - 4][n];

				sxz[m][n] = -sxz[m - 1][n];
				sxz[m + 1][n] = -sxz[m - 2][n];
				sxz[m + 2][n] = -sxz[m - 3][n];
				sxz[m + 3][n] = -sxz[m - 4][n];	

				/******水平求导**********/
				/*应力空间偏导*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx 速度迭代*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*应力空间偏导*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz 速度迭代*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			}
		}
	}

	/*吸收边界条件添加,效果不是很好*/
	if (FW > 0)
	{
		for (int nr = 1; nr <= NX; nr++)
		{
			for (int nc = 1; nc <= NZ; nc++)
			{
				vx[nr][nc] *= absorb_coeff[nr][nc];
				vz[nr][nc] *= absorb_coeff[nr][nc];
				sxx[nr][nc] *= absorb_coeff[nr][nc];
				sxz[nr][nc] *= absorb_coeff[nr][nc];
				szz[nr][nc] *= absorb_coeff[nr][nc];
			}
		}
	}
}