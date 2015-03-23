#include"topoft.h"

void AEA_topo_update_s(int nt, int **flag, float ** vx, float ** vz,
	float ** sxx, float ** szz, float ** sxz,
	float ** lam, float ** u, pml_c pml, pml_array array, int depthDH)
{
	/*全局变量*/
	extern float DT, DX, DZ, FW;
	extern int   NX, NZ;

	/*吸收层数*/
	int ifw_x, ifw_z;

	ifw_x = iround(FW / DX);
	ifw_z = iround(FW / DZ);

	/*中间变量*/
	float vx_dx = 0.0, vz_dz = 0.0;
	float vx_dz = 0.0, vz_dx = 0.0;
	float eff_u_xz = 0.0;
	float eff_lam = 0.0;
	float tractionzz = 1.0;
	float tractionxx = 1.0;
	float tractionxz = 1.0;

	for (int m = ifw_x + 1; m <= NX - ifw_x; m++)
	{
		for (int n = 1; n < depthDH; n++)
		{
			tractionzz = 1.0;
			tractionxx = 1.0;
			tractionxz = 1.0;
			eff_lam = lam[m][n];

			if (flag[m][n] == 0)
			{
				/*有效剪切模量，调和平均*/
				eff_u_xz = 1.0 / ((1.0 / u[m][n] + 1.0 / u[m + 1][n] + 1.0 / u[m][n + 1] + 1.0 / u[m + 1][n + 1]) / 4.0);
			}
			else
			{
				eff_u_xz = u[m][n];
			}

			switch (flag[m][n])
			{
			case -1:
				break;
			case 0:/*内部介质*/
				/*速度空间偏导*/
				break;
			case 1:/*H-point,*/
				eff_lam = 0.0;
				tractionzz = 0.0;
				eff_u_xz = u[m][n];

				if (flag[m + 1][n] == 3)
				{
					tractionxz = 0.0;
				}
				break;
			case 2:/*OL-point,sxx,szz位于真空，处理方式需要检验*/
				eff_lam = 0.0;
				eff_u_xz = u[m][n];
				tractionxz = 0.0; //镜像法，设置strss xz为0
				tractionxx = 0.0; //位于真空之中
				tractionzz = 0.0;
				break;
			case 3:/*OR-point */
				eff_lam = 0.0;
				eff_u_xz = u[m][n];
				tractionxz = 0.0;
				tractionxx = 0.0;
				tractionxz = 0.0;
				break;
			case 4:/*IR-point*/
				eff_lam = 0.0;
				eff_u_xz = u[m][n];
				tractionxx = 0.0;
				tractionzz = 0.0;

				if (flag[m + 1][n] == 3)
				{
					tractionxz = 0.0;
				}
				break;
			case 5:/*IL-point,ok*/
				eff_lam = 0.0;
				eff_u_xz = u[m][n];
				tractionxx = 0.0;
				tractionzz = 0.0;
				break;
			case 6:/*VL-point,ok*/
				eff_lam = 0.0;
				tractionxx = 0.0;
				eff_u_xz = u[m][n];
				break;
			case 7:/*VR-point,ok*/
				eff_lam = 0.0;
				tractionxx = 0.0;
				tractionxz = 0.0;//位于真空之中
				eff_u_xz = u[m][n];
				break;
			}

			/*速度空间偏导*/
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

			/*sxx、szz 正应力迭代*/
			sxx[m][n] += tractionxx*DT*((eff_lam + 2.0 * u[m][n])*vx_dx + eff_lam * vz_dz);
			szz[m][n] += tractionzz*DT*((eff_lam + 2.0 * u[m][n])*vz_dz + eff_lam * vx_dx);

			/*速度空间偏导*/
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

			/*sxz 剪应力迭代*/
			sxz[m][n] += tractionxz*DT*eff_u_xz*(vx_dz + vz_dx);
		}
	}

	/*内部介质迭代*/
	for (int m = ifw_x + 1; m <= NX - ifw_x; m++)
	{
		for (int n = depthDH; n <= NZ - ifw_z; n++)
		{
			/*速度空间偏导*/
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

			/*sxx、szz 正应力迭代*/
			sxx[m][n] += tractionxx*DT*((eff_lam + 2.0 * u[m][n])*vx_dx + eff_lam * vz_dz);
			szz[m][n] += tractionzz*DT*((eff_lam + 2.0 * u[m][n])*vz_dz + eff_lam * vx_dx);

			/*速度空间偏导*/
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

			/*sxz 剪应力迭代*/
			sxz[m][n] += tractionxz*DT*eff_u_xz*(vx_dz + vz_dx);
		}
	}

	/*left PML*/
	for (int m = 1; m < 1 + ifw_x; m++)
	{
		for (int n = 1; n <= NZ - ifw_z; n++)
		{
			tractionzz = 1.0;
			eff_lam = lam[m][n];

			if (flag[m][n] == 0)
			{
				/*有效剪切模量，调和平均*/
				eff_u_xz = 1.0 / ((1.0 / u[m][n] + 1.0 / u[m + 1][n] + 1.0 / u[m][n + 1] + 1.0 / u[m + 1][n + 1]) / 4.0);
			}
			else
			{
				eff_u_xz = u[m][n];
			}

			if (flag[m][n] == 1)
			{
				eff_lam = 0.0;
				tractionzz = 0.0;
				eff_u_xz = u[m][n];
				/* if (flag[m + 1][n] == 3)
				{
					tractionxz = 0.0;
				} */
			}

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
			array.left_v_sxx[m][n] = array.left_v_sxx[m][n] * (1.0 - pml.left[m] * DT / 2.0) / (1.0 + pml.left[m] * DT / 2.0) +
				DT*(eff_lam + 2.0 * u[m][n])*vx_dx / (1.0 + pml.left[m] * DT / 2.0);
			array.left_p_sxx[m][n] += DT*eff_lam * vz_dz;

			sxx[m][n] = array.left_v_sxx[m][n] + array.left_p_sxx[m][n];

			array.left_v_szz[m][n] = array.left_v_szz[m][n] * (1.0 - pml.left[m] * DT / 2.0) / (1.0 + pml.left[m] * DT / 2.0) +
				tractionzz*DT*eff_lam*vx_dx / (1.0 + pml.left[m] * DT / 2.0);
			array.left_p_szz[m][n] += tractionzz*DT*(eff_lam + 2.0 * u[m][n])*vz_dz;

			szz[m][n] = array.left_v_szz[m][n] + array.left_p_szz[m][n];

			/*速度空间偏导*/
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
			array.left_v_sxz[m][n] = array.left_v_sxz[m][n] * (1.0 - pml.left[m] * DT / 2.0) / (1.0 + pml.left[m] * DT / 2.0) +
				DT*eff_u_xz*vz_dx / (1.0 + pml.left[m] * DT / 2.0);
			array.left_p_sxz[m][n] += DT*eff_u_xz*vx_dz;

			sxz[m][n] = array.left_v_sxz[m][n] + array.left_p_sxz[m][n];
		}
	}

	/*right PML*/
	for (int m = NX - ifw_x + 1; m <= NX; m++)
	{
		for (int n = 1; n <= NZ - ifw_z; n++)
		{
			tractionzz = 1.0;
			eff_lam = lam[m][n];

			if (flag[m][n] == 0)
			{
				/*有效剪切模量，调和平均*/
				eff_u_xz = 1.0 / ((1.0 / u[m][n] + 1.0 / u[m + 1][n] + 1.0 / u[m][n + 1] + 1.0 / u[m + 1][n + 1]) / 4.0);
			}
			else
			{
				eff_u_xz = u[m][n];
			}

			if (flag[m][n] == 1)
			{
				eff_lam = 0.0;
				tractionzz = 0.0;
				eff_u_xz = u[m][n];
				/* if (m < NX&&flag[m + 1][n] == 3)
				{
					tractionxz = 0.0;
				} */
			}

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
			array.right_v_sxx[m][n] = array.right_v_sxx[m][n] * (1.0 - pml.right[m] * DT / 2.0) / (1.0 + pml.right[m] * DT / 2.0) +
				DT*(eff_lam + 2.0 * u[m][n])*vx_dx / (1.0 + pml.right[m] * DT / 2.0);
			array.right_p_sxx[m][n] += DT*eff_lam * vz_dz;

			sxx[m][n] = array.right_v_sxx[m][n] + array.right_p_sxx[m][n];

			array.right_v_szz[m][n] = array.right_v_szz[m][n] * (1.0 - pml.right[m] * DT / 2.0) / (1.0 + pml.right[m] * DT / 2.0) +
				tractionzz*DT*eff_lam*vx_dx / (1.0 + pml.right[m] * DT / 2.0);
			array.right_p_szz[m][n] += tractionzz*DT*(eff_lam + 2.0 * u[m][n])*vz_dz;

			szz[m][n] = array.right_v_szz[m][n] + array.right_p_szz[m][n];

			/*速度空间偏导*/
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
			array.right_v_sxz[m][n] = array.right_v_sxz[m][n] * (1.0 - pml.right[m] * DT / 2.0) / (1.0 + pml.right[m] * DT / 2.0) +
				DT*eff_u_xz*vz_dx / (1.0 + pml.right[m] * DT / 2.0);
			array.right_p_sxz[m][n] += DT*eff_u_xz*vx_dz;

			sxz[m][n] = array.right_v_sxz[m][n] + array.right_p_sxz[m][n];
		}
	}

	/*bottom PML*/
	for (int m = 1 + ifw_x; m <= NX - ifw_x; m++)
	{
		for (int n = NZ - ifw_z + 1; n <= NZ; n++)
		{
			eff_lam = lam[m][n];
			if (flag[m][n] == 0)
			{
				/*有效剪切模量，调和平均*/
				eff_u_xz = 1.0 / ((1.0 / u[m][n] + 1.0 / u[m + 1][n] + 1.0 / u[m][n + 1] + 1.0 / u[m + 1][n + 1]) / 4.0);
			}
			else
			{
				eff_u_xz = u[m][n];

			}

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
			array.bottom_v_sxx[m][n] = array.bottom_v_sxx[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT*eff_lam * vz_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_sxx[m][n] += DT*(eff_lam + 2.0 * u[m][n])*vx_dx;

			sxx[m][n] = array.bottom_v_sxx[m][n] + array.bottom_p_sxx[m][n];

			array.bottom_v_szz[m][n] = array.bottom_v_szz[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT*(eff_lam + 2.0 * u[m][n])*vz_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_szz[m][n] += DT*eff_lam*vx_dx;

			szz[m][n] = array.bottom_v_szz[m][n] + array.bottom_p_szz[m][n];

			/*速度空间偏导*/
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
			array.bottom_v_sxz[m][n] = array.bottom_v_sxz[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT* eff_u_xz*vx_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_sxz[m][n] += DT*eff_u_xz*vz_dx;

			sxz[m][n] = array.bottom_v_sxz[m][n] + array.bottom_p_sxz[m][n];
		}
	}

	/*Bottom left corner*/
	for (int m = 1; m < 1 + ifw_x; m++)
	{
		for (int n = NZ - ifw_z + 1; n <= NZ; n++)
		{
			eff_lam = lam[m][n];
			if (flag[m][n] == 0)
			{
				/*有效剪切模量，调和平均*/
				eff_u_xz = 1.0 / ((1.0 / u[m][n] + 1.0 / u[m + 1][n] + 1.0 / u[m][n + 1] + 1.0 / u[m + 1][n + 1]) / 4.0);
			}
			else
			{
				eff_u_xz = u[m][n];
			}

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
			array.bottom_v_sxx[m][n] = array.bottom_v_sxx[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT*eff_lam * vz_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_sxx[m][n] = array.bottom_p_sxx[m][n] * (1.0 - pml.left[m] * DT / 2.0) / (1.0 + pml.left[m] * DT / 2.0) +
				DT*(eff_lam + 2.0 * u[m][n])*vx_dx / (1.0 + pml.left[m] * DT / 2.0);
			sxx[m][n] = array.bottom_v_sxx[m][n] + array.bottom_p_sxx[m][n];

			array.bottom_v_szz[m][n] = array.bottom_v_szz[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT*(eff_lam + 2.0 * u[m][n])*vz_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_szz[m][n] = array.bottom_p_szz[m][n] * (1.0 - pml.left[m] * DT / 2.0) / (1.0 + pml.left[m] * DT / 2.0) +
				DT*eff_lam*vx_dx / (1.0 + pml.left[m] * DT / 2.0);

			szz[m][n] = array.bottom_v_szz[m][n] + array.bottom_p_szz[m][n];

			/*速度空间偏导*/
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
			array.bottom_v_sxz[m][n] = array.bottom_v_sxz[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT* eff_u_xz*vx_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_sxz[m][n] = array.bottom_p_sxz[m][n] * (1.0 - pml.left[m] * DT / 2.0) / (1.0 + pml.left[m] * DT / 2.0) +
				DT*eff_u_xz*vz_dx / (1.0 + pml.left[m] * DT / 2.0);

			sxz[m][n] = array.bottom_v_sxz[m][n] + array.bottom_p_sxz[m][n];
		}
	}

	/*Bottom right corner*/
	for (int m = NX - ifw_x + 1; m <= NX; m++)
	{
		for (int n = NZ - ifw_z + 1; n <= NZ; n++)
		{
			eff_lam = lam[m][n];
			if (flag[m][n] == 0)
			{
				/*有效剪切模量，调和平均*/
				eff_u_xz = 1.0 / ((1.0 / u[m][n] + 1.0 / u[m + 1][n] + 1.0 / u[m][n + 1] + 1.0 / u[m + 1][n + 1]) / 4.0);
			}
			else
			{
				eff_u_xz = u[m][n];
			}

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
			array.bottom_v_sxx[m][n] = array.bottom_v_sxx[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT*eff_lam * vz_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_sxx[m][n] = array.bottom_p_sxx[m][n] * (1.0 - pml.right[m] * DT / 2.0) / (1.0 + pml.right[m] * DT / 2.0) +
				DT*(eff_lam + 2.0 * u[m][n])*vx_dx / (1.0 + pml.right[m] * DT / 2.0);
			sxx[m][n] = array.bottom_v_sxx[m][n] + array.bottom_p_sxx[m][n];

			array.bottom_v_szz[m][n] = array.bottom_v_szz[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT*(eff_lam + 2.0 * u[m][n])*vz_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_szz[m][n] = array.bottom_p_szz[m][n] * (1.0 - pml.right[m] * DT / 2.0) / (1.0 + pml.right[m] * DT / 2.0) +
				DT*eff_lam*vx_dx / (1.0 + pml.right[m] * DT / 2.0);

			szz[m][n] = array.bottom_v_szz[m][n] + array.bottom_p_szz[m][n];

			/*速度空间偏导*/
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
			array.bottom_v_sxz[m][n] = array.bottom_v_sxz[m][n] * (1.0 - pml.bottom[n] * DT / 2.0) / (1.0 + pml.bottom[n] * DT / 2.0) +
				DT* eff_u_xz*vx_dz / (1.0 + pml.bottom[n] * DT / 2.0);
			array.bottom_p_sxz[m][n] = array.bottom_p_sxz[m][n] * (1.0 - pml.right[m] * DT / 2.0) / (1.0 + pml.right[m] * DT / 2.0) +
				DT*eff_u_xz*vz_dx / (1.0 + pml.right[m] * DT / 2.0);

			sxz[m][n] = array.bottom_v_sxz[m][n] + array.bottom_p_sxz[m][n];
		}
	}
}