#include"topoft.h"

void AEA_topo_update_v(int nt, int **flag, float ** vx, float ** vz,
	float ** sxx, float ** szz, float ** sxz,
	float **  rden, pml_c pml, pml_array array, int depthDH)
{
	extern float DT, DX, DZ, FW;
	extern int   NX, NZ;
	
	/*吸收层数*/
	int ifw_x, ifw_z;

	ifw_x = iround(FW / DX);
	ifw_z = iround(FW / DZ);

	/**避免角落网格重叠，对水平方向和垂直方向求导分开进行**/
	float sxx_dx = 0.0, sxz_dz = 0.0;
	float sxz_dx = 0.0, szz_dz = 0.0;
	float eff_rden_x = 0.0, eff_rden_z = 0.0;

	/*速度迭代*/
	for (int m = ifw_x + 1; m <= NX - ifw_x; m++)
	{
		for (int n = 1; n < depthDH; n++)
		{
			/*有效介质参数*/
			if (flag[m][n] == 0)
			{
				eff_rden_x = (rden[m][n] + rden[m + 1][n]) / 2.0;
				eff_rden_z = (rden[m][n] + rden[m][n + 1]) / 2.0;
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
				break;
			case 1:/*H-point,*/
				eff_rden_x = rden[m][n] * 2.0;
				eff_rden_z = rden[m][n];
				
				if (flag[m + 1][n] == 3)
				{
					eff_rden_z = rden[m][n] * 2.0;
				}
				break;
			case 2:/*OL-point*/
				eff_rden_x = rden[m][n] * 2.0;
				eff_rden_z = rden[m][n] * 2.0;
				/*保证稳定性，将OL与IL相邻的点速度设置为0*/

                if (flag[m][n + 1] == 5)
				{
					eff_rden_z = 0.0;
				}
                if (flag[m + 1][n] == 5)
				{
					eff_rden_x = 0.0;
				}
				break;
			case 3:/*OR-point */
				eff_rden_z = rden[m][n] * 2.0;
				//eff_rden_x = rden[m][n] * 2.0;
				eff_rden_x = 0; //位于
				
                if (flag[m][n + 1] == 4)
				{
					eff_rden_z = 0.0;
                }
                if (flag[m - 1][n] == 4)
                {
                    eff_rden_z = 0.0;
                }
				break;
			case 4:/*IR-point*/
				eff_rden_x = rden[m][n] * 2.0;
                eff_rden_z = rden[m][n];
				
				if (flag[m + 1][n] == 3)
				{
					eff_rden_z = rden[m][n] * 2.0;
				}
				break;
			case 5:/*IL-point*/
				eff_rden_x = rden[m][n];
				eff_rden_z = rden[m][n];
				break;
			case 6:/*VL-point*/
				eff_rden_x = rden[m][n];
				eff_rden_z = rden[m][n] * 2.0;
				break;
			case 7:/*VR-point*/
				eff_rden_x = rden[m][n];
				eff_rden_z = rden[m][n] * 2.0;
				break;
			}

			/*应力空间偏导*/
			sxx_dx = (
				a1*(sxx[m + 1][n] - sxx[m][n]) +
				a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
				a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
				a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

			sxz_dz = (
				a1*(sxz[m][n] - sxz[m][n - 1]) +
				a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
				a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
				a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

			/*vx 速度迭代*/
			vx[m][n] += DT * eff_rden_x * (sxx_dx + sxz_dz);
			//
			/*应力空间偏导*/
			sxz_dx = (
				a1*(sxz[m][n] - sxz[m - 1][n]) +
				a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
				a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
				a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

			szz_dz = (
				a1*(szz[m][n + 1] - szz[m][n]) +
				a2*(szz[m][n + 2] - szz[m][n - 1]) +
				a3*(szz[m][n + 3] - szz[m][n - 2]) +
				a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

			/*vz 速度迭代*/
			vz[m][n] += DT * eff_rden_z * (sxz_dx + szz_dz);
		}
	}

	for (int m = ifw_x + 1; m <= NX - ifw_x; m++)
	{
		for (int n = depthDH; n <= NZ - ifw_z; n++)
		{
			/*应力空间偏导*/
			sxx_dx = (
				a1*(sxx[m + 1][n] - sxx[m][n]) +
				a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
				a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
				a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

			sxz_dz = (
				a1*(sxz[m][n] - sxz[m][n - 1]) +
				a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
				a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
				a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

			/*vx 速度迭代*/
			vx[m][n] += DT * eff_rden_x * (sxx_dx + sxz_dz);
			//
			/*应力空间偏导*/
			sxz_dx = (
				a1*(sxz[m][n] - sxz[m - 1][n]) +
				a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
				a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
				a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

			szz_dz = (
				a1*(szz[m][n + 1] - szz[m][n]) +
				a2*(szz[m][n + 2] - szz[m][n - 1]) +
				a3*(szz[m][n + 3] - szz[m][n - 2]) +
				a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

			/*vz 速度迭代*/
			vz[m][n] += DT * eff_rden_z * (sxz_dx + szz_dz);
		}
	}

	/*left PML*/
	for (int m = 1; m < 1 + ifw_x; m++)
	{
		for (int n = 1; n <= NZ - ifw_z; n++)
		{
			/******垂向求导**********/
			/*有效介质参数*/
			if (flag[m][n] == 0)
			{
				eff_rden_x = (rden[m][n] + rden[m + 1][n]) / 2.0;
				eff_rden_z = (rden[m][n] + rden[m][n + 1]) / 2.0;
			}
			else
			{
				eff_rden_x = eff_rden_z = rden[m][n];
			}

			//水平自由地表延伸近PML层
			if (flag[m][n] == 1)
			{
				eff_rden_x = rden[m][n] * 2.0;
				eff_rden_z = rden[m][n];
			}

			sxx_dx = (
				a1*(sxx[m + 1][n] - sxx[m][n]) +
				a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
				a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
				a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;
			array.left_v_vx[m][n] = array.left_v_vx[m][n] * (1.0 - pml.leftv[m] * DT / 2.0) / (1.0 + pml.leftv[m] * DT / 2.0) +
				(DT*eff_rden_x*sxx_dx) / (1.0 + pml.leftv[m] * DT / 2.0);

			sxz_dz = (
				a1*(sxz[m][n] - sxz[m][n - 1]) +
				a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
				a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
				a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;
			array.left_p_vx[m][n] += DT*eff_rden_x*sxz_dz;

			vx[m][n] = array.left_v_vx[m][n] + array.left_p_vx[m][n];

			sxz_dx = (
				a1*(sxz[m][n] - sxz[m - 1][n]) +
				a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
				a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
				a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;
			array.left_v_vz[m][n] = array.left_v_vz[m][n] * (1.0 - pml.leftv[m] * DT / 2.0) / (1.0 + pml.leftv[m] * DT / 2.0) +
				(DT * eff_rden_z*sxz_dx) / (1.0 + pml.leftv[m] * DT / 2.0);

			szz_dz = (
				a1*(szz[m][n + 1] - szz[m][n]) +
				a2*(szz[m][n + 2] - szz[m][n - 1]) +
				a3*(szz[m][n + 3] - szz[m][n - 2]) +
				a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;
			array.left_p_vz[m][n] += DT*eff_rden_z*szz_dz;
			
			vz[m][n] = array.left_v_vz[m][n] + array.left_p_vz[m][n];
		}
	}

	/*right PML*/
	for (int m = NX - ifw_x + 1; m <= NX; m++)
	{
		for (int n = 1; n <= NZ - ifw_z; n++)
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

			if (flag[m][n] == 1)
			{
				eff_rden_x = rden[m][n] * 2.0;
				eff_rden_z = rden[m][n];
			}

			sxx_dx = (
				a1*(sxx[m + 1][n] - sxx[m][n]) +
				a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
				a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
				a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;
			array.right_v_vx[m][n] = array.right_v_vx[m][n] * (1.0 - pml.rightv[m] * DT / 2.0) / (1 + pml.rightv[m] * DT / 2.0) +
				(DT*eff_rden_x*sxx_dx) / (1.0 + pml.rightv[m] * DT / 2.0);

			sxz_dz = (
				a1*(sxz[m][n] - sxz[m][n - 1]) +
				a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
				a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
				a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;
			array.right_p_vx[m][n] += DT*eff_rden_x*sxz_dz;

			vx[m][n] = array.right_v_vx[m][n] + array.right_p_vx[m][n];

			sxz_dx = (
				a1*(sxz[m][n] - sxz[m - 1][n]) +
				a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
				a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
				a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;
			array.right_v_vz[m][n] = array.right_v_vz[m][n] * (1.0 - pml.rightv[m] * DT / 2.0) / (1.0 + pml.rightv[m] * DT / 2.0) +
				(DT * eff_rden_z*sxz_dx) / (1.0 + pml.rightv[m] * DT / 2.0);

			szz_dz = (
				a1*(szz[m][n + 1] - szz[m][n]) +
				a2*(szz[m][n + 2] - szz[m][n - 1]) +
				a3*(szz[m][n + 3] - szz[m][n - 2]) +
				a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;
			array.right_p_vz[m][n] += DT*eff_rden_z*szz_dz;

			vz[m][n] = array.right_v_vz[m][n] + array.right_p_vz[m][n];
		}
	}

	/*bottom PML*/
	for (int m = 1 + ifw_x; m <= NX - ifw_x;m++)
	{
		for (int n = NZ - ifw_z + 1; n <= NZ; n++)
		{
			if (flag[m][n] == 0)
			{
				eff_rden_x = (rden[m][n] + rden[m + 1][n]) / 2.0;
				eff_rden_z = (rden[m][n] + rden[m][n + 1]) / 2.0;
			}
			else
			{
				eff_rden_x = eff_rden_z = rden[m][n];
			}

			sxz_dz = (
				a1*(sxz[m][n] - sxz[m][n - 1]) +
				a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
				a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
				a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;
			array.bottom_v_vx[m][n] = array.bottom_v_vx[m][n] * (1.0 - pml.bottomv[n] * DT / 2.0) / (1.0 + pml.bottomv[n] * DT / 2.0) +
				(DT*eff_rden_x*sxz_dz) / (1.0 + pml.bottomv[n] * DT / 2.0);

			sxx_dx = (
				a1*(sxx[m + 1][n] - sxx[m][n]) +
				a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
				a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
				a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;
			array.bottom_p_vx[m][n] += DT*eff_rden_x*sxx_dx;

			vx[m][n] = array.bottom_v_vx[m][n] + array.bottom_p_vx[m][n];

			szz_dz = (
				a1*(szz[m][n + 1] - szz[m][n]) +
				a2*(szz[m][n + 2] - szz[m][n - 1]) +
				a3*(szz[m][n + 3] - szz[m][n - 2]) +
				a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;
			array.bottom_v_vz[m][n] = array.bottom_v_vz[m][n] * (1.0 - pml.bottomv[n] * DT / 2.0) / (1.0 + pml.bottomv[n] * DT / 2.0) +
				(DT * eff_rden_z*szz_dz) / (1.0 + pml.bottomv[n] * DT / 2.0);

			sxz_dx = (
				a1*(sxz[m][n] - sxz[m - 1][n]) +
				a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
				a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
				a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;
			
			array.bottom_p_vz[m][n] += DT*eff_rden_z*sxz_dx;

			vz[m][n] = array.bottom_v_vz[m][n] + array.bottom_p_vz[m][n];
		}
	}

		/*Bottom left corner*/
	for (int m = 1; m < 1 + ifw_x; m++)
	{
		for (int n = NZ - ifw_z + 1; n <= NZ; n++)
		{
			if (flag[m][n] == 0)
			{
				eff_rden_x = (rden[m][n] + rden[m + 1][n]) / 2.0;
				eff_rden_z = (rden[m][n] + rden[m][n + 1]) / 2.0;
			}
			else
			{
				eff_rden_x = eff_rden_z = rden[m][n];
			}
			
			sxz_dz = (
				a1*(sxz[m][n] - sxz[m][n - 1]) +
				a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
				a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
				a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;
			array.bottom_v_vx[m][n] = array.bottom_v_vx[m][n] * (1.0 - pml.bottomv[n] * DT / 2.0) / (1.0 + pml.bottomv[n] * DT / 2.0) +
				(DT*eff_rden_x*sxz_dz) / (1.0 + pml.bottomv[n] * DT / 2.0);

			sxx_dx = (
				a1*(sxx[m + 1][n] - sxx[m][n]) +
				a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
				a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
				a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;
			array.bottom_p_vx[m][n] = array.bottom_p_vx[m][n] * (1.0 - pml.leftv[m] * DT / 2.0) / (1.0 + pml.leftv[m] * DT / 2.0) +
				(DT*eff_rden_x*sxx_dx) / (1.0 + pml.leftv[m] * DT / 2.0);

			vx[m][n] = array.bottom_v_vx[m][n] + array.bottom_p_vx[m][n];

			szz_dz = (
				a1*(szz[m][n + 1] - szz[m][n]) +
				a2*(szz[m][n + 2] - szz[m][n - 1]) +
				a3*(szz[m][n + 3] - szz[m][n - 2]) +
				a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;
			array.bottom_v_vz[m][n] = array.bottom_v_vz[m][n] * (1.0 - pml.bottomv[n] * DT / 2.0) / (1 + pml.bottomv[n] * DT / 2.0) +
				(DT * eff_rden_z*szz_dz) / (1.0 + pml.bottomv[n] * DT / 2.0);

			sxz_dx = (
				a1*(sxz[m][n] - sxz[m - 1][n]) +
				a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
				a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
				a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

			array.bottom_p_vz[m][n] = array.bottom_p_vz[m][n] * (1.0 - pml.leftv[m] * DT / 2.0) / (1.0 + pml.leftv[m] * DT / 2.0) +
				(DT*eff_rden_z*sxz_dx) / (1.0 + pml.leftv[m] * DT / 2.0);

			vz[m][n] = array.bottom_v_vz[m][n] + array.bottom_p_vz[m][n];
		}
	}

	/*Bottom right corner*/
	for (int m = NX - ifw_x + 1; m <= NX; m++)
	{
		for (int n = NZ - ifw_z + 1; n <= NZ; n++)
		{
			if (flag[m][n] == 0)
			{
				eff_rden_x = (rden[m][n] + rden[m + 1][n]) / 2.0;
				eff_rden_z = (rden[m][n] + rden[m][n + 1]) / 2.0;
			}
			else
			{
				eff_rden_x = eff_rden_z = rden[m][n];
			}
			sxz_dz = (
				a1*(sxz[m][n] - sxz[m][n - 1]) +
				a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
				a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
				a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;
			array.bottom_v_vx[m][n] = array.bottom_v_vx[m][n] * (1.0 - pml.bottomv[n] * DT / 2.0) / (1.0 + pml.bottomv[n] * DT / 2.0) +
				(DT*eff_rden_x*sxz_dz) / (1.0 + pml.bottomv[n] * DT / 2.0);

			sxx_dx = (
				a1*(sxx[m + 1][n] - sxx[m][n]) +
				a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
				a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
				a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;
			array.bottom_p_vx[m][n] = array.bottom_p_vx[m][n] * (1.0 - pml.rightv[m] * DT / 2.0) / (1.0 + pml.rightv[m] * DT / 2.0) +
				(DT*eff_rden_x*sxx_dx) / (1.0 + pml.rightv[m] * DT / 2.0);

			vx[m][n] = array.bottom_v_vx[m][n] + array.bottom_p_vx[m][n];

			szz_dz = (
				a1*(szz[m][n + 1] - szz[m][n]) +
				a2*(szz[m][n + 2] - szz[m][n - 1]) +
				a3*(szz[m][n + 3] - szz[m][n - 2]) +
				a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;
			array.bottom_v_vz[m][n] = array.bottom_v_vz[m][n] * (1.0 - pml.bottomv[n] * DT / 2.0) / (1.0 + pml.bottomv[n] * DT / 2.0) +
				(DT * eff_rden_z*szz_dz) / (1.0 + pml.bottomv[n] * DT / 2.0);

			sxz_dx = (
				a1*(sxz[m][n] - sxz[m - 1][n]) +
				a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
				a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
				a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

			array.bottom_p_vz[m][n] = array.bottom_p_vz[m][n] * (1.0 - pml.rightv[m] * DT / 2.0) / (1.0 + pml.rightv[m] * DT / 2.0) +
				(DT*eff_rden_z*sxz_dx) / (1.0 + pml.rightv[m] * DT / 2.0);

			vz[m][n] = array.bottom_v_vz[m][n] + array.bottom_p_vz[m][n];
		}
	}

}
