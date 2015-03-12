#include"topoft.h"

void AEA_planer_update_v(int nt, float ** vx, float ** vz,
	float ** sxx, float ** szz, float ** sxz,
	float **  rden, float ** absorb_coeff)
{
	extern float DT, DX, DZ, FW;
	extern int   NX, NZ;

	//int nxs, nzs;
	///*震源位置*/
	//nxs = iround(XS / DX);
	//nzs = iround(ZS / DZ);

	/*中间应力导数变量以及有效参数*/
	float sxx_dx = 0.0, sxz_dz = 0.0;
	float sxz_dx = 0.0, szz_dz = 0.0;
	float eff_rden_x = 0.0, eff_rden_z = 0.0;
	for (int m = 1; m <= NX; m++)
	{
		for (int n = 1; n <= NZ; n++)
		{
			if (n == 1)
			{
				eff_rden_x = 2.0*rden[m][n];
				eff_rden_z = rden[m][n];
			}
			else
			{
				/*有效介质参数*/
				eff_rden_x = (rden[m][n] + rden[m + 1][n])*0.5;
				eff_rden_z = (rden[m][n] + rden[m][n + 1])*0.5;
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