#include"topoft.h"

void topo_update_v(int nt, int **flag, float ** vx, float ** vz,
	float ** sxx, float ** szz, float ** sxz,
	float **  rden, float ** absorb_coeff)
{
	extern float DT, DX, DZ, FW;
	extern int   NX, NZ;

	/**������������ص�����ˮƽ����ʹ�ֱ�����󵼷ֿ�����**/
	float sxx_dx = 0, sxz_dz = 0;
	float sxz_dx = 0, szz_dz = 0;
	float eff_rden_x = 0, eff_rden_z = 0;

	/*�ȶ�Z���������*/
	for (int m = 1; m <= NX; m++)
	{
		for (int n = 1; n <= NZ; n++)
		{
			/******������**********/
			/*��Ч���ʲ���*/
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
				/*z����Ӧ���ռ�ƫ��*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z����Ӧ���ռ�ƫ��*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 1:/*H-point,��������,Check*/
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

				/*z����Ӧ���ռ�ƫ��*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z����Ӧ���ռ�ƫ��*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 2:/*OL-point,sxz�����ɱ��洦��Ϊ��,szz,sxz��ֱ����*/
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

				/*z����Ӧ���ռ�ƫ��*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z����Ӧ���ռ�ƫ��*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 3:/*OR-point sxz���ɱ������Ϊ�㣬szz,sxz��ֱ����*/
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

				/*z����Ӧ���ռ�ƫ��*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z����Ӧ���ռ�ƫ��*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
			case 4:/*IR-point sxx,szzλ�����ɱ߽���,��Ϊ0,sxz��ֱ����*/
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
				
				/*z����Ӧ���ռ�ƫ��*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z����Ӧ���ռ�ƫ��*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 5:/*IL-point �޾���ֵ,sxx,szzλ�����ɱ߽���,��Ϊ0*/
				sxx[m][n] = 0;
				szz[m][n] = 0;
				
				/*z����Ӧ���ռ�ƫ��*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z����Ӧ���ռ�ƫ��*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 6:/*VL-point*/
				/*z����Ӧ���ռ�ƫ��*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z����Ӧ���ռ�ƫ��*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			case 7:/*VR-point*/
				/*z����Ӧ���ռ�ƫ��*/
				sxz_dz = (
					a1*(sxz[m][n] - sxz[m][n - 1]) +
					a2*(sxz[m][n + 1] - sxz[m][n - 2]) +
					a3*(sxz[m][n + 2] - sxz[m][n - 3]) +
					a4*(sxz[m][n + 3] - sxz[m][n - 4])) / DZ;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxz_dz;

				/*z����Ӧ���ռ�ƫ��*/
				szz_dz = (
					a1*(szz[m][n + 1] - szz[m][n]) +
					a2*(szz[m][n + 2] - szz[m][n - 1]) +
					a3*(szz[m][n + 3] - szz[m][n - 2]) +
					a4*(szz[m][n + 4] - szz[m][n - 3])) / DZ;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z *szz_dz;
				break;
			}
		}
	}

	/*�ȶ�x���������*/
	for (int m = 1; m <= NX; m++)
	{
		for (int n = 1; n <= NZ; n++)
		{
			/*��Ч���ʲ���*/
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
			case -1:/*����*/
				break;
			case 0:/*�ڲ�����*/
				/******ˮƽ��**********/
				/*Ӧ���ռ�ƫ��*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*Ӧ���ռ�ƫ��*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 1:/*H-point*/
				/******ˮƽ��**********/
				/*Ӧ���ռ�ƫ��*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*Ӧ���ռ�ƫ��*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 2:/*OL-point,sxz�����ɱ��洦��Ϊ��,szz,sxz��ֱ����,ǰ���Ѿ����ù�*/
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
				
				/******ˮƽ��**********/
				/*Ӧ���ռ�ƫ��*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*Ӧ���ռ�ƫ��*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 3:/*OR-point sxz���ɱ������Ϊ��,szz,sxz��ֱ����*/
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
				
				/******ˮƽ��**********/
				/*Ӧ���ռ�ƫ��*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*Ӧ���ռ�ƫ��*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 4:/*IR-point*/
				/******ˮƽ��**********/
				/*Ӧ���ռ�ƫ��*/
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

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*Ӧ���ռ�ƫ��*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 5:/*IL-point*/
								/******ˮƽ��**********/
				/*Ӧ���ռ�ƫ��*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*Ӧ���ռ�ƫ��*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 6:/*VL-point ˮƽ����*/
				sxx[m][n] = 0;
				sxx[m - 1][n] = -sxx[m + 1][n];
				sxx[m - 2][n] = -sxx[m + 2][n];
				sxx[m - 3][n] = -sxx[m + 3][n];

				sxz[m - 1][n] = -sxz[m][n];
				sxz[m - 2][n] = -sxz[m + 1][n];
				sxz[m - 3][n] = -sxz[m + 2][n];
				sxz[m - 4][n] = -sxz[m + 3][n];

				/******ˮƽ��**********/
				/*Ӧ���ռ�ƫ��*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*Ӧ���ռ�ƫ��*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			case 7:/*VR-point ˮƽ����*/
				sxx[m][n] = 0;
				sxx[m + 1][n] = -sxx[m - 1][n];
				sxx[m + 2][n] = -sxx[m - 2][n];
				sxx[m + 3][n] = -sxx[m - 3][n];
				sxx[m + 4][n] = -sxx[m - 4][n];

				sxz[m][n] = -sxz[m - 1][n];
				sxz[m + 1][n] = -sxz[m - 2][n];
				sxz[m + 2][n] = -sxz[m - 3][n];
				sxz[m + 3][n] = -sxz[m - 4][n];	

				/******ˮƽ��**********/
				/*Ӧ���ռ�ƫ��*/
				sxx_dx = (
					a1*(sxx[m + 1][n] - sxx[m][n]) +
					a2*(sxx[m + 2][n] - sxx[m - 1][n]) +
					a3*(sxx[m + 3][n] - sxx[m - 2][n]) +
					a4*(sxx[m + 4][n] - sxx[m - 3][n])) / DX;

				/*vx �ٶȵ���*/
				vx[m][n] += DT * eff_rden_x * sxx_dx;

				/*Ӧ���ռ�ƫ��*/
				sxz_dx = (
					a1*(sxz[m][n] - sxz[m - 1][n]) +
					a2*(sxz[m + 1][n] - sxz[m - 2][n]) +
					a3*(sxz[m + 2][n] - sxz[m - 3][n]) +
					a4*(sxz[m + 3][n] - sxz[m - 4][n])) / DX;

				/*vz �ٶȵ���*/
				vz[m][n] += DT * eff_rden_z * sxz_dx;
				break;
			}
		}
	}

	/*���ձ߽��������,Ч�����Ǻܺ�*/
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