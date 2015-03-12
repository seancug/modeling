#include"topoft.h"

bool checkFD(string &message, float **lam, float **u, float **rden, int NX_f, int NX_e, int NZ_f, int NZ_e)
{
	message.append("\n--------------------------------\n");
	message.append("Check Finite Difference Stabilty\n");
	message.append("--------------------------------\n");
	extern float DT, DX, DZ;
	/*中间变量*/
	float vpmax = 0.0, vp = 0.0;
	float Courant_number = 0.0;
	float suma = 0.0;
	for (int i = NX_f; i <= NX_e; i++)
	{
		for (int j = NZ_f; j <= NZ_e; j++)
		{
			vp = sqrt((lam[i][j] + 2 * u[i][j])*rden[i][j]);
			if (vpmax < vp)
			{
				vpmax = vp; //寻找最大的vp速度
			}
		}
	}
	/*判断稳定性条件*/
	Courant_number = DT*vpmax*sqrt(1.0 / pow(DX, 2) + 1.0 / pow(DZ, 2));
	suma = fabs(a1) + fabs(a2) + fabs(a3) + fabs(a4);
	if (Courant_number <= 1.0 / suma)
	{
		//cout << "the explicit time scheme is stable..." << endl;
		message.append("\tthe explicit scheme is stable...\n");
		float rightDT = 0;
		rightDT = (1.0 / suma) / (vpmax*sqrt(1.0 / pow(DX, 2) + 1.0 / pow(DZ, 2)));
		char temp[252] = "";
		//sprintf(temp, "\tCourant_number: %.3e\n\t1.0/suma: %.3e\n", Courant_number, vpmax);
		//message.append(temp);
		sprintf(temp, "\tthe DT shoule small than %.3e\n", rightDT);
		message.append(temp);
		return true;
	}
	else
	{
		float rightDT = 0;
		rightDT = (1.0 / suma) / (vpmax*sqrt(1.0 / pow(DX, 2) + 1.0 / pow(DZ, 2)));
		//cout << "time step is too large, simulation will be unstable...\n" << endl;
		message.append("\ttime step is too large, simulation will be unstable...\n"); 
		char temp[252] = "";
		sprintf(temp, "\tthe DT shoule small than %.3e\n", rightDT);
		message.append(temp);
		return false;
	}
	
}
