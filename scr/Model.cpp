/*水平双层模型*/

#include"topoft.h"

void model(float **rden, float **lam, float **u, int NX_f, int NX_e, int NZ_f, int NZ_e)
{
	cout << "\nGenerate half space homongenous media..." << endl;
	/* 速度，密度 */
	const float vp0 = 1000.0, vs0 = 300.0, rho0 = 1500.0;
	cout << "P velocity is: " << vp0 << "\tS velocity is: "
		<< vs0 << "\tdensity is: " << rho0 << endl;
	//const double vp1 = 4000.0, vs1 = 3000.0, rh1 = 2200.0;
	for (int i = NX_f; i <= NX_e; i++)
	{
		for (int j = NZ_f; j <= NZ_e; j++)
		{
			rden[i][j] = 1 / rho0; //密度倒数
			lam[i][j] = (pow(vp0, 2) - 2 * pow(vs0, 2))*rho0; //拉梅常数
			u[i][j] = pow(vs0, 2)*rho0;
		}
	}
	cout << "Model generated compelete..." << endl;
}