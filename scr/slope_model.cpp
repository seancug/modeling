#include"topoft.h"

void slope_model(float **rden, float **lam, float **u, float **vp, int NX_f, int NX_e, int NZ_f, int NZ_e)
{
	cout << "\n-------------------------------------" << endl;
	cout << "Generate Half Space Homongenous Media"<<endl;
	cout << "-------------------------------------" << endl;
	/* �ٶȣ��ܶ� */
	extern float DX, DZ, FW;
	int ifw_x, ifw_z;
	ifw_x = iround(FW / DX);
	ifw_z = iround(FW / DZ);

	//const double vp0 = 1000.0, vs0 = 500.0, rho0 = 1500.0;
	const float vp0 = 1000.0, vs0 = 300, rho0 = 1500.0;
	const float vp1 = 1300.0, vs1 = 1300.0 / sqrt(3), rho1 = 1500.0;
	cout << "\tP0 velocity is: " << vp0 << "\tS0 velocity is: "
		<< vs0 << "\tdensity0 is: " << rho0 << endl;
	/*cout << "P1 velocity is: " << vp1 << "\tS1 velocity is: "
		<< vs1 << "\tdensity0 is: " << rho1 << endl;*/
    //const double vp1 = 4000.0, vs1 = 3000.0, rh1 = 2200.0;

	for (int i = NX_f + ifw_x; i <= NX_e - ifw_x; i++)
    {
		for (int j = NZ_f; j <= NZ_e - ifw_z; j++)
		{
			//��бģ��
			//if (j < (500.0 - (i - ifw_x - NX_f) / 4.0))
			//{
			//	rden[i][j] = 0;
			//	vp[i][j] = 0;
			//}
			//else
			//{
			//	rden[i][j] = 1 / rho0; //�ܶȵ���
			//	lam[i][j] = (pow(vp0, 2) - 2 * pow(vs0, 2))*rho0; //��÷����
			//	u[i][j] = pow(vs0, 2)*rho0;
			//	vp[i][j] = vp0;
			//}

			////ˮƽģ��
			//if (j < 20)
			//{
			//	rden[i][j] = 0;
			//	vp[i][j] = 0;
			//}
			//else
			//{
			//	rden[i][j] = 1.0 / rho0; //�ܶȵ���
			//	lam[i][j] = (pow(vp0, 2) - 2.0 * pow(vs0, 2))*rho0; //��÷����
			//	u[i][j] = pow(vs0, 2)*rho0;
			//	vp[i][j] = vp0;
			//}

			////��бģ��,����������б,���30����б����
			//if (j < 500 - (i - ifw_x - NX_f) *(sqrt(3.0) / 3.0) || (j < 10))
			//{
			//	rden[i][j] = 0;
			//	vp[i][j] = 0;
			//}
			//else
			//{
			//	rden[i][j] = 1 / rho0; //�ܶȵ���
			//	lam[i][j] = (pow(vp0, 2) - 2 * pow(vs0, 2))*rho0; //��÷����
			//	u[i][j] = pow(vs0, 2)*rho0;
			//	vp[i][j] = vp0;
			//}

			////��б˫��ģ��,����������б
			//if (j < 400 - (i - ifw_x - NX_f) / 6.0 || (j < 10))
			//{
			//	rden[i][j] = 0;
			//	vp[i][j] = 0;
			//}
			//else
			//{
			//	if (j < 520)
			//	{
			//		rden[i][j] = 1 / rho0; //�ܶȵ���
			//		lam[i][j] = (pow(vp0, 2) - 2 * pow(vs0, 2))*rho0; //��÷����
			//		u[i][j] = pow(vs0, 2)*rho0;
			//		vp[i][j] = vp0;
			//	} 
			//	else
			//	{
			//		rden[i][j] = 1 / rho1; //�ܶȵ���
			//		lam[i][j] = (pow(vp1, 2) - 2 * pow(vs1, 2))*rho1; //��÷����
			//		u[i][j] = pow(vs1, 2)*rho1;
			//		vp[i][j] = vp1;
			//	}
			//	
            //}

			//��бģ��,������б
			if (j < 10 + (i - ifw_x - NX_f) *(sqrt(3.0)) / 3.0 || (j < 10))
			{
				rden[i][j] = 0;
				vp[i][j] = 0;
			}
			else
			{
				rden[i][j] = 1 / rho0; //�ܶȵ���
				lam[i][j] = (pow(vp0, 2) - 2 * pow(vs0, 2))*rho0; //��÷����
				u[i][j] = pow(vs0, 2)*rho0;
				vp[i][j] = vp0;
			}


			////VL ģ��
			//if (i < ifw_x + 10)
			//{
			//	rden[i][j] = 0;
			//	vp[i][j] = 0;
			//}
			//else
			//{
			//	rden[i][j] = 1 / rho0; //�ܶȵ���
			//	lam[i][j] = (pow(vp0, 2) - 2 * pow(vs0, 2))*rho0; //��÷����
			//	u[i][j] = pow(vs0, 2)*rho0;
			//	vp[i][j] = vp0;
			//}
		}
	}

	int ii = 0;
	for (int i = NX_f - 1; i < NX_f + ifw_x; i++)
	{
		ii = NX_e - i + NX_f;
		for (int j = NZ_f; j <= NZ_e - ifw_z; j++)
		{
			vp[i][j] = vp[NX_f + ifw_x][j];
			vp[ii][j] = vp[NX_e - ifw_x][j];

			rden[i][j] = rden[NX_f + ifw_x][j];
			rden[ii][j] = rden[NX_e - ifw_x][j];
			
			lam[i][j] = lam[NX_f + ifw_x][j];
			lam[ii][j] = lam[NX_e - ifw_x][j];

			u[i][j] = u[NX_f + ifw_x][j];
			u[ii][j] = u[NX_e - ifw_x][j];
		}
	}

	for (int j = NZ_e - ifw_z + 1; j <= NZ_e + 1; j++)
	{
		for (int i = NX_f - 1; i <= NX_e + 1; i++)
		{
			vp[i][j] = vp[i][NZ_e - ifw_z];
			rden[i][j] = rden[i][NZ_e - ifw_z];
			lam[i][j] = lam[i][NZ_e - ifw_z];
			u[i][j] = u[i][NZ_e - ifw_z];
		}
    }
}
