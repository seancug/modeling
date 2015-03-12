#include"topoft.h"

void add_sour(int nt, float ** vx, float **vz, float **rden)
{
	extern float  DT, F0, ANGLE;
	extern string SOURCE_WAVELET, SOURCE_TYP;
	extern float DX, DZ, XS, ZS;
	extern float AMP;
	int nxs, nzs; //震源所在的网格点
	float a, t;//时间
	float t0 = 0.08;
	float sour_term;
	float factor;
	factor = (-1.0)*AMP;
	nxs = iround(XS / DX) + 1;
	nzs = iround(ZS / DZ) + 1;
	t = (nt - 1)*DT;

	if (SOURCE_WAVELET == "RICKER")
	{
		a = PI*PI*F0*F0;
		sour_term = factor*(1.0 - 2.0*a*pow(t - t0, 2))*exp(-a*pow(t - t0, 2));
	}
	else if (SOURCE_WAVELET == "FILE")
	{
		//功能未添加，除了雷克子波。还可以增加其余子波

	}

	if (SOURCE_TYP == "EXPLOSION")
	{
		vx[nxs][nzs] -= sour_term*DT*rden[nxs][nzs];
		vz[nxs][nzs] -= sour_term*DT*rden[nxs][nzs];

		vx[nxs + 1][nzs] += sour_term*DT*rden[nxs + 1][nzs];
		vz[nxs + 1][nzs] -= sour_term*DT*rden[nxs + 1][nzs];

		vx[nxs + 1][nzs + 1] += sour_term*DT*rden[nxs + 1][nzs + 1];
		vz[nxs + 1][nzs + 1] += sour_term*DT*rden[nxs + 1][nzs + 1];

		vx[nxs][nzs + 1] -= sour_term*DT*rden[nxs][nzs + 1];
		vz[nxs][nzs + 1] += sour_term*DT*rden[nxs][nzs + 1];
	}
	else if (SOURCE_TYP == "POINT")
	{
		vx[nxs][nzs] += sin(ANGLE)*sour_term*DT*rden[nxs][nzs];
		vz[nxs][nzs] += cos(ANGLE)*sour_term*DT*rden[nxs][nzs];
	}
	else
		nrerror(" Give me right source type! ");
}
