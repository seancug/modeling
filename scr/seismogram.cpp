#include"topoft.h"

/*���沨��*/
void seismogram(int csamp, int ntr, int **recpos, float**sectionvx, float**sectionvz,
	float **vx, float **vz)
{
	extern int NDT;
	extern string SEISMO;
	int ins, nxrec, nzrec;


	ins = csamp / NDT;//���Լ����ʱ�����counting

	for (int itr = 1; itr <= ntr; itr++){
		nxrec = recpos[1][itr];
		nzrec = recpos[2][itr];
		if (SEISMO == "V")//�����ٶȳ�
		{
			sectionvx[itr][ins] = vx[nxrec][nzrec];
			sectionvz[itr][ins] = vz[nxrec][nzrec];
		}
		else if (SEISMO == "S")
		{
			//ˮ����������δ���
			/*sectionvx[itr][ins] = -sxx[nyrec][nxrec][nzrec]
			- syy[nyrec][nxrec][nzrec]
			- szz[nyrec][nxrec][nzrec];*/
		}
		else if (SEISMO == "VS" || SEISMO == "SV")
		{
			
		}
	}
}