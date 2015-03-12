#include"topoft.h"

/*��ü첨��*/
int** receiver(int &ntr, int **flag)
{
	extern float XREC1, XREC2, ZREC1, ZREC2;
	extern float DRX, DRZ, DX, DZ, FW;
	extern string RECTYPE, READREC;
	extern char REC_FILE[STRING_SIZE];
	extern int NZ;
	/*���ղ���*/
	int ifw_x, ifw_z;
	int nxrec1 = 0, nxrec2 = 0, nzrec1 = 0, nzrec2 = 0;
	int drx = 0, drz = 0;
	int ** respos =NULL; //�첨��λ������
	int itr = 0;

	if (READREC == "YES")
	{
		/*���ļ��ж�ȡ������δ���*/
	}
	else if (READREC == "NO")
	{
		/*�ж����ղ㣬����֪�첨��������δ���*/
		ifw_x = iround(FW / DX);
		ifw_z = iround(FW / DZ);
		nxrec1 = iround(XREC1 / DX + 0.5);
		nxrec2 = iround(XREC2 / DX - 0.5);
		nzrec1 = iround(ZREC1 / DZ + 0.5);
		nzrec2 = iround(ZREC2 / DZ - 0.5);
		drx = iround(DRX / DX);
		drz = iround(DRZ / DZ);

		
		if (RECTYPE == "PLANER")
		{
			ntr = iround((nxrec2 - nxrec1) / drx - 0.5) + 1;
			if (ntr != 0)
			{
				respos = imatrix(1, 2, 1, ntr);
			}
			
			itr = 0;
			/*ˮƽ�첨������*/
			for (int i = nxrec1; i <= nxrec2; i += drx)
			{
				itr++;
				respos[1][itr] = i;
				respos[2][itr] = nzrec1;//nzrec1=nzrec2
			}
		}
		else if (RECTYPE == "TOPOGRAPHY")
		{
			/*�̶�x��ʼ����ֹ��Χ��z������������ر�*/
			ntr = iround((nxrec2 - nxrec1) / drx - 0.5) + 1;
			if (ntr != 0)
			{
				respos = imatrix(1, 2, 1, ntr);
			}

			itr = 0;
			for (int i = nxrec1; i <= nxrec2; i += drx)
			{
				itr++;
				respos[1][itr] = i;
				for (int j = 1; j < NZ - ifw_z; j++)
				{
					if (flag[i][j] != -1 && flag[i][j] != 0)
					{
						respos[2][itr] = j + 2;
						break;
					}
				}
			}

			string tmp(REC_FILE);
			if (tmp != "")
			{
				write_repos(respos, ntr);//����첨����Ϣ
			}
		}
		else
		{
			/*�����ͣ�����*/
		}
	}
	
	return respos;
}