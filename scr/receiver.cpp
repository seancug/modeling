#include"topoft.h"

/*获得检波器*/
int** receiver(int &ntr, int **flag)
{
	extern float XREC1, XREC2, ZREC1, ZREC2;
	extern float DRX, DRZ, DX, DZ, FW;
	extern string RECTYPE, READREC;
	extern char REC_FILE[STRING_SIZE];
	extern int NZ;
	/*吸收层数*/
	int ifw_x, ifw_z;
	int nxrec1 = 0, nxrec2 = 0, nzrec1 = 0, nzrec2 = 0;
	int drx = 0, drz = 0;
	int ** respos =NULL; //检波器位置坐标
	int itr = 0;

	if (READREC == "YES")
	{
		/*从文件中读取，功能未添加*/
	}
	else if (READREC == "NO")
	{
		/*判断吸收层，不不知检波器，功能未添加*/
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
			/*水平检波器布置*/
			for (int i = nxrec1; i <= nxrec2; i += drx)
			{
				itr++;
				respos[1][itr] = i;
				respos[2][itr] = nzrec1;//nzrec1=nzrec2
			}
		}
		else if (RECTYPE == "TOPOGRAPHY")
		{
			/*固定x起始和终止范围，z坐标沿着起伏地表*/
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
				write_repos(respos, ntr);//输出检波器信息
			}
		}
		else
		{
			/*无类型，报错*/
		}
	}
	
	return respos;
}