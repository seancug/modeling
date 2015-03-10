#include"topoft.h"
/*
������ر�������з���
�����ǣ�flag ֵ��
	����  -1
	����   0
	H-point  1
	OL-point 2
	OR-point 3
	IR-point 4
	IL-point 5
	VL-point 6
	VR-point 7
vp��ͨ���ݲ��ٶ�������
NX_f, NX_e, NZ_f, NZ_e:ģ�ʹ�С
depth:����ر���͵�����㣬Ĭ��ֵΪz���ȫ����
*/
void flag_mark(int **flag, float **vp, int NX_f, int NX_e, int NZ_f, int NZ_e, int &depth)
{
	cout << "\n-------------" << endl;
	cout << "Mark the Grid" << endl;
	cout << "-------------" << endl;
	
	depth = 1;
	/*������ɵر����*/
	for (int i = NX_f; i <= NX_e; i++)
	{
		for (int j = NZ_f; j <= NZ_e; j++)
		{
			/* H��ʾˮƽ���ɱ���߽���ӣ������Ϸ�Ϊ��գ��·�Ϊ������ʣ����Ҹ�������һ���߽����*/
			if (vp[i][j] > 340 && vp[i][j - 1] <= 340 && vp[i - 1][j] > 340 && vp[i + 1][j] > 340)
			{
				/*H-point*/
				flag[i][j] = 1;
				if (depth < j)
					depth = j;//��¼�������ɵر����
			}
			/*OL��ʾ���λ���󷽵���ǵ�߽���ӣ����ҷ����·���������һ���߽���ӣ����·��������*/
			else if (vp[i][j] > 340 && vp[i][j - 1] <= 340 && vp[i - 1][j] <= 340 &&
				vp[i + 1][j] > 340 && vp[i][j + 1] > 340 && vp[i + 1][j + 1] > 340)
			{
				/*OL-points*/
				flag[i][j] = 2;
				if (depth < j)
					depth = j;
			}
			/*OR��ʾ���λ�����Ϸ�����ǵ�߽���ӣ����󷽺��·���������һ���߽����*/
			else if (vp[i][j] > 340 && vp[i][j - 1] <= 340 && vp[i + 1][j] <= 340 &&
				vp[i - 1][j] > 340 && vp[i][j + 1] > 340 && vp[i - 1][j + 1] > 340)
			{
				/*OR point*/
				flag[i][j] = 3;
				if (depth < j)
					depth = j;
			}
			/* IR��ʾ���λ�����Ϸ����ڽǵ�߽���ӣ����Ϸ����ҷ���������һ���߽����*/
			else if (vp[i + 1][j] > 340 && vp[i][j] > 340 && vp[i - 1][j] > 340
				&& vp[i][j + 1] > 340 && vp[i][j - 1] > 340 && vp[i + 1][j - 1] <= 340)
			{
				/*IR point*/
				flag[i][j] = 4;
				if (depth < j)
					depth = j;
			}
			/* IL��ʾ���λ�����Ϸ����ڽǵ�߽���ӣ����Ϸ����󷽸���һ�߽���ӣ����Ϸ�Ϊ��գ����·�Ϊ�������*/
			else if (vp[i + 1][j] > 340 && vp[i][j] > 340 && vp[i - 1][j] > 340
				&& vp[i][j + 1] > 340 && vp[i][j - 1] > 340 && vp[i - 1][j - 1] <= 340)
			{
				/*IL point*/
				flag[i][j] = 5;
				if (depth < j)
					depth = j;
			}
			/*VL��ʾ�������ߵĴ�ֱ���ɱ���߽���ӣ�������Ϊ��գ����¸�������һ���߽����*/
			else if (vp[i][j] > 340 && vp[i - 1][j] <= 340 && vp[i][j - 1] > 340 && vp[i][j + 1] > 340)
			{
				flag[i][j] = 6;
				if (depth < j)
					depth = j;
			}
			/*VR��ʾ���λ���ҷ��Ĵ�ֱ���ɱ���߽���ӣ������¸�������һ���߽����*/
			else if (vp[i][j] > 340 && vp[i + 1][j] <= 340 && vp[i][j - 1] > 340 && vp[i][j + 1] > 340)
			{
				flag[i][j] = 7;
				if (depth < j)
					depth = j;
			}
			else if (vp[i][j] <= 340)
			{
				flag[i][j] = -1;
			}
		}
	}
	depth += 8;
}