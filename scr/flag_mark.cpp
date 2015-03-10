#include"topoft.h"
/*
对起伏地表网格进行分类
网格标记：flag 值：
	空气  -1
	介质   0
	H-point  1
	OL-point 2
	OR-point 3
	IR-point 4
	IL-point 5
	VL-point 6
	VR-point 7
vp：通过纵波速度来区分
NX_f, NX_e, NZ_f, NZ_e:模型大小
depth:起伏地表，最低的网格点，默认值为z深度全部。
*/
void flag_mark(int **flag, float **vp, int NX_f, int NX_e, int NZ_f, int NZ_e, int &depth)
{
	cout << "\n-------------" << endl;
	cout << "Mark the Grid" << endl;
	cout << "-------------" << endl;
	
	depth = 1;
	/*起伏自由地表分类*/
	for (int i = NX_f; i <= NX_e; i++)
	{
		for (int j = NZ_f; j <= NZ_e; j++)
		{
			/* H表示水平自由表面边界格子，格子上方为真空，下方为固体介质，左右各有至少一个边界格子*/
			if (vp[i][j] > 340 && vp[i][j - 1] <= 340 && vp[i - 1][j] > 340 && vp[i + 1][j] > 340)
			{
				/*H-point*/
				flag[i][j] = 1;
				if (depth < j)
					depth = j;//记录最深自由地表深度
			}
			/*OL表示真空位于左方的外角点边界格子，其右方和下方各至少有一个边界格子，右下方固体介质*/
			else if (vp[i][j] > 340 && vp[i][j - 1] <= 340 && vp[i - 1][j] <= 340 &&
				vp[i + 1][j] > 340 && vp[i][j + 1] > 340 && vp[i + 1][j + 1] > 340)
			{
				/*OL-points*/
				flag[i][j] = 2;
				if (depth < j)
					depth = j;
			}
			/*OR表示真空位于右上方的外角点边界格子，其左方和下方各有至少一个边界格子*/
			else if (vp[i][j] > 340 && vp[i][j - 1] <= 340 && vp[i + 1][j] <= 340 &&
				vp[i - 1][j] > 340 && vp[i][j + 1] > 340 && vp[i - 1][j + 1] > 340)
			{
				/*OR point*/
				flag[i][j] = 3;
				if (depth < j)
					depth = j;
			}
			/* IR表示真空位于右上方的内角点边界格子，其上方和右方各有至少一个边界格子*/
			else if (vp[i + 1][j] > 340 && vp[i][j] > 340 && vp[i - 1][j] > 340
				&& vp[i][j + 1] > 340 && vp[i][j - 1] > 340 && vp[i + 1][j - 1] <= 340)
			{
				/*IR point*/
				flag[i][j] = 4;
				if (depth < j)
					depth = j;
			}
			/* IL表示真空位于左上方的内角点边界格子，其上方和左方各有一边界格子，左上方为真空，右下方为固体介质*/
			else if (vp[i + 1][j] > 340 && vp[i][j] > 340 && vp[i - 1][j] > 340
				&& vp[i][j + 1] > 340 && vp[i][j - 1] > 340 && vp[i - 1][j - 1] <= 340)
			{
				/*IL point*/
				flag[i][j] = 5;
				if (depth < j)
					depth = j;
			}
			/*VL表示真空在左边的垂直自由表面边界格子，格子左方为真空，上下各有至少一个边界格子*/
			else if (vp[i][j] > 340 && vp[i - 1][j] <= 340 && vp[i][j - 1] > 340 && vp[i][j + 1] > 340)
			{
				flag[i][j] = 6;
				if (depth < j)
					depth = j;
			}
			/*VR表示真空位于右方的垂直自由表面边界格子，其上下各有至少一个边界格子*/
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