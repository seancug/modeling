#include"topoft.h"
/*����Լ��첨��λ����Ϣ
**repos  �첨������
  ntr    �첨������*/
void write_repos(int **repos,int ntr)
{
	extern char REC_FILE[STRING_SIZE];
	FILE *pos_fp;
	extern float DX, DZ;
    char tem[124];
    sprintf(tem, "%s%s","../gather/", REC_FILE);
	pos_fp = fopen(REC_FILE, "w");
	fprintf(pos_fp, "trace num\tx position\ty position\n");
	for (int i = 1; i <= ntr;i++)
	{
		fprintf(pos_fp, "%5d\t%8.2f\t%8.2f\n", i, DX*repos[1][i], DZ*repos[2][i]);
	}
	fclose(pos_fp);
}
