#include"topoft.h"

void outputdata(char *s, float **datain, int nrl, int nrh, int ncl, int nch, int prec)
{
	ofstream outfile(s);

	outfile.setf(ios::fixed);
	outfile.setf(ios::showpoint);
	outfile.precision(prec);
	
	
	for (int j = ncl; j <= nch; j++)
	{
		for (int i = nrl; i <= nrh; i++)
		{
			outfile << datain[i][j] << "\t";
		}
		outfile << endl;
	}
	outfile.close();
}

void outputdata(char *s, int **datain, int nrl, int nrh, int ncl, int nch)
{
	/*ofstream outfile(s);

	outfile.setf(ios::fixed);
	outfile.setf(ios::showpoint);
	outfile.right;
	//outfile.precision(prec);

	for (int j = ncl; j <= nch; j++)
	{
		for (int i = nrl; i <= nrh; i++)
		{
			outfile << datain[i][j] << "\t";
		}
		outfile << endl;
	}
	outfile.close();*/
	FILE *fp;
	fp = fopen("mark.temp","w+");
	for (int i = nrl; i <= nrh; i++)
	{
		fprintf(fp, "%4d\t", i);
	}
	fprintf(fp, "\n");
	if (fp != NULL)
	{
		for (int j = ncl; j <= nch; j++)
		{
			for (int i = nrl; i <= nrh; i++)
			{
				fprintf(fp, "%4d\t", datain[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
}