#include"topoft.h"

float* rd_sour(int *nts, float *nt, FILE* fp_source)
{
	/* 局部变量 */
	float *psource;
	int i;
	
	if (fp_source == NULL) nrerror(" Source file could no be opened !");
	fscanf(fp_source,"%i", nts); 
	fscanf(fp_source, "%f", nt);
	psource = dvector(1, *nts);
	for (i = 1; i <= *nts; i++) fscanf(fp_source, "%e", &psource[i]);
	fclose(fp_source);
	return psource;
}