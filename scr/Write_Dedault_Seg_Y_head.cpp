#include"topoft.h"

void Write_Dedault_Seg_Y_head(FILE *fp)
{
	char buf[81];

	sprintf(buf, "%-79.79s\n", "C    China university of geoscience (Wuhan)");
	fwrite(buf, 1, 80, fp);

	sprintf(buf, "%-79.79s\n", "C    Institute of Geophysics & Geomatics");
	fwrite(buf, 1, 80, fp);

	sprintf(buf, "%-79.79s\n", "C    Geophysical Department");
	fwrite(buf, 1, 80, fp);

	sprintf(buf, "%-79.79s\n", "C    Author: Zhang Qiang");
	fwrite(buf, 1, 80, fp);

	sprintf(buf, "%-10.10s %-68.68s\n", "C    Date:", __DATE__);
	fwrite(buf, 1, 80, fp);

	for (int i = 0; i < 35; i++)
	{
		sprintf(buf, "%-79.79s\n", "C");
		fwrite(buf, 1, 80, fp);
	}
}