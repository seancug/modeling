#include"topoft.h"

void writedsk(FILE *fp_out, float amp1, string format)
{
	//char b = ' ';
	if (format == "SU")
		nrerror(" Sorry, SU-format for snapshots not implemented yet. \n");
	else if (format == "ASCII")
		fprintf(fp_out, "%e\t", amp1);
	else if (format == "BINARY")
	{
		fwrite(&amp1, sizeof(float), 1, fp_out);
//		fwrite(&b, sizeof(char), 1, fp_out);
	}
	else if (format == "GRID")
		nrerror(" Sorry, SURFER-format for snapshots not implemented yet. \n");
	else
	{
		/*printf(" Don't know the format for the snapshot-data !\n");
		nrerror(" No output was written. ");*/
	}
}