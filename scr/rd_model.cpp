#include"topoft.h"

/*读入模型*/
float** rd_model(char* filename, int &NX, int &NZ,float &DX,float &DZ)
{
	cout << "--------------------------------------" << endl;
	cout << "Read Model \"" << filename << "\" begin..." << endl;

	//extern float DX, DZ, FW;
	//int ifw_x, ifw_z;
	//ifw_x = iround(FW / DX);
	//ifw_z = iround(FW / DZ);

	FILE *fp;
	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		cout << "Open file \""<<filename<< "\" fails..." << endl;
		return false;
	}
	int ret;		  /* fscanf() return value		 */
	int n1=0;			  /* number of floats per line	   	   */
	int n2=0;           /* number of vector	   	   */
	unsigned long int vnum;   /* how many floats we have filled in the */
	float xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0, zmin = 0.0, zmax = 0.0;
	/*  current vector already */
	unsigned long int vec;    /* number of vector we are reading	*/
	//n1 = NX_e - NX_f + 1;
	char INDEX[6];
	fscanf(fp, "%s", INDEX);
	fscanf(fp, " %d %d ", &n1, &n2);
	fscanf(fp, " %f %f ", &xmin, &xmax);
	fscanf(fp, " %f %f ", &ymin, &ymax);
	fscanf(fp, " %f %f ", &zmin, &zmax);
	if (NX == 0)
		NX = n1;
	else if (NX != n1){
		fscanf(stderr, "NX size don't match the model file size\n");
		return false;
	}
	if (NZ == 0)
		NZ = n2;
	else if (NZ != n2){
		fscanf(stderr, "NZ size don't match the model file size\n");
		return false;
	}

	float dx, dz;
	dx = fabs((xmax - xmin) / (n1 - 1));
	dz = fabs((ymax - ymin) / (n2 - 1));
	if (DX == 0)
		DX = dx;
	else if (DX != dx){
		fscanf(stderr, "NX size don't match the model file size\n");
		return false;
	}
	if (DZ == 0)
		DZ = dz;
	else if (DZ != dz){
		fscanf(stderr, "NZ size don't match the model file size\n");
		return false;
	}
	float **Par;
	Par = dmatrix(0, NX + 1, 1, NZ + 1);//前后多开辟一个
	vnum = 1;		/* offset inside vector, in range 0->n1 */
	vec = n2;		/* vector number, ends up as n2	 */

	unsigned short int p1 = 0;
	while (1) {
		ret = fscanf(fp, " %f ", &Par[vnum][vec]);
		/*  fscanf returns: 0   : characters there, but no conversion (error)
		*		  EOF : eof before conversion
		*		  else: number of conversions
		*/
		if (ret == EOF) {
			if (vnum != 1)
				fprintf(stderr, "\n\tvector #%lu, float #%lu: we encountered an EOF but the vector has not been read in completely", vec + 1, vnum + 1);
			break;  /* else everything is okay: get out of the loop */
		}
		if (ret == 0)
			fprintf(stderr, "\n\tcould not read float in vector #%lu, float #%lu", vec + 1, vnum + 1);
		vnum++;
		if (vnum == (n1 + 1)) {
			vnum = 1;
			static unsigned char w[] = "///-";
			if (p1++ == 3) p1 = 0;
			printf("\r\tread %s :%.1f %c %c", filename, 100.0*(n2 - vec + 1) / n2, '%', w[p1]);
			fflush(stdout);
			vec--;
		}
	}
	printf("\n");
	cout << "--------------------------------------" << endl;
	return Par;
}
