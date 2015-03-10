#include"topoft.h"
#include<cstdlib>
#include<cstdio>
#include<cstddef>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr, "run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}

void parerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr, "Parameter ");
	fprintf(stderr, "%s ", error_text);
	fprintf(stderr, "is not put in!\n");
}


float **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;

	/* allocate pointers to rows */
	m = (float **)malloc((size_t)((nrow + NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (float *)malloc((size_t)((nrow*ncol + NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* initializing matrix */
	for (i = nrl; i <= nrh; i++)
		for (j = ncl; j <= nch; j++) m[i][j] = 0.0;
	/* return pointer to array of pointers to rows */
	return m;
}

void free_dmatrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

float *dvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;
	int i;
	v = (float *)malloc((size_t)((nh - nl + 1 + NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	for (i = 0; i < (nh - nl + 1 + NR_END); i++) v[i] = 0.0;
	return v - nl + NR_END;
}

void free_dvector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

int **imatrix(int nrl, int nrh, int ncl, int nch){
	/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch]
	and intializing the matrix, e.g. m[nrl..nrh][ncl..nch]=0.0 */
	int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	int **m;

	/* allocate pointers to rows */
	m = new int*[nrow + NR_END];
	//m = (int **)malloc((size_t)((nrow + NR_END)*sizeof(int)));
	if (!m) nrerror("allocation failure 1 in function imatrix() ");
	m += NR_END;
	m -= nrl;

	/* allocation rows and set pointers to them */
	m[nrl] = new int[nrow*ncol + NR_END];
	//m[nrl] = (int *)malloc((size_t)((nrow*ncol + NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in function matrix() ");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* initializing matrix */
	for (i = nrl; i <= nrh; i++)
		for (j = ncl; j <= nch; j++) m[i][j] = 0;

	/* return pointer to array of pointer to rows */
	return m;
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch){
	/* free a integer matrix allocated by imatrix() */
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}