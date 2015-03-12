/*  **  Copyright (C) 2014, Geoscience Of China University         **
**  All right reserved                                         **
**  Filename:modeling with topography                          **
**  Current Version:Debug 1.0                                  **
**  Author:Sun BO(Sean)                                        **
**  Finish Date:2014.7.24                                      **
**  For More Information,please contact seancug28@yahoo.com    ** */
#ifndef topoft_h
#define topoft_h
#define _TOPOHT_
#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
#include<cstring>
#include<string>
#include<vector>

/*宏定义*/
#define a1 (1.196289E0)
#define a2 (-7.9752604E-2)
#define a3 (9.5703125E-3)
#define a4 (-6.9754464E-4)
#define iround(x) ((int)(floor)((x)+0.5))
#define min(x,y) ((x<y)?x:y)
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)

#define PI (3.141592653589793)
#define NPAR 30
#define STRING_SIZE 74
#define REQUEST_COUNT 6

using namespace std;

struct pml_c
{
	float *left;
	float *right;
	float *bottom;

	float *leftv;
	float *rightv;
	float *bottomv;
};

struct pml_array
{
	float **left_p_vx;
	float **left_v_vx;
	float **left_p_vz;
	float **left_v_vz;
	float **left_p_sxx;
	float **left_v_sxx;
	float **left_p_szz;
	float **left_v_szz;
	float **left_p_sxz;
	float **left_v_sxz;

	float **right_p_vx;
	float **right_v_vx;
	float **right_p_vz;
	float **right_v_vz;
	float **right_p_sxx;
	float **right_v_sxx;
	float **right_p_szz;
	float **right_v_szz;
	float **right_p_sxz;
	float **right_v_sxz;

	float **bottom_p_vx;
	float **bottom_v_vx;
	float **bottom_p_vz;
	float **bottom_v_vz;
	float **bottom_p_sxx;
	float **bottom_v_sxx;
	float **bottom_p_szz;
	float **bottom_v_szz;
	float **bottom_p_sxz;
	float **bottom_v_sxz;
};
/*函数声明*/
void nrerror(char error_text[]);
void parerror(char error_text[]);
float **dmatrix(long nrl, long nrh, long ncl, long nch);
float *dvector(long nl, long nh);
void free_dmatrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(float *v, long nl, long nh);
string read_Par(ifstream &fp_in);
void check_par();
float* rd_sour(int *nts, float *nt, FILE* fp_source);
void model(float **rden, float **lam, float **u, int NX_f, int NX_e, int NZ_f, int NZ_e);
float** rd_model(char* filename, int &NX, int &NZ, float &DX, float &DZ);
void caculate_coff(float **rden, float **lam, float **u, float **vpIn, float **vsIn, float **rhoIn);
void absorb(float **absorb_coeff, int NX_f, int NX_e, int NZ_f, int NZ_e);
bool checkFD(string &message, float **lam, float **u, float **rden, int NX_f, int NX_e, int NZ_f, int NZ_e);
void update_v(int nt, float ** vx, float ** vz, float ** sxx, float ** szz, float ** sxz,
	float **  rden, float ** absorb_coeff);
void add_sour(int nt, float ** vx, float **vz, float **rden);
void update_s(int nt, float ** vx, float ** vz, float ** sxx, float ** szz, float ** sxz,
	float ** lam, float ** u);
int **imatrix(int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
int** receiver(int &ntr, int** flag);
void write_repos(int **repos, int ntr);
void seismogram(int lsamp, int ntr, int **recpos, float**sectionvx, float**sectionvz,
	float **vx, float **vz);
void writedsk(FILE *fp_out, float amp1, string format);
string snap(int nt, string format, string type,
	float **vx, float **vz, int idx, int idz, int nx1, int nz1, int nx2, int nz2);
string outseis(FILE *fpdata, float **section, int **recpos, int ntr, int ns, string seis_form);
string saveSeis(float **sectionvx, float **sectionvz, int  **recpos, int ntr, int ns);
void surface(int depth, float ** vx, float **vz, float ** sxx, float **sxz, float **szz,
	float ** rden, float **lam, float **u);
void flag_mark(int **flag, float **vp, int NX_f, int NX_e, int NZ_f, int NZ_e, int &depth);
void topo_update_s(int nt, int **flag, float ** vx, float ** vz, float ** sxx, float ** szz, float ** sxz,
	float ** lam, float ** u);
void topo_update_v(int nt, int **flag, float ** vx, float ** vz, float ** sxx, float ** szz, float ** sxz,
	float **  rden, float ** absorb_coeff);
void slope_model(float **rden, float **lam, float **u, float **vp, int NX_f, int NX_e, int NZ_f, int NZ_e);

void AEA_planer_update_s(int nt, float ** vx, float ** vz, float ** sxx, float ** szz, float ** sxz,
	float ** lam, float ** u);
void AEA_planer_update_v(int nt, float ** vx, float ** vz, float ** sxx, float ** szz, float ** sxz,
	float **  rden, float ** absorb_coeff);

void AEA_topo_update_v(int nt, int **flag, float ** vx, float ** vz, float ** sxx, float ** szz, float ** sxz,
	float **  rden, pml_c pml, pml_array array, int depth);
void AEA_topo_update_s(int nt, int **flag, float ** vx, float ** vz, float ** sxx, float ** szz, float ** sxz,
	float ** lam, float ** u, pml_c pml, pml_array array, int depth);
void pml_coff(pml_c pml, float Vs, int NX_f, int NX_e, int NZ_f, int NZ_e);

/*segy function*/
void Write_Dedault_Seg_Y_head(FILE *fpdata);

/*test function*/
void outputdata(char *s, float **datain, int nrl, int nrh, int ncl, int nch, int prec = 3);
void outputdata(char *s, int **datain, int nrl, int nrh, int ncl, int nch);

#endif
