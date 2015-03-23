// modeling with topography.cpp : Defines the entry point for the console application.
//
#include<stdlib.h>
#include"topoft.h"
#include"glo_var.h"

int main(int argc, char* argv[])
{
	//时间计算操作
	clock_t begin, end;
	double timeuse;
	begin = clock();
	cout << "-----------------------------------------------------------\n";
	cout << "------------------------PROGRAM BEGIN----------------------\n" << endl;//开始运行程序

	/* open log-file ? */
	string logfile;
	
	ofstream fpLog;
    ifstream fp("../par/Par.inp");//文件指针，保存参数文本

	if (!fp.fail())
	{
		string message;
		message = read_Par(fp);  
        logfile = "../logs/"+TITLE_NAME + VERBOSE + ".log";
		fpLog.open(logfile.c_str()); 
		cout << message;
		fpLog << message;
		check_par();
	}
	else
    {
        fpLog.open("../logs/fd.log");
		fpLog << "open file error[File maybe doesn't exit]!\n" << endl;
        cout << "open file error[File maybe doesn't exit]!" << endl;
		exit(-1);
	}

	/*模型参数建立或从外部文件输入*/
	float **rden = NULL, **lam = NULL, **u = NULL;
	float **denIn = NULL; float **vsIn = NULL, **vpIn = NULL;
	if (READMOD == "YES")
	{
        char tem[124];
        sprintf(tem, "%s%s","../model/", MFILE_VP);
        vpIn = rd_model(tem, NX, NZ, DX, DZ);
        sprintf(tem, "%s%s","../model/", MFILE_VS);
        vsIn = rd_model(tem, NX, NZ, DX, DZ);
        sprintf(tem, "%s%s","../model/", MFILE_DEN);
        denIn = rd_model(tem, NX, NZ, DX, DZ);

        char temmess[128];
        string message;
        message.append("G R I D\n");
        message.append("=======\n");
        sprintf(temmess, "\tNumber of x Direction Grids.............<numbers> = %d\n", NX);
        message.append(temmess);
        sprintf(temmess, "\tNumber of z Direction Grids.............<numbers> = %d\n", NZ);
        message.append(temmess);
        sprintf(temmess, "\tInterval of X Direction.......................<m> = %.2f\n", DX);
        message.append(temmess);
        sprintf(temmess, "\tInterval of Z Direction.......................<m> = %.2f\n", DZ);
        message.append(temmess);

        cout<<message;
        fpLog<<message;

		if (vpIn == NULL)//读入vp
		{
			cout << "read Vp model fails..." << endl;
			return -1;
		}
		if (vsIn == NULL)//读入vs
		{
			cout << "read Vp model fails..." << endl;
			return -1;
		}
		if (denIn == NULL)//读入vden
		{
			cout << "read Vp model fails..." << endl;
			return -1;
		}
		rden = dmatrix(-1, NX + 2, -1, NZ + 2);
		lam = dmatrix(-1, NX + 2, -1, NZ + 2);
		u = dmatrix(-1, NX + 2, -1, NZ + 2);
		caculate_coff(rden, lam, u, vpIn, vsIn, denIn);
		free_dmatrix(vsIn, 0, NX + 1, 1, NZ + 1);
		free_dmatrix(denIn, 0, NX + 1, 1, NZ + 1);
	}
	else if (READMOD == "NO")
	{
		rden = dmatrix(-1, NX + 2, -1, NZ + 2);
		lam = dmatrix(-1, NX + 2, -1, NZ + 2);
		u = dmatrix(-1, NX + 2, -1, NZ + 2);
		vpIn = dmatrix(-1, NX + 2, -1, NZ + 2);
		//model(rden, lam, u, 1, NX, 1, NZ);
		//outputdata("../../temp/rden.temp", rden, -1, NX + 2, -1, NZ + 2, 4);
		/*outputdata("../../temp/lam.temp", lam, -1, NX + 2, -1, NZ + 2, 4);
		outputdata("../../temp/u.temp", u, -1, NX + 2, -1, NZ + 2, 4);*/

		slope_model(rden, lam, u, vpIn, 1, NX, 1, NZ);
        printf("1");
		//outputdata("..\\..\\temp\\slope_rden.temp", rden, 1, NX, 1, NZ, 4);
		//outputdata("..\\..\\temp\\slope_vp.temp", vp, 1, NX, 1, NZ, 4);/*倾斜模型，确认建立的模型是否正确*/
	}

	/*判断交错网格稳定性*/
	string message;
	if (!checkFD(message, lam, u, rden, 1, NX, 1, NZ))
	{
		fpLog << message << endl;
		cout << message << endl;
		fpLog.close();
		exit(0);
	}
	fpLog << message << endl;
	cout << message << endl;
	message.clear();

	int ns, csnap;
	//计算时间网格数
	NT = iround(TIME / DT);       //时间网格数
	ns = iround(NT / NDT);        //地震采样数
	csnap = iround(TSNAP1 / DT);  //时间网格点，第一个快照
	
	/* 分配速度迭代速度场和应力场内存,从-1，开始考虑到自由地表镜像法 */
	float  ** sxz, ** sxx, ** szz;
	float  ** vx, ** vz;
	vx = dmatrix(-5, NX + 6, -5, NZ + 6);
	vz = dmatrix(-5, NX + 6, -5, NZ + 6);
	sxz = dmatrix(-5, NX + 6, -5, NZ + 6);
	sxx = dmatrix(-5, NX + 6, -5, NZ + 6);
	szz = dmatrix(-5, NX + 6, -5, NZ + 6);
	
	/*确定吸收系数*/
	float  ** absorb_coeff; //吸收系数
	absorb_coeff = dmatrix(1, NX, 1, NZ);
	absorb(absorb_coeff, 1, NX, 1, NZ);
		//outputdata("../../temp/absorb.temp",absorb_coeff, 1, NX, 1, NZ, 5);
	
	/*起伏地表网格分类*/
	int **flag;
	flag = imatrix(1, NX, 1, NZ);
	int freeDH;
	if (FREE_SURF == "IMAGING_TOPOGRAPHY" || FREE_SURF == "AEA_TOPOGRAPHY")
	{
		flag_mark(flag, vpIn, 1, NX, 1, NZ, freeDH);
		outputdata("../model/flag.dat", flag, 1, NX, 1, NZ); //输出标记的网格点，确认没有错误
	}

	if(READMOD == "YES")
		free_dmatrix(vpIn, 0, NX + 1, 1, NZ + 1);
	else if(READMOD == "NO")
		free_dmatrix(vpIn, -1, NX + 2, -1, NZ + 2);
	
	/*确定检波器位置*/
	int ntr; //总接收道数
	int    ** recpos = NULL; //检波器位置,recpos[1][...] x位置,recpos[2][...],Z位置
//	int **test = imatrix(1, NX, 1, NZ);
	recpos = receiver(ntr, flag); //返回检波器位置数组

	//outputdata("../../temp/recpos.temp", recpos, 1, 2, 1, ntr);

	/*输出地震记录,数组分配**  sectionvy,*/
	float  **  sectionvx, **  sectionvz;
	sectionvx = dmatrix(1, ntr, 1, ns);
	sectionvz = dmatrix(1, ntr, 1, ns);

	//float *cp_x, *cp_z;
	//cp_x = dvector(1, ns);
	//cp_z = dvector(1, ns);

	int csamp = NDT;//采样其实时间，以后会增加NDT

	/*PML cofficient*/
	int ifw_x, ifw_z;
	ifw_x = iround(FW / DX);
	ifw_z = iround(FW / DZ);
	pml_c pml;
	/*数组开辟，pml边界*/
	pml.left = dvector(1, ifw_x);
	pml.right = dvector(NX - ifw_x + 1, NX);
	pml.bottom = dvector(NZ - ifw_z + 1, NZ);
	pml.leftv = dvector(1, ifw_x);
	pml.rightv = dvector(NX - ifw_x + 1, NX);
	pml.bottomv = dvector(NZ - ifw_z + 1, NZ);
    pml_coff(pml, 1732.05, 1, NX, 1, NZ);

	pml_array array;

	array.left_p_vx = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_v_vx = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_p_vz = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_v_vz = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_p_sxx = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_v_sxx = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_p_szz = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_v_szz = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_p_sxz = dmatrix(1, ifw_x, 1, NZ - ifw_z);
	array.left_v_sxz = dmatrix(1, ifw_x, 1, NZ - ifw_z);

	array.right_p_vx = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_v_vx = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_p_vz = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_v_vz = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_p_sxx = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_v_sxx = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_p_szz = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_v_szz = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_p_sxz = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	array.right_v_sxz = dmatrix(NX - ifw_x + 1, NX, 1, NZ - ifw_z);

	array.bottom_p_vx = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_v_vx = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_p_vz = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_v_vz = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_p_sxx = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_v_sxx = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_p_szz = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_v_szz = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_p_sxz = dmatrix(1, NX, NZ - ifw_z + 1, NZ);
	array.bottom_v_sxz = dmatrix(1, NX, NZ - ifw_z + 1, NZ);

	
	extern int free_depth;
	unsigned short int p1 = 0;
	printf("\n");

	cout << "--------------------------" << endl;
	cout << "Update Velocity and Stress" << endl;
	cout << "--------------------------" << endl;
	for (int nt = 1; nt <= NT; nt++)
	{
		/*自由边界条件，分类*/
		if (FREE_SURF == "IMAGING_PLANER")/*水平镜像，功能正确*/
		{
			/*更新应力*/
			update_s(nt, vx, vz, sxx, szz, sxz, lam, u);

			/*水平自由地表*/
			surface(free_depth, vx, vz, sxx, sxz, szz, rden, lam, u);

			/*更新速度*/
			update_v(nt, vx, vz, sxx, szz, sxz, rden, absorb_coeff);
		}
		else if (FREE_SURF == "IMAGING_TOPOGRAPHY")/*镜像法起伏地表，需要调试*/
		{
			/*更新应力*/
			topo_update_s(nt, flag, vx, vz, sxx, szz, sxz, lam, u);

			/*更新速度*/
			topo_update_v(nt, flag, vx, vz, sxx, szz, sxz, rden, absorb_coeff);
		}
		else if (FREE_SURF == "AEA_PLANER")/*AEA水平自由地表*/
		{
			/*更新应力*/
			AEA_planer_update_s(nt, vx, vz, sxx, szz, sxz, lam, u);

			/*更新速度*/
			AEA_planer_update_v(nt, vx, vz, sxx, szz, sxz, rden, absorb_coeff);
		}
		else if (FREE_SURF == "AEA_TOPOGRAPHY")/*AEA起伏地表*/
		{
			/*更新应力*/
			AEA_topo_update_s(nt, flag, vx, vz, sxx, szz, sxz, lam, u, pml, array, freeDH);

			/*更新速度*/
			AEA_topo_update_v(nt, flag, vx, vz, sxx, szz, sxz, rden, pml, array, freeDH);
		}

		/*添加震源*/
		add_sour(nt, vx, vz, rden);

		/*保存检波器地震记录*/
		if ((SEISMO != "NO") && ntr > 0 && csamp == nt && nt < NT)
		{
			seismogram(csamp, ntr, recpos, sectionvx, sectionvz, vx, vz);
			csamp += NDT;
		}


		static unsigned char w[] = "///-";
		if (p1++ == 3) p1 = 0;
		printf("\r\tfinish percent:%.1f %c %c", 100.0*(nt) / NT, '%', w[p1]);
		fflush(stdout);

		/*保存波场快照*/
		if ((FILEDS != "NO") && nt == csnap && (nt <= TSNAP2 / DT))
		{
			message = snap(nt, SNAP_FORMAT, FILEDS, vx, vz, 1, 1, 1, 1, NX, NZ);
			fpLog << message << endl;
			//cout << message << endl;
			message.clear();
			csnap += iround(TSNAPINC / DT);
		}
	}

	cout << endl;
	if (ntr > 0 && (SEISMO != "NO"))
	{
		message = saveSeis(sectionvx, sectionvz, recpos, ntr, ns);
		fpLog << message << endl;
		cout << message << endl;
		message.clear();
	}
	//ofstream outcp_x("cp_vx.txt");
	//ofstream outcp_z("cp_vz.txt");
	//for (int j = 1; j <= ns; j++)
	//{
	//	outcp_x << cp_x[j] << endl;
	//	outcp_z << cp_z[j] << endl;
	//}
	//outcp_x.close();
	//outcp_z.close();
	//free_dvector(cp_x, 1, ns);
	//free_dvector(cp_z, 1, ns);
	
	/*释放速度应力*/
	free_dmatrix(vx, -5, NX + 6, -5, NZ + 6);
	free_dmatrix(vz, -5, NX + 6, -5, NZ + 6);
	free_dmatrix(sxx, -5, NX + 6, -5, NZ + 6);
	free_dmatrix(szz, -5, NX + 6, -5, NZ + 6);
	free_dmatrix(sxz, -5, NX + 6, -5, NZ + 6);
	/*释放参数*/
	free_dmatrix(rden, -1, NX + 2, -1, NZ + 2);
	free_dmatrix(lam, -1, NX + 2, -1, NZ + 2);
	free_dmatrix(u, -1, NX + 2, -1, NZ + 2);
	/*释放吸收系数*/
	free_dmatrix(absorb_coeff, 1, NX, 1, NZ);
	/*检波器*/
	free_imatrix(recpos, 1, 2, 1, ntr);
	free_dmatrix(sectionvx, 1, ntr, 1, ns);
	free_dmatrix(sectionvz, 1, ntr, 1, ns);
	/*标记*/
	free_imatrix(flag, 1, NX, 1, NZ);
	/*释放pml系数*/
	free_dvector(pml.left, 1, ifw_x);
	free_dvector(pml.right, NX - ifw_x + 1, NX);
	free_dvector(pml.bottom, NZ - ifw_z + 1, NZ);
	/*释放PML层数组*/
	free_dmatrix(array.left_p_vx, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_v_vx, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_p_vz, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_v_vz, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_p_sxx, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_v_sxx, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_p_szz, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_v_szz, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_p_sxz, 1, ifw_x, 1, NZ - ifw_z);
	free_dmatrix(array.left_v_sxz, 1, ifw_x, 1, NZ - ifw_z);

	free_dmatrix(array.right_p_vx, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_v_vx, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_p_vz, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_v_vz, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_p_sxx, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_v_sxx, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_p_szz, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_v_szz, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_p_sxz, NX - ifw_x + 1, NX, 1, NZ - ifw_z);
	free_dmatrix(array.right_v_sxz, NX - ifw_x + 1, NX, 1, NZ - ifw_z);

	free_dmatrix(array.bottom_p_vx, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_v_vx, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_p_vz, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_v_vz, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_p_sxx, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_v_sxx, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_p_szz, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_v_szz, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_p_sxz, 1, NX, NZ - ifw_z + 1, NZ);
	free_dmatrix(array.bottom_v_sxz, 1, NX, NZ - ifw_z + 1, NZ);
	end = clock();
	timeuse = (double)(end - begin) / CLOCKS_PER_SEC;
	int hour, min, sec;
	hour = (int)(timeuse / 60.0 / 60.0);
	min = (int)((timeuse - 60.0*60.0*hour) / 60.0);
	sec = (int)(timeuse - 60.0*60.0*hour - min * 60.0);
	timeuse = timeuse / 60.0;
	fpLog << "\nCaculate time:" << hour << "h " << min << "min " << sec << "s " << endl;
	cout << "\nCaculate time:" << hour << "h " << min << "min " << sec << "s " << endl;
	fpLog.close();

//	system("pause");
	return 0;
}

