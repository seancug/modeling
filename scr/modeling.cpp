// modeling with topography.cpp : Defines the entry point for the console application.
//
#include<stdlib.h>
#include"topoft.h"
#include"glo_var.h"

int main(int argc, char* argv[])
{
	//ʱ��������
	clock_t begin, end;
	double timeuse;
	begin = clock();
	cout << "-----------------------------------------------------------\n";
	cout << "------------------------PROGRAM BEGIN----------------------\n" << endl;//��ʼ���г���

	/* open log-file ? */
	string logfile;
	
	ofstream fpLog;
    ifstream fp("../par/Par.inp");//�ļ�ָ�룬��������ı�

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

	/*ģ�Ͳ�����������ⲿ�ļ�����*/
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

		if (vpIn == NULL)//����vp
		{
			cout << "read Vp model fails..." << endl;
			return -1;
		}
		if (vsIn == NULL)//����vs
		{
			cout << "read Vp model fails..." << endl;
			return -1;
		}
		if (denIn == NULL)//����vden
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
		//outputdata("..\\..\\temp\\slope_vp.temp", vp, 1, NX, 1, NZ, 4);/*��бģ�ͣ�ȷ�Ͻ�����ģ���Ƿ���ȷ*/
	}

	/*�жϽ��������ȶ���*/
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
	//����ʱ��������
	NT = iround(TIME / DT);       //ʱ��������
	ns = iround(NT / NDT);        //���������
	csnap = iround(TSNAP1 / DT);  //ʱ������㣬��һ������
	
	/* �����ٶȵ����ٶȳ���Ӧ�����ڴ�,��-1����ʼ���ǵ����ɵر��� */
	float  ** sxz, ** sxx, ** szz;
	float  ** vx, ** vz;
	vx = dmatrix(-5, NX + 6, -5, NZ + 6);
	vz = dmatrix(-5, NX + 6, -5, NZ + 6);
	sxz = dmatrix(-5, NX + 6, -5, NZ + 6);
	sxx = dmatrix(-5, NX + 6, -5, NZ + 6);
	szz = dmatrix(-5, NX + 6, -5, NZ + 6);
	
	/*ȷ������ϵ��*/
	float  ** absorb_coeff; //����ϵ��
	absorb_coeff = dmatrix(1, NX, 1, NZ);
	absorb(absorb_coeff, 1, NX, 1, NZ);
		//outputdata("../../temp/absorb.temp",absorb_coeff, 1, NX, 1, NZ, 5);
	
	/*����ر��������*/
	int **flag;
	flag = imatrix(1, NX, 1, NZ);
	int freeDH;
	if (FREE_SURF == "IMAGING_TOPOGRAPHY" || FREE_SURF == "AEA_TOPOGRAPHY")
	{
		flag_mark(flag, vpIn, 1, NX, 1, NZ, freeDH);
		outputdata("../model/flag.dat", flag, 1, NX, 1, NZ); //�����ǵ�����㣬ȷ��û�д���
	}

	if(READMOD == "YES")
		free_dmatrix(vpIn, 0, NX + 1, 1, NZ + 1);
	else if(READMOD == "NO")
		free_dmatrix(vpIn, -1, NX + 2, -1, NZ + 2);
	
	/*ȷ���첨��λ��*/
	int ntr; //�ܽ��յ���
	int    ** recpos = NULL; //�첨��λ��,recpos[1][...] xλ��,recpos[2][...],Zλ��
//	int **test = imatrix(1, NX, 1, NZ);
	recpos = receiver(ntr, flag); //���ؼ첨��λ������

	//outputdata("../../temp/recpos.temp", recpos, 1, 2, 1, ntr);

	/*��������¼,�������**  sectionvy,*/
	float  **  sectionvx, **  sectionvz;
	sectionvx = dmatrix(1, ntr, 1, ns);
	sectionvz = dmatrix(1, ntr, 1, ns);

	//float *cp_x, *cp_z;
	//cp_x = dvector(1, ns);
	//cp_z = dvector(1, ns);

	int csamp = NDT;//������ʵʱ�䣬�Ժ������NDT

	/*PML cofficient*/
	int ifw_x, ifw_z;
	ifw_x = iround(FW / DX);
	ifw_z = iround(FW / DZ);
	pml_c pml;
	/*���鿪�٣�pml�߽�*/
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
		/*���ɱ߽�����������*/
		if (FREE_SURF == "IMAGING_PLANER")/*ˮƽ���񣬹�����ȷ*/
		{
			/*����Ӧ��*/
			update_s(nt, vx, vz, sxx, szz, sxz, lam, u);

			/*ˮƽ���ɵر�*/
			surface(free_depth, vx, vz, sxx, sxz, szz, rden, lam, u);

			/*�����ٶ�*/
			update_v(nt, vx, vz, sxx, szz, sxz, rden, absorb_coeff);
		}
		else if (FREE_SURF == "IMAGING_TOPOGRAPHY")/*��������ر���Ҫ����*/
		{
			/*����Ӧ��*/
			topo_update_s(nt, flag, vx, vz, sxx, szz, sxz, lam, u);

			/*�����ٶ�*/
			topo_update_v(nt, flag, vx, vz, sxx, szz, sxz, rden, absorb_coeff);
		}
		else if (FREE_SURF == "AEA_PLANER")/*AEAˮƽ���ɵر�*/
		{
			/*����Ӧ��*/
			AEA_planer_update_s(nt, vx, vz, sxx, szz, sxz, lam, u);

			/*�����ٶ�*/
			AEA_planer_update_v(nt, vx, vz, sxx, szz, sxz, rden, absorb_coeff);
		}
		else if (FREE_SURF == "AEA_TOPOGRAPHY")/*AEA����ر�*/
		{
			/*����Ӧ��*/
			AEA_topo_update_s(nt, flag, vx, vz, sxx, szz, sxz, lam, u, pml, array, freeDH);

			/*�����ٶ�*/
			AEA_topo_update_v(nt, flag, vx, vz, sxx, szz, sxz, rden, pml, array, freeDH);
		}

		/*�����Դ*/
		add_sour(nt, vx, vz, rden);

		/*����첨�������¼*/
		if ((SEISMO != "NO") && ntr > 0 && csamp == nt && nt < NT)
		{
			seismogram(csamp, ntr, recpos, sectionvx, sectionvz, vx, vz);
			csamp += NDT;
		}


		static unsigned char w[] = "///-";
		if (p1++ == 3) p1 = 0;
		printf("\r\tfinish percent:%.1f %c %c", 100.0*(nt) / NT, '%', w[p1]);
		fflush(stdout);

		/*���沨������*/
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
	
	/*�ͷ��ٶ�Ӧ��*/
	free_dmatrix(vx, -5, NX + 6, -5, NZ + 6);
	free_dmatrix(vz, -5, NX + 6, -5, NZ + 6);
	free_dmatrix(sxx, -5, NX + 6, -5, NZ + 6);
	free_dmatrix(szz, -5, NX + 6, -5, NZ + 6);
	free_dmatrix(sxz, -5, NX + 6, -5, NZ + 6);
	/*�ͷŲ���*/
	free_dmatrix(rden, -1, NX + 2, -1, NZ + 2);
	free_dmatrix(lam, -1, NX + 2, -1, NZ + 2);
	free_dmatrix(u, -1, NX + 2, -1, NZ + 2);
	/*�ͷ�����ϵ��*/
	free_dmatrix(absorb_coeff, 1, NX, 1, NZ);
	/*�첨��*/
	free_imatrix(recpos, 1, 2, 1, ntr);
	free_dmatrix(sectionvx, 1, ntr, 1, ns);
	free_dmatrix(sectionvz, 1, ntr, 1, ns);
	/*���*/
	free_imatrix(flag, 1, NX, 1, NZ);
	/*�ͷ�pmlϵ��*/
	free_dvector(pml.left, 1, ifw_x);
	free_dvector(pml.right, NX - ifw_x + 1, NX);
	free_dvector(pml.bottom, NZ - ifw_z + 1, NZ);
	/*�ͷ�PML������*/
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

