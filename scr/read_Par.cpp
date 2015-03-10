#include <stdlib.h> 
#include<sstream>
#include"topoft.h"

string read_Par(ifstream &fp_in)
{
	extern float DX, DZ, TIME, DT, ANGLE, F0, XS, ZS;
	extern float TSNAP1, TSNAP2, TSNAPINC, FW;
	extern float XREC1, XREC2, ZREC1, ZREC2;
	extern float   DRX, DRZ;
	extern float  DAMPING;
	extern float AMP;
	extern int    NDT;
	extern int   NX, NZ;
	extern string SOURCE_WAVELET, SOURCE_TYP, READMOD, MD_FORMAT, FREE_SURF, FILEDS, SNAP_FORMAT, READREC, SEISMO, SEIS_FORMAT;
	extern char  SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE];
	extern string RECTYPE;
	extern char RECIN[STRING_SIZE];
	extern char   MFILE_VP[STRING_SIZE], MFILE_VS[STRING_SIZE], MFILE_DEN[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char  SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VZ[STRING_SIZE];
	extern string TITLE_NAME, VERBOSE;



	string message;
	
	char tem[252] = "";
	message.append("P U T  I N  P A R A M E T E R S\n");
	message.append("===============================\n");
	//默认值 
	TIME = 0.0;
	ANGLE = 0.0;
	AMP = 1.0;
	TSNAP1 = TSNAP2 = 0.1;
	TSNAPINC = 0.0;
	FW = 5.0;
	DAMPING = 1e-4;
	VERBOSE.assign("YES");
	TITLE_NAME.assign("Test");
	READMOD.assign("NO");
	MD_FORMAT.assign("GRID");
	FILEDS.assign("NO");
	READREC.assign("NO");
	SEISMO.assign("NO");
	SEIS_FORMAT.assign("SEGY");
	FREE_SURF.assign("AEA_TOPOGRAPHY");
	NDT = 1;
	READREC.assign("NO");
	RECTYPE.assign("TOPOGRAPHY");

	string line;
	string modelName;
	string comandLine;
	string parTxt, parNum;
	
	while (getline(fp_in, line, '\n'))
	{
		if (line[0] == '#' || line[0] == '\n' || line == "")//跳过空白注释行
		{
			line.clear();
			continue;
		}
		
		if (line[line.length() - 1] != '/')//命令是否换行
		{
			comandLine.append(line);
			continue;
		}

		comandLine.append(line);
		if (comandLine[0] == '&'&&comandLine[comandLine.size() - 1] == '/')
		{
			int i = 0;
			while (comandLine[i] != '/')
			{
				if (comandLine[i] == ' ') //去掉空格和换行符,制表符
				{
					char a = comandLine[i];
					comandLine.erase(i, 1);
					continue;
				}
				if (comandLine[i] == '\t')
				{
					comandLine.erase(i, 1);
					continue;
				}
				if (comandLine[i] == '\n')
				{
					comandLine.erase(i, 1);
					continue;
				}
				i++;
			}

			comandLine.erase(comandLine.size() - 1, 1);

            stringstream txt(comandLine);
			getline(txt, modelName, ':');
			modelName.erase(0, 1);//去除&

			if (modelName == "GRID")//GRID 参数读入
			{
				message.append("G R I D\n");
				message.append("=======\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "nx")
					{
						getline(txt, parNum, ',');
						NX = atoi(parNum.c_str());
                        sprintf(tem, "\tNumber of x Direction Grids.............<numbers> = %d\n", NX);
                        message.append(tem);
					}
					else if (parTxt == "nz")
					{
						getline(txt, parNum, ',');
						NZ = atoi(parNum.c_str());
                        sprintf(tem, "\tNumber of z Direction Grids.............<numbers> = %d\n", NZ);
                        message.append(tem);
					}
					else if (parTxt == "dx")
					{
						getline(txt, parNum, ',');
						DX = atof(parNum.c_str());
                        sprintf(tem, "\tInterval of X Direction.......................<m> = %.2f\n", DX);
                        message.append(tem);
					}
					else if (parTxt == "dz")
					{
						getline(txt, parNum, ',');
						DZ = atof(parNum.c_str());
                        sprintf(tem, "\tInterval of Z Direction.......................<m> = %.2f\n", DZ);
                        message.append(tem);
					}
				}
			}
			else if (modelName == "GENERAL")
            {
				message.append("\nG E N E R A L\n");
				message.append("==============\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "title")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						TITLE_NAME.assign(parNum);
                        sprintf(tem, "\tThe Title of the Prject.......................<s> = %s\n", TITLE_NAME.c_str());
                        message.append(tem);
					}
					else if (parTxt == "verbose")
                    {
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						VERBOSE.assign(parNum);
                        sprintf(tem, "\tThe VERBOSE of the Prject.....................<s> = %s\n", TITLE_NAME.c_str());
                        message.append(tem);
					}
				}
                 printf("4\n");
			}
			else if (modelName == "TIME")
			{
				message.append("\nT I M E\n");
				message.append("========\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "ti")
					{
						getline(txt, parNum, ',');
						TIME = atof(parNum.c_str());
						sprintf(tem, "\tThe time of Wave Propagation..................<s> = %.2f\n", TIME );
						message.append(tem);
					}
					else if (parTxt == "dt")
					{
						getline(txt, parNum, ',');
						DT = atof(parNum.c_str());
						sprintf(tem, "\tThe time Interval of Wave Propagation.........<s> = %.2e\n", DT);
						message.append(tem);
					}
				}
			}
			else if (modelName == "SOURCE")
			{
				message.append("\nS O U R C E\n");
				message.append("============\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "wavelet")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						SOURCE_WAVELET.assign(parNum);
						sprintf(tem, "\tSource Wavelet.............................<name> = %s\n", SOURCE_WAVELET.c_str());
						message.append(tem);
					}
					else if (parTxt == "fileName")
					{
						getline(txt, parNum, ',');
						strcpy(SOURCE_FILE, parNum.c_str());
						sprintf(tem, "\tIn Put Source Wavelet File.................<name> = %s\n", SOURCE_FILE);
						message.append(tem);

					}
					else if (parTxt == "mechanism")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						SOURCE_TYP.assign(parNum);
						sprintf(tem, "\tThe Mechanism of Source....................<name> = %s\n", SOURCE_TYP.c_str());
						message.append(tem);
					}
					else if (parTxt == "angle")
					{
						getline(txt, parNum, ',');
						ANGLE = atof(parNum.c_str()) / 180.0*PI;
						sprintf(tem, "\tThe Angele of Point source.................<name> = %.2f\n", ANGLE);
						message.append(tem);
					}
					else if (parTxt == "f0")
					{
						getline(txt, parNum, ',');
						F0 = atof(parNum.c_str());
						sprintf(tem, "\tThe Peak frequency of source...............<name> = %.2f\n", F0);
						message.append(tem);
					}
					else if (parTxt == "coorid")
					{
						getline(txt, parNum, ',');
						XS = atof(parNum.c_str());
						getline(txt, parNum, ',');
						ZS = atof(parNum.c_str());
						sprintf(tem, "\tThe coorinate (x,z) of source.................<m> = (%.2f,%.2f)\n", XS, ZS);
						message.append(tem);
					}
					else if (parTxt == "ampli")
					{
						getline(txt, parNum, ',');
						AMP = atof(parNum.c_str());
						sprintf(tem, "\tThe Peak Amplitude of source...............<name> = %.2f\n", AMP);
						message.append(tem);
					}
				}
			}
			else if (modelName == "MODEL")
			{
				message.append("\nM O D E L\n");
				message.append("==========\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "readMod")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						READMOD = parNum;
						sprintf(tem, "\tRead Model from the file......................... = %s\n", READMOD.c_str());
						message.append(tem);
					}
					else if (parTxt == "mfile_vp")
					{
						getline(txt, parNum, ',');
						strcpy(MFILE_VP, parNum.c_str());
						sprintf(tem, "\tThe Path of Vp Model File..................<name> = %s\n", MFILE_VP);
						message.append(tem);
					}
					else if (parTxt == "mfile_vs")
					{
						getline(txt, parNum, ',');
						strcpy(MFILE_VS, parNum.c_str());
						sprintf(tem, "\tThe Path of VS Model File..................<name> = %s\n", MFILE_VS);
						message.append(tem);
					}
					else if (parTxt == "mfile_den")
					{
						getline(txt, parNum, ',');
						strcpy(MFILE_DEN, parNum.c_str());
						sprintf(tem, "\tThe Path of Density Model File.............<name> = %s\n", MFILE_DEN);
						message.append(tem);
					}
					else if (parTxt == "format")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						MD_FORMAT.assign(parNum);
						sprintf(tem, "\tThe Format of Model File...................<name> = %s\n", MD_FORMAT.c_str());
						message.append(tem);
					}
				}
			}
			else if (modelName == "SURFACE")
			{
				message.append("\nS U R F A C E\n");
				message.append("==============\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "free_surf")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						FREE_SURF.assign(parNum);
						sprintf(tem, "\tThe Taple of Free Surface..................<name> = %s\n", FREE_SURF.c_str());
						message.append(tem);
					}
				}
			}
			else if (modelName == "BC_PAR")
			{
				message.append("\nB C _ P A R\n");
				message.append("============\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "fw")
					{
						getline(txt, parNum, ',');
						FW = atof(parNum.c_str());
						if (FW < 0.0) FW = 0.0;
						sprintf(tem, "\tThe distance of PML Layer.....................<m> = %.2f\n", FW);
						message.append(tem);
					}
					else if (parTxt == "damping")
					{
						getline(txt, parNum, ',');
						DAMPING = atof(parNum.c_str());
						sprintf(tem, "\tThe Reflect cofficient........................... = %.2e\n", DAMPING);
						message.append(tem);
					}
				}
			}
			else if (modelName == "SNAP_DEF")
			{
				message.append("\nS N A P _ D E F\n");
				message.append("================\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "fileds")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						FILEDS.assign(parNum);
						sprintf(tem, "\tThe Fileds of Snap shot.......................... = %s\n", FILEDS.c_str());
						message.append(tem);
					}
					else if (parTxt == "format")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						SNAP_FORMAT.assign(parNum);
						sprintf(tem, "\tThe Format of Snap Shot File...............<name> = %s\n", MD_FORMAT.c_str());
						message.append(tem);
					}
					else if (parTxt == "tsnap1")
					{
						getline(txt, parNum, ',');
						TSNAP1 = atof(parNum.c_str());
						sprintf(tem, "\tThe Begin of Snap Shot........................<s> = %.2f\n", TSNAP1);
						message.append(tem);
					}
					else if (parTxt == "tsnap2")
					{
						getline(txt, parNum, ',');
						TSNAP2 = atof(parNum.c_str());
						sprintf(tem, "\tThe End of Snap Shot..........................<s> = %.2f\n", TSNAP2);
						message.append(tem);
					}
					else if (parTxt == "tsnapinc")
					{
						getline(txt, parNum, ',');
						TSNAPINC = atof(parNum.c_str());
						sprintf(tem, "\tThe Interval of Snap Shot.....................<s> = %.2f\n", TSNAPINC);
						message.append(tem);
					}
					else if (parTxt == "snapFile")
					{
						getline(txt, parNum, ',');
						strcpy(SNAP_FILE, parNum.c_str());
						sprintf(tem, "\tThe path of Save Snap Shot File............<name> = %s\n", SNAP_FILE);
						message.append(tem);
					}
				}
			}
			else if (modelName == "REC_DEF")
			{
				message.append("\nR E C _ D E F\n");
				message.append("==============\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "readRec")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						READREC.assign(parNum);
						sprintf(tem, "\tRead the Recevers Position form the File...<name> = %s\n", READREC.c_str());
						message.append(tem);
					}
					else if (parTxt == "rec_file")
					{
						getline(txt, parNum, ',');
						strcpy(RECIN, parNum.c_str());
						sprintf(tem, "\tThe Path of the Recevers Position........ .<name> = %s\n", RECIN);
						message.append(tem);
					}
				}
			}
			else if (modelName == "REC_LINE")
			{
				message.append("\nR E C _ L I N E\n");
				message.append("================\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "rec_type")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						RECTYPE.assign(parNum);
						sprintf(tem, "\tThe Recevers Type................................ = %s\n", RECTYPE.c_str());
						message.append(tem);
					}
					else if (parTxt == "first")
					{
						getline(txt, parNum, ',');
						XREC1 = atof(parNum.c_str());
						getline(txt, parNum, ',');
						ZREC1 = atof(parNum.c_str());
						sprintf(tem, "\tThe First Recever(x,z) of source..............<m> = (%.2f,%.2f)\n", XREC1, ZREC1);
						message.append(tem);
					}
					else if (parTxt == "last")
					{
						getline(txt, parNum, ',');
						XREC2 = atof(parNum.c_str());
						getline(txt, parNum, ',');
						ZREC2 = atof(parNum.c_str());
						sprintf(tem, "\tThe Last Recever(x,z) of source...............<m> = (%.2f,%.2f)\n", XREC2, ZREC2);
						message.append(tem);
					}
					else if (parTxt == "interval")
					{
						getline(txt, parNum, ',');
						DRX = atof(parNum.c_str());
						getline(txt, parNum, ',');
						DRZ = atof(parNum.c_str());
						sprintf(tem, "\tThe Interval of Recever(dx,dz) of source......<m> = (%.2f,%.2f)\n", DRX, DRZ);
						message.append(tem);
					}
					else if (parTxt == "out_rec")
					{
						getline(txt, parNum, ',');
						strcpy(REC_FILE, parNum.c_str());
						sprintf(tem, "\tThe path of Out Put Receivers..............<name> = %s\n", REC_FILE);
						message.append(tem);
					}
				}
			}
			else if (modelName == "GATHER")
			{
				message.append("\nG A T H E R\n");
				message.append("============\n");
				while (!txt.eof())
				{
					getline(txt, parTxt, '=');
					if (parTxt == "seisGather")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						SEISMO.assign(parNum);
						sprintf(tem, "\tThe Fileds of Seismic Gather..................... = %s\n", SEISMO.c_str());
						message.append(tem);
					}
					else if (parTxt == "ndt")
					{
						getline(txt, parNum, ',');
						NDT = atoi(parNum.c_str());
						sprintf(tem, "\tThe resample of time.....................<number> = %d\n", NDT);
						message.append(tem);
					}
					else if (parTxt == "format")
					{
						getline(txt, parNum, ',');
						parNum.erase(0, 1);
						parNum.erase(parNum.length() - 1, 1);
						SEIS_FORMAT.assign(parNum);
						sprintf(tem, "\tThe Format of Seismic Gather...............<name> = %s\n", SEIS_FORMAT.c_str());
						message.append(tem);
					}
					else if (parTxt == "seis_file_vx")
					{
						getline(txt, parNum, ',');
						strcpy(SEIS_FILE_VX, parNum.c_str());
						sprintf(tem, "\tThe path of Save vx File...................<name> = %s\n", SEIS_FILE_VX);
						message.append(tem);
					}
					else if (parTxt == "seis_file_vz")
					{
						getline(txt, parNum, ',');
						strcpy(SEIS_FILE_VZ, parNum.c_str());
						sprintf(tem, "\tThe path of Save vz File...................<name> = %s\n", SEIS_FILE_VZ);
						message.append(tem);
					}
				}
			}
			else
			{
				//message.append("module name: &"+modelName+" does't exit!\n");
				cout << "\n\nmodule name: &" << modelName << " does't exit!\n";
			}
		}
		comandLine.clear();
	}
	line.clear();
	return message;
//	fprintf(fp_log, "Parameters read completed...\n");
}

		
		
		
	
