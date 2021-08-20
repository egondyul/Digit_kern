#include "Elasticity2Da.h"

//------Without PML---------//

Elasticity2Da::Elasticity2Da()
{
	//initializing geometry 
	readData_geometry(); //Nx,Nz, hx, hz, tau, ...


	//memory alloc
	rho_x = new realval[(Nx - 1)*Nz];
	rho_z = new realval[Nx*(Nz - 1)];
	//искомые величины на текущем шаге по времени
	sigma_xx = new realval[Nx*Nz];					// нормальное напряжение упругого материала по x
	sigma_zz = new realval[Nx*Nz];					// нормальное напряжение упругого материала по z
	sigma_xz = new realval[(Nx - 1)*(Nz - 1)];		// касательное напряжение упругого материала
	v_x = new realval[(Nx - 1)*Nz];					// скорость твердх частиц по x
	v_z = new realval[Nx*(Nz - 1)];					// скорость твердых частиц по z

	c11 = new realval[Nx*Nz];
	c13 = new realval[Nx*Nz];
	c33 = new realval[Nx*Nz];
	c55_xz = new realval[(Nx - 1)*(Nz - 1)];


	//for viscoelasticity
	/*R_xx_new = new std::vector<realval>;
	(*R_xx_new).resize(Nx*Nz);
	R_xx_old = new std::vector<realval>;
	(*R_xx_old).resize(Nx*Nz);
	R_zz_new = new std::vector<realval>;
	(*R_zz_new).resize(Nx*Nz);
	R_zz_old = new std::vector<realval>;
	(*R_zz_old).resize(Nx*Nz);
	R_xz_new = new std::vector<realval>;
	(*R_xz_new).resize((Nx-1)*(Nz-1));
	R_xz_old = new std::vector<realval>;
	(*R_xz_old).resize((Nx-1)*(Nz-1));*/

	Rxx = new realval[Nx*Nz];
	Rzz = new realval[Nx*Nz];
	Rxz = new realval[(Nx - 1)*(Nz - 1)];

	/*ofstream fileR;
	strcpy(path1,"R_xz_"); // не для MPI

	std::ofstream outphi;
	outphi.open("R_xz_" + std::to_string(0) + ".txt");
	if (outphi.is_open())
	{
		for (int it = 0; it < Nz-1; it++)
		{
			for (int jt = 0; jt < Nx-1; jt++)
			{
				outphi << (*R_xz_new)[it*(Nx-1)+jt] << " ";
			}
			outphi << std::endl;
		}
	}
	outphi.close();*/


	tau11 = new realval[Nx*Nz];
	tau13 = new realval[Nx*Nz];
	tau33 = new realval[Nx*Nz];
	tau55_xz = new realval[(Nx - 1)*(Nz - 1)];

	src = new realval[Nx*Nz];

	//initialize all parameters: C, tau_tensor, rho
	readData_parameters();

	rx = tau / hx;
	rz = tau / hz;

	aa = tau_sigma / tau;
	bbb = (realval(2.0) * aa - realval(1.0)) / (realval(2.0) * aa + realval(1.0));
	//bbb = (aa - 0.5) / (aa + 0.5);
	//ccc = 1 / (aa + 0.5)/tau;
	ccc = 1 / (tau_sigma + 0.5*tau);

	//for rendering (output)
	coordinate();

	// Initial data (zero)

	int i, j;
#pragma omp parallel num_threads(4)
	{
#pragma omp for private(i,j) schedule(guided)
		for (j = 0; j < Nz; j++)
		{
			for (i = 0; i < Nx; i++)
			{
				sigma_xx[j*Nx + i] = 0.0;
				sigma_zz[j*Nx + i] = 0.0;
				src[j*Nx + i] = 0.0;
				Rxx[j*Nx + i] = 0.0;
				Rzz[j*Nx + i] = 0.0;
				if (i < Nx - 1)
				{
					v_x[j*(Nx - 1) + i] = 0.0;
				}
				if (j < Nz - 1)
				{
					v_z[j*Nx + i] = 0.0;
				}
				if ((j < Nz - 1) && (i < Nx - 1))
				{
					sigma_xz[j*(Nx - 1) + i] = 0.0;
					Rxz[j*(Nx - 1) + i] = 0.0;
				}
			}
		}
	}

	// set source points
	for (i = 0; i < Nx; i++)
	{
		src[inpData.z_src*Nx + i] = tau;
	}

	//mu, lambda
/*#pragma omp parallel
	{
#pragma omp for private(i,j) schedule(guided)
		for (j = 0; j < Nz; j++)
		{
			for (i = 0; i < Nx; i++)
			{
				mu[j*Nx + i] = rho[j*Nx + i] * Vs[j*Nx + i] * Vs[j*Nx + i];
				lam[j*Nx + i] = rho[j*Nx + i] * Vp[j*Nx + i] * Vp[j*Nx + i] - 2 * mu[j*Nx + i];
			}
		}
	}*/

	cur_time = 0.0;
	cur_time_str = 0.0;
}

Elasticity2Da::~Elasticity2Da()
{
	delete[]sigma_xx;
	delete[]sigma_xz;
	delete[]sigma_zz;
	delete[]v_x;
	delete[]v_z;
	delete[]src;

	delete[]rho_x;
	delete[] rho_z;
	delete[] c11;
	delete[] c13;
	delete[] c33;
	delete[] c55_xz;

	delete[] tau11;
	delete[] tau13;
	delete[] tau33;
	delete[] tau55_xz;

	delete[]Rxx;
	delete[]Rzz;
	delete[]Rxz;
}

void Elasticity2Da::readData_geometry()
{
	ifstream parameters;
	strcpy(path1, "INPUT.txt");
	strcpy(FullPath,  path);
	strcat(FullPath,path1);
	parameters.open(FullPath, ios::in);

	if (parameters.is_open())
	{
		double tmp; //smt i don't need
		parameters >> inpData.v0;
		cout << "v0 = " << inpData.v0 << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.Time;
		cout << "Time = " << inpData.Time << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.srcplace;
		cout << "srcplace = " << inpData.srcplace << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.z_src;
		cout << "z_src = " << inpData.z_src << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.z_rec_1;
		cout << "z_rec_1 = " << inpData.z_rec_1 << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.z_rec_2;
		cout << "z_rec_2 = " << inpData.z_rec_2 << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.Nz_PML_l;
		cout << "Nz_PML_l = " << inpData.Nz_PML_l << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.Nz_PML_r;
		cout << "Nz_PML_r = " << inpData.Nz_PML_r << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_xx_rec;
		cout << "sigma_xx_rec = " << inpData.sigma_xx_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_zz_rec;
		cout << "sigma_zz_rec = " << inpData.sigma_zz_rec << "\n";  
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_xz_rec;
		cout << "sigma_xz_rec = " << inpData.sigma_xz_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.v_x_rec;
		cout << "v_x_rec = " << inpData.v_x_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.v_z_rec;
		cout << "v_z_rec = " << inpData.v_z_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_xx_snap;
		cout << "sigma_xx_snap = " << inpData.sigma_xx_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_zz_snap;
		cout << "sigma_zz_snap = " << inpData.sigma_zz_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_xz_snap;
		cout << "sigma_xz_snap = " << inpData.sigma_xz_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.v_x_snap;
		cout << "v_x_snap = " << inpData.v_x_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.v_z_snap;
		cout << "v_z_snap = " << inpData.v_z_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.NSnaps;
		cout << "NSnaps = " << inpData.NSnaps << "\n";
		parameters.ignore(1000, '\n');
		parameters >> tmp;
		parameters.ignore(1000, '\n');
		parameters >> Nx;
		cout << "Nx = " << Nx << "\n";
		parameters.ignore(1000, '\n');
		parameters >> Nz;
		cout << "Nz = " << Nz << "\n";
		parameters.ignore(1000, '\n');
	}
	else
	{
		std::cout << "INPUT.txt failed." << std::endl;
	}
	parameters.close();

	inpData.Time = inpData.Time + 3.0 / inpData.v0; // add time for signal to form

	ifstream grid;
	strcpy(path1,"grid.bin"); // не для MPI
	strcpy(FullPath,  path);
	strcat(FullPath,  path1);
	grid.open(FullPath, ios::binary | ios::in);

	if (grid.is_open())
	{
		//grid.read((char*)&Nx, sizeof(int));
		//grid.read((char*)&Nz, sizeof(int));
		grid.read((char*)&hx, sizeof(realval));
		grid.read((char*)&hz, sizeof(realval));
		grid.read((char*)&tau, sizeof(realval));
	}
	else
	{
		std::cout << "grid.bin failed." << std::endl;
	}
	grid.close();

	cout << "Nx=" << Nx << endl;
	cout << "Nz=" << Nz << endl;

	//fflush(stdout);
}

void Elasticity2Da::readData_parameters()
{
	ifstream Frho_x;
	strcpy(path1,"rho_x.bin");
	strcpy(FullPath, path);
	strcat(FullPath,  path1);
	Frho_x.open(FullPath, ios::binary | ios::in);
	if (!Frho_x.is_open())
	{
		std::cout << "rho_x.bin failed." << std::endl;
	}
	ifstream Frho_z;
	strcpy(path1,  "rho_z.bin");
	strcpy(FullPath,  path);
	strcat(FullPath,  path1);
	Frho_z.open(FullPath, ios::binary | ios::in);
	if (!Frho_z.is_open())
	{
		std::cout << "rho_z.bin failed." << std::endl;
	}

	for (i = 0; i < (Nx - 1)*Nz; i++)
	{
		Frho_x.read((char*)&rho_x[i], sizeof hx);
	}
	for (i = 0; i < Nx*(Nz - 1); i++)
	{
		Frho_z.read((char*)&rho_z[i], sizeof hx);
	}
	Frho_x.close();
	Frho_z.close();

	ifstream Fc11;
	strcpy(path1, "c11.bin");
	strcpy(FullPath, path);
	strcat(FullPath, path1);
	Fc11.open(FullPath, ios::binary | ios::in);
	if (!Fc11.is_open())
	{
		std::cout << "c11.bin failed." << std::endl;
	}
	ifstream Fc13;
	strcpy(path1, "c13.bin");
	strcpy(FullPath,  path);
	strcat(FullPath,  path1);
	Fc13.open(FullPath, ios::binary | ios::in);
	if (!Fc13.is_open())
	{
		std::cout << "c13.bin failed." << std::endl;
	}
	ifstream Fc33;
	strcpy(path1, "c33.bin");
	strcpy(FullPath,  path);
	strcat(FullPath,  path1);
	Fc33.open(FullPath, ios::binary | ios::in);
	if (!Fc33.is_open())
	{
		std::cout << "c33.bin failed." << std::endl;
	}

	ifstream Ftau11;
	strcpy(path1, "tau11.bin");
	strcpy(FullPath, path);
	strcat(FullPath, path1);
	Ftau11.open(FullPath, ios::binary | ios::in);
	if (!Ftau11.is_open())
	{
		std::cout << "tau11.bin failed." << std::endl;
	}
	ifstream Ftau13;
	strcpy(path1, "tau13.bin");
	strcpy(FullPath, path);
	strcat(FullPath,  path1);
	Ftau13.open(FullPath, ios::binary | ios::in);
	if (!Ftau13.is_open())
	{
		std::cout << "tau13.bin failed." << std::endl;
	}
	ifstream Ftau33;
	strcpy(path1, "tau33.bin");
	strcpy(FullPath, path);
	strcat(FullPath, path1);
	Ftau33.open(FullPath, ios::binary | ios::in);
	if (!Ftau33.is_open())
	{
		std::cout << "tau33.bin failed." << std::endl;
	}

	for (i = 0; i < Nx*Nz; i++)
	{
		Fc11.read((char*)&c11[i], sizeof hx);
		Fc13.read((char*)&c13[i], sizeof hx);
		Fc33.read((char*)&c33[i], sizeof hx);
		Ftau11.read((char*)&tau11[i], sizeof hx);
		Ftau13.read((char*)&tau13[i], sizeof hx);
		Ftau33.read((char*)&tau33[i], sizeof hx);
	}

	Fc11.close();
	Fc13.close();
	Fc33.close();
	Ftau11.close();
	Ftau13.close();
	Ftau33.close();

	ifstream Ftau55;
	strcpy(path1, "tau55_xz.bin");
	strcpy(FullPath, path);
	strcat(FullPath, path1);
	Ftau55.open(FullPath, ios::binary | ios::in);
	if (!Ftau55.is_open())
	{
		std::cout << "tau55_xz.bin failed." << std::endl;
	}
	ifstream Fc55;
	strcpy(path1, "c55_xz.bin");
	strcpy(FullPath, path);
	strcat(FullPath, path1);
	Fc55.open(FullPath, ios::binary | ios::in);
	if (!Fc55.is_open())
	{
		std::cout << "c55_xz.bin failed." << std::endl;
	}

	for (i = 0; i < (Nx - 1)*(Nz - 1); i++)
	{
		Ftau55.read((char*)&tau55_xz[i], sizeof hx);
		Fc55.read((char*)&c55_xz[i], sizeof hx);
	}

	Ftau55.close();
	Fc55.close();

	ifstream Ftau_sigma;
	strcpy(path1, "tau_sigma.bin");
	strcpy(FullPath,  path);
	strcat(FullPath, path1);
	Ftau_sigma.open(FullPath, ios::binary | ios::in);
	if (!Ftau_sigma.is_open())
	{
		std::cout << "tau_sigma.bin failed." << std::endl;
	}
	Ftau_sigma.read((char*)&tau_sigma, sizeof(realval));
	Ftau_sigma.close();
}

void Elasticity2Da::filename(char*path, char* name, int i, char* F)
{
	char num[150];
	sprintf(num, "%d", i);
	int size_F = sizeof(F);
	int size_num = sizeof(i);
	strcpy(F, path);
	strcat(F,  name);
	strcat(F,  num);
	const char* s = ".bin";
	strcat(F, s);
}

realval Elasticity2Da::average_x(realval A[], int i, int j, int Nx)
{
	return 0.5*(A[j*Nx + i] + A[j*Nx + i + 1]);
}



realval Elasticity2Da::average_z(realval A[], int i, int j, int Nx)
{
	return 0.5*(A[j*Nx + i] + A[(j + 1)*Nx + i]);
}

realval Elasticity2Da::average_xz(realval A[], int i, int j, int Nx)
{
	realval av;
	if ((A[j*Nx + i] < 10) || (A[j*Nx + i + 1] < 10) || (A[(j + 1)*Nx + i] < 10) || (A[(j + 1)*Nx + i + 1] < 10))
		return 0;
	else
	{
		av = 0.25*(1 / A[j*Nx + i] + 1 / A[j*Nx + i + 1] + 1 / A[(j + 1)*Nx + i] + 1 / A[(j + 1)*Nx + i + 1]);
		return 1 / av;
	}
}

void Elasticity2Da::coordinate()
{
	realval *x_int = new realval[Nx];						// координаты "целых" узлов по x
	realval *z_int = new realval[Nz];						// координаты "целых" узлов по z
	realval *x_med = new realval[Nx - 1];					// координаты "дробных" узлов по x
	realval *z_med = new realval[Nz - 1];

	// Определение координат узлов
	for (i = 0; i < Nx; i++)
	{
		x_int[i] = i * hx;
	}

	for (i = 0; i < Nx - 1; i++)
	{
		x_med[i] = 0.5*hx + i * hx;
	}

	for (i = 0; i < Nz; i++)
	{
		z_int[i] = i * hz;
	}

	for (i = 0; i < Nz - 1; i++)
	{
		z_med[i] = 0.5*hz + i * hz;
	}

	ofstream x_int_res;
	x_int_res.open("x_int.bin", ios::binary | ios::out);

	ofstream z_int_res;
	z_int_res.open("z_int.bin", ios::binary | ios::out);

	ofstream x_med_res;
	x_med_res.open("x_med.bin", ios::binary | ios::out);

	ofstream z_med_res;
	z_med_res.open("z_med.bin", ios::binary | ios::out);

	for (i = 0; i < Nx; i++)
	{
		x_int_res.write((char*)&x_int[i], sizeof rx);
	}

	for (i = 0; i < Nz; i++)
	{
		z_int_res.write((char*)&z_int[i], sizeof rx);
	}

	for (i = 0; i < Nx - 1; i++)
	{
		x_med_res.write((char*)&x_med[i], sizeof rx);
	}

	for (i = 0; i < Nz - 1; i++)
	{
		z_med_res.write((char*)&z_med[i], sizeof rx);
	}

	x_int_res.close();
	z_int_res.close();
	x_med_res.close();
	z_med_res.close();

}

void Elasticity2Da::outputData()
{
	// запись данных в приемниках в файлы

	// создание снэпов
	if (inpData.NSnaps > 0)
	{
		DtPlot = inpData.Time / (realval)inpData.NSnaps;
	}

	if (cur_time == 0 || ((fabs(cur_time_str - counter * DtPlot) <= tau) && (counter <= inpData.NSnaps)))
	{
		ofstream v_x_res;
		ofstream v_z_res;
		ofstream q_x_res;
		ofstream q_z_res;
		ofstream sigma_xx_res;
		ofstream sigma_zz_res;
		ofstream sigma_xz_res;
		ofstream p_res;

		if (inpData.v_x_snap == 1)
		{
			strcpy(path1, "v_x_"); // не для MPI

			filename(path_data, path1, counter, FullPath);
			v_x_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			v_x_res.write((char*)v_x, (sizeof(rx)*(Nx - 1)*Nz));
			//v_x_res.write((char*)&v_x, (sizeof(rx)*(Nx - 1)*Nz));
			v_x_res.close();
		}

		if (inpData.v_z_snap == 1)
		{
			strcpy(path1,"v_z_"); // не для MPI

			filename(path_data, path1, counter, FullPath);
			v_z_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			v_z_res.write((char*)v_z, (sizeof(rx)*(Nz - 1)*Nx));
			//v_z_res.write((char*)&v_z, (sizeof(rx)*(Nz - 1)*Nx));
			v_z_res.close();
		}

		if (inpData.sigma_xx_snap == 1)
		{
			strcpy(path1,  "sigma_xx_"); // не для MPI

			filename(path_data, path1, counter, FullPath);
			sigma_xx_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			sigma_xx_res.write((char*)sigma_xx, (sizeof(rx)*Nx*Nz));
			//sigma_xx_res.write((char*)&sigma_xx, (sizeof(rx)*Nx*Nz));
			sigma_xx_res.close();
		}

		if (inpData.sigma_zz_snap == 1)
		{
			strcpy(path1, "sigma_zz_"); // не для MPI

			filename(path_data, path1, counter, FullPath);
			sigma_zz_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			sigma_zz_res.write((char*)sigma_zz, (sizeof(rx)*Nx*Nz));
			//sigma_zz_res.write((char*)&sigma_zz, (sizeof(rx)*Nx*Nz));
			sigma_zz_res.close();
		}

		if (inpData.sigma_xz_snap == 1)
		{
			strcpy(path1, "sigma_xz_"); // не для MPI

			filename(path_data, path1, counter, FullPath);
			sigma_xz_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			sigma_xz_res.write((char*)sigma_xz, (sizeof(rx)*(Nx - 1)*(Nz - 1)));
			//sigma_xz_res.write((char*)&sigma_xz, (sizeof(rx)*(Nx - 1)*(Nz - 1)));
			sigma_xz_res.close();
		}


		ofstream fileR;
		strcpy(path1,  "R_zz_"); // не для MPI

		std::ofstream outphi;
		outphi.open("R_zz_" + std::to_string(counter) + ".txt");
		if (outphi.is_open())
		{
			for (int it = 0; it < Nz - 1; it++)
			{
				for (int jt = 0; jt < Nx - 1; jt++)
				{
					outphi << Rzz[it*(Nx - 1) + jt] << " ";
				}
				outphi << std::endl;
			}
		}
		outphi.close();

		counter = counter + 1;
	}
}

void Elasticity2Da::Elasticity()
{
	//time_scale output 
	ofstream time_scale_vel;
	strcpy(path1,"time_scale_vel.bin"); // �� ��� MPI 
	strcpy(FullPath,  path);
	strcat(FullPath, path1);
	time_scale_vel.open(FullPath, ios::binary | ios::out);

	time_scale_vel.write((char*)&cur_time, sizeof rx);

	ofstream time_scale_str;
	strcpy(path1, "time_scale_str.bin"); // �� ��� MPI 
	strcpy(FullPath, path);
	strcat(FullPath, path1);
	time_scale_str.open(FullPath, ios::binary | ios::out);


	time_scale_str.write((char*)&cur_time_str, sizeof rx);

	//zero step data output 
	//outputData();

	//_rec ouput  

	//first receiver 
	ofstream v_x_rec_1;
	strcpy(path1, "v_x_rec_1.bin"); // �� ��� MPI 
	strcpy(FullPath,  path_data);
	strcat(FullPath,  path1);
	v_x_rec_1.open(FullPath, ios::binary | ios::out);

	ofstream v_z_rec_1;
	strcpy(path1,  "v_z_rec_1.bin"); // �� ��� MPI 
	strcpy(FullPath, path_data);
	strcat(FullPath, path1);
	v_z_rec_1.open(FullPath, ios::binary | ios::out);

	ofstream sigma_xx_rec_1;
	strcpy(path1, "sigma_xx_rec_1.bin"); // �� ��� MPI 
	strcpy(FullPath,  path_data);
	strcat(FullPath,  path1);
	sigma_xx_rec_1.open(FullPath, ios::binary | ios::out);

	ofstream sigma_zz_rec_1;
	strcpy(path1,  "sigma_zz_rec_1.bin"); // �� ��� MPI 
	strcpy(FullPath, path_data);
	strcat(FullPath,  path1);
	sigma_zz_rec_1.open(FullPath, ios::binary | ios::out);

	ofstream sigma_xz_rec_1;
	strcpy(path1, "sigma_xz_rec_1.bin"); // �� ��� MPI 
	strcpy(FullPath, path_data);
	strcat(FullPath,  path1);
	sigma_xz_rec_1.open(FullPath, ios::binary | ios::out);

	//second receiver 
	ofstream v_x_rec_2;
	strcpy(path1,"v_x_rec_2.bin"); // �� ��� MPI 
	strcpy(FullPath,path_data);
	strcat(FullPath, path1);
	v_x_rec_2.open(FullPath, ios::binary | ios::out);

	ofstream v_z_rec_2;
	strcpy(path1, "v_z_rec_2.bin"); // �� ��� MPI 
	strcpy(FullPath,  path_data);
	strcat(FullPath,path1);
	v_z_rec_2.open(FullPath, ios::binary | ios::out);

	ofstream sigma_xx_rec_2;
	strcpy(path1, "sigma_xx_rec_2.bin"); // �� ��� MPI 
	strcpy(FullPath,  path_data);
	strcat(FullPath,  path1);
	sigma_xx_rec_2.open(FullPath, ios::binary | ios::out);

	ofstream sigma_zz_rec_2;
	strcpy(path1,"sigma_zz_rec_2.bin"); // �� ��� MPI 
	strcpy(FullPath, path_data);
	strcat(FullPath, path1);
	sigma_zz_rec_2.open(FullPath, ios::binary | ios::out);

	ofstream sigma_xz_rec_2;
	strcpy(path1, "sigma_xz_rec_2.bin"); // �� ��� MPI 
	strcpy(FullPath, path_data);
	strcat(FullPath, path1);
	sigma_xz_rec_2.open(FullPath, ios::binary | ios::out);

	if (inpData.v_x_rec == 1)
	{
		for (i = 0; i < Nx - 1; i++)
		{
			v_x_rec_1.write((char*)&v_x[inpData.z_rec_1 * (Nx - 1) + i], sizeof rx);
			v_x_rec_2.write((char*)&v_x[inpData.z_rec_2 * (Nx - 1) + i], sizeof rx);
		}
	}

	if (inpData.v_z_rec == 1)
	{
		for (i = 0; i < Nx; i++)
		{
			v_z_rec_1.write((char*)&v_z[inpData.z_rec_1 * Nx + i], sizeof rx);
			v_z_rec_2.write((char*)&v_z[inpData.z_rec_2 * Nx + i], sizeof rx);
		}
	}

	if (inpData.sigma_xx_rec == 1)
	{
		for (i = 0; i < Nx; i++)
		{
			sigma_xx_rec_1.write((char*)&sigma_xx[inpData.z_rec_1 * Nx + i], sizeof rx);
			sigma_xx_rec_2.write((char*)&sigma_xx[inpData.z_rec_2 * Nx + i], sizeof rx);
		}
	}

	if (inpData.sigma_zz_rec == 1)
	{
		for (i = 0; i < Nx; i++)
		{
			sigma_zz_rec_1.write((char*)&sigma_zz[inpData.z_rec_1 * Nx + i], sizeof rx);
			sigma_zz_rec_2.write((char*)&sigma_zz[inpData.z_rec_2 * Nx + i], sizeof rx);
		}
	}

	if (inpData.sigma_xz_rec == 1)
	{
		for (i = 0; i < Nx - 1; i++)
		{
			sigma_xz_rec_1.write((char*)&sigma_xz[inpData.z_rec_1 * (Nx - 1) + i], sizeof rx);
			sigma_xz_rec_2.write((char*)&sigma_xz[inpData.z_rec_2 * (Nx - 1) + i], sizeof rx);
		}
	}


	realval v0 = inpData.v0;
	int i, j;
	realval f_src, c1, c2, c3, c4;
	realval tmp;

	//std::vector<realval>* swap;

	while (cur_time < inpData.Time)
	{
		/*swap = R_xx_new;
		R_xx_new = R_xx_old;
		R_xx_old = swap;

		swap = R_zz_new;
		R_zz_new = R_zz_old;
		R_zz_old = swap;

		swap = R_xz_new;
		R_xz_new = R_xz_old;
		R_xz_old = swap;*/

		cur_time_str = cur_time + 0.5*tau;
		cur_time = cur_time + tau;

		time_scale_vel.write((char*)&cur_time, sizeof rx);
		time_scale_str.write((char*)&cur_time_str, sizeof rx);

		// Without PML 
		//-----------memory variables------------//
/*#pragma omp parallel
		{
#pragma omp for private(i,j) schedule(guided)
			for (j = 1; j < Nz - 1; j++)
			{
				for (i = 1; i < Nx - 1; i++)
				{
					(*R_xx_new)[j*Nx + i] = (1.0 / (tau_sigma + tau * 0.5))*((tau_sigma - 0.5*tau)*(*R_xx_old)[j*Nx + i] - tau11[j*Nx + i] * rx*(v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + (i - 1)]) - tau13[j*Nx + i] * rz*(v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]));
					(*R_zz_new)[j*Nx + i] = (1.0 / (tau_sigma + tau * 0.5))*((tau_sigma - 0.5*tau)*(*R_zz_old)[j*Nx + i] - tau13[j*Nx + i] * rx*(v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + (i - 1)]) - tau33[j*Nx + i] * rz*(v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]));
				}
			}
		}
		//we don't care about boundary condition except periodic on top wall
#pragma omp parallel
		{
#pragma omp for private(j,f_src) schedule(guided)
			for (j = 1; j < Nz - 1; j++)
			{
				(*R_xx_new)[j*Nx] = (1.0 / (tau_sigma + tau * 0.5))*((tau_sigma - 0.5*tau)*(*R_xx_old)[j*Nx] - tau11[j*Nx] * rx*(v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx-2)]) - tau13[j*Nx] * rz*(v_z[j*Nx] - v_z[(j - 1)*Nx]));
				(*R_zz_new)[j*Nx] = (1 / (tau_sigma + tau * 0.5))*((tau_sigma - 0.5*tau)*(*R_zz_old)[j*Nx] - tau13[j*Nx] * rx*(v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx-2)]) - tau33[j*Nx] * rz*(v_z[j*Nx] - v_z[(j - 1)*Nx]));

				(*R_xx_new)[j*Nx + (Nx - 1)] = (*R_xx_new)[j*Nx];
				(*R_zz_new)[j*Nx + (Nx - 1)] = (*R_zz_new)[j*Nx];

			}
		}
#pragma omp parallel
		{
#pragma omp for private(i,j) schedule(guided)
			for (j = 0; j < Nz - 1; j++)
			{
				for (i = 0; i < Nx - 1; i++)
				{
					(*R_xz_new)[j*(Nx - 1) + i] = (1 / (tau_sigma + 0.5*tau))*((tau_sigma - 0.5*tau)*(*R_xz_old)[j*(Nx - 1) + i] - tau55_xz[j*(Nx - 1) + i] * (rz*(v_x[(j + 1)*(Nx - 1) + i] - v_x[j*(Nx - 1) + i]) + rx * (v_z[j*Nx + i + 1] - v_z[j*Nx + i])));
				}
			}
		}*/

		
		//--------sigma---------//

		//boundary condition for sigma_zz on left wall (j==0)
#pragma omp parallel 
		{
#pragma omp for private(i,j,f_src,tmp) schedule(guided) 
			for (i = 1; i < Nx - 1; i++)
			{
				j = 0;
				f_src = (1 - 2 * (PI*inpData.v0*(cur_time_str - 3.0 / v0))*(PI*v0*(cur_time_str - 3.0 / v0)))*exp(-((PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0))));
				tmp = rz * c33[j*Nx + i] * (1.0 / sqrt(rho_z[j*Nx + i] * c33[j*Nx + i]));
				sigma_zz[j*Nx + i] = (1.0 / (1 + tmp))*((1 - tmp)*sigma_zz[j*Nx + i] + c33[j*Nx + i] * rz * 2 * v_z[j*Nx + i] + rx * c13[j*Nx + i] * (v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + src[j*Nx + i] * f_src);

				//viscocity
				//sigma_zz[j*Nx + i] -= tau_sigma * ((*R_zz_new)[j*Nx + i] - (*R_zz_old)[j*Nx + i]);
			}
		}


#pragma omp parallel 
		{
#pragma omp for private(i,j,f_src) schedule(guided) 
			for (j = 1; j < Nz - 1; j++)
			{
				for (i = 1; i < Nx - 1; i++)
				{
					// sigma_xx, sigma_zz, inside 
					rrxx = rrzz = 0.0;
					for (int k = 0; k < K; k++)
					{
						rr = Rxx[j*Nx + i] * bbb - ccc * (tau11[j*Nx + i] * rx* (v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + tau13[j*Nx + i] * rz*(v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]));
						rrxx =rrxx + tau_sigma * (rr - Rxx[j*Nx + i]);
						/*if (abs(tau_sigma*(rr - Rxx[j*Nx + i])) > 1e-2)
						{
							std::cout << tau_sigma*(rr - Rxx[j*Nx + i]) << "   ";
						}*/
						Rxx[j*Nx + i] = rr;
						rr = Rzz[j*Nx + i] * bbb - ccc * (tau13[j*Nx + i] * rx*(v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + tau33[j*Nx + i] * rz* (v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]));
						rrzz = rrzz + tau_sigma * (rr - Rzz[j*Nx + i]);
						Rzz[j*Nx + i] = rr;

					}
					
					f_src = (1 - 2 * (PI*inpData.v0*(cur_time_str - 3.0 / v0))*(PI*v0*(cur_time_str - 3.0 / v0)))*exp(-((PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0))));
					sigma_xx[j*Nx + i] = sigma_xx[j*Nx + i] + rx * c11[j*Nx + i] * (v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + rz * c13[j*Nx + i] * (v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]) + src[j*Nx + i] * f_src;
					sigma_zz[j*Nx + i] = sigma_zz[j*Nx + i] + rx * c13[j*Nx + i] * (v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + rz * c33[j*Nx + i] * (v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]) + src[j*Nx + i] * f_src;
					//viscocity
					sigma_xx[j*Nx + i] = sigma_xx[j*Nx + i] - rrxx;
					sigma_zz[j*Nx + i] = sigma_zz[j*Nx + i] - rrzz;
					//sigma_zz[j*Nx + i] -= tau_sigma * ((*R_zz_new)[j*Nx + i] - (*R_zz_old)[j*Nx + i]);
					//sigma_xx[j*Nx + i] -= tau_sigma * ((*R_xx_new)[j*Nx + i] - (*R_xx_old)[j*Nx + i]);
					/*if (abs(rrxx)>1e-8)
					{
						std::cout << rrxx << "   ";
					}*/

				}
			}
		}

		//boundary condition for sigma_zz in right wall (j==Nz-1)
#pragma omp parallel 
		{
#pragma omp for private(i,j,f_src,tmp) schedule(guided) 
			for (i = 1; i < Nx - 1; i++)
			{
				j = Nz - 1;
				tmp = rz * c33[j*(Nx)+i] * (1.0 / sqrt(rho_z[(j - 1)*Nx + i] * c33[j*Nx + i]));
				f_src = (1 - 2 * (PI*inpData.v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0)))*exp(-((PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0))));
				sigma_zz[j*Nx + i] = (1.0 / (1.0 + tmp))*((1.0 - tmp)*sigma_zz[j*Nx + i] - c33[j*Nx + i] * rz * 2 * v_z[(j - 1)*Nx + i] + rx * c13[j*Nx + i] * (v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + src[j*Nx + i] * f_src);

				//viscocity
				//sigma_zz[j*Nx + i] -= tau_sigma * ((*R_zz_new)[j*Nx + i] - (*R_zz_old)[j*Nx + i]);

			}
		}


		//������������� ��������� ������� 
		//edge (j=0)
		j = 0;
		i = 0;
		f_src = (1 - 2 * (PI*inpData.v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0)))*exp(-((PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0))));
		tmp = rz * c33[j*Nx + i] * (1 / sqrt(rho_z[j*Nx + i] * c33[j*Nx + i]));
		sigma_zz[j*Nx + i] = (1 / (1 + tmp))*((1 - tmp)*sigma_zz[j*Nx + i] + c33[j*Nx + i] * rz * 2 * v_z[j*Nx + i] + rx * c13[j*Nx + i] * (v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + (Nx - 2)]) + src[j*Nx + i] * f_src);

		sigma_zz[j*Nx + (Nx - 1)] = sigma_zz[j*Nx];

		//viscocity
		//sigma_zz[j*Nx + i] -= tau_sigma * ((*R_zz_new)[j*Nx + i] - (*R_zz_old)[j*Nx + i]);
		//sigma_zz[j*Nx + (Nx-1)] -= tau_sigma * ((*R_zz_new)[j*Nx + (Nx-1)] - (*R_zz_old)[j*Nx + (Nx-1)]);


		//edge (j=Nz-1)
		j = Nz - 1;
		i = 0;
		f_src = (1 - 2 * (PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0)))*exp(-((PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0))));
		tmp = rz * c33[j*Nx + i] * (1 / sqrt(rho_z[(j - 1)*Nx + i] * c33[j*Nx + i]));
		sigma_zz[j*Nx] = (1 / (1 + tmp))*((1 - tmp)*sigma_zz[j*Nx + i] - c33[j*Nx + i] * rz * 2 * v_z[(j - 1)*Nx + i] + rx * c13[j*Nx] * (v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx - 2)]) + src[j*Nx] * f_src);

		sigma_zz[j*Nx + (Nx - 1)] = sigma_zz[j*Nx];
		//viscocity
		//sigma_zz[j*Nx + i] -= tau_sigma * ((*R_zz_new)[j*Nx + i] - (*R_zz_old)[j*Nx + i]);
		//sigma_zz[j*Nx + (Nx - 1)] -= tau_sigma * ((*R_zz_new)[j*Nx + (Nx - 1)] - (*R_zz_old)[j*Nx + (Nx - 1)]);


#pragma omp parallel 
		{
#pragma omp for private(j,f_src) schedule(guided) 
			for (j = 1; j < Nz - 1; j++)
			{
				// top and bottom boundary condition (periodic)
				rrxx = rrzz = 0.0;
				for (int k = 0; k < K; k++)
				{
					rr = Rxx[j*Nx + i] * bbb - ccc * (tau11[j*Nx + i] * rx*(v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + (Nx - 2)]) + tau13[j*Nx + i] * rz* (v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]));
					rrxx = rrxx+ tau_sigma * (rr - Rxx[j*Nx + i]);
					Rxx[j*Nx + i] = rr;
					rr = Rzz[j*Nx + i] * bbb - ccc * (tau13[j*Nx + i] * rx*(v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + (Nx - 2)]) + tau33[j*Nx + i] * rz*(v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]));
					rrzz = rrzz+tau_sigma * (rr - Rzz[j*Nx + i]);
					Rzz[j*Nx + i] = rr;
				}
				
				f_src = (1 - 2 * (PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0)))*exp(-((PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0))));
				sigma_xx[j*Nx] = sigma_xx[j*Nx] + rx * c11[j*Nx] * (v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx - 2)]) + rz * c13[j*Nx] * (v_z[j*Nx] - v_z[(j - 1)*Nx]) + src[j*Nx] * f_src;
				sigma_zz[j*Nx] = sigma_zz[j*Nx] + rx * c13[j*Nx] * (v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx - 2)]) + rz * c33[j*Nx] * (v_z[j*Nx] - v_z[(j - 1)*Nx]) + src[j*Nx] * f_src;

				//viscocity
				sigma_zz[j*Nx] = sigma_zz[j*Nx] - rrzz;
				sigma_xx[j*Nx] = sigma_xx[j*Nx] - rrxx;

				sigma_xx[j*Nx + (Nx - 1)] = sigma_xx[j*Nx];
				sigma_zz[j*Nx + (Nx - 1)] = sigma_zz[j*Nx];


				//sigma_zz[j*Nx] -= tau_sigma * ((*R_zz_new)[j*Nx] - (*R_zz_old)[j*Nx]);
				//sigma_xx[j*Nx] -= tau_sigma * ((*R_xx_new)[j*Nx] - (*R_xx_old)[j*Nx]);
				//sigma_zz[j*Nx+(Nx-1)] -= tau_sigma * ((*R_zz_new)[j*Nx] - (*R_zz_old)[j*Nx]);
				//sigma_xx[j*Nx+(Nx-1)] -= tau_sigma * ((*R_xx_new)[j*Nx] - (*R_xx_old)[j*Nx]);

			}
		}
#pragma omp parallel 
		{
#pragma omp for private(i,j) schedule(guided) 
			for (j = 0; j < Nz - 1; j++)
			{
				for (i = 0; i < Nx - 1; i++)
				{
					rrxz = 0.0;
					for (int k = 0; k < K; k++)
					{
						rr = Rxz[j*(Nx - 1) + i] * bbb - ccc * tau55_xz[j*(Nx - 1) + i] * (rx*(v_z[j*Nx + i + 1] - v_z[j*Nx + i]) + rz * (v_x[(j + 1)*(Nx - 1) + i] - v_x[j*(Nx - 1) + i]));
						rrxz =rrxz+ tau_sigma * (rr - Rxz[j*(Nx - 1) + i]);
						Rxz[j*(Nx - 1) + i] = rr;
					}
					
					sigma_xz[j*(Nx - 1) + i] = sigma_xz[j*(Nx - 1) + i] + c55_xz[j*(Nx - 1) + i] * (rx*(v_z[j*Nx + i + 1] - v_z[j*Nx + i]) + rz * (v_x[(j + 1)*(Nx - 1) + i] - v_x[j*(Nx - 1) + i]));// +src[j*Nx] * f_src;

					//viscocity
					sigma_xz[j*(Nx - 1) + i] = sigma_xz[j*(Nx - 1) + i] - rrxz;
					//sigma_xz[j*(Nx - 1) + i] -= tau_sigma * ((*R_xz_new)[j*(Nx - 1) + i] - (*R_xz_old)[j*(Nx - 1) + i]);
				}
			}
		}
		
		// ---------------u_x, u_z----------------//

		//bc for u_x on left wall (j==0)
#pragma omp parallel 
		{
#pragma omp for private(i,j,tmp) schedule(guided) 
			for (i = 0; i < Nx - 1; i++)
			{
				j = 0;
				tmp = rz * (1 / rho_x[j*(Nx - 1) + i])*(sqrt(rho_x[j*(Nx - 1) + i] * c55_xz[j*(Nx - 1) + i]));
				v_x[j*(Nx - 1) + i] = (1.0 / (1.0 + tmp))*((1.0 - tmp)*v_x[j*(Nx - 1) + i] + (1.0 / rho_x[j*(Nx - 1) + i])*rz * 2 * sigma_xz[j*(Nx - 1) + i] + (1 / rho_x[j*(Nx - 1) + i])*rx*(sigma_xx[j*Nx + i + 1] - sigma_xx[j*Nx + i]));
			}
		}



#pragma omp parallel 
		{
#pragma omp for private(i,j,c1,c2,c3,c4) schedule(guided) 
			for (j = 1; j < Nz - 1; j++)
			{
				for (i = 0; i < Nx - 1; i++)
				{
					v_x[j*(Nx - 1) + i] = v_x[j*(Nx - 1) + i] + (1.0 / rho_x[j*(Nx - 1) + i]) * (rx*(sigma_xx[j*Nx + i + 1] - sigma_xx[j*Nx + i]) + rz * (sigma_xz[j*(Nx - 1) + i] - sigma_xz[(j - 1)*(Nx - 1) + i]));
				}
			}
		}

		//bc on right wall
#pragma omp parallel 
		{
#pragma omp for private(i,j,c1,c2,c3,c4) schedule(guided) 
			for (i = 0; i < Nx - 1; i++)
			{
				j = Nz - 1;
				tmp = (1 / rho_x[j*(Nx - 1) + i])*rz*(sqrt(rho_x[j*(Nx - 1) + i] * c55_xz[(j - 1)*(Nx - 1) + i]));
				v_x[j*(Nx - 1) + i] = (1 / (1 + tmp))*((1 - tmp)*v_x[j*(Nx - 1) + i] - (1 / rho_x[j*(Nx - 1) + i])*rz * 2 * sigma_xz[(j - 1)*(Nx - 1) + i] + (1 / rho_x[j*(Nx - 1) + i]) * (rx*(sigma_xx[j*Nx + i + 1] - sigma_xx[j*Nx + i])));
			}
		}

		// ���������� �������� �������� � ������� ������ �� n+1 (������������ ����������) 
#pragma omp parallel 
		{
#pragma omp for private(i,j,c1,c2,c3,c4) schedule(guided) 
			for (j = 0; j < Nz - 1; j++)
			{
				for (i = 1; i < Nx - 1; i++)
				{
					v_z[j*Nx + i] = v_z[j*Nx + i] + (1.0 / rho_z[j*Nx + i]) * (rx*(sigma_xz[j*(Nx - 1) + i] - sigma_xz[j*(Nx - 1) + i - 1]) + rz * (sigma_zz[(j + 1)*Nx + i] - sigma_zz[j*Nx + i]));
				}
			}
		}
		//������������� ��������� ������� 
#pragma omp parallel 
		{
#pragma omp for private(j,c1,c2,c3,c4) schedule(guided) 
			for (j = 0; j < Nz - 1; j++)
			{
				v_z[j*Nx] = v_z[j*Nx] + (1.0 / rho_z[j*Nx])* (rx*(sigma_xz[j*(Nx - 1)] - sigma_xz[j*(Nx - 1) + (Nx - 2)]) + rz * (sigma_zz[(j + 1)*Nx] - sigma_zz[j*Nx]));
				v_z[j*Nx + (Nx - 1)] = v_z[j*Nx];

			}
		}


		outputData();

		if (inpData.v_x_rec == 1)
		{
			for (i = 0; i < Nx - 1; i++)
			{
				v_x_rec_1.write((char*)&v_x[inpData.z_rec_1 * (Nx - 1) + i], sizeof rx);
				v_x_rec_2.write((char*)&v_x[inpData.z_rec_2 * (Nx - 1) + i], sizeof rx);
			}
		}

		if (inpData.v_z_rec == 1)
		{
			for (i = 0; i < Nx; i++)
			{
				v_z_rec_1.write((char*)&v_z[inpData.z_rec_1 * Nx + i], sizeof rx);
				v_z_rec_2.write((char*)&v_z[inpData.z_rec_2 * Nx + i], sizeof rx);
			}
		}

		if (inpData.sigma_xx_rec == 1)
		{
			for (i = 0; i < Nx; i++)
			{
				sigma_xx_rec_1.write((char*)&sigma_xx[inpData.z_rec_1 * Nx + i], sizeof rx);
				sigma_xx_rec_2.write((char*)&sigma_xx[inpData.z_rec_2 * Nx + i], sizeof rx);
			}
		}

		if (inpData.sigma_zz_rec == 1)
		{
			for (i = 0; i < Nx; i++)
			{
				sigma_zz_rec_1.write((char*)&sigma_zz[inpData.z_rec_1 * Nx + i], sizeof rx);
				sigma_zz_rec_2.write((char*)&sigma_zz[inpData.z_rec_2 * Nx + i], sizeof rx);
			}
		}

		if (inpData.sigma_xz_rec == 1)
		{
			for (i = 0; i < Nx - 1; i++)
			{
				sigma_xz_rec_1.write((char*)&sigma_xz[inpData.z_rec_1 * (Nx - 1) + i], sizeof rx);
				sigma_xz_rec_2.write((char*)&sigma_xz[inpData.z_rec_2 * (Nx - 1) + i], sizeof rx);
			}
		}

	}

	time_scale_vel.close();
	time_scale_str.close();

}
