#pragma once
#include <cmath>
#include <iostream> 
#include <fstream> 
#include <stdio.h>
#include <assert.h>
#include <omp.h>
#include <string.h>
#include <sstream>
#include <time.h>
#include<vector>
#include <chrono>

using namespace std;

typedef double realval;
typedef long myint; //for WINDOWS

class Elasticity2D
{
	struct inputData
	{
		realval v0; //frequency
		realval Time; // time of propagation
		int srcplace; //source placement: 0-fluid; 1-solid
		int z_src; // source z-node
		int z_rec_1; // receiver 1 z-node
		int z_rec_2; // receiver 2 z-node
		int Nz_PML_l; // left PML length
		int Nz_PML_r; // right PML length


		int sigma_xx_rec;  // normal stress in x-direction record lines data: 1-record; 0-do not record
		int sigma_zz_rec;  // normal stress in z-direction record lines data: 1-record; 0-do not record
		int sigma_xz_rec;  // shear stress record lines data: 1-record; 0-do not record
		int v_x_rec;       // x-component of solid velocity record lines data: 1-record; 0-do not record
		int v_z_rec;       // z-component of solid velocity record lines data: 1-record; 0-do not record
		int sigma_xx_snap; // normal stress in x-direction snapshots: 1-make; 0-do not make
		int sigma_zz_snap; // normal stress in z-direction snapshots: 1-make; 0-do not make
		int sigma_xz_snap; // shear stress snapshots: 1-make; 0-do not make
		int v_x_snap;      // x-component of solid velocity snapshots: 1-make; 0-do not make
		int v_z_snap;      // z-component of solid velocity snapshots: 1-make; 0-do not make
		int NSnaps;        // number of snapshots during propagation time (through equal time intervals)

	};

	void filename(char* path, char* name, int i, char* F);
	realval average_x(realval A[], int i, int j, int Nx);
	realval average_z(realval A[], int i, int j, int Nx);
	realval average_xz(realval A[], int i, int j, int Nx);

	void readData_geometry();
	void readData_parameters();
	//for outout files
	void outputData(); //sigma, v and time_scale
	void coordinate();

public:
	Elasticity2D();
	~Elasticity2D();
	void Elasticity();
private:
	inputData inpData;
	int Nx;
	int Nz;
	realval hx;
	realval hz;
	realval tau;
	//char path1[150] = "";
	//char path_data[150] = "10/1/30000/Data/";
	//char path[150] = "10/1/30000/";
	std::string path_data;
	std::string path;
	std::string FullPath;
	std::string path1;
	//char FullPath[150];

private:
	//from input files
	realval* rho_x; //average density by x
	realval* rho_z; //by z
	realval* c11;
	realval* c13;
	realval* c33;
	realval* c55_xz;

	realval *sigma_xx;					// нормальное напряжение упругого материала по x
	realval *sigma_zz;					// нормальное напряжение упругого материала по z
	realval *sigma_xz;		// касательное напряжение упругого материала
	realval *v_x;					// скорость твердх частиц по x
	realval *v_z;

	//for viscoelasticity 
	/*std::vector<realval>*R_xx_old;
	std::vector<realval>*R_zz_old;
	std::vector<realval>*R_xz_old;

	std::vector<realval>*R_xx_new;
	std::vector<realval>*R_zz_new;
	std::vector<realval>*R_xz_new;*/

	realval rr, rrxx, rrzz, rrxz;
	realval aa, bbb, ccc;
	realval*Rxx;
	realval*Rzz;
	realval*Rxz;

	realval tau_sigma;
	realval *tau11;
	realval *tau13;
	realval *tau33;
	realval *tau55_xz;

	int K = 1;

	//coordinate of source
	realval *src;

	//std::vector<realval> src;

	int i, j;
	int counter = 1;
	realval cur_time, cur_time_str;			// время на текущем шаге для скоростей и напряжений соответственно
	realval rx, rz;							// tau/hx, tau/hz
	realval DtPlot;                       // вспомогательные параметры и функция источника и его частота
	const realval PI = 3.14159265;
	int flag = 0;

};

