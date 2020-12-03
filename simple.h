/*************************************************************************
    > File Name: simple.cpp
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年01月15日 星期二 21时59分27秒
 ************************************************************************/

#include<iostream>
#include<string>
#include<vector>
#define Np 21
#define Mp 11
#define Nu 22
#define Mu 11
#define Nv 23
#define Mv 12
#define rho 1.225
#define mu 1.79e-5
#define deltaX 0.00762
#define deltaY 0.0003048
#define deltaT 0.001
#define P_CORR_COUNT 1000
#define PD_TOL 1e-8
#define alpha_p 0.1
#define ITER_MAX 40

using namespace std;

class simple{

public:
	simple();
	simple(double, double, double, double);

	void loop();
	void foutput();

private:

	vector<double> p;
	vector<double> pd;
	vector<double> u;
	vector<double> v;
	vector<double> rhoU;
	vector<double> rhoV;

	vector<vector<double> > p_;
	vector<vector<double> > pd_;
	vector<vector<double> > u_;
	vector<vector<double> > v_;
	vector<vector<double> > rhoU_;
	vector<vector<double> > rhoV_;

	vector<vector<double> > p_N;
	vector<vector<double> > pd_N;
	vector<vector<double> > u_N;
	vector<vector<double> > v_N;
	vector<vector<double> > rhoU_N;
	vector<vector<double> > rhoV_N;

	void initialize(double, double, double, double);
	void momentum_predict();
	void pd_relax_equation();
	void p_correction();
};
