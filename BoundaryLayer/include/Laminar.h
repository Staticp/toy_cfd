/*************************************************************************
    > File Name: Laminar.h
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年03月23日 星期六 22时40分43秒
 ************************************************************************/

#include<iostream>
#include"Mesh.h"
using namespace std;

class Laminar:public Mesh
{
public:
	Laminar(const double, const double, const int, const int);
	virtual void forward();
	virtual	void file_output();

protected:

	vector<double> a_coe;
	vector<double> b_coe;
	vector<double> c_coe;

	vector<vector<double> >A;	//系数阵
	vector<double>b;	//求解矩阵的源

	vector<double> U;	//功能性向量，不参与初始化，下同
	vector<double> V;

	vector<vector<double> >UU;
	vector<vector<double> >VV;

	void initialize();
	void lower_boundary();
	void upper_boundary();
	void assemble_momentum();
	void assemble_continue(const vector<double>&);
	virtual vector<double> calcu_a_coe();
	virtual vector<double> calcu_b_coe(const vector<double>&);
	vector<double> calcu_c_coe();
	vector<double> LU(vector<vector<double> >&, vector<double>&) const;//普通的LU分解
	vector<double> ZG(vector<vector<double> >&, vector<double>&) const;//追赶法
};
