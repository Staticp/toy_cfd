/*************************************************************************
    > File Name: Laminar.h
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年03月23日 星期六 22时40分43秒
 ************************************************************************/

#include"Laminar.h"

using namespace std;

class Turbulence:public Laminar
{
public:
	Turbulence(const double, const double, const int, const int);
	void forward() final; 
	void file_output() final;

protected:

	vector<double> calcu_a_coe() final;
	vector<double> calcu_b_coe(const vector<double>&) final;

private:
	vector<double> calcu_U_y(const vector<double>&);
	vector<double> calcu_nut(const vector<double>&);
	vector<double> calcu_nut_y(const vector<double>&);

	vector<double> nut;	//功能性向量,实际上使用的是1/2节点处的nut
	vector<double> nut_y;	//指nut关于y的偏导数
};
