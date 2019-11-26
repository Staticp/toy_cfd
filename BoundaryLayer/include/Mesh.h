/*************************************************************************
    > File Name: Mesh.h
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年03月23日 星期六 10时30分09秒
 ************************************************************************/

#include<iostream>
#include<vector>
#include<utility>

#define PAIR pair<double,double>
#define e 2.71828

using namespace std;

class Mesh
{
public:
	Mesh();
	Mesh(const int, const int);
	Mesh(const double, const double, const int, const int);
	vector<vector<PAIR> > Domain_com() const
	{
		return Points_com;
	}

	vector<vector<PAIR> > Domain_phy() const
	{
		return Points_phy;
	}

	void phy_fout();
	void com_fout();

protected:

	double X_length;
	double Y_length;
	int M;
	int N;
	double deltaX;
	double deltaEta;
	vector<vector<PAIR> > Points_com;
	vector<vector<PAIR> > Points_phy;
	double Y_(const double);
	double eta_y(const double);
	double eta_yy(const double);

private:
	double Eta();
	double f_eta(const double);
	double p_f_eta(const double);
	double pp_f_eta(const double);
	void calculate_coords();
	double s=4.0;
};
