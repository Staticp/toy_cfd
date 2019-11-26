/*************************************************************************
    > File Name: Mesh.cpp
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年03月23日 星期六 11时37分32秒
 ************************************************************************/

#include "Mesh.h"
#include <cmath>
#include <fstream>
#include <iostream>

Mesh::Mesh():
	X_length(1),Y_length(1),M(10),N(20)
{
	deltaX = X_length/(M-1);
	deltaEta = Eta()/(N-1);
	calculate_coords();
}

Mesh::Mesh(const int m, const int n):
	X_length(1),Y_length(1),M(m),N(n)
{
	deltaX = X_length/(M-1);
	deltaEta = Eta()/(N-1);
	calculate_coords();
}

Mesh::Mesh(const double x, const double y, const int m, const int n):
	X_length(x),Y_length(y),M(m),N(n)
{
	deltaX = X_length/(M-1);
	cout<<"Eta= "<<Eta()<<endl;
	deltaEta = Eta()/(N-1);
	calculate_coords();
}

double Mesh::Eta()
{
	return log((pow(e,s)-pow(e,-s)+sqrt(pow(pow(e,s)-pow(e,-s),2)+4))/2.0)/s;
}

double Mesh::f_eta(const double eta)
{
	return Y_length*sinh(s*eta)/sinh(s);
}

double Mesh::p_f_eta(const double eta)
{
	return s*Y_length*cosh(s*eta)/sinh(s);
}

double Mesh::pp_f_eta(const double eta)
{
	return s*s*Y_length*sinh(s*eta)/sinh(s);
}

double Mesh::Y_(const double eta)
{
	return Y_length*sinh(s*eta)/sinh(s);
}

double Mesh::eta_y(const double eta)
{
	return 1.0/p_f_eta(eta);
}

double Mesh::eta_yy(const double eta)
{
	return -pp_f_eta(eta)/pow(p_f_eta(eta), 3);
}

void Mesh::calculate_coords()
{
	vector<PAIR> tmp_points_com,tmp_points_phy;
	PAIR point_com, point_phy;

	for(int i=0;i!=M;++i){
		for(int j=0;j!=N;++j){
			point_com.first = i*deltaX;
			point_com.second = j*deltaEta;
			point_phy.first = point_com.first;
			point_phy.second = f_eta(point_com.second);

			tmp_points_com.push_back(point_com);	//先存列，这样沿X方向推进快
			tmp_points_phy.push_back(point_phy);
		}
		Points_com.push_back(tmp_points_com);
		Points_phy.push_back(tmp_points_phy);
		tmp_points_com.clear();
		tmp_points_phy.clear();
	}
}

void Mesh::phy_fout()
{
	ofstream out;
	out.open("phy_domain.dat");
	out << "variables=x,y" << endl;

	for(int j=0; j!=N; ++j){
		for(int i=0; i!=M; ++i){
			out << Points_phy[i][j].first << " " 
				<<Points_phy[i][j].second <<endl;
		}
	}
	out.close();
}

void Mesh::com_fout()
{
	ofstream out;
	out.open("com_domain.dat");
	out << "variables=x,y" << endl;

	for(int j=0; j!=N; ++j){
		for(int i=0; i!=M; ++i){
			out << Points_com[i][j].first << " " 
				<<Points_com[i][j].second <<endl;
		}
	}
	out.close();
}

