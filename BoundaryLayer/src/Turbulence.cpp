/*************************************************************************
    > File Name: Laminar.cpp
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年03月23日 星期六 22时49分41秒
 ************************************************************************/

#include <iostream>
#include <cmath>
#include <fstream>
#include "Turbulence.h"

#define kp 0.4
#define Re 160000
#define V_TOL 1e-6
#define U_TOL 1e-6

using namespace std;

Turbulence::Turbulence(const double x,const double y,const int m,const int n):
	Laminar(x,y,m,n)
{
	initialize();	//初始化了第一层的UU和VV,这时UU和VV都只有一层 1*N
	c_coe = calcu_c_coe();
}

vector<double> Turbulence::calcu_U_y(const vector<double>& U_)
{
	vector<double> tmp;
	for(int j=0;j!=N;++j){
		if(j==0)
			tmp.push_back((U_[j+1]-U_[j])/deltaEta*eta_y(j*deltaEta));
		else if(j==N-1)
			tmp.push_back((U_[j]-U_[j-1])/deltaEta*eta_y(j*deltaEta));
		else
			tmp.push_back(0.5*(U_[j+1]-U_[j-1])/deltaEta*eta_y(j*deltaEta));
	}
	return tmp;
}

vector<double> Turbulence::calcu_nut(const vector<double>& U_)	//nut要用什么边界条件呢？
{
	vector<double> tmp;
	vector<double> U_y = calcu_U_y(U_);
	for(int j=0;j!=N;++j){
		if(j==0)
			tmp.push_back(0);
		else
			tmp.push_back(pow(kp*Y_(j*deltaEta),2)*fabs(U_y[j])*Re/X_length);
	}
	return tmp;
}

vector<double> Turbulence::calcu_nut_y(const vector<double>& nut)
{
	vector<double> tmp;
	for(int j=0;j!=N;++j){
		if(j==0){
			tmp.push_back((nut[j+1]-nut[j])/deltaEta*eta_y(j*deltaEta));
			//tmp.push_back(kp*kp*(2*Y_(j*deltaEta)*fabs(U_y[j])
			//			+Y_(j*deltaEta)*Y_(j*deltaEta)*eta_y(j*deltaEta)*(fabs(U_y[j+1])-fabs(U_y[j]))/deltaEta));
		}
		else if(j==N-1){
			tmp.push_back((nut[j]-nut[j-1])/deltaEta*eta_y(j*deltaEta));
			//tmp.push_back(kp*kp*(2*Y_(j*deltaEta)*fabs(U_y[j])
			//			+Y_(j*deltaEta)*Y_(j*deltaEta)*eta_y(j*deltaEta)*(fabs(U_y[j])-fabs(U_y[j-1]))/deltaEta));
		}
		else{
			tmp.push_back(0.5*(nut[j+1]-nut[j-1])/deltaEta*eta_y(j*deltaEta));
//			tmp.push_back(kp*kp*(2*Y_(j*deltaEta)*fabs(U_y[j])
//						+Y_(j*deltaEta)*Y_(j*deltaEta)*eta_y(j*deltaEta)*0.5*(fabs(U_y[j+1])-fabs(U_y[j-1]))/deltaEta));
//			tmp.push_back(kp*kp*(2*Y_(j*deltaEta)*fabs(0.5*(U_[j+1]-U_[j-1])/deltaEta*eta_y(j*deltaEta))
//						+Y_(j*deltaEta)*Y_(j*deltaEta)*fabs(eta_yy(j*deltaEta)*0.5*(U_[j+1]-U_[j-1])/deltaEta
//							+eta_y(j*deltaEta)*eta_y(j*deltaEta)*(U_[j+1]-2*U_[j]+U_[j-1])/(deltaEta*deltaEta)))
//					);
		}
	}
	return tmp;
}

vector<double> Turbulence::calcu_a_coe()	//在计算湍流边界层时，系数A也是随着迭代进行不断变化的
{
	vector<double> tmp;
	for(int j=0;j!=N;++j){
		tmp.push_back((1+nut[j])*pow(eta_y(j*deltaEta), 2)/Re);
	}
	return tmp;
}

vector<double> Turbulence::calcu_b_coe(const vector<double>& V)	//B要根据V来算，所以每一步都要更新
{
	vector<double> tmp;
	for(int j=0;j!=N;++j){
		tmp.push_back((1+nut[j])*eta_yy(j*deltaEta)/Re 
				- V[j] * eta_y(j*deltaEta)
				+ nut_y[j] * eta_y(j*deltaEta)/Re);
	}
	return tmp;
}

void Turbulence::forward()
{
	vector<double> U_tmp, V_tmp;
	double aver_nut = 0.0, aver_b = 0.0;

	for(int i=0; i!=M-1; ++i){	//从第一层开始推进,到倒数第二层
		bool mark = 1;
		int iter=0;

		vector<double> U_N(UU[i].size(), 0.0), V_N(VV[i].size(), 0.0);

		while(mark && iter < 100){
			mark = 0;	// mark用来判断否是收敛
			cout << "目前是第"<<i+1<<"层的第"<<iter+1<<"次迭代！！！！！"<<endl;

			for(int j=0; j!=N; ++j){
				aver_nut = 0.5*(calcu_nut(U)[j] + calcu_nut(U_N)[j]);
				nut.push_back(aver_nut);
			}
			nut_y=calcu_nut_y(nut);

			a_coe = calcu_a_coe();

			for(int j=0; j!=N; ++j){
				aver_b = 0.5*(calcu_b_coe(V)[j] + calcu_b_coe(V_N)[j]);	//用这一步和下一步
				b_coe.push_back(aver_b);
			}
			assemble_momentum();	//装配矩阵,使用了上一步的a_coe和b_coe

			nut.clear();
			nut_y.clear();
			a_coe.clear();
			b_coe.clear();
			U_tmp=ZG(A,b);	//根据动量方程求解出下一层的U_tmp

			A.clear();
			b.clear();	//clear()的原因是之后求解V仍要重复使用
			assemble_continue(U_tmp);
			V_tmp=ZG(A,b);		//根据连续方程求解出下一层的V_N
			A.clear();
			b.clear();

			for(int j=0;j!=N;++j){
				if(fabs(U_tmp[j]-U_N[j]) > U_TOL || fabs(V_tmp[j]-V_N[j]) > V_TOL)//在一层迭代多次，每一个量都收敛就推进一步
					mark=1;
			}
			U_N.swap(U_tmp);
			V_N.swap(V_tmp);
			++iter;
		}
		U.swap(U_N);	//一步收敛后推进到下一步，U是新的值
		V.swap(V_N);
		UU.push_back(U);
		VV.push_back(V);
	}
	U_tmp.clear();
	V_tmp.clear();
	U.clear();
	V.clear();
}

void Turbulence::file_output()
{
	ofstream out;
	out.open("U_tur.dat");
	out<<"variables = x, y, U"<<endl;
	out<<"ZONE I=" << M <<", J=" << N << ",DATAPACKING=POINT"<<endl;
	for(int j=0; j!=N; ++j){
		for(int i=0; i!=M; ++i){
			out << Points_phy[i][j].first << " " 
				<<Points_phy[i][j].second << " " << UU[i][j] << endl;
		}
	}
	out.close();

	out.open("V_tur.dat");
	out<<"variables = x, y, V"<<endl;
	out<<"ZONE I=" << M <<", J=" << N << ",DATAPACKING=POINT"<<endl;
	for(int j=0; j!=N; ++j){
		for(int i=0; i!=M; ++i){
			out << Points_phy[i][j].first << " " 
				<<Points_phy[i][j].second << " " << VV[i][j] << endl;
		}
	}
	out.close();

	out.open("eta_U_025_tur.dat");
	out<< "variables = U/U∞, η" << endl;
    int k = M/4;
	for(int j=0;j!=N;++j){
		out << UU[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_U_050_tur.dat");
	out<< "variables = U/U∞, η" << endl;
	k = M/2;
	for(int j=0;j!=N;++j){
		out << UU[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_U_075_tur.dat");
	out<< "variables = U/U∞, η" << endl;
	k = 3*M/4;
	for(int j=0;j!=N;++j){
		out << UU[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_V_025_tur.dat");
	out<< "variables = V/U∞, η" << endl;
	k = M/4;
	for(int j=0;j!=N;++j){
		out << VV[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_V_050_tur.dat");
	out<< "variables = V/U∞, η" << endl;
	k = M/2;
	for(int j=0;j!=N;++j){
		out << VV[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_V_075_tur.dat");
	out<< "variables = V/U∞, η" << endl;
	k = 3*M/4;
	for(int j=0;j!=N;++j){
		out << VV[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("layer_thickness_approximate_tur.dat");
	out << "variables = X, Thickness" << endl;
	for(double x=0.01; x<X_length; x += X_length/M)
		out << x << " " << x*0.37/pow(Re*x/X_length, 0.2) << endl;
	out.close();

	out.open("layer_thickness_compute_tur.dat");
	out << "variables = X, Thickness" << endl;
	for(int i=0;i!=M;++i){
		for(int j=0;j!=N;++j){
			if(UU[i][j] >= 0.99){
				out << Points_phy[i][j].first << " " << Points_phy[i][j].second << endl;
				break;
			}
		}
	}
	out.close();

	out.open("Cf_approximate_tur.dat");
	out << "variables = X, Thickness" << endl;
	for(double x=0.01; x<X_length; x += X_length/M)
		out << x << " " << 0.074/pow(Re*x/X_length, 0.2) << endl;
	out.close();

	out.open("Cf_compute_tur.dat");
	out << "variables = X, Cf" << endl;
	for(int i=0; i!=M; ++i){
		out << Points_phy[i][0].first << " " << (UU[i][1] - UU[i][0])/deltaEta * eta_y(deltaEta) * X_length/(0.5*Re) << endl;
	}
	out.close();

	out.open("yPlus.dat");
	out << "variables = x, y+" << endl;
	for(int i=0;i!=M;++i)
		out << Points_phy[i][0].first << " " << UU[i][1]*Points_phy[i][1].second*Re/X_length << endl;
	out.close();

	out.open("uPlus1.dat");
	out << "variables = ln(y+), u+" << endl;
	for(int i=0;i!=M;++i)
		out << UU[i][1]*Points_phy[i][1].second*Re/X_length << " " <<  UU[i][1]*Points_phy[i][1].second*Re/X_length << endl;
	out.close();

	out.open("uPlus2.dat");
	out << "variables = ln(y+), u+" << endl;
	for(int i=0;i!=M;++i)
		out << UU[i][1]*Points_phy[i][1].second*Re/X_length << " " <<  2.5*log(UU[i][1]*Points_phy[i][1].second*Re/X_length)+5 << endl;
	out.close();
}
