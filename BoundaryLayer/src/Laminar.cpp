/*************************************************************************
    > File Name: Laminar.cpp
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年03月23日 星期六 22时49分41秒
 ************************************************************************/

#include <iostream>
#include <cmath>
#include <fstream>
#include "Laminar.h"

#define Re 160000
#define U_TOL 1E-6
#define V_TOL 1E-6

using namespace std;

Laminar::Laminar(const double x,const double y,const int m,const int n):
	Mesh(x,y,m,n)
{
	initialize();	//初始化了第一层的UU和VV,这时UU和VV都只有一层 1*N
	a_coe = calcu_a_coe();
	c_coe = calcu_c_coe();
}

void Laminar::initialize()
{
	lower_boundary();

	for(int j=1;j!=N-1;++j){
		U.push_back(1.0);
		V.push_back(0.0);
	}

	upper_boundary();
	UU.push_back(U);
	VV.push_back(V);	//先把第一层初始化的加进去，之后每迭代一步就加一层,边求解边添加元素
}

void Laminar::lower_boundary()
{
	U.push_back(0.0);
	V.push_back(0.0);
}

void Laminar::upper_boundary()
{
	U.push_back(1.0);
	V.push_back(0.0);
}

vector<double> Laminar::calcu_a_coe()
{
	vector<double> tmp;
	for(int j=0;j!=N;++j){
		tmp.push_back(pow(eta_y(j*deltaEta), 2)/Re);
	}
	return tmp;
}

vector<double> Laminar::calcu_b_coe(const vector<double>& V)	//B要根据V来算，所以每一步都要更新
{
	vector<double> tmp;
	for(int j=0;j!=N;++j){
		tmp.push_back(eta_yy(j*deltaEta)/Re- V[j] *eta_y(j*deltaEta));
	}
	return tmp;
}

vector<double> Laminar::calcu_c_coe()
{
	vector<double> tmp;
	for(int j=0;j!=N;++j){
		tmp.push_back(eta_y(j*deltaEta));
	}
	return tmp;
}

void Laminar::forward()
{
	vector<double> U_tmp, V_tmp;
	double aver_b = 0.0;

	for(int i=0; i!=M-1; ++i){	//从第一层开始推进,到倒数第二层
		bool mark = 1;
		int iter=0;

		vector<double> U_N(UU[i].size(), 0.0), V_N(VV[i].size(), 0.0);

		while(mark && iter < 100){
			mark = 0;	// mark用来判断否是收敛
			cout << "目前是第"<<i+1<<"层的第"<<iter+1<<"次迭代！！！！！"<<endl;

			for(int j=0; j!=N; ++j){
				aver_b = 0.5*(calcu_b_coe(V)[j] + calcu_b_coe(V_N)[j]);	//用这一步和下一步
				b_coe.push_back(aver_b);
			}
			assemble_momentum();	//装配矩阵,使用了上一步的a_coe和b_coe

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

void Laminar::assemble_momentum()	// 构建了A和b用于之后的矩阵求解
{
	vector<double> a;

	for(int j=0; j!=N; ++j){
		if(j==0){	//首先装配下边界矩阵
			a.push_back(1.0);
			for(int k=1;k!=N;++k)
				a.push_back(0.0);
			A.push_back(a);
			a.clear();
			b.push_back(0.0);
		}

		else if(j!=0 && j!=N-1){		//装配非边界的矩阵,使用了当前步的U和下一步的V
			for(int k=0;k!=N;++k){
				if(k == j-1)
					a.push_back(-a_coe[j]/(deltaEta*deltaEta) + b_coe[j]/(2*deltaEta));
				else if(k == j)
					a.push_back(2*U[j]/deltaX + 2*a_coe[j]/(deltaEta*deltaEta));//线性化
				else if(k == j+1)
					a.push_back(-a_coe[j]/(deltaEta*deltaEta) - b_coe[j]/(2*deltaEta));
				else
					a.push_back(0.0);
			}
			A.push_back(a);
			a.clear();
			b.push_back(2*U[j]*U[j]/deltaX
					+a_coe[j]*(U[j+1]-2*U[j]+U[j-1])/(deltaEta*deltaEta)
					+b_coe[j]*(U[j+1]-U[j-1])/(2*deltaEta));
		}

		else if(j==N-1){	//装配上边界的矩阵
			for(int k=0;k!=N-2;++k)
				a.push_back(0.0);
			a.push_back(-1.0);
			a.push_back(1.0);
			A.push_back(a);
			a.clear();
			b.push_back(0);
		}
	}
}

void Laminar::assemble_continue	//option2
(const vector<double>& U_N)
{
	vector<double> a;
	for(int j=0;j!=N;++j){
		if(j==0){
			a.push_back(1.0);
			for(int k=1;k!=N;++k)
				a.push_back(0);
			A.push_back(a);
			a.clear();
			b.push_back(0);
		}

		else{
			for(int k=0;k!=N;++k){
				if(k == j-1)
					a.push_back(-0.5*(c_coe[j]+c_coe[j-1])/deltaEta);
				else if(k == j)
					a.push_back(0.5*(c_coe[j]+c_coe[j-1])/deltaEta);
				else
					a.push_back(0.0);
			}
			A.push_back(a);
			a.clear();
			b.push_back(0.5*(U[j]-U_N[j]+U[j-1]-U_N[j-1])/deltaX);
		}
	}
}
/*
void Laminar::assemble_continue	//连续方程option3
(const vector<double>& U_N)
{
	vector<double> a;
	for(int j=0;j!=N;++j){
		if(j==0){
			a.push_back(1.0);
			for(int k=1;k!=N;++k)
				a.push_back(0);
			A.push_back(a);
			a.clear();
			b.push_back(0);
		}

		else {
			for(int k=0;k!=N;++k){
				if(k == j-1)
					a.push_back(-0.5*(c_coe[j]+c_coe[j-1])/deltaEta);
				else if(k == j)
					a.push_back(0.5*(c_coe[j]+c_coe[j-1])/deltaEta);
				else
					a.push_back(0.0);
			}
			A.push_back(a);
			a.clear();
			b.push_back(0.5*(U[j-1]+U[j]-U_N[j-1]-U_N[j])/deltaX-0.5*(c_coe[j]+c_coe[j-1])*(V[j]-V[j-1])/deltaEta);
		}
	}
}*/

vector<double> Laminar::LU(vector<vector<double> >& A, vector<double>& b) const
{
	const int n = b.size();
	double upper[n][n] = {0};
	double lower[n][n] = {0};
	double y[n] = {0};
	double x[n] = {0};
	vector<double> X;

	for(int j=0; j<n; ++j)
		upper[0][j] = A[0][j];
	for(int i=1; i<n ;++i)
		lower[i][0] = A[i][0]/upper[0][0];

	for(int k=1; k<n; ++k){		//Doolittle 分解
        for(int j=k; j<n; ++j){
            double sum = 0;
            for(int m = 0; m<k; ++m){
                sum += lower[k][m] * upper[m][j];
            }
            upper[k][j] = A[k][j] - sum;
        }

		for(int i=k+1;i<n;++i)
		{
			double sum = 0;
			for (int m=0; m<k; ++m){
				sum += lower[i][m] * upper[m][k];
			}
			lower[i][k] = (A[i][k] - sum) / upper[k][k];
		}
	}

	y[0] = b[0];

	for(int k=0; k<n; ++k){
		double sum = 0;
		for(int m=0; m<k; ++m)
			sum += lower[k][m] * y[m];
		y[k] = b[k] - sum;
	}

	x[n-1] = y[n-1] / upper[n-1][n-1];
	for(int k=n-2; k>-1; --k){
		double sum = 0;
		for(int m=k+1; m<n; ++m)
			sum += upper[k][m] * x[m];
		x[k] = (y[k]-sum) / upper[k][k];
	}

	for(int i=0; i<n; ++i)
		X.push_back(x[i]);

	return X;
}

vector<double> Laminar::ZG(vector<vector<double> >& A,vector<double>& b) const
{	
	const int n = b.size();
	vector<double> X;
	double u[n]={0};
	double l[n]={0};
	double x[n]={0};
	double y[n]={0};

	u[0] = A[0][0];
	for(int i=1;i<n;++i){
		l[i]=A[i][i-1]/u[i-1];
		u[i]=A[i][i]-l[i]*A[i-1][i];
	}

	y[0]=b[0];
	for(int i=1;i<n;++i)
		y[i]=b[i]-l[i]*y[i-1];

	x[n-1]=y[n-1]/u[n-1];
	for(int i=n-2;i>-1;--i)
		x[i]=(y[i]-A[i][i+1]*x[i+1])/u[i];

	for(int i=0;i<n;++i)
		X.push_back(x[i]);

	return X;
}

void Laminar::file_output()
{
	ofstream out;
	out.open("U_lam.dat");
	out<<"variables = x, y, U"<<endl;
	out<<"ZONE I=" << M <<", J=" << N << ",DATAPACKING=POINT"<<endl;
	for(int j=0; j!=N; ++j){
		for(int i=0; i!=M; ++i){
			out << Points_phy[i][j].first << " " 
				<<Points_phy[i][j].second << " " << UU[i][j] << endl;
		}
	}
	out.close();

	out.open("V_lam.dat");
	out<<"variables = x, y, V"<<endl;
	out<<"ZONE I=" << M <<", J=" << N << ",DATAPACKING=POINT"<<endl;
	for(int j=0; j!=N; ++j){
		for(int i=0; i!=M; ++i){
			out << Points_phy[i][j].first << " " 
				<<Points_phy[i][j].second << " " << VV[i][j] << endl;
		}
	}
	out.close();

	out.open("eta_U_025_lam.dat");
	out<< "variables = U/U∞, η" << endl;
    int k = M/4;
	for(int j=0;j!=N;++j){
		out << UU[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_U_050_lam.dat");
	out<< "variables = U/U∞, η" << endl;
	k = M/2;
	for(int j=0;j!=N;++j){
		out << UU[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_U_075_lam.dat");
	out<< "variables = U/U∞, η" << endl;
	k = 3*M/4;
	for(int j=0;j!=N;++j){
		out << UU[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_V_025_lam_analystic.dat");	//公式错了啊，是eta*f'-f 不是eta*f'-f''
	out<< "variables = V/U∞, η" << endl;
	out<< (0*0-0)/(2*sqrt(0.25*Re)) <<" "<<0<<endl;
	out<< (0.5*0.1659-0.0415)/(2*sqrt(0.25*Re)) <<" "<<0.5<<endl;
	out<< (1.0*0.3298-0.1656)/(2*sqrt(0.25*Re)) <<" "<<1.0<<endl;
	out<< (1.5*0.4868-0.3701)/(2*sqrt(0.25*Re)) <<" "<<1.5<<endl;
	out<< (2.0*0.6298-0.6500)/(2*sqrt(0.25*Re)) <<" "<<2.0<<endl;
	out<< (2.5*0.7513-0.9963)/(2*sqrt(0.25*Re)) <<" "<<2.5<<endl;
	out<< (3.0*0.8460-1.3968)/(2*sqrt(0.25*Re)) <<" "<<3.0<<endl;
	out<< (3.5*0.9130-1.8377)/(2*sqrt(0.25*Re)) <<" "<<3.5<<endl;
	out<< (4.0*0.9555-2.3057)/(2*sqrt(0.25*Re)) <<" "<<4.0<<endl;
	out<< (4.5*0.9795-2.7901)/(2*sqrt(0.25*Re)) <<" "<<4.5<<endl;
	out<< (5.0*0.9915-3.2833)/(2*sqrt(0.25*Re)) <<" "<<5.0<<endl;
	out<< (5.5*0.9969-3.7806)/(2*sqrt(0.25*Re)) <<" "<<5.5<<endl;
	out<< (6.0*0.9990-4.2796)/(2*sqrt(0.25*Re)) <<" "<<6.0<<endl;
	out<< (6.5*0.9997-4.7793)/(2*sqrt(0.25*Re)) <<" "<<6.5<<endl;
	out<< (7.0*0.9999-5.2792)/(2*sqrt(0.25*Re)) <<" "<<7.0<<endl;
	out<< (7.5*1.0000-5.7792)/(2*sqrt(0.25*Re)) <<" "<<7.5<<endl;
	out<< (8.0*1.0000-6.2792)/(2*sqrt(0.25*Re)) <<" "<<8.0<<endl;
	out.close();

	out.open("eta_V_050_lam_analystic.dat");	//公式错了啊，是eta*f'-f 不是eta*f'-f''
	out<< "variables = V/U∞, η" << endl;
	out<< (0*0-0)/(2*sqrt(0.5*Re)) <<" "<<0<<endl;
	out<< (0.5*0.1659-0.0415)/(2*sqrt(0.5*Re)) <<" "<<0.5<<endl;
	out<< (1.0*0.3298-0.1656)/(2*sqrt(0.5*Re)) <<" "<<1.0<<endl;
	out<< (1.5*0.4868-0.3701)/(2*sqrt(0.5*Re)) <<" "<<1.5<<endl;
	out<< (2.0*0.6298-0.6500)/(2*sqrt(0.5*Re)) <<" "<<2.0<<endl;
	out<< (2.5*0.7513-0.9963)/(2*sqrt(0.5*Re)) <<" "<<2.5<<endl;
	out<< (3.0*0.8460-1.3968)/(2*sqrt(0.5*Re)) <<" "<<3.0<<endl;
	out<< (3.5*0.9130-1.8377)/(2*sqrt(0.5*Re)) <<" "<<3.5<<endl;
	out<< (4.0*0.9555-2.3057)/(2*sqrt(0.5*Re)) <<" "<<4.0<<endl;
	out<< (4.5*0.9795-2.7901)/(2*sqrt(0.5*Re)) <<" "<<4.5<<endl;
	out<< (5.0*0.9915-3.2833)/(2*sqrt(0.5*Re)) <<" "<<5.0<<endl;
	out<< (5.5*0.9969-3.7806)/(2*sqrt(0.5*Re)) <<" "<<5.5<<endl;
	out<< (6.0*0.9990-4.2796)/(2*sqrt(0.5*Re)) <<" "<<6.0<<endl;
	out<< (6.5*0.9997-4.7793)/(2*sqrt(0.5*Re)) <<" "<<6.5<<endl;
	out<< (7.0*0.9999-5.2792)/(2*sqrt(0.5*Re)) <<" "<<7.0<<endl;
	out<< (7.5*1.0000-5.7792)/(2*sqrt(0.5*Re)) <<" "<<7.5<<endl;
	out<< (8.0*1.0000-6.2792)/(2*sqrt(0.5*Re)) <<" "<<8.0<<endl;
	out.close();

	out.open("eta_V_075_lam_analystic.dat");	//公式错了啊，是eta*f'-f 不是eta*f'-f''
	out<< "variables = V/U∞, η" << endl;
	out<< (0*0-0)/(2*sqrt(0.75*Re)) <<" "<<0<<endl;
	out<< (0.5*0.1659-0.0415)/(2*sqrt(0.75*Re)) <<" "<<0.5<<endl;
	out<< (1.0*0.3298-0.1656)/(2*sqrt(0.75*Re)) <<" "<<1.0<<endl;
	out<< (1.5*0.4868-0.3701)/(2*sqrt(0.75*Re)) <<" "<<1.5<<endl;
	out<< (2.0*0.6298-0.6500)/(2*sqrt(0.75*Re)) <<" "<<2.0<<endl;
	out<< (2.5*0.7513-0.9963)/(2*sqrt(0.75*Re)) <<" "<<2.5<<endl;
	out<< (3.0*0.8460-1.3968)/(2*sqrt(0.75*Re)) <<" "<<3.0<<endl;
	out<< (3.5*0.9130-1.8377)/(2*sqrt(0.75*Re)) <<" "<<3.5<<endl;
	out<< (4.0*0.9555-2.3057)/(2*sqrt(0.75*Re)) <<" "<<4.0<<endl;
	out<< (4.5*0.9795-2.7901)/(2*sqrt(0.75*Re)) <<" "<<4.5<<endl;
	out<< (5.0*0.9915-3.2833)/(2*sqrt(0.75*Re)) <<" "<<5.0<<endl;
	out<< (5.5*0.9969-3.7806)/(2*sqrt(0.75*Re)) <<" "<<5.5<<endl;
	out<< (6.0*0.9990-4.2796)/(2*sqrt(0.75*Re)) <<" "<<6.0<<endl;
	out<< (6.5*0.9997-4.7793)/(2*sqrt(0.75*Re)) <<" "<<6.5<<endl;
	out<< (7.0*0.9999-5.2792)/(2*sqrt(0.75*Re)) <<" "<<7.0<<endl;
	out<< (7.5*1.0000-5.7792)/(2*sqrt(0.75*Re)) <<" "<<7.5<<endl;
	out<< (8.0*1.0000-6.2792)/(2*sqrt(0.75*Re)) <<" "<<8.0<<endl;
	out.close();

	out.open("eta_V_025_lam.dat");
	out<< "variables = V/U∞, η" << endl;
	k = M/4;
	for(int j=0;j!=N;++j){
		out << VV[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_V_050_lam.dat");
	out<< "variables = V/U∞, η" << endl;
	k = M/2;
	for(int j=0;j!=N;++j){
		out << VV[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("eta_V_075_lam.dat");
	out<< "variables = V/U∞, η" << endl;
	k = 3*M/4;
	for(int j=0;j!=N;++j){
		out << VV[k][j] <<" "<< Points_phy[k][j].second 
			* sqrt(X_length * Re/Points_phy[k][j].first) << endl;
	}
	out.close();

	out.open("layer_thickness_compute_lam.dat");
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

	out.open("layer_thickness_analystic_lam.dat");
	out << "variables = X, Thickness" << endl;
	for(double x=0.01; x<X_length; x += X_length/M)
		out << x << " " << x*5.0/sqrt(Re*x/X_length) << endl;
	out.close();

	out.open("Cf_compute_lam.dat");
	out << "variables = X, Cf" << endl;
	for(int i=0; i!=M; ++i){
		out << Points_phy[i][0].first << " " << (UU[i][1] - UU[i][0])/deltaEta * eta_y(deltaEta) * X_length/(0.5*Re) << endl;
	}
	out.close();

	out.open("Cf_analystic_lam.dat");
	out << "variables = X, Cf" << endl;
	for(double x=0.01; x<X_length; x += 2*X_length/M)
		out << x << " " << 0.664/sqrt(Re*x/X_length) << endl;
	out.close();
}
