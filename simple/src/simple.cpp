/*************************************************************************
    > File Name: simple.cpp
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年01月15日 星期二 23时24分22秒
 ************************************************************************/

#include <cmath>
#include <fstream>
#include "simple.h"

void simple::initialize(double p0=0, double pd0=0, double u0=0, double v0=0)
{
	for(int i=0; i!=Mp; ++i){
		p.clear();
		for(int j=0; j!=Np; ++j){
			p.push_back(p0);
		}
		p_.push_back(p);
	}

	for(int i=0; i!=Mp; ++i){
		pd.clear();
		for(int j=0; j!=Np; ++j){
			pd.push_back(pd0);
		}
		pd_.push_back(pd);
	}

	for(int i=0; i!=Mu; ++i){
		u.clear();
		rhoU.clear();
		for(int j=0; j!=Nu; ++j){
			if(i==Mu-1){
				u.push_back(0.3048);
				rhoU.push_back(rho*0.3048);
			}
			else{
				u.push_back(u0);
				rhoU.push_back(rho*u0);
			}
		}
		u_.push_back(u);
		rhoU_.push_back(rhoU);
	}

	for(int i=0; i!=Mv; ++i){
		v.clear();
		rhoV.clear();
		for(int j=0; j!=Nv; ++j){
			if(i==4&&j==14){
				v.push_back(0.1524);
				rhoV.push_back(rho*0.1524);
			}
			else{
				v.push_back(v0);
				rhoV.push_back(rho*v0);
			}
		}
		v_.push_back(v);
		rhoV_.push_back(rhoV);
	}

	cout <<"initializal p field is:"<<endl;
	for(int i=0; i!=Mp; ++i){
		for(int j=0; j!=Np; ++j){
			cout << p_[i][j]<<" ";
		}
		cout <<endl;
	}

	cout <<"initializal pd field is:"<<endl;
	for(int i=0; i!=Mp; ++i){
		for(int j=0; j!=Np; ++j){
			cout << pd_[i][j]<<" ";
		}
		cout <<endl;
	}

	cout <<"initializal u field is:"<<endl;
	for(int i=0; i!=Mu; ++i){
		for(int j=0; j!=Nu; ++j){
			cout << u_[i][j]<<" ";
		}
		cout <<endl;
	}

	cout <<"initializal v field is:"<<endl;
	for(int i=0; i!=Mv; ++i){
		for(int j=0; j!=Nv; ++j){
			cout << v_[i][j]<<" ";
		}
		cout <<endl;
	}
}

simple::simple()
{
	initialize();
	cout << "0 parametre Constructor is used!" << endl;
}

simple::simple(double p0, double pd0, double u0, double v0)
{
	initialize(p0, pd0, u0, v0);
	cout << "4 parametres Constructor is used!" << endl;
}

void simple::momentum_predict()
{
	double rhoU_star = 0;
	double rhoV_star = 0;
	double u_star = 0;
	double v_star = 0;

	for(int i=0; i!=Mu; ++i){
		rhoU.clear();
		u.clear();
		for(int j=0; j!=Nu; ++j){
			if(i==0){
				rhoU_star=rho*0;
			}
			else if(i==Mu-1){
				rhoU_star=rho*0.3048;
			}
			else if(i!=0 && i!=Mu-1 && j!=0 && j!=Nu-1){
				rhoU_star=rhoU_[i][j]
				+(
					-0.5*(rhoU_[i][j+1]*u_[i][j+1] - rhoU_[i][j-1]*u_[i][j-1])/deltaX
					-0.25*(rhoU_[i+1][j]*(v_[i+1][j]+v_[i+1][j+1]) - rhoU_[i-1][j]*(v_[i][j]+v_[i][j+1]))/deltaY
					+mu*((u_[i][j+1]-2*u_[i][j]+u_[i][j-1])/(deltaX*deltaX)
					+(u_[i+1][j]-2*u_[i][j]+u_[i-1][j])/(deltaY*deltaY))
				)*deltaT
				-deltaT/deltaX*(p_[i][j]-p_[i][j-1]);
			}
			else if(j==0){
				rhoU_star=rhoU_[i][j+1]
				+(
					-0.5*(rhoU_[i][j+2]*u_[i][j+2] - rhoU_[i][j]*u_[i][j])/deltaX
					-0.25*(rhoU_[i+1][j+1]*(v_[i+1][j+1]+v_[i+1][j+2]) - rhoU_[i-1][j+1]*(v_[i][j+1]+v_[i][j+2]))/deltaY
					+mu*((u_[i][j+2]-2*u_[i][j+1]+u_[i][j])/(deltaX*deltaX)
					+(u_[i+1][j+1]-2*u_[i][j+1]+u_[i-1][j+1])/(deltaY*deltaY))
				)*deltaT
				-deltaT/deltaX*(p_[i][j+1]-p_[i][j]);
			}
			else if(j==Nu-1){
				rhoU_star=rhoU_[i][j-1]
				+(
					-0.5*(rhoU_[i][j]*u_[i][j] - rhoU_[i][j-2]*u_[i][j-2])/deltaX
					-0.25*(rhoU_[i+1][j-1]*(v_[i+1][j-1]+v_[i+1][j]) - rhoU_[i-1][j-1]*(v_[i][j-1]+v_[i][j]))/deltaY
					+mu*((u_[i][j]-2*u_[i][j-1]+u_[i][j-2])/(deltaX*deltaX)
					+(u_[i+1][j-1]-2*u_[i][j-1]+u_[i-1][j-1])/(deltaY*deltaY))
				)*deltaT
				-deltaT/deltaX*(p_[i][j-1]-p_[i][j-2]);
			}

			rhoU.push_back(rhoU_star);
			u.push_back(rhoU_star/rho);
		}
		rhoU_N.push_back(rhoU);
		u_N.push_back(u);
	}

	for(int i=0; i!=Mv; ++i){
		rhoV.clear();
		v.clear();
		for(int j=0; j!=Nv; ++j){
			if(i==0){
				rhoV_star=rho*0;
			}
			else if(i==Mv-1){
				rhoV_star=rho*0;
			}
			else if(i!=0 && i!=Mv-1 && j!=0 && j!=Nv-1){
				rhoV_star=rhoV_[i][j]+
				+(
					-0.25*(rhoV_[i][j+1]*(u_[i-1][j]+u_[i][j]) - rhoV_[i][j-1]*(u_[i-1][j-1]+u_[i][j-1]))/deltaX
					-0.5*(rhoV_[i+1][j]*v_[i+1][j]-rhoV_[i-1][j]*v_[i-1][j])/deltaY
					+mu*((v_[i][j+1]-2*v_[i][j]+v_[i][j-1])/(deltaX*deltaX)
					+(v_[i+1][j]-2*v_[i][j]+v_[i-1][j])/(deltaY*deltaY))
				)*deltaT
				-deltaT/deltaY*(p_[i][j-1]-p_[i-1][j-1]);
			}
			else if(j==0){
				rhoV_star=rho*0;
			}
			else if(j==Nv-1){
				rhoV_star=rhoV_[i][j-1]+
				+(
					-0.25*(rhoV_[i][j]*(u_[i-1][j-1]+u_[i][j-1]) - rhoV_[i][j-2]*(u_[i-1][j-2]+u_[i][j-2]))/deltaX
					-0.5*(rhoV_[i+1][j-1]*v_[i+1][j-1]-rhoV_[i-1][j-1]*v_[i-1][j-1])/deltaY
					+mu*((v_[i][j]-2*v_[i][j-1]+v_[i][j-2])/(deltaX*deltaX)
					+(v_[i+1][j-1]-2*v_[i][j-1]+v_[i-1][j-1])/(deltaY*deltaY))
				)*deltaT
				-deltaT/deltaY*(p_[i][j-2]-p_[i-1][j-2]);
			}
			rhoV.push_back(rhoV_star);
			v.push_back(rhoV_star/rho);
		}
		rhoV_N.push_back(rhoV);
		v_N.push_back(v);
	}

	rhoU_=rhoU_N;
	rhoV_=rhoV_N;
	u_=u_N;
	v_=v_N;
	u_N.clear();
	v_N.clear();
	rhoU_N.clear();
	rhoV_N.clear();
}

void simple::pd_relax_equation()
{
	double pd_star = 0;
	double b = -deltaT/(deltaX*deltaX);
	double c = -deltaT/(deltaY*deltaY);
	double a = -2*(b + c);
	bool state = true;
	int iter_num = 1;

	while(state && iter_num < P_CORR_COUNT){
		state = false;
		for(int i=0; i!=Mp; ++i){
			pd.clear();
			for(int j=0; j!=Np; ++j){
				if(i!=0 && i!=Mp-1 && j!=0 && j!=Np-1){
					pd_star=-1/a * (b*(pd_[i][j+1]+pd_[i][j-1])+c*(pd_[i+1][j]+pd_[i-1][j])
						+ (rhoU_[i][j+1]-rhoU_[i][j])/deltaX + (rhoV_[i+1][j+1]-rhoV_[i][j+1])/deltaY);
				}
				else{
					pd_star=0;
				}

				if(fabs(pd_star - pd_[i][j]) > PD_TOL){
					state = true;
				}
				pd.push_back(pd_star);
			}
			pd_N.push_back(pd);
		}
		pd_=pd_N;
		++iter_num;
	}
	pd_N.clear();
}

void simple::p_correction()
{
	double p_corr = 0;

	for(int i=0; i!=Mp; ++i){
		p.clear();
		for(int j=0; j!=Np; ++j){
			p_corr = p_[i][j] + alpha_p * pd_[i][j];
			p.push_back(p_corr);
		}
		p_N.push_back(p);
	}
	p_=p_N;
	p_N.clear();
}

void simple::loop()
{
	int iter_num=1;
	while(iter_num < ITER_MAX){
		cout << "begin!" << endl;
		momentum_predict();
		cout << "momentum predict complete!" << endl;

	cout <<" rhou field is:"<<endl;
	for(int i=0; i!=Mu; ++i){
		for(int j=0; j!=Nu; ++j){
			cout << rhoU_[i][j]<<" ";
		}
		cout <<endl;
	}

	cout <<" rhov field is:"<<endl;
	for(int i=0; i!=Mv; ++i){
		for(int j=0; j!=Nv; ++j){
			cout << rhoV_[i][j]<<" ";
		}
		cout <<endl;
	}

		pd_relax_equation();
		cout << "p'relax equation compute completed!" << endl;
		p_correction();
		cout << "pressure_correction complete!" << endl;
		++iter_num;
		cout << "loop count is " << iter_num << endl;
	}
}

void simple::foutput()
{
	ofstream out;
	out.open("result.dat");
	out << "variables=x,y,p"<< endl;
	for(int i=0; i!=Mp; ++i){
		for(int j=0; j!=Np; ++j){
			out << j*deltaX << "  " << i*deltaY << "  " << u_[i][j] << endl;
		}
	}
	out.close();
}
