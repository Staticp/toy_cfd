/*************************************************************************
    > File Name: main.cpp
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年03月23日 星期六 20时50分42秒
 ************************************************************************/

#include<iostream>
#include "Turbulence.h"

using namespace std;

int main()
{
	cout << "laminar boundary layer computing begin!" << endl;

	Laminar lam(1.0,0.2,200,200);
	lam.phy_fout();
	lam.com_fout();
	lam.forward();
	lam.file_output();

	cout << "turbulent boundary layer computing begin!" << endl;

	Turbulence tur(1.0,0.2,200,200);
	tur.phy_fout();
	tur.com_fout();
	tur.forward();
	tur.file_output();
	return 0;
}
