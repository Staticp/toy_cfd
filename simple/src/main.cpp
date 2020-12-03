/*************************************************************************
    > File Name: main.cpp
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年01月17日 星期四 22时47分11秒
 ************************************************************************/

#include <iostream>
#include "simple.h"

using namespace std;

int main()
{
	simple Couette(0,0,0,0);
	Couette.loop();
	Couette.foutput();
	return 0;
}
