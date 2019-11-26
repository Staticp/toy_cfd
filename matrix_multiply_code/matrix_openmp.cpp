/*************************************************************************
    > File Name: matrix.cpp
    > Author: static
    > Mail: wxtchina@foxmail.com 
    > Created Time: 2019年04月07日 星期日 19时05分16秒
 ************************************************************************/

#include<iostream>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<ctime>
#include<random>
#include<omp.h>

using namespace std;

int main()
{
	int i,j,k;
	int M,N,L;
    double t_start, t_end;
	uniform_int_distribution<int> u(-1000000,1000000);
	default_random_engine e;

	cout << "calculate A*B" << endl;
	cout << "please input dimensions of A(M*L)" <<endl;
	cin >> M >> L;
	cout << "please input dimensions of B(L*N)" <<endl;
	cin >> L >> N;

	double **A = new double* [M];
	for(i=0;i<M;++i)
		A[i] = new double[L];

	double **B = new double* [L];
	for(i=0;i<L;++i)
		B[i] = new double[N];

	double **BT = new double* [N];
	for(i=0;i<N;++i)
		BT[i] = new double[L];

	double **C = new double* [M];
	for(i=0;i<M;++i)
		C[i] = new double[N];

//	srand((unsigned)time(NULL));
	for(i=0;i<M;++i)
		for(j=0;j<L;++j)
			A[i][j] = double(u(e)/100000.0);
//			A[i][j] = rand()/(double)(RAND_MAX/100)-50;

	for(i=0;i<L;++i){
		for(j=0;j<N;++j)
			B[i][j] = double(u(e)/100000.0);
//			B[i][j] = rand()/(double)(RAND_MAX/100)-50;
	}

	for(i=0;i<N;++i)
		for(j=0;j<L;++j)
			BT[i][j] = B[j][i];

	for(i=0;i<M;++i)
		for(j=0;j<N;++j)
			C[i][j] = 0;
	
	omp_set_num_threads(8);

	t_start = omp_get_wtime();

	#pragma omp parallel for shared(A,BT,C) private(i,j,k)
	for(i=0;i<M;++i)
		for(j=0;j<N;++j){
			for(k=0;k<L;++k)
				C[i][j] += A[i][k]*BT[j][k];
		}
	t_end = omp_get_wtime();

	cout <<"求解时间为"<< (double)(t_end-t_start) << "s"<<endl;

	ofstream out;
	out.open("openmp_result.dat");

	out <<"求解时间为"<< (double)(t_end-t_start) << "s"<<endl;
	out<<"A矩阵为"<<endl;
	for(i=0;i<M;++i){
		for(j=0;j<L;++j)
			out <<A[i][j]<<" ";
		out<<endl;
	}

	out<<"B矩阵为"<<endl;
	for(i=0;i<L;++i){
		for(j=0;j<N;++j)
			out <<B[i][j]<<" ";
		out<<endl;
	}

	out<<"相乘得C为"<<endl;
	for(i=0;i<M;++i){
		for(j=0;j<N;++j)
			out <<C[i][j]<<" ";
		out <<endl;
	}
	out.close();

	for(i=0;i<M;++i)
		delete []A[i];
	delete [] A;

	for(i=0;i<L;++i)
		delete []B[i];
	delete [] B;

	for(i=0;i<N;++i)
		delete []BT[i];
	delete [] BT;

	for(i=0;i<M;++i)
		delete []C[i];
	delete [] C;

	return 0;
}
