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
#include<mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
	MPI_Status status;
	MPI_Init(&argc, &argv);

	int numprocs, myid;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	int i,j,k;
	int M,N,L;
	clock_t t_start, t_end;
	uniform_int_distribution<int> u(-1000000,1000000);
	default_random_engine e;
	
	if(myid==0)
	{	
		cout << "calculate A*B" << endl;
		cout << "please input dimensions of A(M*L)" <<endl;
		cin >> M >> L;
		cout << "please input dimensions of B(L*N)" <<endl;
		cin >> L >> N;

	}
	MPI_Bcast(&M,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&L,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&N,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	double *A = new double [M*L];
	double *B = new double [L*N];
	double *BT = new double [N*L];
	double *C = new double [M*N];

	double *a, *c;

//	srand((unsigned)time(NULL));
	for(i=0;i<M*L;++i)
		A[i] = double(u(e)/100000.0);

	for(i=0;i<L*N;++i)
		B[i] = double(u(e)/100000.0);

	for(i=0;i<N;++i)
        for(j=0;j<L;++j)
    		BT[i*L+j] = B[j*N+i];

	for(i=0;i<M*N;++i)
		C[i] = 0;


	int m=M/numprocs;

	a = new double[m*L];
	c = new double[m*N];

if(myid==0)
	t_start = clock();
	MPI_Scatter(A,m*L,MPI_DOUBLE,a,m*L,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(B,L*N,MPI_DOUBLE,0,MPI_COMM_WORLD);

	for(i=0;i<m;++i)
		for(j=0;j<N;++j){
			double sum=0;
			for(k=0;k<L;++k){
				sum += a[i*L+k]*BT[j*L+k];
			}
			c[i*N+j]=sum;
		}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(c,m*N,MPI_DOUBLE,C,m*N,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int remain = m*numprocs;
	
	if(myid==0){
		for(i=remain;i<M;++i)
			for(j=0;j<N;++j){
				double sum=0;
				for(k=0;k<L;++k){
					sum += A[i*L+k]*BT[j*L+k];
				}
				C[i*N+j]=sum;
			}
	t_end = clock();
	}


	if(myid==0){
		cout <<"求解时间为"<< (double)(t_end-t_start)/CLOCKS_PER_SEC << "s"<<endl;

		ofstream out;
		out.open("mpi_result.dat");
		out <<"求解时间为"<< (double)(t_end-t_start)/CLOCKS_PER_SEC << "s"<<endl;
		out<<"A矩阵为"<<endl;
		for(i=0;i<M;++i){
			for(j=0;j<L;++j)
				out <<A[i*L+j]<<" ";
			out<<endl;
		}

		out<<"B矩阵为"<<endl;
		for(i=0;i<L;++i){
			for(j=0;j<N;++j)
				out <<B[i*N+j]<<" ";
			out<<endl;
		}

		out<<"相乘得C为"<<endl;
		for(i=0;i<M;++i){
			for(j=0;j<N;++j)
				out <<C[i*N+j]<<" ";
			out<<endl;
		}

		out.close();

	}
	
	delete [] A;
	delete [] B;
	delete [] BT;
	delete [] C;
	delete [] a;
	delete [] c;

	MPI_Finalize();
	return 0;
}
