#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <fstream>

using namespace std;

void construct_U(vector<vector<double> >& U, const vector<vector<double> >& S)
{
	int n = S.size(), m = S[0].size(), k;
	int q = m * (m - 1) / 2;
	int p = (m + 1)*(m + 2) / 2;

	vector<double> u;

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < p; ++j) {
			if (j == 0)
				u.push_back(1);
			else if (j < m + 1)
				u.push_back(S[i][j - 1]);
			else if (j < 1 + m + q) {
				k = j - 1 - m;
				for (int d = 1; d < m - k; ++d) {
					u.push_back(S[i][k] * S[i][k + d]);
				}
			}
			else
				u.push_back(S[i][j - 1 - m - q] * S[i][j - 1 - m - q]);
		}
		U.push_back(u);
		u.clear();
	}
}

void trans_U(vector<vector<double> > &UT, const vector<vector<double> >& U)
{
	int n = U.size(), p = U[0].size();

	vector<double> ut;

	for (int i = 0; i < p; ++i) {
		for (int j = 0; j < n; ++j) {
			ut.push_back(U[j][i]);
		}
		UT.push_back(ut);
		ut.clear();
	}
}

vector<vector<double> > multiply(const vector<vector<double> >&A, const vector<vector<double> >&B)
{
	double sum = 0;
	vector<double> times;
	vector<vector<double> > TIMES;

	int m = A.size(), n = B[0].size(), l=A[0].size();

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			sum = 0;
			for (int k = 0; k < l; ++k) {
				sum += A[i][k] * B[k][j];
			}
			times.push_back(sum);
		}
		TIMES.push_back(times);
		times.clear();
	}
	return TIMES;
}

vector<double> calcu_b(const vector<vector<double> >& UT, const vector<double>&ys)
{
	vector<double> b;
	double sum = 0;

	int m = UT.size(), n=ys.size();

	for (int i = 0; i < m; ++i) {
		sum = 0;
		for (int k = 0; k < n; ++k) {
			sum += UT[i][k] * ys[k];
		}
		b.push_back(sum);
	}
	return b;
}

vector<double> LU(vector<vector<double> >& A, vector<double>& b)
{
	const int n = b.size();

	cout << "求解方程组的阶为" << n << endl;

	double **upper = new double *[n];
	double **lower = new double *[n];
	double *x = new double[n];
	double *y = new double[n];

	for (int i = 0; i < n; ++i) {
		upper[i] = new double[n];
		lower[i] = new double[n];
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			upper[i][j] = 0;
			lower[i][j] = 0;
		}
	}

	cout << "upper矩阵为" << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << upper[i][j] << " ";
		}
		cout << endl;
	}
	cout << "lower矩阵为" << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << lower[i][j] << " ";
		}
		cout << endl;
	}

	vector<double> X;

	for (int j = 0; j<n; ++j)
		upper[0][j] = A[0][j];
	for (int i = 1; i<n; ++i)
		lower[i][0] = A[i][0] / upper[0][0];

	for (int k = 1; k<n; ++k) {		//Doolittle 分解
		for (int j = k; j<n; ++j) {
			double sum = 0;
			for (int m = 0; m<k; ++m) {
				sum += lower[k][m] * upper[m][j];
			}
			upper[k][j] = A[k][j] - sum;
		}

		for (int i = k + 1; i<n; ++i)
		{
			double sum = 0;
			for (int m = 0; m<k; ++m) {
				sum += lower[i][m] * upper[m][k];
			}
			lower[i][k] = (A[i][k] - sum) / upper[k][k];
		}
	}

	cout << "upper矩阵为" << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << upper[i][j] << " ";
		}
		cout << endl;
	}
	cout << "lower矩阵为" << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << lower[i][j] << " ";
		}
		cout << endl;
	}

	y[0] = b[0];

	for (int k = 0; k<n; ++k) {
		double sum = 0;
		for (int m = 0; m<k; ++m)
			sum += lower[k][m] * y[m];
		y[k] = b[k] - sum;
	}

	x[n - 1] = y[n - 1] / upper[n - 1][n - 1];
	for (int k = n - 2; k>-1; --k) {
		double sum = 0;
		for (int m = k + 1; m<n; ++m)
			sum += upper[k][m] * x[m];
		x[k] = (y[k] - sum) / upper[k][k];
	}

	for (int i = 0; i<n; ++i)
		X.push_back(x[i]);

	for (int i = 0; i < n; ++i) {
		delete[]upper[i];
		delete[]lower[i];
	}
	delete[]lower;
	delete[]upper;
	delete[]x;
	delete[]y;

	return X;
}

double respond(const vector<double>& X, const vector<double>& beta)
{
	double y=0, sum1=0, sum2=0, sum3=0;
	int m = X.size(), k;
	int q = m * (m - 1) / 2;
	int p = (m + 1)*(m + 2) / 2;

	//求响应值的第一种方法
	//vector<double> u;

	//for (int j = 0; j < p; ++j) {
	//	if (j == 0)
	//		u.push_back(1);
	//	else if (j < m + 1)
	//		u.push_back(X[j - 1]);
	//	else if (j < 1 + m + q) {
	//		k = j - 1 - m;
	//		for (int d = 1; d < m - k; ++d) {
	//			u.push_back(X[k] * X[k + d]);
	//		}
	//	}
	//	else
	//		u.push_back(X[j - 1 - m - q] * X[j - 1 - m - q]);
	//}

	//for (int i = 0; i < p; ++i) {
	//	y += u[i] * beta[i];
	//}

	//求响应值的第二种方法, 其实都一样
	for (int i = 0; i < m; ++i) {
		sum1 += beta[i+1] * X[i];
		for (int j = 1; j < m - i; ++j) {
			sum2 += beta[i + 1 + m] * X[i]*X[i+j];
		}
		sum3 += beta[i + 1 + m + q] * X[i] * X[i];
	}
	y = beta[0] + sum1 + sum2 + sum3;

	return y;
}

int main()
{
	int m=2, n=100;	//m个设计变量,n个样本点
	vector<vector<double> > S, U, UT, U_times;
	vector<double> X, beta, ys, B, s_temp;

	ofstream out;

	uniform_int_distribution<int> u(-2000000, 2000000);
	default_random_engine e;

	double y_respond, sum=0;

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			s_temp.push_back(double(u(e)/1000000.0));
		}
		S.push_back(s_temp);
		s_temp.clear();
	}

	cout << "样本点为:" << endl;
	for (const auto &a : S) {
		for (const auto b : a) {
			cout << b << " ";
		}
		cout << endl;
	}

	for (int i = 0; i < n; ++i) {
		//ys.push_back(double(u(e) / 1000000.0));
		//sum = 0;
		//for (int j = 0; j < m; ++j) {
		//	sum += pow(S[i][j], 2);
		//}
		//ys.push_back(sum);

		ys.push_back(pow(1 - S[i][0], 2) + 100 * pow(S[i][1] - pow(S[i][0], 2), 2));
	}

	out.open("random.dat");
	for (int i = 0; i < int(S.size()); ++i) {
		for (int j = 0; j < int(S[0].size()); ++j) {
			out << S[i][j] << " ";
		}
		out << ys[i] << endl;
	}
	out.close();

	cout << "样本点响应值为:" << endl;
	for (const auto a : ys)
		cout << a << endl;

	construct_U(U, S);

	cout << "UT为: " << endl;
	for (const auto &a : U) {
		for (const auto b : a) {
			cout << b << " ";
		}
		cout << endl;
	}

	trans_U(UT,U);

	cout << "UT为: " << endl;
	for (const auto &a : UT) {
		for (const auto b : a) {
			cout << b << " ";
		}
		cout << endl;
	}

	U_times = multiply(UT,U);

	cout << "UT*U，即A为: " << endl;
	for (const auto &a : U_times) {
		for (const auto b : a) {
			cout << b << " ";
		}
		cout << endl;
	}

	B = calcu_b(UT, ys);

	cout << "组建方程的B为：" << endl;
	for (const auto x : B)
		cout << x << " ";
	cout << endl;

	beta = LU(U_times, B);

	cout << "求得beta为：" << endl;
	for (const auto x : beta)
		cout << x << " ";
	cout << endl;

	vector<vector<vector<double> > > SPACE;
	vector<vector<double> >Space;
	vector<double>space;

	//三维测试
	//for (double x1 = -1.0; x1 <= 1.01; x1 += 0.02) {
	//	for (double x2 = -1.0; x2 <= 1.01; x2 += 0.02) {
	//		for (double x3 = -1.0; x3 <= 1.01; x3 += 0.02) {
	//			space.push_back(x1);
	//			space.push_back(x2);
	//			space.push_back(x3);
	//			Space.push_back(space);
	//			space.clear();
	//		}
	//	}
	//	SPACE.push_back(Space);
	//	Space.clear();
	//}

	for (double x1 = -2.0; x1 <= 2.01; x1 += 0.05) {
		for (double x2 = -2.0; x2 <= 2.01; x2 += 0.05) {
			space.push_back(x1);
			space.push_back(x2);
			Space.push_back(space);
			space.clear();
		}
		SPACE.push_back(Space);
		Space.clear();
	}

	for (auto &c : SPACE) {
		for (auto &d : c) {
			cout << "<";
			for (auto e : d) {
				cout << e << " ";
			}
			cout << ">" << endl;
		}
		cout << endl;
	}

	//单样本测试
	//for (int i = 0; i < m; ++i)
	//	X.push_back(double(u(e) / 1000000.0));

	//y_respond = respond(X, beta);

	//cout << y_respond << endl;

	int nn = SPACE.size();
	int mm = SPACE[0].size();

	out.open("QRSM.dat");
	out << "variables=x1, x2, respond" << endl;
	for (int i = 0; i < nn; ++i) {
		for (int j = 0; j < mm; ++j) {
			for (int k = 0; k < m; ++k) {
				out << SPACE[i][j][k] << " ";
			}
			out << respond(SPACE[i][j], beta) << endl;
		}
	}
	out.close();

	out.open("real.dat");
	out << "variables=x1, x2, respond" << endl;
	for (int i = 0; i < nn; ++i) {
		for (int j = 0; j < mm; ++j) {
			//sum = 0.0;
			for (int k = 0; k < m; ++k) {
				out << SPACE[i][j][k] << " ";
				/*sum += pow(SPACE[i][j][k], 2);*/
			}
			sum = pow(1 - SPACE[i][j][0], 2) + 100 * pow(SPACE[i][j][1] - pow(SPACE[i][j][0], 2), 2);
			out << sum << endl;
		}
	}
	out.close();

	system("pause");

	return 0;
}