#include<iostream>
#include<cmath>
#define esp 1e-9
#define N 3  
using namespace std;

void Relaxation(double a[][N], double y[], int n);
double norm_inf(double x[], int n);
void Jacobi(double a[][N], double y[], int n);
void Gauss_Seidel(double a[N][N], double y[N], int n);

int main()
{	
	double a[N][N] = { { 8, -3, 2 }, { 4, 11, -1 }, { 6, 3, 12 } };
	double y[N] = { 20, 33, 36 };
	Jacobi(a, y, 3);
	cout<<endl;
	Gauss_Seidel(a, y, 3);
	cout<<endl;
	Relaxation(a, y, 3);
	cout<<endl;
	return 0;
}

double norm_inf(double x[], int n)  
{
	double max = 0;
	for (int i = 0; i < n; i++)
	{
		double tem = 0;
		for (int j = 0; j < n; j++)
			tem = tem + x[i*n + j];
		if (tem>max)
			max = tem;
	}
	return max;
}

void Jacobi(double a[N][N], double y[N], int n)
{
	double B[N][N]; double y_[N]; double xk[N], xk1[N], temp; 
	int k = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i==j)
				B[i][j] = 0;
			else
				B[i][j] = -(a[i][j] / a[i][i]);
	for (int i = 0; i < n; i++)
	{
		y_[i] = y[i] / a[i][i];
		xk[i] = 0;
		xk1[i] = 1;
	}
	while (fabs(norm_inf(xk1,n)-norm_inf(xk,n))>esp)
	{
		k++;
		for (int i = 0; i < N; i++)
		{
			xk[i] = xk1[i];
		}
		for (int i = 0; i < N; i++)
		{
			temp = 0;
			for (int j = 0; j < N; j++)
			{
				temp+= B[i][j] * xk[j];
			}
			xk1[i] = temp+y_[i];
		}
	}
	cout << "Jacobi迭代法解："<<endl;
	for (int i = 0; i < N; i++)
		cout << "x" << i + 1 <<"="<< xk1[i] << endl;
	cout << "迭代次数：" << k << endl;
}

void Gauss_Seidel(double a[N][N], double y[N], int n)
{
	double B[N][N]; double y_[N]; double xk[N], xk1[N];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				B[i][j] = 0;
			else
				B[i][j] = -(a[i][j] / a[i][i]);
	for (int i = 0; i < n; i++)
	{
		y_[i] = y[i] / a[i][i];
		xk[i] = 0;
		xk1[i] = 1;
	}
	double temp; int k = 0;
	while (fabs(norm_inf(xk1, n) - norm_inf(xk, n))>esp)
	{
		k++;
		for (int i = 0; i < N; i++)
			xk[i] = xk1[i];
		for (int i = 0; i < N; i++)
		{
			temp = 0;
			for (int j = 0; j < i; j++)
				temp += B[i][j] * xk1[j];							
			for (int j = i; j < N; j++)
				temp += B[i][j] * xk[j];
			xk1[i] = temp + y_[i];
		}
	}
	cout << "高斯赛德尔迭代法解："<<endl;
	for (int i = 0; i < N; i++)
		cout << "x" << i + 1 << "=" << xk1[i] << endl;
	cout << "迭代次数：" << k << endl;
	cout<<endl;
}


void Relaxation(double a[N][N], double y[N], int n) 
{
	double B[N][N]; double y_[N]; double xk[N], xk1[N];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				B[i][j] = 0;
			else
				B[i][j] = -(a[i][j] / a[i][i]);
	for (int i = 0; i < n; i++)
	{
		y_[i] = y[i] / a[i][i];
		xk[i] = 0;
		xk1[i] = 1;
	}
	double temp; int k = 0; double w = 0.99;
	while (fabs(norm_inf(xk1, n) - norm_inf(xk, n))>esp)
	{
		k++;
		for (int i = 0; i < N; i++)
			xk[i] = xk1[i];
		for (int i = 0; i < N; i++)
		{
			temp = 0;
			for (int j = 0; j < i; j++)
				temp += B[i][j] * xk1[j];							
			for (int j = i; j < N; j++)
				temp += B[i][j] * xk[j];
			xk1[i] = w*(temp + y_[i])+(1-w)*xk[i];
		}
	}
	cout << "松弛迭代法解为："<<endl;
	for (int i = 0; i < N; i++)
		cout << "x" << i + 1 << "=" << xk1[i] << endl;
	cout << "迭代次数：" << k << endl;
}

