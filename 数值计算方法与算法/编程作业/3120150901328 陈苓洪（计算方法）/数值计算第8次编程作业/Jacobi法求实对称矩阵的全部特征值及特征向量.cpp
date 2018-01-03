#include <iostream>
#include <iomanip>
using namespace std;
int main()
{
	double p[3][4],x[3],z[3][4],eps,a, b, c, s = 0, q ;
	int i, j;
	p[0][0] = 2; p[1][0] = p[2][0] = p[2][1] = 1; p[1][2] = p[0][1] = p[0][2] = -1;
	p[1][1] = 5; p[0][3] = -5; p[1][3] = 8; p[2][3] = 11; p[2][2] = 10;
	cout << "该线性方程的系数增广矩阵为：" << endl;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 4; j++)
			cout << p[i][j] << "   ";
		cout << endl;
	}
	x[0] = x[1] = x[2] = 1;
	cout << "初值：" << endl;
	cout << "x1=" << x[0] << "   x2=" << x[1] << "   x3=" << x[2] << endl;
	
	cout << "请输入误差限：";
	cin >> eps;
	q = 1 + eps;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 4; j++)
			z[i][j] = p[i][j];
	for (i = 0; i < 3; i++)
		for (j = 0; j < 4; j++)
			if (i!=j)
				z[i][j] = z[i][j] / z[i][i];
	for (i = 0; i < 3;i++)
		for (j = 0; j < 3;j++)
			if (i == j)
				z[i][j] = 1;
			else
				z[i][j] = -z[i][j];
	for (i = 0; i < 3; i++)
	{	
		for (j = 0; j<4; j++)
			cout <<setw(5)<< z[i][j] << "    ";
		cout << endl;
	}
	while (q>eps)
	{
		a = x[0];
		b = x[1];
		c = x[2];
		for (i = 0; i < 3; i++)
		{
			s = 0;
			for (j = 0; j < 3; j++)
				if (i!=j)
					s = s + z[i][j] * x[j];
			x[i] = s+z[i][3];
		}
		q = pow(pow(x[0] - a, 2) + pow(x[1] - b, 2) + pow(x[2] - c, 2), 0.5);
	}
	cout << "线性方程组的解:" << endl;
	cout << "x1=" << x[0] << "   x2=" << x[1] << "   x3=" << x[2]<<endl;
	return 0;
}
