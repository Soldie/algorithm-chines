#include <iostream>
using namespace std;
int main()
{
	double a[2][2]={0,1,1,1,};
	int i, j, k = 0;
	double x[] = { 1, 1 },x1[2];
	cout << "�þ���" << endl;
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
			cout << a[i][j] << " ";
		cout << endl;
	}
	cout << "��ʼֵ��x1=" << x[0] << "   x2=" << x[1] << endl;
	double eps;
	cout << "����������ޣ�";
	cin >> eps;
	double p;
	p = eps + 1;
	while (p>eps)
	{
		x1[0] = x[0];
		x1[1] = x[1];
		for (i = 0; i < 2; i++)
		{
			k = 0;
			for (j = 0; j < 2; j++)
				k = k + a[i][j] * x1[j];
			x[i] = k;
		}
		p = fabs(x[0] / x1[0] - x[1] / x1[1]);
	}
	cout<<endl;
	double y = x[0] /x1[0];
	cout << "����ֵ:" <<y<< endl;
	cout << "����������" << x[0] << "  " << x[1] << endl;
	return 0;
}
