#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
int main ()
{
	double * Gauss(double a[4][4],double b[4],int n=4);	//高斯列主元法
	double * Doolittle(double a[4][4],double b[4],int n);	//LU分解法
	void Guss_Jordan(int n,double *a);	//高斯约当
	double a[4][4]={{7.2,2.3,-4.4,0.5},{1.3,6.3,-3.5,2.8},{5.6,0.9,8.1,-1.3,},{1.5,0.4,3.7,5.9}};
	double b[4]={15.1,1.8,16.6,36.9};

	cout<<"高斯列主元法："<<endl;
	double *p;
	p=Gauss(a,b,4);
	for(int i=0;i<4;i++)
		cout<<setiosflags(ios::fixed)<<setprecision(7)<<setw(15)<<*(p+i);
	cout<<endl<<endl;
	cout<<"LU分解法："<<endl;
	p=Doolittle(a,b,4);
	for(int i=0;i<4;i++)
		cout<<setiosflags(ios::fixed)<<setprecision(7)<<setw(15)<<*(p+i);
	cout<<endl<<endl;

	cout<<"高斯-约当消元法："<<endl;
	double arr[2][2] = { {1,3 }, {2,4 } };
	Guss_Jordan(2, arr[0]);
	cout<<endl;

	return 0;
}

//高斯列主元法
double * Gauss(double a[4][4],double b[4],int n=4)
{
	
	for(int i=1;i<n;i++)
	{
		int k=i;
		for(int j=i+1;j<n;j++)
			if(fabs(a[k][i])<fabs(a[j][i]))
				k=j;
		for(int j=i;j<n;j++)
		{
			double t=a[i][j];
			a[i][j]=a[k][j];
			a[k][j]=t;
		}
		double t=b[i];
		b[i]=b[k];
		b[k]=t;
		for(int j=i+1;j<n;j++)
		{
			a[j][i]=a[j][i]/a[i][i];
			for(int k=i+1;k<n;k++)
			{
				a[j][k]=a[j][k]-a[j][i]*a[i][k];
			}
			b[j]=b[j]-a[j][i]*b[i];
		}
	}
	for(int i=n-1;i>=0;i--)
	{
		for(int j=i+1;j<n;j++)
		{
			b[i]=b[i]-a[i][j]*b[j];
			b[i]=b[i]/a[i][i];
		}
	}
	return b;
}

//LU分解法
double * Doolittle(double a[4][4],double b[4],int n=4)
{
	double l[4][4],u[4][4],y[4],x[4];
	for(int k=0;k<n;k++)
	{
		//计算U的第k行元素
		for(int j=k;j<n;j++)
		{
			u[k][j]=a[k][j];
			for(int r=0;r<k-1;r++)
			{
				u[k][j]-=l[k][r]*u[r][j];
			}
		}
		//计算L的第k列元素
		for(int i=k+1;i<n;i++)
		{
			l[i][k]=a[i][k];
			for(int r=0;r<k-1;r++)
			{
				l[i][k]-=l[i][r]*u[r][k];
			}
			l[i][k]/=u[k][k];
		}
	}
	//解方程组LY=b
	for(int i=0;i<n;i++)
	{
		y[i]=b[i];
		for(int j=0;j<i-1;j++)
		{
			y[i]-=l[i][j]*y[j];
		}
	}
	//解方程组UX=Y
	for(int i=n-1;i>=0;i--)
	{
		x[i]=y[i];
		if(i!=n-1)
		{
			for(int j=i+1;j<n;j++)
			{
				x[i]-=u[i][j]*x[j];
			}
			x[i]/=u[i][i];
		}
	}
	return x;
}

void Guss_Jordan(int n,double *a)  
{
	double*a_temp = new double[n*n * 2];
	double *temp = new double[n*n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i==j)
			{
				*(temp + i*n + j) = 1;
			}
			else
			{
				*(temp + i*n + j) = 0;
			}
		}
	}
	for (int i = 0; i < n; i++)//把a与单位矩阵组合
	{
		for (int j = 0; j < 2*n; j++)
		{
			if (j>n-1)
			{
				*(a_temp+i*n*2+j) = *(temp + i*n + j - n);
			}
			else
			{
				*(a_temp + i*n * 2 + j) = *(a+i *n + j);
			}
		}
	}
	delete[]temp;
	for (int i = 0; i < n; i++)
	{
		double m = *(a_temp + i*n * 2 + i);
		if (m==0)
		{
			cout << "程序失败"; return;
		}
		for (int j = 0; j < 2*n; j++)
		{
			*(a_temp + i*n * 2 + j) = *(a_temp + i*n * 2 + j) / m;
		}
		for (int k = i+1; k < n; k++)
		{
			double Multiples = *(a_temp + k * 2 * n + 0);
			for (int j = 0; j < 2 * n; j++)
			{
				*(a_temp + k * 2 * n + j) = *(a_temp + k * 2 * n + j) - Multiples*(*(a_temp + (k - 1) * 2 * n + j));
			}
		}
		for (int k = i; k >0; k--)
		{
			double Multiples =*(a_temp+(k-1)*n*2+k) /(*(a_temp + k * 2 * n + k));
			for (int j = 0; j < 2 * n; j++)
			{
				*(a_temp + (k-1) * 2 * n + j) = *(a_temp + (k-1) * 2 * n + j) - Multiples*(*(a_temp + k * 2 * n + j));
			}
		}
	}
	
	cout << "逆矩阵为：\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2 * n; j++)
		{
			if (j>n - 1)
			{
				cout << a_temp[i * 2 * n + j]<<'\t';
			}

		}
		cout << endl;
	}
	delete[]a_temp;

}