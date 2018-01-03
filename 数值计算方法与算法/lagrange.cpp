#include <iostream>
using namespace std;
int main ()
{
	double function(double);   //求解函数 
	long double x,f;
	cout<<"请输入待计算x值："<<endl;
	cin>>x;
	if(x<0.10||x>1.05) 
	{
		cout<<"输入值无效，请重新输入："<<endl;
		cin>>x;
	}
	else
		f=function(x);
	cout<<"L(x)="<<f<<endl;
	return 0;
} 
double function(double x)
{
	struct point
	{
		double x0;
		double y0;
	};
	point poi[20]={{0.10,0.1103329},{0.15,0.1736223},{0.20,0.2426552},{0.25,0.3176729},
				{0.30,0.3989105},{0.35,0.4865951},{0.40,0.5809439},{0.45,0.6821617},
				{0.50,0.7904390},{0.55,0.9059492},{0.60,1.0288456},{0.65,1.1592592},
				{0.70,1.2972951},{0.75,1.4430292},{0.80,1.5965053},{0.85,1.7577308},
				{0.90,1.9266733},{0.95,2.1032563},{1.00,2.2873552},{1.05,2.4787929}};
	//初始化			
				
	int i,j,k,a,b;
	long double f=0,f1=1;
	if(x<=0.25) 
	{
		i=0; j=7;
	}
	else if(x>=0.90)
	{
		i=11; j=18;
	}
	else 
	{
		i=0; j=19;
		int m;
		while(1)           //二分法求x的位置 
		{
			m=(i+j)/2;
			if(poi[m].x0-x>0)
			{
				j=m;
			}
			else i=m+1;
			if(i>=j)
			{
				m=(i+j)/2;
				break;
			}
		}
		i=m-3;
		j=m+4; 
	}
	a=i; b=j;
	for(i=a;i<=b;i++)
	{
		for(j=a;j<=b;j++)
		{
			if(j==i) continue;
			else
			{
				f1*=(x-poi[j].x0)/(poi[i].x0-poi[j].x0);
			}
		}
		f+=(poi[i].y0)*f1;
		f1=1;
	}
	return f;
}
