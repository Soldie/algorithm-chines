#include <iostream>
#include <iomanip>
using namespace std;
struct point
{
	double x0;
	double y0;
};

int main ()
{
	void function1(double m,struct point *p,int *pi,int *pj); //二分法求i,j 
	double newton(double m,struct point *p,int i,int j);     //newton插值函数
	double m,f;
	int i,j,*pi=&i,*pj=&j; 
	point poi[20]={{0.10,0.1103329},{0.15,0.1736223},{0.20,0.2426552},{0.25,0.3176729},
				{0.30,0.3989105},{0.35,0.4865951},{0.40,0.5809439},{0.45,0.6821617},
				{0.50,0.7904390},{0.55,0.9059492},{0.60,1.0288456},{0.65,1.1592592},
				{0.70,1.2972951},{0.75,1.4430292},{0.80,1.5965053},{0.85,1.7577308},
				{0.90,1.9266733},{0.95,2.1032563},{1.00,2.2873552},{1.05,2.4787929}};
	//初始化
	/*for(int i=0;i<20;i++)
	{
		cin>>poi[i].x0>>poi.y0;
	}*/      //标准输入20组数据 
	cout<<"Please enter date:"<<endl;
	cin>>m;
	if(m<0.10 || m>1.05)
	{
		cout<<"Invalid data! Please re-enter the date:"<<endl;
		cin>>m; 
	 } 
	function1(m,poi,pi,pj);
	f=newton(m,poi,i,j);
	cout<<"f("<<m<<")= "<<setiosflags(ios::fixed)<<setprecision(7)<<f<<endl;
	return 0;
}

double newton(double m,struct point *p,int i,int j)    //newton插值函数
{
	double function2(struct point *p,int i,int j);
	int k,ai,aj,a,b;
	double Nx,Nx0=1;
	Nx=p[i].y0;
	for(a=i+1;a<=j;a++)
	{
		for(b=0;b<a-i;b++)
		{
			Nx0*=m-p[i+b].x0;
		}
		ai=i;
		aj=i+b;
		Nx+=Nx0*function2(p,ai,aj);
	}
	return Nx;
}

void function1(double m,struct point *p,int *pi,int *pj) //二分法求i,j 
{
	int i,j;
	if(m<=0.25)
	{
		i=0;
		j=7;	
	}		
	else if(m>=0.90)
	{
		i=12;
		j=19;
	}	
	else 
	{
		int l=0,r=19,x;
		while(1)
		{
			x=(l+r)/2;
			if(p[x].x0==m)
			{
				i=x-3;
				j=x+4;
				break;
			}
			else if(p[x].x0>m)
			{
				r=x-1;
			}
			else 
			{
				l=x+1;
			}
			if(l>=r)
			{
				x=(l+r)/2;
				i=x-3;
				j=x+4;
				break;
			}
		}	
	}
	*pi=i;
	*pj=j;	
}

double function2(struct point *p,int i,int j)     //差商函数 
{
	double f;
	if(i==j)
	{
		f=p[i].y0;
	}
	else
	{
		f=(function2(p,i+1,j)-function2(p,i,j-1))/(p[j].x0-p[i].x0);
	}
	return f;
}