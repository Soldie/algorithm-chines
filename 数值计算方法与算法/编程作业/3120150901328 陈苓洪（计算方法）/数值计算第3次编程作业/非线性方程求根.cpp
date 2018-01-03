#include <iostream>
#include <iomanip>
#define EPS 0.5e-6
using namespace std;
int main()
{
	double erfen(double (*fun)(double x0),double a,double b);	//二分法
	double fun(double x);   //原函数
	double fun_1(double x); //原函数的一阶导函数
	double fun_2(double x); //原函数的二阶导函数
	double newton(double (*f)(double),double (*f_1)(double),double (*f_2)(double x),double m,double n);  //牛顿迭代法函数
	double xianjiefa(double (*fun)(double x),double m,double n);	//弦截法
	double biliqiugen(double (*fun)(double x),double m,double n);	//比例求根

	double m,n,a,b;
	cout<<"请输入求根区间："<<endl; 
	cin>>m>>n;
	cout<<endl;

	a=fun(m)*fun(n);
	if(a>0)
	{
		cout<<"输入区间无解!";
		return 0;
	}
	a=newton(fun,fun_1,fun_2,m,n);
	cout<<"Newton迭代法求得零点： ";
	cout<<setiosflags(ios::fixed)<<setprecision(7)<<a<<endl<<endl;

	a=erfen(fun,m,n);
	cout<<"二分法求解：";
	cout<<setiosflags(ios::fixed)<<setprecision(7)<<a<<endl<<endl;

	a=xianjiefa(fun,m,n);
	cout<<"弦截法求解：";
	cout<<setiosflags(ios::fixed)<<setprecision(7)<<a<<endl<<endl;

	a=biliqiugen(fun,m,n);
	cout<<"比例求根：";
	cout<<setiosflags(ios::fixed)<<setprecision(7)<<a<<endl<<endl;


	return 0; 
}

double erfen(double (*fun)(double x0),double a,double b)
{
	double x;
	while(b-a>EPS || a-b>EPS)
	{
		x=(a+b)/2;
		if(fun(x)>0)
			b=x;
		else a=x;
	}
	return x;
}

double newton(double (*f)(double x),double (*f_1)(double),double (*f_2)(double x),double m,double n)
{
	double a,a_2,i,x0,x1;
	a=f(m);
	a_2=f_2(m);
	if(a*a_2>=0) x1=x0=m;
	else  x1=x0=n;
	do{
		x0=x1;
		x1=x0-f(x0)/f_1(x0);
	}while(fabs(x1-x0)>EPS);
	return x1; 
}

double xianjiefa(double (*fun)(double x),double m,double n)
{
	double temp;
	while(m-n>EPS || n-m>EPS)
	{
		temp=n-(fun(n)*(n-m))/(fun(n)-fun(m));
		if(fun(temp)*fun(m)>0)
			m=temp;
		else n=temp;
	}
	return temp;
}

double biliqiugen(double (*fun)(double x),double m,double n)
{
	double temp=1.0;
	while(1)
	{
		temp=m-fun(m)*(m-n)/(fun(m)-fun(n));
		if(fun(temp)*fun(m)>0)
			m=temp;
		else n=temp;
		if(sqrt(fun(temp))<EPS)
			break;
	}
	return temp;
}

double fun(double x)
{
	return x*x*x-7.7*x*x+19.2*x-15.3;
}
double fun_1(double x)
{
	return 3*x*x-7.7*2*x+19.2;
}
double fun_2(double x)
{
	return 3*2*x-7.7*2;
}