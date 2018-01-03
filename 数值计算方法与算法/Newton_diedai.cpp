#include <iostream>
#include <cmath>
#include <iomanip>
#define EPS 1e-6
using namespace std;
int main ()
{
	double fun(double x);   //原函数
	double fun_1(double x); //原函数的一阶导函数
	double fun_2(double x); //原函数的二阶导函数
	double newton(double (*f)(double),double (*f_1)(double),double (*f_2)(double x),double m,double n);  //牛顿迭代法函数
	double m,n,a,b;
	cout<<"请输入求根区间："<<endl; 
	cin>>m>>n;
	a=fun(m)*fun(n);
	if(a>0)
	{
		cout<<"输入区间无解!";
		return 0;
	}
	a=newton(fun,fun_1,fun_2,m,n);
	cout<<"Newton迭代法求得零点： ";
	cout<<setiosflags(ios::fixed)<<setprecision(7)<<a<<endl;
	return 0; 
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
