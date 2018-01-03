#include<iostream>
#include<cmath>
#define f(x) (1.0/(1+x))    //积分函数
#define epsilon 0.0001    
#define MAXREPT 10          //迭代次数
using namespace std;

int main()
{
	double Romberg(double , double );
    double a,b;
    cout<<"Romberg积分,请输入积分范围:"<<endl;
    cin>>a>>b;
    cout<<"积分结果:"<<Romberg(a, b)<<endl;
    system("pause");
    return 0;
}

double Romberg(double a, double b)
{ 
    int m, n;
    double h, x, s, q,ep,p; 
    double *y = new double[MAXREPT];
    h = b - a;
    y[0] = h*(f(a) + f(b))/2.0;
    m = 1;
    n = 1;
    ep = epsilon + 1.0;
    while ((ep >= epsilon) && (m < MAXREPT))
    {
        p = 0.0;
        for (int i=0; i<n; i++)
        {
            x = a + (i+0.5)*h;
            p = p + f(x);
        }
        p = (y[0] + h*p)/2.0;      
        s = 1.0;
        for (int k=1; k<=m; k++)
        {
            s = 4.0*s;
            q = (s*p - y[k-1])/(s - 1.0);
            y[k-1] = p;
            p = q;
        }
        p = fabs(q - y[m-1]);
        m = m + 1;
        y[m-1] = q;
        n = n + n; h = h/2.0;
    }
    return (q);
}
