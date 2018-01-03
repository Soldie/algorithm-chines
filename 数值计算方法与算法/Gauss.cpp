#include <iostream>
#include <iomanip>
using namespace std;
int main ()
{
	void Gauss(double *A,int m,int n);  //Gauss函数声明，m:行，n:列
	void Show(double *A,int m,int n);   //标准输出矩阵 
	double A[]={-1,3,2,4,0,1,1,2,3,2,0,5};
	int m=3,n=4;
	Gauss(A,3,4);
	Show(A,3,4);
	return 0;
 } 
 
void Gauss(double *A,int m,int n)
{
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			if(j>i){
				for(int k=1;k<n;k++)
				{
					A[j*n+k]-=A[i*n+k]*(A[j*n+i]/A[i*n+i]);
				}
				A[j*n+i]=0;
			}
		}
	}
}

void Show(double *A,int m,int n)
{
	
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			cout<<setw(3)<<A[i*n+j]<<" ";
		}
		cout<<endl;
	}
	
}
