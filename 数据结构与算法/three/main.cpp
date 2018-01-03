#include <iostream>
#include "linked.h"
using namespace std;
int main ()
{
	int arr[10],location=5,data=8,n=10;
	for(int i=0;i<10;i++)
		arr[i]=i+1;
	Linked<int> linked(arr,n);
	linked.display();	//输出链表
	cout<<endl<<"插入操作："<<endl;
	linked.insert(location,data);	//在location插入data
	linked.display();
	cout<<endl<<"删除操作："<<endl;
	linked.remove(location+2);	//删除位置location+2的元素
	linked.display();
	return 0;
}