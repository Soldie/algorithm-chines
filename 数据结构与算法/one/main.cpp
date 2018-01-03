#include <iostream>
#include "SeqList.h"
using namespace std;
int main ()
{
	int arr[10],location=5,data=8,length=10;
	for(int i=0;i<length;i++)	//初始化数组
		arr[i]=i+1;
	SeqList<int> list(arr,10);
	list.display();	//输出线性表
	list.insert(location,data);	//在表中location位置插入data
	list.display();
	return 0;
}