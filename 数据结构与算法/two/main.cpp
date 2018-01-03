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
	list.remove(location+1);	//删除表中location+1位置的元素
	list.display();
	//list.remove(data-5);	//删除表中元素值为data-5的元素
	//list.display();
	return 0;
}