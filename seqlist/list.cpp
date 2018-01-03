#include <iostream>
#include "seqlistClass.h"
using namespace std;

//顺序表的各种操作 
int main ()
{
	int arr[10]={1,2,3,4,5,6,7,8,9,0};
	SeqList<int> myList;
	SeqList<int> myList_1(arr,10);
	myList.PrintList();
	myList_1.PrintList();
	return 0;
}
