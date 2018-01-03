#include "seqlistClass.h"
#include <iostream>
using namespace std;

//无参构造函数，建立一个空的顺序表
template <class DataType> 
SeqList<DataType>::SeqList()
{
	length=0;
} 

//有参构造函数，建立一个长度为 n 的顺序表
template <class DataType>
SeqList<DataType>::SeqList(DataType arr[],int n) 
{
	if(n<1 || n>MaxSize)
		throw "参数非法";
	for(int i=0;i<n;i++)
		data[i]=arr[i];
	length=n;
}

//析构函数
template <class DataType>
SeqList<DataType>::~SeqList(){ }

//获取线性表的长度
template <class DataType>
int SeqList<DataType>::Length()   
{
	return length;
}

//获取元素值为 x 的位置
template <class DataType>
int SeqList<DataType>::Location(DataType x) 
{
	for(int i=0;i<length;i++)
		if(data[i]==x) return i+1;
	return -1;
}

//获取 i 位置的元素值
template <class DataType>
DataType SeqList<DataType>::Get(int i)
{
	if(i<0 || i>length) throw "数据非法";
	else return data[i];
}

//插入操作，在位置 i 处插入值为 x 的元素
template <class DataType>
void SeqList<DataType>::Insert(int i,DataType x)
{
	if(length==MaxSize) throw "上溢";
	if(i<0 || i>length+1) throw "数据非法";
	for(int j=length;j>=i;j--)
		data[j]=data[j-1];
	data[i-1]=x;
	length++;
}

//删除操作，删除位置 i 处的元素并返回位置 i 的元素值
template <class DataType>
DataType SeqList<DataType>::Delete(int i)
{
	if(length==0) throw "下溢";
	if(i<0 || i>length) throw "数据非法";
	DataType x;
	x=data[i-1];
	for(int j=i-1;j<length;j++)
		data[j]=data[j+1];
	length--;
	return x;
}

//输出顺序表全部元素
template <class DataType>
void SeqList<DataType>::PrintList()
{
	if(length==0)
	{
		cout<<"线性表为空、、";
		return ;
	}
	for(int i=0;i<length;i++)
	{
		cout<<data[i]<<" ";
	}
}

//输出元素表中位置 m 到位置 n 的元素(要求 m <= n)
template <class DataType>
void SeqList<DataType>::PrintList(int m,int n)
{
	if((m<1 || m>length) && (n<1 || n>length) && m>n)
		throw "数据非法";
	for(int i=m-1;i<n;i++)
		cout<<data[i]<<" "; 
}

//更改位置 i 的元素值为 x
template <class DataType>
void SeqList<DataType>::Change(int i,DataType x)
{
	if(i<1 || i>length)
		throw "数据非法";
	data[i-1]=x;
}