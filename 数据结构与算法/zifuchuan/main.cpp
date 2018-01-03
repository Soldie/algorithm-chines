#include <iostream>
#include <string>
#define MaxSize 100
using namespace std;

class String
{
public:
	String(string str,int _length);	//初始化
	bool insertInto(char data,int location);	//在location位置插入数据data
	void display();	//打印全部数据
private:
	char arr[MaxSize];
	int length;
};

String::String(string str,int _length)
{
	for(int i=0;i<_length;i++)
		arr[i]=str.at(i);
	length=_length;
}

bool String::insertInto(char data,int location)
{//在location位置插入数据data
	if(location<1 || location>length+1 || length==MaxSize)
		return false;
	for(int i=length;i>=location;i--)
		arr[i]=arr[i-1];
	arr[location-1]=data;
	length++;
	return true;
}

void String::display()
{
	for(int i=0;i<length;i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

int main ()
{
	string str;
	cout<<"请输入一个字符串："<<endl<<endl; 
	cin>>str;
	String test(str,str.size());
	cout<<"输出字符串顺序结构的字符串："<<endl<<endl;
	test.display();
	test.insertInto('A',5);	//在位置5插入字符‘A’
	cout<<"执行插入操作后的字符串："<<endl<<endl;
	test.display();
	return 0;
}
