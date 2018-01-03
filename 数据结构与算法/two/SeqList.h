#define MAXSIZE 100
template <class dataType> 
class SeqList
{
public:
	SeqList(){};
	SeqList(dataType arr[],int _length);	//有参构造函数，建立一个长度为 n 的顺序表
	bool remove(int location);	//删除位置location的元素
	//bool remove(dataType _data);	//删除值为_data的元素
	void display();	//打印数据
private:
	dataType data[MAXSIZE];
	int length;
};

template <class dataType>
SeqList<dataType>::SeqList(dataType arr[],int _length)
{
	if(_length<1 || _length>MAXSIZE)
		return;
	for(int i=0;i<_length;i++)
		data[i]=arr[i];
	length=_length;
}


//删除位置location的元素
template <class dataType>
bool SeqList<dataType>::remove(int location)
{
	if(location<1 || location>length)
		return false;
	for(int i=location-1;i<length-2;i++)
		data[i]=data[i+1];
	length--;
	return true;
}

//删除值为_data的元素
//template <class dataType>
//bool SeqList<dataType>::remove(dataType _data)
//{
//	bool f=false;
//	for(int i=0;i<length;i++)
//	{
//		if(data[i]==_data)
//		{
//			f=true;
//			for(int j=i;i<length-2;i++)
//				data[j]=data[j+1];
//		}
//	}
//	return f;
//}

template <class dataType>
void SeqList<dataType>::display()
{
	cout<<"顺序表数据："<<endl;
	for(int i=0;i<length;i++)
	{
		cout<<data[i]<<" ";
	}
	cout<<endl<<"线性表长度："<<length<<endl;
}