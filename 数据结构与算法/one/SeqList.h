#define MAXSIZE 100
template <class dataType> 
class SeqList
{
public:
	SeqList(){};
	SeqList(dataType arr[],int _length);	//有参构造函数，建立一个长度为 n 的顺序表
	bool insert(int location,dataType _date);	//在location位置插入数据_data
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

//在location位置插入数据_data
template <class dataType>
bool SeqList<dataType>::insert(int location,dataType _data)	
{	
	if(MAXSIZE==length)
		return false;
	else if(location<1 || location>length)
		return false;
	for(int i=length;i>=location;i--)
		data[i]=data[i-1];
	data[location-1]=_data;
	length++;
	return true;
}

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