const int MaxSize=100;
template <class DataType>
class SeqList
{
public:
	SeqList();      //无参构造函数，建立一个空的顺序表
	SeqList(DataType arr[],int n);   //有参构造函数，建立一个长度为 n 的顺序表
	~SeqList();
	int Length();   //获取线性表的长度
	int Location(DataType x);   //获取元素值为 x 的位置
	DataType Get(int i);  //获取 i 位置的元素值
	void Insert(int i, DataType x);  //插入操作，在位置 i 处插入值为 x 的元素
	DataType Delete(int i);   //删除操作，删除位置 i 处的元素并返回位置 i 的元素值
	void PrintList();  //输出顺序表全部元素
	void PrintList(int m,int n);   //输出元素表中位置 m 到位置 n 的元素 (m <= n)
	void Change(int i,DataType x);  //更改位置 i 的元素值为 x

private:
	DataType data[MaxSize];
	int length;
};