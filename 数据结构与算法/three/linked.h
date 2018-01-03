template <class dataType>
struct Node
{
	dataType data;
	Node<dataType> *next;
};

template <class dataType>
class Linked
{
public:
	Linked();
	Linked(dataType arr[],int _length);	
	~Linked();
	bool insert(int location,dataType _data);	//在location插入_data
	bool remove(int location);	//删除location位置的元素
	//bool remove(dataType _data);	//删除元素值为_data的元素
	void display();	//输出链表
private:
	Node<dataType> *first;
};

template<class dataType>
Linked<dataType>::Linked()
{
	first=new Node;
	first->=NULL;
}

//构造函数(头插法)
template<class dataType>
Linked<dataType>::Linked(dataType arr[],int _length)
{
	first=new Node<dataType>;
	first->next=NULL;
	for(int i=0;i<_length;i++)
	{
		Node<dataType> *f;
		f=new Node<dataType>;
		f->data=arr[i];
		f->next=first->next;
		first->next=f;
	}
}

//析构函数
template<class dataType>
Linked<dataType>::~Linked()
{
	while(first!=NULL)
	{
		Node<dataType> *f;
		f=first;
		first=first->next;
		delete f;
	}
}

//在location插入_data
template <class dataType>
bool Linked<dataType>::insert(int location,dataType _data)
{
	Node<dataType> *f=new Node<dataType>;
	f=first->next;
	int m=1;
	while(m!=location && f!=NULL)
	{
		if(m==location-1)
		{
			Node<dataType> *a=new Node<dataType>;
			a->data=_data;
			a->next=f->next;
			f->next=a;
			return true;
		}
		f=f->next;
		m++;
	}
	return false;
}

//删除location位置的元素
template <class dataType>
bool Linked<dataType>::remove(int location)
{
	Node<dataType> *f;
	f=first;
	int m=0;
	while(m!=location-1 && f!=NULL)
	{
		f=f->next;
		m++;
		if(m==location-2)
		{
			Node<dataType> *a;
			a=f->next;
			f->next=(f->next)->next;
			delete a;
			return true;
		}
	}
	return false;
}

//输出链表
template<class dataType>
void Linked<dataType>::display()
{
	if(first==NULL) return;
	Node<dataType> *f;
	f=first->next;
	cout<<"链表数据为："<<endl;
	while(f!=NULL)
	{
		cout<<f->data<<" ";
		f=f->next;
	}
	cout<<endl;
}