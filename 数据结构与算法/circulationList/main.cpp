#include <iostream>
using namespace std;

struct Node
{
	long data;
	Node *left,*right;
};

class List
{
public:
	//~List();	//析构函数，释放内存
	List(long a[],long length);
	long getLength();	//获取链表长度
	bool insertInto(long _data,long location);	//在location插入_data
	bool remove(long location);	//删除位置location的元素
	void display();	//打印链表
private:
	Node *head;
};

//List::~List()
//{
//	while(head->right!=head)
//	{
//		Node *p=head;
//		head=head->right;
//		((p->right)->right)->left=p;
//		p->right=(p->right)->right;
//		delete p;
//	}
//}


List::List(long a[],long length)
{
	//初始化链表
	Node *p=new Node;
	if(p==NULL)
	{
		cout<<"内存不足!";
		return ;
	}
	head=p;
	p->data=a[0];
	p->left=p;
	p->right=p;
	for(long i=1;i<length;i++)
	{
		Node *node=new Node;
		node->data=a[i];
		node->right=p;
		node->left=p->left->right;
		p->left->right=node;
		p->left=node;
	}
}

long List::getLength()
{	//获取链表长度
	Node *p=head->right;
	if(p==head) return 0;
	long count=1;
	while(p!=head)
	{
		count++;
		p=p->right;
	}
	return count;
}

bool List::insertInto(long _data,long location)
{//在location插入_data
	Node *p=head->right;
	long count=0;
	while(p!=head && count<=location-1)
	{
		p=p->right;
		count++;
	}
	if(p==head) return false;
	Node *node=new Node;
	node->data=_data;
	node->right=p->right;
	node->left=(p->right)->left;
	(p->right)->left=node;
	p->right=node;
	return true;
}

bool List::remove(long location)
{//删除位置location的元素
	Node *p=head->right;
	long count=1;
	while(p!=head && count<=location-1)
	{
		p=p->right;
		count++;
	}
	if(p==head) return false;
	((p->right)->right)->left=p;
	p->right=(p->right)->right;
	return true;
}

void List::display()
{	//打印链表
	Node *p=head->right;
	while(p!=head)
	{
		cout<<p->data<<" ";
		p=p->right;
	}
	cout<<endl;
}

int main ()
{
	long a[10]={1,2,3,4,5,6,7,8,9,10};
	List list(a,10);	//定义List类的对象
	//for(long i=0;i<10;i++)
	//	list.insertInto(a[i],i+1);	//将数组元素全部插入链表中
	cout<<"初始化链表："<<endl<<endl; 
	list.display();	//打印链表，方便查看操作是否成功
	list.insertInto(11,5);	//在位置5插入11
	cout<<"插入操作："<<endl<<endl;
	list.display();	//打印链表，方便查看操作是否成功
	long location=5;
	list.remove(location);	//删除location位置的元素
	cout<<"删除操作："<<endl<<endl;
	list.display();	//打印链表，方便查看操作是否成功
	return 0;
}
