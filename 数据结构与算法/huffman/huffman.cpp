#include <iostream>
#include <algorithm>
using namespace std;

struct Node
{
	int data;
	Node *lChild,*rChild;
};

class HuffmanTree
{
public:
	HuffmanTree(int arr[],int count);	//构造函数,count表示结点个数,arr数组元素即为结点元素值
	void display();	//先序遍历哈夫曼树
private:
	Node *base;	//树根
	int nodeCount;	//结点个数
};

bool cmp(Node a,Node b)
{
	return a.data>b.data;
}

HuffmanTree::HuffmanTree(int arr[],int count)
{//构造函数,count表示结点个数,arr数组元素即为结点元素值
	nodeCount=count;
	Node *p=new Node;
	if(p==NULL)
	{
		cout<<"内存不足";
		return ;
	}
	delete p;
	sort(arr,arr+count);	//对数组进行排升序
	Node tree[100];
	int len=count;	//数组tree的长度
	for(int i=0;i<len;i++)
		(tree+i)->data=arr[i];
	while(len>1)
	{
		sort(tree,tree+len,cmp);
		Node *p=new Node;
		p->lChild=(tree+len-2);
		p->rChild=(tree+len-1);
		p->data=(tree+len-2)->data+(tree+len-1)->data;

		cout<<p->data<<"     ";	//输出元素值，验证哈夫曼编码
		cout<<"lChild:"<<p->lChild->data<<"  ";
		cout<<"rChild:"<<p->rChild->data<<endl<<endl;

		(tree+len-2)->data=p->data;
		(tree+len-2)->lChild=p->lChild;
		(tree+len-2)->rChild=p->rChild;
		len--;	
	}
	base=&tree[0];
}

//void HuffmanTree::display()
//{	//先序遍历
//	Node *p=base;
//	Node *stack[100]; //创建栈
//	int top=-1;	//栈顶
//	while(p!=NULL || top>=0)
//	{
//		while(p!=NULL)
//		{
//			cout<<(p->data)<<endl;
//			//cout<<"**"<<endl;
//			if(p->rChild)
//				stack[++top]=p->rChild;
//			p=p->lChild;
//		}
//		if(top>=0) p=stack[top--];
//	}
//}

int main ()
{
	int arr[10];
	for(int i=0;i<10;i++)
		arr[i]=i;
	HuffmanTree myHuffman(arr,10);
	//myHuffman.display();
	return 0;
}