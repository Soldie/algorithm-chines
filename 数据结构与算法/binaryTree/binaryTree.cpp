#include <iostream>
using namespace std;

struct Node
{
	int data;
	Node *lChild,*rChild;
};

class BinaryTree
{
public:
	BinaryTree(int arr[],int count);	//用数组arr的元素初始化二叉树（生成的是完全二叉树）
	void display();	//	先序遍历
private:
	Node *base;		//树根
	int nodeCount;		//节点个数
};

BinaryTree::BinaryTree(int arr[],int count)
{//用数组arr的元素初始化二叉树（生成的是完全二叉树）
	Node *first=new Node;
	if(first==NULL)
	{
		cout<<"内存不足！";
		return ;
	}
	delete first;
	nodeCount=count;
	Node *tree=new Node[nodeCount];
	for(int i=1;i<=nodeCount;i++)
	{
		(tree+i-1)->data=arr[i-1];
		if(i*2<=nodeCount)
			(tree+i-1)->lChild=(tree+2*i-1);
		else (tree+i-1)->lChild=NULL;
		if(i*2+1<=nodeCount)
			(tree+i-1)->rChild=(tree+2*i);
		else (tree+i-1)->rChild=NULL;
	}
	base=tree;
}

void BinaryTree::display()
{	//先序遍历
	Node *p=base;
	Node *stack[100]; //创建栈
	int top=-1;	//栈顶
	while(p || top>=0)
	{
		while(p)
		{
			cout<<p->data<<endl;
			if(p->rChild)
				stack[++top]=p->rChild;
			p=p->lChild;
		}
		if(top>=0) p=stack[top--];
	}
}

int main()
{
	int arr[10];
	cout<<"用于初始化二叉树的数组元素："<<endl<<endl;
	for(int i=0;i<10;i++)
	{
		arr[i]=i*i;
		cout<<arr[i]<<"  ";
    }
	BinaryTree myTree(arr,10);
	cout<<endl<<"先序遍历二叉树："<<endl<<endl;
	myTree.display();
	return 0;
}
