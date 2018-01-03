#include <iostream>
#include <stdio.h>
using namespace std;
#define MAX_QUEUE_LENGTH 64 	// 队列最大长度
#define MAX_VERTEX_NUM 50	// 顶点最大数量
// 环形队列
typedef struct Queue{
	int buffer[MAX_QUEUE_LENGTH];	// 队列缓冲区
	int begin;		// 开始位置
	int end;		// 结束位置
	int length;		// 队列长度
}Queue;
void InitQueue(Queue* pQ);
int EnQueue(Queue* pQ, int Elem);
int DeQueue(Queue* pQ);
int QueueEmpty(Queue* pQ);
// 边（弧）
typedef struct ArcNode {
	int vertex;				// 顶点域。即顶点在顶点数组中的下标。
	struct ArcNode* next;	// 链域。指向下一条边（弧）。
}ArcNode;

// 顶点
typedef struct VexNode {
	const char* name;			// 使用一个字符串作为顶点保存的值。也可以认为是顶点的名字。
	struct ArcNode* firstarc;	// 指向第一条边（弧）。
}VexNode;

// 图。使用邻接表存储
typedef struct Graph {
	VexNode vexlist[MAX_VERTEX_NUM];	// 顶点表
	int length;							// 顶点数量
}Graph;

void BreadthFirstSearch(Graph* pGraph);
void InitGraph(Graph* pGraph);
void DeleteGraph(Graph* pGraph);
int visited[MAX_VERTEX_NUM];	// 顶点访问标志数组。0 表示未访问；大于 0 表示访问的先后次序。
Queue queue;	// 队列

int main()
{
	Graph graph;
	int i;
	// 初始化图
	InitGraph(&graph);
	// 初始化队列
	InitQueue(&queue);
	// 初始化访问标志数组
	for(i=0; i<graph.length; i++)
		visited[i] = 0;
	// 广度优先搜索
	BreadthFirstSearch(&graph);
	// 销毁图
	DeleteGraph(&graph);
	return 0;
}

/*
功能：
	广度优先搜索。

参数：
	pGraph -- 图指针
*/
void BreadthFirstSearch(Graph* pGraph)
{
	ArcNode* pArcNode;		// 边（弧）节点指针
	int i, v;				// 顶点序号。从 0 开始计数。
	int nVisitCount = 0;	// 访问计数器。
	ArcNode* Stack[MAX_VERTEX_NUM];
	int top = 0;
	v = 0;
	for (i = 0;nVisitCount < 7;) {
		if (visited[i] == 0) {
			printf("%s ", pGraph->vexlist[i].name);
			visited[i] = 1;
			nVisitCount++;
			pArcNode = pGraph->vexlist[i].firstarc;
			while (pArcNode != NULL) {
				if (visited[pArcNode->vertex] == 0) {
					Stack[top] = pArcNode;
					visited[pArcNode->vertex] = 1;
					top++;
					pArcNode = pArcNode->next;
				}
			}
		}else{
			pArcNode = Stack[v];
			v++;
			printf("%s ", pGraph->vexlist[pArcNode->vertex].name);
			visited[pArcNode->vertex] = 1;
			nVisitCount++;
			pArcNode = pGraph->vexlist[pArcNode->vertex].firstarc;
			while (pArcNode != NULL) {
				if (visited[pArcNode->vertex] == 0) {
					Stack[top] = pArcNode;
					visited[pArcNode->vertex] = 1;
					top++;					
				}
				pArcNode = pArcNode->next;
			}
		}
	}
	return;
}

/*
功能：
	使用给定的数据初始化图的邻接表

参数：
	pGraph -- 图指针
*/
typedef struct VertexArrayEntry {
	const char* name;				// 顶点名称。NULL 表示顶点序列结束。
	int VexIndex[MAX_VERTEX_NUM];	// 与该顶点邻接的顶点序列。-1 表示序列结束。
}VertexArrayEntry;
const VertexArrayEntry VertexArray[] = {
	{ "V0", {2, 5, 1, -1}		},
	{ "V1", {3, 6, 2, -1}		},
	{ "V2", {0, 3, 6, 1, -1}	},
	{ "V3", {1, 2, 4, -1}		},
	{ "V4", {1, 3, -1}			},
	{ "V5", {0, 6, -1}			},
	{ "V6", {2, 5, 0, -1}		},
	{ NULL 						}	// 结束标志
};
void InitGraph(Graph* pGraph)
{
	int i, j;
	ArcNode** pArcNodePtr;	// 指向指针的指针	
	// 重置图中的数据
	pGraph->length = 0;
	// 使用给定的数据初始化图的邻接表
	for(i=0; i<MAX_VERTEX_NUM ;i++)
	{
		if(NULL == VertexArray[i].name)
			break;		
		pGraph->vexlist[i].name = VertexArray[i].name;
		pArcNodePtr = &pGraph->vexlist[i].firstarc;
		for(j=0; j<MAX_VERTEX_NUM; j++)
		{
			if(-1 == VertexArray[i].VexIndex[j])
			{
				*pArcNodePtr = NULL;
				break;
			}		
			*pArcNodePtr = (ArcNode*)malloc(sizeof(ArcNode));
			(*pArcNodePtr)->vertex = VertexArray[i].VexIndex[j];
			pArcNodePtr = &(*pArcNodePtr)->next;
		}	
		pGraph->length++;
	}
}
/*
功能：
	销毁图

参数：
	pGraph -- 图指针
*/
void DeleteGraph(Graph* pGraph)
{
	int i;
	ArcNode* pArcNode;
	
	for(i=0; i<pGraph->length; i++)
	{
		while(pGraph->vexlist[i].firstarc != NULL)
		{
			pArcNode = pGraph->vexlist[i].firstarc;
			pGraph->vexlist[i].firstarc = pArcNode->next;
			free(pArcNode);
		}
	}
	pGraph->length = 0;
}

/*
功能：
	初始化队列。

参数：
	pQ -- 队列指针
*/
void InitQueue(Queue* pQ)
{
	pQ->begin = pQ->end = 0;
	pQ->length = 0;
}
/*
功能：
	将元素插入队尾。

参数：
	pQ -- 队列指针
	Elem -- 入队的元素
  
返回值：
	如果插入成功返回入队元素的值。
	如果插入失败返回 -1。
*/
int EnQueue(Queue* pQ, int Elem)
{
	// 队列满，入队失败。
	if(MAX_QUEUE_LENGTH == pQ->length)
		return -1;
		
	pQ->buffer[pQ->end] = Elem;
	pQ->end = (pQ->end + 1) % MAX_QUEUE_LENGTH;
	pQ->length++;
	
	return Elem;
}

/*
功能：
	将队首元素出队

参数：
	pQ -- 队列指针
  
返回值：
	如果出队成功返回出队元素的值。
	如果出队失败返回 -1。
*/
int DeQueue(Queue* pQ)
{
	int Elem;
	// 队列空，出队失败
	if(QueueEmpty(pQ))
		return -1;
	Elem = pQ->buffer[pQ->begin];
	pQ->begin = (pQ->begin + 1) % MAX_QUEUE_LENGTH;
	pQ->length--;
	return Elem;
}

/*
功能：
	判断队列是否为空。

参数：
	pQ -- 队列指针
  
返回值：
	如果队列空返回 1（真）
	如果队列非空返回 0（假）
*/
int QueueEmpty(Queue* pQ)
{
	return 0 == pQ->length ? 1 : 0;
}


