#include <iostream>
using namespace std;
#define MAX_VERTEX_NUM 50	// 顶点最大数量
// 边（弧）
typedef struct ArcNode {
	int vertex;				// 顶点域。即顶点在顶点数组中的下标。
	struct ArcNode* next;	// 链域。指向下一条边（弧）。
}ArcNode;
// 顶点
typedef struct VexNode {
	const char* name;			// 使用一个字符串作为顶点保存的值。也可以认为是顶点的名字。
	struct ArcNode* firstarc;	// 指向第一条边（弧）
}VexNode;
// 图。使用邻接表存储
typedef struct Graph {
	VexNode vexlist[MAX_VERTEX_NUM];	// 顶点表
	int length;							// 顶点数量
}Graph;
void DepthFirstSearch(Graph* pGraph);
void InitGraph(Graph* pGraph);
void DeleteGraph(Graph* pGraph);
int visited[MAX_VERTEX_NUM];	// 顶点访问标志数组。0 表示未访问；大于 0 表示访问的先后次序

int main()
{
	Graph graph;
	int i;
	InitGraph(&graph); // 初始化图
	// 初始化访问标志数组
	for(i=0; i<graph.length; i++)
		visited[i] = 0;
	// 深度优先搜索
	DepthFirstSearch(&graph);
	// 销毁图
	DeleteGraph(&graph);
	return 0;
}
void DepthFirstSearch(Graph* pGraph)
{
	ArcNode* Stack[MAX_VERTEX_NUM];		// 栈。
	int top = 0;						// 栈顶。0 表示空栈

	ArcNode* pArcNode;		// 边（弧）节点指针
	int i, v;				// 顶点序号。从 0 开始计数。
	int nVisitCount = 0;	// 访问计数器。
	int f[7] = {0,0,0,0,0,0,0};
	for(i=0;nVisitCount<7;){
		v=i;
		if(f[v] == 0) {
			printf("%s ",pGraph->vexlist[v].name);
			f[v]=1;	
		}else{
			pArcNode = Stack[top-1];
			top--;
			printf("%s ",pGraph->vexlist[pArcNode->vertex].name);
			pArcNode = pGraph->vexlist[pArcNode->vertex].firstarc;
			while (pArcNode != NULL) {
				if (f[pArcNode->vertex] != 1) {
					Stack[top] = pArcNode;
					top++;v++;
					f[pArcNode->vertex] = 1;					
				}
				pArcNode = pArcNode->next;
			}
		}			
		pArcNode = pGraph->vexlist[v].firstarc;
		while(pArcNode != NULL && f[pArcNode->vertex]!=1){
			Stack[top] = pArcNode;
			top++;v++;
			f[pArcNode->vertex]=1;
			pArcNode = pArcNode->next;
		}
		nVisitCount++;
	}	
	return;
}
typedef struct VertexArrayEntry {
	const char* name;				// 顶点名称。NULL 表示顶点序列结束。
	int VexIndex[MAX_VERTEX_NUM];	// 与该顶点邻接的顶点序列。-1 表示序列结束。
}VertexArrayEntry;
const VertexArrayEntry VertexArray[] = {
	{ "V0", {2, 5, 6, 1, -1}	},
	{ "V1", {3, 6, 4, 2, -1}	},
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
功能：销毁图
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
