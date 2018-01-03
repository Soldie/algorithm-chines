
// 房屋喷漆Dlg.cpp : 实现文件
//

#include "stdafx.h"
#include "房屋喷漆.h"
#include "房屋喷漆Dlg.h"
#include "afxdialogex.h"

//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// C房屋喷漆Dlg 对话框




C房屋喷漆Dlg::C房屋喷漆Dlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(C房屋喷漆Dlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void C房屋喷漆Dlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(C房屋喷漆Dlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON1, &C房屋喷漆Dlg::OnBnClickedButton1)
	ON_COMMAND(ID_32777, &C房屋喷漆Dlg::On32777)
	ON_COMMAND(ID_32778, &C房屋喷漆Dlg::On32778)
	ON_COMMAND(ID_step1, &C房屋喷漆Dlg::Onstep1)
	ON_COMMAND(ID_32772, &C房屋喷漆Dlg::On32772)
	ON_COMMAND(ID_32773, &C房屋喷漆Dlg::On32773)
	ON_COMMAND(ID_32774, &C房屋喷漆Dlg::On32774)
	ON_COMMAND(ID_32775, &C房屋喷漆Dlg::On32775)
	ON_COMMAND(ID_32776, &C房屋喷漆Dlg::On32776)
	ON_COMMAND(ID_32780, &C房屋喷漆Dlg::On32780)
END_MESSAGE_MAP()


// C房屋喷漆Dlg 消息处理程序

BOOL C房屋喷漆Dlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: 在此添加额外的初始化代码
	picDC=GetDlgItem(IDC_Pic)->GetDC();	//获得图片控件的DC
	CRect rect;
	GetDlgItem(IDC_Pic)->GetClientRect(&rect);	//获得图片控件客户区大小
	picClient_W=rect.Width();
	picClient_H=rect.Height();
	CenterX=rect.CenterPoint().x;
	CenterY=rect.CenterPoint().y;
	//graph=new Graphics(picDC->m_hDC);	//创建图形对象，大家公用
	graph=Graphics::FromHDC(picDC->m_hDC);
	int ColorVal=GetSysColor(COLOR_BTNFACE);	//颜色值int 
	char red=GetRValue(ColorVal);	//获得颜色的R分量值
	char green=GetGValue(ColorVal);
	BYTE blue=GetBValue(ColorVal);
	int alpha=255;	//不透明，0全透明
	BkColor=Color(red,green,blue);		//颜色对象，透明度可略



	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void C房屋喷漆Dlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void C房屋喷漆Dlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR C房屋喷漆Dlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void C房屋喷漆Dlg::OnBnClickedButton1()
{
	// TODO: 在此添加控件通知处理程序代码
	exit(0);
}


void C房屋喷漆Dlg::On32777()
{
	// TODO: 在此添加命令处理程序代码
	ProblemDlg problemDlg;
	problemDlg.DoModal();
}


void C房屋喷漆Dlg::On32778()
{
	// TODO: 在此添加命令处理程序代码
	//绘制外框
	graph->Clear(BkColor);
	Pen myPen(Color::Black,1);
	Rect rect(120,10,400,400);
	graph->DrawRectangle(&myPen,rect);
	Sleep(100);

	//绘制格子线
	/*Point point[]={
	};*/
	for(int i=120;i<=520;i+=50)
	{
		graph->DrawLine(&myPen,Point(i,10),Point(i,410));
		Sleep(100);
	}
	for(int i=10;i<=410;i+=50)
	{
		graph->DrawLine(&myPen,Point(120,i),Point(520,i));
		Sleep(100);
	}
}


void C房屋喷漆Dlg::Onstep1()
{
	// TODO: 在此添加命令处理程序代码
	int alpha=255;
	Pen myPen=(Color::Red);	//笔对象
	////竖
	//for(int i=85;i<=385;i+=50)
	//	graph->DrawLine(&myPen,Point(60+i,i),Point(60+i,i+50));
	////横
	//for(int i=95;i<=445;i+=50)
	//	graph->DrawLine(&myPen,Point(i,i-10),Point(i+50,i-10));

	for(int k=95;k<=445;k+=50)
	{
		graph->DrawLine(&myPen,Point(k,k-10),Point(k+50,k-10));
		Sleep(100);
		graph->DrawLine(&myPen,Point(50+k,k-10),Point(50+k,k+40));
		Sleep(100);
	}
}


void C房屋喷漆Dlg::On32772()
{
	// TODO: 在此添加命令处理程序代码
	int alpha=255;
	Pen myPen=(Color::Blue);	//笔对象
	for(int k=95;k<=395;k+=50)
	{
		//横
		graph->DrawLine(&myPen,Point(k,k+40),Point(k+50,k+40));
		Sleep(100);
		//竖
		graph->DrawLine(&myPen,Point(50+k,k+40),Point(50+k,k+90));
		Sleep(100);
	}
}


void C房屋喷漆Dlg::On32773()
{
	// TODO: 在此添加命令处理程序代码
	Pen myPen=(Color::Green);	//笔对象
	graph->DrawLine(&myPen,Point(95,280),Point(145,280));
	Sleep(500);
	graph->DrawLine(&myPen,Point(95,290),Point(145,280));
}


void C房屋喷漆Dlg::On32774()
{
	// TODO: 在此添加命令处理程序代码
	Pen myPen=(Color::Green);	//笔对象
	graph->DrawLine(&myPen,Point(140,380),Point(140,430));
	Sleep(500);
	graph->DrawLine(&myPen,Point(140,380),Point(150,430));
}


void C房屋喷漆Dlg::On32775()
{
	// TODO: 在此添加命令处理程序代码
	Pen myPen=(Color::Green);	//笔对象
	graph->DrawLine(&myPen,Point(190,420),Point(190,370));
	Sleep(500);
	graph->DrawLine(&myPen,Point(190,370),Point(190,320));
	Sleep(500);
	graph->DrawLine(&myPen,Point(190,320),Point(195,370));
	Sleep(500);
	graph->DrawLine(&myPen,Point(195,370),Point(200,420));
}


void C房屋喷漆Dlg::On32776()
{
	// TODO: 在此添加命令处理程序代码
	Pen myPen=(Color::Green);	//笔对象
	graph->DrawLine(&myPen,Point(240,380),Point(240,430));
	Sleep(500);
	graph->DrawLine(&myPen,Point(240,380),Point(250,430));
}


void C房屋喷漆Dlg::On32780()
{
	// TODO: 在此添加命令处理程序代码
	//第一步
	int alpha=255;
	Pen myPen1=(Color::Red);	//笔对象
	for(int k=95;k<=445;k+=50)
	{
		graph->DrawLine(&myPen1,Point(k,k-10),Point(k+50,k-10));
		Sleep(100);
		graph->DrawLine(&myPen1,Point(50+k,k-10),Point(50+k,k+40));
		Sleep(100);
	}
	Sleep(1000);
	//第二步
	Pen myPen2=(Color::Blue);	//笔对象
	for(int k=95;k<=395;k+=50)
	{
		//横
		graph->DrawLine(&myPen2,Point(k,k+40),Point(k+50,k+40));
		Sleep(200);
		//竖
		graph->DrawLine(&myPen2,Point(50+k,k+40),Point(50+k,k+90));
		Sleep(200);
	}
	Sleep(1000);
	//第三步
	Pen myPen=(Color::Green);	//笔对象
	graph->DrawLine(&myPen,Point(95,280),Point(145,280));
	Sleep(500);
	graph->DrawLine(&myPen,Point(95,290),Point(145,280));
	Sleep(1000);
	//第四步
	graph->DrawLine(&myPen,Point(140,380),Point(140,430));
	Sleep(500);
	graph->DrawLine(&myPen,Point(140,380),Point(150,430));
	Sleep(1000);
	//第五步
	graph->DrawLine(&myPen,Point(190,420),Point(190,370));
	Sleep(500);
	graph->DrawLine(&myPen,Point(190,370),Point(190,320));
	Sleep(500);
	graph->DrawLine(&myPen,Point(190,320),Point(195,370));
	Sleep(500);
	graph->DrawLine(&myPen,Point(195,370),Point(200,420));
	Sleep(1000);
	//第六步
	graph->DrawLine(&myPen,Point(240,380),Point(240,430));
	Sleep(500);
	graph->DrawLine(&myPen,Point(240,380),Point(250,430));
}
