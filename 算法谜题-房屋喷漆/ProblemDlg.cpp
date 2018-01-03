// ProblemDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "房屋喷漆.h"
#include "ProblemDlg.h"
#include "afxdialogex.h"


// ProblemDlg 对话框

IMPLEMENT_DYNAMIC(ProblemDlg, CDialogEx)

ProblemDlg::ProblemDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(ProblemDlg::IDD, pParent)
{

}

ProblemDlg::~ProblemDlg()
{
}

void ProblemDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(ProblemDlg, CDialogEx)
	ON_BN_CLICKED(IDOK, &ProblemDlg::OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, &ProblemDlg::OnBnClickedCancel)
END_MESSAGE_MAP()


// ProblemDlg 消息处理程序


void ProblemDlg::OnBnClickedOk()
{
	// TODO: 在此添加控件通知处理程序代码
	C房屋喷漆Dlg homePage;
	homePage.DoModal();
}


void ProblemDlg::OnBnClickedCancel()
{
	// TODO: 在此添加控件通知处理程序代码
	exit(0);
}
