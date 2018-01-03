
// 房屋喷漆Dlg.h : 头文件
//

#pragma once


// C房屋喷漆Dlg 对话框
class C房屋喷漆Dlg : public CDialogEx
{
// 构造
public:
	C房屋喷漆Dlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_MY_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()

	CDC *picDC;	//图片控件DC（设备描述表）
	int picClient_W;	//宽度：图片控件客户区（绘制图形、文字、图像）
	int picClient_H;	//高度：图片控件客户区
	int CenterX,CenterY;	//图片控件的客户区中心
	Graphics *graph;	//图形对象指针
	Color BkColor;	//图片控件客户区背景色（用于清除客户区内容)
	Bitmap *src;	//原图像，Bitmap类用于位图
	Bitmap *dest;	//结果图像（显示）
	int w,h;	//结果图像dest的宽和高
	BYTE *R_val;	//图像红色分量值指针（new分配）
	BYTE *G_val;	//图像绿色分量值指针，
	BYTE *B_val;	//图像蓝色分量值指针	

public:
	afx_msg void OnBnClickedButton1();
	afx_msg void On32777();
	afx_msg void On32778();
	afx_msg void Onstep1();
	afx_msg void On32772();
	afx_msg void On32773();
	afx_msg void On32774();
	afx_msg void On32775();
	afx_msg void On32776();
	afx_msg void On32780();
};

#include "ProblemDlg.h"
