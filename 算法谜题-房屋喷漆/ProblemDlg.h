#pragma once


// ProblemDlg 对话框

class ProblemDlg : public CDialogEx
{
	DECLARE_DYNAMIC(ProblemDlg)

public:
	ProblemDlg(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~ProblemDlg();

// 对话框数据
	enum { IDD = IDD_Problem_Dlg };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
};

#include "房屋喷漆Dlg.h"
