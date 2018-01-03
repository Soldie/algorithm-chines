
// 房屋喷漆.h : PROJECT_NAME 应用程序的主头文件
//

#pragma once

#ifndef __AFXWIN_H__
	#error "在包含此文件之前包含“stdafx.h”以生成 PCH 文件"
#endif

#include "resource.h"		// 主符号
#include "gdiplusinit.h"


// C房屋喷漆App:
// 有关此类的实现，请参阅 房屋喷漆.cpp
//

class C房屋喷漆App : public CWinApp
{
public:
	C房屋喷漆App();

// 重写
public:
	virtual BOOL InitInstance();

// 实现

	DECLARE_MESSAGE_MAP()
	GdiplusStartupInput GdiInput;
	ULONG_PTR GdiToken;
	int ExitInstance(void);
};

extern C房屋喷漆App theApp;