
// mosaicGUITest.h : mosaicGUITest Ӧ�ó������ͷ�ļ�
//
#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"       // ������


// CmosaicGUITestApp:
// �йش����ʵ�֣������ mosaicGUITest.cpp
//

class CmosaicGUITestApp : public CWinAppEx
{
public:
	CmosaicGUITestApp();


// ��д
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();

// ʵ��
	BOOL  m_bHiColorIcons;

	virtual void PreLoadState();
	virtual void LoadCustomState();
	virtual void SaveCustomState();

	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CmosaicGUITestApp theApp;
