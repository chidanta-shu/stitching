
// mosaicGUITestView.h : CmosaicGUITestView 类的接口
//

#pragma once


class CmosaicGUITestView : public CScrollView
{
protected: // 仅从序列化创建
	CmosaicGUITestView();
	DECLARE_DYNCREATE(CmosaicGUITestView)

// 特性
public:
	CmosaicGUITestDoc* GetDocument() const;

// 操作
public:

// 重写
public:
	virtual void OnDraw(CDC* pDC);  // 重写以绘制该视图
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual void OnInitialUpdate(); // 构造后第一次调用
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

// 实现
public:
	virtual ~CmosaicGUITestView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// 生成的消息映射函数
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnCalibTool();
};

#ifndef _DEBUG  // mosaicGUITestView.cpp 中的调试版本
inline CmosaicGUITestDoc* CmosaicGUITestView::GetDocument() const
   { return reinterpret_cast<CmosaicGUITestDoc*>(m_pDocument); }
#endif

