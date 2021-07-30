
// mosaicGUITestView.h : CmosaicGUITestView ��Ľӿ�
//

#pragma once


class CmosaicGUITestView : public CScrollView
{
protected: // �������л�����
	CmosaicGUITestView();
	DECLARE_DYNCREATE(CmosaicGUITestView)

// ����
public:
	CmosaicGUITestDoc* GetDocument() const;

// ����
public:

// ��д
public:
	virtual void OnDraw(CDC* pDC);  // ��д�Ի��Ƹ���ͼ
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual void OnInitialUpdate(); // ������һ�ε���
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

// ʵ��
public:
	virtual ~CmosaicGUITestView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// ���ɵ���Ϣӳ�亯��
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnCalibTool();
};

#ifndef _DEBUG  // mosaicGUITestView.cpp �еĵ��԰汾
inline CmosaicGUITestDoc* CmosaicGUITestView::GetDocument() const
   { return reinterpret_cast<CmosaicGUITestDoc*>(m_pDocument); }
#endif

