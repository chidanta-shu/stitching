
// mosaicGUITestView.cpp : CmosaicGUITestView ���ʵ��
//

#include "stdafx.h"
#include <CameraCalib.h>
// SHARED_HANDLERS ������ʵ��Ԥ��������ͼ������ɸѡ�������
// ATL ��Ŀ�н��ж��壬�����������Ŀ�����ĵ����롣
#ifndef SHARED_HANDLERS
#include "mosaicGUITest.h"
#endif

#include "mosaicGUITestDoc.h"
#include "mosaicGUITestView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CmosaicGUITestView

IMPLEMENT_DYNCREATE(CmosaicGUITestView, CScrollView)

BEGIN_MESSAGE_MAP(CmosaicGUITestView, CScrollView)
	// ��׼��ӡ����
	ON_COMMAND(ID_FILE_PRINT, &CScrollView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CScrollView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CmosaicGUITestView::OnFilePrintPreview)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_COMMAND(ID_CALIB_TOOL, &CmosaicGUITestView::OnCalibTool)
END_MESSAGE_MAP()

// CmosaicGUITestView ����/����

CmosaicGUITestView::CmosaicGUITestView()
{
	// TODO:  �ڴ˴���ӹ������

}

CmosaicGUITestView::~CmosaicGUITestView()
{
}

BOOL CmosaicGUITestView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO:  �ڴ˴�ͨ���޸�
	//  CREATESTRUCT cs ���޸Ĵ��������ʽ

	return CScrollView::PreCreateWindow(cs);
}

// CmosaicGUITestView ����

void CmosaicGUITestView::OnDraw(CDC* /*pDC*/)
{
	CmosaicGUITestDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO:  �ڴ˴�Ϊ����������ӻ��ƴ���
}

void CmosaicGUITestView::OnInitialUpdate()
{
	CScrollView::OnInitialUpdate();

	CSize sizeTotal;
	// TODO:  �������ͼ�ĺϼƴ�С
	sizeTotal.cx = sizeTotal.cy = 100;
	SetScrollSizes(MM_TEXT, sizeTotal);
}


// CmosaicGUITestView ��ӡ


void CmosaicGUITestView::OnFilePrintPreview()
{
#ifndef SHARED_HANDLERS
	AFXPrintPreview(this);
#endif
}

BOOL CmosaicGUITestView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// Ĭ��׼��
	return DoPreparePrinting(pInfo);
}

void CmosaicGUITestView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  ��Ӷ���Ĵ�ӡǰ���еĳ�ʼ������
}

void CmosaicGUITestView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  ��Ӵ�ӡ����е��������
}

void CmosaicGUITestView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CmosaicGUITestView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CmosaicGUITestView ���

#ifdef _DEBUG
void CmosaicGUITestView::AssertValid() const
{
	CScrollView::AssertValid();
}

void CmosaicGUITestView::Dump(CDumpContext& dc) const
{
	CScrollView::Dump(dc);
}

CmosaicGUITestDoc* CmosaicGUITestView::GetDocument() const // �ǵ��԰汾��������
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CmosaicGUITestDoc)));
	return (CmosaicGUITestDoc*)m_pDocument;
}
#endif //_DEBUG


// CmosaicGUITestView ��Ϣ�������


void CmosaicGUITestView::OnCalibTool()
{
	// TODO: �ڴ���������������
	CameraCalib::calibmain("test");
}
