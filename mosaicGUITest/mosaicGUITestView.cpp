
// mosaicGUITestView.cpp : CmosaicGUITestView 类的实现
//

#include "stdafx.h"
#include <CameraCalib.h>
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
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
	// 标准打印命令
	ON_COMMAND(ID_FILE_PRINT, &CScrollView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CScrollView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CmosaicGUITestView::OnFilePrintPreview)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_COMMAND(ID_CALIB_TOOL, &CmosaicGUITestView::OnCalibTool)
END_MESSAGE_MAP()

// CmosaicGUITestView 构造/析构

CmosaicGUITestView::CmosaicGUITestView()
{
	// TODO:  在此处添加构造代码

}

CmosaicGUITestView::~CmosaicGUITestView()
{
}

BOOL CmosaicGUITestView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO:  在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式

	return CScrollView::PreCreateWindow(cs);
}

// CmosaicGUITestView 绘制

void CmosaicGUITestView::OnDraw(CDC* /*pDC*/)
{
	CmosaicGUITestDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO:  在此处为本机数据添加绘制代码
}

void CmosaicGUITestView::OnInitialUpdate()
{
	CScrollView::OnInitialUpdate();

	CSize sizeTotal;
	// TODO:  计算此视图的合计大小
	sizeTotal.cx = sizeTotal.cy = 100;
	SetScrollSizes(MM_TEXT, sizeTotal);
}


// CmosaicGUITestView 打印


void CmosaicGUITestView::OnFilePrintPreview()
{
#ifndef SHARED_HANDLERS
	AFXPrintPreview(this);
#endif
}

BOOL CmosaicGUITestView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// 默认准备
	return DoPreparePrinting(pInfo);
}

void CmosaicGUITestView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  添加额外的打印前进行的初始化过程
}

void CmosaicGUITestView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO:  添加打印后进行的清理过程
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


// CmosaicGUITestView 诊断

#ifdef _DEBUG
void CmosaicGUITestView::AssertValid() const
{
	CScrollView::AssertValid();
}

void CmosaicGUITestView::Dump(CDumpContext& dc) const
{
	CScrollView::Dump(dc);
}

CmosaicGUITestDoc* CmosaicGUITestView::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CmosaicGUITestDoc)));
	return (CmosaicGUITestDoc*)m_pDocument;
}
#endif //_DEBUG


// CmosaicGUITestView 消息处理程序


void CmosaicGUITestView::OnCalibTool()
{
	// TODO: 在此添加命令处理程序代码
	CameraCalib::calibmain("test");
}
