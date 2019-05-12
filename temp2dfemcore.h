#ifndef TEMP2DFEMCORE_H
#define TEMP2DFEMCORE_H
#include "datatype.h"
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include "armadillo"
#include "widget.h"
#include <QVector>

using namespace arma;
class CTemp2DFEMCore
{
    const double TEMP = 273.15;
public:
    CTemp2DFEMCore(Widget *parent = nullptr, const char *fn = nullptr);
    ~CTemp2DFEMCore();
    int Load2DMeshCOMSOL();   //加载COMSOL分网文件
    int preCalculation();   //计算基本几何参数
    int setCondition(Demo showWhat); //根据模型的不同设置边界条件和负载情况
    int StaticAxisAssemble(); //二维轴对称温度场装配
    int DirectSolve();  //直接法求解
    int DirectSolve1();
    int PostProcessing();
    int GenerateMetisMesh(int partition);
    int drawBDR();  //绘制边界
    int DDTLMSolve();   //通过DDTLM求解

private:
    QVector<double> SparseSolve(umat &loc, vec &val, vec F1, int size);
    int m_num_pts;
    int m_num_VtxEle;
    int m_num_EdgEle;
    int m_num_TriEle;
    C2DNode *mp_2DNode;
    CVtxElement *mp_VtxEle;
    CEdgElement *mp_EdgEle;
    CTriElement *mp_TriEle;
    mat S;
    vec F;
    Widget *thePlot;
    QCustomPlot *customplot;
    QCPGraph *graph1;   //绘制边界，分网，负载区域
    QList<QCPCurve*> *mesh;
    QList<QCPCurve*> *mesh1;
    const char *meshfile;
    char metismesh[256];
    int m_num_part{0};  //分区个数
    int *epartTable;    //保存三角形单元在第几个分区
    int *npartTable;    //保存节点在第几个分区
    sp_mat* X;
};

#endif // TEMP2DFEMCORE_H
