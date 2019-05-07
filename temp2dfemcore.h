#ifndef TEMP2DFEMCORE_H
#define TEMP2DFEMCORE_H
#include "datatype.h"
#include "armadillo"
#include <QWidget>

using namespace arma;
class CTemp2DFEMCore
{
    const double TEMP = 273.15;
public:
    CTemp2DFEMCore();
    ~CTemp2DFEMCore();
    int Load2DMeshCOMSOL(const char *fn);   //加载COMSOL分网文件
    int preCalculation();   //计算基本几何参数
    int setCondition(); //设置边界条件和负载情况，此处参考的模型为二维轴对称模型
    int StaticAxisAssemble(); //二维轴对称温度场装配
    int DirectSolve();  //直接法求解
    int PostProcessing(QWidget *parent = nullptr);

private:
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
};

#endif // TEMP2DFEMCORE_H
