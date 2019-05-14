#ifndef TEMP3DFEMCORE_H
#define TEMP3DFEMCORE_H
#include <string>
#include "widget.h"
#include "datatype.h"
class CTemp3DFEMCore
{
public:
    CTemp3DFEMCore(Widget *parent = nullptr, const char *fn = nullptr );
    ~CTemp3DFEMCore();
    int Load3DFEMCOMSOL();
    int preCalculation();   //计算基本几何参数
    int setCondition(Demo showWhat); //根据模型的不同设置边界条件和负载情况
    int StaticAxisAssemble(); //二维轴对称温度场装配
    int DirectSolve();  //直接法求解
    int PostProcessing();
    int GenerateMetisMesh(int partition);
    int drawBDR();  //绘制边界
    int DDTLMSolve();   //通过DDTLM求解

private:
    Widget *m_Widget;
    const char *m_COMSOLMesh;
    int m_num_pts{0};
    int m_num_VtxEle{0};
    int m_num_EdgEle{0};
    int m_num_TriEle{0};
    int m_num_TetEle{0};
    C3DNode *mp_3DNode;
    CVtxElement *mp_VtxEle;
    CEdgElement *mp_EdgEle;
    CTriElement *mp_TriEle;
    CTetElement *mp_TetEle;
};

#endif // TEMP3DFEMCORE_H
