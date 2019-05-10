#ifndef DATATYPE_H
#define DATATYPE_H

struct C2DNode{
    double x{0}, y{0};//the co
    double V{0};//the solution
    int bdr{0};//boundary type，根据COMSOL中的设置来定义 0无影响，1为强迫边界条件，用于求解，2和3分别为第二类边界条件和第三类边界条件，用于区域分解
};

struct CVtxElement{
    int n{0};
    int domain{0};
};

struct CEdgElement{
    int n[2]{0};
    double x[2]{0};
    double y[2]{0};
    int domain{0};
    double h{0};   //传热系数
    double Text{0};    //外部温度
    double heatflux{0};   //热通量
    double d{0};  //线单元的长度
    double xavg{0};  //线单元的平均半径
    int bdr{0};    //边界条件
};

struct CTriElement{
    int n[3]{0};// ni, nj, nk;//
    int domain{0};
    double x[3]{0};
    double y[3]{0};
    double q[3]{0};// Qi, Qj, Qk;
    double r[3]{0};
    double cond{0};
    double Area{0};//单元的面积
    double xavg{0};//轴对称模型时，单元的平均半径
    double source{0};
    bool LinearFlag{0};//定义逻辑变量LinearFlag，用来判断具体单元是否处于线性区域
    int Material{0};
};

struct CInterfacePoint{
    double Y0{0};   //传输线导纳
    double Vi{0};   //入射电压或者反射电压
};

class ResistMarix{
public:
    double Y11;
    double Y12;
    double Y13;
    double Y22;
    double Y23;
    double Y33;
};

enum Demo{
    SOLVESIMPLE,    //直接法求解简单模型
    SOLVECONTACTOR, //直接法求解接触器模型
    METISTEST,  //区域分解API测试
    BDRTEST, //查看负载和边界区域是否正确
    DDTLM   //DDTLM求解
};
#endif // DATATYPE_H
