#ifndef DATATYPE_H
#define DATATYPE_H

struct C2DNode{
    double x{0}, y{0};//the co
    double V{0};//the solution
    int bdr{0};//boundary type，根据COMSOL中的设置来定义 1为强迫边界条件，0无影响
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
    double heatflux{0};   //热通量
    double d{0};  //线单元的长度
    double ravg{0};  //线单元的平均半径
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
    double ravg{0};//轴对称模型时，单元的平均半径
    bool LinearFlag{0};//定义逻辑变量LinearFlag，用来判断具体单元是否处于线性区域
};
#endif // DATATYPE_H
