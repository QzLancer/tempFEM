#include "widget.h"
#include <QApplication>
#include "temp2dfemcore.h"


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Widget w;
    w.show();

    CTemp2DFEMCore temp;
    temp.Load2DMeshCOMSOL("D:\\tempFEM\\tempFEM0\\tempFEM\\mesh_heatexcg.mphtxt");
    temp.preCalculation();
    temp.setCondition();
    temp.StaticAxisAssemble();
    temp.DirectSolve();
    return a.exec();
}
