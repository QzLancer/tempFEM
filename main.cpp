#include "widget.h"
#include <QApplication>

int Drawnephogram();

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Widget w;
    w.show();

    CTemp2DFEMCore temp;
    temp.Load2DMeshCOMSOL("D:\\tempFEM\\tempFEM0\\tempFEM\\mesh_contactor.mphtxt");

    return a.exec();
}
