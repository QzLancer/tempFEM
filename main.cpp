#include <QApplication>
#include <QDebug>
#include "temp2dfemcore.h"

void solvesimple(Widget *parent);
void solvecontactor(Widget *parent);
void metistest(Widget *parent, int part);
void bdrtest(Widget *parent);

int main(int argc, char *argv[])
{

    Demo showWhat = METISTEST;
    QApplication a(argc, argv);
    Widget w;
    w.show();
    switch(showWhat){
        case SOLVESIMPLE:
            solvesimple(&w);
        break;
        case SOLVECONTACTOR:
            solvecontactor(&w);
        break;
        case METISTEST:
            metistest(&w, 4);
        break;
        case BDRTEST:
            bdrtest(&w);
        break;


    }

    return a.exec();

}

void solvesimple(Widget *parent){
   CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_heatexcg.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->preCalculation();
    temp->setCondition(SOLVESIMPLE);
    temp->StaticAxisAssemble();
    temp->DirectSolve();
//    temp.PostProcessing(parent);
}

void solvecontactor(Widget *parent){

}

void metistest(Widget *parent, int part){
    CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->preCalculation();
    temp->setCondition(METISTEST);
    temp->StaticAxisAssemble();
    temp->GenerateMetisMesh(part);
    temp->drawBDR();
}

void bdrtest(Widget *parent){
    CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->preCalculation();
    temp->setCondition(METISTEST);
//    temp->drawBdr();
//    temp->drawLoad();
}
