#include <QApplication>
#include <QDebug>
#include "temp2dfemcore.h"
#include "temp3dfemcore.h"
void solvesimple(Widget *parent);
void solvecontactor(Widget *parent);
void metistest(Widget *parent, int part);
void bdrtest(Widget *parent);
void ddtlm(Widget *parent, int part);
void solve3dcontactor(Widget *parent);

int main(int argc, char *argv[])
{

    Demo showWhat = SOLVE3DCONTACTOR;
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
        case DDTLM:
            ddtlm(&w, 8);
        case SOLVE3DCONTACTOR:
            solve3dcontactor(&w);
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
    CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->preCalculation();
    temp->setCondition(SOLVECONTACTOR);
    temp->StaticAxisAssemble();
    temp->DirectSolve1();
}

void metistest(Widget *parent, int part){
    CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->preCalculation();
    temp->setCondition(SOLVECONTACTOR);
    temp->StaticAxisAssemble();
    temp->GenerateMetisMesh(part);
    temp->drawBDR();
}

void bdrtest(Widget *parent){
    CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->preCalculation();
    temp->setCondition(SOLVECONTACTOR);
    temp->drawBDR();
//    temp->drawLoad();
}

void ddtlm(Widget *parent, int part){
    CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->preCalculation();
    temp->setCondition(SOLVECONTACTOR);
    temp->GenerateMetisMesh(part);
    temp->DDTLMSolve();
}

void solve3dcontactor(Widget *parent){
    CTemp3DFEMCore *temp = new CTemp3DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor3D.mphtxt");
    temp->Load3DFEMCOMSOL();
    temp->preCalculation();
}
