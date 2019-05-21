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
void ddtlm3d(Widget *parent, int part);
void solvecontactorNR(Widget *parent);
void solve3dcontactorNR(Widget *parent);
void NRddtlm(Widget *parent, int part);
void NRddtlm3d(Widget *parent, int part);

int main(int argc, char *argv[])
{

    Demo showWhat = NRDDTLM3D;
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
            ddtlm(&w, 4);
        break;
        case SOLVE3DCONTACTOR:
            solve3dcontactor(&w);
        break;
        case DDTLM3D:
            ddtlm3d(&w, 4);
        break;
        case SOLVECONTACTORNR:
            solvecontactorNR(&w);
        break;
        case SOLVE3DCONTACTORNR:
            solve3dcontactorNR(&w);
        break;
        case NRDDTLM:
            NRddtlm(&w, 4);
        break;
        case NRDDTLM3D:
            NRddtlm3d(&w, 4);
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
    temp->setCondition();
    temp->Static3DAssemble();
}

void ddtlm3d(Widget *parent, int part)
{
    CTemp3DFEMCore *temp = new CTemp3DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor3D.mphtxt");
    temp->Load3DFEMCOMSOL();
    temp->GenerateMetisMesh(part);
    temp->preCalculation();
    temp->setCondition();
    temp->DDTLM3DSolve1();
}

void solvecontactorNR(Widget *parent){
    CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->preCalculation();
    temp->setCondition(SOLVECONTACTORNR);
    temp->NRSolve();
}

void solve3dcontactorNR(Widget *parent)
{
    CTemp3DFEMCore *temp = new CTemp3DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor3D.mphtxt");
    temp->Load3DFEMCOMSOL();
    temp->preCalculation();
    temp->setCondition();
    temp->NRSolve();
}

void NRddtlm(Widget *parent, int part){
    CTemp2DFEMCore *temp = new CTemp2DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor.mphtxt");
    temp->Load2DMeshCOMSOL();
    temp->GenerateMetisMesh(part);
    temp->preCalculation();
    temp->setCondition(SOLVECONTACTORNR);
    temp->DRDDTLMSolve();
}

void NRddtlm3d(Widget *parent, int part){
    CTemp3DFEMCore *temp = new CTemp3DFEMCore(parent, "..\\tempFEM\\model\\mesh_contactor3D.mphtxt");
    temp->Load3DFEMCOMSOL();
    temp->GenerateMetisMesh(part);
    temp->preCalculation();
    temp->setCondition();
    temp->NRDDTLMSolve();
}
