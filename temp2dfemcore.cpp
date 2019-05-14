#include "Temp2DFEMCore.h"
#include <stdio.h>
#include <stdlib.h>
#include <QDebug>
#include <math.h>
#include <iostream>
#include <vector>
#include "mpmetis.h"
#include <QPen>
#include "slu_ddefs.h"
#include <iomanip>
#include <omp.h>
#include <QVector>

CTemp2DFEMCore::CTemp2DFEMCore(Widget *parent, const char *fn):mp_2DNode(nullptr),
    mp_VtxEle(nullptr),
    mp_EdgEle(nullptr),
    mp_TriEle(nullptr),
    thePlot(parent),
    customplot(thePlot->customplot),
    mesh(new QList<QCPCurve*>),
    mesh1(new QList<QCPCurve*>),
    meshfile(fn),
    epartTable(new int(0)),
    npartTable(new int(0))
{

}

CTemp2DFEMCore::~CTemp2DFEMCore()
{
    if(!mp_2DNode){
        delete mp_2DNode;
        mp_2DNode = nullptr;
    }
    if(!mp_VtxEle){
        delete mp_VtxEle;
        mp_VtxEle = nullptr;
    }
    if(!mp_EdgEle){
        delete mp_EdgEle;
        mp_EdgEle = nullptr;
    }
    if(!mp_TriEle){
        delete mp_TriEle;
        mp_TriEle = nullptr;
    }
}

int CTemp2DFEMCore::Load2DMeshCOMSOL()
{
    FILE *fp = nullptr;
    int err;
    char ch[256];
    err = fopen_s(&fp, meshfile, "r");
    if (!fp) {
        qDebug() << "Error: openning file!";
        return 1;
    }

    //--------------Read the head-----------------------------
    for (int i = 0; i < 18; i++) {
        fgets(ch, 256, fp);
    }
    //-----------------mesh point-----------------------------
    if (fscanf_s(fp, "%d # number of mesh points\n", &m_num_pts) != 1) {
        qDebug() << "Error: reading num_bdr_ns!";
        return 1;
    }
    else qDebug() << m_num_pts << "number of mesh points.";
    mp_2DNode = new C2DNode[m_num_pts];
    int pts_ind;//the beginning of the points index

    if (fscanf_s(fp, "%d # lowest mesh point index\n", &pts_ind) != 1) {
        qDebug() << "Error: reading pts_ind!";
        return 1;
    }
    fgets(ch, 256, fp);

    for (int i = pts_ind; i < m_num_pts; i++) {
        //读取x,y坐标
        if (fscanf_s(fp, "%lf %lf \n", &(mp_2DNode[i].x), &(mp_2DNode[i].y)) != 2) {
            qDebug() << "Error: reading mesh point!";
            return 1;
        }
        //        else{
        //            qDebug() << mp_2DNode[i].x << " " << mp_2DNode[i].y << "\n";
        //        }
    }
    //---------------vertexnode-------------------------------
    for (int i = 0; i < 7; i++)
        fgets(ch, 256, fp);
    int num_vtx_ns;

    if (fscanf_s(fp, "%d # number of nodes per element\n", &num_vtx_ns) != 1) {
        qDebug() << "Error: reading num_vtx_ns!";
        return 1;
    }

    if (fscanf_s(fp, "%d # number of elements\n", &m_num_VtxEle) != 1) {
        qDebug() << "Error: reading m_num_VtxEle!";
        return 1;
    }
    else qDebug() << m_num_VtxEle <<  "number of Vertex elements.";
    fgets(ch, 256, fp);
//    qDebug() << m_num_VtxEle;
    mp_VtxEle = new CVtxElement[m_num_VtxEle];
    for (int i = 0; i < m_num_VtxEle; i++) {
        if (fscanf_s(fp, "%d \n", &((mp_VtxEle + i)->n)) != 1) {
            qDebug() << "Error: reading vertex element points!";
            return 1;
        }
    }
    //---------------vertexdomain-------------------------------
    for (int i = 0; i < 2; i++)
        fgets(ch, 256, fp);
    for (int i = 0; i < m_num_VtxEle; i++){
        if (fscanf_s(fp, "%d \n", &(mp_VtxEle[i].domain)) != 1){
            qDebug() << "Error: reading vertex domain!";
            return 1;
        }
        else{
            mp_VtxEle[i].domain++;
        }
    }
    //----------------edgnode-----------------------------------
    for (int i = 0; i < 6; i++){
        fgets(ch, 256, fp);
    }
    if(fscanf_s(fp, "%d # number of elements\n", &m_num_EdgEle) != 1){
        qDebug() <<  "Error: reading m_num_EdgEle";
        return 1;
    }
    else qDebug() << m_num_EdgEle << "number of Edge elements.";
    mp_EdgEle = new CEdgElement[m_num_EdgEle];
    fgets(ch, 256, fp);
    for (int i = 0; i < m_num_EdgEle; i++){
        if(fscanf_s(fp, "%d %d\n", &mp_EdgEle[i].n[0], &mp_EdgEle[i].n[1]) != 2){
            qDebug() << "Error: reading edg element points!" ;
            return 1;
        }
    }
    //----------------edgdomain---------------------------------
    for(int i = 0; i < 2; i++){
        fgets(ch, 256, fp);
    }
    for (int i = 0; i < m_num_EdgEle; i++){
        if(fscanf_s(fp, "%d \n", &mp_EdgEle[i].domain) != 1){
            qDebug() << "Error: reading edgdomain!";
            return 1;
        }
        else{
            mp_EdgEle[i].domain++;
//            qDebug() << "mp_EdgEle: " << mp_EdgEle[i].domain;
        }
    }
    //----------------trinode-----------------------------------
    for (int i = 0; i < 6; i++){
        fgets(ch, 256, fp);
    }
    if(fscanf_s(fp, "%d # number of elements\n", &m_num_TriEle) != 1){
        qDebug() << "Error: reading m_num_TriEle!";
        return 1;
    }
    else qDebug() << m_num_TriEle << "number of Triangle elements.";
    fgets(ch, 256, fp);
    mp_TriEle = new CTriElement[m_num_TriEle];
    for (int i = 0; i < m_num_TriEle; i++){
        if (fscanf_s(fp, "%d %d %d \n", &mp_TriEle[i].n[0], &mp_TriEle[i].n[1], &mp_TriEle[i].n[2]) != 3) {
            qDebug() << "Error: reading elements points!";
            return 1;
        }
    }
    //----------------tridomain---------------------------------
    for (int i = 0; i < 2; i++){
        fgets(ch, 256, fp);
    }
    for (int i = 0; i < m_num_TriEle; i++){
        if(fscanf_s(fp, "%d \n", &mp_TriEle[i].domain) != 1){
            qDebug() << "Error: reading tridomain!";
            return 1;
        }
        else{
//            mp_TriEle[i].domain++;
        }
    }
    fclose(fp);
    S = zeros<mat>(m_num_pts, m_num_pts);
    F = zeros<vec>(m_num_pts);

    qDebug() << "Load2DMeshCOMSOL successfully!";
    return 0;
}

int CTemp2DFEMCore::preCalculation()
{
    //计算所有三角形单元的q,r,面积Area,平均半径xavg(轴对称模型使用)
    for(int i = 0; i < m_num_TriEle; i++){
        for(int j = 0; j < 3; j++){
            mp_TriEle[i].x[j] = mp_2DNode[mp_TriEle[i].n[j]].x;
            mp_TriEle[i].y[j] = mp_2DNode[mp_TriEle[i].n[j]].y;
        }
        mp_TriEle[i].q[0] = mp_TriEle[i].y[1]-mp_TriEle[i].y[2];
        mp_TriEle[i].q[1] = mp_TriEle[i].y[2]-mp_TriEle[i].y[0];
        mp_TriEle[i].q[2] = mp_TriEle[i].y[0]-mp_TriEle[i].y[1];
        mp_TriEle[i].r[0] = mp_TriEle[i].x[2]-mp_TriEle[i].x[1];
        mp_TriEle[i].r[1] = mp_TriEle[i].x[0]-mp_TriEle[i].x[2];
        mp_TriEle[i].r[2] = mp_TriEle[i].x[1]-mp_TriEle[i].x[0];
        mp_TriEle[i].Area = abs(mp_TriEle[i].q[0]*mp_TriEle[i].r[1]-mp_TriEle[i].q[1]*mp_TriEle[i].r[0])/2.;
        mp_TriEle[i].xavg = (mp_TriEle[i].x[0]+mp_TriEle[i].x[1]+mp_TriEle[i].x[2])/3;
    }
    //计算所有线单元的长度d，平均半径xavg(轴对称模型使用)
    for(int i = 0; i < m_num_EdgEle; i++){
        for(int j = 0; j < 2; j++){
            mp_EdgEle[i].x[j] = mp_2DNode[mp_EdgEle[i].n[j]].x;
            mp_EdgEle[i].y[j] = mp_2DNode[mp_EdgEle[i].n[j]].y;
        }
        mp_EdgEle[i].d = std::sqrt(pow(mp_EdgEle[i].x[0]-mp_EdgEle[i].x[1],2)+
                pow(mp_EdgEle[i].y[0]-mp_EdgEle[i].y[1],2));
        mp_EdgEle[i].xavg = (mp_EdgEle[i].x[0]+mp_EdgEle[i].x[1])/2;
    }

   //观察求解的值是否正确
//    std::cout << "mp_TriEle.Area:" << "\n";
//    for(int i = 0; i< m_num_TriEle; i++){
//        for(int j = 0; j < 3; j++){
//            std::cout << mp_TriEle[i].Area << " ";
//        }
//        std::cout << "\n";
//    }
    return 0;
}

//不同模型主要改写的是这一部分
int CTemp2DFEMCore::setCondition(Demo showWhat)
{
    if(showWhat == SOLVESIMPLE){
        qDebug() << "setCondition:SOLVESIMPLE.";
        //本模型中只有第一类和第二类边界条件，通过检索线单元来处理
        for (int i = 0; i < m_num_EdgEle; i++){
            if(mp_EdgEle[i].domain == 1 || mp_EdgEle[i].domain == 4){
                mp_EdgEle[i].bdr = 0;
            }
            else if(mp_EdgEle[i].domain == 2 || mp_EdgEle[i].domain == 5 || mp_EdgEle[i].domain == 6){
                mp_EdgEle[i].bdr = 1;
                for(int j = 0; j < 2; j++){
                    mp_2DNode[mp_EdgEle[i].n[j]].bdr = 1;
                    mp_2DNode[mp_EdgEle[i].n[j]].V = 273.15;
                }
            }
            else if(mp_EdgEle[i].domain == 3){
                mp_EdgEle[i].bdr = 2;
                mp_EdgEle[i].heatflux = 500000;
            }
        }
        //三角形单元参数
        for (int i = 0; i < m_num_TriEle; i++)
            mp_TriEle[i].cond = 52;

    //    //调试
    //    for (int i = 0; i < m_num_pts; i++)
    //        cout << mp_2DNode[i].bdr << endl;
    }
    else if(showWhat == SOLVECONTACTOR){
        qDebug() << "setcondition:SOLVEDDTLM.";
        //热源设置
        for(int i = 0; i < m_num_TriEle; i++){
            if((mp_TriEle[i].domain == 10) | (mp_TriEle[i].domain == 11)){
                mp_TriEle[i].source = 500000;
            }
            else{
                mp_TriEle[i].source = 0;
            }
        }
        //热导率设置
        for(int i = 0; i < m_num_TriEle; i++){
            if((mp_TriEle[i].domain == 10) | (mp_TriEle[i].domain == 11)){
                mp_TriEle[i].cond = 400;
                mp_TriEle[i].Material = 0;
//                qDebug() << "Domain1: " << i;
            }
            else if((mp_TriEle[i].domain == 2) | (mp_TriEle[i].domain == 3) | (mp_TriEle[i].domain == 4) | (mp_TriEle[i].domain == 6) | (mp_TriEle[i].domain == 8)){
                mp_TriEle[i].cond = 76.2;
                mp_TriEle[i].Material = 1;
//                qDebug() << "Domain2: " << i;
            }
            else if((mp_TriEle[i].domain == 1) | (mp_TriEle[i].domain == 7) | (mp_TriEle[i].domain == 9)){
                mp_TriEle[i].cond = 0.26;
                mp_TriEle[i].Material = 2;
//                qDebug() << "Domain3: " << i;
            }
            else if((mp_TriEle[i].domain == 5) | (mp_TriEle[i].domain == 12) | (mp_TriEle[i].domain == 13)){
                mp_TriEle[i].cond = 0.03;
                mp_TriEle[i].Material = 3;
//                qDebug() << "Domain4: " << i;
            }   //实际上这里应该是空气，暂时用尼龙的热导率代替
        }
        //第三类边界条件设置
        for(int i = 0; i< m_num_EdgEle; i++){
            if((mp_EdgEle[i].domain == 2) | (mp_EdgEle[i].domain == 9) | (mp_EdgEle[i].domain == 49)){
                mp_EdgEle[i].bdr = 3;
                mp_EdgEle[i].h = 20;
                mp_EdgEle[i].Text = 293.15;
//                qDebug() << i;
            }
//            cout << mp_EdgEle[i].bdr << endl;
        }
    }
    return 0;
}

using namespace arma;

int CTemp2DFEMCore::StaticAxisAssemble()
{
    int pos = 0;
    int numbdr = 0;
    for(int i = 0; i < m_num_EdgEle; i++){
        if(mp_EdgEle[i].bdr == 3) numbdr++;
    }
    umat locs(2, 9*m_num_TriEle+4*numbdr);
    mat vals(1, 9*m_num_TriEle+4*numbdr);
    std::ofstream mycoutS("..\\tempFEM\\test\\S.txt");
    std::ofstream mycoutF("..\\tempFEM\\test\\F.txt");
    std::ofstream mycouth("..\\tempFEM\\test\\h.txt");
    std::ofstream mycoutText("..\\tempFEM\\test\\Text.txt");
    std::ofstream mycoutd("..\\tempFEM\\test\\d.txt");
    std::ofstream mycoutFl("..\\tempFEM\\test\\Fl.txt");
    std::ofstream mycoutxavg("..\\tempFEM\\test\\xavg.txt");
    const double PI = 3.14159265358979323846;
    //三角形单元装配过程
//    mat Se = zeros<mat>(3,3);
    double Se;
    double Fe;
    double Sl;
//    double Fl;
    vec Fl = zeros<vec>(2);
    for(int k = 0; k < m_num_TriEle; k++){
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                Se = (PI*mp_TriEle[k].cond*mp_TriEle[k].xavg*(mp_TriEle[k].r[i]*mp_TriEle[k].r[j]+mp_TriEle[k].q[i]*mp_TriEle[k].q[j]))/(2*mp_TriEle[k].Area);
                S(mp_TriEle[k].n[i],mp_TriEle[k].n[j]) = S(mp_TriEle[k].n[i],mp_TriEle[k].n[j]) + Se;
                locs(0, pos) = mp_TriEle[k].n[i];
                locs(1, pos) = mp_TriEle[k].n[j];
                vals(0, pos) = Se;
                ++pos;
            }
            Fe = PI*mp_TriEle[k].source*mp_TriEle[k].Area*(mp_TriEle[k].x[i]+3.*mp_TriEle[k].xavg)/6.;
            F(mp_TriEle[k].n[i]) = F(mp_TriEle[k].n[i]) + Fe;
        }
    }
    qDebug() << rank(S)<<S.size();



    //线单元装配过程
    for(int k = 0; k < m_num_EdgEle; k++){
        for(int i = 0; i < 2; i++){
            if(mp_EdgEle[k].bdr == 3){
                for(int j = 0; j < 2; j++){
                    if(i == j){
                        Sl = PI*mp_EdgEle[k].h*mp_EdgEle[k].d*(2*mp_EdgEle[k].xavg+2*mp_EdgEle[k].x[i])/6.;
                    }
                    else{
                        Sl = PI*mp_EdgEle[k].h*mp_EdgEle[k].d*(2*mp_EdgEle[k].xavg)/6.;
                    }
//                    cout << Sl << endl;
                    S(mp_EdgEle[k].n[i], mp_EdgEle[k].n[j]) = S(mp_EdgEle[k].n[i], mp_EdgEle[k].n[j]) + Sl;
                    locs(0, pos) = mp_EdgEle[k].n[i];
                    locs(1, pos) = mp_EdgEle[k].n[j];
                    vals(0, pos) = Sl;
                    ++pos;
                }
                Fl(i) = PI*mp_EdgEle[k].h*mp_EdgEle[k].Text*mp_EdgEle[k].d*(2*mp_EdgEle[k].xavg+mp_EdgEle[k].x[i])/3.;
                F(mp_EdgEle[k].n[i]) = F(mp_EdgEle[k].n[i])+Fl(i);
            }else if(mp_EdgEle[k].bdr == 2){
                Fl(i) = PI*mp_EdgEle[k].heatflux*mp_EdgEle[k].d*(mp_EdgEle[k].x[i]+2.*mp_EdgEle[k].xavg)/3.;
                F(mp_EdgEle[k].n[i]) = F(mp_EdgEle[k].n[i])+Fl(i);
            }
        }
    }
    qDebug() << "Assemble successful!";
//    cout << locs;
//    cout << vals;
    X = new sp_mat(true, locs, vals, m_num_pts, m_num_pts, true, true);

//    mycoutS << S;
//    qDebug()<< F(2)<<aaa;

//     printf("%.10lf",S(10,10));
//    mycoutFl << Fl << endl;
//    qDebug() << rank(S);
//    FILE* fp;
//    fopen_s(&fp, "D:\\tempFEM\\tempFEM0\\tempFEM\\S.txt", "w");;
//    for(int i=0;i < m_num_pts;++i){
//        for(int j=0;j < m_num_pts;++j){
//            fprintf_s(fp,"%.10lf\t",S(i,j));
//        }
//        fprintf_s(fp,"\n");
//    }
//    fclose(fp);
    return 0;
}

//通过直接法来求解包含各类边界条件的温度场问题
int CTemp2DFEMCore::DirectSolve()
{
    vec A = zeros<vec>(m_num_pts);
    std::vector<int> v1;
    //处理第一类边界条件
    for(int i = 0; i < m_num_pts; i++){
        if(mp_2DNode[i].bdr == 1){  //如果节点处于边界上，那么该节点为一个固定值
            A(i) = TEMP;
        }
        else{
            v1.push_back(i);
        }
    }
    //删掉第一类边界条件包含的行，重构S和F
    int j = v1.size();
    vec A1 = zeros<vec>(j);
    qDebug() << "test";
    mat S1 = zeros<mat>(j, m_num_pts);
    mat S2 = zeros<mat>(j, j);
    vec F1 = zeros<vec>(j);
    vec F2 = zeros<vec>(j);
    for(int i = 0; i < j; i++){
        S1.row(i) = S.row(v1[i]);
        F1(i) = F(v1[i]);
    }
    F2 = S1*A;
    F1 = F1 - F2;
    for(int i = 0; i < j; i++){
        S2.col(i) = S1.col(v1[i]);
    }
    A1 = solve(S2, F1);

    for(int i = 0; i < j; i++){
        A(v1[i]) = A1[i];
        mp_2DNode[v1[i]].V = A1[i];
    }
    for(int i = 0; i < m_num_pts; i++){
        cout << mp_2DNode[i].V << endl;
    }

    return 0;
}

int CTemp2DFEMCore::DirectSolve1()
{
    SuperMatrix sluA;
    NCformat *Astore;
    double   *a;
    int      *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SuperMatrix U;      /* factor U */
    SuperMatrix B;
    int      nrhs, ldx, info, m, n, nnz;
    double   *rhs;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;

    set_default_options(&options);

    /* create matrix A in Harwell-Boeing format.*/
    m = m_num_pts; n = m_num_pts; nnz = X->n_nonzero;
    a = const_cast<double *>(X->values);
    asub = (int*)const_cast<unsigned int*>(X->row_indices);
    xa = (int*)const_cast<unsigned int*>(X->col_ptrs);
    dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = (NCformat *)sluA.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", sluA.nrow, sluA.ncol, Astore->nnz);

    nrhs = 1;
    if (!(rhs = doubleMalloc(m * nrhs))) ABORT("Malloc fails for rhs[].");
    //将内存拷贝过来
    //memmove(rhs, unknown_b, 5*sizeof(double));
    for (int i = 0; i < m; i++){
        rhs[i] = F(i);
    }
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");
    if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);

    dgssv(&options, &sluA, perm_c, perm_r, &L, &U, &B, &stat, &info);
    if (info == 0) {
        qDebug()<<"Ok.";
        /* This is how you could access the solution matrix. */
        double *sol = (double*)((DNformat*)B.Store)->nzval;
        for(int i = 0; i < m_num_pts; ++i){
            mp_2DNode[i].V = sol[i];
            qDebug() << mp_2DNode[i].V;
        }
    }else {
        qDebug() << "info = " << info;
    }

    SUPERLU_FREE(rhs);
//    SUPERLU_FREE(xact);
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
//    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    return 0;
}

int CTemp2DFEMCore::PostProcessing()
{
//    QCustomPlot *customPlot = new QCustomPlot(parent);
//    QVBoxLayout *layout = new QVBoxLayout(parent);
//    layout->addWidget(customPlot);
//    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
//    customPlot->axisRect()->setupFullAxesBox(true);
//    customPlot->xAxis->setLabel("x");
//    customPlot->yAxis->setLabel("y");
//    customPlot->xAxis->setRange(0.02,0.1);
//    customPlot->yAxis->setRange(0,0.14);
//    customPlot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
//    customPlot->axisRect()->setupFullAxesBox(true);
//    customPlot->xAxis->setLabel("x");
//    customPlot->yAxis->setLabel("y");
//    //二维插值


//    // set up the QCPColorMap:
//    QCPColorMap *colorMap = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
//    int nx = 200;
//    int ny = 200;
//    colorMap->data()->setSize(nx, ny); // we want the color map to have nx * ny data points
//    colorMap->data()->setRange(QCPRange(-4, 4), QCPRange(-4, 4)); // and span the coordinate range -4..4 in both key (x) and value (y) dimensions
//    // now we assign some data, by accessing the QCPColorMapData instance of the color map:
//    double x, y, z;
//    for (int xIndex=0; xIndex<nx; ++xIndex)
//    {
//      for (int yIndex=0; yIndex<ny; ++yIndex)
//      {
//        colorMap->data()->cellToCoord(xIndex, yIndex, &x, &y);
//        cout << x << " " << y << endl;
//        double r = 3*qSqrt(x*x+y*y)+1e-2;
//        z = 2*x*(qCos(r+2)/r-qSin(r+2)/r); // the B field strength of dipole radiation (modulo physical constants)
//        colorMap->data()->setCell(xIndex, yIndex, z);
//      }
//    }

//    // add a color scale:
//    QCPColorScale *colorScale = new QCPColorScale(customPlot);
//    customPlot->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
//    colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
//    colorMap->setColorScale(colorScale); // associate the color map with the color scale
//    colorScale->axis()->setLabel("Magnetic Field Strength");

//    // set the color gradient of the color map to one of the presets:
//    colorMap->setGradient(QCPColorGradient::gpPolar);
//    // we could have also created a QCPColorGradient instance and added own colors to
//    // the gradient, see the documentation of QCPColorGradient for what's possible.

//    // rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient:
//    colorMap->rescaleDataRange();

//    // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
//    QCPMarginGroup *marginGroup = new QCPMarginGroup(customPlot);
//    customPlot->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
//    colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

//    // rescale the key (x) and value (y) axes so the whole color map is visible:
//    customPlot->rescaleAxes();
    return 0;
}

int CTemp2DFEMCore::GenerateMetisMesh(int partition)
{
    m_num_part = partition;
    if(meshfile == nullptr){
        return 1;
    }
    FILE *fp = nullptr;

//    QVBoxLayout *layout = new QVBoxLayout(thePlot);
//    layout->addWidget(customplot);
//    customplot->xAxis->setLabel("x");
//    customplot->xAxis->setRange(0, 0.09);
//    customplot->xAxis->setAutoTickStep(false);
//    customplot->xAxis->setTicks(false);
//    customplot->yAxis->setLabel("y");
//    customplot->yAxis->setRange(-0.09, 0.09);
//    customplot->xAxis2->setTicks(false);
    customplot->yAxis->setScaleRatio(customplot->xAxis, 1.0);

//    customplot->yAxis->setAutoTickStep(false);
//    customplot->yAxis->setAutoTickLabels(false);
//    customplot->yAxis->setTicks(false);
//    customplot->yAxis->grid()->setVisible(false);
//    customplot->yAxis->grid()->setZeroLinePen(Qt::NoPen);
//    customplot->yAxis2->setTicks(false);
//    customplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    //生成metis的输入分网文件
    strcpy(metismesh,meshfile);
    strcat(metismesh, ".metis");
    fp = fopen(metismesh, "w+");
    if(!fp){
        return 1;
    }
    //写入单元数目和单元节点编号（从1开始）
    fprintf(fp, "%d\n", m_num_TriEle);
    for(int i = 0; i < m_num_TriEle; i++){
        fprintf(fp, "%d %d %d \n", mp_TriEle[i].n[0]+1, mp_TriEle[i].n[1]+1, mp_TriEle[i].n[2]+1);
    }
    fclose(fp);
    //调用metis进行分区
    char part[3];
    sprintf(part, "%d", partition);
    char a[] = "QzLancer";
    char * myargv[] = {a, metismesh, part};
    mpmetis(3,myargv);
    //读入单元分区表
    char epartname[256];
    sprintf(epartname,"%s.epart.%d",metismesh,partition);
    fp = fopen(epartname,"r");
    if(!fp){
        qDebug()<<"Error in open epart file!";
        return 1;
    }
    if(epartTable){
        free(epartTable);//释放之前申请的空间
    }
    epartTable = (int*)malloc(m_num_TriEle*sizeof(int));
    for(int i=0;i < m_num_TriEle;++i){
        fscanf(fp,"%d\n",epartTable+i);
    }
    fclose(fp);
    //读入节点分区表
    char npartname[256];
    sprintf(npartname,"%s.npart.%d",metismesh,partition);
    fp = fopen(npartname,"r");
    if(!fp){
        qDebug()<<"Error in open npart file!";
        return 1;
    }
    if(npartTable){
        free(npartTable);//释放之前申请的空间
    }
    npartTable = (int*)malloc(m_num_pts*sizeof(int));
    for(int i=0;i < m_num_pts;++i){
        fscanf(fp,"%d\n",npartTable+i);
    }
    fclose(fp);
    //绘制分区之后的分网
    QColor cc[7];
    cc[0] = QColor(150, 0, 0);
    cc[1] = QColor(0, 150, 0);
    cc[2] = QColor(0, 0, 150);
    cc[3] = QColor(150, 150, 0);
    cc[4] = QColor(0, 150, 150);
    cc[5] = QColor(150, 0, 150);
    cc[6] = QColor(150, 150, 150);

    for(int i=0;i < m_num_TriEle;++i){
        QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
        mesh->push_back(newCurve);
        QVector <double> x1(4);
        QVector <double> y1(4);
        for (int j = 0; j < 3; j++) {
            x1[j] = mp_2DNode[mp_TriEle[i].n[j]].x;
            y1[j] = mp_2DNode[mp_TriEle[i].n[j]].y;
        }
        x1[3] = x1[0];
        y1[3] = y1[0];
        newCurve->setBrush(cc[epartTable[i]]);
        newCurve->setData(x1,y1);
    }
    customplot->replot();
    thePlot->setWindowTitle("区域分解");
//    qApp->processEvents();//强制刷新界面
//    for(int i = 0; i < m_num_TriEle; ++i){
//        qDebug() << epartTable[i];
//    }
    return 0;
}

int CTemp2DFEMCore::drawBDR()
{
    Widget *w1 = new Widget();
    QCustomPlot *customplot1 = w1->customplot;
    w1->setWindowTitle("边界单元和负载域");
    w1->show();
    customplot1->yAxis->setScaleRatio(customplot->xAxis, 1.0);

    QColor cc[7];
    cc[0] = QColor(150, 0, 0);
    cc[1] = QColor(0, 150, 0);
    cc[2] = QColor(0, 0, 150);
    cc[3] = QColor(150, 150, 0);
    cc[4] = QColor(0, 150, 150);
    cc[5] = QColor(150, 0, 150);
    cc[6] = QColor(150, 150, 150);


    //    不同材质区域绘制
        for(int i = 0; i < m_num_TriEle; i++){
            QCPCurve *curve2 = new QCPCurve(customplot1->xAxis, customplot1->yAxis); //curve2负责负载域绘制
            curve2->setLineStyle(QCPCurve::lsNone );
            mesh1->push_back(curve2);
            QVector <double> x1(4);
            QVector <double> y1(4);
            for (int j = 0; j < 3; j++) {
                x1[j] = mp_2DNode[mp_TriEle[i].n[j]].x;
                y1[j] = mp_2DNode[mp_TriEle[i].n[j]].y;
            }
            x1[3] = x1[0];
            y1[3] = y1[0];
            curve2->setBrush(cc[mp_TriEle[i].Material]);
            curve2->setData(x1,y1);
        }

//边界绘制
//    qDebug() << m_num_EdgEle;
    for(int i = 0; i < m_num_EdgEle; i++){
         QCPGraph *graph1 = customplot1->addGraph(); //graph1负责边界单元绘制
         QVector <double> x1(2);
         QVector <double> y1(2);
         for(int j = 0; j < 2; j++){
             x1[j] = mp_2DNode[mp_EdgEle[i].n[j]].x;
             y1[j] = mp_2DNode[mp_EdgEle[i].n[j]].y;
         }
         QPen graph1pen;
         if(mp_EdgEle[i].bdr == 0){
             graph1pen.setColor(Qt::white);
         }
         else graph1pen.setColor(Qt::black);
         graph1pen.setWidthF(3);
         graph1->setPen(graph1pen);
         graph1->setData(x1,y1);
    }



    customplot1->replot();

    return 0;
}

int CTemp2DFEMCore::DDTLMSolve()
{
    //PART I:处理分区信息
    //A :查找各个分区上的交界点，并且去除掉边界点
    //1.创建partition个大小为num_pts的数组，初始化为-1
    //-1表示该位置的节点不在该分区内，方便交界点过来询问
    int * TriEle_num_part = (int*)calloc(m_num_part, sizeof(int));
    int ** npart = (int**)malloc(m_num_part * sizeof(int*));
    for(int i = 0;i < m_num_part;++i){
        npart[i] = (int*)malloc(m_num_pts * sizeof(int));
        if(npart[i]){
            for(int j = 0;j < m_num_pts;++j){
                npart[i][j] = -1;
            }
        }
    }

    //2.遍历所有单元，得到交界点列表，保存的是原始编号
    QList<int> interfacePoints;
    int IFPOINT = 1000;//标记为交界点
    for(int i = 0;i < m_num_TriEle;++i){
//        qDebug() << "epartTable[" << i << "]= " << epartTable[i];
        if(epartTable[i] != npartTable[mp_TriEle[i].n[0]]){
            npartTable[mp_TriEle[i].n[0]] = IFPOINT;

        }
        if(epartTable[i] != npartTable[mp_TriEle[i].n[1]]){
            npartTable[mp_TriEle[i].n[1]] = IFPOINT;
        }
        if(epartTable[i] != npartTable[mp_TriEle[i].n[2]]){
            npartTable[mp_TriEle[i].n[2]] = IFPOINT;
        }

        //计算每个分区内的单元数目
        TriEle_num_part[epartTable[i]]++;
    }
    for(int i = 0;i < m_num_pts;++i){
        if(npartTable[i] == IFPOINT){
            interfacePoints.push_back(i);
            //qDebug()<<i;
        }
    }
    qDebug()<<"interface points size: "<<interfacePoints.size();
    std::ofstream myinter("../tempFEM/test/interfacePoints.txt");
    for(int i = 0; i < interfacePoints.size(); i++){
        myinter << interfacePoints.at(i) << endl;
    }
    //3.交界点上的电压
    double *inter_voltage = (double *)calloc(interfacePoints.size(), sizeof(double));
    for(int i = 0;i < interfacePoints.size();++i){
        inter_voltage[i] = mp_2DNode[interfacePoints.at(i)].V;
    }
    //4.第i个节点数组对应第i个分区，保存节点k在分区内的编号
    for(int i = 0;i < m_num_TriEle;++i){
        int k = mp_TriEle[i].n[0];
        int m = mp_TriEle[i].n[1];
        int n = mp_TriEle[i].n[2];
        npart[epartTable[i]][k] = IFPOINT;
        npart[epartTable[i]][m] = IFPOINT;
        npart[epartTable[i]][n] = IFPOINT;
    }
    //5.线单元归入到各自的区域中，全局到局部检索
    int * EdgEle_num_part = (int*)calloc(m_num_part, sizeof(int));
    int **Edgpart = (int**)malloc(m_num_part * sizeof(int*));
    for(int i = 0; i < m_num_part; ++i){
        int order = 0;
        Edgpart[i] = (int*)malloc(m_num_EdgEle * sizeof(int));
        for(int j = 0; j < m_num_EdgEle; j++){
            int n0 = mp_EdgEle[j].n[0];
            int n1 = mp_EdgEle[j].n[1];
            int Node1Part = npart[i][n0];
            int Node2Part = npart[i][n1];
            if((Node1Part == IFPOINT) & (Node2Part == IFPOINT)){
                Edgpart[i][j] = order++;
                ++EdgEle_num_part[i];
            }   else Edgpart[i][j] = 0;
        }
        qDebug() << "Number of EdgElements in partition "<<i<<" is "<<EdgEle_num_part[i];
    }
    //6.分区内节点重新编号，全局到局部检索
    int * freenodepart = (int*)malloc(m_num_part * sizeof(int));
    for(int i = 0;i < m_num_part;++i){
        int order = 0;
        for(int j = 0;j < m_num_pts;++j){
            if(npart[i][j] == IFPOINT){
                npart[i][j] = order++;
            }
        }
        freenodepart[i] = order;
        qDebug()<<"Partition "<<i<<" free nodes: "<<order;
    }

//    std::ofstream mycoutnpart("..\\tempFEM\\test\\npart.txt");
//    for(int i = 0; i < m_num_pts; i++){
//        for(int j = 0; j < m_num_part; j++){
//            mycoutnpart << npart[j][i] << " ";
//        }
//        mycoutnpart << endl;
//    }
    //7.仿照epartTable构造lpartTable，这里只包含第三类边界条件
    int *lpartTable = (int*)malloc(m_num_EdgEle*sizeof(int));
    QVector<int> numbdr(m_num_part);
    for(int i = 0; i < m_num_EdgEle; i++){
        if(mp_EdgEle[i].bdr == 3){
            for(int j = 0; j < m_num_part; j++){
                if(Edgpart[j][i] != 0){
                    ++numbdr[j];
//                    qDebug() << "i = " << i;
                    lpartTable[i] = j;
                    break;
               }
            }
        }
        else lpartTable[i] = -1;
//        if(lpartTable[i] != -1)
//        qDebug() << "lpartTable" << i << " = " << lpartTable[i];
    }

    //8.定义传输线，初始化定义入射过程中节点的电压
    CInterfacePoint ** tl= (CInterfacePoint**) malloc(m_num_part * sizeof(CInterfacePoint*));
    for(int i = 0;i < m_num_part;++i){
        tl[i] = (CInterfacePoint *)malloc(interfacePoints.size() * sizeof(CInterfacePoint));
        for(int j = 0;j < interfacePoints.size();++j){
            //tl[i][j].Vi = pmeshnode[interfacePoints.at(j)].A;
           // fscanf(fp,"%lf \n",&(tl[i][j].Vi));
            tl[i][j].Vi = 0;//tl[i][j].Vi;//pmeshnode[interfacePoints.at(j)].A -
            tl[i][j].Y0 = 0.05;
        }
        //输出每个分区单元数目
        qDebug()<<"Number of TriElements in partition "<<i<<" is "<<TriEle_num_part[i];
    }
//    qDebug() << "partinterfacepoints 0" << " = " << partinterfacepoints[0];


    //PART II:DDTLM迭代
    const int MAX_ITER = 200;
    omp_set_num_threads(m_num_part);
    QVector <umat>locs;
    QVector <mat>vals;
    QVector <int>pos(m_num_part);
    QVector <vec>F;
    QVector <vec>Va;
    QVector <vec>Va_old;
    QVector <mat>partS;

    //构造稀疏矩阵
    for(int part = 0;part < m_num_part;++part){
        umat loc(2, 9*TriEle_num_part[part]+4*numbdr[part]+interfacePoints.size());
        locs.push_back(std::move(loc));
        mat val(1, 9*TriEle_num_part[part]+4*numbdr[part]+interfacePoints.size());
        vals.push_back(std::move(val));
        vec F1 = zeros<vec>(freenodepart[part]);
        F.push_back(std::move(F1));
        vec V1 = zeros<vec>(freenodepart[part]);
        Va.push_back(V1);
        Va_old.push_back(V1);
        pos[part] = 0;
        mat partS1 = zeros<mat>(freenodepart[part],freenodepart[part]);
        partS.push_back(partS1);
    }
    //单元装配到自己的part中
    //三角形单元装配
    const double PI = 3.1415926535;
    double Se, Fe, Sl, Fl;
    for(int k = 0; k < m_num_TriEle; k++){
        int ePart = epartTable[k];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                Se = (PI*mp_TriEle[k].cond*mp_TriEle[k].xavg*(mp_TriEle[k].r[i]*mp_TriEle[k].r[j]+mp_TriEle[k].q[i]*mp_TriEle[k].q[j]))/(2*mp_TriEle[k].Area);
                locs[ePart](0,pos[ePart]) = npart[ePart][mp_TriEle[k].n[i]];
                locs[ePart](1,pos[ePart]) = npart[ePart][mp_TriEle[k].n[j]];
//                qDebug() << "npart[ePart][mp_TriEle[k].n[i]] = " << npart[ePart][mp_TriEle[k].n[i]];
//                qDebug() << "npart[ePart][mp_TriEle[k].n[j]] = " << npart[ePart][mp_TriEle[k].n[j]];
                vals[ePart](0,pos[ePart]) = Se;
                partS[ePart](npart[ePart][mp_TriEle[k].n[i]], npart[ePart][mp_TriEle[k].n[j]]) = partS[ePart](npart[ePart][mp_TriEle[k].n[i]], npart[ePart][mp_TriEle[k].n[j]]) + Se;
                ++pos[ePart];
            }
            Fe = PI*mp_TriEle[k].source*mp_TriEle[k].Area*(mp_TriEle[k].x[i]+3*mp_TriEle[k].xavg)/6;
            F[ePart](npart[ePart][mp_TriEle[k].n[i]]) = F[ePart](npart[ePart][mp_TriEle[k].n[i]]) + Fe;
        }
    }

    //线单元装配
    for(int k = 0; k < m_num_EdgEle; k++){
        int lPart = lpartTable[k];
        if(mp_EdgEle[k].bdr == 3){
//            qDebug() << lPart;
            for(int i = 0; i < 2; i++){
                for(int j = 0; j < 2; j++){
                    if(i == j){
                        Sl = PI*mp_EdgEle[k].h*mp_EdgEle[k].d*(2*mp_EdgEle[k].xavg+2*mp_EdgEle[k].x[i])/6.;
                    }
                    else{
                        Sl = PI*mp_EdgEle[k].h*mp_EdgEle[k].d*(2*mp_EdgEle[k].xavg)/6;
                    }
                    locs[lPart](0,pos[lPart]) = npart[lPart][mp_EdgEle[k].n[i]];
                    locs[lPart](1,pos[lPart]) = npart[lPart][mp_EdgEle[k].n[j]];
                    vals[lPart](0,pos[lPart]) = Sl;
                    partS[lPart](npart[lPart][mp_EdgEle[k].n[i]], npart[lPart][mp_EdgEle[k].n[j]]) = partS[lPart](npart[lPart][mp_EdgEle[k].n[i]], npart[lPart][mp_EdgEle[k].n[j]]) + Sl;
                    ++pos[lPart];
                }
                Fl = PI*mp_EdgEle[k].h*mp_EdgEle[k].Text*mp_EdgEle[k].d*(2*mp_EdgEle[k].xavg+mp_EdgEle[k].x[i])/3;
                F[lPart](npart[lPart][mp_EdgEle[k].n[i]]) = F[lPart](npart[lPart][mp_EdgEle[k].n[i]]) + Fl;
            }
        }
    }
    //查看每个单元的系数矩阵
    std::ofstream mypartS("../tempFEM/test/partS.txt", ios::ate);
    for(int part = 0; part < m_num_part; ++part){
        mypartS << partS[part] << endl << endl;
    }

    std::ofstream mypartF("../tempFEM/test/partF.txt", ios::ate);
    for(int part = 0; part < m_num_part; ++part){
        mypartF << F[part] << endl << endl;
    }

    for(int i = 0; i < m_num_part; i++){
        qDebug() << "pos" << i << " = " << pos[i];
//        qDebug() << "locs" << i << "size = " << 9*TriEle_num_part[i]+4*numbdr[i];
    }
    //求解过程
    int time = 0;
    QVector<int> pos1 = pos;    //记录线单元装配之后的位置
    QVector<vec> F1 = F;
    while(time++ < MAX_ITER ){
        pos = pos1;
        for(int part = 0; part < m_num_part; ++part){
#pragma omp parallel for
            //入射过程求解，每个部分装配上交界点，然后解算
            //叠加每个part的右侧列向量和系数矩阵
            for(int i = 0; i < interfacePoints.size(); i++){
                int n = npart[part][interfacePoints.at(i)];
                if(n != -1){
                    tl[part][i].Vi = inter_voltage[i] - tl[part][i].Vi;
//                    qDebug() << "interfacePoints.at(i) = " << interfacePoints.at(i);
//                    qDebug() <<"n = " << n;
                    double current = 2*tl[part][i].Vi*tl[part][i].Y0;
                    locs[part](0,pos[part]) = n;
                    locs[part](1,pos[part]) = n;
                    vals[part](0,pos[part]) = tl[part][i].Y0;
                    F[part](n) = F1[part](n) + current;
//                    qDebug() << "part = " << part;
//                    qDebug() << "i = " << i;
                }
                else{
                    locs[part](0,pos[part]) = 0;
                    locs[part](1,pos[part]) = 0;
                    vals[part](0,pos[part]) = 0;
                }
                ++pos[part];
            }

            //求解过程
//            SparseSolve(locs[part], vals[part], F[part], freenodepart[part]);
            sp_mat X1(true, locs[part], vals[part], freenodepart[part], freenodepart[part], true, true);
            SuperMatrix sluA;
            NCformat *Astore;
            double   *a;
            int      *asub, *xa;
            int      *perm_c; /* column permutation vector */
            int      *perm_r; /* row permutations from partial pivoting */
            SuperMatrix L;      /* factor L */
            SuperMatrix U;      /* factor U */
            SuperMatrix B;
            int      nrhs, ldx, info, m, n, nnz;
            double   *rhs;
            mem_usage_t   mem_usage;
            superlu_options_t options;
            SuperLUStat_t stat;

            set_default_options(&options);

            /* create matrix A in Harwell-Boeing format.*/
            m = freenodepart[part]; n = freenodepart[part]; nnz = X1.n_nonzero;
            a = const_cast<double *>(X1.values);
            asub = (int*)const_cast<unsigned int*>(X1.row_indices);
            xa = (int*)const_cast<unsigned int*>(X1.col_ptrs);
            dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
            Astore = (NCformat *)sluA.Store;
            printf("Dimension %dx%d; # nonzeros %d\n", sluA.nrow, sluA.ncol, Astore->nnz);

            nrhs = 1;
            if (!(rhs = doubleMalloc(m * nrhs))) ABORT("Malloc fails for rhs[].");
            //将内存拷贝过来
            //memmove(rhs, unknown_b, 5*sizeof(double));
            for (int i = 0; i < m; i++){
                rhs[i] = F[part](i);
            }
            dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

            if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");
            if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");

            /* Initialize the statistics variables. */
            StatInit(&stat);
            dgssv(&options, &sluA, perm_c, perm_r, &L, &U, &B, &stat, &info);
            if (info == 0) {
//                qDebug()<<"Ok.";
                /* This is how you could access the solution matrix. */
                double *sol = (double*)((DNformat*)B.Store)->nzval;
                for(int i = 0; i < freenodepart[part]; ++i){
                    Va[part](i) = sol[i];
//                    qDebug() << "Va[" << part << "](" << i << ") = " << Va[part](i);
                }
            }else {
                qDebug() << "info = " << info;
            }

            SUPERLU_FREE(rhs);
        //    SUPERLU_FREE(xact);
            SUPERLU_FREE(perm_r);
            SUPERLU_FREE(perm_c);
        //    Destroy_CompCol_Matrix(&A);
            Destroy_SuperMatrix_Store(&B);
            Destroy_SuperNode_Matrix(&L);
            Destroy_CompCol_Matrix(&U);
        }
        //反射过程求解
        for(int i = 0; i < interfacePoints.size(); ++i){
            double I = 0;
            double Y = 0;
            for(int j = 0; j < m_num_part; ++j){
                if(npart[j][interfacePoints.at(i)] != -1){
                    Y += tl[j][i].Y0;
                    tl[j][i].Vi = Va[j][npart[j][interfacePoints.at(i)]] - tl[j][i].Vi;
                    I += 2*tl[j][i].Vi*tl[j][i].Y0;
                }
            }
            inter_voltage[i] = I / Y;
        }
        //判断误差
        double outter_error = 1;
        double a = 0, b = 0;
        for(int i = 0; i < m_num_part; ++i){
            for(int j = 0; j < freenodepart[i]; ++j){
                a += (Va_old[i][j] - Va[i][j])*(Va_old[i][j] - Va[i][j]);
                b += Va[i][j] * Va[i][j];
            }
        }
        outter_error = sqrt(a) / sqrt(b);
        qDebug() << "outter_error = " << outter_error;
        if(outter_error < Precision){
            break;
        }
        else{
            for(int i = 0; i < m_num_part; ++i){
                for(int j = 0; j < freenodepart[i]; ++j){
                    Va_old[i][j] = Va[i][j];
                }
            }
        }
    }
    //输出当前的Va
    std::ofstream mypartV("../tempFEM/test/partV.txt", ios::ate);
    for(int part = 0; part < m_num_part; ++part){
        mypartV << Va[part] << endl << endl;
    }
    //整合结果
    std::ofstream mytemp("../tempFEM/test/temp.txt");
    double temp[751];
    for(int part = 0; part < m_num_part; part++){
        for(int i = 0; i < m_num_pts; i++){
            int n = npart[part][i];
            if(n != -1){
                temp[i] = Va[part][n];
            }
        }
    }
    for(int i = 0; i < m_num_pts; i++){
        mytemp << temp[i] << endl;
    }

    return 0;
}
