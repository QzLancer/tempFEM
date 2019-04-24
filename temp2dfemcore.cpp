#include "Temp2DFEMCore.h"
#include <stdio.h>
#include <stdlib.h>
#include <QDebug>
#include <math.h>
#include <iostream>
#include <vector>

CTemp2DFEMCore::CTemp2DFEMCore():mp_2DNode(nullptr),
    mp_VtxEle(nullptr),
    mp_EdgEle(nullptr),
    mp_TriEle(nullptr)
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

int CTemp2DFEMCore::Load2DMeshCOMSOL(const char *fn)
{
    FILE *fp = nullptr;
    int err;
    char ch[256];
    err = fopen_s(&fp, fn, "r");
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
        if(fscanf_s(fp, "%d \n", &mp_TriEle[0].domain) != 1){
            qDebug() << "Error: reading tridomain!";
            return 1;
        }
        else{
            mp_TriEle[0].domain++;
        }
    }
    fclose(fp);
    S = zeros<mat>(m_num_pts, m_num_pts);
    F = zeros<vec>(m_num_pts);
    return 0;
}

int CTemp2DFEMCore::preCalculation()
{
    //计算所有三角形单元的q,r,面积Area,平均半径TriRadius(轴对称模型使用)
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
        mp_TriEle[i].Area = (mp_TriEle[i].q[0]*mp_TriEle[i].r[1]-mp_TriEle[i].q[1]*mp_TriEle[i].r[0])/2;
        mp_TriEle[i].ravg = (mp_TriEle[i].x[0]+mp_TriEle[i].x[1]+mp_TriEle[i].x[2])/3;
    }
    //计算所有线单元的长度d，平均半径EdgRadius(轴对称模型使用)
    for(int i = 0; i < m_num_EdgEle; i++){
        for(int j = 0; j < 2; j++){
            mp_EdgEle[i].x[j] = mp_2DNode[mp_EdgEle[i].n[j]].x;
            mp_EdgEle[i].y[j] = mp_2DNode[mp_EdgEle[i].n[j]].y;
        }
        mp_EdgEle[i].d = std::sqrt(pow(mp_EdgEle[i].x[0]-mp_EdgEle[i].x[1],2)+
                pow(mp_EdgEle[i].y[0]-mp_EdgEle[i].y[1],2));
        mp_EdgEle[i].ravg = (mp_EdgEle[i].x[0]+mp_EdgEle[i].x[1])/2;
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
int CTemp2DFEMCore::setCondition()
{
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
    return 0;
}

using namespace arma;

int CTemp2DFEMCore::StaticAxisAssemble()
{
    const double PI = 3.1415926535;
    //三角形单元装配过程
    mat Se = zeros<mat>(3,3);
    double Fl;
    for(int k = 0; k < m_num_TriEle; k++){
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                Se(i,j) = (PI*mp_TriEle[k].cond*mp_TriEle[k].ravg*(mp_TriEle[k].r[i]*mp_TriEle[k].r[j]+mp_TriEle[k].q[i]*mp_TriEle[k].q[j]))/(2*mp_TriEle[k].Area);
                S(mp_TriEle[k].n[i],mp_TriEle[k].n[j]) = S(mp_TriEle[k].n[i],mp_TriEle[k].n[j]) + Se(i,j);
            }
        }
//        cout << Se << "\n";
    }
    //线单元装配过程
    for(int k = 0; k < m_num_EdgEle; k++){
        if(mp_EdgEle[k].bdr == 2){
            for(int i = 0; i < 2; i++){
                Fl = PI*mp_EdgEle[k].heatflux*mp_EdgEle[k].d*(mp_EdgEle[k].x[i]+2*mp_EdgEle[k].ravg)/3;
                F(mp_EdgEle[k].n[i]) = F(mp_EdgEle[k].n[i])+Fl;
            }
        }
    }
//    cout << F;
    return 0;
}

//通过直接法来求解包含第二类边界条件的温度场问题
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
    mat S1 = zeros<mat>(j, m_num_pts);
    mat S2 = zeros<mat>(j, j);
    vec F1 = zeros<vec>(j);
    vec F2 = zeros<vec>(j);
    vec A1 = zeros<vec>(j);
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
    }
    cout << A;
    return 0;
}

int CTemp2DFEMCore::PostProcessing()
{

    return 0;
}
