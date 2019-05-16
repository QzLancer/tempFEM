#include "temp3dfemcore.h"
#include <math.h>
#include "mpmetis.h"
#include <QVector>
#pragma execution_character_set("utf-8")

using namespace arma;

CTemp3DFEMCore::CTemp3DFEMCore(Widget *parent, const char *fn):
    m_COMSOLMesh(fn),
    mp_3DNode(nullptr),
    mp_VtxEle(nullptr),
    mp_EdgEle(nullptr),
    mp_TriEle(nullptr),
    mp_TetEle(nullptr),
    m_tpartTable(new int(0)),
    m_npartTable(new int(0))
{

}



int CTemp3DFEMCore::Load3DFEMCOMSOL(){
    FILE *fp = nullptr;
    int err;
    char ch[256];
    err = fopen_s(&fp, m_COMSOLMesh, "r");
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
    mp_3DNode = new C3DNode[m_num_pts];
    int pts_ind;//the beginning of the points index

    if (fscanf_s(fp, "%d # lowest mesh point index\n", &pts_ind) != 1) {
        qDebug() << "Error: reading pts_ind!";
        return 1;
    }
    fgets(ch, 256, fp);

    for (int i = pts_ind; i < m_num_pts; i++) {
        //读取x,y坐标
        if (fscanf_s(fp, "%lf %lf %lf\n", &(mp_3DNode[i].x), &(mp_3DNode[i].y), &(mp_3DNode[i].z)) != 3) {
            qDebug() << "Error: reading mesh point!";
            return 1;
        }
        //                else{
        //                    qDebug() << mp_3DNode[i].x << " " << mp_3DNode[i].y << " " << mp_3DNode[i].z << "\n";
        //                }
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
            mp_TriEle[i].domain++;
        }
    }
    //----------------tetnode-----------------------------------
    for (int i = 0; i < 6; i++){
        fgets(ch, 256, fp);
    }
    if(fscanf_s(fp, "%d # number of elements\n", &m_num_TetEle) != 1){
        qDebug() << "Error: reading m_num_TriEle!";
        return 1;
    }
    else qDebug() << m_num_TetEle << "number of Tetrahedron elements.";
    fgets(ch, 256, fp);
    mp_TetEle = new CTetElement[m_num_TetEle];
    for (int i = 0; i < m_num_TetEle; i++){
        if (fscanf_s(fp, "%d %d %d %d\n", &mp_TetEle[i].n[0], &mp_TetEle[i].n[1], &mp_TetEle[i].n[2], &mp_TetEle[i].n[3]) != 4) {
            qDebug() << "Error: reading Tet elements points!";
            return 1;
        }
    }
    //----------------tetdomain-----------------------------------

    for (int i = 0; i < 2; i++){
        fgets(ch, 256, fp);
    }
    for (int i = 0; i < m_num_TetEle; i++){
        if(fscanf_s(fp, "%d \n", &mp_TetEle[i].domain) != 1){
            qDebug() << "Error: reading tetdomain!";
            return 1;
        }
        else{
            //            mp_TetEle[i].domain++;
        }
    }

    fclose(fp);

    return 0;
}

int CTemp3DFEMCore::preCalculation(){
    //I:计算所有四面体单元中的p,q,r,s,volume
    mat p0 = zeros<mat>(3,3);
    //    mat q0 = zeros<mat>(3,3);
    //    mat r0 = zeros<mat>(3,3);
    //    mat s0 = zeros<mat>(3,3);
    mat v0 = zeros<mat>(4,4);
    for(int i = 0; i < m_num_TetEle; ++i){
        for(int j = 0; j < 4; ++j){
            mp_TetEle[i].x[j] = mp_3DNode[mp_TetEle[i].n[j]].x;
            mp_TetEle[i].y[j] = mp_3DNode[mp_TetEle[i].n[j]].y;
            mp_TetEle[i].z[j] = mp_3DNode[mp_TetEle[i].n[j]].z;
        }
        //求解p
        p0 = {{mp_TetEle[i].x[1], mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].p[0] = det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].p[1] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[1], mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].p[2] = det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[1], mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], mp_TetEle[i].z[2]}};
        mp_TetEle[i].p[3] = -det(p0);

        //求解q
        p0 = {{1, mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {1, mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {1, mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].q[0] = -det(p0);
        p0 = {{1, mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {1, mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {1, mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].q[1] = det(p0);
        p0 = {{1, mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {1, mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {1, mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].q[2] = -det(p0);
        p0 = {{1, mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {1, mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {1, mp_TetEle[i].y[2], mp_TetEle[i].z[2]}};
        mp_TetEle[i].q[3] = det(p0);

        //求解r
        p0 = {{mp_TetEle[i].x[1], 1, mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[2], 1, mp_TetEle[i].z[2]},
              {mp_TetEle[i].x[3], 1, mp_TetEle[i].z[3]}};
        mp_TetEle[i].r[0] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], 1, mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[2], 1, mp_TetEle[i].z[2]},
              {mp_TetEle[i].x[3], 1, mp_TetEle[i].z[3]}};
        mp_TetEle[i].r[1] = det(p0);
        p0 = {{mp_TetEle[i].x[0], 1, mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[1], 1, mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[3], 1, mp_TetEle[i].z[3]}};
        mp_TetEle[i].r[2] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], 1, mp_TetEle[i].z[0]},
              {mp_TetEle[i].x[1], 1, mp_TetEle[i].z[1]},
              {mp_TetEle[i].x[2], 1, mp_TetEle[i].z[2]}};
        mp_TetEle[i].r[3] = det(p0);

        //求解s
        p0 = {{mp_TetEle[i].x[1], mp_TetEle[i].y[1], 1},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], 1},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], 1}};
        mp_TetEle[i].s[0] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], 1},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], 1},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], 1}};
        mp_TetEle[i].s[1] = det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], 1},
              {mp_TetEle[i].x[1], mp_TetEle[i].y[1], 1},
              {mp_TetEle[i].x[3], mp_TetEle[i].y[3], 1}};
        mp_TetEle[i].s[2] = -det(p0);
        p0 = {{mp_TetEle[i].x[0], mp_TetEle[i].y[0], 1},
              {mp_TetEle[i].x[1], mp_TetEle[i].y[1], 1},
              {mp_TetEle[i].x[2], mp_TetEle[i].y[2], 1}};
        mp_TetEle[i].s[3] = det(p0);

        //求解volume
        v0 = {{1, mp_TetEle[i].x[0], mp_TetEle[i].y[0], mp_TetEle[i].z[0]},
              {1, mp_TetEle[i].x[1], mp_TetEle[i].y[1], mp_TetEle[i].z[1]},
              {1, mp_TetEle[i].x[2], mp_TetEle[i].y[2], mp_TetEle[i].z[2]},
              {1, mp_TetEle[i].x[3], mp_TetEle[i].y[3], mp_TetEle[i].z[3]}};
        mp_TetEle[i].Volume = det(v0)/6;
    }

    //B:计算所有三角形单元的Area
    for(int i = 0; i < m_num_TriEle; ++i){
        for(int j = 0; j < 3; ++j){
            mp_TriEle[i].x[j] = mp_3DNode[mp_TriEle[i].n[j]].x;
            mp_TriEle[i].y[j] = mp_3DNode[mp_TriEle[i].n[j]].y;
            mp_TriEle[i].z[j] = mp_3DNode[mp_TriEle[i].n[j]].z;
            //                qDebug() << "mp_TriEle " << i << " x " << j << " = " << mp_TriEle[i].x[j];
            //                qDebug() << "mp_TriEle " << i << " y " << j << " = " << mp_TriEle[i].y[j];
            //                qDebug() << "mp_TriEle " << i << " z " << j << " = " << mp_TriEle[i].z[j];
        }
        vec ab = zeros<vec>(3);
        vec ac = zeros<vec>(3);
        vec s0 = zeros<vec>(3);
        ab = {mp_TriEle[i].x[1]-mp_TriEle[i].x[0], mp_TriEle[i].y[1]-mp_TriEle[i].y[0], mp_TriEle[i].z[1]-mp_TriEle[i].z[0]};
        ac = {mp_TriEle[i].x[2]-mp_TriEle[i].x[0], mp_TriEle[i].y[2]-mp_TriEle[i].y[0], mp_TriEle[i].z[2]-mp_TriEle[i].z[0]};
        s0 = cross(ab, ac);
        mp_TriEle[i].Area = sqrt(s0(0)*s0(0) + s0(1)*s0(1) + s0(2)*s0(2))/2;
        //            qDebug() << "mp_TriEle" << i << "area = " << mp_TriEle[i].Area;
        //           ab.print("ab = ");
        //           ac.print("ac = ");
    }

    return 0;
}

//不同模型改写这一部分
int CTemp3DFEMCore::setCondition()
{
    //热源设置
    //        std::ofstream mytetsource("../tempFEM/test/tetsource.txt");
    for(int i = 0; i < m_num_TetEle; ++i){
        if((mp_TetEle[i].domain == 7) | (mp_TetEle[i].domain == 9)){
            mp_TetEle[i].source = 500000;
        }else{
            mp_TetEle[i].source = 0;
        }
        //         mytetsource << "mp_TetEle " << i << " domain = " << mp_TetEle[i].source << endl;
    }
    //热导率设置
    for(int i = 0; i < m_num_TetEle; ++i){
        if((mp_TetEle[i].domain == 7) | (mp_TetEle[i].domain == 9)){
            mp_TetEle[i].cond = 400;
            mp_TetEle[i].Material = 0;
            mp_TetEle[i].LinerFlag = 1;
        }
        else if((mp_TetEle[i].domain == 2) | (mp_TetEle[i].domain == 4) | (mp_TetEle[i].domain == 10) | (mp_TetEle[i].domain == 12)| (mp_TetEle[i].domain == 13)){
            mp_TetEle[i].cond = 76.2;
            mp_TetEle[i].Material = 1;
            mp_TetEle[i].LinerFlag = 1;
        }
        else if((mp_TetEle[i].domain == 1) | (mp_TetEle[i].domain == 6) | (mp_TetEle[i].domain == 8)){
            mp_TetEle[i].cond = 0.26;
            mp_TetEle[i].Material = 2;
            mp_TetEle[i].LinerFlag = 1;
        }
        else if((mp_TetEle[i].domain == 3) | (mp_TetEle[i].domain == 5) | (mp_TetEle[i].domain == 11)){
            mp_TetEle[i].cond = 0.03;
            mp_TetEle[i].Material = 3;
            mp_TetEle[i].LinerFlag = 1;
        }
    }
    //第三类边界条件设置
    for(int i = 0; i < m_num_TriEle; ++i){
        if((mp_TriEle[i].domain == 1) | (mp_TriEle[i].domain == 2) | (mp_TriEle[i].domain == 3) | (mp_TriEle[i].domain == 4) |
                (mp_TriEle[i].domain == 5) | (mp_TriEle[i].domain == 6) | (mp_TriEle[i].domain == 91) | (mp_TriEle[i].domain == 92) |
                (mp_TriEle[i].domain == 93) | (mp_TriEle[i].domain == 137) | (mp_TriEle[i].domain == 143) | (mp_TriEle[i].domain == 183)){
            mp_TriEle[i].bdr = 3;
            mp_TriEle[i].h = 20;
            mp_TriEle[i].Text = 293.15;
        }
    }

    return 0;
}

int CTemp3DFEMCore::Static3DAssemble()
{
    const double PI = 3.14159265358979323846;
    int pos = 0;
    int numbdr = 0;
    for(int i = 0; i < m_num_TriEle; ++i){
        if(mp_TriEle[i].bdr == 3){
            //            qDebug() << "i = " << i;
            ++numbdr;
        }
    }
    umat locs(2, 16*m_num_TetEle+9*numbdr);
    mat vals(1, 16*m_num_TetEle+9*numbdr);
    vec F = zeros<vec>(m_num_pts);
    double St;
    double Ft;
    double Se;
    double Fe;
//        mat S = zeros<mat>(m_num_pts, m_num_pts);
    //四面体单元装配
    //    std::ofstream myFt("../tempFEM/test/Ft.txt");
    for(int k = 0; k < m_num_TetEle; ++k){
        for(int i = 0; i < 4; ++i){
            for(int j = 0; j < 4; ++j){
                St = mp_TetEle[k].cond*(mp_TetEle[k].q[i]*mp_TetEle[k].q[j]+mp_TetEle[k].r[i]*mp_TetEle[k].r[j]+mp_TetEle[k].s[i]*mp_TetEle[k].s[j])/(36*mp_TetEle[k].Volume);
                locs(0, pos) = mp_TetEle[k].n[i];
                locs(1, pos) = mp_TetEle[k].n[j];
                vals(0, pos) = St;
//                                S(mp_TetEle[k].n[i], mp_TetEle[k].n[j]) = S(mp_TetEle[k].n[i], mp_TetEle[k].n[j]) + St;
                ++pos;
            }
            Ft = mp_TetEle[k].source*mp_TetEle[k].Volume/4;
            F(mp_TetEle[k].n[i]) = F(mp_TetEle[k].n[i]) + Ft;
        }
        //        myFt << "Ft " << k << " = " << Ft << endl;
    }

    //三角形单元装配
    //    std::ofstream myFe("../tempFEM/test/Fe.txt");
    for(int k = 0; k < m_num_TriEle; ++k){
        for(int i = 0; i < 3; ++i){
            if(mp_TriEle[k].bdr == 3){
                for(int j = 0; j < 3; ++j){
                    if(i == j){
                        Se = mp_TriEle[k].h*mp_TriEle[k].Area/6;
                    }else{
                        Se = mp_TriEle[k].h*mp_TriEle[k].Area/12;
                    }
//                                        S(mp_TetEle[k].n[i], mp_TetEle[k].n[j]) = S(mp_TetEle[k].n[i], mp_TetEle[k].n[j]) + Se;
                    locs(0, pos) = mp_TriEle[k].n[i];
                    locs(1, pos) = mp_TriEle[k].n[j];
                    //                    if(Se == 0){
                    //                        qDebug() << "k = " << k << ", i = " << i << ", j = " << j << ", pos = " << pos;
                    //                          qDebug() << "mp_TriEle[k].h = " << mp_TriEle[k].h;
                    //                          qDebug() << "mp_TriEle[k].Area = " << mp_TriEle[k].Area;
                    //                    }
                    vals(0, pos) = Se;
                    ++pos;
                    //                    qDebug() << "k = " << k;
                    //                    qDebug() << "i = " << i;
                    //                    qDebug() << "j = " << j;
                    //                    qDebug() << "pos = " << pos;
                }
                Fe = mp_TriEle[k].h*mp_TriEle[k].Text*mp_TriEle[k].Area/3;
                F(mp_TriEle[k].n[i]) = F(mp_TriEle[k].n[i]) + Fe;
                //                qDebug()<< mp_TriEle[k].Area;
            }
            else if(mp_TriEle[k].bdr == 2){
                Fe = mp_TriEle[k].heatflux*mp_TriEle[k].Area/3;
                F(mp_TriEle[k].n[i]) = F(mp_TriEle[k].n[i]) + Fe;
            }else Fe = 0;
        }
        //        myFe << "Fe " << k << " = " << Fe << endl;
    }

    //    std::ofstream mylocs("../tempFEM/test/locs.txt");
    //    mylocs << locs.t();
    //    std::ofstream myvals("../tempFEM/test/vals.txt");
    //    myvals << vals.t();
        std::ofstream myF3D("../tempFEM/test/F3D.txt");
        myF3D << F;
//        std::ofstream myS3D("../tempFEM/test/S3D.txt");
//        myS3D << S;

    sp_mat X(true, locs, vals, m_num_pts, m_num_pts, true, true);
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
    m = m_num_pts; n = m_num_pts; nnz = X.n_nonzero;
    a = const_cast<double *>(X.values);
    asub = (int*)const_cast<unsigned int*>(X.row_indices);
    xa = (int*)const_cast<unsigned int*>(X.col_ptrs);
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
        std::ofstream myTemp3D("../tempFEM/test/Temp3D.txt");
        /* This is how you could access the solution matrix. */
        double *sol = (double*)((DNformat*)B.Store)->nzval;
        for(int i = 0; i < m_num_pts; ++i){
            mp_3DNode[i].V = sol[i];
            myTemp3D << mp_3DNode[i].V << endl;
        }
        myTemp3D.close();
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

int CTemp3DFEMCore::GenerateMetisMesh(int partition)
{
    m_num_Part = partition;
    if(m_COMSOLMesh == nullptr){
        return 1;
    }
    FILE *fp = nullptr;

    //生成metis的输入分网文件
    strcpy(m_METISMesh,m_COMSOLMesh);
    strcat(m_METISMesh, ".metis");
    fp = fopen(m_METISMesh, "w+");
    if(!fp){
        return 1;
    }
    //写入单元数目和单元节点编号（从1开始）
    fprintf(fp, "%d\n", m_num_TetEle);
    for(int i = 0; i < m_num_TetEle; i++){
        fprintf(fp, "%d %d %d %d\n", mp_TetEle[i].n[0]+1, mp_TetEle[i].n[1]+1, mp_TetEle[i].n[2]+1, mp_TetEle[i].n[3]+1);
    }
    fclose(fp);
    //调用metis进行分区
    char part[3];
    sprintf(part, "%d", partition);
    char a[] = "QzLancer";
    char * myargv[] = {a, m_METISMesh, part};
    mpmetis(3,myargv);
    //读入单元分区表
    char epartname[256];
    sprintf(epartname,"%s.epart.%d",m_METISMesh,m_num_Part);
    fp = fopen(epartname,"r");
    if(!fp){
        qDebug()<<"Error in open epart file!";
        return 1;
    }
    if(m_tpartTable){
        free(m_tpartTable);//释放之前申请的空间
    }
    m_tpartTable = (int*)malloc(m_num_TetEle*sizeof(int));
    for(int i=0;i < m_num_TetEle;++i){
        fscanf(fp,"%d\n",m_tpartTable+i);
    }
    fclose(fp);
    //读入节点分区表
    char npartname[256];
    sprintf(npartname,"%s.npart.%d",m_METISMesh,partition);
    fp = fopen(npartname,"r");
    if(!fp){
        qDebug()<<"Error in open npart file!";
        return 1;
    }
    if(m_npartTable){
        free(m_npartTable);//释放之前申请的空间
    }
    m_npartTable = (int*)malloc(m_num_pts*sizeof(int));
    for(int i=0;i < m_num_pts;++i){
        fscanf(fp,"%d\n",m_npartTable+i);
    }
    fclose(fp);

    return 0;
}

int CTemp3DFEMCore::DDTLM3DSolve()
{
    //PART I:处理分区信息
    //A :查找各个分区上的交界点，并且去除掉边界点
    //1.创建partition个大小为num_pts的数组，初始化为-1
    //-1表示该位置的节点不在该分区内，方便交界点过来询问
    int * num_Tet_part = (int*)calloc(m_num_Part, sizeof(int));
    int ** npart = (int**)malloc(m_num_Part * sizeof(int*));
    for(int i = 0; i < m_num_Part; ++i){
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
    for(int i = 0; i < m_num_TetEle; ++i){
        if((m_tpartTable[i] != m_npartTable[mp_TetEle[i].n[0]])){
            m_npartTable[mp_TetEle[i].n[0]] = IFPOINT;
        }
        if((m_tpartTable[i] != m_npartTable[mp_TetEle[i].n[1]])){
            m_npartTable[mp_TetEle[i].n[1]] = IFPOINT;
        }
        if((m_tpartTable[i] != m_npartTable[mp_TetEle[i].n[2]])){
            m_npartTable[mp_TetEle[i].n[2]] = IFPOINT;
        }
        if((m_tpartTable[i] != m_npartTable[mp_TetEle[i].n[3]])){
            m_npartTable[mp_TetEle[i].n[3]] = IFPOINT;
        }
        //计算每个分区内的单元数目
        num_Tet_part[m_tpartTable[i]]++;
    }
    for(int i = 0;i < m_num_pts;++i){
        if(m_npartTable[i] == IFPOINT){
            interfacePoints.push_back(i);
            //qDebug()<<i;
        }
    }
    qDebug()<<"interface points size: "<<interfacePoints.size();
    //3.定义传输线，有一部分是没有用到的
    CInterfacePoint ** tl= (CInterfacePoint**) malloc(m_num_Part * sizeof(CInterfacePoint*));
    for(int i = 0;i < m_num_Part;++i){
        tl[i] = (CInterfacePoint *)malloc(interfacePoints.size() * sizeof(CInterfacePoint));
        for(int j = 0;j < interfacePoints.size();++j){
            tl[i][j].Vi = 0;
            tl[i][j].Y0 = 1;
        }
        //输出每个分区单元数目
        qDebug()<<"Number of tet elements in partition "<<i<<" is "<<num_Tet_part[i];
    }
    //4.交界点上的电压
    double *inter_voltage = (double *)calloc(interfacePoints.size(), sizeof(double));
    //5.第i个节点数组对应第i个分区，保存节点k在分区内的编号
    for(int i = 0; i < m_num_TetEle; ++i){
        for(int j = 0; j < 4; ++j){
            npart[m_tpartTable[i]][mp_TetEle[i].n[j]] = IFPOINT;
            //            qDebug() << "npart row: " << mp_TetEle[i].n[j] << " col: " << j << " = " << m_tpartTable[i];
        }
    }
    //6.分区内节点重新编号
    int * freenodepart = (int*)malloc(m_num_Part * sizeof(int));
    for(int i = 0;i < m_num_Part;++i){
        int order = 0;
        for(int j = 0; j < m_num_pts; ++j){
            if(npart[i][j] != -1){
                npart[i][j] = order++;
            }
        }
        freenodepart[i] = order;
        qDebug()<<"Partition "<<i<<" free nodes: "<<order;
    }
    //7.为自由节点值分配空间
    double **Va = (double **)malloc(m_num_Part * sizeof(double*));
    double **Va_old = (double **)malloc(m_num_Part * sizeof(double*));
    for(int i = 0;i < m_num_Part;++i){
        Va[i] = (double *)calloc(freenodepart[i],  sizeof(double));
        Va_old[i] = (double *)calloc(freenodepart[i], sizeof(double));
    }
    //8.将三角形单元归入到各自的区域中
    int *num_Tri_part = (int*)calloc(m_num_Part, sizeof(int));
    int *epartTable = (int*)calloc(m_num_TriEle, sizeof(int));
    for(int i = 0; i < m_num_TriEle; ++i){
        int n0 = mp_TriEle[i].n[0];
        int n1 = mp_TriEle[i].n[1];
        int n2 = mp_TriEle[i].n[2];
        for(int part = 0; part < m_num_Part; ++part){
            if((npart[part][n0] != -1) & (npart[part][n1] != -1) & (npart[part][n2] != -1)){
                epartTable[i] = part;
                ++num_Tri_part[part];
                //                qDebug() << "epartTable " << i << " = " << epartTable[i];
                break;
            }
        }
    }
    for(int part = 0; part < m_num_Part; ++part){
        qDebug()<<"Number of tet elements in partition "<<part<<" is "<<num_Tri_part[part];
    }
    //    //B: 计算迭代中的一些不变量
    //    CTetResisMatrix *TetRM = (CTetResisMatrix*)malloc(m_num_TetEle * sizeof(CTetResisMatrix));
    //    CTriResistMarix *TriRM = (CTriResistMarix*)malloc(m_num_TetEle * sizeof(CTriResistMarix));
    //    //Tet单元上三角矩阵计算
    //    for(int i = 0; i < m_num_TetEle; ++i){
    //        TetRM[i].Y11 = mp_TetEle[i].q[0]*mp_TetEle[i].q[0] + mp_TetEle[i].r[0]*mp_TetEle[i].r[0] + mp_TetEle[i].s[0]*mp_TetEle[i].s[0];
    //        TetRM[i].Y12 = mp_TetEle[i].q[0]*mp_TetEle[i].q[1] + mp_TetEle[i].r[0]*mp_TetEle[i].r[1] + mp_TetEle[i].s[0]*mp_TetEle[i].s[1];
    //        TetRM[i].Y13 = mp_TetEle[i].q[0]*mp_TetEle[i].q[2] + mp_TetEle[i].r[0]*mp_TetEle[i].r[2] + mp_TetEle[i].s[0]*mp_TetEle[i].s[2];
    //        TetRM[i].Y14 = mp_TetEle[i].q[0]*mp_TetEle[i].q[3] + mp_TetEle[i].r[0]*mp_TetEle[i].r[3] + mp_TetEle[i].s[0]*mp_TetEle[i].s[3];
    //        TetRM[i].Y22 = mp_TetEle[i].q[1]*mp_TetEle[i].q[1] + mp_TetEle[i].r[1]*mp_TetEle[i].r[1] + mp_TetEle[i].s[1]*mp_TetEle[i].s[1];
    //        TetRM[i].Y23 = mp_TetEle[i].q[1]*mp_TetEle[i].q[2] + mp_TetEle[i].r[1]*mp_TetEle[i].r[2] + mp_TetEle[i].s[1]*mp_TetEle[i].s[2];
    //        TetRM[i].Y24 = mp_TetEle[i].q[1]*mp_TetEle[i].q[3] + mp_TetEle[i].r[1]*mp_TetEle[i].r[3] + mp_TetEle[i].s[1]*mp_TetEle[i].s[3];
    //        TetRM[i].Y33 = mp_TetEle[i].q[2]*mp_TetEle[i].q[2] + mp_TetEle[i].r[2]*mp_TetEle[i].r[2] + mp_TetEle[i].s[2]*mp_TetEle[i].s[2];
    //        TetRM[i].Y34 = mp_TetEle[i].q[2]*mp_TetEle[i].q[3] + mp_TetEle[i].r[2]*mp_TetEle[i].r[3] + mp_TetEle[i].s[2]*mp_TetEle[i].s[3];
    //        TetRM[i].Y44 = mp_TetEle[i].q[3]*mp_TetEle[i].q[3] + mp_TetEle[i].r[3]*mp_TetEle[i].r[3] + mp_TetEle[i].s[3]*mp_TetEle[i].s[3];

    //        TetRM[i].Y11 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y12 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y13 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y14 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y22 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y23 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y24 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y33 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y34 /= 36*mp_TetEle[i].Volume;
    //        TetRM[i].Y44 /= 36*mp_TetEle[i].Volume;

    //    }

    //    //Tri单元上三角矩阵计算
    //    for(int i = 0; i < m_num_TriEle; ++i){
    //        TriRM[i].Y11 = mp_TriEle[i].Area*2;
    //        TriRM[i].Y12 = mp_TriEle[i].Area;
    //        TriRM[i].Y13 = mp_TriEle[i].Area;
    //        TriRM[i].Y22 = mp_TriEle[i].Area*2;
    //        TriRM[i].Y23 = mp_TriEle[i].Area;
    //        TriRM[i].Y33 = mp_TriEle[i].Area*2;

    //        TriRM[i].Y11 /= 12;
    //        TriRM[i].Y12 /= 12;
    //        TriRM[i].Y13 /= 12;
    //        TriRM[i].Y22 /= 12;
    //        TriRM[i].Y23 /= 12;
    //        TriRM[i].Y33 /= 12;
    //    }
    //输出测试
    std::ofstream mynpart("../tempFEM/test/npart.txt");
    for(int i = 0; i < m_num_pts; ++i){
        for(int part = 0; part < m_num_Part; ++part){
            mynpart << npart[part][i] << " ";
        }
        mynpart << endl;
    }
    mynpart.close();

    std::ofstream myepartTable("../tempFEM/test/epartTable.txt");
    for(int i = 0; i < m_num_TriEle; ++i){
        myepartTable << epartTable[i] << endl;
    }
    myepartTable.close();

    //PART II:DDTLM迭代
    int MAX_ITER = 200;
    omp_set_num_threads(m_num_Part);
    for(int iter = 0; iter < MAX_ITER; ++iter){
#pragma omp parallel for
        for(int part = 0; part < m_num_Part; ++part){
            double St;
            double Se;
            double Ft;
            double Fe;
            vec F = zeros<vec>(m_num_pts);
            int pos = 0;
            umat locs(2, 16*m_num_TetEle+9*m_num_TriEle+interfacePoints.size());
            locs.zeros();
            mat vals(1, 16*m_num_TetEle+9*m_num_TriEle+interfacePoints.size());
            //每个区域的四面体单元装配
            for(int k = 0; k < m_num_TetEle; ++k){
                int tetpart = m_tpartTable[k];
                if(tetpart == part){
                    for(int i = 0; i < 4; ++i){
                        int n0 = mp_TetEle[k].n[i];
                        for(int j = 0; j < 4; ++j){
                            St = mp_TetEle[k].cond*(mp_TetEle[k].q[i]*mp_TetEle[k].q[j]+mp_TetEle[k].r[i]*mp_TetEle[k].r[j]+mp_TetEle[k].s[i]*mp_TetEle[k].s[j])/(36*mp_TetEle[k].Volume);
                            int n1 = mp_TetEle[k].n[j];
                            locs(0, pos) = npart[part][n0];
                            locs(1, pos) = npart[part][n1];
                            vals(0, pos) = St;
                            ++pos;
                        }
                        Ft = mp_TetEle[k].source*mp_TetEle[k].Volume/4;
                        F(npart[part][n0]) = F(npart[part][n0]) + Ft;
                    }
                }
            }
//            qDebug() << "pos0 " << part << " = " << pos;
            //每个区域的三角形单元装配
            for(int k = 0; k < m_num_TriEle; ++k){
                int tripart = epartTable[k];
                if(tripart == part){
                    for(int i = 0; i < 3; ++i){
                        int n0 = mp_TriEle[k].n[i];
                        for(int j = 0; j < 3; ++j){
                            int n1 = mp_TriEle[k].n[j];
                            if(i == j){
                                Se = mp_TriEle[k].h*mp_TriEle[k].Area/6;
                            }else{
                                Se = mp_TriEle[k].h*mp_TriEle[k].Area/12;
                            }
                            locs(0, pos) = npart[part][n0];
                            locs(1, pos) = npart[part][n1];
                            vals(0, pos) = Se;
                            ++pos;
                        }
                        Fe = mp_TriEle[k].h*mp_TriEle[k].Text*mp_TriEle[k].Area/3;
                        //qDebug() << "k = " << k << ", i = " << i << endl;
                        F(npart[part][n0]) = F(npart[part][n0]) + Fe;
                    }
                }
            }
            //入射装配
            for(int i = 0; i < interfacePoints.size(); ++i){
                int n = npart[part][interfacePoints.at(i)];
                if(n != -1){
                    tl[part][i].Vi = inter_voltage[i] - tl[part][i].Vi;
                    double current = 2*tl[part][i].Vi*tl[part][i].Y0;
                    locs(0, pos) = n;
                    locs(1, pos) = n;
                    vals(0, pos) = tl[part][i].Y0;
                    F[n] = F[n] + current;
                }
                else{
                    locs(0, pos) = 0;
                    locs(0, pos) = 0;
                    vals(0, pos) = 0;
                }
                ++pos;
            }
            //入射过程求解
            sp_mat X1(true, locs, vals, freenodepart[part], freenodepart[part], true, true);
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
                rhs[i] = F(i);
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
                    Va[part][i] = sol[i];
//                    qDebug() << "Va[" << part << "](" << i << ") = " << Va[part][i];
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
            for(int j = 0; j < m_num_Part; ++j){
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
        for(int i = 0; i < m_num_Part; ++i){
            for(int j = 0; j < freenodepart[i]; ++j){
                a += (Va_old[i][j] - Va[i][j])*(Va_old[i][j] - Va[i][j]);
                b += Va[i][j] * Va[i][j];
            }
        }
        outter_error = sqrt(a) / sqrt(b);
        qDebug() << "outter_error = " << outter_error;
        if(outter_error < Precision){
            for(int part = 0; part < m_num_Part; ++part){
                for(int i = 0; i < m_num_pts; ++i){
                    int n = npart[part][i];
                    if(n != -1){
                        mp_3DNode[i].V = Va[part][n];
                    }
                }
            }
            break;
        }
        else{
            for(int i = 0; i < m_num_Part; ++i){
                for(int j = 0; j < freenodepart[i]; ++j){
                    Va_old[i][j] = Va[i][j];
                }
            }
        }
    }

    std::ofstream myTemp3DDDTLM("../tempFEM/test/Temp3DDDTLM.txt");
    for(int i = 0; i < m_num_pts; ++i){
        myTemp3DDDTLM << mp_3DNode[i].V << endl;
    }
    //求解结束，回收空间
    free(num_Tet_part);
    for(int i = 0; i < m_num_Part; ++i){
        free(npart[i]);
    }
    free(npart);
    for(int i = 0; i < m_num_Part; ++i){
        free(tl[i]);
    }
    free(tl);
    free(inter_voltage);
    free(freenodepart);
    for(int i = 0; i < m_num_Part; ++i){
        free(Va[i]);
        free(Va_old[i]);
    }
    free(Va);
    free(Va_old);
    free(num_Tri_part);
    free(epartTable);
    //    free(TetRM);
    //    free(TriRM);
    return 0;
}

int CTemp3DFEMCore::DDTLM3DSolve1()
{
    //PART I:处理分区信息
    //A :查找各个分区上的交界点，并且去除掉边界点
    //1.创建partition个大小为num_pts的数组，初始化为-1
    //-1表示该位置的节点不在该分区内，方便交界点过来询问
    int * num_Tet_part = (int*)calloc(m_num_Part, sizeof(int));
    int ** npart = (int**)malloc(m_num_Part * sizeof(int*));
    for(int i = 0; i < m_num_Part; ++i){
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
    for(int i = 0; i < m_num_TetEle; ++i){
        if((m_tpartTable[i] != m_npartTable[mp_TetEle[i].n[0]])){
            m_npartTable[mp_TetEle[i].n[0]] = IFPOINT;
        }
        if((m_tpartTable[i] != m_npartTable[mp_TetEle[i].n[1]])){
            m_npartTable[mp_TetEle[i].n[1]] = IFPOINT;
        }
        if((m_tpartTable[i] != m_npartTable[mp_TetEle[i].n[2]])){
            m_npartTable[mp_TetEle[i].n[2]] = IFPOINT;
        }
        if((m_tpartTable[i] != m_npartTable[mp_TetEle[i].n[3]])){
            m_npartTable[mp_TetEle[i].n[3]] = IFPOINT;
        }
        //计算每个分区内的单元数目
        num_Tet_part[m_tpartTable[i]]++;
    }
    for(int i = 0;i < m_num_pts;++i){
        if(m_npartTable[i] == IFPOINT){
            interfacePoints.push_back(i);
            //qDebug()<<i;
        }
    }
    qDebug()<<"interface points size: "<<interfacePoints.size();
    std::ofstream myinterfacepoint("../tempFEM/test/interfacepoint.txt");
    for(int i = 0; i < interfacePoints.size(); ++i){
        myinterfacepoint << interfacePoints.at(i) << endl;
    }
    myinterfacepoint.close();
    //3.定义传输线，有一部分是没有用到的
    CInterfacePoint ** tl= (CInterfacePoint**) malloc(m_num_Part * sizeof(CInterfacePoint*));
    for(int i = 0;i < m_num_Part;++i){
        tl[i] = (CInterfacePoint *)malloc(interfacePoints.size() * sizeof(CInterfacePoint));
        for(int j = 0;j < interfacePoints.size();++j){
            tl[i][j].Vi = 0;
            tl[i][j].Y0 = 0.01;
        }
        //输出每个分区单元数目
        qDebug()<<"Number of tet elements in partition "<<i<<" is "<<num_Tet_part[i];
    }
    //4.交界点上的电压
    double *inter_voltage = (double *)calloc(interfacePoints.size(), sizeof(double));
    //5.第i个节点数组对应第i个分区，保存节点k在分区内的编号
    for(int i = 0; i < m_num_TetEle; ++i){
        for(int j = 0; j < 4; ++j){
            npart[m_tpartTable[i]][mp_TetEle[i].n[j]] = IFPOINT;
            //            qDebug() << "npart row: " << mp_TetEle[i].n[j] << " col: " << j << " = " << m_tpartTable[i];
        }
    }
    //6.分区内节点重新编号
    int * freenodepart = (int*)malloc(m_num_Part * sizeof(int));
    for(int i = 0;i < m_num_Part;++i){
        int order = 0;
        for(int j = 0; j < m_num_pts; ++j){
            if(npart[i][j] != -1){
                npart[i][j] = order++;
            }
        }
        freenodepart[i] = order;
        qDebug()<<"Partition "<<i<<" free nodes: "<<order;
    }
    std::ofstream mynpart("../tempFEM/test/npart.txt");
    for(int i = 0; i < m_num_pts; ++i){
        for(int part = 0; part < m_num_Part; ++part){
            mynpart << npart[part][i] << " ";
        }
        mynpart << endl;
    }
    mynpart.close();
//    //7.为自由节点值分配空间
//    double **Va = (double **)malloc(m_num_Part * sizeof(double*));
//    double **Va_old = (double **)malloc(m_num_Part * sizeof(double*));
//    for(int i = 0;i < m_num_Part;++i){
//        Va[i] = (double *)calloc(freenodepart[i],  sizeof(double));
//        Va_old[i] = (double *)calloc(freenodepart[i], sizeof(double));
//    }
    //8.将三角形单元归入到各自的区域中
    int *num_Tri_part = (int*)calloc(m_num_Part, sizeof(int));
    int *epartTable = (int*)calloc(m_num_TriEle, sizeof(int));
    int *numbdr = (int*)calloc(m_num_Part, sizeof(int));
    for(int i = 0; i < m_num_TriEle; ++i){
        int n0 = mp_TriEle[i].n[0];
        int n1 = mp_TriEle[i].n[1];
        int n2 = mp_TriEle[i].n[2];
        for(int part = 0; part < m_num_Part; ++part){
            if((npart[part][n0] != -1) & (npart[part][n1] != -1) & (npart[part][n2] != -1)){
                epartTable[i] = part;
                ++num_Tri_part[part];
                if(mp_TriEle[i].bdr == 3) ++numbdr[part];
                break;
            }
        }
    }
    for(int part = 0; part < m_num_Part; ++part){
        qDebug()<<"Number of tri elements in partition "<<part<<" is "<<num_Tri_part[part];
        qDebug()<<"Number of bdr in partition "<<part<<" is "<<numbdr[part];
    }

    //PART II:DDTLM迭代
    const int MAX_ITER = 500;
    omp_set_num_threads(m_num_Part);
    QVector <umat>locs;
    QVector <mat>vals;
    QVector <int>pos(m_num_Part);
    QVector <vec>F;
    QVector <vec>Va;
    QVector <vec>Va_old;
    QVector <mat>partS;
//    umat locs1 = zeros<umat>(2, 16*m_num_TetEle+9*(numbdr[0]+numbdr[1]+numbdr[2]+numbdr[3]));
//    mat vals1 = zeros<mat>(1, 16*m_num_TetEle+9*(numbdr[0]+numbdr[1]+numbdr[2]+numbdr[3]));
//    int pos1 = 0;
    //构造稀疏矩阵
    for(int part = 0;part < m_num_Part;++part){
        umat loc(2, 16*num_Tet_part[part]+9*numbdr[part]+interfacePoints.size());
        locs.push_back(std::move(loc));
        mat val(1, 16*num_Tet_part[part]+9*numbdr[part]+interfacePoints.size());
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
    //四面体单元装配
    double St, Ft, Se, Fe;
    for(int k = 0; k < m_num_TetEle; k++){
        int tPart = m_tpartTable[k];
//         qDebug() << "k = " << k;
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                St = mp_TetEle[k].cond*(mp_TetEle[k].q[i]*mp_TetEle[k].q[j]+mp_TetEle[k].r[i]*mp_TetEle[k].r[j]+mp_TetEle[k].s[i]*mp_TetEle[k].s[j])/(36*mp_TetEle[k].Volume);
                locs[tPart](0,pos[tPart]) = npart[tPart][mp_TetEle[k].n[i]];
                locs[tPart](1,pos[tPart]) = npart[tPart][mp_TetEle[k].n[j]];
                vals[tPart](0,pos[tPart]) = St;
                partS[tPart](npart[tPart][mp_TetEle[k].n[i]], npart[tPart][mp_TetEle[k].n[j]]) = partS[tPart](npart[tPart][mp_TetEle[k].n[i]], npart[tPart][mp_TetEle[k].n[j]]) + St;
                ++pos[tPart];
//                locs1(0, pos1) = mp_TetEle[k].n[i];
//                locs1(1, pos1) = mp_TetEle[k].n[j];
//                vals1(0, pos1) = St;
//                ++pos1;

            }
            Ft = mp_TetEle[k].source*mp_TetEle[k].Volume/4;
            F[tPart](npart[tPart][mp_TetEle[k].n[i]]) = F[tPart](npart[tPart][mp_TetEle[k].n[i]]) + Ft;
        }
    }

    //三角形单元装配
    for(int k = 0; k < m_num_TriEle; k++){
        int ePart = epartTable[k];
        if(mp_TriEle[k].bdr == 3){
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                    if(i == j){
                        Se = mp_TriEle[k].h*mp_TriEle[k].Area/6;
                    }
                    else{
                        Se = mp_TriEle[k].h*mp_TriEle[k].Area/12;
                    }
                    locs[ePart](0,pos[ePart]) = npart[ePart][mp_TriEle[k].n[i]];
                    locs[ePart](1,pos[ePart]) = npart[ePart][mp_TriEle[k].n[j]];
                    vals[ePart](0,pos[ePart]) = Se;
                    partS[ePart](npart[ePart][mp_TriEle[k].n[i]], npart[ePart][mp_TriEle[k].n[j]]) = partS[ePart](npart[ePart][mp_TriEle[k].n[i]], npart[ePart][mp_TriEle[k].n[j]]) + Se;
                    ++pos[ePart];
//                    locs1(0, pos1) = mp_TriEle[k].n[i];
//                    locs1(1, pos1) = mp_TriEle[k].n[j];
//                    vals1(0, pos1) = Se;
//                    ++pos1;
                }
                Fe = mp_TriEle[k].h*mp_TriEle[k].Text*mp_TriEle[k].Area/3;
                F[ePart](npart[ePart][mp_TriEle[k].n[i]]) = F[ePart](npart[ePart][mp_TriEle[k].n[i]]) + Fe;
            }
        }
    }
    qDebug() << "assemble finish";
    //查看每个part的右侧列向量
//    std::ofstream mypartF("../tempFEM/test/partF.txt", ios::ate);
//    for(int i = 0; i < m_num_pts; ++i){
//        for(int part = 0; part < m_num_Part; ++part){
//            int n =npart[part][i];
//            if(n == -1){
//                mypartF << 0 << " ";
//            }
//            else{
//                mypartF << F[part][n] << " ";
//            }
//        }
//         mypartF << endl;
//    }
//    mypartF.close();

//    组合F，与直接法得到的F进行对比
//    std::ofstream myasmF("../tempFEM/test/asmF.txt", ios::ate);
//    vec asmF = zeros<vec>(m_num_pts);
//    for(int part = 0; part < m_num_Part; ++part){
//        for(int i = 0; i < m_num_pts; ++i){
//            if(npart[part][i] != -1){
//                asmF(i) += F[part][npart[part][i]];
//            }
//        }
//    }
//    myasmF << asmF;
////    构建节点从part到全局的检索
//    umat part2global = zeros<umat>(m_num_Part, m_num_pts);
//    for(int part = 0; part < m_num_Part; ++part){
//        int j = 0;
//        for(int i = 0; i < m_num_pts; ++i){
//            if(npart[part][i] != -1){
//                part2global(part, j) = i;
//                ++j;
//            }
//        }
//    }
//    std::ofstream mypart2global("../tempFEM/test/part2global.txt");
//    mypart2global << part2global.t();
////    组合S,通过组合每个part后直接求解来判断S装配的结果是否正确
//    umat asmlocs = zeros<umat>(2, 16*(m_num_TetEle)+9*(numbdr[0]+numbdr[1]+numbdr[2]+numbdr[3]));
//    mat asmvals = zeros<mat>(1, 16*(m_num_TetEle)+9*(numbdr[0]+numbdr[1]+numbdr[2]+numbdr[3]));
//    int asmpos = 0;
//    for(int part = 0; part < m_num_Part; ++part){
//        for(int i = 0; i < pos[part]; ++i){
//            asmlocs(0, asmpos) = part2global(part, locs[part](0,i));
//            asmlocs(1, asmpos) = part2global(part, locs[part](1,i));
//            asmvals(0, asmpos) = vals[part](0,i);
//            ++asmpos;
//        }
//    }

////    求解过程
//    sp_mat X(true, asmlocs, asmvals, m_num_pts, m_num_pts, true, true);
//    SuperMatrix sluA;
//    NCformat *Astore;
//    double   *a;
//    int      *asub, *xa;
//    int      *perm_c; /* column permutation vector */
//    int      *perm_r; /* row permutations from partial pivoting */
//    SuperMatrix L;      /* factor L */
//    SuperMatrix U;      /* factor U */
//    SuperMatrix B;
//    int      nrhs, ldx, info, m, n, nnz;
//    double   *rhs;
//    mem_usage_t   mem_usage;
//    superlu_options_t options;
//    SuperLUStat_t stat;

//    set_default_options(&options);

//    /* create matrix A in Harwell-Boeing format.*/
//    m = m_num_pts; n = m_num_pts; nnz = X.n_nonzero;
//    a = const_cast<double *>(X.values);
//    asub = (int*)const_cast<unsigned int*>(X.row_indices);
//    xa = (int*)const_cast<unsigned int*>(X.col_ptrs);
//    dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
//    Astore = (NCformat *)sluA.Store;
//    printf("Dimension %dx%d; # nonzeros %d\n", sluA.nrow, sluA.ncol, Astore->nnz);

//    nrhs = 1;
//    if (!(rhs = doubleMalloc(m * nrhs))) ABORT("Malloc fails for rhs[].");
//    //将内存拷贝过来
//    //memmove(rhs, unknown_b, 5*sizeof(double));
//    for (int i = 0; i < m; i++){
//        rhs[i] = asmF(i);
//    }
//    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

//    if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");
//    if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");

//    /* Initialize the statistics variables. */
//    StatInit(&stat);
//    dgssv(&options, &sluA, perm_c, perm_r, &L, &U, &B, &stat, &info);
//    if (info == 0) {
////                qDebug()<<"Ok.";
//        /* This is how you could access the solution matrix. */
//        double *sol = (double*)((DNformat*)B.Store)->nzval;
//        for(int i = 0; i < m_num_pts; ++i){
//           mp_3DNode[i].V = sol[i];
//           qDebug() << "V " << i << " " << mp_3DNode[i].V;
//        }
//    }else {
//        qDebug() << "info = " << info;
//    }

//    SUPERLU_FREE(rhs);
////    SUPERLU_FREE(xact);
//    SUPERLU_FREE(perm_r);
//    SUPERLU_FREE(perm_c);
////    Destroy_CompCol_Matrix(&A);
//    Destroy_SuperMatrix_Store(&B);
//    Destroy_SuperNode_Matrix(&L);
//    Destroy_CompCol_Matrix(&U);
    for(int i = 0; i < m_num_Part; i++){
        qDebug() << "pos" << i << " = " << pos[i];
//        qDebug() << "locs" << i << "size = " << 9*TriEle_num_part[i]+4*numbdr[i];
    }

    //求解过程
    int time = 0;
    QVector<int> pos2 = pos;    //记录线单元装配之后的位置
    QVector<vec> F1 = F;
    while(time++ < MAX_ITER ){
        pos = pos2;
        #pragma omp parallel for
        for(int part = 0; part < m_num_Part; ++part){
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
//            printf("Dimension %dx%d; # nonzeros %d\n", sluA.nrow, sluA.ncol, Astore->nnz);

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
            for(int j = 0; j < m_num_Part; ++j){
                if(npart[j][interfacePoints.at(i)] != -1){
                    Y += tl[j][i].Y0;
                    tl[j][i].Vi = Va[j][npart[j][interfacePoints.at(i)]] - tl[j][i].Vi;
                    I += 2*tl[j][i].Vi*tl[j][i].Y0;
                }
            }
            inter_voltage[i] = I / Y;
//            qDebug() << "inter_voltage " << i << " = " << inter_voltage[i];
        }
        //判断误差
        double outter_error = 1;
        double a = 0, b = 0;
        for(int i = 0; i < m_num_Part; ++i){
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
            for(int i = 0; i < m_num_Part; ++i){
                for(int j = 0; j < freenodepart[i]; ++j){
                    Va_old[i][j] = Va[i][j];
                }
            }
        }
    }
    //输出当前的Va
    std::ofstream mypartV("../tempFEM/test/partV.txt", ios::ate);
    for(int part = 0; part < m_num_Part; ++part){
        mypartV << Va[part] << endl << endl;
    }
    //整合结果
    std::ofstream mytemp("../tempFEM/test/temp.txt");
    double temp[15076];
    for(int part = 0; part < m_num_Part; part++){
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

    //输出测试


    std::ofstream myepartTable("../tempFEM/test/epartTable.txt");
    for(int i = 0; i < m_num_TriEle; ++i){
        myepartTable << epartTable[i] << endl;
    }
    myepartTable.close();

    //求解结束，回收空间
    free(num_Tet_part);
    for(int i = 0; i < m_num_Part; ++i){
        free(npart[i]);
    }
    free(npart);
    for(int i = 0; i < m_num_Part; ++i){
        free(tl[i]);
    }
    free(tl);
    free(inter_voltage);
    free(freenodepart);
    free(num_Tri_part);
    free(epartTable);
    //    free(TetRM);
    //    free(TriRM);
    return 0;
}
