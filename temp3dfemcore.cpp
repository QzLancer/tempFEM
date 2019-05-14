#include "temp3dfemcore.h"



CTemp3DFEMCore::CTemp3DFEMCore(Widget *parent, const char *fn):
    m_COMSOLMesh(fn),
    mp_3DNode(nullptr),
    mp_VtxEle(nullptr),
    mp_EdgEle(nullptr),
    mp_TriEle(nullptr),
    mp_TetEle(nullptr)
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

int CTemp3DFEMCore::preCalculation()
{


    return 0;
}
