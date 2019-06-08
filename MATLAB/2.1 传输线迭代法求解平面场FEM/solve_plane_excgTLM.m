%简单模型二维平面热传导温度场的传输线法求解，引入热交换边界条件（第二类边界条件）
%材料的热导率为52W/(m.K)，边界的温度为273.15K，边界上热通量为5e5
%参考颜威利《电气工程电磁场数值分析》 P180 温度场有限元分析
%该书关于温度场的第二类边界条件和第三类边界条件的叙述有误，但是方法是没错的
%引入线单元
%2019/2/26
%By QzLancer

%-------------------------------读取分网文件
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_heatexcg.mphtxt');
%-------------------------------初始化参数
Cond = 52;
BoundTemp = 273.15;
HeatFlux = 5e5;
e = 0.001;
Y0 = 1;
%-------------------------------根据插值要求求解各单元的几何参数，目前只考虑内部单元
%导出所有三角形单元RZ轴的坐标
R = Coor(:,1);
Z = Coor(:,2);
TriR = R(TriElement);
TriZ = Z(TriElement);
%计算出所有三角形单元的p,q,r和Area
p(:,1) = TriR(:,2).*TriZ(:,3) - TriZ(:,2).*TriR(:,3);
p(:,2) = TriR(:,3).*TriZ(:,1) - TriZ(:,3).*TriR(:,1);
p(:,3) = TriR(:,1).*TriZ(:,2) - TriZ(:,1).*TriR(:,2);
q(:,1) = TriZ(:,2) - TriZ(:,3);
q(:,2) = TriZ(:,3) - TriZ(:,1);
q(:,3) = TriZ(:,1) - TriZ(:,2);
r(:,1) = TriR(:,3) - TriR(:,2);
r(:,2) = TriR(:,1) - TriR(:,3);
r(:,3) = TriR(:,2) - TriR(:,1);
Area = (q(:,1).*r(:,2) - q(:,2).*r(:,1))/2;
%计算出所有线单元的坐标
EdgR = R(EdgElement);
EdgZ = Z(EdgElement);
%过滤掉不需要的线单元，只保留有热交换的线单元，保存其节点和坐标
pHeatExcgBond = find(EdgR(:, 1) == 0.02 & EdgR(:, 2) == 0.02 & ...
                      EdgZ(:, 1) >= 0.04 & EdgZ(:, 2) >= 0.04 & ...
                      EdgZ(:, 1) <= 0.1 & EdgZ(:, 2) <= 0.1);
EdgR1 = EdgR(:, 1);
HeatExcgElement = EdgElement(pHeatExcgBond, :);
HeatExcgR = EdgR(pHeatExcgBond, :);
HeatExcgZ = EdgZ(pHeatExcgBond, :);
%计算所有线单元的d
d = sqrt((HeatExcgR(:, 1) - HeatExcgR(:, 2)).^2 + (HeatExcgZ(:, 1) - HeatExcgZ(:, 2)).^2);
%------------------------------系数矩阵装配，面单元F装配（该模型不存在面单元F）
Se = zeros(3,3,length(TriElement));
F = zeros(length(Coor), 1);
Ym = zeros(length(Coor),length(Coor));
for k = 1:length(TriElement)
    for i = 1:3
        for j = 1:3
            if i == j
                Ym(TriElement(k,i),TriElement(k,j)) = Ym(TriElement(k,i),TriElement(k,j))+2*Y0;
            else
                Ym(TriElement(k,i),TriElement(k,j)) = Ym(TriElement(k,i),TriElement(k,j))-Y0;
            end
            Se(i,j,k) = Cond.*(r(k,i)*r(k,j)+q(k,i)*q(k,j))/(4*Area(k));
        end
    end
end
%------------------------------线单元F装配
for k = 1:length(HeatExcgElement)
    for i = 1:2
        Fl = HeatFlux*d(k)/2;
        F(HeatExcgElement(k, i)) = F(HeatExcgElement(k, i)) + Fl;
    end
end
%------------------------------求解
%检索出边界点
% EquTemp = zeros(length(Coor),1);
Temp = zeros(length(Coor),1);
Boundary = find(Z==0 | Z==0.14 | R==0.1);
FreeNodes = find(~(Z==0 | Z==0.14 | R==0.1));
%入射和反射过程
Va0 = ones(length(Coor),1);
Va1 = zeros(length(Coor),1);
Vr = zeros(3,3,length(TriElement));
Vc = zeros(3,3,length(TriElement));                              
Vi = zeros(3,3,length(TriElement));
I = zeros(length(Coor),1);
m = 1;
while norm(Va1-Va0)/norm(Va0) >= e
    Va0 = Va1;
    Va1 = zeros(length(Coor),1);
    Va1(Boundary) = BoundTemp;
    F1 = Ym*Va1;
    F2 = F+I-F1;
    Va1(FreeNodes) = Ym(FreeNodes,FreeNodes)\F2(FreeNodes);
    I = zeros(length(Coor),1);
    for k = 1:length(TriElement)
        for i = 1:3
            for j = 1:3
                Vr(i,j,k) = Va1(TriElement(k,i))-Va1(TriElement(k,j))-Vi(i,j,k);
                Vc(i,j,k) = 2*Vr(i,j,k)*Y0/(Y0-Se(i,j,k));
                Vi(i,j,k) = Vc(i,j,k) - Vr(i,j,k);
                I(TriElement(k,i)) = I(TriElement(k,i))+2*Y0*Vi(i,j,k);                
            end
        end
    end
    m = m+1;
end
%------------------------------后处理
Interp1 = scatteredInterpolant(R,Z,Va1);
tx = 0.02:1e-3:0.1;
ty = 0:1e-3:0.14;
[qx,qy] = meshgrid(tx,ty);
qz = Interp1(qx,qy);
% mesh(qx,qy,qz);
subplot(1,2,2);
contourf(qx,qy,qz,20);colorbar;
hold on;
EdgNode1 = EdgElement(:, 1);
EdgNode2 = EdgElement(:, 2);
EdgCoor1 = Coor(EdgNode1, :);
EdgCoor2 = Coor(EdgNode2, :);
% hold on;
% plot(EdgCoor1(:, 1), EdgCoor1(:, 2), 'x');
% hold on;
% plot(EdgCoor2(:, 1), EdgCoor2(:, 2), 'x');
Edglength = length(EdgElement); 
for i = 1:Edglength
    plot([EdgCoor1(i, 1), EdgCoor2(i, 1)], [EdgCoor1(i, 2), EdgCoor2(i, 2)],'Color','k','Linewidth',1);
    hold on;
end
%-------------------------读取COMSOL计算出来的结果并进行后处理
[fileID, Errmessage] = fopen('solution_heatexcg.txt', 'r');
if fileID == -1
    disp(Errmessage);
end
for i=1:9
    fgetl(fileID);
end
comsoldata = fscanf(fileID,'%lf %lf %lf\n',[3,length(Coor)]);
comsoldata = comsoldata';
fclose(fileID);
Interp2 = scatteredInterpolant(comsoldata(:,1),comsoldata(:,2),comsoldata(:,3));
qz = Interp2(qx,qy);
subplot(1,2,1);
contourf(qx,qy,qz,20);colorbar;