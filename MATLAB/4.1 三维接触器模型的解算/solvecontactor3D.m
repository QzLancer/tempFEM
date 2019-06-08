%绘制三维接触器COMSOL求解结果和MATLAB求解结果的云图，并且绘图分析
%By QzLancer
%2019/5/27
clear all;
%-------------------------------读取分网文件
clear all;
[Coor,VtxElement,VtxEntity,TriElement,TriEntity,TetElement,TetEntity] = read3Dmesh('../../model/mesh_contactor3D_87779.mphtxt');
tic;
%-------------------------------初始化参数
BoundTemp = 273.15;
HeatFlux = 5e3;
h = 20;
T0 = 293.15;
%-------------------------------根据插值要求求解四面体单元的几何参数
%导出四面体单元XYZ的坐标
X = Coor(:, 1);
Y = Coor(:, 2);
Z = Coor(:, 3);
TetX = X(TetElement);
TetY = Y(TetElement);
TetZ = Z(TetElement);
TetLen = length(TetElement);
%求解各单元体积Volume
Volume = zeros(TetLen, 1);
for i = 1:TetLen
    Ve = det([ones(4, 1), TetX(i, :)', TetY(i, :)', TetZ(i, :)']) / 6;
    Volume(i, 1) = Ve;
end
clear Ve;
%求解p,q,r,s
p = zeros(TetLen, 4);
q = zeros(TetLen, 4);
r = zeros(TetLen, 4);
s = zeros(TetLen, 4);
%求解p
for i = 1:TetLen
    p(i, 1) = det([TetX(i, [2,3,4])', TetY(i, [2,3,4])', TetZ(i, [2,3,4])']);
    p(i, 2) = -det([TetX(i, [1,3,4])', TetY(i, [1,3,4])', TetZ(i, [1,3,4])']);
    p(i, 3) = det([TetX(i, [1,2,4])', TetY(i, [1,2,4])', TetZ(i, [1,2,4])']);
    p(i, 4) = -det([TetX(i, [1,2,3])', TetY(i, [1,2,3])', TetZ(i, [1,2,3])']);
end
%求解q
for i = 1:TetLen
    q(i, 1) = -det([ones(3, 1), TetY(i, [2,3,4])', TetZ(i, [2,3,4])']);
    q(i, 2) = det([ones(3, 1), TetY(i, [1,3,4])', TetZ(i, [1,3,4])']);
    q(i, 3) = -det([ones(3, 1), TetY(i, [1,2,4])', TetZ(i, [1,2,4])']);
    q(i, 4) = det([ones(3, 1), TetY(i, [1,2,3])', TetZ(i, [1,2,3])']);
end
% 求解r
for i = 1:TetLen
    r(i,1) = -det([TetX(i, [2,3,4])', ones(3, 1), TetZ(i, [2,3,4])']);
    r(i,2) = det([TetX(i, [1,3,4])', ones(3, 1), TetZ(i, [1,3,4])']);
    r(i,3) = -det([TetX(i, [1,2,4])', ones(3, 1), TetZ(i, [1,2,4])']);
    r(i,4) = det([TetX(i, [1,2,3])', ones(3, 1), TetZ(i, [1,2,3])']);
end
% 求解s
for i = 1:TetLen
    s(i,1) = -det([TetX(i, [2,3,4])', TetY(i, [2,3,4])', ones(3, 1)]);
    s(i,2) = det([TetX(i, [1,3,4])', TetY(i, [1,3,4])', ones(3, 1)]);
    s(i,3) = -det([TetX(i, [1,2,4])', TetY(i, [1,2,4])', ones(3, 1)]);
    s(i,4) = det([TetX(i, [1,2,3])', TetY(i, [1,2,3])', ones(3, 1)]);
end
%-------------------------------求解边界三角形单元的几何参数
%找出热交换区域节点编号和坐标
pHeatExcgBond = find((TriEntity == 1) | (TriEntity == 2) | (TriEntity == 3) | (TriEntity == 4) |...
(TriEntity == 5) | (TriEntity == 6) | (TriEntity == 91) | (TriEntity == 92) |...
                     (TriEntity == 93) | (TriEntity == 137) | (TriEntity == 143) | (TriEntity == 183));
TriExcgElement = TriElement(pHeatExcgBond, :);
TriExcgX = X(TriExcgElement);
TriExcgY = Y(TriExcgElement);
TriExcgZ = Z(TriExcgElement);
%设置热导率
cond = zeros(length(TetElement),1);
source = zeros(length(TetElement),1);
pEntity = find((TetEntity == 7) | (TetEntity == 9));
cond(pEntity) = 400;
source(pEntity) = 500000;
pEntity = find((TetEntity == 2) | (TetEntity == 4) | (TetEntity == 10) | (TetEntity == 12) | (TetEntity == 13));
cond(pEntity) = 76.2;
pEntity = find((TetEntity == 1) | (TetEntity == 6) | (TetEntity == 8));
cond(pEntity) = 0.26;
pEntity = find((TetEntity == 3) | (TetEntity == 5) | (TetEntity == 11));
cond(pEntity) = 0.03;
%求解两条向量，叉乘的模和三维三角形的面积
vec1(:, 1) = TriExcgX(:, 2) - TriExcgX(:, 1);
vec1(:, 2) = TriExcgY(:, 2) - TriExcgY(:, 1);
vec1(:, 3) = TriExcgZ(:, 2) - TriExcgZ(:, 1);
vec2(:, 1) = TriExcgX(:, 3) - TriExcgX(:, 1);
vec2(:, 2) = TriExcgY(:, 3) - TriExcgY(:, 1);
vec2(:, 3) = TriExcgZ(:, 3) - TriExcgZ(:, 1);
vec = cross(vec1,vec2);
Area = sqrt(vec(:, 1).^2 + vec(:, 2).^2 + vec(:, 3).^2) / 2;
%-------------------------------四面体单元分析和装配
S = zeros(length(Coor));
F = zeros(length(Coor), 1);
for k = 1:TetLen
    for i = 1:4
        for j = 1:4
            Se = cond(k)*(q(k,i)*q(k,j)+r(k,i)*r(k,j)+s(k,i)*s(k,j))/(36*Volume(k));
            S(TetElement(k,i), TetElement(k,j)) = S(TetElement(k,i), TetElement(k,j)) + Se;
        end
        Ft = source(k)*Volume(k)/4;
        F(TetElement(k,i)) = F(TetElement(k,i)) + Ft;
    end
end
%-------------------------------边界三角形单元分析和装配
for k = 1:length(TriExcgElement)
    for i = 1:3
        for j = 1:3
            if(i==j)
                Se = h*Area(k)/6;
            else
                Se = h*Area(k)/12;
            end
            S(TriExcgElement(k,i),TriExcgElement(k,j)) =  S(TriExcgElement(k,i),TriExcgElement(k,j)) + Se;
        end
        Fe = h*T0*Area(k)/3;
        F(TriExcgElement(k,i)) = F(TriExcgElement(k,i)) + Fe;
    end
end
%-------------------------------采用直接法求解
Temp3 = S\F;
toc;


%----------------------------------COMSOL结果读取
[fileID, Errmessage] = fopen('solve3D.txt', 'r');
if fileID == -1
    disp(Errmessage);
end
for i=1:9
    fgetl(fileID);
end
comsoldata = fscanf(fileID,'%lf %lf %lf %lf\n',[4,999999]);
comsoldata = comsoldata';
tx = -0.025:0.001:0.025;
ty = -0.025:0.001:0.025;
tz = -0.025:0.001:0.025;
Interp1 = scatteredInterpolant(comsoldata(:,1),comsoldata(:,2),comsoldata(:,3),comsoldata(:,4));
[qx, qy, qz] = meshgrid(tx, ty, tz);
Temp1 = Interp1(qx, qy, qz);
py = squeeze(qy(:,28,:));
pz = squeeze(qz(:,28,:));
pTemp1 = squeeze(Temp1(:,28,:));
subplot(1,3,1);
contourf(py,pz,pTemp1,20);colorbar;
xlabel('y','FontName','Times New Roman','FontSize',15);
ylabel('z','FontName','Times New Roman','FontSize',15);
title('COMSOL','FontName','Times New Roman','FontSize',15);
axis equal;
%----------------------------------C++结果读取
[fileID, Errmessage] = fopen('../../test/Temp3D_contactor_87779.txt', 'r');
if fileID == -1
    disp(Errmessage);
end
cppdata = fscanf(fileID,'%lf %lf %lf %lf\n',[4,999999]);
cppdata = cppdata';
% Interp2 = scatteredInterpolant(cppdata(:,1),cppdata(:,2),cppdata(:,3),cppdata(:,4));
Interp2 = scatteredInterpolant(cppdata(:,1),cppdata(:,2),cppdata(:,3),cppdata(:,4));
[qx, qy, qz] = meshgrid(tx, ty, tz);
Temp2 = Interp2(qx, qy, qz);
py = squeeze(qy(:,28,:));
pz = squeeze(qz(:,28,:));
pTemp2 = squeeze(Temp2(:,28,:));
subplot(1,3,3);
contourf(py,pz,pTemp2,20);colorbar;
xlabel('y','FontName','Times New Roman','FontSize',15);
ylabel('z','FontName','Times New Roman','FontSize',15);
title('C++','FontName','Times New Roman','FontSize',15);
axis equal;
%------------------------------MATLAB结果
Interp3 = scatteredInterpolant(cppdata(:,1),cppdata(:,2),cppdata(:,3),Temp3);
[qx, qy, qz] = meshgrid(tx, ty, tz);
Temp3 = Interp3(qx, qy, qz);
py = squeeze(qy(:,28,:));
pz = squeeze(qz(:,28,:));
pTemp3 = squeeze(Temp3(:,28,:));
subplot(1,3,2);
contourf(py,pz,pTemp3,20);colorbar;
xlabel('y','FontName','Times New Roman','FontSize',15);
ylabel('z','FontName','Times New Roman','FontSize',15);
title('MATLAB','FontName','Times New Roman','FontSize',15);
axis equal;
%--------------------------绘制节点的误差比较图
figure(2);
ppz = pz(:,26);
p = [1,2,3,4,5,10,15,20,25,30,35,40,45,46,47,48,49,50,51];
ppy = py(p,26);
ppTemp1 = pTemp1(p,26);
ppTemp2 = pTemp3(p,26);
ppTemp3 = pTemp2(p,26);
scatter(ppy,ppTemp1,'filled','MarkerFaceColor',[.5 0 0]);
hold on;
scatter(ppy,ppTemp2,'filled','^','MarkerFaceColor',[0 .5 0]);
hold on;
scatter(ppy,ppTemp3,'filled','d','MarkerFaceColor',[0 0 .5]);
grid on;
err1 = 0;
err2 = 0
for i = 1:length(p)
    err1 = err1 + abs(ppTemp1(i)-ppTemp2(i))/(ppTemp1(i)*length(p));
    err2 = err2 + abs(ppTemp2(i)-ppTemp3(i))/(ppTemp2(i)*length(p));
end