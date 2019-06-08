%----------------------------------COMSOL结果读取
[fileID, Errmessage] = fopen('COMSOL3DNR_87779.txt', 'r');
if fileID == -1
    disp(Errmessage);
end
for i=1:9
    fgetl(fileID);
end
figure('color',[1 1 1]);
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
[fileID, Errmessage] = fopen('../../test/Temp3DDDTLM_87779.txt', 'r');
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
title('DDTLM','FontName','Times New Roman','FontSize',15);
axis equal;
%--------------------------绘制节点的误差比较图
figure('color',[1 1 1]);
ppz = pz(:,26);
p = [1,2,3,4,5,10,15,20,25,30,35,40,45,46,47,48,49,50,51];
ppy = py(p,26);
ppTemp1 = pTemp1(p,26);
ppTemp2 = pTemp2(p,26);
scatter(ppy,ppTemp1,'filled','^','MarkerFaceColor',[.5 0 0]);
hold on;
scatter(ppy,ppTemp2,'filled','MarkerFaceColor',[0 0 .5]);
hold on;
grid on;
err1 = 0;
for i = 1:length(p)
    err1 = err1 + abs(ppTemp1(i)-ppTemp2(i))/(ppTemp1(i)*length(p));
end