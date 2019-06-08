%-----------------------------------读取C++结果
[fileID, Errmessage] = fopen('../../test/TempDDTLM_13666.txt', 'r');
tx = 0:1e-4:0.026;
ty = -0.025:1e-4:0.024;
[qx,qy] = meshgrid(tx,ty);
if fileID == -1
    disp(Errmessage);
end
cppdata = fscanf(fileID, '%lf %lf %lf\n', [3, 99999]);
cppdata = cppdata';
fclose(fileID);
interp3 = scatteredInterpolant(cppdata(:,1), cppdata(:,2), cppdata(:,3));
figure('color',[1 1 1]);
qz3 = interp3(qx, qy);
subplot(1,3,3);
contourf(qx,qy,qz3,20);colorbar;
hold on;
title('DDTLM','FontName','Times New Roman','FontSize',15);
axis equal;
%-------------------------读取COMSOL计算出来的结果并进行后处理
[fileID, Errmessage] = fopen('COMSOLNR_13666.txt', 'r');
if fileID == -1
    disp(Errmessage);
end
for i=1:9
    fgetl(fileID);
end
comsoldata = fscanf(fileID,'%lf %lf %lf\n',[3,99999]);
comsoldata = comsoldata';
fclose(fileID);
Interp2 = scatteredInterpolant(comsoldata(:,1),comsoldata(:,2),comsoldata(:,3));
qz1 = Interp2(qx,qy);
subplot(1,3,1);
contourf(qx,qy,qz1,20);colorbar;
title('COMSOL','FontName','Times New Roman','FontSize',15);
axis equal;
%-----------------------绘制节点的误差比较图
px = qx(251,:);
p = [51,101,151,201,221,231,236,241,246,251];
px = px(p);
pz1 = qz1(251,p);
pz3 = qz3(251,p);
figure('color',[1 1 1]);
scatter(px,pz1,'filled','^','MarkerFaceColor',[.5 0 0]);
hold on;
scatter(px,pz3,'filled','MarkerFaceColor',[0 0 .5]);
grid on;
err1 = 0;
for i = 1:length(p)
    err1 = err1 + abs(pz1(i)-pz3(i))/(pz1(i)*length(p));
end