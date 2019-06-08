%将二维COMSOL分网文件通过mpmetis.exe分解
%需要导出的数据有：
%By QzLancer
%2019/4/11
%------------------------------------调用readcomsol，读取COMSOL分网文件
clear all;
close all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_source.mphtxt');
fp=fopen('mesh_source.mpmetis','w');
a=length(TriElement);
fprintf(fp,'%d\n',a);
fprintf(fp,'%d %d %d\n',TriElement');
system('mpmetis.exe mesh_source.mpmetis 2');