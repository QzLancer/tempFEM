%����άCOMSOL�����ļ�ͨ��mpmetis.exe�ֽ�
%��Ҫ�����������У�
%By QzLancer
%2019/4/11
%------------------------------------����readcomsol����ȡCOMSOL�����ļ�
clear all;
close all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_source.mphtxt');
fp=fopen('mesh_source.mpmetis','w');
a=length(TriElement);
fprintf(fp,'%d\n',a);
fprintf(fp,'%d %d %d\n',TriElement');
system('mpmetis.exe mesh_source.mpmetis 2');