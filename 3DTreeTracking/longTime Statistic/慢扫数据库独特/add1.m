function add1(p0)
% p0���ļ�·���������������������ɨmat�ļ������ó���Ƭ������
dt=dir(p0);
p=[dt.isdir];pt=dt(p);
n=length(pt);
for k=3:n-1
    m=strfind(pt(k).name,'PA');
    name{k}=pt(k).name(m:m+5);
end
for i=1:61
    dirFile=strcat(p0,'\allResult\allResult',num2str(i),'.mat');
    load(dirFile);
    allResult.name=name{1,i+2};
    newname=strcat(p0,'\newallResult\allResult',num2str(i));
    save(newname,'allResult');
    clear allResult
end
end