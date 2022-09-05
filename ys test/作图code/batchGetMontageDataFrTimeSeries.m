function [dataAll,dataOri] = batchGetMontageDataFrTimeSeries()
%cxy ������������� �����յ�rsmY reporter���
dirAll='F:\GacS\cxy20181205';
dataAll={};
dataOri={};
AllList=dir(dirAll);
for ifile=1:numel(AllList)-2 %cxy�ļ��п�ͷ��#�����Դӵ�һ����ʼ
    dirFile=[dirAll,'\',AllList(ifile+2).name];
    disp(AllList(ifile+2).name);
    [~,~,~,dataHist_sf,dataHist_mS] = getMultiFieldImageForcxy(dirFile);
    
    temp=dataHist_mS.intensityCollect;
    temp1=dataHist_sf.intensityCollect;
    
    dataOri{ifile}(:,1)=temp1';
    dataOri{ifile}(:,2)=temp';
    
%     b=temp1<=2000 & temp>3 & temp<=500;%ֱ�����߼����� ����  
    
    b=temp1<=3000 & temp>3 & temp<=1000;%ֱ�����߼����� ����  
    a=[];
    a(1,:)=dataHist_sf.intensityCollect(b);
    a(2,:)=dataHist_mS.intensityCollect(b);
    a=a';
    dataAll{ifile}=a;
    clear dataHist_sf dataHist_mS 
end
save([dirAll,'\dataAll.mat'],'dataAll');
save([dirAll,'\dataOri.mat'],'dataOri');
end

