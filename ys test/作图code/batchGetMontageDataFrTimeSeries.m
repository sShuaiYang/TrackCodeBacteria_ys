function [dataAll,dataOri] = batchGetMontageDataFrTimeSeries()
%cxy 数据批处理程序 蓝光诱导rsmY reporter表达
dirAll='F:\GacS\cxy20181205';
dataAll={};
dataOri={};
AllList=dir(dirAll);
for ifile=1:numel(AllList)-2 %cxy文件有开头是#，所以从第一个开始
    dirFile=[dirAll,'\',AllList(ifile+2).name];
    disp(AllList(ifile+2).name);
    [~,~,~,dataHist_sf,dataHist_mS] = getMultiFieldImageForcxy(dirFile);
    
    temp=dataHist_mS.intensityCollect;
    temp1=dataHist_sf.intensityCollect;
    
    dataOri{ifile}(:,1)=temp1';
    dataOri{ifile}(:,2)=temp';
    
%     b=temp1<=2000 & temp>3 & temp<=500;%直接用逻辑矩阵 方便  
    
    b=temp1<=3000 & temp>3 & temp<=1000;%直接用逻辑矩阵 方便  
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

