function [a,b]=test_getFinalResult(dirFile,iMap)
colorAll=colormap(jet(8));
% dirFile='F:\data\2015-11-18 EMCCD HJD1_B test\2015-11-18 jzy test1\allResult\';
nameList=dir(dirFile);
gamaAll=[];
fpProduceAll=[];
miuAll=[];
GRratioAll=[];
gfpAll=[];
rfpAll=[];
deltaAll=[];
for i=1:numel(nameList)-2
% for i=1:400
    load(strcat(dirFile,'\',nameList(i+2).name,'\miu&Presult\3_只计算第一个点_new_1点.mat'))
%     load(strcat(dirFile,'\',nameList(i+2).name,'\miu&Presult\1.mat'))
    gamaAll=[gamaAll;gama];
    fpProduceAll=[fpProduceAll;fpProduce];
    miuAll=[miuAll;miu];
    GRratioAll=[GRratioAll;GRratio];
%     gfpAll=[gfpAll;gfp];
%     rfpAll=[rfpAll;rfp];
%     deltaAll=[deltaAll;deltaRFP];
end
result.gfpAll=gfpAll;
result.rfpAll=rfpAll;
result.gamaAll=gamaAll;
result.deltaAll=deltaAll;
% return
% hold on;plot(rfpAll./gfpAll,gamaAll,'lineStyle','none','marker','.')
% result.gamaAll=gamaAll;   
result.fpProduceAll=fpProduceAll;   
result.miuAll=miuAll;
result.GRratioAll=GRratioAll;

a=result.fpProduceAll;
a(a>=300)=[];
a(a<=-100)=[];
[b,a]=hist(a,linspace(-100,300,400));
hold on;plot(a,b./sum(b),'color',colorAll(iMap,:));
a=mean(result.fpProduceAll);

% a=result.gamaAll;
% a(a>=0.05)=[];
% a(a<=-0.02)=[];
% [b,a]=hist(a,linspace(-0.02,0.05,400));
% hold on;plot(a,b./sum(b),'color',colorAll(iMap,:));
% % a=mean(result.gamaAll);
% a=gaussFit(a,b);
return                           
b1=result.miuAll;
% b1=result.miuAll;
u=mean(b1);
disp(u);
% b1=b1*0.9176+0.00163;
% hold on;plot(b,a,'lineStyle','none','Marker','.')
a(isnan(a))=[];
b1(isnan(b1))=[];
% [b2,a2]=hist(b,linspace(0,10,100));
[b2,a2]=hist(b1,linspace(0.001,0.05,100));
[b1,a1]=hist(a,linspace(0.001,0.05,100));
a1(end)=[];
b1(end)=[];
a1(1)=[];
b1(1)=[];
a2(end)=[];
b2(end)=[];
a2(1)=[];
b2(1)=[];
hold on;
% subplot(2,3,iMap);
% plot(result.GRratioAll,result.rfpAll.*(log(2)/62*(result.GRratioAll-1)-result.gamaAll),'linestyle','none','marker','.','markerSize',2);hold on;plot(-0.02:0.01:0.02,-0.02:0.01:0.02,'r');
plot(result.fpProduceAll,result.miuAll,'linestyle','none','marker','.','markerSize',2);hold on;plot(-0.02:0.01:0.02,-0.02:0.01:0.02,'r');
xlim([0,0.02])
ylim([0,0.02])

% plot(a1,b1./sum(b1));hold on;
% plot(a2,b2./sum(b2),'g')
% hold on
% close all
% result=[gaussFit(a1,b1),gaussFit(a2,b2)];
% gama=result.gamaAll;
% fpProduce=result.fpProduceAll;
% miu=result.miuAll;
end
function a
b=result.GRratioAll;
a=result.gamaAll;
a1=linspace(0,0.02,51);
clear b1
clear c1
for i=1:100
b1(i)=mean(b(a>a1(i)&a<=a1(i+1)));
c1(i)=mean([a1(i),a1(i+1)]);
end
c1(isnan(b1))=[];
b1(isnan(b1))=[];
end
function b=gaussFit(xData,yData)
ft = fittype( 'gauss1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [-Inf -Inf 0];
u=xData(yData==max(yData));
opts.StartPoint = [max(yData),u(1),0.001];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData', yData', ft, opts );
b=fitresult.b1;
end