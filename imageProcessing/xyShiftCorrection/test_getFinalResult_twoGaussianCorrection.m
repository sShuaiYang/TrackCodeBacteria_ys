function [result,coef]=test_getFinalResult_twoGaussianCorrection(dirFile,index)
% 算出来的稀释率有可能出现双峰，可以预见，
% 较低的峰对应的就是分解率可以忽略条件下的生长率
% 较高的峰对应的就是分解率较高时候的情况，程序需要进行判定、
clc
disp(dirFile)
nameList=dir(dirFile);
gamaAll=[];
fpProduceAll=[];
miuAll=[];
GRratioAll=[];
isDeadAll=[];
gfpAll=[];
rfpAll=[];
for i=1:numel(nameList)-2
    load(strcat(dirFile,'\',nameList(i+2).name,'\miu&Presult\3_只计算第一个点_new_1点.mat'))
%     load(strcat(dirFile,'\',nameList(i+2).name,'\miu&Presult\1.mat'))
    gamaAll=[gamaAll;gama];
    fpProduceAll=[fpProduceAll;fpProduce];
    miuAll=[miuAll;miu];
    GRratioAll=[GRratioAll;GRratio];
    isDeadAll=[isDeadAll;isDead];
    gfpAll=[gfpAll;gfp];
    rfpAll=[rfpAll;rfp];
end
result.gamaAll=gamaAll;   
result.fpProduceAll=fpProduceAll;   
result.miuAll=miuAll;
result.GRratioAll=GRratioAll;
result.isDead=isDeadAll;
result.gfp=gfpAll;
result.rfp=rfpAll;
coef=1;
% b=result.GRratioAll;
a=result.gamaAll;
% a=result.fpProduceAll;
b1=result.miuAll;

hold on;subplot(2,1,2);h=plot(a,b1,'.','linestyle','none','markerSize',2);
xlim([0,0.03]);
ylim([0,0.03]);
hold on;plot(0:0.001:0.03,0:0.001:0.03,'r')

% % b1=b1*0.9176+0.00163;
% % hold on;plot(b,a,'lineStyle','none','Marker','.')
% a(isnan(a))=[];
% b1(isnan(b1))=[];
% b1(abs(b1)==inf)=[];
% 
% % [b2,a2]=hist(b,linspace(0,10,100));
% [b2,a2]=hist(b1,linspace(0,0.03,60));
% [b1,a1]=hist(a,linspace(0,0.03,60));
% % [b2,a2]=hist(b1,linspace(-50,50,300));
% % [b1,a1]=hist(a,linspace(-50,50,300));
% a1(end)=[];
% b1(end)=[];
% a1(1)=[];
% b1(1)=[];
% a2(end)=[];
% b2(end)=[];
% a2(1)=[];
% b2(1)=[];
% 
% % a1,b1为真实生长率   a2,b2为计算出的稀释率
% % [yb2,coef]=detailAnalysis(a2,b2);
% yb2=b2;
% coef=1;
% if coef==1
%     plot(a1,b1./sum(b1));hold on;plot(a2,b2./sum(b2),'g');plot(a2,yb2./sum(yb2),'r');
% %     result=[gaussFit(a1,b1),gaussFit(a2,b2)];
% else
%     plot(a1,b1./sum(b1));hold on;plot(a2,b2./sum(b2),'g');plot(a2,yb2./sum(yb2),'r');
%     result=[gaussFit(a1,b1),gaussFit(a2,yb2)];
% end
% close all;
end
function b=gaussFit(xData,yData)
ft = fittype( 'gauss1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
u=xData(yData==max(yData));
opts.StartPoint = [max(yData),u(1),u(1)/10];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData', yData', ft, opts );
gof1.adjrsquare=0;
n=0;
while gof1.adjrsquare~=gof.adjrsquare && gof.adjrsquare>=0.99
    n=n+1;
    if n==10
        break
    end
    gof1=gof;
    opts.StartPoint=[fitresult.a1,fitresult.b1,fitresult.c1];
    [fitresult, gof] = fit( xData', yData', ft, opts );
end
b=fitresult.b1;
end
function [yb2,coef]=detailAnalysis(a2,b2)
ft = fittype( 'gauss1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
u=a2(b2==max(b2));
opts.StartPoint = [max(b2),u(1),u(1)/10];
opts.Upper = [Inf Inf Inf];
opts.Robust = 'Bisquare';
% Fit model to data.
[fitresult, gof] = fit( a2',b2', ft, opts );
if gof.adjrsquare<=0.995
    coef=2;
else
    coef=1;
    yb2=fitresult.a1*exp(-((a2-fitresult.b1)/fitresult.c1).^2);
    return
end
clear opts
ft = fittype( 'gauss2' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
% Fit model to data.
opts.Upper = [2*max(b2) 0.05 Inf 2*max(b2) 0.05 Inf];
opts.Lower = [0 0 0 0 0 0];
[fitresult, gof] = fit( a2', b2', ft, opts );
gof1.adjrsquare=0;
n=0;
while gof1.adjrsquare~=gof.adjrsquare && gof.adjrsquare>=0.9
    n=n+1;
    if n==10
        break
    end
    gof1=gof;
    opts.StartPoint=[fitresult.a1,fitresult.b1,fitresult.c1,fitresult.a2,fitresult.b2,fitresult.c2];
    [fitresult, gof] = fit( a2', b2', ft, opts );
end
if fitresult.b1<fitresult.b2
    yb2=fitresult.a1*exp(-((a2-fitresult.b1)/fitresult.c1).^2);
else
    yb2=fitresult.a2*exp(-((a2-fitresult.b2)/fitresult.c2).^2);
end
end