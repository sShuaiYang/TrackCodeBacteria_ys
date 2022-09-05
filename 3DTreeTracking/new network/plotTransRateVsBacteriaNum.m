% function [bacteriaNum,expNum]=plotTransRateVsBacteriaNum(resultAll)
% for iResult=14:size(resultAll,2)
%     distMatrix=full(resultAll{iResult}.finalMatrix);
%     bacteriaNum(iResult)=size(distMatrix,1);
%     degreeNum=[];
%     for i=1:size(distMatrix,1)
%         iLine=distMatrix(i,:);
%         iLine=iLine/160000*20/10/(10.^(-4))*10^(-5);
%         iLine=1-iLine;
%         degreeNum(i)=1-prod(iLine);
%     end
%     regionHist=[];
%     maxDegree=min(0.05,(max(degreeNum)+0.0001));
%     regionHist(1,:)=0:0.000001:maxDegree;
%     for i=1:numel(regionHist)
%         regionHist(2,i)=numel(degreeNum(degreeNum>=regionHist(1,i)));
%     end
%     hold on
%     if mod(iResult,3)==0
%         plot(regionHist(1,:),regionHist(2,:),'r')
%     end
%     [expNum(iResult),a]=prepareCurveData(regionHist(1,:),regionHist(2,:));
%     if mod(iResult,3)==0
%         plot(regionHist(1,:),a*exp(-regionHist(1,:)*expNum(iResult)),'b')
%     end
%     warning off %#ok<WNOFF>
% end
% end
function [bacteriaNum,expNum]=plotTransRateVsBacteriaNum(result)
% 对直连的网络的统计
result=result{1};
for iResult=20:size(result,2)
    distMatrix=full(result(iResult).finalMatrix);
    bacteriaNum(iResult)=size(distMatrix,1);
    degreeNum=[];
    for i=1:size(distMatrix,1)
        iLine=distMatrix(i,:);
        degreeNum(i)=numel(iLine(iLine==1));
    end
    regionHist=[];
    maxDegree=max(degreeNum);
    regionHist(1,:)=1:1:maxDegree;
    for i=1:numel(regionHist)
        regionHist(2,i)=numel(degreeNum(degreeNum>=regionHist(1,i)));
    end
    hold on
    if mod(iResult,3)==0
        plot(regionHist(1,:),regionHist(2,:),'r')
    end
    [expNum(iResult),a]=prepareCurveData(regionHist(1,:),regionHist(2,:));
    if mod(iResult,3)==0
        plot(regionHist(1,:),a*exp(-regionHist(1,:)*expNum(iResult)),'b')
    end
    warning off %#ok<WNOFF>
end
end
function [b,a]=prepareCurveData(xData,yData)

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
opts.StartPoint = [NaN -Inf];
opts.Upper = [Inf Inf];
opts.Robust = 'Bisquare';

% Fit model to data.
[fitresult, gof] = fit( xData', yData', ft, opts );

% Plot fit with data.
% disp(gof.rsquare);
b=abs(fitresult.b);
a=abs(fitresult.a);
% if gof.rsquare<=0.95
%     b=0;
% end
end