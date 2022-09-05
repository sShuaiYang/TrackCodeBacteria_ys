function resultImage=easyHist(originalData,plotColor)
% p=[1:9,ceil(logspace(1,log10(max(originalData)),20))];
% p=ceil(linspace(1,400,400));
% p=1:max(originalData);
% pResult=[];
% for iNum=1:(numel(p)-1)
%     if iNum<=10
%         pResult=[pResult;iNum,numel(originalData(originalData==iNum))];
%     else
%         pResult=[pResult;(p(iNum)+p(iNum+1))/2,numel(originalData(originalData>=p(iNum) & originalData<p(iNum+1)))/(p(iNum+1)-p(iNum))];
%     end
% end
% pResult(:,2)=pResult(:,2)/max(pResult(:,2));
resultImage=oneDimGaussian(originalData,2,2);   % for ratio
% resultImage=oneDimGaussian(originalData,0.3,0.3);
hold on;plot(resultImage(end/2:end,1),resultImage(end/2:end,2),'Color',plotColor)
% % hold on;plot(pResult(:,1),pResult(:,2),'Color',plotColor,'LineStyle','none','Marker','o')
end
function resultImage=oneDimGaussian(degreeResult,interval,sigma)
xMin=0;
xMax=max(degreeResult);
rowPos=-xMax:interval:xMax;
rowPos=rowPos';
resultImage=zeros(numel(rowPos),2);
for i=1:numel(degreeResult)
%     disp(i)
    rowGauss=1/(2*pi)^0.5/sigma*exp(-(rowPos-degreeResult(i)).^2/2/sigma^2)+1/(2*pi)^0.5/sigma*exp(-(rowPos+degreeResult(i)).^2/2/sigma^2);
    resultImage(:,2)=resultImage(:,2)+rowGauss;
end
resultImage(:,1)=rowPos;
sumAll=sum(resultImage(:,2)*interval);
resultImage(:,2)=resultImage(:,2)/sumAll;
end