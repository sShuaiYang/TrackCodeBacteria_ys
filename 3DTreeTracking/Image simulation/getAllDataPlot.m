function getAllDataPlot()
dirFile=uigetdir();
nameList=dir(dirFile);
colorAll=colormap(spring(11));
colorAll=[colorAll;0.48,0.06,0.89;0.48,0.06,0.89];
for i=1:numel(nameList)-2
    if i<=11
        dirFrameInfo=strcat(dirFile,'\',nameList(i+2).name,'\bacteriaFrameInfo');
        load(dirFrameInfo)
        [degreeResult] = mixNumPlot(bacteriaFrameInfo);
        load(strcat(dirFile,'\',nameList(i+2).name,'\cosSitaAndLength'))
        b=curveFitting(cosSitaAndLength(:,1),cosSitaAndLength(:,2));
        hold on;plot(degreeResult(:,1),degreeResult(:,2),'DisplayName',num2str(b),'Color',colorAll(i,:),'marker','none','lineStyle','--','lineWidth',3);
    else
        dirFrameInfo=strcat(dirFile,'\',nameList(i+2).name,'\bacteriaFrameInfo');
        load(dirFrameInfo)
        [degreeResult] = mixNumPlot(bacteriaFrameInfo);
        if i==12
            degreeResult(degreeResult(:,1)>35,:)=[];
        end
        if i==13
            degreeResult(degreeResult(:,1)>48,:)=[];
        end
        hold on;plot(degreeResult(:,1),degreeResult(:,2),'DisplayName',num2str(i),'Color',colorAll(i,:),'marker','.','lineWidth',3);
    end
end
end
function b=curveFitting(xdata,ydata)

% Set up fittype and options.
ft = fittype( 'exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = -Inf;
opts.StartPoint = 0.970592781760616;
opts.Upper = Inf;

% Fit model to data.
fitresult= fit(xdata,ydata,ft,opts);
b=fitresult.b;
end
function [degreeResult] = mixNumPlot(bacteriaFrameInfo)
% generate direct link matrix
% aveLen=40;
multi=3;
aimFrame=numel(bacteriaFrameInfo);
finalInfo=bacteriaFrameInfo{aimFrame};
aveLen=50;
distMatrix=pdist2(finalInfo.centroidInfo,finalInfo.centroidInfo);

distMatrix(distMatrix<=aveLen*multi & distMatrix~=0)=1;
distMatrix(distMatrix~=1)=0;
distMatrix=logical(distMatrix);

% branch mix
degree=[];
branchIndex=bacteriaFrameInfo{aimFrame}.bacteriaInfo(:,5);
for i=1:size(distMatrix,2)
    iLine=distMatrix(i,:);
    degree(i,1)=numel(iLine(iLine==1));
    attachBranch=branchIndex(iLine);
    attachBranch=unique(attachBranch);
    attachNum(i,1)=numel(attachBranch);
end
attachNum=attachNum./degree;
degreeResult=[];
for i=1:max(degree)
    degreeResult=[degreeResult;[i,mean(attachNum(degree==i))]];
end
end

