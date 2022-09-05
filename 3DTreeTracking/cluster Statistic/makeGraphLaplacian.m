function [pResult,laplacianMatrix]=makeGraphLaplacian(adjacentMatrix)
adjacentMatrix=full(adjacentMatrix);
laplacianMatrix=zeros(size(adjacentMatrix));
laplacianMatrix(adjacentMatrix==1)=-1;
for i=1:size(laplacianMatrix,1)
    iLine=adjacentMatrix(i,:);
    laplacianMatrix(i,i)=numel(iLine(iLine==1));
end
laplacianMatrix=normalizeLaplacian(laplacianMatrix);
eigenvalue=eig(double(laplacianMatrix));
[x,xResult]=gaussianConvolution(eigenvalue);
pResult=cat(2,x',xResult');
plot(x,xResult);
set(gca,'xLim',[0,2]);
end
function laplacianMatrixNew=normalizeLaplacian(laplacianMatrix)
laplacianMatrixNew=zeros(size(laplacianMatrix));
matrixSize=size(laplacianMatrix,1);
for i=1:matrixSize
    for j=1:matrixSize
        if i==j && laplacianMatrix(i,i)~=0
            laplacianMatrixNew(i,j)=1;
        else
            if i~=j && laplacianMatrix(i,i)~=0 && laplacianMatrix(j,j)~=0 && laplacianMatrix(i,j)==-1
            laplacianMatrixNew(i,j)=-(laplacianMatrix(i,i)*laplacianMatrix(j,j))^-(0.5);
            end
        end
    end
end
end
function createSpectralDensity(laplacianMatrix)
eigenvalue=eig(double(laplacianMatrix));
p=linspace(0,2,201);
interval=p(3)-p(1);
pResult=[];
for iNum=2:2:numel(p)-1
    pHist=numel(eigenvalue(eigenvalue>=p(iNum-1)&eigenvalue<p(iNum+1)));
    pResult=[pResult;p(iNum),pHist/interval];
end
plot(pResult(:,1),pResult(:,2))
end
function pResult=countHistNormal(laplacianMatrix)
% allNum=unique(originalData);
% pResult=[];
% for iNum=1:numel(allNum)
%     pHist=numel(originalData(originalData==allNum(iNum)));
%     pResult=[pResult;allNum(iNum),pHist];
% end
p=linspace(-4*constFix,6*constFix,401);
interval=p(3)-p(1);
pResult=[];
for iNum=2:2:numel(p)-1
    pHist=0;
    for i=1:numel(originalData)
        if originalData(i)>=p(iNum-1) && originalData(i)<p(iNum+1)
            pHist=pHist+1;
        end
    end
    pResult=[pResult;p(iNum),pHist/interval];
end
end
function [x,xResult]=gaussianConvolution(eigenvalue)
x=-1:0.001:3;
parfor i=1:numel(eigenvalue)
    sigma=0.0075;
   gaussianResult{i}=1/((2*pi)^0.5*sigma)*exp(-(x-eigenvalue(i)).^2/2/sigma^2);
end
xResult=[];
for i=1:numel(eigenvalue)
    if i==1
        xResult=gaussianResult{i};
    else
    xResult=xResult+gaussianResult{i};
    end
end
xResult=xResult/numel(eigenvalue);
end
    