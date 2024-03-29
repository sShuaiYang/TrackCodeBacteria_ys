function pResult=makeAdjacentMatrix(distMatrix)
eigenvalue=eig(double(full(distMatrix)));
nodeNum=size(distMatrix,1);
possibleLink=(numel(distMatrix(distMatrix==1))/2)/(nodeNum*(nodeNum-1)/2);
constFix=(possibleLink*(1-possibleLink)*nodeNum).^0.5;
pResult=gaussianConvolution(eigenvalue,constFix);
standardX=-4*constFix:0.01:6*constFix;
standardY=(4*nodeNum*possibleLink*(1-possibleLink)-standardX.^2).^0.5/2/pi/(nodeNum*possibleLink*(1-possibleLink));
h=SpectralDensity(pResult(:,1)/constFix,pResult(:,2)*constFix,standardX/constFix,standardY*constFix);
pResult(:,1)=pResult(:,1)/constFix;
pResult(:,2)=pResult(:,2)*constFix;
end
function pResult=gaussianConvolution(eigenvalue,constFix)
x=-4*constFix:0.001*constFix:6*constFix;
parfor i=1:numel(eigenvalue)
    sigma=0.025;
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
pResult=cat(2,x',xResult');
end
function h=SpectralDensity(X1, Y1, X2, Y2)
%CREATEFIGURE(X1,Y1,X2,Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 30-Nov-2012 22:39:08

% Create figure
figure1 = figure();

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'Color',[0 0 1]);

% Create xlabel
xlabel('��*��Np(1-p)��^0.5','FontSize',16);

% Create ylabel
ylabel('��*��Np(1-p)��^0.5','FontSize',16);

% Create plot
plot(X2,Y2,'LineStyle','--','Color',[1 0 0]);
h=gcf;
end