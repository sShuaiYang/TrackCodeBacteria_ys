function [xNew,yNew]=curveGuassionConvert(x,y)
xNew=-min(x):0.01:max(x);
parfor i=1:numel(x)
    sigma=0.1;
   gaussianResult{i}=y(i)/((2*pi)^0.5*sigma)*exp(-(xNew-x(i)).^2/2/sigma^2);
end
yNew=[];
for i=1:numel(x)
    if i==1
        yNew=gaussianResult{i};
    else
    yNew=yNew+gaussianResult{i};
    end
end
yNew=yNew/numel(x)^0.5;
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