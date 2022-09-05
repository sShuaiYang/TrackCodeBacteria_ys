function finalMatrix=paraMeterGet(bioTree,gfpImage,rfpImage)
alphaG=0.4722*(1:1:10);   % for SFGFP
alphaR=0.0534*(1:1:10);
for i=1:numel(alphaG)
    for j=1:numel(alphaR)
        paraG(i,j)=alphaG(i);
        paraR(i,j)=alphaR(j);
        [gamaAll,miuAll,fpProduceAll]=calculateMiuAndProduceRate(bioTree,gfpImage,rfpImage,alphaG(i),alphaR(j));
        [slope,rsquare] = prepareCurveData( gamaAll', miuAll' );
        slopeAll(i,j)=slope;
        rsquareAll(i,j)=rsquare;
    end
end
finalMatrix.paraG=paraG;
finalMatrix.paraR=paraR;
finalMatrix.slopeAll=slopeAll;
finalMatrix.rsquareAll=rsquareAll;
end

function [slope,rsquare] = prepareCurveData( gamaAll, miuAll )

% Set up fittype and options.
ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = -Inf;
opts.Robust = 'Bisquare';
opts.StartPoint = 0.950222048838355;
opts.Upper = Inf;

% Fit model to data.
[fitresult, gof] = fit( gamaAll, miuAll, ft, opts );

slope=fitresult.a;
rsquare=gof.rsquare;
end


