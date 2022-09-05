function plotCorrelationMap(allData)
clc
try
    nameList=dir(dirFile);
catch error
    disp('please run 1.getExperimentFile');
    return
end
disp('please input treePlot type')
disp('1.CyOFP');
disp('2.sfGFP');
disp('3.mScaletI');
disp('4.iRFP');
treeType=input('____:');
try
    switch treeType
        case 1
            result=allData.correlationDataCYOFP;
        case 2
            result=allData.correlationDataGFP;
        case 3
            result=allData.correlationDatamScalet;
        case 4
            result=allData.correlationDataRFP;
    end
    result(isnan(result(:,2)),:)=[];
    generation=unique(result(:,1));
    for i=1:numel(generation)
        finalResult(i,1)=generation(i);
        corrNum=result(result(:,1)==generation(i),2);
        finalResult(i,2)=mean(corrNum);
        finalResult(i,3)=var(corrNum);
    end
    errorbar(finalResult(:,1),finalResult(:,2),finalResult(:,3));
    xlim(gca,[0 max(finalResult(:,1))+1]);
    prepareCurveFitData(finalResult(:,1),finalResult(:,2))
catch err
    disp('do not have this fluo color')
end
end
function [xData, yData] = prepareCurveFitData(a,b)

% Set up fittype and options.
ft = fittype( 'exp(-x/a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = -Inf;
opts.Robust = 'Bisquare';
opts.StartPoint = 3;
opts.Upper = Inf;

% Fit model to data.
[fitresult, gof] = fit(a,b, ft, opts );

% Plot fit with data.
hold on
h = plot( fitresult,a,b);
text(6,0.9,['generation=',num2str(fitresult.a)]);
end