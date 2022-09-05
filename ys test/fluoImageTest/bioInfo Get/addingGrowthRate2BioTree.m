function [bioTree] = addingGrowthRate2BioTree(bioTree)
%将每个细菌的生长率添加到bioTree中
%处理bioTree中的每个root和node
% Shuai Yang  2020.05.21
for iFrame = 1:numel(bioTree)
    
    for iRoot = 1:numel(bioTree{iFrame}.root)% 如果root为空，numel为0 for循环不执行
        bioTree{iFrame}.root{iRoot}.gR = ones(1,2)*NaN;
        traceInfo = bioTree{iFrame}.root{iRoot}.traceInfo; 
        traceNum = numel(traceInfo.measurment);
        timer = bioTree{1, 1}.bioTreeTimer(iFrame:iFrame+traceNum-1);
        if traceNum >= 3
            cellArea = zeros(1,traceNum);
            cellLength = zeros(1,traceNum);
            for iTrace = 1: traceNum
                cellArea(iTrace) = traceInfo.measurment{iTrace}.FilledArea;
                cellLength(iTrace) = traceInfo.measurment{iTrace}.MajorAxisLength;
            end 
             [gR,~,fitgof] = growthRatefit_poly1(timer, cellLength);
             if gR>0.001 && gR<0.1 && fitgof.rsquare >= 0.8 % gR threshold 0.003或者0.001
                 bioTree{iFrame}.root{iRoot}.gR(1) = gR;%第一列生长 min-1
                 bioTree{iFrame}.root{iRoot}.gR(2) = fitgof.rsquare;%第二列拟合的rsquare
             end
        end
    end
    
    for iNode = 1:numel(bioTree{iFrame}.node)
        for iOut = 1:numel(bioTree{iFrame}.node{iNode}.Out)
            bioTree{iFrame}.node{iNode}.Out{iOut}.gR = ones(1,2)*NaN;
            traceInfo = bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo;
            traceNum = numel(traceInfo.measurment);
            timer = bioTree{1, 1}.bioTreeTimer(iFrame:iFrame+traceNum-1);
            
            if traceNum >=3
                cellArea = zeros(1,traceNum);
                cellLength = zeros(1,traceNum);
                for iTrace = 1: traceNum
                    cellArea(iTrace) = traceInfo.measurment{iTrace}.FilledArea;
                    cellLength(iTrace) = traceInfo.measurment{iTrace}.MajorAxisLength;
                end
                [gR,~,fitgof] = growthRatefit_poly1(timer, cellLength);
                if gR>0.003 && gR<0.1 && fitgof.rsquare>=0.8
                    bioTree{iFrame}.node{iNode}.Out{iOut}.gR(1) = gR;
                    bioTree{iFrame}.node{iNode}.Out{iOut}.gR(2) = fitgof.rsquare;
                end
            end
            
        end
        
    end
end

end

%% growth rate fit
function [gR,fitperiodata,gof] = growthRatefit_poly1( t, dataPickI )
% Shuai Yang 2021/01/10
% ln(L) = ln(L0) + ut;
% 用 y = p1*x+p2 进行拟合生长率 u = p1
% 初始t不一定需要为0 影响p1值
% dL/L = ut;L = L0*e^(ut);
% ln(L) = ln(L0) + ut;
% (L2-L1)/L1 = e^(ut2-ut1)-1;
% e^u(t2-t1) = L2/L1;u = ln(L2/L1)/(t2-t1);
% doubling time DT = log(2)/u;
% dataPickI 为细菌长度或者细菌面积随时间的变化
% timeSeries 为时间序列 hr /min /s
% 也可参考 growthRatefit_exp1
% row number
t = t(:);
dataPickI = dataPickI(:);
% 取对数后 smooth
dataPickI = log(dataPickI);
dataPickI = smooth(dataPickI,'lowess');
dataPickI = smooth(dataPickI,'rlowess');

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Robust = 'Bisquare';
opts.Upper = [Inf Inf];
warning('off');
% Fit model to data.
[fitperiodata,gof] = fit(t,dataPickI,ft,opts);
gR = fitperiodata.p1;
end