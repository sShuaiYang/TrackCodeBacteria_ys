function wells96_eachPrmtPlot(prmt_data,dirSave)
% 对与wells96显微镜获取的prmtdata 进行每个样品作图 
% 与eachSampleIntRatioPlot2.mlx类似
% 目前针对显微镜实验
% Shuai Yang 2022/10/10

fieldsName = fieldnames(prmt_data);
inFluoTF = contains(fieldsName,'int');
intFluo = fieldsName(inFluoTF);
chNum = sum(inFluoTF);
baseCh = intFluo{1}(4:end);% 内参荧光通道 一般默认第一个，可修改
expChs = cellfun(@(x) x(4:end),intFluo,'UniformOutput',false);

prmtNames = [prmt_data.promoterName];
prmtNum = numel(prmtNames);
for iPrmt = 1:prmtNum
    sampName = char(prmt_data(iPrmt).promoterName);
    fname = ['intChange','_',sampName];
    f1 = figure('Name',fname,'Position',[1632 790 1146 781]);  
    lightColor  = [0,0.33,0.74];
    darkColor  = [0.5,0.5,0.5];

    if any(prmt_data(iPrmt).wellsIllum == 0)
        pltcolor = darkColor;
    else
        pltcolor = lightColor;
    end

    pltNames = "sample"+string(prmt_data(iPrmt).wellsPosIdx);

    % fluoInt 作图
    for iCh = 1:chNum
        ax = nexttile;
        fluoInt = prmt_data(iPrmt).(intFluo{iCh});
        plt1 = plot(prmt_data(iPrmt).time_lag_c,fluoInt, ...
            'MarkerSize',4,'Marker','o','Color',pltcolor, ...
            'LineWidth',1);

        set(ax,'FontName','Arial','FontSize',10);
        box(ax,'off');
        title(ax,[sampName,'-',intFluo{iCh}],'Interpreter','none');
        xlabel(ax,'time(min)')
        ylabel(ax,[intFluo{iCh}, ' (a.u.)'])

        if iCh ==1
            legend(ax,'show');
            for iSamp = 1:numel(pltNames)
                set(plt1(iSamp),'DisplayName',pltNames(iSamp));
            end

        end        

    end

    % fluoInt Ratio 作图
    Ch1 = baseCh;
    for iCh = 1:chNum
        if strcmp(expChs{iCh}, baseCh)
            continue
        end
        Ch2 = expChs{iCh};
        intRatio = prmt_data(iPrmt).([Ch2,Ch1,'Ratio']);

        ax = nexttile;
        plt1 = plot(prmt_data(iPrmt).time_lag_c,intRatio, ...
            'MarkerSize',4,'Marker','o','Color',pltcolor, ...
            'LineWidth',1);

        set(ax,'FontName','Arial','FontSize',10);
        box(ax,'off');
        title(ax,[sampName,'-',Ch2,Ch1,'Ratio'],'Interpreter','none');
        xlabel(ax,'time(min)')
        ylabel(ax,[Ch2,Ch1,'Ratio', ' (a.u.)'])
    end
    saveas(f1,[dirSave,'\',fname,'.fig']); 
    close(f1)
  
end

end