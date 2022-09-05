function  genealogicalTreePlot_MultiType(pedigreeAll,dirSave,type)
%生成pedigreeAll后 用户自定义选择做那种类型的图
% Shuai Yang 2020.05.12
%type = 1， 3D 
%type = 2， line
%type = 3， arc %默认
%type = 0， 3种都做 

if nargin < 3 || isempty(type)
    type = 3;
end

switch type
    case 1
        treePlot_xyt(pedigreeAll,dirSave);
    case 2
        treePlot_lineLink(pedigreeAll,dirSave);
    case 3
        treePlot_arcLink(pedigreeAll,dirSave);
    case 0
        treePlot_xyt(pedigreeAll,dirSave);
        treePlot_lineLink(pedigreeAll,dirSave);
        treePlot_arcLink(pedigreeAll,dirSave) ;
end

end

%% xyt tree plot
function treePlot_xyt(pedigreeAll,dirSave)
for i = 1:numel(pedigreeAll)
    fname = strcat('pedigree3D_#',num2str(pedigreeAll(i).frameIdx),'-',num2str(pedigreeAll(i).rootIdx));
    close all
    f1 = figure('Name',fname);
    % set(f1,'visible','off')
    axes1 = axes('Parent',f1);
    % hold(axes1,'on');
    for iLeaf = 1:pedigreeAll(i).leafNum
        xyt = pedigreeAll(i).tree(iLeaf).data_3D;
        %有时候bioTree 用1,1,2,2...这样的重复image生成的，每一frame都重复了
        if isequal(pedigreeAll(i).tree(1).data_3D(1,3),pedigreeAll(i).tree(1).data_3D(2,3))
            xyt = xyt(1:2:end,:);
        end
        plot3(xyt(:,1),xyt(:,2),xyt(:,3),'Color',[0,0.45,0.74],'LineWidth',1)
        hold on
    end
    ylabel({'YPos'});
    xlabel({'XPos'});
    zlabel({'time (min)'});
    title(fname,'Interpreter','none');
    
    % hold(axes1,'off')
    hold off
    set(axes1,'FontSize',12);
    
    saveas(f1,[dirSave,'\',fname,'.fig'])
    close all
end

end
%% line link plot function
function treePlot_lineLink(pedigreeAll,dirSave)
for i = 1:numel(pedigreeAll)
    fname = strcat('pedigreeLine_#',num2str(pedigreeAll(i).frameIdx),'-',num2str(pedigreeAll(i).rootIdx));
    close all
    f1 = figure('Name',fname);
    % set(f1,'visible','off')
    axes1 = axes('Parent',f1,'Position',[0.13 0.11 0.5 0.8]);
    % hold(axes1,'on');
    for iLeaf = 1:pedigreeAll(i).leafNum
        linePts = pedigreeAll(i).tree(iLeaf).linePts;
        line(linePts(:,1),linePts(:,2),'Color',[0,0.45,0.74],'LineWidth',1)
        hold on
    end
    ylabel({'time (min)'});
    title(fname,'Interpreter','none');
    
    % hold(axes1,'off')
    hold off
    set(axes1,'FontSize',12,'XTick',zeros(1,0));
    
    saveas(f1,[dirSave,'\',fname,'.fig'])
    close all
end


end
%% arc link plot function
function treePlot_arcLink(pedigreeAll,dirSave)
for i = 1:numel(pedigreeAll)
    fname = strcat('pedigreeArc_#',num2str(pedigreeAll(i).frameIdx),'-',num2str(pedigreeAll(i).rootIdx));
    close all
    
    f1 = figure('Name',fname);
    % set(f1,'visible','off')
    axes1 = axes('Parent',f1,'Position',[0.13 0.11 0.5 0.8]);
    % hold(axes1,'on');
    for iLeaf = 1:pedigreeAll(i).leafNum
        arcLinkPts = pedigreeAll(i).tree(iLeaf).arcLinkPts;
        line(arcLinkPts(:,1),arcLinkPts(:,2),'Color',[0,0.45,0.74],'LineWidth',1)
        hold on
    end
    ylabel({'time (min)'});
    title(fname,'Interpreter','none');
    % hold(axes1,'off')
    hold off
    set(axes1,'FontSize',12,'XTick',zeros(1,0));
    saveas(f1,[dirSave,'\',fname,'.fig'])
    close all
end
end