function [pedigreeAll] = genealogicalTreePlot_3D(pedigreeAll)
%xy是细菌质心在视野中的真实位置，z轴为时间
% shuai Yang 2020.05.10
% dirSave = 'C:\Users\XJY\Desktop\genealogicalTreeTest';
for i = 1:numel(pedigreeAll)
    [leafLabels] = leafLabelsGeneration(pedigreeAll,i);
    
    for iLeaf = 1:pedigreeAll(i).leafNum
        if pedigreeAll(i).tree(iLeaf).label==1 % 如果label等于1 root标签
            data_3D = zeros(numel(pedigreeAll(i).tree(iLeaf).timer),3);
            data_3D(:,3)= pedigreeAll(i).tree(iLeaf).timer;%z data
            for iTime = 1:numel(pedigreeAll(i).tree(iLeaf).timer)
                data_3D(iTime,1:2) = ...
                    pedigreeAll(i).tree(iLeaf).traceInfo.measurment{iTime}.Centroid;
                %质心的xypoisiton
            end
        else
            %不等于1，data_3D第一个列是mothercell形成node前的最后位置，加进去 方便作图时连接
            data_3D = zeros(numel(pedigreeAll(i).tree(iLeaf).timer)+1,3);
            data_3D(2:end,3)= pedigreeAll(i).tree(iLeaf).timer;%z data
            for iTime = 1:numel(pedigreeAll(i).tree(iLeaf).timer)
                data_3D(iTime+1,1:2) = ...
                    pedigreeAll(i).tree(iLeaf).traceInfo.measurment{iTime}.Centroid;
                %质心的xypoisiton
            end
            motherLabel = fix(pedigreeAll(i).tree(iLeaf).label/2);%找到其母代细菌的label
            posIdx = find(leafLabels == motherLabel);%根据label找位置
            data_3D(1,:) = pedigreeAll(i).tree(posIdx).data_3D(end,:);
        end
        pedigreeAll(i).tree(iLeaf).data_3D =data_3D;
    end
    
end
% treePlot_xyt(pedigreeAll,dirSave);
end

%%
function [leafLabels] = leafLabelsGeneration(pedigreeAll,i)
leafLabels = zeros(pedigreeAll(i).leafNum,1);
for iLeaf = 1:pedigreeAll(i).leafNum
    leafLabels(iLeaf)=pedigreeAll(i).tree(iLeaf).label;
end
end
%xyt tree plot
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