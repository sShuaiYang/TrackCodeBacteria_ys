function pairCellDistanceHistPlot(gRPos_class,dirSave)
%函数gRAndDistanceAnalysis获得gRPos_class后计算
%距离不同距离的两两细菌（pair）的数目 对 r作图
% plot function
for i=1:numel(gRPos_class)
    if ~isfield (gRPos_class{i},'plotData')
        continue
    end
    distCut=50;%distance less than 50um
    distLogical=gRPos_class{i}.plotData(:,4)<distCut;
    if sum(distLogical)>50 %在小于50um的前提下判断数据点数目
        
        xData = gRPos_class{i}.plotData(distLogical,4);%distance

        nbin = 100;%distCut/nbin = 0.5； 相当于0.5um的距离间隔
        
        figure,
%         [N,edges] = histcounts(xData,nbin);
        z = zeros(2,nbin-1);
        %         z(2,:) = N;
        edges = linspace(min(xData),max(xData),nbin);
        
        for iBin = 1:numel(edges)-1
            z(1,iBin) = (edges(iBin)+edges(iBin+1))/2;
            z(2,iBin) = sum(xData>=edges(iBin)&xData<edges(iBin+1));
        end
        plot(z(1,:),z(2,:),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','--','Color',[1 0 0])
        xlabel('r (um)');
        ylabel('cell pairs');
        title(['Generation',num2str(i)]);
        fname1 = ['g',num2str(i),'_distHist'];
        set(gca,'FontSize',12,'YScale','linear');
%         set(gca,'FontSize',12,'YScale','log');
        saveas(gcf,[dirSave,'\',fname1,'.fig']);
        saveas(gcf,[dirSave,'\',fname1,'.tif']);
    end
end
close all
end
