templogic=zeros(1,length(bacInfo));
k=1;
bacInfo_new={};
for i=1:length(bacInfo)

    if isfield(bacInfo{1,i},'meanGFP')
        bacInfo_new{1,k}=bacInfo{1,i};
        k=k+1;
    end

end

figure,
for i=1:length(bacInfo_new)
     line('XData',repmat(bacInfo_new{1, i}.time,1, length(bacInfo_new{1, i}.meanCyOFP))/60,'YData',bacInfo_new{1, i}.meanCyOFP,...
         'MarkerSize',3,'Marker','o','LineStyle','none','Color',[0.850980392156863 0.329411764705882 0.101960784313725]);
     xlim ([0,12])
     hold on
end

figure,
for i=1:length(bacInfo_new)
     line('XData',bacInfo_new{1, i}.time/60,'YData',mean(bacInfo_new{1, i}.meanGFP),...
         'MarkerSize',3,'Marker','o','LineStyle','none','Color',[0.850980392156863 0.329411764705882 0.101960784313725]);
     xlim ([0,12])
     hold on
end

figure,
sfGFP=[];
for i=1:length(bacInfo_new)
    sfGFP(i,1)=bacInfo_new{1, i}.time/60;
    sfGFP(i,2)=mean(bacInfo_new{1, i}.meanGFP);
    sfGFP(i,3)=std(bacInfo_new{1, i}.meanGFP);
end
% sfGFP=sfGFP';
% fill([sfGFP(1,:), fliplr(sfGFP(1,:))],[sfGFP(2,:)+sfGFP(3,:), fliplr(sfGFP(2,:)-sfGFP(3,:))], [0.850980392156863 0.329411764705882 0.101960784313725], 'FaceAlpha', 0.3,'linestyle','none');
errorbar(sfGFP(:,1),sfGFP(:,2),sfGFP(:,3))

mScarlet=[];
for i=1:length(bacInfo_new)
    mScarlet(i,1)=bacInfo_new{1, i}.time/60;
    mScarlet(i,2)=mean(bacInfo_new{1, i}.meanmScalet);
    mScarlet(i,3)=std(bacInfo_new{1, i}.meanmScalet);
end
% mScarlet=msmScarlet(1:25,:);
a3=mScarlet(1:25,2);

CyOFP=[];
for i=1:length(bacInfo_new)
    CyOFP(i,1)=bacInfo_new{1, i}.time/60;
    CyOFP(i,2)=mean(bacInfo_new{1, i}.meanCyOFP);
    CyOFP(i,3)=std(bacInfo_new{1, i}.meanCyOFP);
end
% mScarlet=msmScarlet(1:25,:);
a4=CyOFP(1:25,2);
Venus=[];
for i=1:length(bacInfo_new)
    Venus(i,1)=bacInfo_new{1, i}.time/60;
    Venus(i,2)=mean(bacInfo_new{1, i}.meanVenus);
    Venus(i,3)=std(bacInfo_new{1, i}.meanVenus);
end
% mScarlet=msmScarlet(1:25,:);
a5=Venus(1:25,2);

figure,
 for row = 1:4
    for col = 1:4
        nm = 4*(row-1)+col;
        subplot(4,4,nm)
        stem(lgs/2,cr(:,nm),'o')
        title(sprintf('c_{%d%d}',row,col))
        ylim([0 1])
    end
 end

figure,
 for row = 1:4
    for col = 1:4
        nm = 4*(row-1)+col;
        subplot(4,4,nm)
        stem(lgs/2,cr(:,nm),'o')
        title(sprintf('c_{%d%d}',row,col))
        ylim([0 1])
    end
 end

%  figure
%  for row = 1:4
%     for col = 1:4
%         nm = 4*(row-1)+col;
%         subplot(4,4,nm)
%         stem(lgs1,cr1(:,nm),'LineStyle','-','MarkerSize',5)
%         title(sprintf('c_{%d%d}',row,col))
%         ylim([0 1])
%         xlim([-12,12])
%         box off
%     end
%  end
% 
% figure
% for row = 1:4
%     for col = 1:4
%         nm = 4*(row-1)+col;
%         subplot(4,4,nm)
%         line([0,0],[0,1],'color',[0.6,0.6,0.6])
%         hold on,
%         line(lgs2/2,cr2(:,nm),'MarkerSize',7,'Marker','o','LineStyle','none')
%         title(sprintf('c_{%d%d}',row,col))
%                 ylim([0 1])
%         xlim([-8,8])
%     end
% end



