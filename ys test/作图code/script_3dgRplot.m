%  
% label={'PAO1','GacA','RetS','RsmA','RsmYZ'};
% 
% figure,scatter3(gr2(1:5,1),gr2(1:5,2),gr2(1:5,3),'filled')
% hold on,scatter3(gr2(6:10,1),gr2(6:10,2),gr2(6:10,3),'filled')
% hold on,scatter3(gr2(11:15,1),gr2(11:15,2),gr2(11:15,3),'filled')
% for j=1:5
% for i=0:2
%     text(gr2(j+5*i,1),gr2(j+5*i,2),gr2(j+5*i,3),label{j})
% end
% end

 
label={'PAO1','GacA','RetS','RsmA','RsmYZ'};

hold on,scatter3(gr2(1:5,1),gr2(1:5,2),gr2(1:5,3),'filled')
hold on,scatter3(gr2(6:10,1),gr2(6:10,2),gr2(6:10,3),'filled')
hold on,scatter3(gr2(11:15,1),gr2(11:15,2),gr2(11:15,3),'filled')
for j=1:5
for i=0:2
    text(gr2(j+5*i,1),gr2(j+5*i,2),gr2(j+5*i,3),label{j})
end
end