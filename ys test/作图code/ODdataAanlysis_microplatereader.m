function [amean,astd] = ODdataAanlysis_microplatereader(ODdata)
%ø����������ͼ������ODdata��ʱ�������зֲ�
F=10/60;
sam_rep=3;%ÿ����Ʒ�Ĳ��Դ��� sam ����sample
a=ODdata;
[xdim,ydim]=size(a);
num_sam=ydim/sam_rep;%��Ʒ����Ŀ
amean=zeros(xdim,num_sam);
astd=zeros(xdim,num_sam);
for isam=1:num_sam
    amean(:,isam)=mean(a(:,(isam-1)*sam_rep+1:isam*sam_rep),2);
    astd(:,isam)=std(a(:,(isam-1)*sam_rep+1:isam*sam_rep),0,2);
end

label={'PAO1','GacA','RetS','RsmA','RsmYZ'};

batchstdshade_ys(amean,astd,F,label);
figure,
h=line(0:F:(size(amean,1)-1)*F, amean,'LineWidth',2);
for i=1:size(amean,2)
    h(i,1).DisplayName=label{i};
end
ylabel('OD600');

% ���� xlabel
xlabel('Time (h)');

% ���� title
title('FAB 2.5mM Ca');


end


function batchstdshade_ys(amean,astd,F,label)
%���������ݵ�Ƚ���ʱ��������(std����sem)������Ӱ��,�д���ʵ���ظ�����
%����ԭ����stdshade_ys�޸�
%amean��һ��Ϊblank����
% amean=amean(:,2:end)-amean(:,1);%��ȥblank����
% astd=astd(:,2:end);
% label=label(2:end);

% amean=amean(:,1:end)-amean(:,1);
F=(0:size(amean,1)-1)*F; 
F=F';
figure;
h=line(F,amean,'lineWidth',2);

for i=1:size(amean,2)
    h(i,1).DisplayName=label{i};
    hold on;
    fill([F; flipud(F)],[amean(:,i)+astd(:,i); flipud(amean(:,i)-astd(:,i))],h(i,1).Color,'FaceAlpha', 0.4,'linestyle','none');
end

end