function batchscatterhistForcxydata(dataAll)
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
dirSave='E:\cxy20181113У����ɾѡ\scatterhist';
dataNum=size(dataAll,2);
for i=1:dataNum
    figure,h=scatterhist(dataAll{i}(:,1),dataAll{i}(:,2),'Kernel','on');ylabel('mScarlet');xlabel('sfGFP');
    saveas(h,strcat(dirSave,'\',num2str(i,'%02d'),'.fig'));
end
close all;
end

