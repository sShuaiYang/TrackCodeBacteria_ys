function batchscatterhistForcxydata(dataAll)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
dirSave='E:\cxy20181113校正并删选\scatterhist';
dataNum=size(dataAll,2);
for i=1:dataNum
    figure,h=scatterhist(dataAll{i}(:,1),dataAll{i}(:,2),'Kernel','on');ylabel('mScarlet');xlabel('sfGFP');
    saveas(h,strcat(dirSave,'\',num2str(i,'%02d'),'.fig'));
end
close all;
end

