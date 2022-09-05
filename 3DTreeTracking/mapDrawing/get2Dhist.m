function [countMap,xbin,ybin]=get2Dhist(data1,data2,bin1,bin2)
% 二维统计时获得surface颜色统计图的程序
xbin=bin1(2:2:end);
ybin=bin2(2:2:end);
countMap=zeros((size(bin2,2)-1)/2,(size(bin1,2)-1)/2);
for i=2:2:size(bin1,2)
    dataCol=data2( data1>=bin1(i-1)& data1<bin1(i+1));
    for j=2:2:size(bin2,2)
        countMap(j/2,i/2)=numel(dataCol(dataCol>=bin2(j-1)&dataCol<bin2(j+1)));
    end
end
end