function [processedImages] = phaseContrastImageProcessing_ys(phaseImages)
%phase Contrast image processing code
% phase contrast 16bit图像 转化为8bit图像
%图像为black图像 转化为white图像
% phase Contrast图像比正常明场图像细菌胖一些
%phaseContrastImageProcessingTest2ys 进行测试效果
%Shuai Yang 2020.04.09
%phaseContrastImageProcessingTest3ys 重新进行测试效果 2020.04.17

convertedImages = uint8 (zeros (size(phaseImages)));
parfor iframe = 1:size(phaseImages,3)
    I0 = phaseImages(:,:, iframe);
    I0 = uint8(double(I0-min(I0(:)))/double(max(I0(:))-min(I0(:)))*255);
    I = imcomplement(I0);% equal to I = 255 - I0;
    I = imsubtract(I,I0);
    convertedImages (:,:,iframe) = I;
end
phaseImages = convertedImages;

processedImages = false(size(phaseImages));

areaThreshold = 250;
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3) = 24;%created edgefilter

parfor iframe = 1:size(phaseImages,3)
    I = phaseImages (:,:,iframe);
    I = imgaussfilt(I, 2);
    I = imfilter(I,edgeFilter);
    % edgefilter后会产生一个2个pixel亮的边框
    I_border = ones(size(I))*255; I_border(3:size(I,1)-2,3:size(I,1)-2) = 0;
    I_border = uint8(I_border);
    I = imsubtract(I,I_border);
    I = imgaussfilt(I,2);
    [Iobrcbr] = cleanupImageUsingReconstuction_PhCImages(I);
    
    tempBW = edge((Iobrcbr),'log', 1/255, 1);%1/255,threshold;1,sigma.
    tempBW = bwareaopen( tempBW, 20);% 删除小的线段
    % tempBW = bwmorph(tempBW,'bridge');% 不要 会出现很多小的错误   
    [tempBW] = circleLeakingEdgeContours(tempBW);
    
    tempBW = imclearborder(tempBW);
    tempBW = imfill( tempBW,'holes');
    tempBW = imopen( tempBW, strel('disk',2));
    tempBW = bwareaopen( tempBW, areaThreshold);
    % imshow(labeloverlay(I1, tempBW))
    % title('Iobrcbr mask')
%     [tempBW] = bwModification(tempBW);
    %     imshow(labeloverlay(I1, tempBW))
    %     title('Iobrcbr mask edge')
    tempBW = bwareaopen( tempBW, areaThreshold);
    processedImages(:,:,iframe) = tempBW;
end

end
%%

function [bwSegNew]= bwModification(tempBW)
%对细菌的mask进行一些形态学处理
% 腐蚀 骨架 删除endpoints 
%再用骨架进行重构
tempBW = imopen(tempBW,strel('disk',2));
tempBW = imopen(tempBW,ones(3));
tempBW = bwmorph(tempBW,'open');

tempBW = imerode(tempBW,ones(4));%腐蚀进行分割
tempBW = bwareaopen(tempBW, 5);%删除小的零星的pixel
skel = bwmorph(tempBW,'skel',Inf);%获取骨架
E = bwmorph(skel, 'endpoints');
skel = skel&(~E); %删除endpoint 要不膨胀会连在一起
E = bwmorph(skel, 'endpoints');
skel = skel&(~E); %删除两次
E = bwmorph(skel, 'endpoints');
skel = skel&(~E); %删除三次
% E = bwmorph(skel, 'endpoints');
% skel = skel&(~E); %删除四次
%重组细菌
se = strel('disk',5);
bwSegNew = imdilate(skel,se);
bwSegNew=imfill( bwSegNew,'holes');

% bwSegNew=bwmorph(bwSegNew,'hbreak');
% bwSegNew=imopen(bwSegNew,strel('disk',2));
% bwSegNew=imopen(bwSegNew,ones(3));
% bwSegNew=bwmorph(bwSegNew,'open');

end
%%
function [BW] = circleLeakingEdgeContours(BW)
%对edge的图像进行处理 将距离<2pixeld的endpoints 连起来形成闭合的图形
distThresh = 5;%大于3个pixel的对角线距离
E = bwmorph(BW, 'endpoints');
[x,y] = find(E);
C = [x,y];%endpoints 的xy位置
tempLogic=C(:,1)>5&C(:,1)<2043&C(:,2)>5&C(:,2)<2043;%不管视野边框周围的细菌[5,2043]
C=C(tempLogic,:);
if size(C,1)<2 %如果endpoint的数目小于2不再计算
    return
end
D=pdist(C);
Z = squareform(D);
pointIdx = true(1,size(C,1));
for iPoint = 1: size(C,1)
    if pointIdx(iPoint) == true
        [~, IX] = sort(Z(:,iPoint));
        k = IX(2);% 找到第二小的距离值的索引,第一是自己距离为0
        if Z(k,iPoint) <= distThresh && pointIdx(k) == true
            pt1 = C(iPoint,:);pt2 = C(k,:);
            E = connectTwopointsInImage (E,pt1,pt2);
            pointIdx(k) = false;
            %一个点只连接距离最近的那个点 连接后 这两个点都不在进入库进行下一次判断
        end
        pointIdx(iPoint) = false; %已经循环过的points不再重复判断
    end
    
end
BW = BW|E;
end
function [E] = connectTwopointsInImage (E,pt1,pt2)

xbin = linspace(pt1(1),pt2(1),abs(pt1(1)-pt2(1))+1);%判断两个点x位置的间隔
ybin = linspace(pt1(2),pt2(2),abs(pt1(2)-pt2(2))+1);%判断两个点y位置的间隔

if numel(xbin)>=numel(ybin) %哪个间隔大按哪个来
    ybin = linspace(pt1(2),pt2(2),abs(pt1(1)-pt2(1))+1);
    ybin = round(ybin);
else
    xbin = linspace(pt1(1),pt2(1),abs(pt1(2)-pt2(2))+1);
    xbin = round(xbin);
end
linkpts = zeros (numel(xbin),2);
linkpts (:,1) = xbin;%建立要连接点的像素位置
linkpts (:,2) = ybin;
for iLink = 1:numel(xbin)
    E(linkpts (iLink,1),linkpts (iLink,2)) = true;
end

end