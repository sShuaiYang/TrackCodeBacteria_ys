function [processedImages] = phaseContrastImageProcessing_ys(phaseImages)
%phase Contrast image processing code
% phase contrast 16bitͼ�� ת��Ϊ8bitͼ��
%ͼ��Ϊblackͼ�� ת��Ϊwhiteͼ��
% phase Contrastͼ�����������ͼ��ϸ����һЩ
%phaseContrastImageProcessingTest2ys ���в���Ч��
%Shuai Yang 2020.04.09
%phaseContrastImageProcessingTest3ys ���½��в���Ч�� 2020.04.17

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
    % edgefilter������һ��2��pixel���ı߿�
    I_border = ones(size(I))*255; I_border(3:size(I,1)-2,3:size(I,1)-2) = 0;
    I_border = uint8(I_border);
    I = imsubtract(I,I_border);
    I = imgaussfilt(I,2);
    [Iobrcbr] = cleanupImageUsingReconstuction_PhCImages(I);
    
    tempBW = edge((Iobrcbr),'log', 1/255, 1);%1/255,threshold;1,sigma.
    tempBW = bwareaopen( tempBW, 20);% ɾ��С���߶�
    % tempBW = bwmorph(tempBW,'bridge');% ��Ҫ ����ֺܶ�С�Ĵ���   
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
%��ϸ����mask����һЩ��̬ѧ����
% ��ʴ �Ǽ� ɾ��endpoints 
%���ùǼܽ����ع�
tempBW = imopen(tempBW,strel('disk',2));
tempBW = imopen(tempBW,ones(3));
tempBW = bwmorph(tempBW,'open');

tempBW = imerode(tempBW,ones(4));%��ʴ���зָ�
tempBW = bwareaopen(tempBW, 5);%ɾ��С�����ǵ�pixel
skel = bwmorph(tempBW,'skel',Inf);%��ȡ�Ǽ�
E = bwmorph(skel, 'endpoints');
skel = skel&(~E); %ɾ��endpoint Ҫ�����ͻ�����һ��
E = bwmorph(skel, 'endpoints');
skel = skel&(~E); %ɾ������
E = bwmorph(skel, 'endpoints');
skel = skel&(~E); %ɾ������
% E = bwmorph(skel, 'endpoints');
% skel = skel&(~E); %ɾ���Ĵ�
%����ϸ��
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
%��edge��ͼ����д��� ������<2pixeld��endpoints �������γɱպϵ�ͼ��
distThresh = 5;%����3��pixel�ĶԽ��߾���
E = bwmorph(BW, 'endpoints');
[x,y] = find(E);
C = [x,y];%endpoints ��xyλ��
tempLogic=C(:,1)>5&C(:,1)<2043&C(:,2)>5&C(:,2)<2043;%������Ұ�߿���Χ��ϸ��[5,2043]
C=C(tempLogic,:);
if size(C,1)<2 %���endpoint����ĿС��2���ټ���
    return
end
D=pdist(C);
Z = squareform(D);
pointIdx = true(1,size(C,1));
for iPoint = 1: size(C,1)
    if pointIdx(iPoint) == true
        [~, IX] = sort(Z(:,iPoint));
        k = IX(2);% �ҵ��ڶ�С�ľ���ֵ������,��һ���Լ�����Ϊ0
        if Z(k,iPoint) <= distThresh && pointIdx(k) == true
            pt1 = C(iPoint,:);pt2 = C(k,:);
            E = connectTwopointsInImage (E,pt1,pt2);
            pointIdx(k) = false;
            %һ����ֻ���Ӿ���������Ǹ��� ���Ӻ� �������㶼���ڽ���������һ���ж�
        end
        pointIdx(iPoint) = false; %�Ѿ�ѭ������points�����ظ��ж�
    end
    
end
BW = BW|E;
end
function [E] = connectTwopointsInImage (E,pt1,pt2)

xbin = linspace(pt1(1),pt2(1),abs(pt1(1)-pt2(1))+1);%�ж�������xλ�õļ��
ybin = linspace(pt1(2),pt2(2),abs(pt1(2)-pt2(2))+1);%�ж�������yλ�õļ��

if numel(xbin)>=numel(ybin) %�ĸ�������ĸ���
    ybin = linspace(pt1(2),pt2(2),abs(pt1(1)-pt2(1))+1);
    ybin = round(ybin);
else
    xbin = linspace(pt1(1),pt2(1),abs(pt1(2)-pt2(2))+1);
    xbin = round(xbin);
end
linkpts = zeros (numel(xbin),2);
linkpts (:,1) = xbin;%����Ҫ���ӵ������λ��
linkpts (:,2) = ybin;
for iLink = 1:numel(xbin)
    E(linkpts (iLink,1),linkpts (iLink,2)) = true;
end

end