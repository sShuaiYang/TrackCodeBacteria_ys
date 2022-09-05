function outim = imshowColorlabel_ys(L,bwscreen)
% IMSHOWLABEL is used to display an integer image 
%给connected区域随机标记颜色
%Originla from Schnitzcells tacking codes
% L = bwlabel(bw1)
% Shuai Yang 2020.09.29
% L 为根据mask获得的标记connected区域
% bwscreen 为要重叠图像 作为mask的背景


if islogical(L)
%     L = double(L);%二值图像单一颜色
    L = bwlabel(L);
else
    L = logical(double(L));
    L = bwlabel(L);
end

% L2 has every non-background blob in range [2,256] and 
% sets background to one, corresp. to first entry in mymap
L2 = mod(L,255)+2;
L2(L == 0) = 1;

% M is the maximum color table entry, at most 256 colors
M = min(max(L(:))+2,256);
% create a color map
mymap = hsv(M);
% explicitly set the colormap's first entry to black for background
mymap(1,:)=[0 0 0];
% get sequence of random integers in range [1,maxcolors-1]
[~,I] = sort(rand(M-1,1));  
% randomly reorder mymap color entries [2,maxcolors]
mymap(2:end,:) = mymap(I+1,:);

if nargin>=2
	rgb = 0.5 * ind2rgb(L2,mymap);
	bwscreen = double(bwscreen);
	bwscreen = 0.5 * bwscreen / double(max(bwscreen(:)));
	rgb(:,:,1) = rgb(:,:,1) + bwscreen;
	rgb(:,:,2) = rgb(:,:,2) + bwscreen;
	rgb(:,:,3) = rgb(:,:,3) + bwscreen;
	figure;imshow(rgb);
	outim = rgb;
else
	figure; imshow(L2, mymap);	
    outim = L2;
end
end