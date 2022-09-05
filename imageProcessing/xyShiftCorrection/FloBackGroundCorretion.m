function imageGFP=FloBackGroundCorretion(imageGFP,backGround)
backGround=double(backGround);
fluoBack=backGroundGet(imageGFP);
averageBack=mean(mean(backGround));
for iframe=1:size(imageGFP,3)
%     dispFrame(iframe);
    imageTemp=double(imageGFP(:,:,iframe));
%     averageGFP=mean(mean(imageTemp));
    imageGFP(:,:,iframe)=(imageTemp-fluoBack)./backGround*averageBack;    % change by jzy 10.14
end
% fprintf('\n');
imageGFP=uint16(imageGFP);
end
function dispFrame(iConnect)
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end
function  [gBack]=backGroundGet(gfpImage)
%首先要对红色荧光和绿色荧光图像去除背景
    gImage=gfpImage;
    gIndex=gImage(:);
    gIndex=sort(gIndex);
    gBack=mean(gIndex(1:round(end/6)));
end