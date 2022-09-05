function [imageStack,backGroundPara]=backGroundCorrection(imageStack,backGround,imageType)
if strcmp(imageType,'8bit')
costantHighIn=1.2; %16 bit image use 1.25; 8bit image use 1.5
costanLowIn=0.7; %16 bit image use 0.85; 8bit image use 0.6
ampIntensity=2; %16 bit image use 2.8; 8bit image use 1
c=0.5; %16 bit image use 0.6; 8bit image use 0.3
end
if strcmp(imageType,'16bit')
    %really 11 bit
%     costantHighIn=1.6; %16 bit image use 1.25; 8bit image use 1.5
%     costanLowIn=0.85; %16 bit image use 0.85; 8bit image use 0.6
%     ampIntensity=4.5; %16 bit image use 2.8; 8bit image use 1
%     c=0.6; %16 bit image use 0.6; 8bit image use 0.3
%     
  costantHighIn=1.6; %16 bit image use 1.25; 8bit image use 1.5
    costanLowIn=0.85; %16 bit image use 0.85; 8bit image use 0.6
    ampIntensity=3.0; %16 bit image use 2.8; 8bit image use 1
    c=0.6; %16 bit image use 0.6; 8bit image use 0.3
    
    %really 16 bit
%     costantHighIn=1.6; %16 bit image use 1.25; 8bit image use 1.5
%     costanLowIn=0.85; %16 bit image use 0.85; 8bit image use 0.6
%     ampIntensity=4; %16 bit image use 2.8; 8bit image use 1
%     c=0.5;    
end
fprintf('\n');
backImage=double(backGround);
averageBack=mean(mean(backGround));
for iframe=1:size(imageStack,3)
    dispFrame(iframe);
    imageTemp=double(imageStack(:,:,iframe));
    averageI=mean(mean(imageTemp));
    ratio= averageI/averageBack;
    highIn=averageBack*costantHighIn;
    lowIn=averageBack*costanLowIn;
    slope=1/(highIn-lowIn);
    b=-lowIn/(highIn-lowIn);
    imageStack(:,:,iframe)=((((imageTemp./ratio)./backImage).*averageBack).*slope+b).*255; % backGround Normilzed
    imageStack(:,:,iframe)=(imageStack(:,:,iframe)-c*mean(mean( imageStack(:,:,iframe)))).*ampIntensity; % Intensity AMP
end
fprintf('\n');
imageStack=uint8(imageStack);
backGroundPara.costantHighIn=costantHighIn;
backGroundPara.costanLowIn=costanLowIn;
backGroundPara.ampIntensity=ampIntensity;
backGroundPara.c=c;
end
function dispFrame(iConnect)
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end
