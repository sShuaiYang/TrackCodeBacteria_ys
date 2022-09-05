function batchconvertJob(dirOpen,dirSaveImage,backGround,picType)
% disp('select the original tiff image folder');
% dirOpen=uigetdir();%put the convert folder
% disp('select the converted tiff matrix folder');
% dirSaveImage=uigetdir();
addpath(dirOpen);
nameList=dir(dirOpen);

% add 10.3
% backGround=uint16(double(backGround)/65535*2048);

for i=1:(length(nameList)-2)
    disp(strcat('Strart Job',num2str(i)));
    disp(nameList(i+2).name);
    fileName=strcat(dirOpen,'\',nameList(i+2).name);
    if strcmp(fileName(end-2:end), 'tif')
        %     imageStack=convertifftoMatrixJob(fileName); % this is for 8bit image
        tic;imageStack=import_tiff_stack(fileName);toc;
        
        % add 10.3
%         xcrop=1:2160;
%         ycrop=1280:2560;
%         imageStack=imageStack(xcrop,ycrop,:);
%         imageStack=uint16(double(imageStack)/65535*2048);
        
        imageStack=backGroundCorrection(imageStack,backGround,'16bit');
        n_str=length(nameList(i+2).name);
        disp('save matrix images file')
        tic;
        if strcmp(picType,'b')
            imageStack=255-imageStack;
        end
        save(strcat(dirSaveImage,'\',nameList(i+2).name(1:n_str-4),'.mat'),'imageStack');
        clear imageStack;
        toc;
    end
end
rmpath(dirOpen);
end