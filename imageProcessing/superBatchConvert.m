function dirAll=superBatchConvert(picType)
clc;
dirAll=uigetdir();
% nameList=dir(dirAll);
cd(dirAll);
% for i=1:(length(nameList)-2)
%     disp(strcat('AllJob',num2str(i)));
%     disp(nameList(i+2).name);
%     dirOpen=strcat(dirAll,'\',nameList(i+2).name,'\','original data');
%     dirSave=strcat(dirAll,'\',nameList(i+2).name,'\','tiff2matlab');
%     cd(nameList(i+2).name);
dirOpen=strcat(dirAll,'\','original data');
dirSave=strcat(dirAll,'\','tiff2matlab');
tempFileName=dir();
backGround=import_tiff_stack(tempFileName(3).name);

% xcrop=1:2160;
% ycrop=1280:2560;
% backGround=backGround(xcrop,ycrop,:);

mkdir('tiff2matlab');
cd('original data');
batchconvertJob(dirOpen,dirSave,backGround,picType);
cd('..');
cd('..');
% end
end