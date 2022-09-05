function batchFluoColonyImageResultGetProtocal(dirFolder,fixedChannel)
%体式显微镜数据处理 多个文件夹批处理
%Shuai Yang 2021.03.29
fileList = dir(dirFolder);
for iFile = 1 : numel(fileList)
    
    if(isequal(fileList(iFile).name,'.')||... % 去除系统自带的两个隐文件夹
            isequal(fileList(iFile).name,'..')||...
            ~fileList(iFile).isdir)                % 去除遍历中不是文件夹的
        continue;
    end
    disp(fileList(iFile).name);
    dirFile = [dirFolder,filesep,fileList(iFile).name];
    fluoColonyImageResultGetProtocal(dirFile,fixedChannel);
    
end

end