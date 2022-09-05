function timeLapseMultiFieldBioInfoReshape(dirFile)
%load 每个视野的bioInfo，将其转化为时间点的bioInfo
%每一个荧光拍摄的时间点，一个bioInfo
% shuai Yang 2020.05.25
disp('Multifield bioInfo reshape to time-points bioInfo')
fieldList = dir([dirFile,filesep,'field*']);
fieldList = fieldListClean (fieldList);%只保留field的文件和两个系统文件
allBioInfo = cell(1, length(fieldList));
k = 1;
for iField= 1:(length(fieldList))
    if ~strcmp(fieldList(iField).name(1:5),'field')
        continue
    end
    disp(fieldList(iField).name);
    dirField = strcat(dirFile,'\',fieldList(iField).name);
    
    if isfile(strcat(dirField,'\','bioInfo.mat'))
        load([dirField,'\bioInfo.mat']);
        allBioInfo {k} = bioInfo;
        k = k+1;
    end
end

dirSave = strcat(dirFile,'\timeBioInfo');
mkdir(dirSave);

bioInfo = struct();
fields = fieldnames(allBioInfo{1});%结构体的fields
for iFds = 1 :numel(fields)
    bioInfo.(fields{iFds}) = [];
end

timePoints =  size(allBioInfo{1},2);
for t = 1:timePoints
    for iField = 1:numel(allBioInfo)
        bioInfo(iField)= allBioInfo{iField}(t);
    end
    fname = strcat(dirSave,'\','bioInfo_t',num2str(t-1,'%04d'),'.mat');
    save (fname,'bioInfo');
end

end
%%
function fieldList = fieldListClean (fieldList)
% 只保留field 文件和两个系统文件
templogic = false (1, numel (fieldList));
for iField = 1: numel (fieldList)
    
    if(isequal(fieldList(iField).name,'.')||... % 系统自带的两个隐文件夹
            isequal(fieldList(iField).name,'..')||...
            strcmp(fieldList(iField).name(1:5),'field')) 
        templogic(iField) = true;
    end
 
end
fieldList = fieldList(templogic);
end