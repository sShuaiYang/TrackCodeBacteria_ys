function  cameraCorrectionForMultiFluoChannels(dirFile,fixedChannel,isCheck)
%���ڽ�����ͬӫ��ͨ��ͼ��֮���Ư�ƣ�ʹ��ͬӫ��ͨ����ϸ��ͼ���غ�
% isCheck = 0 or 1�������� default 0���Ƿ��ֶ�ȷ�Ͻ��
% 2020.04.02 Shuai Yang
% ��׼channel �����Զ���
% ����֮һ
% allChannels = {'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};

channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1�����ͨ�� �̲���
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2�����ͨ��

warning('on')
%%
fieldList =  dir([dirFile,filesep,'field*']);
dirField_check = [dirFile,filesep,fieldList(1).name];

% �ж��Ƿ�ֻ��cam 1 �����channel
chTF_cam1 = false (1,numel(channels_cam1));
for iCh_cam1 = 1:numel(channels_cam1)
    chName_cam1 = channels_cam1{iCh_cam1};
    if isfolder([dirField_check,filesep,chName_cam1])
        chTF_cam1(iCh_cam1) = true;
    end
end
if sum(chTF_cam1) == 0
    warning('All channels were captured by camera 2');
    warning('No fluo channel thus need camera correction');
    disp('All channels were captured by camera 2');
    disp('No fluo channel thus need camera correction');
    return
end

% �ж��Ƿ���ҪУ�� �Ƿ���cam2 �����channels
chTF_cam2 = false (1,numel(channels_cam2));
for iCh_cam2 = 1:numel(channels_cam2)
    chName_cam2 = channels_cam2{iCh_cam2};
    if isfolder([dirField_check,filesep,chName_cam2])
        chTF_cam2(iCh_cam2) = true;
    end
end


if sum(chTF_cam2) == 0
    warning('All channels were captured by camera 1');
    warning('No fluo channel thus need camera correction');
    disp('All channels were captured by camera 1');
    disp('No fluo channel thus need camera correction');
    return
end

% �ж�ͼ���Ƿ�У���� ����chan cam2 ��ͼ���Ƿ���СֵΪ0
for iCh_cam2 = 1:numel(channels_cam2)
    if chTF_cam2(iCh_cam2) == 1
        chName_cam2 = channels_cam2{iCh_cam2};
        dirImage = [dirField_check,filesep,chName_cam2];
        imageList = dir([dirImage,filesep,'*.tif']);
        image_Name = [dirImage,filesep,imageList(1).name];
        image_cam2 = import_tiff_stack(image_Name);
        if min(double(image_cam2(:))) == 0
            warning('Fluo Images captured by camera 2 have ALREADY been corrected.');
            disp('Fluo Images captured by camera 2 have ALREADY been corrected.');
            return
        else
            break
        end        
    end    
end
%%
if isfile(strcat(dirFile,'\transformInfo.mat'))
    load([dirFile,'\transformInfo.mat'])
    existedCameraTformCheck(dirFile,channels_cam1,channels_cam2,transformInfo);
else   
    [Im_cam1,Im_cam2] = autoGetTwoCamImages(channels_cam1,channels_cam2,dirFile);
    try
        [transformInfo] = autoTrackGetTwoCamImtransformInfo(Im_cam1,Im_cam2,dirFile);
    catch
        [transformInfo] = getTwoCamImtransformInfo(Im_cam1,Im_cam2,dirFile);
        %   msg = 'Correction failed. transformInfo does not exist';
        %   error(msg)
    end
end
%%
if nargin == 2 || isempty(isCheck)|| isCheck == 0
    % % % %�����ֶ� ֱ�ӽ���У��
    close all
    transformTwoCamImages_ys(transformInfo,dirFile);
    
elseif isCheck == 1
    % % ���ֶ�ȷ���Ƿ���ȷ
    
    prompt = 'Please check the correction result��is it OK? Y/N [Y]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
    end
    
    if upper(str) == 'Y'
        close all;
        transformTwoCamImages_ys(transformInfo,dirFile);
    else
        close all;
        msg = 'Error occurred. correction failed ';
        error(msg)
    end
    
end

end
%%
function existedCameraTformCheck(dirFile,channels_cam1,channels_cam2,transformInfo)
% check ���е�tform�Ƿ�����
% shuai yang 2020.05.28
checkImIdx = 1; %check ��frame���� ��һ��frame,default = 1;
fieldList =  dir([dirFile,filesep,'field*']);
dirField_check = [dirFile,filesep,fieldList(1).name];

for iCh_cam1 = 2:numel(channels_cam1)
    chName_cam1 = channels_cam1{iCh_cam1};
    if isfolder([dirField_check,filesep,chName_cam1])
        dirImage = [dirField_check,filesep,chName_cam1];
        imageList = dir([dirImage,filesep,'*.tif']);
        image_Name = [dirImage,filesep,imageList(checkImIdx).name];
        Im_cam1 = import_tiff_stack(image_Name);
        break
    end

end
% ����ͼimportӫ��ͼ�����������ӫ��ͨ��ͼ������ͼimport����ͼ��
if ~exist('Im_cam1','var')
    chName_cam1 = channels_cam1{1};
    dirImage = [dirField_check,filesep,chName_cam1];
    imageList = dir([dirImage,filesep,'*.tif']);
    image_Name = [dirImage,filesep,imageList(checkImIdx).name];    
    Im_cam1 = import_tiff_stack(image_Name);
end

for iCh_cam2 = 1:numel(channels_cam2)
    chName_cam2 = channels_cam2{iCh_cam2};
    if isfolder([dirField_check,filesep,chName_cam2])
        dirImage = [dirField_check,filesep,chName_cam2];
        imageList = dir([dirImage,filesep,'*.tif']);
        image_Name = [dirImage,filesep,imageList(checkImIdx).name];
        Im_cam2 = import_tiff_stack(image_Name);
        break
    end
end

if ~isa(transformInfo.(chName_cam2),'double')
    tform = transformInfo.(chName_cam2);
    Im_cam2 = flip(Im_cam2,1);
    Im_cam2 = imwarp(Im_cam2,tform,'OutputView',imref2d(size(Im_cam2)));
else
    tform = transformInfo.(chName_cam1);
    Im_cam1 = flip(Im_cam1,1);
    Im_cam1 = imwarp(Im_cam1,tform,'OutputView',imref2d(size(Im_cam1)));
end
if strcmp(chName_cam1,'BF1')
    I0 = Im_cam1;
    I0 = uint8(double(I0-min(I0(:)))/double(max(I0(:))-min(I0(:)))*255);
    I = imcomplement(I0);% equal to I = 255 - I0;
    I = imsubtract(I,I0);
    Im_cam1 = I;
else
    I = Im_cam1-100;
    Im_cam1 = uint8(double(I-min(I(:)))/double(max(I(:))-min(I(:)))*255);
end

I = Im_cam2-100;
Im_cam2 = uint8(double(I-min(I(:)))/double(max(I(:))-min(I(:)))*255);

% C = imfuse(image_cam1,image_cam2,'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);
C = cat(3,Im_cam2,Im_cam1,Im_cam1*0);
figure;
imshow(C,[])

end
%% 
function  [Im_cam1,Im_cam2] = autoGetTwoCamImages(channels_cam1,channels_cam2,dirFile)
% �Զ���ȡĳ��ʵ����������������ͼ��
% Shuai Yang 2022/6/24
checkImIdx = 1; %check ��frame���� ��һ��frame,default = 1;
fieldList =  dir([dirFile,filesep,'field*']);
dirField_check = [dirFile,filesep,fieldList(1).name];

for iCh_cam1 = 2:numel(channels_cam1)
    chName_cam1 = channels_cam1{iCh_cam1};
    if isfolder([dirField_check,filesep,chName_cam1])
        dirImage = [dirField_check,filesep,chName_cam1];
        imageList = dir([dirImage,filesep,'*.tif']);
        image_Name = [dirImage,filesep,imageList(checkImIdx).name];
        Im_cam1 = import_tiff_stack(image_Name);
        break
    end
end
% ����ͼimportӫ��ͼ�����������ӫ��ͨ��ͼ������ͼimport����ͼ��
if ~exist('Im_cam1','var')
    chName_cam1 = channels_cam1{1};
    dirImage = [dirField_check,filesep,chName_cam1];
    imageList = dir([dirImage,filesep,'*.tif']);
    image_Name = [dirImage,filesep,imageList(checkImIdx).name];    
    Im_cam1 = imread(image_Name);
end

for iCh_cam2 = 1:numel(channels_cam2)
    chName_cam2 = channels_cam2{iCh_cam2};
    if isfolder([dirField_check,filesep,chName_cam2])
        dirImage = [dirField_check,filesep,chName_cam2];
        imageList = dir([dirImage,filesep,'*.tif']);
        image_Name = [dirImage,filesep,imageList(checkImIdx).name];
        Im_cam2 = imread(image_Name);
        break
    end
end

if ~exist('Im_cam1','var') || ~exist('Im_cam2','var')
    Im_cam1 = [];
    Im_cam2 = [];
    disp('No image files or wrong folder ')
end

end