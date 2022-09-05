function  multiFiledPhCImCellSegTrack_supperSegger(dirFile) 
% supperSegger for phase contrast image segmentaion and tracking
%Shuai Yang 202.09.24
%Pa PhC image test suppersegger.mlx  test
supperSeggerImagesNamingConvention(dirFile);% supperSegger 格式转化

% Set constants
% Load the constants and set the desired values
CONST = loadConstants ('100XPa',1, 1);
% fit up to 5 foci in each cell;[5]
% Max number of foci to fit in each fluorescence channel (default = [0 0])
CONST.trackLoci.numSpots = [0,0];
% find the neighbors
CONST.trackOpti.NEIGHBOR_FLAG = false;
% not verbose state
CONST.parallel.verbose = 0;
% Setting clean flag to true to resgment data
clean_flag = 0;
%The segmentation is performed every skip files
skip = 1;
startEnd = [2 10];% 
showWarnings = 1;
close all; % close all figures
for iField = 1:numel(fieldList)
    dirField = [dirFile,filesep,fieldList(iField).name];
    disp (fieldList(iField).name);
    dirPhaseIm = [dirField,filesep,'BF1'];
    dirsupperSeg = [dirPhaseIm,filesep,'supperSegger'];
    BatchSuperSeggerOpti(dirsupperSeg,skip,clean_flag,CONST,startEnd,showWarnings);
    close all
    % clc
end

%  deleteSupperSeggerFiles(dirFile);

end
%% 
function  supperSeggerImagesNamingConvention(dirFile)
fieldList = dir([dirFile,filesep,'field*']);
for iField = 1:numel(fieldList)
    dirField = [dirFile,filesep,fieldList(iField).name];
    disp (fieldList(iField).name);
    dirPhaseIm = [dirField,filesep,'BF1'];
    phaseImList = dir([dirPhaseIm,filesep,'*.tif']);
    dirsupperSeg = [dirPhaseIm,filesep,'supperSegger'];
    if ~isfolder(dirsupperSeg)
        mkdir(dirsupperSeg)
    else
        [status] = rmdir(dirsupperSeg, 's');
        mkdir(dirsupperSeg);
        if ~status
            warning('Folder supperSeg already exists, cannot delete ')
        end
    end
    dirraw_im = [dirsupperSeg,filesep,'raw_im'];
    if ~isfolder(dirraw_im)
        mkdir(dirraw_im)
    end
    basename = 'BF1';
    elementsTime = 't';
    elementsXY = 'xy';
    currentXY = str2double(fieldList(iField).name(end-3:end));
    elementsPhase = 'c';
    c = 1;% phase image
    
    dirPhaseNew = [dirsupperSeg,filesep,elementsXY,num2str(currentXY),filesep,'phase'];
    if ~isfolder(dirPhaseNew)
        mkdir(dirPhaseNew)
    end
    
    parfor iImage = 1:numel(phaseImList)
        currentTime = iImage;
        newImageName = [basename,elementsTime,sprintf('%05d',currentTime),elementsXY,sprintf('%03d',currentXY),elementsPhase,num2str(c)];
        currentIm = imread([dirPhaseIm,filesep,phaseImList(iImage).name]);
        currentIm = rescale(double(currentIm))*255;
        currentIm = uint8(currentIm);
        imwrite(currentIm,[dirraw_im,filesep,newImageName,'.tif']);
        imwrite(currentIm,[dirPhaseNew,filesep,newImageName,'.tif']);
        %         copyfile([dirPhaseIm,filesep,phaseImList(iImage).name],[dirraw_im,filesep,newImageName,'.tif']);
        %         copyfile([dirPhaseIm,filesep,phaseImList(iImage).name],[dirPhaseNew,filesep,newImageName,'.tif']);
    end
    
end
end
%%
function deleteSupperSeggerFiles(dirFile)

fieldList = dir([dirFile,filesep,'field*']);
for iField = 1:numel(fieldList)
    dirField = [dirFile,filesep,fieldList(iField).name];
    disp (fieldList(iField).name);
    dirPhaseIm = [dirField,filesep,'BF1'];
    dirsupperSeg = [dirPhaseIm,filesep,'supperSegger'];
    [status] = rmdir(dirsupperSeg, 's');
    if ~status
        warning('Delete supperSegger files failed')
    end 
end
end