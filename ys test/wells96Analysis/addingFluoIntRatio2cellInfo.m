function cellInfo_struct = addingFluoIntRatio2cellInfo(dirAll,cellInfo_struct,Ch2,Ch1,dataFolder)
% Shuai Yang 2021.09.09
% add intenstiy ratio to 96 wells cellInfo
% Exampel：Ch2 = 'sfGFP';Ch1 = 'mScarletI' ; default ch1 为 baseCh
% intRatio 结构是 单个细菌的平均值；方差；
for iSamp = 1:numel(cellInfo_struct)
    
    % case for empty samples
    if isempty(cellInfo_struct(iSamp).(['int',Ch2]) )
        continue
    end
    
    % load allDataCollect or dataCollect
    dirSample = [dirAll,filesep,cellInfo_struct(iSamp).sampIdx];
    dirLoad = [dirSample,filesep,dataFolder];%result_basic2/result_basic 
    
    if isfile(strcat(dirLoad,'\allDataCollect.mat'))
        load([dirLoad,filesep,'allDataCollect.mat']);
    elseif isfile(strcat(dirLoad,'\dataCollect.mat'))
        load([dirLoad,filesep,'dataCollect.mat']);
        allDataCollect{1} = dataCollect;
    else
         disp('allDataCollect or dataCollect does not exist')
    end 
    
    intRatio = NaN(numel(allDataCollect),2);
    for iTime = 1:numel(allDataCollect)
        try
            intCh1 = allDataCollect{iTime}.(['int',Ch2]);
            intCh2 = allDataCollect{iTime}.(['int',Ch1]);
            intRatio(iTime,1) = mean(intCh1./intCh2,'omitnan');
            intRatio(iTime,2) = std(intCh1./intCh2,'omitnan');
        catch
            msg = ['No cell was dectected in frame ',num2str(iTime),';',...
                ' Specify the ',Ch2,Ch1,'Ratio',' value as NaN.'];
            warning(msg);
            disp(msg);
            intRatio(iTime,1) = NaN;
            intRatio(iTime,2) = NaN;
        end
    end
    
    cellInfo_struct(iSamp).([Ch2,Ch1,'Ratio']) = intRatio;
end

end
