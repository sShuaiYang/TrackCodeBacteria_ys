function [rockFrame,allDataOri]=findRockFrameAndDelete(allDataOri)
rockFrame=[];
for iFrame=1:19999
    Vel=[];
    for jCell=1:size(allDataOri,2)
        %        cellFrame=allDataOri{jCell}.velocityData(:,1);
        %         if max(cellFrame==iFrame)==1
        %             Vel=[Vel,allDataOri{jCell}.velocityData(cellFrame==iFrame,14)];
        %         end
        if allDataOri{jCell}.velocityData(1,1)<=iFrame&&allDataOri{jCell}.velocityData(end,1)>=iFrame
            frameNum=find(allDataOri{jCell}.velocityData(:,1)==iFrame);
            Vel=[Vel,allDataOri{jCell}.velocityData(frameNum,14)];
        end
    end
    if min(Vel)>0.2&&size(Vel,2)>5
       rockFrame=[rockFrame,iFrame];
    end
end
allDataOri=rockToZero2(allDataOri,rockFrame);
end