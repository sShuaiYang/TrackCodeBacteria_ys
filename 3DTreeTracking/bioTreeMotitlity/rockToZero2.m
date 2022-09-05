function allDataNew=rockToZero2(allDataOri,rockFrame)
cellNum=size(allDataOri,2);
  for i=1:size(rockFrame,2)
      for k=1:cellNum
          if  allDataOri{k}.velocityData(1,1)<=rockFrame(i)&&allDataOri{k}.velocityData(end,1)>=rockFrame(i)
              FrameNum=find(allDataOri{k}.velocityData(:,1)==rockFrame(i));              
              rockValue=allDataOri{k}.velocityData(FrameNum,[2,3,7,8,12,13]);
              allDataOri{k}.velocityData(FrameNum,[2,3,4,7,8,9,12,13,14])=0;
%               if  i==size(rockFrame,2)
                  tempDataSize=size(allDataOri{k}.denosieData,1)-FrameNum;
                  allDataOri{k}.denosieData(FrameNum+1:end,[2,3,4,5,6,7])=allDataOri{k}.denosieData(FrameNum+1:end,[2,3,4,5,6,7])-ones(tempDataSize,1)*rockValue;
%                   continue;
%               end
%               if allDataOri{k}.velocityData(1,1)<=rockFrame(i+1)&&allDataOri{k}.velocityData(end,1)>=rockFrame(i+1)
%                  NextFrameNum=allDataOri{k}.velocityData(:,1)==rockFrame(i+1); 
%                  tempDataSize=NextFrameNum-FrameNum;
%                  allDataOri{k}.denosieData(FrameNum+1:NextFrameNum,[2,3,4,5,6,7])=allDataOri{k}.denosieData(FrameNum+1:NextFrameNum,[2,3,4,5,6,7])-ones(tempDataSize,1)*rockValue;   
%               else
%                  tempDataSize=size(allDataOri{k}.denosieData,1)-FrameNum; 
%                  allDataOri{k}.denosieData(FrameNum+1:end,[2,3,4,5,6,7])=allDataOri{k}.denosieData(FrameNum+1:end,[2,3,4,5,6,7])-ones(tempDataSize,1)*rockValue;
%               end
          end
      end
  end
  allDataNew=allDataOri;   
end