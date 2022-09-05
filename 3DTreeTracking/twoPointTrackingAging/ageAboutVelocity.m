function ageAboutVelocity(bioTree)
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        bioTree=bioTreeVelocity(bioTree,[iframe,iRoot,0]);
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            bioTree=bioTreeVelocity(bioTree,[iframe,iRoot,0]);
        end
    end
end
end
function bioTree=bioTreeVelocity(bioTree,aimInfo)
frameFrequency=1/3;
if aimInfo(3)==0
    for i=1:size(bioTree{aimInfo(1)}.root{aimInfo(2)}.traceInfo.measurment,2)
        positionInfo(i,:)=bioTree{aimInfo(1)}.root{aimInfo(2)}.traceInfo.measurment{i}.Centroid;
    end
    positionInfo(:,1)=wden(positionInfo(:,1),'sqtwolog','s','sln',3,'sym10');
    positionInfo(:,2)=wden(positionInfo(:,2),'sqtwolog','s','sln',3,'sym10');
    bioTree{aimInfo(1)}.root{aimInfo(2)}.positionInfo=positionInfo;
    bioTree{aimInfo(1)}.root{aimInfo(2)}.velocityInfo=diff(positionInfo);
    bioTree{aimInfo(1)}.root{aimInfo(2)}.meanVelocity=mean(sum(bioTree{aimInfo(1)}.root{aimInfo(2)}.velocityInfo.^2,2).^(0.5));
    timeDelay=logspace(1,log10(2400),100);
    [msd,tao]=getMSD(positionInfo,timeDelay,frameFrequency);
    p=polyfit(log10(tao(~isnan(msd))),log10(msd(~isnan(msd))),1);
    bioTree{aimInfo(1)}.root{aimInfo(2)}.MSDslope=p(1);
end
if aimInfo(3)~=0
    for i=1:size(bioTree{aimInfo(1)}.node{aimInfo(2)}.Out{aimInfo(3)}.traceInfo.measurment,2)
        positionInfo(i,:)=bioTree{aimInfo(1)}.node{aimInfo(2)}.Out{aimInfo(3)}.traceInfo.measurment{i}.Centroid;
    end
    positionInfo(:,1)=wden(positionInfo(:,1),'sqtwolog','s','sln',3,'sym10');
    positionInfo(:,2)=wden(positionInfo(:,2),'sqtwolog','s','sln',3,'sym10');
    bioTree{aimInfo(1)}.node{aimInfo(2)}.Out{aimInfo(3)}.positionInfo=positionInfo;
    bioTree{aimInfo(1)}.node{aimInfo(2)}.Out{aimInfo(3)}.velocityInfo=diff(velocityInfo);
    bioTree{aimInfo(1)}.node{aimInfo(2)}.Out{aimInfo(3)}.meanVelocity=mean(sum(bioTree{aimInfo(1)}.root{aimInfo(2)}.velocityInfo.^2,2).^(0.5));
     timeDelay=logspace(1,log10(2400),100);
    [msd,tao]=getMSD(positionInfo,timeDelay,frameFrequency);
    p=polyfit(log10(tao(~isnan(msd))),log10(msd(~isnan(msd))),1);
    bioTree{aimInfo(1)}.node{aimInfo(2)}.Out{aimInfo(3)}.MSDslope=p(1);
end
end
function [msd,tao]=getMSD(position,timeDelay,frameFrequency)
tao=[1,2,3,4,5,6,7,8,9,fix(timeDelay(timeDelay>=10))];
msd=[];
for iTao=1:size(tao,2)
    pos_pre=position(1:end-tao(iTao),:);
    pos_next=position(1+tao(iTao):end,:);
    dataTemp=(pos_next-pos_pre).^2;
    msd=[msd,mean(dataTemp(:,1)+dataTemp(:,2))];
end
% change tao unit to second
tao = tao./frameFrequency;
end