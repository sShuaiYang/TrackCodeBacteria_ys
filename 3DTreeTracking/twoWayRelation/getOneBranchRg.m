function [degreeRadius]=getOneBranchRg(bacteriaFrameInfo,ibranch)
tTime=0;
for ibranch=1:7  
    disp(ibranch)
    pointCentroid=[];
    for iTime=1:5:numel(bacteriaFrameInfo)
        branchIndex=bacteriaFrameInfo{iTime}.bacteriaInfo(:,5);
        orderNum=1:size(branchIndex,1);
        aimOrder=orderNum(branchIndex==ibranch);
        if ~isempty(aimOrder) && isempty(pointCentroid)
            pointCentroid=bacteriaFrameInfo{iTime}.centroidInfo(aimOrder,:);
        end
        if ~ isempty(pointCentroid) && ~isempty(aimOrder)
            tTime=tTime+1;
            degreeRadius(tTime,1)=iTime;
            degreeRadius(tTime,2)=numel(aimOrder);
            icentroid=bacteriaFrameInfo{iTime}.centroidInfo(aimOrder,:);
            properCentroid=pointCentroid;
            %         properCentroid=[mean(icentroid(:,1)),mean(icentroid(:,2))];
            degreeRadius(tTime,3)=(mean(pdist2(icentroid,properCentroid).^2)).^0.5;
        end
    end
end
end
% result=[];
% for i=1:5:numel(bacteriaFrameInfo)
%     meanR=mean(b(a==i));
%     if meanR~=0
%         result=[result;i,meanR];
%     end
% end