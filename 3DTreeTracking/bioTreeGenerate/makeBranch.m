function treeBranch=makeBranch(maskImage,firstFrame,endFrame) %This function can genrate CC array in the images Stacks
treeBranch=[];
connectMask=connectMatrix(maskImage,firstFrame,endFrame);
clear maskImage;
for iBranch=firstFrame:endFrame-1
    CCbranch=bwconncomp(connectMask(:,:,:,iBranch),26);
    treeBranch=[treeBranch,CCbranch];
end
end
function connectMask=connectMatrix(maskImage,firstFrame,endFrame)
connectMask=(false(size(maskImage,1),size(maskImage,2),2,endFrame-1));
for iframe=firstFrame:endFrame-1
    connectMask(:,:,1,iframe)=(maskImage(:,:,iframe));
    connectMask(:,:,2,iframe)=(maskImage(:,:,iframe+1));
%     connectMask(:,:,3,iframe)=im2uint8(false(size(maskImage(:,:,1))));
end
end