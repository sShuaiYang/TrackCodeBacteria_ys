function stctureFactor=intergrateMask(fftMaskImage)
stctureFactor=zeros(2,size(fftMaskImage,1));
Radius=0:1:size(fftMaskImage,1)-1;
for iX=1:size(fftMaskImage(:,1))
    for iY=1:size(fftMaskImage(:,2))
        iRadius=sqrt((iX-1)^2+(iY-1)^2);
        stctureFactor(1,Radius<=iRadius & Radius>iRadius-1)=stctureFactor(1,Radius<=iRadius & Radius>iRadius-1)+abs(fftMaskImage(iX,iY));
        stctureFactor(2,Radius<=iRadius & Radius>iRadius-1)=stctureFactor(2,Radius<=iRadius & Radius>iRadius-1)+1;
    end
end
% stctureFactor(2,:)=1./(stctureFactor(1,:)./stctureFactor(2,:));
stctureFactor(2,:)=1./(stctureFactor(1,:));
stctureFactor(1,:)=Radius.^2;
end