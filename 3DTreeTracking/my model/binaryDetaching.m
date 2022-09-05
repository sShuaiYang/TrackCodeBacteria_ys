function bacPre=binaryDetaching(bNum,gNum)
bacPre=ones(bNum,1);
detachingRate=0.5;
for gTime=1:gNum
    bacNew=[];
    for i=1:numel(bacPre)
        bacNew=[bacNew;bacPre(i)+1;1];
    end
    bacNeed=false(numel(bacNew),1);
    for i=1:numel(bacNew)
        if rand<detachingRate
            bacNeed(i)=false;
        else
            bacNeed(i)=true;
        end
    end
    bacPre=bacNew(bacNeed);
end
end
