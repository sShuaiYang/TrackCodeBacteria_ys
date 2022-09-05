function allTreeSave(bioTree,saveFile)
saveFile=strcat(saveFile,'\finalBioTree');
bytesTemp=whos('bioTree');
bioTree1=bioTree(1);
bytesTemp1=whos('bioTree1');
bytes=bytesTemp.bytes-bytesTemp1.bytes;
treeNum=fix(bytes/2000000000)+2;
savePoint=1;
savePoint=findSavePoint(bioTree,treeNum,savePoint);
startFrame=1;
for iTree=1:treeNum
    if treeNum==1
        save(saveFile,'bioTree');
        return;
    end
    bioTreeStack=bioTree(startFrame:savePoint(iTree));
    saveFile1=strcat(saveFile,num2str(iTree));
    if savePoint(iTree)==1
        save(saveFile1,'bioTreeStack','-v7.3');
    else
        save(saveFile1,'bioTreeStack');
    end
    startFrame=savePoint(iTree)+1;
    bioTreeStack=[];
end
end
function savePoint=findSavePoint(bioTree,treeNum,savePoint)
if treeNum==1
    return;
end
if ~isempty(savePoint)
    startFrame=2;
else
    startFrame=1;
end
if size(savePoint,2)==treeNum-1
    savePoint=[savePoint,size(bioTree,2)];
    return;
end
for iframe=100:100:size(bioTree,2)
    bioTreeTmep=bioTree(startFrame:iframe);
    bytesTemp=whos('bioTreeTmep');
    bytes=bytesTemp.bytes;
    pointNum=fix(bytes/2000000000);
    if pointNum==1
        savePoint=[savePoint,iframe-100];
        startFrame=iframe;
        if size(savePoint,2)==treeNum-1
            savePoint=[savePoint,size(bioTree,2)];
            return;
        end
    end
end
end