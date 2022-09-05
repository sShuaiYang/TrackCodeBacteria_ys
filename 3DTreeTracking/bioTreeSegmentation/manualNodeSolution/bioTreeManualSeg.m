function bioTree=bioTreeManualSeg(bioTree)
bioTree=autoTidyBioTree(bioTree);
nodeSortPre=[];
[nodeSort,~]=findNode(bioTree);
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
%     bioTree=type1NodeReduction(bioTree,1,25);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,2);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
% nodeSortPre=[];
% while ~isequal(nodeSortPre,nodeSort)
%     nodeSortPre=nodeSort;
%     bioTree=type_InNodeReduction(bioTree,1);
%     bioTree=type2NodeReduction(bioTree,1);
%     bioTree=mixMatchReducion(bioTree,2);
%     bioTree=type3NodeReduction(bioTree,15);
%     bioTree=type4NodeReduction(bioTree);
%     bioTree=autoTidyBioTree(bioTree);
%     [nodeSort,~]=findNode(bioTree);
% end
end