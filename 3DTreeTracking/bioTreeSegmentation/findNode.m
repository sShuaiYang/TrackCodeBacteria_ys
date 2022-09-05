function [nodeSort,nodeList]=findNode(bioTree)
nodeList=[];
% nodeSort=zeros(1000,1000);
nodeSort=zeros(4000,4000);%20200109ys

for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            try
                bioTree{iframe}.node{iNode}.In;
                bioTree{iframe}.node{iNode}.Out;
            catch
                continue
            end                
            nodeSort(size(bioTree{iframe}.node{iNode}.In,2),size(bioTree{iframe}.node{iNode}.Out,2))=nodeSort(size(bioTree{iframe}.node{iNode}.In,2),size(bioTree{iframe}.node{iNode}.Out,2))+1;
            if ~(size(bioTree{iframe}.node{iNode}.In,2)==1 && size(bioTree{iframe}.node{iNode}.Out,2)==2)
%                 disp([iframe,iNode,size(bioTree{iframe}.node{iNode}.In,2),size(bioTree{iframe}.node{iNode}.Out,2)])
                nodeList=[nodeList;[iframe,iNode]];
            end
%             if size(bioTree{iframe}.node{iNode}.In,2)==3 && size(bioTree{iframe}.node{iNode}.Out,2)==3
%                 nodeList=[nodeList;[iframe,iNode]];
%             end
%             if ~((size(bioTree{iframe}.node{iNode}.In,2)==1 && size(bioTree{iframe}.node{iNode}.Out,2)==2) ...
%                     ||(size(bioTree{iframe}.node{iNode}.In,2)==1 && size(bioTree{iframe}.node{iNode}.Out,2)==1)) ...
%                     ||(size(bioTree{iframe}.node{iNode}.In,2)==2 && size(bioTree{iframe}.node{iNode}.Out,2)==2))
%                 nodeList=[nodeList;[iframe,iNode]];
%             end
        end
    end
end
end