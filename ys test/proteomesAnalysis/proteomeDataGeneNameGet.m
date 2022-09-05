function [proteomeData,geneSymbols] = proteomeDataGeneNameGet(proteomeData)
% gene name get
geneSymbols = cell(1,size(proteomeData,1));
for iProtein = 1:size(proteomeData,1)
    FastaHeaders = proteomeData(iProtein).FastaHeaders;%从FastaHeaders中获得GeneSymbol
    k1 = strfind(FastaHeaders,'GN=');
    k2 = strfind(FastaHeaders,'PE=');
    if ~isempty(k1)&&~isempty(k2)
        proteomeData(iProtein).GeneSymbol = FastaHeaders(k1+3:k2-2);
        geneSymbols{iProtein} = proteomeData(iProtein).GeneSymbol;
    end
end
end
