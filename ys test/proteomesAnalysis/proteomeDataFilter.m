function proteomeData =  proteomeDataFilter(proteomeData)
fliterTF = true(size(proteomeData,1),1);
% filter "potentialContaminant"
if isfield(proteomeData,'PotentialContaminant')
    for iProtein = 1:size(proteomeData,1)
        if strcmp(proteomeData(iProtein).PotentialContaminant,'+')
            fliterTF(iProtein) = false;
        end
    end
end
% filter "Reverse"
if isfield(proteomeData,'Reverse')
    for iProtein = 1:size(proteomeData,1)
        if strcmp(proteomeData(iProtein).Reverse,'+')
            fliterTF(iProtein) = false;
        end
    end
end
% filter "Reverse"
if isfield(proteomeData,'OnlyIdentifiedBySite')
    for iProtein = 1:size(proteomeData,1)
        if strcmp(proteomeData(iProtein).OnlyIdentifiedBySite,'+')
            fliterTF(iProtein) = false;
        end
    end
end
proteomeData = proteomeData(fliterTF);
end