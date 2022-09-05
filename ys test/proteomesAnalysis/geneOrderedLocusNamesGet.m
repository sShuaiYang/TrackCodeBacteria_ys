function proteomeData = geneOrderedLocusNamesGet(proteomeData)
% ys 2020.08.28
parfor iProtein = 1:size(proteomeData,1)
    FastaHeaders = proteomeData(iProtein).FastaHeaders;%从FastaHeaders中获得GeneSymbol
    k1 = strfind(FastaHeaders,'GN=');
    k2 = strfind(FastaHeaders,'PE=');
    if ~isempty(k1)&&~isempty(k2)
        geneName = FastaHeaders(k1+3:k2-2);
        locusName = geneOrderedLocusNamesGetFrGN(geneName);
        proteomeData(iProtein).LocusName = locusName;
    end
end

end


function locusName = geneOrderedLocusNamesGetFrGN(geneName)
filein = 'D:\proteomecomponent_Chromosome\uniprot-proteome_UP000002438+AND+(proteomecomponent_Chromosome).txt';
fileID = fopen(filein,'r');
nline = 0;
while ~feof(fileID) % 判断是否为文件末尾
    line_text = fgetl(fileID); % 从文件读行
    nline = nline+1;
    if strcmp(line_text(1:2),'GN') && contains(line_text,geneName)
        k = strfind(line_text,'OrderedLocusNames=');
        
        if ~isempty(k)
            locusName = line_text(k+18:k+24);
            break
        else%有可能OrderedLocusNames 在下面一行或两行
            lineStart = nline;
            while isempty(k)
                line_text = fgetl(fileID);
                k = strfind(line_text,'OrderedLocusNames=');
                nline = nline+1;
            end
            locusName = line_text(k+18:k+24);
            if (nline -lineStart) > 3
                warning(strcat('Ordered Locus Name of',geneName,' may be wrong'));
                disp(strcat('lineNum','[',numestr(lineStart),',',numestr(nline),']'));
            end
            break
            
        end
        
    end
end
fclose(fileID);
end
