function caculateFluoImage(bioTree,bacteriaFrameInfo,gfpImage,rfpImage)
stackInterval=100;

% gain gfpImage & rfpImage data
for iframe=1:stackInterval:numel(bacteriaFrameInfo)
    bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
    for iBac=1:size(bacteriaInfo,1)
    