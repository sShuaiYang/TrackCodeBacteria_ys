function image=gpi
% getPixelIdxFromImageJ
MIJ.run('getSelectionCoordinates');
pause(0.1)
image=MIJ.getCurrentImage;
image=im2bw(uint8(image));
end