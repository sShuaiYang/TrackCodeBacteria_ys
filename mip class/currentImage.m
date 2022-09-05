classdef currentImage
    
    properties
        imageSize
        mosaicWidth
        mosaicHeight
        imageOri
        imageBF
        maskBF
        maskFluo
        mosaicWhiteImage
        mosaicBlueMask
        mosaicGreenMask
        mosaicRedMask
        displayImage;
        focusImage;
        
    end
    
    methods
        function obj=initializeCurrentImage(obj)
                        obj.imageSize=[2048,2048];
%             obj.imageSize=[1392,1040];
            obj.mosaicWidth=800;
            obj.mosaicHeight=600;
            obj.imageOri=zeros(obj.imageSize(1),obj.imageSize(2),2,'uint16');
            obj.imageBF=zeros(obj.imageSize(1),obj.imageSize(2),2,'uint8');
            obj.maskBF=false(obj.imageSize(1),obj.imageSize(2));
            obj.maskFluo=false(obj.imageSize(1),obj.imageSize(2));
            obj.mosaicWhiteImage=true(obj.mosaicHeight,obj.mosaicWidth);
            obj.mosaicBlueMask=false(obj.mosaicHeight,obj.mosaicWidth);
            obj.mosaicGreenMask=false(obj.mosaicHeight,obj.mosaicWidth);
            obj.mosaicRedMask=false(obj.mosaicHeight,obj.mosaicWidth);
            
        end
    end
    
end

