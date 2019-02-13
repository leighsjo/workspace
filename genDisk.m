function circ = genDisk(canvas,radius, x0,y0)
%This function is used to generate a filled circle.
%   Detailed explanation goes here
[imageSizeX, imageSizeY] = size(canvas);
[cIdx, rIdx] = meshgrid(1:imageSizeX, 1:imageSizeY);
circ = (rIdx - y0).^2 + (cIdx - x0).^2 <= radius.^2;

end

