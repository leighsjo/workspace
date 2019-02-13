function resizeFlag(imName,imSize)
%resizeFlag is used to resize flag if no additional pattern needs to be
%added.
%   imName      - image name
%   imSize      - desired image size
%
if nargin <2
    imSize = [225, 360];
end
path = '/Users/leixu/Documents/Administration/GroupMembers/Website/';
im = double(imread([path,imName]));
imNew = uint8(imresize(im, imSize, 'nearest'));
figure,imshow(imNew)
filename = [path,imName(1:end-4), '_old.png'];
imwrite(imNew,filename,'png');
end

