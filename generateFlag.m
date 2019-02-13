function generateFlag(imName)
%generateFlag is used to generate flag with additional pattern and resize
%the image to predetermined size.
%   imName      - image name

path = '/Users/leixu/Documents/Administration/GroupMembers/Website/';
im = double(imread([path,imName]));
[m,n,~] = size(im);
filter = 1-0.15*BrownianField(0.3,[m,n]);
filter(m+1:size(filter,1),:) = [];
filter(:,n+1:size(filter,2)) = [];
filter = (filter - min(filter(:)))/(max(filter(:))-min(filter(:)));
hsv = rgb2hsv(im);
hsv(:,:,3) = hsv(:,:,3).*filter;
imNew = uint8(imresize(hsv2rgb(hsv), [225 360], 'nearest'));
figure,imshow(imNew)
filename = [path,imName(1:end-4), '_old.png'];
imwrite(imNew,filename,'png');
end

