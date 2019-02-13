function [detectedFace,BB]= faceDetect(im,showImage)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin <2
    showImage = 'false';
end

faceDetector=vision.CascadeObjectDetector('ProfileFace'); 
im=rgb2gray(im); % convert to gray
BB=step(faceDetector,im); % Detect faces
iimg = insertObjectAnnotation(im, 'rectangle', BB, 'Face'); %Annotate detected faces.
if strcmpi(showImage, 'true')
    figure
    imshow(iimg); 
    title('Detected face');
    hold on,
    for i=1:size(BB,1)
        rectangle('position',BB(i,:),'Linewidth',2,'Linestyle','-','Edgecolor','w');
    end
    hold off
    N=size(BB,1);
    detectedFace = cell(N,1);
    for i=1:N
        detectedFace{i}=imcrop(im,BB(i,:));
    figure
    imshow(detectedFace{i}); 
    title('Outcropped face');
    end
end
end

