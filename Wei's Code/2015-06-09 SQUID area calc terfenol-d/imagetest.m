%Use this code to determine appropriate threshold value, validate by visual
%inspection yourself.

image1 = imread('sample.bmp'); % image at t = 0
image2 = imread('blank.bmp'); % testing image at 0 < t < t_f

threshold = 0.25; % value to be determined, use this value in main code

% MAIN CODE
difference = imabsdiff(image1,image2);

imwrite(difference,'difference.bmp','bmp')

filtered = im2bw(difference, threshold);
all = im2bw(image2, threshold);

imwrite(filtered,'filtered.bmp','bmp')
imwrite(all,'all.bmp','bmp')

% The number of white pixels is simply the sum of all the image pixel values since each white pixel has value 1.
    % If the white pixels have value 255 then divide the sum by 255.
    sampleArea = sum(filtered(:));
    % The number of black pixels is simply the total number of pixels in the image minus the number of white pixels.
    totalArea = sum(all(:));
    % Now calculate the ratio of white/black, and display in MATLAB window.
    areaRatio = sampleArea/totalArea;
    
wholeArea = (2.54*0.5)^2; % 0.5 in^2, expressed in cm^2
sampleArea = wholeArea*areaRatio
    