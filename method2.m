clear all;
close all;
clc;
disp('')
disp('Prosze Czekac - trwa przetwarzanie obrazu...')
Image = im2double(imread('color1.jpg'));

Dilation = dilationFunction(Image);
Erosion = erosionFunction(Image);
ImgOpen = dilationFunction(Erosion);
ImgClose = erosionFunction(Dilation);
ContrastUpgrade = contrastFunction(Image,Dilation,Erosion);
MorphGrad = topHatGradFunction(Dilation-Erosion);
TopHatLight = topHatGradFunction(Image - ImgOpen);
TopHatDark = topHatGradFunction(ImgClose - Image);
DylRecon = dilatRecoFunction(Image, Dilation);
EroRecon = eroRecoFunction(Image, Erosion);
OpenRecon = dilatRecoFunction(Image, Erosion);
CloseRecon = eroRecoFunction(Image, ImgOpen);
GeoEro = geoEroFunction(Image, Erosion);
GeoDyl = geoDilaFunction(Image, Dilation);

disp('Zakonczono przetwarzanie. Wyswietlam obrazy...')

figure('Name','Przeksztalcenia morfologiczne'), set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,5,1)
imshow(Image);
title('Obraz oryginalny');
subplot(2,5,2)
imshow(Dilation);
title('Dylacja');
subplot(2,5,3)
imshow(Erosion);
title('Erozja');
subplot(2,5,4)
imshow(ImgOpen);
title('Otwarcie');
subplot(2,5,5)
imshow(ImgClose);
title('Zamkniecie');
subplot(2,5,6)
imshow(ContrastUpgrade);
title('Poprawa kontrastu');
subplot(2,5,7)
imshow(MorphGrad);
title('Gradient morfologiczny');
subplot(2,5,8)
imshow(TopHatLight);
title('Top-hat jasny');
subplot(2,5,9)
imshow(TopHatDark);
title('Top-hat ciemny');

figure('Name','Przeksztalcenia geodezyjne'), set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,4,1)
imshow(Image);
title('Obraz oryginalny');
subplot(2,4,2)
imshow(GeoEro);
title('Erozja Geodezyjna');
subplot(2,4,3)
imshow(GeoDyl);
title('Dylacja Geodezyjna');
subplot(2,4,4)
imshow(DylRecon);
title('Rekonstrukcja przez dylacje');
subplot(2,4,5)
imshow(EroRecon);
title('Rekonstrukcja przez erozje');
subplot(2,4,6)
imshow(OpenRecon);
title('Otwarcie przez rekonstrukcje');
subplot(2,4,7)
imshow(CloseRecon);
title('Zamkniecie przez rekonstrukcje');
disp('Wyswietlono obrazy')





function [result] = pixelValueComparatorFunction(p1, p2)

result = 0;
if p1(3) == p2(3) 
    if p1(2) == p2(2) 
        if p1(1) == p2(1) 
            result = 0;
        elseif p1(1) > p2(1)
            result = 1;
        else
            result = -1;
        end
    elseif p1(2) > p2(2)
        result = 1;
    else
        result = -1;
    end
elseif p1(3) > p2(3)
    result = 1;
else
    result = -1;
end

end






function [ PROCESSED_IMAGE ] = dilationFunction( image )
modelImage = rgb2hsv(image);

Comparator = zeros(9,3);
length1 = size(image,1);
length2 = size(image,2);

PROCESSED_IMAGE = modelImage;

for i=2:length1-2
    for j=2:length2-2
        element = 1;
        for temp=-1:1
            for l=-1:1
                Comparator(element,: ) = modelImage(i+temp,j+l,:);
                element=element+1;
            end
        end
        maximalElement = 0;
        z = 2;
        maximalPixelsVector = [];
        for x=1:9
            for l=z:9         
                euklides = sqrt((Comparator(x,1)-Comparator(l,1)).^2 + (Comparator(x,2)-Comparator(l,2)).^2 + (Comparator(x,3)-Comparator(l,3)).^2);
                
                if (euklides > maximalElement)
                    maximalElement = euklides;
                    maximalPixelsVector = [Comparator(x,:); Comparator(l,:)];
                elseif euklides == maximalElement
                    maximalPixelsVector = [maximalPixelsVector; Comparator(x,:); Comparator(l,:)];
                end
            end
            z = z + 1;
        end
        
        px = maximalPixelsVector(1,:);
        n = size(maximalPixelsVector,1);
        for p=1:n-1
            if pixelValueComparatorFunction(px, maximalPixelsVector(p+1,:)) < 0
                px = maximalPixelsVector(p+1,:);
            end
        end
        
        PROCESSED_IMAGE(i,j,:) = px;
    end
end
PROCESSED_IMAGE = hsv2rgb(PROCESSED_IMAGE);
end


function [ PROCESSED_IMAGE ] = erosionFunction( image )
modelImage = rgb2hsv(image);

Comparator = zeros(9,3);
length1 = size(image,1);
length2 = size(image,2);

PROCESSED_IMAGE = zeros(size(modelImage));

for i=2:length1-2
    for j=2:length2-2
        element = 1;
        for temp=-1:1
            for l=-1:1
                Comparator(element,: ) = modelImage(i+temp,j+l,:);
                element=element+1;
            end
        end
        
        maximalElement = 0;
        z = 2;
        maximalPixelsVector = [];
        for x=1:9
            for z=z:9
                euklides = sqrt((Comparator(x,1)-Comparator(z,1)).^2 + (Comparator(x,2)-Comparator(z,2)).^2 + (Comparator(x,3)-Comparator(z,3)).^2);
                
                if (euklides > maximalElement)
                    maximalElement = euklides;
                    maximalPixelsVector = [Comparator(x,:); Comparator(z,:)];
                elseif euklides == maximalElement
                    maximalPixelsVector = [maximalPixelsVector; Comparator(x,:); Comparator(z,:)];
                end
            end
            z = z + 1;
        end
        
        px = maximalPixelsVector(1,:);
        n = size(maximalPixelsVector,1);
        for p=1:n-1
            if pixelValueComparatorFunction(px, maximalPixelsVector(p+1,:)) > 0
                px = maximalPixelsVector(p+1,:);
            end
        end
        
        PROCESSED_IMAGE(i,j,:) = px;
    end
end
PROCESSED_IMAGE = hsv2rgb(PROCESSED_IMAGE);
end


function [PROCESSED_IMAGE] = topHatGradFunction(image)
PROCESSED_IMAGE = image;
for i = 1:size(image, 1)
    for j = 1:size(image, 2)
        temp = norm(image(i,j));
        PROCESSED_IMAGE(i, j, 1) = temp;
        PROCESSED_IMAGE(i, j, 2) = temp;
        PROCESSED_IMAGE(i, j, 3) = temp;
    end
end
end

function [ PROCESSED_IMAGE ] = eroRecoFunction( image, erosion )
PROCESSED_IMAGE = erosion;
for i=1:8
    PROCESSED_IMAGE = geoEroFunction(image,PROCESSED_IMAGE); 
    if PROCESSED_IMAGE == image
        break
   
    end
end
end

function [ PROCESSED_IMAGE ] = dilatRecoFunction( image, dilation )
PROCESSED_IMAGE = dilation;
for i=1:8
    PROCESSED_IMAGE = geoDilaFunction(image,PROCESSED_IMAGE); 
    if PROCESSED_IMAGE == image
        break
   
    end
end
end

function [PROCESSED_IMAGE] = contrastFunction(image, dilation, erosion)
PROCESSED_IMAGE = dilation;
for i=1:size(PROCESSED_IMAGE)
    for j=size(PROCESSED_IMAGE)      
        if (norm(image(i,j)-dilation(i,j)) <= norm(image(i,j)-erosion(i,j)))
            PROCESSED_IMAGE(i,j) = dilation(i,j);
        else
            PROCESSED_IMAGE(i,j) = erosion(i,j);
        end       
    end
end
end

function [ PROCESSED_IMAGE ] = geoEroFunction(image, erosion)
modelImage = rgb2hsv(image);
image_erosion = rgb2hsv(erosion);
PROCESSED_IMAGE = modelImage;
length1 = size(modelImage,1);
length2 = size(modelImage,2);
for i=1:length1
    for j=1:length2
        euclid1 = sqrt(image_erosion(i,j,1).^2 + image_erosion(i,j,2).^2 + image_erosion(i,j,3).^2);
        euclid2 = sqrt(modelImage(i,j,1).^2 + modelImage(i,j,2).^2 + modelImage(i,j,3).^2);
        if (euclid1 == euclid2)
            if pixelValueComparatorFunction(image_erosion(i,j,:), modelImage(i,j,:)) < 0
                PROCESSED_IMAGE(i,j,:) = image_erosion(i,j);
            end
        elseif euclid1 > euclid2
            PROCESSED_IMAGE(i,j,:) = image_erosion(i,j, :);
        end
    end
end

PROCESSED_IMAGE = hsv2rgb(PROCESSED_IMAGE);

end

function [ PROCESSED_IMAGE ] = geoDilaFunction(image, dilation)
modelImage = rgb2hsv(image);
image_dilation = rgb2hsv(dilation);
PROCESSED_IMAGE = modelImage;
length1 = size(modelImage,1);
length2 = size(modelImage,2);
for i=1:length1
    for j=1:length2
        euclid1 = sqrt(image_dilation(i,j,1).^2 + image_dilation(i,j,2).^2 + image_dilation(i,j,3).^2);
        euclid2 = sqrt(modelImage(i,j,1).^2 + modelImage(i,j,2).^2 + modelImage(i,j,3).^2);
        if (euclid1 == euclid2)
            if pixelValueComparatorFunction(image_dilation(i,j,:), modelImage(i,j,:)) > 0
                PROCESSED_IMAGE(i,j,:) = image_dilation(i,j);
            end
        elseif euclid1 < euclid2
            PROCESSED_IMAGE(i,j,:) = image_dilation(i,j, :);
        end
    end
end

PROCESSED_IMAGE = hsv2rgb(PROCESSED_IMAGE);

end
