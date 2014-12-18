
% HW-4: High Dynamic Range Imaging and Tone-mapping
% EECS 395
%
% Akshat Thirani
%
% 3rd November, 2014

% Part 2: Camera response curves 

% We need to model the brightness Z(i,j) = E(i) * B(j)

%           inverse of camera response curve: g(Z(i,j) = ln(E(i)) + ln(B(j))

% [g,lE] = gsolve(Z,B,l) where

%  Z(i,j) is the pixel values of pixel location number i in image j
%  B(j)   is the log delta t, or log shutter speed, for image j
%  l      is lamdba, the constant that determines the amount of smoothness



eB = [1/385 1/189 1/94 1/47];
B = zeros(4,1);

for i=1:4
    B(i) = log(eB(i));
    % B = log shutter speed
end

% Read the files

    folder = '/Users/Akshat/Desktop/EECS 395/HWs/HW-4/ImageStream/ImageSet2/';
    
    i1 = 'hw4_1.jpg';
    i2 = 'hw4_2.jpg';
    i3 = 'hw4_3.jpg';
    i4 = 'hw4_4.jpg';
    
    rows = 768;
    cols = 1024;    % Not sure

    ji1 = imread([folder i1]);
    ji2 = imread([folder i2]);
    ji3 = imread([folder i3]);
    ji4 = imread([folder i4]);
    
    I = cat(4,cat(4,cat(4,cat(4,ji1,ji2),ji3),ji4));
    
    
    for i = 1:size(I,4)            
        tR = (I(:,:,1,i));
        zr1(:,i) = reshape(tR(randsample(length(tR(:)),1000)), 1000,1);

        tG = (I(:,:,2,i));
        zg1(:,i) = reshape(tG(randsample(length(tG(:)),1000)),1000,1);

        tB = (I(:,:,3,i));
        zb1(:,i) = reshape(tB(randsample(length(tB(:)),1000)),1000,1);
    end
        
    g = zeros(1000,1);
    lE = zeros(1000,1);
   
    %%
    
    l = 5;

    [g_r, lE_r] = gsolve(zr1, B, l);
    [g_g, lE_g] = gsolve(zg1, B, l);
    [g_b, lE_b] = gsolve(zb1, B, l);
    
%         g_r = response curve for red
%         g_g = response curve for green
%         g_b = response curve for blue
% 
%         lE_r = log radiance for red
%         lE_g = log radiance for green
%         lE_b = log radiance for blue
    

x = 0:255;


for i=1:255
    for j=1:4
        Xr(i,j)=lE_r(i)+log(eB(j));
        Xg(i,j)=lE_g(i)+log(eB(j));
        Xb(i,j)=lE_b(i)+log(eB(j));
    end
end



plot(x, g_b);
title('Blue: g vs color values')
xlabel('[0,255]');
ylabel('g(Z)');
hold on;

plot(Xb,'.');
hold off;

%% _________________________________________________________________

% initialize output image to zeros, same size as one of the input images
 
outIndex_r = zeros(rows,cols);
outIndex_g = zeros(rows,cols);
outIndex_b = zeros(rows,cols);

for p = 1:256
    for k=1:4
        index = find(I(:, :, 1, k) == p);
        radiance = g_r(p) + log(eB(k));
        radiance = radiance/4;
        outIndex_r(index) = outIndex_r(index) + real(exp((radiance)));
        
        index = find(I(:, :, 2, k) == p);
        radiance = g_g(p) + log(eB(k));
        radiance = radiance/4;
        outIndex_g(index) = outIndex_g(index) + real(exp((radiance)));

        index = find(I(:, :, 3, k) == p);
        radiance = g_b(p) + log(eB(k));
        radiance = radiance/4;
        outIndex_b(index) = outIndex_b(index) + real(exp((radiance)));

    end
end

outIndex = cat(3, cat(3, outIndex_r, outIndex_g), outIndex_b);

% Dynamic Range:

max(exp(outIndex_r(:)))-min(exp(outIndex_r(:))) / mean2(exp(outIndex_r))

% _________________________________________________________________

% Normalizing the color channel 

maxE=0;
minE=1000;

for i=1:rows
    for j=1:cols
        temp1 = outIndex_r(i,j);
        temp2 = outIndex_g(i,j);
        temp3 = outIndex_b(i,j);
        
        if (max([temp1 temp2 temp3])>maxE)
            maxE=real(max([temp1 temp2 temp3]));
        end
        
        if (min([temp1 temp2 temp3])<minE)
            minE=real(min([temp1 temp2 temp3]));
        end
        
    end
end

% Result: maxE and minE

normE_r = zeros(rows,cols);
normE_g = zeros(rows,cols);
normE_b = zeros(rows,cols);

for i=1:rows
    for j=1:cols
        normE_r(i,j) = (outIndex_r(i,j)-minE)/(maxE-minE);
        normE_g(i,j) = (outIndex_g(i,j)-minE)/(maxE-minE);
        normE_b(i,j) = (outIndex_b(i,j)-minE)/(maxE-minE);
    end
end

% Apply a gamma curve to the image

gamma_r = zeros(rows,cols);
gamma_g = zeros(rows,cols);
gamma_b = zeros(rows,cols);

g=4;

for i=1:rows
    for j=1:cols
        gamma_r(i,j) = (normE_r(i,j))^g;
        gamma_g(i,j) = (normE_g(i,j))^g;
        gamma_b(i,j) = (normE_b(i,j))^g;
    end
end

G = cat(3,cat(3,cat(3,gamma_r,gamma_g),gamma_b));

% imagesc(gamma_r)

% I liked g = 3

% Global tone mapping operator from Reinhard '02

normE = cat(3,cat(3,cat(3,normE_r,normE_g),normE_b));

l = rgb2gray(normE);

% Log average exposure

lAvg=0;

for i=1:rows
    for j=1:cols
        lAvg = lAvg + log(l(i,j));
    end
end

lAvg = lAvg/(rows*cols);

lAvg = exp(lAvg);

% Scale the image

a = .15;

t = zeros(rows,cols);

for i=1:rows
    for j=1:cols
        t(i,j) = (a/lAvg)*(l(i,j));
    end
end

tmax = max(max(t));

% Reinhard tone-mapping operator

lTone = zeros(rows,cols);

for i=1:rows
    for j=1:cols
        lTone(i,j) = (t(i,j)*(1+(t(i,j)/(tmax^2))))/(t(i,j)+1);
    end
end

% Scaling operator

M = zeros(rows,cols);

for i=1:rows
    for j=1:cols
        M(i,j) = lTone(i,j)/l(i,j);
    end
end

Rnew = zeros(rows,cols);
Gnew = zeros(rows,cols);
Bnew = zeros(rows,cols);

% for i=1:rows
%     for j=1:cols
%         Rnew(i,j) = M(i,j) * normE(i,j,1);
%         Gnew(i,j) = M(i,j) * normE(i,j,2);
%         Bnew(i,j) = M(i,j) * normE(i,j,3);
%     end
% end
%
for i=1:rows
    for j=1:cols
        Rnew(i,j) = M(i,j) * G(i,j,1);
        Gnew(i,j) = M(i,j) * G(i,j,2);
        Bnew(i,j) = M(i,j) * G(i,j,3);
    end
end

RGBImage = cat(3, cat(3, Rnew, Gnew), Bnew);

imshow(RGBImage)




