function main()
%MAIN minimal implementation of GLRT
%   now for the implementation to give a better feeling I will first build
%   2 example images, one containing only noise and another containing a
%   particle in the middle. For the particle I just use a gaussian
%   approximation to the PSF.

% lets agree on a ROI size of 13 x 13
ROIsize = [13,13];

% ROI image containing noise
imRand = uint16(ones(ROIsize).*20);
imRand = imnoise(imRand,'poisson');

% ROI image containing a gaussian spot on its center with some noise. For
% this I first build a gaussian with appropriate dimentions - G. For a
% water immersion obj the best-fit gasuss to the PSF has: sigma = 0.25
% wavelength / NA
% emission wavelength
emWave   = 600; % in nm
% objective NA
objNA    = 1.2;
% sigma of gaussian in nm
sigma_nm = 0.25 * emWave/objNA;
% typical pixel size for SR microscopy
pxSizeNm = 100;
% sigma of gaussian in pixels
sigma_px = sigma_nm / pxSizeNm;
% get input for my 2D gaussian function generator
xid = 1:ROIsize(1);
yid = 1:ROIsize(2);
% gaussian PSF will be in the middle of the ROI
pos = [mean(xid),mean(yid)];
% gaussian PSF, see end of file for the function definition
[G] = gaus2D(pos, [sigma_px, sigma_px], xid,yid,1);
% create image and add noise
imParticle = G.*150+20;
imParticle = uint16(imParticle);
imParticle = imnoise(imParticle,'poisson');
%%
testROI = cat(3,imRand,imParticle);
LRT = zeros(1,2);
for i = 1:2
    im = testROI(:,:,i);
    im = double(im);
    assert(size(im,1)==size(im,2),'ROI must be squared')
    % get size of ROI
    ws = size(im,1);
    N   = ws^2;
    % calculation of G_bar and sGbar2.
    % get model gaussian, a bit redundant but good to have it here also
    % get input for my 2D gaussian function generator
    xid = 1:ws;
    yid = 1:ws;
    % gaussian PSF will be in the middle of the ROI
    pos = [mean(xid),mean(yid)];
    % gaussian PSF -  note that sigma_px is a input given to the GLRT
    [G] = gaus2D(pos, [sigma_px, sigma_px], xid,yid,1);
    %       Mean of the gaussian PSF
    G_mean  = mean(G(:));
    %       G bar: gaussian PSF minus its mean
    G_bar = G - G_mean;
    %       sum of G bar squared
    sGbar2 = sum(sum(G_bar.^2));
    %       convolution of the image with G_bar
    imConv = sum(sum(G_bar.*im));
    %       calculation of I_hat
    I_hat  = imConv/sGbar2;
    %       denominator, which depends on how good it fits
    %       to random noise
    den = var(im(:))*N;
    %   calculation of the difference betwen the likelihood
    %   of having the data explained by random noise or by
    %   a gaussian. L(H_0)-L(H_1)
    LRT(i) = ( (N)/2 ) * (log(1-((I_hat^2 * sGbar2)/den)));
                   
end

figure(1)
subplot(1,2,1)
imagesc(imRand)
axis image
title({'no particle example',['GLRT: ' num2str(LRT(1),3)]})

subplot(1,2,2)
imagesc(imParticle)
axis image
title({'particle example',['GLRT: ' num2str(LRT(2),3)]})

end

function [G] = gaus2D(pos, sig, xid,yid,maxCount)
    %Simple 2D gaussian, just to make the code easier to read
    switch nargin
        case 4
            A = 100;
        case 5
            A = maxCount;
        otherwise
            error('Not enough input');
    end

    % due to the way meshgrid works y comes first and x second. this because I
    % want x to be the first dimention of my array and y the second
    [y,x] = meshgrid(xid,yid);
    sigX = sig(1);
    sigY = sig(2);
    x0 = pos(1);
    y0 = pos(2);

    xPart = ((x-x0).^2) ./ (2*sigX^2);
    yPart = ((y-y0).^2) ./ (2*sigY^2);


    G = A.*exp( -(xPart + yPart));

end