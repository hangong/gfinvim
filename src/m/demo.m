% Demostration of Graham Finlayson's Invariant Image Derivation

%   References:
%      Graham Finlayson et al. "Entropy Minimization for Shadow Removal".
%      IJCV, 2009.

% Copyright 2014 Han Gong, University of Bath

I = imread('../data/4.tif'); % read an image
[I1D,IL1] = gfinvim(I,'entropy','shannon','demo',true);
figure;
subplot(1,3,1); imshow(I); title('original');
subplot(1,3,2); imshow(I1D); title('1D invariant');
subplot(1,3,3); imshow(IL1); title('L1 chromaticity');
