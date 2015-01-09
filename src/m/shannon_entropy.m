function retval=shannon_entropy(x,bw)
%SHANNON_ENTROPY computes Shannon's entropy.
%   SHANNON_ENTROPY(X,BW) computes the entropy of an input vector X using
%   bandwidth BW.

%   Copyright 2014 Han Gong, University of Bath

nbins = round((max(x)-min(x))/bw); % number of bins
P = hist(x(:),nbins,true); % Compute histogram
P(P==0) = []; % remove zeros-entries
P = P./numel(x);
retval = -sum(P.*log2(P));
