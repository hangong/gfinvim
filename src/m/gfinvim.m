function [I1D,IL1] = gfinvim(I,varargin)
%GFINVIM computes Finlayson's 1D invariant image and its L1 chromaticity
%image. 
%   GFINVIM(I) returns a 1D invariant image of colour image I.
%
%   [I1D,IL1] = GFINVIM(I) returns both 1D invariant image I1D and its L1
%   chromaticity image IL1.
%
%   GFINVIM(... ,'Param1', value1, 'Param2', value2, ...) uses 
%   parameter-value pairs to control the algorithm.
%       Parameter name   Value
%       --------------   -----
%       'nd'             a positive integer specifies the number of angles
%                        to traverse (default: 180).
%       'entropy'        a string which is either 'shannon' or 'quadratic'
%                        that indicates the entropy computation method. If
%                        'shannon', Shannon's entropy (default) is used.
%                        Otherwise, a quadratic entropy (slow) is used.
%       'demo'           a BOOL value. If true, it shows the entropy plot
%                        and the corresponding 1D invariant image of each
%                        projection angle. Otherwise (default), it runs
%                        quietly.

%   Copyright 2014 Han Gong <gong@fedoraproject.org>, University of Bath

%   References:
%      Graham Finlayson et al. "Entropy Minimization for Shadow Removal".
%      IJCV, 2009.

para = inputparse(I,varargin{:}); % parse input parameters
I = max(im2double(para.I),eps); % to avoid Inf
RM = (I(:,:,1).*I(:,:,2).*I(:,:,3)).^(1/3); % geometric mean
sz = size(RM); nel = sz(1)*sz(2); % size and number of element of image
C = bsxfun(@rdivide,I,RM); % chromaticity image
rho = reshape(log(C),[],3)'; % log chromaticity on a plane

U = [1/sqrt(2), -1/sqrt(2), 0; 1/sqrt(6), 1/sqrt(6), -2/sqrt(6)]; % eigens
X = U*rho; % 2D points on a plane orthogonal to [1,1,1]

ang = linspace(0,180,para.nd); % angles to traverse 
E = zeros(para.nd,1); % entropy
if para.demo, figure; h = imshow(zeros(sz)); tt = title(''); end
for a = 1:para.nd % traverse each angle
    i = ang(a);
    e_t = [cos(i*pi/180);sin(i*pi/180)]; % projection vector
    mX = X'*e_t; % projection on e_t
    if strcmp(para.entropy,'shannon') % Shannon's entropy
        mX_std = std(mX);
        bw = 3.5*mX_std*nel^(-1/3); % bin width
        mXR = mean(mX)+3*[-mX_std,mX_std]; % range of data for entropy
        E(a) = shannon_entropy(mX(mX>mXR(1)&mX<mXR(2)),bw); 
    else % quadratic entropy using FGT (very slow) 
        X_var = var(X,0,2); % variance of X
        s = sqrt(X_var'*e_t); s = abs(1.06*s*(nel^(-1/5)));
        w = ones(1,nel); % weights
        sigma = 2*s;
        [xc,Ak] = fgt_model(mX',w,sigma,10,sqrt(nel),6);
        E(a) = sum(fgt_predict(mX',xc,Ak,sigma,6))/(4*pi*s*s);
    end
    if para.demo
        set(h,'CDATA',reshape(mat2gray(mX),sz)); drawnow;
        set(tt,'String',['Angle: ',num2str(i)]);
    end
end

if para.demo, plot(E); title('entropy plot'); end % plot entropy
[~,ma] = min(E); mi = ang(ma); % get minimum entropy
e_t = [cos(i*mi/180);sin(i*mi/180)]; % projection vector
e = [-sin(i*mi/180);cos(i*mi/180)]; % illumination vector
I1D = reshape(mat2gray(exp(X'*e_t)),sz); % 1D invariant image

if nargout>1
    p_th = e_t*e_t'; % project 2D points onto e_t
    X_th = p_th*X;   % X_th = e_t*(X*e_t')' <= slower but clear 
                     % to understand the projection
    mX = X'*e; mX_th = X_th'*e; % illumination brightness
    [~,bidx] = sort(mX,'descend'); bidx = bidx(1:ceil(0.01*nel)); % the %1
    X_E = (median(mX(bidx))-median(mX_th(bidx)))*e; % extra illumination
    X_th = bsxfun(@plus,X_th,X_E); % add back extra illumination
    rho_ti = U'*X_th; % convert to 3D log chromaticity
    c_ti = exp(rho_ti)'; % convert to 3D chromaticity
    r_ti = bsxfun(@rdivide,c_ti,sum(c_ti,2)); % L1 chromaticity
    IL1 = reshape(r_ti,[sz,3]); % L1 chromaticity image
end

function para = inputparse(I,varargin)
% parse input parameters
defaultnd = 180;
defaultentropy = 'shannon';
defaultdemo = false;
expectedentropy = {'shannon','quadratic'};

p = inputParser;
VER = version; % software version
if any(regexp(VER,'[a-z,A-Z]')) % it's MATLAB
    p.addRequired('I',@(x) isnumeric(x) && size(x,3)==3);
    p.addParamValue('nd',defaultnd,@(x) isnumeric(x) && x>1);
    p.addParamValue('entropy',defaultentropy,...
                 @(x) any(validatestring(x,expectedentropy)));
    p.addParamValue('demo',defaultdemo,@(x) isa(x,'logical'));
    p.parse(I,varargin{:});
else % otherwise, it should be GNU Octave
    p = p.addRequired('I',@(x) isnumeric(x) && size(x,3)==3);
    p = p.addParamValue('nd',defaultnd,@(x) isnumeric(x) && x>1);
    p = p.addParamValue('entropy',defaultentropy,...
                 @(x) any(validatestring(x,expectedentropy)));
    p = p.addParamValue('demo',defaultdemo,@(x) isa(x,'logical'));
    p = p.parse(I,varargin{:});
end
para = p.Results; % get parameters
