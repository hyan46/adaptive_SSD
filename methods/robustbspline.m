function [ yhat,lambda,gamma] = robustbspline(y,B,W,lambda,gamma,maxIter,Schangetol)
%% Smooth decomposition methods
%
%
%   Parameters:
%       -----------------
%       y      :    Signal to be smooth
%       B      :    Basis for Background
%       lambda :    Smoothing Parameter
%       W      :    Weight Coeff
%       maxIter:    Maximum Iteration
%       -----------------
%  Output:
%       -----------------
%       S      :    Defect
%       BetaS  :    Coef for Defect
%       Ye     :    Noise part
%       yhat   :    Background
%       lambda :    Smoothing Parameter
%       gamma  :    Defect Parameter
%
%
%% Test & prepare the variables
%---

if nargin<7
    Schangetol = 1e-2;
    if nargin<6
        maxIter = 10;
    end
end


isAutoLambda = isempty(lambda);
isAutoGamma = isempty(gamma);

addpath(genpath(corrPATH('/Dropbox/MATLAB/toolbox/tensor_toolbox_2.5/')))

plus0 = @(x) (x>0).*x;
constructD =@(n) diff(eye(n),1);
sizey = size(y);
ndim = length(sizey);

Lbs = 2;
isOrth = 1;
BetaS = y*0;


if numel(lambda) == 1
    lambda = ones(ndim,1)*lambda;
end

Ysize = size(y);


SChange = 10;
H = cell(ndim,1);
vec = @(x) x(:);
S = zeros(size(y));
[ L,D,U,C,V,Z ] = psplineinitial( y,B );


iIter = 0;
t = 1;

%

while SChange > Schangetol && iIter < maxIter
    iIter = iIter + 1;
    Sold = S;
    BetaSold = BetaS;
    told = t;
    gcv = @(x) splinegcv(x,y,C,Z,0,[]);
    
    if isAutoLambda && iIter==1
        lambda = fminbnd(gcv,1e-2,1e3);
        lambda = ones(ndim,1)*lambda;
    end
    Err = 1e-4;
    % %
    yhat = bsplineweight( (y-S),C,Z,W,lambda,400,Err);
    Ye = (y - yhat);
    maxYe = max(abs(Ye(:)));
    if isAutoGamma
        gamma = 1.4828*median(abs(Ye(:)-median(Ye(:))));
    end
    S = wthresh(Ye,'s',gamma/2);
    %S = (Ye>gamma/2).*Ye;
    %S = sign(Ye) .* plus0(abs(Ye)- gamma/2);
    SChange = sum(sum(sum((S - Sold).^2)));
    
       
end


end






