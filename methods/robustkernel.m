function [yhat,lambda,gamma,Beta] = robustkernel(y,B,lambda,gamma,maxIter,Schangetol)
%% Smooth decomposition methods
%
%
%   Parameters:
%       -----------------
%       y      :    Signal to be smooth
%       B      :    Basis for Background
%       lambda :    Smoothing Parameter
%       gamma  :    Defect Parameter
%       isOrth :    Whether you have orthogonal basis for Defect
%       Initial:    Initial value for lambda and gamma
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

if ~iscolumn(y)
    y = y';
end

if nargin<7
    Schangetol = 1e-2;
    if nargin<6
        maxIter = 10;
    end
end

isAutoLambda = isempty(lambda);
isAutoGamma = isempty(gamma);

addpath(genpath('tensor_toolbox/'))


plus0 = @(x) (x>0).*x;
sizey = size(y);
ndim = length(sizey);


SChange = 10;
vec = @(x) x(:);
S = zeros(sizey);

iIter = 0;
t = 1;

H = B/(B + lambda * eye(size(B)));

%
while SChange > Schangetol && iIter < maxIter
    
    iIter = iIter + 1;
    Sold = S;
    yhat = H*(y-S);
    
    Ye = (y - yhat);
    if isAutoGamma
        gamma = 2*norminv(1-0.1)*1.4828*median(abs(Ye(:)-median(Ye(:))))*1;
    end
    S = wthresh(Ye,'s',gamma/2); 
    SChange = sum(sum(sum((S - Sold).^2)));

end
Beta = (B + lambda * eye(size(B)))\(y-S);
end






