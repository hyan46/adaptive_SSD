function [ S,BetaS,Ye,yhat,Beta,lambda,gamma] = smoothDecompautokernel(y,B,Bs,lambda,gamma,maxIter,Schangetol)
%% Smooth decomposition methods
%
%
%   Parameters:
%       -----------------
%       y      :    Signal to be smooth
%       B      :    Basis for Background
%       Bs     :    Basis for Defect
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




plus0 = @(x) (x>0).*x;
sizey = size(y);
ndim = length(sizey);

if isempty(Bs)
    Lbs = 2;
    isOrth = 1;
    BetaS = y*0;
    Bs = eye(length(y));
    X = BetaS;
else
    Lbs = 2*(normest(Bs,1e-2)*1.1).^2;
    X = zeros(size(Bs,2),1);
    BetaS = zeros(size(Bs,2),1);
end

Ysize = size(y);


SChange = 10;
vec = @(x) x(:);
S = zeros(size(y));


iIter = 0;
t = 1;

H = B/(B + lambda * eye(size(B)));

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
    yhat = H*(y-S);
    
    BetaSe = X + 2/Lbs*Bs'*(y -Bs*X - yhat);

    maxYe = max(abs(BetaSe(:)));
    if isAutoGamma && mod(iIter,5)==1
        gamma = graythresh(abs(BetaSe)/maxYe)*maxYe*Lbs;
    end
    
    BetaS = wthresh(BetaSe,'h',gamma/Lbs);
    S = Bs * BetaS;
    t = (1+sqrt(1+4*told^2))/2;
    if iIter==1
        X = BetaS;
    else
        X = BetaS+(told-1)/t*(BetaS-BetaSold);
    end
    
    SChange = sum(sum(sum((S - Sold).^2)));
    
    
end

Ye = y - yhat;
Beta = (B + lambda * eye(size(B)))\(y-S);
end






