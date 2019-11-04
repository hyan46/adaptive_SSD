function [ S,BetaS,Ye,yhat,lambda,gamma] = bsplineSmoothDecompauto(y,B,Bs,W,lambda,gamma,maxIter,Schangetol)
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
% isOrth = @(x) norm(x'*x - eye(size(x,2))) < 100*eps;
constructD =@(n) diff(eye(n),1);
% delta = 0.5;
sizey = size(y);
ndim = length(sizey);
% isBsOrthognal = all(cellfun(isOrthogonal,Bs));

if isempty(Bs)
    Lbs = 2;
    isOrth = 1;
    BetaS = y*0;
else
    Lbs = 2*norm(Bs{1})^2*norm(Bs{2})^2;
    isOrth = 2;
    X = zeros(size(Bs{1},2),size(Bs{2},2));
    BetaS = zeros(size(Bs{1},2),size(Bs{2},2));
end

if isempty(B)
    for i = 1:ndim
        B{i} = eye(size(y,i));
    end
end

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
    if isempty(W)
        for idim = 1:ndim
            if lambda(idim) >= 0
                H{idim} = Z{idim}*diag(1./(ones(size(C{idim},1),1) + lambda(idim)*diag(C{idim})))*Z{idim}';
            else
                H{idim} = 1/sizey(idim)*ones(sizey(idim));
            end
        end
        
        if ndim == 1
            yhat = H{1}*(y-S);
        elseif ndim == 2
            yhat = H{1}*(y-S)*H{2};
        elseif ndim >= 3
            yhat = double(ttm(tensor(y-S),H));
        end
    else
           yhat = bsplineweight( (y-S).*W,C,Z,W,lambda,400,Err); 
    end
    
    if isOrth == 2   % non-orthogonal basis
        if isempty(W)
            BetaSe = X + 2/Lbs*Bs{1}'*(y -Bs{1}*X*Bs{2}' - yhat)*Bs{2};
        else
            BetaSe = X + 2/Lbs*Bs{1}'*(W.*(y -Bs{1}*X*Bs{2}' - yhat))*Bs{2};
        end
        
        %         gamma = fminbnd(@lassogcv,1e-3,1);
        maxYe = max(abs(BetaSe(:)));
        if isAutoGamma && mod(iIter,5)==1
            gamma = graythresh(abs(BetaSe)/maxYe)*maxYe*Lbs;
        end
        
        BetaS = wthresh(BetaSe,'h',gamma/Lbs);
        % sign(BetaSe) .* plus0(abs(BetaSe)- gamma/Lbs);
        S = Bs{1} *BetaS* Bs{2}';
        t = (1+sqrt(1+4*told^2))/2;
        if iIter==1
            X = BetaS;
        else
            X = BetaS+(told-1)/t*(BetaS-BetaSold);
        end
        
    elseif isOrth == 1  % identity matrix
        Ye = (y - yhat); 
        
        maxYe = max(abs(Ye(:)));
        if isAutoGamma
            gamma = 2*graythresh(abs(Ye)/maxYe)*maxYe;
        end
%         S = wthresh(Ye,'s',gamma/2);
        S = (Ye>gamma/2).*Ye;
        %S = sign(Ye) .* plus0(abs(Ye)- gamma/2);
        BetaS = S;
        
    end
    
%      SChange = norm(tensor(double(BetaS~=0) - double(BetaSold~=0)))/(norm(tensor(double(BetaS~=0)))+1e-10)
     SChange = sum(sum(sum((S - Sold).^2)));
%     SChange = norm(tensor(double(S~=0) - double(Sold~=0)))/(norm(tensor(double(S~=0)))+1e-10)


       
end


% Beta = (B{1}'*B{1}+lambda*D{1}'*D{1})\B{1}'*(yhat)*B{2}/(B{2}'*B{2}+lambda*D{2}'*D{2});
Ye = y - yhat;
% BetaS = [];
end






