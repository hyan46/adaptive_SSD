function [ Yk,lambda ] = bsplineweight( Y,C,Z,W,lambda,maxIter,Err)
%

if isempty(Z)
   [ L,D,U,C,V,Z ] = psplineinitial( Y,C ); 
end
if nargin<5
    Err = 1e-5;
    if nargin<4
        maxIter = 50;
    end
end


sizey = size(Y);
ndim = length(sizey);
H = cell(ndim,1);
for idim = 1:ndim
    H{idim} = Z{idim}*diag(1./(ones(size(C{idim},1),1) + lambda(idim)*diag(C{idim})))*Z{idim}';
end


Yk = Y.*0;
Yknew = 0;

if isempty(W)
    if ndim == 1
        Yk = H{1}*(Y);
    elseif ndim == 2
        Yk = H{1}*(Y)*H{2};
    elseif ndim >= 3
        Yk = double(ttm(tensor(Y),H));
    end
    
    
else
    for iIter = 1:maxIter
        iIter ;
        Yknew = H{1}*((Y - Yk).*W + Yk)*H{2};
        diffnow = sum(sum(abs(Yk-Yknew)))/sum(sum(abs(Yk)));
        if diffnow<Err
            str = ['error bound achieved in Iteration' int2str(iIter)];
            %disp(str)
            break;
        end
        Yk = Yknew;
    end
end



end

