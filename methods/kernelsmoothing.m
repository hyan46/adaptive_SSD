function [ H] = kernelsmoothing( x,xtest,K,lambda)
    H = K(xtest,x)/(K(x,x) + lambda * eye(size(K(x,x))));
end

