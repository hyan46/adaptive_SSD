function [ K ] = kernelmatrix(type,c)
% 
% if ~iscolumn(x)
%         x = x';
% end
if strcmp(type,'gaussian') 
    K = @(x,y) exp(-pdist2(x,y).^2 ./ (2*c^2));
elseif strcmp(type,'expexp') 
    K = @(x,y) exp(exp(-pdist2(x,y).^2 ./ (2*c^2)));
elseif strcmp(type,'uniform')
    K = @(x,y) 1/2/c*double(pdist2(x,y)<c);
%      K = bsxfun(@rdivide, K, sum(abs(K)));
elseif strcmp(type,'polynomial')
    K = @(x,y) (x'*y+c(1))^c(2);
elseif strcmp(type,'Epanechnikov')
    K = @(x,y) (1-pdist2(x,y)./c).^2.*(pdist2(x,y)<c);
elseif strcmp(type,'cauchy')
    K = @(x,y) 1./(1+(pdist2(x,y)./c).^2 )/pi;
end

% Normalize


end

