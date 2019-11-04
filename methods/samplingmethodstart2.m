function [y,A,par,d,dres] = samplingmethodstart(Shat,Bs,S,W,type,Lambda,e,Aold,parold)
Anomaly=im2bw(Shat,0.006);
% Sampling methods

ini = 0;
ndim = sum(size(W)~=1);
n = size(W,1);
nW = sum(sum(W));
if ndim ==1
    xw = find(W~=0);
    xallw = find(W+1>0);
    x=xw/n;
    xall = xallw/n;
elseif ndim ==2
    [xx,yy]=find(W);
    [xxall,yyall]=find(W+1>0);
    [xxnew,yynew]=find(W==max(W(:)));
    x = [xx,yy]/n;
    xnew = [xxnew,yynew]/n;
    xall = [xxall,yyall]/n;
end

lW = logical(W);

if isempty(Aold)
    [xxnew,yynew]=find(W~=0);
    xnew = [xxnew,yynew]/n;
    A=min(pdist2(xall,xnew),[],2);
    if ndim==2
        A = reshape(A,[n,n]);
    end
else
    Anew = (pdist2(xall,xnew));
    if ndim==2
        Anew = reshape(Anew,[n,n]);
    end
    
    A = min(Anew,Aold);
end
AA = A;
if strcmp(type,'grid') ||strcmp(type,'random') 
    d = 1/15;
else
    d = min(max(AA(Anomaly)),0.1);
end
dres = min(max(AA(Anomaly)),0.1);

if isempty(d) || d==0
    d = 0.1;
end


if strcmp(type,'MED')
    Shat = Bs*S*Bs;
    B = real((e +Shat/sum(sum(W))).^(Lambda));
    if sum(sum(lW)) < ini
        y = A;
    else
        y = A.*B;
    end
    
    AA = A;
    AA(B<4*e) = 0;
    d = min(max(AA(:)),0.1);
    par =[];
elseif strcmp(type,'grid')
    xpoint = round(linspace(1,n,Lambda));
    xarr = zeros(1,n);
    xarr(xpoint) = 1;
    yy = xarr'*xarr;
    yy(yy~=0)=1:Lambda^2;
    y = yy;
    
    par = parold;
    if sum(lW(:))== Lambda^2
        Znew = gridsearch(n,Lambda,Shat>0.001,e);
        y = Znew;
        par = y;
    end
    if sum(lW(:)) > Lambda^2
        y = par;
    end
    B = [];
    
end

y(W==1) = 0;

end

