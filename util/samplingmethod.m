function [y,A,B] = samplingmethod(Bs,S,W,Wnew,k,e,type,Aold)

% Sampling methods

if nargin<6
    Aold = [];
end

    
ndim = sum(size(W)~=1);
n = size(W,1);
nW = sum(sum(W));
if ndim ==1
    xw = find(W);
    xallw = find(W+1);
    x=xw/n;
    xall = xallw/n;
    Shat = Bs*S;
    What = Bs*W;
elseif ndim ==2
    [xx,yy]=find(W);
    [xxall,yyall]=find(W+1);
    [xxnew,yynew]=find(Wnew);
    x = [xx,yy]/n;
    xnew = [xxnew,yynew]/n;
    xall = [xxall,yyall]/n;
end

if strcmp(type,'MEDsimple')
    Shat = Bs*S*Bs;
    What = Bs*W*Bs;
    a = (e+Shat'/sum(sum(W)));
    if isempty(Aold)
        A = pdist2(xall,x);AA = A;
        A = min(bsxfun(@times,A,a(W==1)'.^k),[],2);
        if ndim==2
            A = reshape(A,[n,n]);
        end
    else
        AAnew = (pdist2(xall,x));
        A = min(AAnew,Aold);AA = A;
        A = min(bsxfun(@times,A,a(W==1)'.^k),[],2);
    end
    %B = real(a'.^k);
    if ndim==2
        A = reshape(A,[n,n]);
    end
    %y = A.*B;
    y  =A;
    A = AA;
    B = A;
    
elseif strcmp(type,'MED2')
    Shat = Bs*S*Bs;
    What = Bs*W*Bs;
    if isempty(Aold)
        Ad=min(pdist2(xall,x),[],2);
        if ndim==2
            Ad = reshape(Ad,[n,n]);
        end
        A = (1./Ad).^k;
    else
        Anew = (pdist2(xall,xnew));
        if ndim==2
            Anew = reshape(Anew,[n,n]);
        end
        
        A = Aold + (1./Anew).^k;
    end
    
    B = real((e +Shat/sum(sum(W))));
    y = (1./A).^(1/k).*B.^f;
elseif strcmp(type,'simple')
    Shat = Bs*S*Bs;
    %What = Bs0*W*Bs0;
    if isempty(Aold)
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
    B = real((e +Shat/sum(sum(W))).^(k));
    y = A.*B;
elseif strcmp(type,'Exp')
    Shat = Bs*S*Bs;
    What = Bs*W*Bs;
    A= exp((-k)*What);
    B = (e+Shat);
    y= A.*B;
    
end

y(W==1) = 0;

end

