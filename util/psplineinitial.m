function [ L,D,U,C,V,Z ] = psplineinitial( Y,B )
%     Initial the Bspline estimation with L,D,U,C,V,Z
ndim = sum(size(Y)~=1);
L = cell(ndim,1);
D = cell(ndim,1);
U = cell(ndim,1);
C = cell(ndim,1);
V = cell(ndim,1);
Z = cell(ndim,1);


for idim = 1:ndim
    L{idim} = sqrtm(B{idim}'*B{idim});
    L{idim} = L{idim} + 1e-8*eye(size(L{idim}));
    D{idim} = diff(eye(size(B{idim},2)),1);
    [U{idim},C{idim},V{idim}] =  svd((L{idim}')\(D{idim}'*D{idim})/(L{idim}));
    Z{idim} = B{idim}/(L{idim}')*U{idim};
end

end

