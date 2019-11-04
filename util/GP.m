c = 0.1



K = @(x,y,ox,oy) exp(-pdist2(x-ox,y-oy).^2 ./ (2*c^2));

S = @(x,y) K(x,y,0.1,0.4)+K(x,y,0.2,0.7)+K(x,y,0.1,0.2);




xx = 0:0.01:1
X = meshgrid(xx,xx);
[i,j]=find(X~=0)
 
for ii = 1:10201

A(ii) = S(i(ii)/100,j(ii)/100);
end

%%
mA=reshape(A,[101,101]);
imagesc(mA)