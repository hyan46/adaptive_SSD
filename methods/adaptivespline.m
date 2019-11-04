function [ Yk,Yks,W,yhat ] = adaptivespline( Y, W, B, Bs ,neach)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


if nargin<6
    maxIter = 10;
end

Ys = Y.*W;
[S,BetaS,ye,yhat,l,g] = bsplineSmoothDecompauto(Ys,B,Bs,W,0.1,0.2,5,1e-3);
Sw = S.*W;
Yk = bsplineweight( Sw+ye,Bs,[],W,[0.01 0.01],2000,1e-3);
Yks = bsplineweight( Sw,Bs,[],W,[0.01 0.01],2000,1e-3);

e = abs(Yk);
Pe = e/sum(sum(e));
x = discretesample(Pe,neach);
W(x)=1;


end

