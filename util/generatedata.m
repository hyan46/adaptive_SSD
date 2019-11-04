function [ Y,Ytrue,Strue ] = generatedata(n,sigma,delta,datacase,R)
% Generate Data for simulation study



if datacase == 1
    % Multiple defect
    load('circSim2.mat')
    kk = exp(-((1:n)-n/2).^2/n^2*4);
    Ytrue = kk'*kk;
    BetaS(1:2,:)=0;
    BetaS(:,1:2)=0;
    BetaS(:,(end-1):end)=0;
    BetaS((end-1):end,:)=0;
    Strue = Bs{1} * BetaS *Bs{2}';
    Strue(Strue < 0.4) = 0;
    Strue = imresize(Strue,[n,n]);
    Strue(156:200,1:159) = 0;

elseif datacase == 2
    load('circSim2.mat')
    [XX,YY] = meshgrid(1:n,1:n);
    Strue = (XX-n/2).^2+(YY-n/2).^2<R^2;
elseif datacase == 3
    load('circSim2.mat')
    [XX,YY] = meshgrid(1:n,1:n);
    Strue = (XX-3*n/4).^2+(YY-n/4).^2<(R/2)^2;
    Strue((3*n/4-R/2):(3*n/4+R/2),(n/4-R/2):(n/4+R/2)) = 1;
end
    %Ytrue = B{1}*Beta*B{2}';
    Ytrue = imresize(Ytrue,[n,n]);
    Y = Ytrue + delta*Strue + normrnd(0,sigma,size(Strue,1),size(Strue,2));

end

