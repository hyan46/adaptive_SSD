%% Generate Data
clc
clear all
addpath('util')
addpath('methods')
%%da
n = 200;
%delta = 0.3;
sigma = 0.05;
datacase = 1;
R = 30;
e = 0.4;

%[ Y,Ytrue,Strue ] = generatedata(n,sigma,delta,datacase,R);
kk = exp(-((1:n)-n/2).^2/n^2*4);
Ytrue = kk'*kk;

%
k = 13;
[ Bsg] = bsplineBasis(n,k,3);
mBs = rand(size(Bsg,2),size(Bsg,2))<0.03;
Strue = (Bsg*mBs*Bsg');
Strue(Strue<0.3) = 0;
imagesc(Strue)
delta = 0.5
Y = Ytrue + delta*Strue + normrnd(0,sigma,size(Strue,1),size(Strue,2));


%%
nx = size(Y,1);
ny = size(Y,2);
%
% Initialization
p = 0
F = kernelmatrix('gaussian',0.7);


B0 =  F((1:n)'/n,(1:n)'/n);
%%


gamma = 0.8;
maxIter = 100;
nrep = 1;
samplesize = 0;
W1 = Y*0+1;


eta = 40;
expeng = 1.5;
yhatold = 0;
yhat = 1.3;
mY = Y*0;
mYold = Y*0+2;
%
c = 2;



trueIntry = find(Strue~=0);
%

Aold = [];
ipoint = 1;

d = 0.2;
alltype = {'random','MED','grid'};

Fs = kernelmatrix('gaussian',d/10);
Bs0 = Fs((1:n)'/n,(1:n)'/n);

%%
nmethod = length(alltype);
nrep = 10;
npoint = 500;
allprecision = zeros(nmethod,npoint,nrep);
allrecall = zeros(nmethod,npoint,nrep);
comptime = zeros(nmethod,npoint,nrep);
samplesize = zeros(nmethod,npoint,nrep);
bkgmse = zeros(nmethod,npoint,nrep);
Anomaly250 = cell(nmethod,1);
Anomaly400 = cell(nmethod,1);
%%
for jrep = 1:405
    mBs = mBs*0;
    mBs(1,:) = 1;
    mBs(end,:) = 1;
    mBs(:,1) = 1;
    mBs(:,end) = 1;
    idxBs = find(mBs==0);
    mBs = mBs*0;
    mBs(idxBs(randperm(numel(idxBs),7))) = 1;
    Strue = (Bsg*mBs*Bsg');
    Strue(Strue<0.3) = 0;
    imagesc(Strue)
    delta = 0.5;
    Y = Ytrue + delta*Strue + normrnd(0,sigma,size(Strue,1),size(Strue,2));
    trueIntry = find(Strue~=0);
    %%
    
    for imethod = 2
        %%
        tic
        W = Y*0;
        W(end/2,end/2) = 1;
        %         if strcmp(type,'exp')||strcmp(type,'random')||strcmp(type,'MED')
        %             W(randperm(n^2,40)) = 1;
        %         end
        
        Aold = [];
        %
        type = alltype{imethod};
        if strcmp(type,'random')
            e = 0;
            Lambda = 0;
        elseif strcmp(type,'MED')
            e = 1e-9;
            Lambda = 1/10;
        elseif strcmp(type,'grid')
            e = 4;
            Lambda = 15;
        end
        parold = [];
        %
        
        for ipoint = 1:npoint
            %
            %d
            %
            
            Fs = kernelmatrix('gaussian',d/6);
            Bs_estimate0 = Fs((1:n)'/n,(1:n)'/n);
            %
            
            [xx,yy]=find(W~=0);
            [xxall,yyall]=find(W+1>0);
            x = [xx,yy]/200;
            xall = [xxall,yyall]/200;
            lW = logical(W);
            y = Y(lW);
            %
            B = F(x,x);
            Bs_estimate = Fs(x,x);         L = 2*(normest(Bs_estimate,1e-2)*1.1).^2;
            lambda = 0.001;
            nw = sum(sum(lW));
            
            %
            % Estimation
            %
            %mse(mYold - mY)
            
            if mse(mYold - mY)>1e-10 && sum(sum(lW))>30
                tic
                [ yhat,lambda,gamma,Beta] = robustkernel(y,B,lambda,[],maxIter);
                Ye = y - yhat;
                S = wthresh(Ye,'s',gamma/2);
                mBeta = W*0;
                mBeta(lW) = Beta;
                mYold = mY;
                mY = B0*mBeta*B0;
                shat = 1.4828*median(abs(Ye(:)-median(Ye(:))));
                
                gamma = 2*shat*sqrt(mean(sum(Bs_estimate.^2,2)))*norminv(1-0.05);
                
                [ S1,BetaS,lambda,gamma] = sparsekernel(y-yhat,Bs_estimate,lambda,gamma);
                
            else
                
                yhat = mY(lW);
                Ye = y - yhat;
                shat = 1.4828*median(abs(Ye(:)-median(Ye(:))));
                %gamma =  2*norminv(1-0.1)*shat;
                
                S = wthresh(Ye,'s',gamma/2);
                
                gamma = 2*shat*sqrt(mean(sum(Bs_estimate.^2,2)))*norminv(1-0.05);
                
                [ S1,BetaS,lambda,gamma] = sparsekernel(y-yhat,Bs_estimate,lambda,gamma);
                
            end
            
            %
            normalizedS = 2*(1-normcdf(S)).*(S>0);
            
            % estimated anomaly
            mBetaS = W*0;
            mBetaS(lW) = BetaS;
            Shat = ((Bs_estimate0*mBetaS*Bs_estimate0));
            levels = graythresh(Shat);
            Anomaly=im2bw(Shat,0.006);
            mnormalizedS = W*0;
            mnormalizedS(lW) = normalizedS;
            
            mS = W*0;
            mS(lW) = normalizedS;
            
            tic
            [miny,Aold,parold,d,dres] = samplingmethodstart2(Shat,Bs0,mnormalizedS,W,type,Lambda,e,Aold,parold);
            if ipoint == 250 || ipoint == 400
                toc
            end
            
            if isempty(dres)
                dres = 0;
            end
            miny(lW) = 0;
            miny(miny<0)  = 0;
            [ii,jj]= find(miny==max(max(miny)),1,'first');
            W(ii,jj) = max(max(W))+1;
            samplesize(imethod,ipoint,jrep) = sum(sum(lW))-1;
            bkgmse(imethod,ipoint,jrep) = mse(mY-Ytrue);
            dall(imethod,ipoint,jrep) = d;
            dresall(imethod,ipoint,jrep) = dres;
            %
            predIntry = find(Anomaly~=0);
            
            [ ~,~,allprecision(imethod,ipoint,jrep),allrecall(imethod,ipoint,jrep),~,~,~,~ ] = evaluateprec(predIntry,trueIntry,numel(mS));
            
            comptime(imethod,ipoint,jrep) = toc;
            
            
            if ipoint == 250
                Anomaly250{imethod} = Shat;
            elseif ipoint == 400
                Anomaly400{imethod} = Shat;
            end
%             
                        subplot(2,2,1)
                        imagesc(lW+Anomaly+Strue)
                        subplot(2,2,2)
                        imagesc(miny/mean(mean(miny)))
                        subplot(2,2,3)
                        plot(allrecall(imethod,:,jrep))
                        Fall = allprecision.*allrecall./(allprecision+allrecall);
                        %backgroundmse(i)= mse(mY - Ytrue);
                        subplot(2,2,4)
                        plot(Fall(imethod,:,jrep))

                        pause(0.01)
            %if j == 1
            %                 canomaly{i} = Anomaly;
            %                 cW{i} = W;
            %             end
            %title(sum(sum(lW)))
            
        end
        
    end
    
    
    
    %if mod(jrep,5) ==0 || jrep ==1
    %    save simu_res.mat allprecision allrecall samplesize comptime bkgmse
    %end
end
%
%save res.mat
%%
F = meanprecision.*meanrecall./(meanprecision+meanrecall);

plot(F')
legend(alltype{1},alltype{2},alltype{3})
%save onerun.mat
