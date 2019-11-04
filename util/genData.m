function [ Y,Ytrue,Strue,Bt,Bst,Beta] = genData( iData,delta)
% 1:  theta data from RMS of the Guided Wave Test 
% 2:  the true data of theta of Loaded Circel 
% 3:  Simulated Loaded Circle theta with Scattered Defect
% 4:  Simulated loaded Circle with clustered Defect;
% 5:  Multicrystalline Bending Beam
% 6:  Monocrystalline Bending Beam
% 7:  Indentation data
% 8/9/10/11:  GTW, simulation datasets

if nargin<2
    delta = 0.3;
end
% if nargout<3
% Ytrue = [];
% Strue = [];
%     B = [];
Bst = [];
Beta = [];
% Strue = [];
Bt = [];
% end

if iData == 1 
    dataPATH = corrPATH('\Dropbox\MATLAB\DamagedComposites\AverageRMS.mat');
    load(dataPATH);
    Y = AverageRMS(1:270,11:100);
    Ytrue = [];
    
elseif iData == 2      % true loaded circle
%     dataPATH = corrPATH('\Dropbox\MATLAB\Polariscope samples\Simulation Data\disc');
%     load([dataPATH, '/theta.dat'])
%     thetatrue = theta;
%     [ thetawhole,~,~,~,~] = stressmap( dataPATH,0);
%     rec = [184.5100  119.5100  350.9800  343.9800];
    % [Y] = imcrop(thetawhole,rec);
%     thetatruecrop = imcrop(thetatrue,rec);
%     [Y,~,~,~,~] = stressmap( dataPATH,rec);
%     Ytrue = thetatruecrop/180*pi;
    load('circSim3.mat')
     
%     K = 300;
%     Strue(randsample(prod(size(Ytrue)),K)) = delta;
%     sigma = 0.1;

    Y = Ytrue + delta*Strue + normrnd(0,0.05,size(Strue,1),size(Strue,2));


elseif iData == 3     % Simulated loaded Circle with scattered Defect;
    load( 'circSim1.mat')
%     delta = 0.1;
    Ytrue = B{1}*Beta*B{2}';
%     S(S~=0) = 0.3;
    S(S~=0) = delta;
    Strue = S;
    Y = Ytrue + Strue + normrnd(0,0.05,size(S,1),size(S,2));
    Bt = B;
    Bst = Bs;
    
elseif iData == 4     % Simulated loaded Circle with clustered Defect;
    %
    load('circSim2.mat')
    
    Ytrue = B{1}*Beta*B{2}';
    BetaS(1:2,:)=0;
    BetaS(:,1:2)=0;
    BetaS(:,(end-1):end)=0;
    BetaS((end-1):end,:)=0;
    Strue = Bs{1} * BetaS *Bs{2}';
    Strue(Strue < 0.4) = 0;
    
    sigma = 0.05;
%     delta = 1;
    Y = Ytrue + delta*Strue + normrnd(0,sigma,size(S,1),size(S,2));
    Bt = B;
    Bst = Bs;
    
    
    
    %     %     S(S~=0) = 0.3;
    %     %     Strue = S;
    %     %     Y = Ytrue + Strue + normrnd(0,0.05,size(S,1),size(S,2));
    %     nx = size(Ytrue,1);
    %     ny = size(Ytrue,2);
    %
    %     snkx = 20;
    %     snky = 20;
    %     skx = round(nx/snkx);
    %     sky = round(ny/snky);
    %     sd = 3;
    %     ssd = 3;
    %     Bsx = bsplineBasis(nx,skx,ssd,3);
    %     Bsy = bsplineBasis(ny,sky,ssd,3);
    %     Bs = {Bsx,Bsy};
    %     BetaS = zeros(size(Bs{1},2),size(Bs{2},2));
    %     BetaS(randsample(prod(size(BetaS)),20)) = 1;
    %
    %     Bt = B;
    %     Bst = Bs;
    
elseif iData == 5    % load Multicrystalline Bending Beam
    %
    dataPATH = corrPATH('/Dropbox/MATLAB/Smoothing/Code/paper/Case Study/Data');
    load([dataPATH, '/4 pt multi.mat'])
    Y = Stat.shear;
    rY = imrotate(Y,-0.5);
    rec = [100    89.5    2909    199];
    Y= imcrop(rY,rec);
    imagesc(Y)
    %     title('MultiCrystalline Bending Beam','FontSize',20)
    set(gcf,'color','w')
    set(gca,'FontSize',12)
    colorbar
    Ytrue=0; Strue = 0;
elseif iData == 6    % load Multicrystalline mono Bending Beam
    %
    dataPATH = corrPATH('/Dropbox/MATLAB/Smoothing/Code/paper/Case Study/Data');
    load([dataPATH, '/4 pt mono.mat'])
    
    Y = Data.TauMax;
    %     rY = imrotate(Y,-0.5);
    %     rec = [100    89.5    2909    199];
    
    %     Y= imcrop(rY,rec);
    imagesc(Y)
    %     title('MultiCrystalline Bending Beam','FontSize',20)
    set(gcf,'color','w')
    set(gca,'FontSize',12)
    colorbar
    
elseif iData == 7    % load the indentation data
    %
    dataPATH = corrPATH('/Dropbox/MATLAB/Smoothing/Code/paper/Case Study/Data');
    load([dataPATH, '/indents_top_row.mat'])
    
    Yfull = Stat.shear;
    Y = imcrop(Yfull,[70,140,549,89]);
    imagesc(Y)
    %     title('MultiCrystalline Bending Beam','FontSize',20)
    set(gcf,'color','w')
    set(gca,'FontSize',12)
    colorbar
    Ytrue=0;
    Strue=0;
    
elseif ismember(iData,[8,9,10,11]) % GTW
    
    dataPATH = corrPATH('/Dropbox/MATLAB/Wavenumber_to_Effectivethickness');
    load([dataPATH, '/Wavenumber_Damage',int2str(iData - 7),'_NotSmooth.mat'])
    Y = IW; 
    imagesc(Y)
    set(gcf,'color','w')
    set(gca,'FontSize',12)
    colorbar
    
elseif iData == 12
    dataPATH = corrPATH('/Dropbox/MATLAB/Wavenumber_to_Effectivethickness');
    load([dataPATH, '/Wavenumber_Experiment_NotSmooth.mat'])
    Y = IW; 
    imagesc(Y)
    set(gcf,'color','w')
    set(gca,'FontSize',12)
    colorbar

    
end


end





