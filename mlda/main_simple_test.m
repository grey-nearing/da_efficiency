%% --- Init Workspace -----------------------------------------------------
clear all; close all; clc;

restoredefaultpath
addpath('tools');

%% --- Experimental Controls ----------------------------------------------

% dimensions
nSites = 129; % number of sites
nMin = 1000; % minimum number of points per site

% regression parameters
kFold = 10; % number of k-fold validation splits
verbose = 0; % regression screen report

% regression types
modelTypes = [{'linearRegression'},{'neuralNetwork'}];%,{'randomForest'}];%,{'gaussianProcess'}];
modelTypesShort = [{'MLR'},{'ANN'}];%,{'TBG'}];%,{'GPR'}];
nModels = length(modelTypes);

% mutual information bins
infoBinEdges = -1:0.01:2;

% load new data
newData = 0;

%% --- Load Data ----------------------------------------------------------

% if loading new data
if newData
    
    % find number of time steps
    temp = load('open_loop/open_1.out');
    nTimes = length(temp);
    
    % init storage
    lprm = zeros(nTimes,nSites)./0;
    open = zeros(nTimes,nSites)./0;
    scan = zeros(nTimes,nSites)./0;
    
    % screen report
    fprintf('Loading data ... \n');
    fprintf(strcat(repmat('.',[1,nSites]),'\n')); tic;
    
    for s = 1:nSites
        
        % screen report
        fprintf('.');
        
        try
            % load open loop
            fname = strcat('open_loop/open_',num2str(s),'.out');
            readData = load(fname);
            assert(length(readData) == nTimes);
            open(:,s) = readData(:,1);
            
            % load lprm
            fname = strcat('lprm_data/lprm_',num2str(s),'.txt');
            readData = load(fname);
            assert(length(readData) == nTimes);
            lprm(:,s) = readData(:,5);
            
            % load open loop
            fname = strcat('scan_data/scan_',num2str(s),'.txt');
            readData = load(fname);
            assert(length(readData) == nTimes);
            scan(:,s) = readData(:,5);
            
        catch
            open(:,s) = 0/0;
            lprm(:,s) = 0/0;
            scan(:,s) = 0/0;
        end
        
    end % s-loop
    fprintf('. finished; time = %f[s] \n',toc);
    
    % convert missing values
    open(open<-9990) = 0/0;
    lprm(lprm<-9990) = 0/0;
    scan(scan<-9990) = 0/0;
    
    % save final data file
    fname = 'training_data.mat';
    save(fname,'open','lprm','scan');
    
else % if not loading new data
    fname = 'training_data.mat';
    load(fname)
end
    
% sample an equal number of points from each site
enoughData = zeros(nSites,1);
allData = zeros(nMin,nSites,3)./0;

for s = 1:nSites
    
    siteData = [open(:,s),lprm(:,s),scan(:,s)];
    iMissing = find(any(isnan(siteData')));
    siteData(iMissing,:) = [];
    nSiteData = length(siteData);
    
    if nSiteData >= nMin
        enoughData(s) = 1;
        iUse = randperm(nSiteData,nMin);
        allData(:,s,1) = siteData(iUse,1);
        allData(:,s,2) = siteData(iUse,2);
        allData(:,s,3) = siteData(iUse,3);
    end
    
end % s-loop

% remove sites without enough data
allData(:,enoughData==0,:) = [];

% dimensions
nSites = size(allData,2);

%% --- Calclate Global and Site-Specific Statistics -----------------------

% global stats
o = allData(:,:,1); o = o(:);
l = allData(:,:,2); l = l(:);
s = allData(:,:,3); s = s(:);

nseGlobal(1) = nse(o,s);
rhoGlobal(1) = rho(o,s);
nmiGlobal(1) = nmi(o,s,infoBinEdges,infoBinEdges);

nseGlobal(2) = nse(l,s);
rhoGlobal(2) = rho(l,s);
nmiGlobal(2) = nmi(l,s,infoBinEdges,infoBinEdges);

% site-specific
for s = 1:nSites
       
    % open loop
    nseSites(s,1) = nse(squeeze(allData(:,s,1)),squeeze(allData(:,s,end)));
    rhoSites(s,1) = rho(squeeze(allData(:,s,1)),squeeze(allData(:,s,end)));
    nmiSites(s,1) = nmi(squeeze(allData(:,s,1)),squeeze(allData(:,s,end)),infoBinEdges,infoBinEdges);
    
    % remote sensing obs
    nseSites(s,2) = nse(squeeze(allData(:,s,2)),squeeze(allData(:,s,end)));
    rhoSites(s,2) = rho(squeeze(allData(:,s,2)),squeeze(allData(:,s,end)));
    nmiSites(s,2) = nmi(squeeze(allData(:,s,2)),squeeze(allData(:,s,end)),infoBinEdges,infoBinEdges);
     
end % s-loop

%% --- Train Global Models ------------------------------------------------

% screen report
fprintf('--- Training Global Models -------------------------- \n');

% init storage - assimilation
assimGlobalPreds = zeros(nMin,nModels,nSites)./0;
assimGlobalStats = cell(nModels,nSites);
assimGlobalSens  = zeros(2,nModels,nSites)./0;
assimGlobalNSE   = zeros(nModels,nSites)./0;
assimGlobalNMI   = zeros(nModels,nSites)./0;

% init storage - open loop only
openGlobalPreds = zeros(nMin,nModels,nSites)./0;
openGlobalStats = cell(nModels,nSites);
openGlobalSens  = zeros(1,nModels,nSites)./0;
openGlobalNSE   = zeros(nModels,nSites)./0;
openGlobalNMI   = zeros(nModels,nSites)./0;

% init storage - remote sensing only
lprmGlobalPreds = zeros(nMin,nModels,nSites)./0;
lprmGlobalStats = cell(nModels,nSites);
lprmGlobalSens  = zeros(1,nModels,nSites)./0;
lprmGlobalNSE   = zeros(nModels,nSites)./0;
lprmGlobalNMI   = zeros(nModels,nSites)./0;


% run the analysis
for s = 1:nSites
    
    % screen report
    fprintf('Training at site %d of %d ... ',s,nSites); tic;
    
    for m = 1:nModels
        
        % extract datam
        Is = 1:nSites; Is(s) = [];
        Dtrain = reshape(allData(:,Is,:),[nMin*(nSites-1),3]);
        Xtrain = Dtrain(:,[1,2]);
        Ytrain = Dtrain(:,end); 
        
        Xtest = squeeze(allData(:,s,[1,2]));
        Ytest = squeeze(allData(:,s,end));
        
        % train k-fold models
        [assimGlobalPreds(:,m,s),assimGlobalStats{m,s},assimGlobalSens(:,m,s)] = ...
            calval_benchmark(Xtrain,Ytrain,Xtest,Ytest,modelTypes{m},infoBinEdges,infoBinEdges,verbose);
        [openGlobalPreds(:,m,s),openGlobalStats{m,s},openGlobalSens(:,m,s)] = ...
            calval_benchmark(Xtrain(:,1),Ytrain,Xtest(:,1),Ytest,modelTypes{m},infoBinEdges,infoBinEdges,verbose);
        [lprmGlobalPreds(:,m,s),lprmGlobalStats{m,s},lprmGlobalSens(:,m,s)] = ...
            calval_benchmark(Xtrain(:,2),Ytrain,Xtest(:,2),Ytest,modelTypes{m},infoBinEdges,infoBinEdges,verbose);
        
        assimGlobalNSE(m,s) = assimGlobalStats{m,s}.test.nse;
        openGlobalNSE(m,s) = openGlobalStats{m,s}.test.nse;
        lprmGlobalNSE(m,s) = lprmGlobalStats{m,s}.test.nse;
        
        assimGlobalNMI(m,s) = assimGlobalStats{m,s}.test.nmi;
        openGlobalNMI(m,s) = openGlobalStats{m,s}.test.nmi;
        lprmGlobalNMI(m,s) = lprmGlobalStats{m,s}.test.nmi;
        
    end % m-loop
    
    % screen report
    fprintf('. finished; time = %f[s]. \n',toc);

    formatstr = strcat({'--------------  '},repmat('%s\t ',[1,nModels]),{' -------------- '},repmat('%s\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        modelTypesShort{:},modelTypesShort{:});

    formatstr = strcat({'--- NSE Raw  = '},repmat('%5.3f\t ',[1,nModels]),{' --- NMI Raw  = '},repmat('%5.3f\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        repmat(nseSites(s,1),[1,nModels]),repmat(nmiSites(s,1),[1,nModels]));

    formatstr = strcat({'--- NSE Both = '},repmat('%5.3f\t ',[1,nModels]),{' --- NMI Both = '},repmat('%5.3f\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        assimGlobalNSE(:,s),assimGlobalNMI(:,s));
    
    formatstr = strcat({'--- NSE Open = '},repmat('%5.3f\t ',[1,nModels]),{' --- NMI Open = '},repmat('%5.3f\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        openGlobalNSE(:,s),openGlobalNMI(:,s));
    
    formatstr = strcat({'--- NSE LPRM = '},repmat('%5.3f\t ',[1,nModels]),{' --- NMI LPRM = '},repmat('%5.3f\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        lprmGlobalNSE(:,s),lprmGlobalNMI(:,s));
    
    fprintf('------------------------------------------------------------------------------------\n');

end % s-loop

% screen report
fprintf('--- All Global Models Finished ---------------------- \n');

%% --- Train Site-Specific Models -----------------------------------------

% screen report
fprintf('--- Training Site-Specific Models ------------------- \n');

% init storage - assimilation
assimSitePreds = zeros(nMin,nModels,nSites)./0;
assimSiteStats = cell(nModels,nSites);
assimSiteSens  = zeros(2,kFold,nModels,nSites)./0;
assimSiteNSE   = zeros(nModels,nSites)./0;
assimSiteNMI   = zeros(nModels,nSites)./0;

% init storage - open loop only
openSitePreds = zeros(nMin,nModels,nSites)./0;
openSiteStats = cell(nModels,nSites);
openSiteSens  = zeros(1,kFold,nModels,nSites)./0;
openSiteNSE   = zeros(nModels,nSites)./0;
openSiteNMI   = zeros(nModels,nSites)./0;

% init storage - remote sensing only
lprmSitePreds = zeros(nMin,nModels,nSites)./0;
lprmSiteStats = cell(nModels,nSites);
lprmSiteSens  = zeros(1,kFold,nModels,nSites)./0;
lprmSiteNSE   = zeros(nModels,nSites)./0;
lprmSiteNMI   = zeros(nModels,nSites)./0;

% run the analysis
for s = 1:nSites
    
    % screen report
    fprintf('Training at site %d of %d ... ',s,nSites); tic;
    
    for m = 1:nModels
        
        % extract data
        X = squeeze(allData(:,s,[1,2]));
        Y = squeeze(allData(:,s,end));
        
        % train k-fold models
        [assimSitePreds(:,m,s),assimSiteStats{m,s},assimSiteSens(:,:,m,s)] = ...
            kfold_benchmark(X,Y,kFold,modelTypes{m},infoBinEdges,infoBinEdges,verbose);
        [openSitePreds(:,m,s),openSiteStats{m,s},openSiteSens(:,:,m,s)] = ...
            kfold_benchmark(X(:,1),Y,kFold,modelTypes{m},infoBinEdges,infoBinEdges,verbose);
        [lprmSitePreds(:,m,s),lprmSiteStats{m,s},lprmSiteSens(:,:,m,s)] = ...
            kfold_benchmark(X(:,2),Y,kFold,modelTypes{m},infoBinEdges,infoBinEdges,verbose);
        
        assimSiteNSE(m,s) = assimSiteStats{m,s}.all.nse;
        openSiteNSE(m,s) = openSiteStats{m,s}.all.nse;
        lprmSiteNSE(m,s) = lprmSiteStats{m,s}.all.nse;
        
        assimSiteNMI(m,s) = assimSiteStats{m,s}.all.nmi;
        openSiteNMI(m,s) = openSiteStats{m,s}.all.nmi;
        lprmSiteNMI(m,s) = lprmSiteStats{m,s}.all.nmi;
        
    end % m-loop
    
    % screen report
    fprintf('. finished; time = %f[s]. \n',toc);

    formatstr = strcat({'--------------  '},repmat('%s\t ',[1,nModels]),{' -------------- '},repmat('%s\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        modelTypesShort{:},modelTypesShort{:});

    formatstr = strcat({'--- NSE Raw  = '},repmat('%5.3f\t ',[1,nModels]),{' --- NMI Raw  = '},repmat('%5.3f\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        repmat(nseSites(s,1),[1,nModels]),repmat(nmiSites(s,1),[1,nModels]));

    formatstr = strcat({'--- NSE Both = '},repmat('%5.3f\t ',[1,nModels]),{' --- NMI Both = '},repmat('%5.3f\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        assimSiteNSE(:,s),assimSiteNMI(:,s));
    
    formatstr = strcat({'--- NSE Open = '},repmat('%5.3f\t ',[1,nModels]),{' --- NMI Open = '},repmat('%5.3f\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        openSiteNSE(:,s),openSiteNMI(:,s));
    
    formatstr = strcat({'--- NSE LPRM = '},repmat('%5.3f\t ',[1,nModels]),{' --- NMI LPRM = '},repmat('%5.3f\t ',[1,nModels]),{' ---\n'});
    fprintf(formatstr{1},...
        lprmSiteNSE(:,s),lprmSiteNMI(:,s));
    
    fprintf('------------------------------------------------------------------------------------\n');

end % s-loop

% screen report
fprintf('--- All Site-Specific Models Finished --------------- \n');

%% --- Plot Site-Specific Stats -------------------------------------------

O = allData(:,:,end);
M = allData(:,:,1);
L = allData(:,:,2);

% calculate global stats
for m = 1:nModels
    
    A = squeeze(assimGlobalPreds(:,m,:));

    nseGlobal(m+2) = nse(A(:),O(:));
    rhoGlobal(m+2) = rho(A(:),O(:));
    nmiGlobal(m+2) = nmi(A(:),O(:),infoBinEdges,infoBinEdges);
    
end % m-loop

nseGlobal
rhoGlobal
nmiGlobal

% % % calcualte non-redundant info
% % dxs =  0.01:0.001:0.10;
% % for dx = 1:length(dxs)
% %     infoBinEdges = -1:dxs(dx):2;
% %     
% %     O = allData(:,:,end);
% %     M = allData(:,:,1);
% %     L = allData(:,:,2);
% %     [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(O,M,L,infoBinEdges,infoBinEdges,infoBinEdges);
% %     sumInfo(dx) = (Ixyz);
% %     totalInfo(dx) = (Ixy + Ixz - Ixyz) / Hx
% %     
% % end
% % 
% % plot(dxs,sumInfo)
% 
% O = allData(:,:,end);
% M = allData(:,:,1);
% L = allData(:,:,2);
% 
% % calculate global stats
% for m = 1:nModels
%     
%     A = squeeze(assimGlobalPreds(:,m,:));
% 
%     nseGlobal(m+2) = nse(A(:),O(:));
%     rhoGlobal(m+2) = rho(A(:),O(:));
%     
% %     [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(O,A,L,infoBinEdges,infoBinEdges,infoBinEdges);
% %     [Ixy/Hx,Ixz/Hx]
% % 
% %     [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(O,A,M,infoBinEdges,infoBinEdges,infoBinEdges);
% %     [Ixy/Hx,Ixz/Hx]
% 
% %     nmiGlobal(m+2) = Ixy/Hx; %nmi(A(:),O(:),infoBinEdges,infoBinEdges);
%     I = randperm(numel(A),1e4);
%     nmiGlobal(m+2) = knn_info(A(I)',O(I)',5);
%     nmiGlobal(1)   = knn_info(M(I)',O(I)',5);
% %     efficiency(m) = (nmiGlobal(m+2) - Ixz/Hx) / (totalInfo - - Ixz/Hx);
%     
% end % m-loop
% 
% nseGlobal
% rhoGlobal
% nmiGlobal


%% --- END SCRIPT ---------------------------------------------------------