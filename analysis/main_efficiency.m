clear all; close all; clc
restoredefaultpath; addpath(genpath(pwd))

%% --- Experiment Parameters ----------------------------------------------

% bins - based on the convergence analysis
B = [0:0.03:1.03];

% # of data points to use from each site
Nmin = 1200;

%% --- Load Raw Data ------------------------------------------------------

% load data
load('./data/site_data_1.mat');

% dimensions
Is = find(outdata.Ndata(:,1,1) >= Nmin);
Ndata = outdata.Ndata(Is,:,:);
[Nsites,Nobs,Nmod] = size(Ndata);

% init storage
val = zeros(Nmin,Nsites); 
mod = zeros(Nmin,Nsites); 
obs = zeros(Nmin,Nsites); 
cdf = zeros(Nmin,Nsites); 
fil = zeros(Nmin,Nobs,Nmod,Nsites); 

% loop through sites
for s = 1:Nsites; tic

    % load data from sites
    load(strcat('./data/site_data_',num2str(Is(s)),'.mat'));
            
    % pull sites with enough data
    val(:,s)     = outdata.val  (1:Nmin);
    mod(:,s)     = outdata.Mopen(1:Nmin);
    obs(:,s)     = outdata.obs  (1:Nmin);
    cdf(:,s)     = outdata.cdf  (1:Nmin);
    fil(:,:,:,s) = outdata.Menkf(1:Nmin,:,:);
    
    % make sure there are no missing values
    assert(min(val(:,s))>0); assert(max(val(:,s))<=1); assert(~any(isnan(val(:,s))));
    assert(min(mod(:,s))>0); assert(max(mod(:,s))<=1); assert(~any(isnan(mod(:,s))));
    assert(min(obs(:,s))>0); assert(max(obs(:,s))<=1); assert(~any(isnan(obs(:,s))));
    assert(min(cdf(:,s))>0); assert(max(cdf(:,s))<=1); assert(~any(isnan(cdf(:,s))));
%     assert(min(fil(:,:,:,s))>0); assert(max(fil(:,:,:,s))<=1); assert(~any(isnan(fil(:,:,:,s))));
    
    % screen report
    fprintf('loaded data from site %d of %d; time = %f \n',s,Nsites,toc);
    
end

fil = permute(fil,[1,4,2,3]);

% make sure there are no missing values
assert(min(val(:))>0); assert(max(val(:))<=1); assert(~any(isnan(val(:))));
assert(min(mod(:))>0); assert(max(mod(:))<=1); assert(~any(isnan(mod(:))));
assert(min(obs(:))>0); assert(max(obs(:))<=1); assert(~any(isnan(obs(:))));
assert(min(cdf(:))>0); assert(max(cdf(:))<=1); assert(~any(isnan(cdf(:))));

% find missing enkf runs
missing = zeros(Nobs,Nmod);
for o = 1:Nobs
    for m = 1:Nmod
        assert(min(fil(:))>0); 
        assert(max(fil(:))<=1);
        a = isnan(fil(:,:,o,m));
        if any(a(:))
            missing(o,m) = 1;
        end
    end
end

%% --- Calculate Efficiency Metrics ---------------------------------------

% --- calculate static parts ---------

% pull a sample
X = mod;
Y = obs;
C = cdf;
Z = val;

% total info stats (raw obs)
[Ixyz,~,Ixz,Iyz,~,~,Hz] = mutual_info_3(X(:),Y(:),Z(:),B,B,B);
Imod = Ixz/Hz;                         % info in model
Iobs = Iyz/Hz;                         % info in raw obs
Icon = (Iyz-Ixyz)/Hz;                  % conditional obs info
Itot = (Ixz+Iyz-Ixyz)/Hz;              % total info

% total info stats (cdf obs)
[Ixyz,~,~,Iyz,~,~,Hz] = mutual_info_3(X(:),C(:),Z(:),B,B,B);
Iobs_cdf = Iyz/Hz/Hz;                  % info cdf obs
Icon_cdf = (Iyz-Ixyz)/Hz;              % info cdf obs

% information loss via cdf matching
Ecdf     = (Iobs - Iobs_cdf)/Iobs;
Ecdf_con = (Icon - Icon_cdf)/Icon;

% rmse stats
rmseX = sqrt(mean((X(:)-Z(:)).^2));
rmseY = sqrt(mean((Y(:)-Z(:)).^2));
rmseC = sqrt(mean((C(:)-Z(:)).^2));

% corr stats
cc = corrcoef(X(:)-mean(X(:)),Z(:)-mean(Z(:))); acorX = cc(2);
cc = corrcoef(Y(:)-mean(X(:)),Z(:)-mean(Z(:))); acorY = cc(2);
cc = corrcoef(C(:)-mean(X(:)),Z(:)-mean(Z(:))); acorC = cc(2);

% --- calclate parts that change with EnKF variables

% timing variable
etime = 0;
t = 0;

% loop through model and obs covartiance experiments
for m = 1:Nmod
    for o = 1:Nobs
        
        % dont run missing enkf configurations
        if missing(o,m); continue; end;
        
        % screen report
        tic;
        fprintf('running m = %d, o = %d ...',m,o);
        
        % pull a sample
        F = fil(:,:,o,m);
        
        % total info stats (enkf)
        Ifil(o,m) = info(F(:),Z(:),B,B,2);          % info in enkf
        
        % rmse stats
        rmseF(o,m) = sqrt(mean((F(:)-Z(:)).^2));
        
        % corr stats
        cc = corrcoef(F(:)-mean(F(:)),Z(:)-mean(Z(:))); acorF(o,m) = cc(2);

        % screen report
        time = toc;
        etime = etime + time;
        atime = etime/t;
        rtime = (Nobs*Nmod-t)*atime;
        fprintf('. time = %f \n',time);
        
    end
end

% efficiency metrics
Ey = (Ifil - Imod) ./ Icon;
Eda  = Ifil ./ Itot;

%% --- Save Results -------------------------------------------------------

% % store in output structure
% Eresults.Imod = Imod;
% Eresults.Iobs = Iobs;
% Eresults.Icdf = Icdf;
% Eresults.Icon = Icon;
% Eresults.Itot = Itot;
% Eresults.Ifil = Ifil;
% Eresults.Eda  = Eda;
% Eresults.Ey   = Ey;
% Eresults.Ecdf = Ecdf;
% 
% % save to file
% save('./data/effficiency_results.mat','Eresults','-v7.3');

%% --- Plot Results -------------------------------------------------------

% init figure counter
fignum  = 0;

% load model and obs covariacnes for axis labels
Ocov = load('../enkf/obs_perts_master.txt');
Mcov = load('../enkf/state_perts_master.txt')'; 
Mcov(13:end) = [];

% init figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[295         704        3222         794]);

% plot Eda results
subplot(1,3,1);
surf(Ocov,Mcov,Eda');
set(gca,'YScale','log');
set(gca,'XScale','log');
view([1.5850e+02   3.0800e+01]);
set(gca,'fontsize',16);
set(gca,'xlim',[min(Ocov),max(Ocov)],'ylim',[min(Mcov),max(Mcov)]);
xlabel('Obs. Cov. Scale Factor','fontsize',18)
ylabel('Model State Cov.','fontsize',18)
title('E_D_A Results','fontsize',26);
colorbar

% plot Ey results
subplot(1,3,2);
surf(Ocov,Mcov,Ey');
set(gca,'YScale','log');
set(gca,'XScale','log');
set(gca,'xlim',[min(Ocov),max(Ocov)],'ylim',[min(Mcov),max(Mcov)]);%,'zlim',[0,0.02]);
view([1.5850e+02   3.0800e+01]);
set(gca,'fontsize',16);
xlabel('Obs. Cov. Scale Factor','fontsize',18)
ylabel('Model State Cov.','fontsize',18)
title('E_Y Results','fontsize',26);
colorbar

% plot RMSE results
subplot(1,3,3);
%surf(Ocov,Mcov,rmseF./rmseX);
surf(Ocov,Mcov,(acorF-acorX)'./acorX');
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'xlim',[min(Ocov),max(Ocov)],'ylim',[min(Mcov),max(Mcov)]);%,'zlim',[-1,1]);
view([1.5850e+02   3.0800e+01]);
set(gca,'fontsize',16);
xlabel('Obs. Cov. Scale Factor','fontsize',18)
ylabel('Model State Cov.','fontsize',18)
title('Anomaly Corr. Results','fontsize',26);
colorbar

% save plot
fname = strcat('figures/Efficiency_Sensitivity');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- END SCRIPT ---------------------------------------------------------

