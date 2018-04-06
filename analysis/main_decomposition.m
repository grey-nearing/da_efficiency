clear all; close all; clc
% clearvars -except outdata; close all clc
restoredefaultpath; addpath(genpath(pwd))
fignum = 0;

%% ************************************************************************
%% --- Experiment Parameters ----------------------------------------------
%% ************************************************************************

% minimum # of data points to use from each site
Nmin = 1200;

% load model and obs covariacnes for axis labels
Ocov = load('../enkf/obs_perts_master.txt');
Mcov = load('../enkf/state_perts_master.txt');

% bins for informmation calculations
Bw = 0.03;                  % bin width
B = -1:Bw:2;                % bin edges
Bc = B + 0.5*(B(2)-B(1));   % bin centers

%% ************************************************************************
%% --- Load Raw Data ------------------------------------------------------
%% ************************************************************************

% load data
load('./data/site_data_1.mat');

% dimensions
Isites = find(outdata.Ndata(:,1,1) >= Nmin);
Ndata = outdata.Ndata(Isites,:,:);
[Nsites,Nobs,Nmod] = size(Ndata);

% init storage
val = zeros(Nmin,Nsites); 
obs = zeros(Nmin,Nsites); 
cdf = zeros(Nmin,Nsites); 
mod = zeros(Nmin,Nobs,Nmod,Nsites); 
ens = zeros(Nmin,Nobs,Nmod,Nsites); 

% loop through sites
for s = 1:Nsites; tic

    % load data from sites
    load(strcat('./data/site_data_',num2str(Isites(s)),'.mat'));
            
    % pull sites with enough data
    val(:,s)       = outdata.val  (1:Nmin);
    obs(:,s)       = outdata.obs  (1:Nmin);
    cdf(:,s)       = outdata.cdf  (1:Nmin);
    mod(:,:,:,s)   = outdata.Mback(1:Nmin,:,:);
    ens(:,:,:,s)   = outdata.Sback(1:Nmin,:,:);
    
    % screen report
    fprintf('loaded data from site %d of %d; time = %f \n',s,Nsites,toc);
    
end

mod = permute(mod,[1,4,2,3]);
ens = permute(ens,[1,4,2,3]);

% find missing enkf runs
missing = zeros(Nobs,Nmod);
for o = 1:Nobs
    for m = 1:Nmod
        assert(min(mod(:))>0); 
        assert(max(mod(:))<=1);
        a = isnan(mod(:,:,o,m));
        if any(a(:))
            missing(o,m) = 1;
        end
    end
end
i = find(all(missing)); mod(:,:,:,i) = []; ens(:,:,:,i) = []; Mcov(i) = [];
i = find(all(missing')); mod(:,:,i,:) = []; ens(:,:,i,:) = []; Ocov(i) = [];
        
% dimensions
[Ntimes,Nsites,Nobs,Nmod] = size(mod);

% make sure there are no missing data
assert(all(~isnan(mod(:))))
assert(all(~isnan(ens(:))))
assert(all(~isnan(val(:))))
assert(all(~isnan(obs(:))))
assert(all(~isnan(cdf(:))))

% make sure all data is in a feasible range
assert(min(val(:))>0); assert(max(val(:))<=1);
assert(min(obs(:))>0); assert(max(obs(:))<=1);
assert(min(cdf(:))>0); assert(max(cdf(:))<=1);
assert(min(mod(:))>0); assert(max(mod(:))<=1);
assert(min(ens(:))>0); assert(max(ens(:))<=1);

% make sure there are no missing values
assert(min(val(:))>0); assert(max(val(:))<=1); assert(~any(isnan(val(:))));
assert(min(mod(:))>0); assert(max(mod(:))<=1); assert(~any(isnan(mod(:))));
assert(min(ens(:))>0); assert(max(ens(:))<=1); assert(~any(isnan(ens(:))));
assert(min(obs(:))>0); assert(max(obs(:))<=1); assert(~any(isnan(obs(:))));
assert(min(cdf(:))>0); assert(max(cdf(:))<=1); assert(~any(isnan(cdf(:))));

% load cdf covariances
for s = 1:Nsites
    fname = strcat('../data/lprm/site_data/cdf_sig_',num2str(Isites(s)),'.txt');
    Ccov(s) = load(fname);
end 
Ccov = repmat(Ccov,[Nmin,1]);
Ccov = Ccov(:);

% turn into single vectors over all sites
% Y  = obs(:); % observation data
Y  = cdf(:);   % cdf-matched obs
Z  = val(:);   % validation data

% number of samples
N = length(Y);
assert(length(Z) ==N)

% correction factor
cf = 1/N;

% init storage
Dd_store  = zeros(Nobs,Nmod)./0;
Dm_store  = zeros(Nobs,Nmod)./0;
Df_store  = zeros(Nobs,Nmod)./0;
Dp_store  = zeros(Nobs,Nmod)./0;
Dpm_store = zeros(Nobs,Nmod)./0;
Dpv_store = zeros(Nobs,Nmod)./0;
Dl_store  = zeros(Nobs,Nmod)./0;
Dlm_store = zeros(Nobs,Nmod)./0;
Dlv_store = zeros(Nobs,Nmod)./0;

%% ************************************************************************
%% --- Calculate Data Histograms ------------------------------------------
%% ************************************************************************

for o = 16:Nobs
    for m = 1:Nmod
        
        % screen report
        fprintf('Calculating divergences @ obs %d/%d, mod %d/%d ...',o,Nobs,m,Nmod); tic;
        if missing(o,m); continue; end;
        
        % pull data
        Xb = mod(:,:,o,m); Xb = Xb(:);  % posterior mean
        Qb = ens(:,:,o,m); Qb = Qb(:);  % posterior stdev
        
        % data histograms
        [Pxyz,Pxy,Pxz,Pyz,Px,Py,Pz] = hist3(Xb,Y,Z,B,B,B);
        %         Pxyz = max(Pxyz,cf); Pxyz = Pxyz / sum(Pxyz(:));
        %         Pxy  = max(Pxy,cf);  Pxy  = Pxy  / sum(Pyz (:));
        %         Pxz  = max(Pxz,cf);  Pxz  = Pxz  / sum(Pxz (:));
        %         Pyz  = max(Pyz,cf);  Pyz  = Pyz  / sum(Pyz (:));
        %         Px   = max(Px,cf);   Px   = Px   / sum(Px  (:));
        %         Py   = max(Py,cf);   Py   = Py   / sum(Py  (:));
        %         Pz   = max(Pz,cf);   Pz   = Pz   / sum(Pz  (:));
        
        % data prior
        Pz_x = zeros(length(B))./0;
        for x = 1:length(B)
            for z = 1:length(B)
                if Pxz(x,z) == 0
                    Pz_x(x,z) = 0;
                elseif Px(x) == 0
                    error('Pxz ~= 0 && Px == 0');
                else
                    Pz_x(x,z) = Pxz(x,z)/Px(x);
                end
            end
        end
        assert(isempty(find(Pz_x(:) > 1,1)));
        assert(isempty(find(Pz_x(:) < 0,1)));
        assert(isempty(find(abs(sum(Pz_x,1)-1)<1e-15 & abs(sum(Pz_x,1))<1e-15,1)));
        
        % data likelihood
        Py_xz = zeros(length(B),length(B),length(B))./0;
        for x = 1:length(B)
            for y = 1:length(B)
                for z = 1:length(B)
                    if Pxyz(x,y,z) == 0
                        Py_xz(x,y,z) = 0;
                    elseif Pxz(x,z) == 0
                        error('Pxyz ~= 0 && Pxz == 0');
                    else
                        Py_xz(x,y,z) = Pxyz(x,y,z)/Pxz(x,z);
                    end
                end
            end
        end
        assert(isempty(find(Py_xz(:) > 1,1)));
        assert(isempty(find(Py_xz(:) < 0,1)));
        assert(isempty(find(abs(sum(Py_xz,2)-1)<1e-15 & abs(sum(Py_xz,2))<1e-15,1)));
        
        % Markov likelihood
        Py_z = zeros(length(B))./0;
        for y = 1:length(B)
            for z = 1:length(B)
                if Pyz(y,z) == 0
                    Py_z(y,z) = 0;
                elseif Pz(z) == 0
                    error('Pxy ~= 0 && Px == 0');
                else
                    Py_z(y,z) = Pyz(y,z)/Pz(z);
                end
            end
        end
        assert(isempty(find(Py_z(:) > 1,1)));
        assert(isempty(find(Py_z(:) < 0,1)));
        assert(isempty(find(abs(sum(Py_z,2)-1)<1e-15 & abs(sum(Py_z,2))<1e-15,1)));
        
        % data posterior
        Pz_xy = zeros(length(B),length(B),length(B))./0;
        for x = 1:length(B)
            for y = 1:length(B)
                for z = 1:length(B)
                    if Pxyz(x,y,z) == 0
                        Pz_xy(x,y,z) = 0;
                    elseif Pxy(x,y) == 0
                        error('Pxyz ~= 0 && Pxy == 0');
                    else
                        Pz_xy(x,y,z) = Pxyz(x,y,z)/Pxy(x,y);
                    end
                end
            end
        end
        assert(isempty(find(Pz_xy(:) > 1,1)));
        assert(isempty(find(Pz_xy(:) < 0,1)));
        assert(isempty(find(abs(sum(Pz_xy,3)-1)<1e-15 & abs(sum(Pz_xy,3))<1e-15,1)));
        
        % init storage
        Dd  = zeros(N,length(B))./0;
        Dm  = zeros(N,length(B))./0;
        Df  = zeros(N,length(B))./0;
        Dp  = zeros(N,length(B))./0;
        Dpm = zeros(N,length(B))./0;
        Dpv = zeros(N,length(B))./0;
        Dl  = zeros(N,length(B))./0;
        Dlm = zeros(N,length(B))./0;
        Dlv = zeros(N,length(B))./0;
        
        % MC integration over divergences for an info metric:
        % This integration strategy is necessary becaseu the background
        % variances are not constant.
        for n = 1:N
            
            % current data point
            x = find(Xb(n)>=B,1,'last');      % index into data posterior
            y = find(Y(n) >=B,1,'last');      % index into data retrieval
            
            % --- pure data pdfs ------------------------------------------
            % full data posterior by computation
            pd = Pz_x(x,:).*squeeze(Py_xz(x,y,:))';             % grab posterior
            pd = pd./nansum(pd);                                % normalize
            assert(isempty(find(isnan(pd),1)));
            assert(max(abs(pd-squeeze(Pz_xy(x,y,:))'))<1e-15);  % check agreement
            
            % Markov assumption
            pm = Pz_x(x,:).*Py_z(y,:);                          % grab posterior
            pm = pm./nansum(pm);                                % normalize
            assert(isempty(find(isnan(pm),1)));
            
            % pure enkf
            Fz_x = max(cf,pdf('Normal',Bc,Xb(x),Qb(n)));
            Fz_x = Fz_x./sum(Fz_x);
            Fy_z = max(cf,pdf('Normal',Y(y),Bc,Ocov(o))*Ccov(n));
            Fy_z = Fy_z./sum(Fy_z);
            pf = Fz_x.*Fy_z;                                    % grab posterior
            pf = pf./nansum(pf);                                % normalize
            assert(isempty(find(isnan(pf),1)));
            
            % --- enkf priors ---------------------------------------------
            % calculate the sample mean and varaince
            pmu  = nansum(     Bc   .*squeeze(Pz_x(x,:)));
            psig = sqrt(nansum((Bc-pmu).^2.*squeeze(Pz_x(x,:))));
            
            % full enkf prior
            Fz_x = max(cf,pdf('Normal',Bc,Xb(x),Qb(n)));
            Fz_x = Fz_x./sum(Fz_x);
            pp = Fz_x.*Py_z(y,:);                               % grab posterior
            pp = pp./nansum(pp);                                % normalize
            assert(isempty(find(isnan(pp),1)));
            
            % enkf prior mean
            Fz_x = max(cf,pdf('Normal',Bc,Xb(x),psig));
            Fz_x = Fz_x./sum(Fz_x);
            ppm = Fz_x.*Py_z(y,:);                              % grab posterior
            ppm = ppm./nansum(ppm);                             % normalize
            assert(isempty(find(isnan(ppm),1)));
            
            % enkf prior var
            Fz_x = max(cf,pdf('Normal',Bc,pmu,Qb(n)));
            Fz_x = Fz_x./sum(Fz_x);
            ppv = Fz_x.*Py_z(y,:);                              % grab posterior
            ppv = ppv./nansum(ppv);                             % normalize
            assert(isempty(find(isnan(ppv),1)));
            
            % --- enkf likelihoods ----------------------------------------
            % calculate the sample mean and varaince
            lmu = nansum(repmat(Bc',[1,length(B)]).*Py_z);
            lsig = sqrt(nansum(repmat((Bc-lmu(z))'.^2,[1,length(B)]).*Py_z));
            
            % full enkf likelihood
            Fy_z = max(cf,pdf('Normal',Y(y),Bc,Ocov(o))*Ccov(n));
            Fy_z = Fy_z./sum(Fy_z);
            pl = Pz_x(x,:).*Fy_z;                               % grab posterior
            pl = pl./nansum(pl);                                % normalize
            assert(isempty(find(isnan(pl),1)));
            
            % enkf likelihood mean
            Fy_z = max(cf,pdf('Normal',Y(y),Bc,lsig));
            Fy_z = Fy_z./sum(Fy_z);
            plm = Pz_x(x,:).*Fy_z;                              % grab posterior
            plm = plm./nansum(plm);                             % normalize
            assert(isempty(find(isnan(plm),1)));
            
            % enkf likelihood variance
            Fy_z = max(cf,pdf('Normal',Y(y),lmu,Ocov(o))*Ccov(n));
            Fy_z = Fy_z./sum(Fy_z);
            plv = Pz_x(x,:).*Fy_z;                              % grab posterior
            plv = plv./nansum(plv);                             % normalize
            assert(isempty(find(isnan(plv),1)));
            
            % --- divergences ---------------------------------------------
            % calculate divergences
            Dd (n,:) = pd.*log(pd./squeeze(Pz_xy(x,y,:))');
            Dm (n,:) = pd.*log(pd./pm);
            Df (n,:) = pd.*log(pd./pf);
            Dp (n,:) = pd.*log(pd./pp);
            Dpm(n,:) = pd.*log(pd./ppm);
            Dpv(n,:) = pd.*log(pd./ppv);
            Dl (n,:) = pd.*log(pd./pl);
            Dlm(n,:) = pd.*log(pd./plm);
            Dlv(n,:) = pd.*log(pd./plv);
            
            % deal with grandmas
            Dd (n,pd==0) = 0;
            Dm (n,pd==0) = 0;
            Df (n,pd==0) = 0;
            Dp (n,pd==0) = 0;
            Dpm(n,pd==0) = 0;
            Dpv(n,pd==0) = 0;
            Dl (n,pd==0) = 0;
            Dlm(n,pd==0) = 0;
            Dlv(n,pd==0) = 0;
            
            % deal with grandmas
            Dd (n,pd==0) = 0;
            Dm (n,pm==0) = 0;
            Df (n,pf==0) = 0;
            Dp (n,pp==0) = 0;
            Dpm(n,ppm==0) = 0;
            Dpv(n,ppv==0) = 0;
            Dl (n,pl==0) = 0;
            Dlm(n,plm==0) = 0;
            Dlv(n,plv==0) = 0;

            % error checking
            assert(isempty(find(isnan(Dd (n,:)),1)));
            assert(isempty(find(isnan(Dm (n,:)),1)));
            assert(isempty(find(isnan(Df (n,:)),1)));
            assert(isempty(find(isnan(Dp (n,:)),1)));
            assert(isempty(find(isnan(Dpm(n,:)),1)));
            assert(isempty(find(isnan(Dpv(n,:)),1)));
            assert(isempty(find(isnan(Dl (n,:)),1)));
            assert(isempty(find(isnan(Dlm(n,:)),1)));
            assert(isempty(find(isnan(Dlv(n,:)),1)));
            
        end % monte carlo loop
        
        % integrate divergences
        Dd_store(o,m)  = mean(Dd(:) );
        Dm_store(o,m)  = mean(Dm(:) );
        Df_store(o,m)  = mean(Df(:) );
        Dp_store(o,m)  = mean(Dp(:) );
        Dpm_store(o,m) = mean(Dpm(:));
        Dpv_store(o,m) = mean(Dpv(:));
        Dl_store(o,m)  = mean(Dl(:) );
        Dlm_store(o,m) = mean(Dlm(:));
        Dlv_store(o,m) = mean(Dlv(:));
        
        % screen report
        t = toc; fprintf('. finished; t = %f \n',t);
        
        disp([Dd_store(o,m),Dm_store(o,m) ,Df_store(o,m) ]./Df_store(o,m));
        disp([Dp_store(o,m),Dpm_store(o,m),Dpv_store(o,m)]./Df_store(o,m));
        disp([Dl_store(o,m),Dlm_store(o,m),Dlv_store(o,m)]./Df_store(o,m));
        
    end
end

%% ************************************************************************
%% --- save results -------------------------------------------------------
%% ************************************************************************

% load efficiency metrics
% load('./data/effficiency_results.mat');

% gather everything into a matrix
% allResults(:,:,1) = Eresults.Itot-Eresults.Ifil;
allResults(:,:,2) = Df_store;% ./ allResults(:,:,1);
allResults(:,:,3) = Dm_store;% ./ allResults(:,:,1);%./Df_store;

allResults(:,:,4) = Dpm_store;% ./allResults(:,:,1);%./Df_store;
allResults(:,:,5) = Dpv_store;% ./allResults(:,:,1);%./Df_store;
allResults(:,:,6) = Dp_store ;% ./allResults(:,:,1);%./Df_store;

allResults(:,:,7) = Dlm_store;% ./allResults(:,:,1);%./Df_store;
allResults(:,:,8) = Dlv_store;% ./allResults(:,:,1);%./Df_store;
allResults(:,:,9) = Dl_store ;% ./allResults(:,:,1);%./Df_store;

% allResults = allResults./numel(Dd);
% allResults(:,:,1) = Eresults.Itot-Eresults.Ifil;
% allResults = allResults./repmat(allResults(:,:,1),[1,1,9]);
% allResults(:,:,1) = Eresults.Itot-Eresults.Ifil;


% save
save('./data/DecompResults.mat','allResults','-v7.3');

%% ************************************************************************
%% --- plot results -------------------------------------------------------
%% ************************************************************************

% figure names
figNames = [{'E_Y'}         ,{'D_E_n_K_F'}              ,{'D_M_a_r_k_o_v'},...
    {'D_p_r_i_o_r_,_\mu'}   ,{'D_p_r_i_o_r_,_\sigma'}   ,{'D_p_r_i_o_r'} ,...
    {'D_l_i_k_e_,_\mu'}     ,{'D_l_i_k_e_,_\sigma'}     ,{'D_l_i_k_e'}];

% init figure
% fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[295          15        2197        1483]);

for s = 2:9
    
    subplot(3,3,s)
    surf(Mcov,Ocov,allResults(:,:,s));
    set(gca,'XScale','log','YScale','log')
    view([1.5850e+02   3.0800e+01]);
    set(gca,'fontsize',16);
    set(gca,'ylim',[min(Ocov),max(Ocov)],'xlim',[min(Mcov),max(Mcov)]);
    ylabel('Obs. Cov. Scale Factor','fontsize',18)
    xlabel('Model State Cov.','fontsize',18)
    title(figNames{s},'fontsize',26);
%     set(gca,'zlim',[0,2e6]);
%     caxis([0,2e6])
%     view(2)
%     if s > 2; caxis([0,1]); end;
%     if s == 1; colormap gray; end
%     if s > 1; colormap parula; end
    c = colorbar;
    c.YLabel.String = '[nats]';
    c.YLabel.FontSize = 16;
    
end


%% ************************************************************************
%% --- END PROGRAM --------------------------------------------------------
%% ************************************************************************

