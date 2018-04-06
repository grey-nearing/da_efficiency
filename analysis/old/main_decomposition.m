clear all
close all
clc
restoredefaultpath
addpath(genpath(pwd))

if 0
%% *************************************************************************************
%% --- runtime parameters
%% *************************************************************************************

 M = 1000;
 R = 0.04;
 B = [-1:0.02:2];

%% *************************************************************************************
%% --- load data 
%% *************************************************************************************

 % screen report
 fprintf('Loading data ... \n'); tic;

 % load tabular data from previous m-file
 load('data/all_data.mat');

 Ndata = outdata.Ndata;
 Mopen = outdata.Mopen;
 Menks = outdata.Menks;
 Senks = outdata.Senks;
 Mback = outdata.Mback;
 Sback = outdata.Sback;
 Val   = outdata.Val;
 Obs   = outdata.Obs;
 CDF   = outdata.CDF;

 % find sties withe enough data
 I = find(Ndata<M);

 % ensure all sites have the same number of data points
 Mopen(:,I) = []; Mopen(M+1:end,:) = [];
 Menks(:,I) = []; Menks(M+1:end,:) = [];
 Senks(:,I) = []; Senks(M+1:end,:) = [];
 Mback(:,I) = []; Mback(M+1:end,:) = [];
 Sback(:,I) = []; Sback(M+1:end,:) = [];
   Val(:,I) = [];   Val(M+1:end,:) = [];
   Obs(:,I) = [];   Obs(M+1:end,:) = [];
   CDF(:,I) = [];   CDF(M+1:end,:) = [];

 % turn into single vectors over all sites
 X  = Mopen(:); % prior mean
 Y  = Obs(:);   % observation data
 C  = CDF(:);   % cdf-matched obs  
 Z  = Val(:);   % validation data   
 Xa = Menks(:); % posterior mean
 Qa = Senks(:); % posterior var
 Xb = Mback(:); % posterior mean
 Qb = Sback(:); % posterior var

 % number of samples
 N = length(X);
 assert(length(Y)==N)
 assert(length(Z)==N)
 assert(length(Xa)==N)
 assert(length(Qa)==N)
 assert(length(Xb)==N)
 assert(length(Qb)==N)

 % correction factor
 cf = 1/N;

 % screen report
 t = toc; fprintf('Finished loading data: %d \n',t);

%% *************************************************************************************
%% --- total information analysis using the raw data 
%% *************************************************************************************

 % screen report
 fprintf('Calculating point metrics ... \n'); tic;

 % calculate the point metrics
 [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(Xb,Y,Z,B,B,B); % in data
 Itot = (Ixz+Iyz-Ixyz)/Hz;                                  % total available
 Imod = Ixz/Hz;                                             % from model
 Iobs = Iyz/Hz;                                             % from observation
 Isyn = Ixyz/Hz;                                            % synergistic

 % calculate the point metrics
 [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(Xb,C,Z,B,B,B); % in data
 Itot_cdf = (Ixz+Iyz-Ixyz)/Hz;                              % total available
 Iobs_cdf = Iyz/Hz;                                         % from observation
 Isyn_cdf = Ixyz/Hz;                                        % synergistic
 Icdf_loss = (Itot-Itot_cdf)/Itot;

 [Ixz,Hx,Hz] = mutual_info(Xa,Z,B,B);                       % in enkf
 Ienks     = Ixz/Hz;                                        % from enkf
 Eenks     = Ienks/Itot;                                    % efficiency of enkf
 Eenks_cdf = Ienks/Itot_cdf;                                % efficiency of enkf

 % fraction of *obs* info used by the assimilation algortihm
 Eenks_obs     = (Ienks-Imod)/(Iobs-Isyn);
 Eenks_obs_cdf = (Ienks-Imod)/(Iobs_cdf-Isyn_cdf);

 % screen report
 obs_info  = [Imod,Iobs    ,Isyn    ,Itot    ,Ienks]
 cdf_info  = [Imod,Iobs_cdf,Isyn_cdf,Itot_cdf,Icdf_loss]
 tot_effy  = [Eenks,Eenks_obs,Eenks_cdf,Eenks_obs_cdf]

 % screen report
 t = toc; fprintf('Finished calculating data histograms: %d \n',t);

%% *************************************************************************************
%% --- calculate data histograms 
%% *************************************************************************************

 % screen report
 fprintf('Calculating data histograms ... \n'); tic;

 % data histograms
 [Pxyz,Pxy,Pxz,Pyz,Px,Py,Pz] = hist3(Xb,Y,Z,B,B,B);   

 % data posterior
 for x = 1:length(B)
  for y = 1:length(B)
   for z = 1:length(B)
    Pz_xy(x,y,z) = Pxyz(x,y,z)/Pxy(x,y);
%    if isnan(Pz_xy(x,y,z)); Pz_xy(x,y,z) = 0; end;
   end
   Pz_xy(x,y,:) = Pz_xy(x,y,:)/nansum(Pz_xy(x,y,:));
  end
 end

 % data prior and likelihood
 for x = 1:length(B)
  for z = 1:length(B)
   Pz_x(x,z) = Pxz(x,z)/Px(x);              
%   if isnan(Pz_x(x,z)); Pz_x(x,z) = 0; end;
   for y = 1:length(B)
    Py_xz(x,y,z) = Pxyz(x,y,z)/Pxz(x,z);
%    if isnan(Py_xz(x,y,z)); Py_xz(x,y,z) = 0; end;
   end
   Py_xz(x,:,z) = Py_xz(x,:,z)/nansum(Py_xz(x,:,z));
  end
  Pz_x(x,:) = Pz_x(x,:)/nansum(Pz_x(x,:));
 end

 % get rid of missing values 
% Py_xz(isnan(Py_xz)) = cf;
% Pz_xy(isnan(Pz_xy)) = cf;
% Pz_x (isnan(Pz_x) ) = cf;
% Py_xz(Py_xz<cf) = cf;
% Pz_xy(Pz_xy<cf) = cf;
% Pz_x (Pz_x<cf)  = cf;
% for x = 1:length(B)
%  for z = 1:length(B)
%   Py_xz(x,:,z) = Py_xz(x,:,z)/nansum(Py_xz(x,:,z));
%  end
%  for y = 1:length(B)
%   Pz_xy(x,:,z) = Pz_xy(x,y,:)/nansum(Pz_xy(x,y,:));
%  end
%  Pz_x(x,:) = Pz_x(x,:)/nansum(Pz_x(x,:));
% end

 % check that the prior, likelihood and posterior all agree
 for x = 1:length(B)
  for y = 1:length(B)
%   if any(isnan(Pz_xy(x,y,:))); Pz_xy(x,y,:) = 0; end 
   p1(x,y,:) = Pz_x(x,:).*squeeze(Py_xz(x,y,:))';
   p1(x,y,:) = p1(x,y,:)./nansum(p1(x,y,:));
  end
 end
 p2 = Pz_xy(:); p1 = p1(:);
 p1(isnan(p1)) = 0;
 assert(max(abs(p2-p1))<1e-10)

% % data histograms
% [Pxcz,Pxc,~,Pcz,~,Pc,~] = hist3(Xb,C,Z,B,B,B);   
%
% % data prior and likelihood
% for x = 1:length(B)
%  for z = 1:length(B)
%   for y = 1:length(B)
%    Pc_xz(x,y,z) = Pxcz(x,y,z)/Pxz(x,z);
%    if isnan(Pc_xz(x,y,z)); Pc_xz(x,y,z) = 0; end;
%   end
%   Pc_xz(x,:,z) = Pc_xz(x,:,z)/nansum(Pc_xz(x,:,z));
%  end
% end
%
% % get rid of missing values 
% Pc_xz(isnan(Pc_xz)) = 0;

 % screen report
 t = toc; fprintf('Finished calculating data histograms: %d \n',t); 

 % to avoid having to run the above loops every time
 fprintf('Saving data ... \n'); tic;
 save('./data/DecompProcessedData.mat');
 t = toc; fprintf('Finished saving data: %d \n',t);

else

 % load data
 fprintf('Loading data ... \n'); tic;
 load('./data/DecompProcessedData.mat');
 t = toc; fprintf('Finished loading data: %d \n',t);

 % report data stats
 obs_info
 cdf_info
 tot_effy

end

%% *************************************************************************************
%% --- subset data 
%% *************************************************************************************

% !! TEMPORARY !! Make runtimes shorter for coding/debugging
Nn = N;1e4; 
fprintf('\n Original # of data points: %d \n',N);
fprintf(' Reduced # of data points: %d \n \n',Nn);
I = randperm(N,Nn); N = Nn;
X = X(I); Y = Y(I); Z = Z(I); C = C(I);
Xa = Xa(I); Qa = Qa(I);
Xb = Xb(I); Qb = Qb(I);

%% *************************************************************************************
%% --- information analysis using KL-divs on the enkf distributions
%% *************************************************************************************

% divergences to be integrated (via MC) into an information measure
KLdata = zeros(N,1);   % total info: from cond. to prior (==Itot)
KLenkf = zeros(N,1);   % enkf loss: from anal. to cond.
KLbaye = zeros(N,1);   % same as above but with Bayes
KLlyss = zeros(N,1);   % same as above but with data prior
KLlysp = zeros(N,1);   % same as above but with data prior
KLlyps = zeros(N,1);   % same as above but with data prior
KLlypp = zeros(N,1);   % same as above but with data prior
%KLlcss = zeros(N,1);   % same as above but with data prior
%KLlcsp = zeros(N,1);   % same as above but with data prior
%KLlcps = zeros(N,1);   % same as above but with data prior
%KLlcpp = zeros(N,1);   % same as above but with data prior
KLpppp = zeros(N,1);   % same as above but with data lieklihood
KLppss = zeros(N,1);   % same as above but with data lieklihood
KLppsp = zeros(N,1);   % same as above but with data lieklihood
KLppps = zeros(N,1);   % same as above but with data lieklihood

% progress bar
bar = ceil(linspace(1,N,50));
fprintf('\n Calculating KL-div of EnKF analysis ... \n'); tic;
fprintf(repmat('.',1,50)); fprintf('\n');       % progress bar

% bin centers
Bc = B-0.5*(B(2)-B(1));

% loop through samples (this is an MC integration over divergences for an info metric)
for n = 1:N

 % progress bar
 if ~isempty(find(bar==n,1,'first')); fprintf('.'); end;

 % data prior
 pz = Pz;

% --------------
 % data conditional 
 x = find(Xb(n)>=B,1,'last');                      % index into data posterior
 y = find(Y(n) >=B,1,'last');                      % index into data posterior
 pd = squeeze(Pz_xy(x,y,:))';                      % data posterior 

 % data bayes
 pc = Pz_x(x,:).*squeeze(Py_xz(x,y,:))';           % grab posterior 
 pc = pc./nansum(pc);                                 % normalize
 assert(max(abs(pd-pc))<1e-10);                    % must match the one from hist3

% --------------
 % enkf analysis
 r = randn(N,1)*Qa(n) + Xa(n);                     % sample enkf posterior 
 pa = histc(r,B)';                                 % turn into a histogram
 pa = max(pa,1/N/1000);
 pa = pa./nansum(pa);                                 % normalize

 % enkf Bayes
 if Qb(n)<1e-2; Qb(n) = 1e-2; end;                 % Probably shouldn't have to do this.
 px = pdf('Normal',Bc,Xb(n),Qb(n));                % calculate prior 
 py = pdf('Normal',C(n),Bc,R);                     % calculate likelihood
 px = max(px,1/N/1000);
 py = max(py,1/N/1000);
 px = px./nansum(px); py = py./nansum(py);               % normalize
 pb = px.*py;                                      % calculate posterior
 pb = pb./nansum(pb);                                 % normalize

%% --------------
% % calculate the sample mean and varaince
% lmu  = nansum(repmat(Bc,[length(B),1]).*squeeze(Pc_xz(x,:,:)),2)';
% lsig = sqrt(sum(repmat(Bc.^2,[length(B),1]).*squeeze(Pc_xz(x,:,:)),2))';
%
% % likelihood (data prior, enkf likelihood): sample mean and sample variance
% pyss = pdf('Normal',C(n),lmu,lsig);               % calculate likelihood
% pyss(isnan(pyss)) = 0;                            % deal with nans
% pyss = pyss./nansum(pyss);                           % normalize likelihood 
% plssc = Pz_x(x,:).*pyss;                          % calculate posterior
% plssc = plssc./nansum(plssc);                        % normalize posterior
%
% % likelihood (data prior, enkf likelihood): sample mean and prescribed variance
% pysp = pdf('Normal',C(n),lmu,R);                  % calculate likelihood
% pysp(isnan(pysp)) = 0;                            % deal with nana
% pysp = pysp./nansum(pysp);                           % normalize likelihood 
% plspc = Pz_x(x,:).*pysp;                          % calculate posterior
% plspc = plspc./nansum(plspc);                        % normalize posterior
%
% % likelihood (data prior, enkf likelihood): prescribed mean and sample variance
% pyps = pdf('Normal',C(n),Bc,lsig);                % calculate likelihood 
% pyps(isnan(pyps)) = 0;                            % deal with grandma
% pyps = pyps./nansum(pyps);                           % normalize likelihood 
% plpsc = Pz_x(x,:).*pyps;                          % calculate posterior
% plpsc = plpsc./nansum(plpsc);                        % normalize posterior
%
% % likelihood (data prior, enkf likelihood): prescribed mean and prescribed variance
% pypp = pdf('Normal',C(n),Bc,R);                   % calculate likelihood
% pypp(isnan(pypp)) = 0;                            % deal with grandma
% pypp = pypp./nansum(pypp);                           % normalize likelihood 
% plppc = Pz_x(x,:).*pypp;                          % calculate posterior                   
% plppc = plppc./nansum(plppc);                        % normalize posterior

% --------------
 % calculate the sample mean and varaince
 lmu  = nansum(repmat(Bc,[length(B),1]).*squeeze(Py_xz(x,:,:)),2)';
 lsig = sqrt(nansum(repmat(Bc.^2,[length(B),1]).*squeeze(Py_xz(x,:,:)),2))';

 % likelihood (data prior, enkf likelihood): sample mean and sample variance
 pyss = max(1/N/1000,pdf('Normal',Y(n),lmu,lsig));               % calculate likelihood
 pyss(isnan(pyss)) = 1/N/1000;                            % deal with nans
 pyss = pyss./nansum(pyss);                           % normalize likelihood 
 plssy = Pz_x(x,:).*pyss;                          % calculate posterior
 plssy = plssy./nansum(plssy);                        % normalize posterior

 % likelihood (data prior, enkf likelihood): sample mean and prescribed variance
 pysp = max(1/N/1000,pdf('Normal',Y(n),lmu,R));                  % calculate likelihood
 pysp(isnan(pysp)) = 1/N/1000;                            % deal with nana
 pysp = pysp./nansum(pysp);                           % normalize likelihood 
 plspy = Pz_x(x,:).*pysp;                          % calculate posterior
 plspy = plspy./nansum(plspy);                        % normalize posterior

 % likelihood (data prior, enkf likelihood): prescribed mean and sample variance
 pyps = max(1/N/1000,pdf('Normal',Y(n),Bc,lsig));                % calculate likelihood
 pyps(isnan(pyps)) = 1/N/1000;                     % deal with grandma
 pyps = pyps./nansum(pyps);                           % normalize likelihood 
 plpsy = Pz_x(x,:).*pyps;                          % calculate posterior
 plpsy = plpsy./nansum(plpsy);                        % normalize posterior

 % likelihood (data prior, enkf likelihood): prescribed mean and prescribed variance
 pypp = max(1/N/1000,pdf('Normal',Y(n),Bc,R));     % calculate likelihood
 pypp(isnan(pypp)) = 1/N/1000;                     % deal with grandma
 pypp = pypp./nansum(pypp);                           % normalize likelihood 
 plppy = Pz_x(x,:).*pypp;                          % calculate posterior                   
 plppy = plppy./nansum(plppy);                        % normalize posterior

% --------------
 % calculate the sample mean and varaince
 pmu  = nansum(     Bc   .*squeeze(Pz_x(x,:)));
 psig = sqrt(nansum(Bc.^2.*squeeze(Pz_x(x,:))));

 % Model Distribution Error: switch prior (enkf prior, data likelihood)
 pppp = pdf('Normal',Bc,Xb(n),Qb(n));                % calculate prior 
 pppp = pppp.*squeeze(Py_xz(x,y,:))';                  % grab posterior 
 pppp = max(pppp,1/N/1000);
 pppp = pppp./nansum(pppp);                                 % normalize

 % prior (data likelihood, enkf prior): sample mean and variance
 pyss = max(1/N/1000,pdf('Normal',Bc,pmu,psig));               % calculate likelihood
 ppss = pyss.*squeeze(Py_xz(x,y,:))';
 ppss = max(ppss,1/N/1000);
 ppss = ppss/nansum(ppss);

 % prior (data likelihood, enkf prior): sample mean and variance
 pysp = max(1/N/1000,pdf('Normal',Bc,pmu,Qb(n)));               % calculate likelihood
 ppsp = pysp.*squeeze(Py_xz(x,y,:))';
 ppsp = max(ppsp,1/N/1000);
 ppsp = ppsp/nansum(ppsp);

 % prior (data likelihood, enkf prior): sample mean and variance
 pyps = max(1/N/1000,pdf('Normal',Bc,Xb(n),psig));               % calculate likelihood
 ppps = pyps.*squeeze(Py_xz(x,y,:))';
 ppps = max(ppps,1/N/1000);
 ppps = ppps/nansum(ppps);

% --------------
% --------------
% % deal with zero probabilities (this is not unambiguous and it affects the results)
% pa(pa<cf) = cf; pa = pa./nansum(pa);
% pb(pb<cf) = cf; pb = pb./nansum(pb);
% pc(pc<cf) = cf; pc = pc./nansum(pc);
% pd(pd<cf) = cf; pd = pd./nansum(pd);
% pz(pz<cf) = cf; pz = pz./nansum(pz);
% pp(pp<cf) = cf; pp = pp./nansum(pp);
%
% plssy(plssy<cf) = cf; plssy = plssy./nansum(plssy);
% plspy(plspy<cf) = cf; plspy = plspy./nansum(plspy);
% plpsy(plpsy<cf) = cf; plpsy = plpsy./nansum(plpsy);
% plppy(plppy<cf) = cf; plppy = plppy./nansum(plppy);

% plssc(plssc<cf) = cf; plssc = plssc./nansum(plssc);
% plspc(plspc<cf) = cf; plspc = plspc./nansum(plspc);
% plpsc(plpsc<cf) = cf; plpsc = plpsc./nansum(plpsc);
% plppc(plppc<cf) = cf; plppc = plppc./nansum(plppc);

 % force 0*log(~)=0
 I = find(pc>0);
 pa = pa(I); pb = pb(I); pc = pc(I); pd = pd(I); pz = pz(I); 
 plssy = plssy(I); plpsy = plpsy(I); plspy = plspy(I); plppy = plppy(I);
% plssc = plssc(I); plpsc = plpsc(I); plspc = plspc(I); plppc = plppc(I);
 pppp = pppp(I); ppss = ppss(I); ppsp = ppsp(I); ppps = ppps(I);

 % divergence integrations
 KLdata(n) = pc*log(pc./pz)';     % total info: from cond. to prior (==Itot)
 KLenkf(n) = pc*log(pc./pa)';     % enkf loss: from anal. to cond.
 KLbaye(n) = pc*log(pc./pb)';     % same as above but with Bayes

 KLlyss(n) = pc*log(pc./plssy)';  % same as above but with data prior
 KLlysp(n) = pc*log(pc./plspy)';  % same as above but with data prior
 KLlyps(n) = pc*log(pc./plpsy)';  % same as above but with data prior
 KLlypp(n) = pc*log(pc./plppy)';  % same as above but with data prior

% KLlcss(n) = pc*log(pc./plssc)';  % same as above but with data prior
% KLlcsp(n) = pc*log(pc./plspc)';  % same as above but with data prior
% KLlcps(n) = pc*log(pc./plpsc)';  % same as above but with data prior
% KLlcpp(n) = pc*log(pc./plppc)';  % same as above but with data prior

 KLpppp(n) = pc*log(pc./pppp)';    % same as above but with data lieklihood
 KLppss(n) = pc*log(pc./ppss)';
 KLppsp(n) = pc*log(pc./ppsp)';
 KLppps(n) = pc*log(pc./ppps)';

 if any(isnan([KLdata(n),KLenkf(n),KLbaye(n)... 
              ,KLlyss(n),KLlysp(n),KLlyps(n),KLlypp(n)...
              ,KLpppp(n),KLppss(n),KLppsp(n),KLppps(n)]))
%              ,KLlcss(n),KLlcsp(n),KLlcps(n),KLlcpp(n)...
  plot(pa,'b'); hold on; plot(pb,'c'); plot(pd,'g'); plot(pz,'k'); plot(plppy,'r'); plot(pppp,'m'); hold off; 
  keyboard
 end

end; fprintf('\n');

% screen report
t = toc; fprintf('Finished calculating divergences: %d \n',t);

% remove nanas
I = find(any(isnan(                         ...
              [KLdata,KLenkf,KLbaye         ...
              ,KLlyss,KLlysp,KLlyps,KLlypp  ...
              ,KLpppp,KLppss,KLppsp,KLppps  ...
              ]),2)); 
%              ,KLlcss,KLlcsp,KLlcps,KLlcpp  ...
fprintf('\n Number of grandmas: %d \n',length(I))

KLdata(I) = []; KLbaye(I) = []; 
KLlyss(I) = []; KLlysp(I) = []; KLlyps(I) = []; KLlypp(I) = [];
%KLlcss(I) = []; KLlcsp(I) = []; KLlcps(I) = []; KLlcpp(I) = [];
KLpppp(I) = []; KLppss(I) = []; KLppsp(I) = []; KLppps(I) = [];

% expected KL divergences
iKLdata = sum(KLdata)/N/Hz;
iKLenkf = sum(KLenkf)/N/Hz;
iKLbaye = sum(KLbaye)/N/Hz;

iKLlyss = sum(KLlyss)/N/Hz;
iKLlysp = sum(KLlysp)/N/Hz;
iKLlyps = sum(KLlyps)/N/Hz;
iKLlypp = sum(KLlypp)/N/Hz;

%iKLlcss = sum(KLlcss)/N/Hz;
%iKLlcsp = sum(KLlcsp)/N/Hz;
%iKLlcps = sum(KLlcps)/N/Hz;
%iKLlcpp = sum(KLlcpp)/N/Hz;

iKLpppp = sum(KLpppp)/N/Hz;
iKLppss = sum(KLppss)/N/Hz;
iKLppsp = sum(KLppsp)/N/Hz;
iKLppps = sum(KLppps)/N/Hz;

% screen report
info = [iKLdata]
filt = [iKLenkf,iKLbaye]
liky = [iKLlypp,iKLlyss,iKLlysp,iKLlyps]
%likc = [iKLlcss,iKLlcsp,iKLlcps,iKLlcpp]
prir = [iKLpppp,iKLppss,iKLppsp,iKLppps]

% -------------------------
% save
save('./data/DecompResults.mat');


