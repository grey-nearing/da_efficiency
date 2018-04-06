clear all
close all
clc

% desired length of dataset
Nens = 50;
Nsites = 129; % for SCAN_Site_List.txt, 383 for List_Sites_SCAN_CRN_COSMOS.txt

% parameters of distributions over radiation, temp, and precip
R = [ 1.0,-0.5,-0.8;
     -0.5, 1.0, 0.5;
     -0.8, 0.5, 1.0];

s1 = 0.3;
s2 = 50;
s3 = 0.5;

%A = zeros(N);
%for n1 = 1:N
%  for n2 = max(1,n1-16):min(N,n1+16)
%    [n1,n2]
%    A(n1,n2) = exp(1)^-abs(n1-n2);
%  end
%end
%
%C = [A*R(1,1),A*R(1,2),A*R(1,3);
%     A*R(2,1),A*R(2,2),A*R(2,3);
%     A*R(3,1),A*R(3,2),A*R(3,3)];
%
%C = chol(C);

% loop through sites
for s = 1:Nsites

 % load true forcing
 try 
  fname = strcat('./site_data/forcing_',num2str(s),'.txt');
  fid = fopen(fname);
  U = fscanf(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
  fclose(fid);
  N = size(U,1)/12;
  U = reshape(U,[12,N])';
 catch
  continue
 end

 % temporarily hold date stamps
 Dates = U(:,1:4);
 U(:,1:4) = [];

 % loop through ensembles
 for u = 1:Nens

  fprintf('Site: %d of %d - Ensemble %d of %d \n',s,Nsites,u,Nens);

  Pert = randn(N,3)*chol(R);
%  Pert = randn(1,N*3)*C;
%  Pert = reshape(Pert,N,3);
%  corrcoef(Pert)
%  test1 = Pert(1:end-1,:);
%  test2 = Pert(2:end,:);
%  corrcoef(test1(:,1),test2(:,1))
%  corrcoef(test1(:,2),test2(:,2))
%  corrcoef(test1(:,3),test2(:,3))

  Uhat = U;
  Uhat(:,6) = max(0,Uhat(:,6).*Pert(:,1)*s1);
  Uhat(:,7) = max(0,Uhat(:,7)+Pert(:,2)*s2);
  Uhat(:,8) = max(0,Uhat(:,8).*Pert(:,3)*s3);

  Uhat = [Dates,Uhat];

  fname = strcat('./site_data/forcing_',num2str(s),'_',num2str(u),'.txt');
  save(fname,'Uhat','-ascii');

%  fid = fopen(fname,'w');
%  for t = 1:N
%   fprintf(fid,'%4.4d %2.2d %2.2d %2.2d %17.10f %17.10f %17.10f %17.10f %17.10f %17.10f %17.10f %17.10f \n',Uhat(t,:));
%  end
%  fclose(fid);

  %save(fname,'Uhat','-ascii');
 end % ensemble
end % sites 



