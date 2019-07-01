function [Ixy,Hx,Hy] = knn_info(X,Y,k)

% which algorithm to use?
algo = 1; % [1,2];

% dimensions of input samples
[Nx,Dx] = size(X);
[Ny,Dy] = size(X);
assert(Ny == Nx); N = Nx;

% init storages of neighborhood counts
ex = zeros(N,1)./0;
ey = zeros(N,1)./0;
Nx = zeros(N,1)./0;
Ny = zeros(N,1)./0;

% choose a distance metric - ENTROPIES DON'T WORK FOR NON-EUCLIDEAN NORMS
distanceMetric = 'euclidean';
% distanceMetric = 'minkowski';
% distanceMetric = 'squaredeuclidean';
% distanceMetric = 'mahalanobis';
% distanceMetric = 'cityblock';
% distanceMetric = 'chebychev';
% distanceMetric = 'cosine';
% distanceMetric = 'hamming';

% compute distances according to the chosen metric
% distances = squareform(pdist([X,Y],distanceMetric));
distancesX = squareform(pdist(X,distanceMetric));
distancesY = squareform(pdist(Y,distanceMetric));
distances = max(distancesX,distancesY);

% neighborhood census
for n = 1:N
    
    % find kth smallest distance, not counting the identity
    [~,distSort] = sort(distances(n,:));
    kthSmallest = distSort(k+1);
    
    % neighborhood radius as maximum norm under r.v. partition
    exi = distancesX(n,kthSmallest);
    eyi = distancesY(n,kthSmallest);
    if algo == 1
        ep = max(exi,eyi);
    end
    
    % neighborhood counts
    if algo == 1
        Nx(n) = sum(distancesX(n,:) < ep);
        Ny(n) = sum(distancesY(n,:) < ep);
    elseif algo == 2
        Nx(n) = sum(distancesX(n,:) <= exi);
        Ny(n) = sum(distancesY(n,:) <= eyi);
    end
    
    % find kth smallest distance along single partitions
    [~,distSortX] = sort(distancesX(n,:));
    kthSmallestX = distSortX(k+1);
    ex(n) = distancesX(n,kthSmallestX);
    
    [~,distSortY] = sort(distancesY(n,:));
    kthSmallestY = distSortY(k+1);
    ey(n) = distancesY(n,kthSmallestY);
    
end % n-loop

% mutual information calculation
if algo == 1
    Ixy = psi(k) - mean(psi(Nx)+psi(Ny)) + psi(N);
elseif algo == 2
    Ixy = psi(k) - 1/k  - mean(psi(Nx-1)+psi(Ny-1)) + psi(N);
end

% entropy calculations (these only work for the euclidean norm
cx = pi^(Dx/2) / gamma(Dx/2+1) / 2^Dx; % unit sphere in x-space
cy = pi^(Dy/2) / gamma(Dy/2+1) / 2^Dy; % unit sphere in y-space
Hx = -psi(k) + psi(N) + log(cx) + Dx * mean(log(ex*2));
Hy = -psi(k) + psi(N) + log(cy) + Dy * mean(log(ey*2));







