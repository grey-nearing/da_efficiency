function kIndexes = kfold_partitioning(N,K)
    
extraSamples = ceil(N/K)*K - N;
randIndexes = randperm(N+extraSamples);
randIndexes(randIndexes > N) = 0/0;
kIndexes = reshape(randIndexes,[K,ceil(N/K)]);


