function [kfold_preds,stats,sens] = kfold_benchmark(inputs,targets,kFold,modelType,Xedges,Yedges,verbose)

try
    
% dimensions
nData = length(targets);
assert(nData == size(inputs,1));

% kfold partitioning
kIndexes = kfold_partitioning(nData,kFold);

% out-of-sample predictions
kfold_preds = zeros(size(targets))./0;

% Screen report: Initialize progress bar
if verbose
    fprintf('Training k-Fold models - %s ... \n',modelType);
    fprintf(strcat(repmat('.',[1,2*kFold+1]),'\n')); tic;
end

% standardize
[trainInputs,muInputs,sgInputs] = zscore(trainInputs);
testInputs = (testInputs-muInputs)./sgInputs;

[trainTargets,muTargets,sgTargets] = zscore(trainTargets);
testTargets = (testTargets-muTargets)./sgTargets;

% k-fold loop:
for k = 1:kFold
    
    % Screen report:
    if verbose; fprintf('..'); end
    
    % pull random sample
    iTest = kIndexes(k,:);
    iTrain = find(~ismember(1:nData,iTest));
    iTrain(isnan(iTrain)) = [];
    iTest(isnan(iTest)) = [];
    
    switch modelType
        case 'linearRegression'
            
            % train a benchmark
            model = regress(targets(iTrain),[ones(length(iTrain),1),inputs(iTrain,:)]);
            
            % make predictions
            predictions = [ones(nData,1),inputs] * model;

            %  native sensitivity analysis
            sens(:,k) = model(2:end);

        case 'neuralNetwork'
            
            % get model parameters
            modelParms  = set_ann_parms;
            
            % train a benchmark
            model = trainANN(inputs(iTrain,:),targets(iTrain),modelParms);
            
            % make predictions
            predictions = model(inputs')';
            
            % no native sensitivity analysis
            sens(:,k) = zeros(size(inputs,2),1)./0;
            
        case 'randomForest'
            
            % get model parameters
            modelParms  = set_tbg_parms;
            
            % train a benchmark
            model = trainTBG(inputs(iTrain,:),targets(iTrain),modelParms);
            
            % make predictions
            predictions = predict(model.RegressionEnsemble,inputs);

            %  native sensitivity analysis
            sens(:,k) = oobPermutedPredictorImportance(model.RegressionEnsemble);
            
        case 'gaussianProcess'
            
            % get model parameters
            modelParms  = set_gpr_parms;
            
            % train a benchmark
            model = trainGPR(inputs(iTrain,:),targets(iTrain),modelParms);
            
            % make predictions
            predictions = predict(model.RegressionGP,inputs);
            
            %  native sensitivity analysis
            sens(:,k) = -log(model.RegressionGP.KernelInformation.KernelParameters(1:end-1,1));
            
        otherwise
            
            % invalid model
            error('Invalid model type: %s',modelType);
            
    end % select model
    
    % unstandardize
    predictions = (predictions .* sgTargets) + muTargets;
    targets = (targets .* sgTargets) + muTargets;
    
    % stats over train/test periods
    stats.kfold.train.nse = nse(predictions(iTrain),targets(iTrain));
    stats.kfold.test.nse  = nse(predictions(iTest), targets(iTest));
    stats.kfold.train.rho = rho(predictions(iTrain),targets(iTrain));
    stats.kfold.test.rho  = rho(predictions(iTest), targets(iTest));
    stats.kfold.train.nmi = nmi(predictions(iTrain),targets(iTrain),Xedges,Yedges);
    stats.kfold.test.nmi  = nmi(predictions(iTest),targets(iTest),Xedges,Yedges);
    
    % pull only out-of-sample preds
    kfold_preds(iTest) = predictions(iTest);

    % re-standardize
    targets = (targets - muTargets) ./ sgTargets;
    
end % k-loop

% stats over train/test periods
stats.all.nse = nse(kfold_preds,targets);
stats.all.rho = rho(kfold_preds,targets);
stats.all.nmi = nmi(kfold_preds,targets,Xedges,Yedges);

% Screen report: Initialize progress bar
if verbose
    fprintf('. finished; NSE = %f; time = %f[s] \n',stats.all.nse,toc);
end

catch
    keyboard
end
