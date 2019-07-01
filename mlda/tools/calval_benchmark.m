function [testPreds,stats,sens] = calval_benchmark(...
    trainInputs,trainTargets,...
    testInputs,testTargets,...
    modelType,...
    Xedges,Yedges,...
    verbose)

% dimensions
nTrain = length(trainTargets);
nTest = length(testTargets);
assert(nTrain == size(trainInputs,1));
assert(nTest == size(testInputs,1));

% standardize
[trainInputs,muInputs,sgInputs] = zscore(trainInputs);
testInputs = (testInputs-muInputs)./sgInputs;

[trainTargets,muTargets,sgTargets] = zscore(trainTargets);
testTargets = (testTargets-muTargets)./sgTargets;

switch modelType
    case 'linearRegression'
        
        % train a benchmark
        model = regress(trainTargets,[ones(nTrain,1),trainInputs]);
        
        % make predictions
        try
        trainPreds = [ones(nTrain,1),trainInputs] * model;
        testPreds = [ones(nTest,1),testInputs] * model;
        catch
            keyboard
        end
        
        %  native sensitivity analysis
        sens = model(2:end);
        
    case 'neuralNetwork'
        
        % get model parameters
        modelParms  = set_ann_parms;
        
        % train a benchmark
        model = trainANN(trainInputs,trainTargets,modelParms);
        
        % make predictions
        trainPreds = model(trainInputs')';
        testPreds = model(testInputs')';
        
        % no native sensitivity analysis
        sens = zeros(size(trainInputs,2),1)./0;
        
    case 'randomForest'
        
        % get model parameters
        modelParms  = set_tbg_parms;
        
        % train a benchmark
        model = trainTBG(trainInputs,trainTargets,modelParms);
        
        % make predictions
        trainPreds = predict(model.RegressionEnsemble,trainInputs);
        testPreds = predict(model.RegressionEnsemble,testInputs);
        
        %  native sensitivity analysis
        sens = oobPermutedPredictorImportance(model.RegressionEnsemble);
        
    case 'gaussianProcess'
        
        % get model parameters
        modelParms  = set_gpr_parms;
        
        % train a benchmark
        model = trainGPR(trainInputs,trainTargets,modelParms);
        
        % make predictions
        trainPreds = predict(model.RegressionGP,trainInputs);
        testPreds = predict(model.RegressionGP,testInputs);
        
        %  native sensitivity analysis
        sens(:,k) = -log(model.RegressionGP.KernelInformation.KernelParameters(1:end-1,1));
        
    otherwise
        
        % invalid model
        error('Invalid model type: %s',modelType);
        
end % select model

% unstandardize
trainPreds = (trainPreds .* sgTargets) + muTargets;
testPreds = (testPreds .* sgTargets) + muTargets;
trainTargets = (trainTargets .* sgTargets) + muTargets;
testTargets = (testTargets .* sgTargets) + muTargets;

% stats over train/test periods
stats.train.nse = nse(trainPreds,trainTargets);
stats.test.nse  = nse(testPreds,testTargets);
stats.train.rho = rho(trainPreds,trainTargets);
stats.test.rho  = rho(testPreds,testTargets);
stats.train.nmi = nmi(trainPreds,trainTargets,Xedges,Yedges);
stats.test.nmi  = nmi(testPreds,testTargets,Xedges,Yedges);

% Screen report: Initialize progress bar
if verbose
    fprintf('. finished; NSE = %f; time = %f[s] \n',stats.test.nse,toc);
end
