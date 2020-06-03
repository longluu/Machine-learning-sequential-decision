%%%%%%%%%%%%% Run the full model fitting program %%%%%%%%%%%%%
clear 
nLoops = 30;
nCore = 1;
negLLHMinLoop = [];
fitParameterAll = cell(nLoops/nCore, nCore);
negLLH = NaN(nLoops/nCore, nCore);

%% Choose the parameter for the model fit
% Model type 1: Self-consistent Bayes
%            2: Optimal Bayes
%            3: Gamma distribution
% Optimization algorithm 1: Nelder-Mead simplex
%                        2: Sequential quadratic programming 
%                        3: Simulated annneaing
%                        4: Genetic algorithm
% Set starting point for optimization 0: chooose randomly
%                                     1: set starting point
% Use binned or unbinned likelihood: 0: unbinned
%                                    1: binned
modelType = 2;
optimizationAlgorithm = 1;
SetStartPoint = 0;
initialValueAll = [4.72      6.27      8.95           0.0000       22.24       5.52     1.00      2.50];
                %[stdNoiseInitial                lapseRateInitial  priorRangeInitial stdDelayInitial smoothFactorInitial stdMotor]
binnedLH = 0;
plotFitProgress = 0;
plotEstFit = 0;
subjectID = {'xfl', 'aj', 'zw', 'skb'};
expNumber = [2 3];
includeIncongruentTrials = ''; % empty if not include incongruent trials Incongruent
fixMotorNoise = 1;
fixLapseRate = 1;

for ll = 1 : length(subjectID)
    initialValue = initialValueAll;
    
    %% Extract the data
    [binaryDecision, percentCW, nTrialsPerCondition, estimateDataExp2, angleDiff, stdMotor] = dataForFitting(subjectID(ll), expNumber(1), includeIncongruentTrials);
    [~, ~, ~, estimateDataExp3, ~] = dataForFitting(subjectID(ll), expNumber(2), includeIncongruentTrials);
        
    save('dataAll', 'angleDiff', 'estimateDataExp2', 'estimateDataExp3',...
             'stdMotor', 'nTrialsPerCondition', 'percentCW')
    if fixMotorNoise
        initialValue(end) = mean(stdMotor);
    end

    %% Start the fitting
    currentFolder = pwd;
    fileNameRoot = 'FitResult';
    fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
    fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
    if isempty(subjectID{ll})
        subjectID{ll} = 'Average';
    end
    subjectName = ['Subject: ' subjectID{ll}];
    expID = ['Exp: ' num2str(expNumber)];
    fprintf(fileID, '%27s %27s \r\n', subjectName, expID); 
    fprintf(fileID,'%9s %8s %9s %27s %20s  %8s %12s %12s %12s \r\n', '-LLHTotal', '-LLH1', '-LLH2', ...
                    'stdNoiseLevel', 'LapseRate', 'Prior', 'MemNoise', 'Smoothness', 'MotorNoise');
    counter = 1;
    for kk = 1 : nLoops/nCore
        if nCore ==1
            [tempFitParameter, tempNegLLH] = modelFitBayesJoint_New(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, '', fixMotorNoise, includeIncongruentTrials, fixLapseRate);
            fitParameterAll{kk} = tempFitParameter;
            negLLH(kk) = tempNegLLH;
        else  
            parfor ii = 1 : nCore
                [tempFitParameter, tempNegLLH] = modelFitBayes_New(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, '', fixMotorNoise, includeIncongruentTrials);
                fitParameterAll{kk, ii} = tempFitParameter;
                negLLH(kk, ii) = tempNegLLH;
            end
        end

        % Save the result
        for ii = 1 : nCore
            fitParameter = fitParameterAll{kk, ii};
            fprintf(fileID, ['Iteration' '-' num2str(counter) ' \r\n']);
            fprintf(fileID,'%9.1f %9.4f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f\r\n', negLLH(kk, ii),...
                fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8));
            counter = counter + 1;
        end
    end

    % Print the best params
    negLLH = negLLH(:);
    [negLLHSort,indSort] = sort(negLLH, 'ascend')
    fitParameterAll = fitParameterAll(:);
    fitParameter = fitParameterAll{indSort(1)};
    fprintf(fileID, ['Best params' ' \r\n']);
    fprintf(fileID,'%9.1f %9.4f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f\r\n', negLLH(indSort(1)),...
        fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8));
    fclose(fileID);
end