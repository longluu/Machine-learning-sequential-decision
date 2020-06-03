%%%%%%%%%%%%% Run the full model fitting program %%%%%%%%%%%%%
clear 
nLoops = 20;
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
modelType = 1;
optimizationAlgorithm = 1;
SetStartPoint = 0;
initialValue = [2.7355    4.9799    7.1056           0.0000     27   20    1   0.9678    2.6327];
                %[stdNoiseInitial                lapseRateInitial  priorRangeInitial stdDelayInitial smoothFactorInitial stdMotor]
binnedLH = 0;
plotFitProgress = 0;
plotEstFit = 0;
subjectID = {'ll'};
expNumber = [1 2];
includeIncongruentTrials = ''; % empty if not include incongruent trials Incongruent
fixMotorNoise = 1;
fixLapseRate = 1;

for tt = 1 : length(subjectID)
    %% Extract the data
    [~, percentCW1, nTrialsPerCondition1, estimateDataExp1, ~, ~] = dataForFitting(subjectID(tt), expNumber(1), includeIncongruentTrials);    
    [~, percentCW2, nTrialsPerCondition2, estimateDataExp2, angleDiff, stdMotor] = dataForFitting(subjectID(tt), expNumber(2), includeIncongruentTrials);
        
    save('dataAll', 'angleDiff', 'estimateDataExp1', 'estimateDataExp2', ...
             'stdMotor', 'nTrialsPerCondition1', 'percentCW1', 'nTrialsPerCondition2', 'percentCW2')
    if fixMotorNoise
        initialValue(end) = mean(stdMotor);
    end

    %% Start the fitting
    currentFolder = pwd;
    fileNameRoot = 'FitResult';
    fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
    fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
    if isempty(subjectID{tt})
        subjectID{tt} = 'Average';
    end
    subjectName = ['Subject: ' subjectID{tt}];
    expID = ['Exp: ' num2str(expNumber)];
    fprintf(fileID, '%27s %27s \r\n', subjectName, expID); 
    fprintf(fileID,'%9s %8s %9s %27s %20s  %8s %8s %12s %12s %12s \r\n', '-LLHTotal', '-LLH1', '-LLH2', ...
                    'stdNoiseLevel', 'LapseRate', 'Prior1', 'Prior2', 'MemNoise', 'Smoothness', 'MotorNoise');
    counter = 1;
    for kk = 1 : nLoops/nCore
        if nCore ==1
            [tempFitParameter, tempNegLLH] = modelFitBayesJoint12(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, '', fixMotorNoise, includeIncongruentTrials, fixLapseRate);
            fitParameterAll{kk} = tempFitParameter;
            negLLH(kk) = tempNegLLH;
        else  
        end

        % Save the result
        for ii = 1 : nCore
            fitParameter = fitParameterAll{kk, ii};
            fprintf(fileID, ['Iteration' '-' num2str(counter) ' \r\n']);
            fprintf(fileID,'%9.1f %9.4f %9.4f %9.4f %16.4f  %10.4f %10.4f %10.4f %8.4f %9.4f\r\n', negLLH(kk, ii),...
                fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8), fitParameter(9));
            counter = counter + 1;
        end
    end

    % Print the best params
    negLLH = negLLH(:);
    [negLLHSort,indSort] = sort(negLLH, 'ascend')
    fitParameterAll = fitParameterAll(:);
    fitParameter = fitParameterAll{indSort(1)};
    fprintf(fileID, ['Best params' ' \r\n']);
    fprintf(fileID,'%9.1f %9.4f %9.4f %9.4f %16.4f  %10.4f %10.4f %10.4f %8.4f %9.4f\r\n', negLLH(indSort(1)),...
        fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8), fitParameter(9));
    fclose(fileID);
end